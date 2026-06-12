#!/bin/bash
#SBATCH --job-name=fargo_corner_val
#SBATCH --account=teyssier
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:a100:4
#SBATCH --mem=16G
#SBATCH --time=00:20:00
#SBATCH --output=val/validate_corners_%j.out

# Validate the CUDA-aware MPI tag fix for nested-mesh ghost exchange.
#
# The pre-fix GPU path used comm->facedim (only 0/1/2) as the MPI tag. When two
# inter-CPU communicators between the same CPU pair share a (facedim, level) --
# which happens e.g. with a periodic dimension split across 2 ranks, where a
# rank's regular neighbour and its periodic-wrap neighbour are the SAME peer --
# both messages carry the same tag and MPI can only tell them apart by posting
# order. Mismatched buffers then corrupt corner/edge ghost cells. The fix gives
# each communicator a unique comm->tag (as the CPU path already does).
#
# To force the collision we enable PeriodicZ (dynfric's nested levels are all
# centred in Z), which makes the per-rank collision multiplicity == 2
# (confirmed with the FARGO_TAGDIAG diagnostic on the CPU build).
#
# Binaries (identical except the MPI tag in J_comm_device.c):
#   fargo3d_gpu_fixed  -> unique comm->tag   (the fix)
#   fargo3d_gpu_buggy  -> comm->facedim tag  (pre-fix behaviour)
# Reference: val/cpu_pz  (CPU build, same PeriodicZ 4-rank decomposition).

set -u
cd "$SLURM_SUBMIT_DIR"

source /usr/share/Modules/init/bash 2>/dev/null
module purge
module load cudatoolkit/12.6 openmpi/cuda-12.6/gcc/4.1.6

echo "=================================================================="
echo "Host: $(hostname)"
echo "CUDA_VISIBLE_DEVICES=${CUDA_VISIBLE_DEVICES:-unset}"
nvidia-smi --query-gpu=index,name,compute_cap --format=csv 2>&1 | head
echo "=================================================================="

PAR=setups/dynfric/dynfric.par
RUN="mpirun -np 4 --bind-to none"
OVERRIDES=(-o "PeriodicZ=yes")

for tag in fixed buggy; do
  bin=./fargo3d_gpu_${tag}
  out=val/gpu_pz_${tag}
  echo ""
  echo ">>> Running GPU (${tag}) 4-rank MPICUDA, PeriodicZ -> ${out}"
  rm -rf "$out"; mkdir -p "$out"
  $RUN "$bin" -k "${OVERRIDES[@]}" -o "OutputDir=${out}" "$PAR" > val/run_pz_${tag}.log 2>&1
  echo "exit=$? ; full log -> val/run_pz_${tag}.log ; last 8 lines:"
  tail -8 val/run_pz_${tag}.log
done

echo ""
echo "##################################################################"
echo "## TEST 1 (decisive): GPU fixed vs GPU buggy  -- identical binaries"
echo "## except the MPI tag, so ANY ghost difference is caused by the tag."
echo "##################################################################"
python scripts/compare_outputs.py val/gpu_pz_fixed val/gpu_pz_buggy --snap 2
fvb_rc=$?

echo ""
echo "##################################################################"
echo "## TEST 2: CPU reference vs GPU FIXED   (expect AGREE)"
echo "##################################################################"
python scripts/compare_outputs.py val/cpu_pz val/gpu_pz_fixed --snap 2
fixed_rc=$?

echo ""
echo "##################################################################"
echo "## TEST 3: CPU reference vs GPU BUGGY   (expect CORNER MISMATCH)"
echo "##################################################################"
python scripts/compare_outputs.py val/cpu_pz val/gpu_pz_buggy --snap 2
buggy_rc=$?

echo ""
echo "##################################################################"
echo "fixed-vs-buggy exit code : $fvb_rc    (1 = they differ -> tag matters)"
echo "cpu-vs-fixed   exit code : $fixed_rc  (0 = fixed agrees with CPU)"
echo "cpu-vs-buggy   exit code : $buggy_rc  (1 = buggy corrupts corners)"
if [ "$fixed_rc" -eq 0 ] && [ "$buggy_rc" -eq 1 ] && [ "$fvb_rc" -eq 1 ]; then
  echo "VALIDATION PASSED: the unique-tag fix corrects corner/edge ghost exchange."
else
  echo "VALIDATION INCONCLUSIVE: inspect the per-field diffs above."
fi
echo "##################################################################"
