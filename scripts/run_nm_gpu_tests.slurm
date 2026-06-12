#!/bin/bash
#SBATCH --job-name=fargo_nm_gpu
#SBATCH --account=teyssier
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:a100:4
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --output=test_out/nm_gpu_%j.out

# GPU tier of the nested-mesh test suite (test_suite/nm_gpu.py).
# Compute nodes have no compiler wrappers (mpicc), so build the binaries on
# the login node first, then submit from the repository root:
#   python3 test_suite/nm_gpu.py --build-only
#   sbatch scripts/run_nm_gpu_tests.slurm

set -u
cd "$SLURM_SUBMIT_DIR"
mkdir -p test_out

source /usr/share/Modules/init/bash 2>/dev/null
module purge
module load cudatoolkit/12.6 openmpi/cuda-12.6/gcc/4.1.6

export FARGO_ARCH=DELLA
export CUDA=$CUDA_HOME

echo "Host: $(hostname)"
echo "CUDA_VISIBLE_DEVICES=${CUDA_VISIBLE_DEVICES:-unset}"
nvidia-smi -L || true
echo "=================================================================="

python3 test_suite/nm_gpu.py
exit $?
