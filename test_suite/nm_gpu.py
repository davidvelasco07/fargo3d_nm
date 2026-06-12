#!/usr/bin/env python3
"""GPU tier (requires CUDA devices; run via scripts/run_nm_gpu_tests.slurm).

Checks, on the 2D nested setup with PeriodicZ (the configuration with the
highest MPI-tag collision multiplicity):
  1. gpu np1 vs cpu np4   : same physics within floating-point tolerance;
  2. gpu np4 (staged MPI)  vs gpu np1 : bit-exact decomposition invariance;
  3. gpu np4 (CUDA-aware MPI, MPICUDA=1) vs gpu np1 : bit-exact, plus the
     FARGO_TAGDIAG diagnostic must confirm tag-collision multiplicity > 1,
     i.e. the unique-tag fix is actually being exercised.

Compute nodes have the MPI runtime but no compiler wrappers (mpicc), so
binaries must be built beforehand on the login node:
    python3 test_suite/nm_gpu.py --build-only
populates the binary cache that the SLURM job then uses."""
import os
import re
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.realpath(__file__)),
                                "..", "scripts"))
import nm_test as T

PAR = "setups/nmtest/nmtest.par"
OVR = {"PeriodicZ": "yes"}
# GPU FMA/contraction noise after 2 steps: relative bound plus an absolute
# floor for near-zero velocities (measured max abs diff is ~1e-15).
GPU_CPU_TOL_REL = 1e-10
GPU_CPU_TOL_ABS = 1e-13

t = T.NestedMeshTest("nm_gpu")
try:
    cpu = t.build("nmtest")
    gpu = t.build("nmtest", GPU=1)
    gpucuda = t.build("nmtest", GPU=1, MPICUDA=1)
    if "--build-only" in sys.argv:
        print("binaries cached; submit scripts/run_nm_gpu_tests.slurm")
        sys.exit(0)

    cpu4 = t.run(cpu, PAR, "cpu_np4", np=4, overrides=OVR)
    gpu1 = t.run(gpu, PAR, "gpu_np1", np=1, overrides=OVR)
    t.compare(cpu4, gpu1, snap=2, tol_rel=GPU_CPU_TOL_REL,
              tol_abs=GPU_CPU_TOL_ABS, label="gpu np1 vs cpu np4 (tolerance)")

    gpu4 = t.run(gpu, PAR, "gpu_np4", np=4, overrides=OVR)
    t.compare(gpu1, gpu4, snap=2, label="gpu np4 (staged) vs gpu np1")

    gc4 = t.run(gpucuda, PAR, "gpucuda_np4", np=4, overrides=OVR,
                env={"FARGO_TAGDIAG": "1"})
    t.compare(gpu1, gc4, snap=2, label="gpu np4 (CUDA-aware) vs gpu np1")

    log = t.runlog("gpucuda_np4")
    mults = [int(m) for m in
             re.findall(r"worst multiplicity = (\d+)", log)]
    t.check("tag-collision condition exercised (multiplicity > 1)",
            mults and max(mults) > 1,
            "TAGDIAG multiplicities: %s" % mults)
except T.TestFailure as exc:
    t.check("execution", False, str(exc))
t.finish()
