# fargo3d_nm

FARGO3D + Nested Meshes.

This repository contains the FARGO3D hydrodynamics code (legacy version) with a
JUPITER-style nested-mesh implementation, plus fixes and tooling for running the
nested-mesh communications on GPUs with CUDA-aware MPI:

- Unique MPI tags in the GPU ghost-exchange path (`src/J_comm_device.c`),
  matching the CPU path and making message matching independent of posting order.
- Fix for a GPU crash (`invalid pitch argument in cudaMemcpy2D`) when copying
  the overlap mask to the device in 2D setups without the X dimension
  (`src/J_overlap.c`).
- Optional MPI tag-collision diagnostic, enabled with the `FARGO_TAGDIAG`
  environment variable (`src/J_comm_init.c`).
- Build scripts ported from Python 2 to Python 3 (`scripts/`).
- `DELLA` target platform (NVIDIA A100, `sm_80`, CUDA-aware OpenMPI) in
  `src/makefile`.
- A nested-mesh test suite (`test_suite/nm_*.py`, run with `make testnm`)
  checking decomposition invariance, restarts, posterior refinement, mass
  budget, 3D corner exchanges, and GPU/CUDA-aware MPI runs on the physical
  cells of every refinement level.

## Documentation

See [docs/NESTED_MESHES.md](docs/NESTED_MESHES.md) for how to set up and use
mesh refinement (the `#GRIDINFO` directive, `RELNESTING`, build options,
running on CPU/GPU, restart and posterior refinement, output format, and
the test suite).

## Building

Example build of a nested-mesh setup with GPU and CUDA-aware MPI:

```bash
make SETUP=dynfric GPU=1 PARALLEL=1 MPICUDA=1 FARGO_ARCH=DELLA
```

See `README` for the original FARGO3D documentation.

## License

GPL-3.0 — see `license.txt`.
