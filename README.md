# fargo3d_nm

FARGO3D + Nested Meshes.

This repository contains the FARGO3D hydrodynamics code (legacy version) with a
JUPITER-style nested-mesh implementation, plus fixes and tooling for running the
nested-mesh communications on GPUs with CUDA-aware MPI:

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
