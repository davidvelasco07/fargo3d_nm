# Nested meshes in fargo3d_nm

This document describes the nested-mesh (JUPITER-style) machinery in this
repository: how to define refined patches, how to build and run a nested
setup on CPU and GPU, what the code does internally, and how outputs are
organized.

## Overview

A simulation consists of a base mesh (level 0, defined by the usual
`Nx/Ny/Nz`, `Xmin/Xmax`, ... parameters) plus any number of refined
rectangular patches at levels 1, 2, 3, ... Each level refines the previous
one by a **factor 2 per dimension** (refinement of individual dimensions can
be disabled, see `INHIBREFDIM` below). Patches of level `l+1` must be fully
contained in a patch of level `l`.

At runtime:

- Every patch is split among MPI ranks like an ordinary FARGO3D mesh.
- Ghost cells of each patch are filled by communicators built at startup:
  same-level exchange (`GHOST`), coarse-to-fine interpolation (also `GHOST`,
  level `l-1` to `l`), fine-to-coarse restriction (`MEAN`), and conservative
  flux correction at coarse-fine interfaces (`FLUX`).
- Levels are advanced recursively with subcycling: a fine level may take 2
  substeps per coarse step (see "Time stepping" below).

## Setting up refinement

### The `#GRIDINFO` directive

Refined patches are declared directly in the parameter (`.par`) file, after
a line starting with `#GRIDINFO`. Every subsequent non-comment line defines
one patch with `2*NDIM + 1` numbers:

```text
<min coords (NDIM values)>  <max coords (NDIM values)>  <level>
```

where `NDIM` is the number of active dimensions selected in the `.opt` file
(`-DX`, `-DY`, `-DZ`). Lines starting with `#` are comments, so patches can
be toggled by commenting them out.

Example for a 2D (Y,Z) setup, [setups/dynfric/dynfric.par](../setups/dynfric/dynfric.par):

```text
#GRIDINFO
-1.0 -0.5   0.0 0.5   1
-1.0 -0.25 -0.5 0.25  2
-1.0 -0.125 -0.75 0.125 3
```

Example for a 3D setup, [setups/td3d/td3d.par](../setups/td3d/td3d.par):

```text
#GRIDINFO
-175.0 -175.0 -175.0  175.0 175.0 175.0  1
-87.5  -87.5  -87.5   87.5  87.5  87.5   2
```

### Absolute vs. relative coordinates (`RELNESTING`)

The interpretation of the patch coordinates depends on the `.opt` flag
`-DRELNESTING`:

- **Without** `RELNESTING`: coordinates are absolute physical coordinates,
  in the same units as `Xmin/Xmax`, etc.
- **With** `RELNESTING`: coordinates are relative to the full domain, where
  the domain center maps to `0` and the half-extent maps to `1`. A value
  `c` in dimension `d` means the physical position

  ```text
  x = c * (Dmax - Dmin)/2 + (Dmax + Dmin)/2
  ```

  so `-1 ... 1` spans the whole domain. In the `dynfric` example above
  (Y in `[0, 10]`, Z in `[-10, 10]`), the level-1 line `-1.0 -0.5 0.0 0.5`
  means Y in `[0, 5]` and Z in `[-5, 5]`. The advantage is that the
  refinement layout is independent of the physical domain size.

The parsed (always absolute) patch layout is written at startup to
`<OutputDir>/nested_meshes.txt` — check it to confirm the geometry you
intended.

### Rules and limits

- Patch bounds are snapped to cell boundaries of the parent level; the
  actual extent is what appears in `nested_meshes.txt` and the output
  descriptors.
- A patch of level `l` must lie inside a level `l-1` patch (the code prints
  "Ill-ordered refinement information" and aborts otherwise). List patches
  in increasing level order.
- At most `MAXGRIDS` patches (1000, see [src/J_def.h](../src/J_def.h)).
- The maximum level actually used can be capped at runtime with the
  `-j <level>` command-line flag; `#GRIDINFO` lines above that level are
  discarded. This lets you keep one parameter file and run coarse previews
  with e.g. `-j 1`.
- Refinement of a given dimension can be disabled entirely (patches then
  refine only the other dimensions) via the `INHIBREFDIM1/2/3` switches;
  dimensions absent from the build are inhibited automatically.
- On periodic dimensions a patch may not be larger than the period.

### Relevant `.opt` flags

| Flag | Meaning |
|---|---|
| `-DRELNESTING` | `#GRIDINFO` coordinates are relative (see above). |
| `-DCOMM=INTERNAL` / `EXTERNAL` / `ASYMETRIC` | Convention used for staggered velocity components at coarse-fine interfaces. The setups in this repo use `INTERNAL`. |
| `-DNONMONOTONIC` | Coarse-to-fine ghost interpolation uses plain (unlimited) linear interpolation instead of van Leer slope limiting. |
| `-DNOGHOSTSX` | No ghost zones in X (e.g. axisymmetric 2D setups). |

See [setups/dynfric/dynfric.opt](../setups/dynfric/dynfric.opt) for a complete
2D nested example and `setups/td3d`, `setups/p3dthd`, `setups/heatforce` for
3D ones.

## Building

The build scripts are Python 3. A nested-mesh setup is built like any other
FARGO3D setup:

```bash
# CPU, serial
make SETUP=dynfric

# CPU, MPI-parallel
make SETUP=dynfric PARALLEL=1

# GPU + CUDA-aware MPI (device buffers passed directly to MPI)
make SETUP=dynfric GPU=1 PARALLEL=1 MPICUDA=1 FARGO_ARCH=DELLA
```

`FARGO_ARCH=DELLA` selects the Princeton Della platform block in
[src/makefile](../src/makefile) (A100, `sm_80`, CUDA-aware OpenMPI); add an
analogous block for your own cluster. Load the matching modules first and
export `CUDA=$CUDA_HOME`, e.g.:

```bash
module load cudatoolkit/12.6 openmpi/cuda-12.6/gcc/4.1.6
export CUDA=$CUDA_HOME
```

By default outputs contain only the physical (interior) cells of each
patch. `FULLDEBUG=1` builds with `-g` and additionally enables
`WRITEGHOSTS`, which folds the ghost zones into the output files — a
low-level debugging aid only; such runs must use the `-k` flag and their
outputs cannot be merged or restarted from.

## Running

```bash
# 4 MPI ranks
mpirun -np 4 ./fargo3d setups/dynfric/dynfric.par

# Override parameters on the command line
mpirun -np 4 ./fargo3d -o "OutputDir=run1" -o "PeriodicZ=yes" setups/dynfric/dynfric.par
```

Command-line flags most relevant to nested runs (full list: run without
arguments):

| Flag | Meaning |
|---|---|
| `-m` / `-k` | Merge / do **not** merge per-process output files. `WRITEGHOSTS` builds require `-k`. |
| `-o "Param=value"` | Override any `.par` parameter. |
| `-j <level>` | Cap the maximum refinement level (discard deeper `#GRIDINFO` lines). |
| `-s <n>` | Restart from output `n` (split, i.e. per-process outputs). |
| `-S <n>` | Restart from output `n` (merged outputs). |
| `-r "<out> <grid> <xmin> <ymin> <zmin> <xmax> <ymax> <zmax>"` | Posterior refinement: when restarting, insert a new sub-patch refining grid number `<grid>` of output `<out>` over the given physical extent (clamped to the parent patch). |
| `-0` | Write the initial (or restart) output and exit. |
| `-D <n>` | Manually select GPU device `n`. |
| `-C` | Force all functions on the CPU in a GPU build. |

Note that on **restart** the patch layout is read from the chosen output's
descriptor files, not from `#GRIDINFO` — editing `#GRIDINFO` does not affect
restarted runs; use `-r` to add refinement a posteriori.

### GPU runs with CUDA-aware MPI

With `GPU=1 PARALLEL=1 MPICUDA=1`, all nested-mesh ghost exchanges pass GPU
device pointers directly to MPI ("Communications are done directly on GPU"
is printed at startup). Requirements:

- A genuinely CUDA-aware MPI library (e.g. OpenMPI built with CUDA support).
- One GPU per rank; ranks are bound to devices using the local-rank
  environment variable configured per platform in the makefile
  (`OMPI_COMM_WORLD_LOCAL_RANK` for OpenMPI).

Each communicator carries a unique MPI tag (`comm->tag`), on both the CPU
and GPU paths, so message matching does not depend on posting order. To
audit a given setup/decomposition for messages that would have shared a tag
under the legacy facedim-based scheme, export `FARGO_TAGDIAG=1`; each rank
then reports its tag-collision multiplicity at startup.

## Time stepping and subcycling

Levels are advanced by recursive iteration ([src/J_iter.c](../src/J_iter.c)):
each coarse step advances the finer level by `TimeStepRatio` substeps (1 or
2), then applies the conservative flux correction (`FLUX` communicators),
restricts the fine solution onto the coarse patch (`MEAN`), refills
boundaries, and interpolates coarse data to the fine ghost zones. The CFL
condition is evaluated globally across levels. Subcycling is enabled
(ratio 2 per level) by default.

## Outputs

Each snapshot `N` produces a directory `<OutputDir>/outputNNNNN/`
containing:

- `DescriptorN.dat` — the patch/CPU layout (levels, sizes, ghost count,
  fields) needed to interpret the data files and to restart.
- One file per field, per CPU patch:
  `<fluid><field><N>_<patch>_<level>.dat`, e.g. `gasdensity2_9_1.dat` is
  the density of patch 9 (a level-1 patch) at snapshot 2. Velocity files
  contain `NDIM` contiguous components. Files are raw double-precision
  arrays holding the physical (interior) cells of the patch — outputs never
  contain ghost cells unless the build enables the `WRITEGHOSTS` debug
  option (`FULLDEBUG=1`), which is meant for low-level debugging only.

Note that nested-mesh outputs are always split per CPU patch, so two runs
with different rank counts cannot be compared file by file; use
`scripts/nm_compare.py` (below), which reassembles each refinement level
onto its full canvas before comparing.

In addition, the run directory contains `nested_meshes.txt` (the parsed
refinement layout), `variables.par`, `info.log`, and the monitoring time
series (`mass.dat`, `force.dat`, and the `*_all_levels.dat` variants summed
over all levels).

## Test suite

The nested-mesh test suite works exclusively on physical (interior) cells
of standard builds — ghost cells are never written or compared. Errors in
the ghost exchange contaminate interior cells within a couple of timesteps,
so decomposition-invariance checks on the interiors catch them without any
debug output options.

### CPU tier

```bash
make testnm        # runs all CPU-tier tests
make testlist      # lists every available test (legacy + nm_*)
```

Run from the repository root with `mpirun` in the `PATH` (override the
launcher with `NM_MPIRUN` if needed). The tests, in
[test_suite/](../test_suite/):

| Test | What it checks |
|------|----------------|
| `nm_serial_vs_mpi` | 1-rank vs 4-rank 2D run, bit-exact interiors on all levels |
| `nm_periodic` | same with `PeriodicZ`, the configuration with the highest MPI-tag collision multiplicity |
| `nm_restart` | restart reproduces an uninterrupted run to machine precision and is bit-exact across decompositions |
| `nm_subpatch` | posterior refinement (`-r`) adds a patch, the result restarts and keeps evolving |
| `nm_conservation` | base-level mass budget stays within a calibrated bound and is bit-identical across decompositions |
| `nm_3d` | 1-rank vs 8-rank 3D run (2x2x2, corner-type rank adjacency), bit-exact |

### GPU tier

```bash
sbatch scripts/run_nm_gpu_tests.slurm   # or: python3 test_suite/nm_gpu.py on a GPU node
```

`test_suite/nm_gpu.py` compares a GPU run against a CPU reference (within
floating-point tolerance), then checks that 4-rank GPU runs — both with
staged MPI and with CUDA-aware MPI (`MPICUDA=1`) — are bit-identical to the
1-rank GPU run, and uses the `FARGO_TAGDIAG` diagnostic to confirm the
tag-collision condition is actually exercised.

### Supporting tools

- [scripts/nm_test.py](../scripts/nm_test.py) — the test framework
  (cached builds, `mpirun` invocation, comparisons, reporting).
- [scripts/nm_compare.py](../scripts/nm_compare.py) — standalone comparer:
  reassembles every refinement level from the per-patch files of two runs
  and reports per-level/per-field maximum differences. Works across
  different MPI decompositions:

  ```bash
  python3 scripts/nm_compare.py REF_DIR TEST_DIR --snap 2 [--tol-rel X --tol-abs Y]
  ```

- [scripts/compare_outputs.py](../scripts/compare_outputs.py) — low-level
  debug tool for `WRITEGHOSTS` builds: reports differences separately for
  interior, edge-ghost and corner-ghost regions of every patch file
  (requires `--ngh` since `WRITEGHOSTS` descriptors report a zero ghost
  width). Used by
  [scripts/validate_corners.slurm](../scripts/validate_corners.slurm), the
  original GPU ghost-exchange validation job. Not used by the test suite.

## Troubleshooting

- "Cannot merge outputs when dumping ghost values": the build has
  `WRITEGHOSTS` (e.g. `FULLDEBUG=1`); run with `-k`.
- "Ill-ordered refinement information": a `#GRIDINFO` patch is not
  contained in a coarser patch, or levels are not listed in increasing
  order.
- "Too many grids specified": increase `MAXGRIDS` in
  [src/J_def.h](../src/J_def.h) and recompile.
- A startup check (`CheckBC`) verifies that every ghost zone of every patch
  is covered by a communicator or a physical boundary condition and aborts
  with a diagnostic if there is a hole.
