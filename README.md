# FARGO3D #

#### A versatile MULTIFLUID HD/MHD code that runs on clusters of CPUs or GPUs, with special emphasis on protoplanetary disks. 

This fork adds support for **statically nested meshes** (mesh refinement), on CPU and GPU, including multi-rank runs with CUDA-aware MPI. See [Nested meshes](#nested-meshes) below.

------------------------

##### Website: [fargo.in2p3.fr](http://fargo.in2p3.fr)

##### [Documentation](https://fargo3d.bitbucket.io/)

------------------------

### First run

#### Sequential CPU

``` 
make SETUP=fargo PARALLEL=0 GPU=0
./fargo3d in/fargo.par
```

#### Parallel CPU

```
make SETUP=fargo PARALLEL=1 GPU=0
mpirun -np 4 ./fargo3d in/fargo.par
```

#### Sequential GPU

```
make SETUP=fargo PARALLEL=0 GPU=1
./fargo3d in/fargo.par
```

#### Parallel GPU

```
make SETUP=fargo PARALLEL=1 GPU=1
mpirun -np 2 ./fargo3d in/fargo.par
```

------------------------

### Nested meshes

Nested meshes let you refine selected regions of the domain with finer grids, in the spirit of the JUPITER code. Each refinement level halves the cell size (in every dimension for which refinement is enabled), and levels are evolved with time subcycling. The feature works in serial, with MPI, on GPUs, and with CUDA-aware MPI (device-to-device halo exchanges between patches).

#### Defining the patches

Refined patches are declared directly in the parameter file, in a `#GRIDINFO` block placed at the end. Each line defines one patch with its extent and refinement level:

```
#GRIDINFO
#xmin ymin zmin   xmax ymax zmax   level
-0.5  0.8  1.52   0.5  1.2  1.5708  1
-0.25 0.9  1.55   0.25 1.1  1.5708  2
```

- Coordinates are absolute by default, in the order of the active dimensions (for a 2D run, a line has `xmin ymin xmax ymax level`). If the code is compiled with `-DRELNESTING` (see the `.opt` file), coordinates are instead given relative to the domain center and normalized by the half-width of each dimension (so the full domain spans -1 to 1).
- `level` is the refinement level: level 1 is twice as fine as the base mesh, level 2 four times, etc. A level-N patch must lie inside level-(N-1) refined territory.
- Patch edges are snapped outward to cell boundaries of the parent level. The grids actually used are reported in `nested_meshes.txt` in the output directory — check it after the first run.
- Refinement along individual dimensions can be disabled with the parameters `INHIBREFDIM1`, `INHIBREFDIM2`, `INHIBREFDIM3` (set to `Yes` in the `.par` file to keep that dimension unrefined on all levels).

#### Build options

See `setups/p3dthd/p3dthd.opt` for a working example. The relevant block is:

```
COMM = INTERNAL
FARGO_OPT += -DCOMMADAPT
FARGO_OPT += -DCOMMGPU      # GPU builds: inter-level comms filled on the device
#FARGO_OPT += -DRELNESTING  # optional: relative GRIDINFO coordinates
```

Then build and run as usual; nested meshes are activated by the presence of the `#GRIDINFO` block:

```
make SETUP=p3dthd PARALLEL=1 GPU=1 MPICUDA=1
mpirun -np 2 ./fargo3d setups/p3dthd/p3dthd.par
```

With MPI, every level is domain-decomposed across all ranks. Multi-rank results are identical to single-rank results (bitwise, in our tests), on both CPU and GPU.

#### Outputs and restarts

Field files are written per refinement level: `gasdens<n>_grid<level>.dat`, together with a `Descriptor<n>.dat` file that records the grid hierarchy and decomposition. Restarts read the grid hierarchy back from the descriptors, so the `#GRIDINFO` block is only parsed on fresh starts.

A small nested-mesh verification setup is provided in `setups/p3dthd/verify_cudacomm.par` (two refinement levels on a coarse 3D spherical disk), which can be used to validate a build by comparing single-rank against multi-rank outputs.

------------------------

### Description of subdirectories:

* planets: library of planetary systems.

* scripts: python scripts needed to build the code.

* setups: this is where all the custom setup definition are stored. The name of the setups correspond to the names of the directories found in setups/

* src: this is where all the source files of the code are found. Some of them may be redefined in the setups/ subdirectory (the makefile uses the VPATH variable, which behaves much like the PATH variable of the shell, as it allows to decide in which order a given source file is sought within different directories).

* std: this is where all the files that contain some standard definitions (everything that is not   a source file, not a script, and that the user is not supposed to modify). This includes, for   instance, the definition of the some standard boundary conditions, the units (scaling rules) of   the code parameters, etc.

* test_suite: contains python scripts that are used to test various features of the code. They are invoked as follows. We take the example of the permut.py script (which tests that the output of the Orszag-Tang vortex is independent of the choice of coordinates, XY, YZ or XZ). In the main directory (parent directory of test_suite/), simply issue: make testpermut The rule is therefore simply to issue: make test[name of python script without extension] for any script found in this subdirectory. All these scripts should use the 'test' python module found in scripts/

* utils: contains some utilities to post-process the data.