#!/usr/bin/env python3
"""Nested meshes: a 1-rank and a 4-rank run of the same 2D problem must
produce bit-identical interior cells on every refinement level.

The 4-rank decomposition makes the nested patches straddle rank boundaries,
exercising the inter-rank ghost communications (faces, edges and same-level
neighbour exchanges); any error there contaminates interior cells within a
couple of timesteps and is caught by the bit-exact comparison."""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.realpath(__file__)),
                                "..", "scripts"))
import nm_test as T

t = T.NestedMeshTest("nm_serial_vs_mpi")
try:
    cpu = t.build("nmtest")
    ref = t.run(cpu, "setups/nmtest/nmtest.par", "np1", np=1)
    par = t.run(cpu, "setups/nmtest/nmtest.par", "np4", np=4)
    t.compare(ref, par, snap=2)
except T.TestFailure as exc:
    t.check("execution", False, str(exc))
t.finish()
