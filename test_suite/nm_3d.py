#!/usr/bin/env python3
"""Nested meshes in 3D: 1 rank vs 8 ranks, bit-exact interiors.

The 2x2x2 decomposition of the 3D setup creates corner-type rank
adjacencies that do not exist in 2D, so this exercises the full set of
face/edge/corner ghost communications of the nested-mesh machinery."""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.realpath(__file__)),
                                "..", "scripts"))
import nm_test as T

PAR = "setups/nmtest3d/nmtest3d.par"

t = T.NestedMeshTest("nm_3d")
try:
    cpu = t.build("nmtest3d")
    ref = t.run(cpu, PAR, "np1", np=1)
    par = t.run(cpu, PAR, "np8", np=8)
    t.compare(ref, par, snap=2)
except T.TestFailure as exc:
    t.check("execution", False, str(exc))
t.finish()
