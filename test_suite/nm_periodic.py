#!/usr/bin/env python3
"""Nested meshes with a periodic direction: 1 rank vs 4 ranks, bit-exact.

With PeriodicZ enabled, ranks exchange ghost zones across the wrap-around
boundary as well, which multiplies the number of same-(peer,facedim,level)
communicator pairs. This is exactly the configuration where non-unique MPI
tags (the facedim-as-tag bug fixed in J_comm_device.c) would be most
dangerous, so it deserves its own decomposition-invariance check."""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.realpath(__file__)),
                                "..", "scripts"))
import nm_test as T

OVR = {"PeriodicZ": "yes"}

t = T.NestedMeshTest("nm_periodic")
try:
    cpu = t.build("nmtest")
    ref = t.run(cpu, "setups/nmtest/nmtest.par", "np1", np=1, overrides=OVR)
    par = t.run(cpu, "setups/nmtest/nmtest.par", "np4", np=4, overrides=OVR)
    t.compare(ref, par, snap=2)
except T.TestFailure as exc:
    t.check("execution", False, str(exc))
t.finish()
