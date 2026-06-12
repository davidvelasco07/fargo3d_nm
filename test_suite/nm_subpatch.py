#!/usr/bin/env python3
"""Posterior refinement (-r): add a subpatch to an existing output, then
restart from it and keep evolving.

A 1-rank run advances to snapshot 2. The -r tool then inserts a new level-1
patch into a region of the base grid not covered by the original GRIDINFO
patches (it interpolates the parent data, rewrites the descriptor and
exits). Restarting from the modified output must work, and the descriptor
of the final snapshot must contain the extra grid."""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.realpath(__file__)),
                                "..", "scripts"))
import nm_test as T
import nm_compare

PAR = "setups/nmtest/nmtest.par"
# Refine the base grid at output 2 (grids are numbered finest-first, so the
# base grid is #2). The box y in [5,10], z in [-5,5] (third dim unused in 2D)
# shares its y=5 face with the original level-1 patch: CheckBC rejects
# isolated interior boxes on this setup, and the adjacency additionally
# exercises same-level ghost exchange between old and new patches.
SUBPATCH = '-r "2 2 5.0 -5.0 0.0 10.0 5.0 0.0"'

t = T.NestedMeshTest("nm_subpatch")
try:
    cpu = t.build("nmtest")
    t.run(cpu, PAR, "base", np=1, overrides={"Ntot": 2})

    work = t.clone_output("base", "refined")
    # The -r pass reads snapshot 2, adds the subpatch and exits cleanly.
    t.run(cpu, PAR, "refined", np=1, extra_args="-s 2 " + SUBPATCH)
    t.check("-r reports completion",
            "Subpatch creation completed" in t.runlog("refined"),
            "completion message missing from log")

    desc = nm_compare.parse_descriptor(
        os.path.join(work, "output00002", "Descriptor2.dat"))
    t.check("descriptor gains a grid after -r", len(desc["grids"]) == 4,
            "expected 4 grids, found %d" % len(desc["grids"]))

    # The refined output must be restartable and evolve further.
    t.run(cpu, PAR, "refined", np=1, overrides={"Ntot": 4}, extra_args="-s 2")
    desc4 = nm_compare.parse_descriptor(
        os.path.join(work, "output00004", "Descriptor4.dat"))
    t.check("subpatch survives evolution to snapshot 4",
            len(desc4["grids"]) == 4,
            "expected 4 grids at snap 4, found %d" % len(desc4["grids"]))
except T.TestFailure as exc:
    t.check("execution", False, str(exc))
t.finish()
