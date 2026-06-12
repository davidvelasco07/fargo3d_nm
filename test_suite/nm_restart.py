#!/usr/bin/env python3
"""Nested-mesh restarts must reproduce an uninterrupted run bit-exactly.

A reference 1-rank run advances to snapshot 4. A second 1-rank run stops at
snapshot 2 and is restarted (-s 2) to snapshot 4. A restart rebuilds the
fine-level time-interpolation state from the dump instead of carrying it
over, so agreement with the uninterrupted run is to machine precision
(absolute tolerance 1e-13), not bit-exact. The same restart on 4 ranks
(nested restarts read sequential-type descriptors, so the source output
comes from the 1-rank run) must be bit-identical to the serial restart:
restart is required to be deterministic and decomposition-invariant."""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.realpath(__file__)),
                                "..", "scripts"))
import nm_test as T

PAR = "setups/nmtest/nmtest.par"

t = T.NestedMeshTest("nm_restart")
try:
    cpu = t.build("nmtest")
    ref = t.run(cpu, PAR, "ref", np=1, overrides={"Ntot": 4})
    half = t.run(cpu, PAR, "half", np=1, overrides={"Ntot": 2})

    # Serial restart: continue the truncated run to snapshot 4.
    rst1 = t.clone_output("half", "restart1")
    t.run(cpu, PAR, "restart1", np=1, overrides={"Ntot": 4}, extra_args="-s 2")
    t.compare(ref, rst1, snap=4, tol_abs=1e-13,
              label="serial restart vs straight run (machine precision)")

    # Parallel restart from the sequential (1-rank) snapshot.
    rst4 = t.clone_output("half", "restart4")
    t.run(cpu, PAR, "restart4", np=4, overrides={"Ntot": 4}, extra_args="-s 2")
    t.compare(ref, rst4, snap=4, tol_abs=1e-13,
              label="4-rank restart vs straight run (machine precision)")
    t.compare(rst1, rst4, snap=4,
              label="4-rank restart vs serial restart (bit-exact)")
except T.TestFailure as exc:
    t.check("execution", False, str(exc))
t.finish()
