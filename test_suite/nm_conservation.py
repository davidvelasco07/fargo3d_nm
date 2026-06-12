#!/usr/bin/env python3
"""Mass budget of the base level in a (nearly) closed nested-mesh domain.

With PeriodicZ enabled and reflecting Y boundaries, the only mass-budget
changes come from the small numerical wall flux that this scheme exhibits
under gravity even WITHOUT refinement (measured at ~6e-6 relative over 4
steps with the body centered; placing patches changes it only marginally).
Exact conservation is therefore not an invariant we can demand; instead we
check that
  1. the relative drift stays below a calibrated bound (1e-4), an order of
     magnitude above the unrefined baseline, so a gross spurious source at
     the coarse-fine interfaces would fail the test, and
  2. the drift is bit-identical between 1-rank and 4-rank runs, which is
     the sharp decomposition-invariance signal for the nested comms.
The grid is uniform Cartesian, so total mass is proportional to the plain
sum of the base-level density canvas (cells under refined patches hold the
restriction of the fine data)."""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.realpath(__file__)),
                                "..", "scripts"))
import numpy as np

import nm_test as T
import nm_compare

PAR = "setups/nmtest/nmtest.par"
# Ybody=5 keeps the gravitating body away from the Ymin wall, where it
# would drive a much larger (purely physical-setup) wall flux.
OVR = {"PeriodicZ": "yes", "Ntot": 4, "Ybody": 5.0}
TOL = 1e-4  # relative drift over 4 steps; unrefined baseline is ~6e-6


def base_mass(outdir, snap):
    desc = nm_compare.parse_descriptor(
        os.path.join(outdir, "output%05d" % snap, "Descriptor%d.dat" % snap))
    canv = nm_compare.assemble(outdir, snap, desc, "gasdensity", 1)
    base = min(desc["grids"], key=lambda g: g["level"])
    return float(np.sum(canv[base["number"]]))


t = T.NestedMeshTest("nm_conservation")
try:
    cpu = t.build("nmtest")
    final = {}
    for tag, np_ in (("np1", 1), ("np4", 4)):
        out = t.run(cpu, PAR, tag, np=np_, overrides=OVR)
        m0 = base_mass(out, 0)
        m4 = base_mass(out, 4)
        final[tag] = m4
        drift = abs(m4 - m0) / abs(m0)
        t.check("mass drift %s = %.3e < %g" % (tag, drift, TOL), drift < TOL,
                "m0=%.15e m4=%.15e" % (m0, m4))
    t.check("drift identical across decompositions",
            final["np1"] == final["np4"],
            "np1=%.17e np4=%.17e" % (final["np1"], final["np4"]))
except T.TestFailure as exc:
    t.check("execution", False, str(exc))
t.finish()
