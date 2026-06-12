#!/usr/bin/env python3
"""
Compare two nested-mesh FARGO3D runs on their physical (interior) cells.

Outputs contain one file per CPU patch, so two runs with different MPI
decompositions cannot be compared file by file. This tool reassembles, for
every refinement-level grid, the full interior data from the per-CPU patch
files (using the placement information in DescriptorN.dat) and compares the
resulting canvases between the two runs. Ghost cells are never involved:
outputs are expected to come from standard builds (no WRITEGHOSTS).

It is strict: missing files, size mismatches, or incomplete coverage of a
level canvas are errors, not silent skips.

Usage:
  nm_compare.py REF_DIR TEST_DIR [--snap N] [--tol-rel X] [--tol-abs Y]

Exit code: 0 if all fields on all levels agree within tolerance, 1 otherwise
(2 on structural errors).
"""
import argparse
import os
import sys

import numpy as np

FIELDS = {"gasdensity": 1, "gasvelocity": None}  # None -> NDIM components


class DescriptorError(Exception):
    pass


def parse_descriptor(path):
    """Parse DescriptorN.dat -> dict with ndim and grid/patch layout."""
    if not os.path.exists(path):
        raise DescriptorError("missing descriptor %s" % path)
    with open(path) as fh:
        lines = fh.readlines()
    ndim = int(lines[1].split()[0])
    ngrid, ncpu, ngh, nvar = (int(x) for x in lines[2].split()[:4])
    grids = []
    l = 3 + nvar
    for _ in range(ngrid):
        number, level, cpus = (int(x) for x in lines[l].split()[:3])
        sizes = tuple(int(x) for x in lines[l + 1].split()[:3])
        patches = []
        for c in range(cpus):
            base = l + 8 + 3 * c
            psize = [int(x) for x in lines[base].split()]
            ppos = [int(x) for x in lines[base + 1].split()]
            patches.append(dict(size=tuple(psize[0:3]), cpu=psize[3],
                                gid=psize[4],
                                pcmin=tuple(ppos[0:3]), pcmax=tuple(ppos[3:6])))
        grids.append(dict(number=number, level=level, sizes=sizes,
                          patches=patches))
        l += 8 + 3 * cpus
    return dict(ndim=ndim, ngh=ngh, grids=grids)


def assemble(outdir, snap, desc, field, ncomp):
    """Reassemble one field of every level grid onto its full canvas.

    Returns {grid_number: ndarray of shape (ncomp, n2, n1, n0)}.
    """
    sub = os.path.join(outdir, "output%05d" % snap)
    canvases = {}
    for g in desc["grids"]:
        n0, n1, n2 = g["sizes"]
        canvas = np.full((ncomp, n2, n1, n0), np.nan)
        for p in g["patches"]:
            s0, s1, s2 = p["size"]
            fn = os.path.join(sub, "%s%d_%d_%d.dat"
                              % (field, snap, p["gid"], g["level"]))
            if not os.path.exists(fn):
                raise DescriptorError("missing data file %s" % fn)
            data = np.fromfile(fn)
            if data.size != ncomp * s0 * s1 * s2:
                raise DescriptorError(
                    "size mismatch %s: %d doubles, expected %d (=%dx%dx%dx%d)"
                    % (fn, data.size, ncomp * s0 * s1 * s2, ncomp, s2, s1, s0))
            arr = data.reshape(ncomp, s2, s1, s0)
            i0, i1, i2 = p["pcmin"]
            canvas[:, i2:i2 + s2, i1:i1 + s1, i0:i0 + s0] = arr
        if np.isnan(canvas).any():
            raise DescriptorError(
                "incomplete coverage of grid %d (level %d) for %s in %s"
                % (g["number"], g["level"], field, sub))
        canvases[g["number"]] = canvas
    return canvases


def compare(ref_dir, test_dir, snap, tol_rel, tol_abs, verbose=True):
    """Return True if runs agree within tolerance on all interior cells."""
    dref = parse_descriptor(os.path.join(ref_dir, "output%05d" % snap,
                                         "Descriptor%d.dat" % snap))
    dtst = parse_descriptor(os.path.join(test_dir, "output%05d" % snap,
                                         "Descriptor%d.dat" % snap))
    if len(dref["grids"]) != len(dtst["grids"]):
        raise DescriptorError("different number of level grids (%d vs %d)"
                              % (len(dref["grids"]), len(dtst["grids"])))
    ok = True
    for field, nc in FIELDS.items():
        ncomp = nc if nc else dref["ndim"]
        cref = assemble(ref_dir, snap, dref, field, ncomp)
        ctst = assemble(test_dir, snap, dtst, field, ncomp)
        for g in dref["grids"]:
            num, lev = g["number"], g["level"]
            if num not in ctst:
                raise DescriptorError("grid %d missing in test run" % num)
            a, b = cref[num], ctst[num]
            if a.shape != b.shape:
                raise DescriptorError("grid %d shape mismatch %s vs %s"
                                      % (num, a.shape, b.shape))
            absd = np.abs(a - b)
            maxabs = absd.max()
            denom = np.maximum(np.abs(a), 1e-30)
            maxrel = (absd / denom).max()
            viol = int((absd > tol_abs + tol_rel * np.abs(a)).sum())
            good = viol == 0
            ok = ok and good
            if verbose:
                print("  %-12s level %d (grid %d): max abs = %.3e  "
                      "max rel = %.3e  %s"
                      % (field, lev, num, maxabs, maxrel,
                         "OK" if good else "VIOLATIONS=%d" % viol))
    return ok


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("ref_dir")
    ap.add_argument("test_dir")
    ap.add_argument("--snap", type=int, default=2)
    ap.add_argument("--tol-rel", type=float, default=0.0)
    ap.add_argument("--tol-abs", type=float, default=0.0)
    args = ap.parse_args()
    print("=" * 70)
    print("REF : %s\nTEST: %s\nsnapshot #%d   tol_rel=%g tol_abs=%g"
          % (args.ref_dir, args.test_dir, args.snap, args.tol_rel,
             args.tol_abs))
    print("=" * 70)
    try:
        ok = compare(args.ref_dir, args.test_dir, args.snap,
                     args.tol_rel, args.tol_abs)
    except DescriptorError as exc:
        print("STRUCTURAL ERROR: %s" % exc)
        sys.exit(2)
    print("RESULT: %s" % ("runs agree within tolerance" if ok
                          else "RUNS DIFFER beyond tolerance"))
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
