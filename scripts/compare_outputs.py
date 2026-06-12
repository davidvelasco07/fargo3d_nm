#!/usr/bin/env python3
"""
Compare two FARGO3D nested-mesh split outputs (run with -k, FULLDEBUG/WRITEGHOSTS).

For every per-patch field file (gasdensity / gasvelocity) it reports the maximum
relative difference between the reference run and the test run, broken down into
  - interior cells
  - edge ghost cells (one ghost dimension)
  - corner ghost cells (two ghost dimensions)
so the correctness of corner/edge halo exchange can be assessed directly.

Usage:
  compare_outputs.py REF_DIR TEST_DIR [--snap N] [--ngh 3]
"""
import os
import sys
import glob
import argparse
import numpy as np


def parse_descriptor(path):
    """Return {gid: (level, nx, ny, nz)} and Nghost from a Descriptor file."""
    with open(path) as f:
        lines = f.readlines()
    ndim, coordtype = (int(x) for x in lines[1].split()[:2])
    parts = lines[2].split()
    Ngrid, Ncpu, Nghost, Nvar = (int(parts[i]) for i in range(4))
    dims = {}
    l = 3 + Nvar
    for _ in range(Ngrid):
        head = lines[l].split()
        level = int(head[1])
        cpus = int(head[2])
        for cpu in range(cpus):
            cpuN = [int(x) for x in lines[l + 8 + cpu * 3].split()]
            nx, ny, nz, gid = cpuN[0], cpuN[1], cpuN[2], cpuN[4]
            dims[gid] = (level, nx, ny, nz)
        l += 8 + 3 * cpus
    return dims, Nghost, ndim


def rel(a, b):
    denom = np.maximum(np.abs(a), 1e-30)
    return np.abs(a - b) / denom


def region_masks(shape, ngh):
    """shape is (nz, ny, nx); ngh is (ngz, ngy, ngx). Returns interior, edge, corner masks."""
    nz, ny, nx = shape
    ngz, ngy, ngx = ngh
    # per-axis "is ghost" boolean
    def ax(n, g):
        m = np.zeros(n, bool)
        if g > 0:
            m[:g] = True
            m[-g:] = True
        return m
    gz = ax(nz, ngz)[:, None, None]
    gy = ax(ny, ngy)[None, :, None]
    gx = ax(nx, ngx)[None, None, :]
    nghost_dirs = gz.astype(int) + gy.astype(int) + gx.astype(int)
    interior = nghost_dirs == 0
    edge = nghost_dirs == 1
    corner = nghost_dirs >= 2
    return interior, edge, corner


def summarize(ref_dir, test_dir, snap, default_ngh):
    sub = "output%05d" % snap
    desc = os.path.join(ref_dir, sub, "Descriptor%d.dat" % snap)
    dims, Nghost, ndim = parse_descriptor(desc)
    if default_ngh is not None:
        Nghost = default_ngh
    if Nghost == 0:
        # WRITEGHOSTS builds write ngh=0 in the descriptor (the ghosts are
        # folded into the patch sizes). Without the true width the region
        # masks would be empty and every ghost diff vacuously zero.
        sys.exit("ERROR: descriptor reports Nghost=0; pass the real ghost "
                 "width explicitly with --ngh (3 for standard builds).")

    fields = {"gasdensity": 1, "gasvelocity": ndim}
    agg = {}  # field -> {region -> max rel}
    worst = {}
    for field, ncomp in fields.items():
        pattern = os.path.join(ref_dir, sub, "%s%d_*.dat" % (field, snap))
        files = sorted(glob.glob(pattern))
        a_int = a_edge = a_cor = 0.0
        a_int_abs = a_edge_abs = a_cor_abs = 0.0
        ncells_corner = 0
        nbad_corner = 0
        wf = None
        for rf in files:
            tf = os.path.join(test_dir, sub, os.path.basename(rf))
            if not os.path.exists(tf):
                print("  MISSING in test: %s" % os.path.basename(rf))
                continue
            base = os.path.basename(rf)
            gid = int(base.split("_")[1])
            ref = np.fromfile(rf)
            tst = np.fromfile(tf)
            if ref.size != tst.size or ref.size == 0:
                print("  SIZE MISMATCH %s (%d vs %d)" % (base, ref.size, tst.size))
                continue
            r = rel(ref, tst)
            absd = np.abs(ref - tst)
            level, nx, ny, nz = dims.get(gid, (None, None, None, None))
            shaped = (nx is not None) and (ncomp * nx * ny * nz == ref.size)
            if shaped:
                arr_r = r.reshape(ncomp, nz, ny, nx) if ncomp > 1 else r.reshape(nz, ny, nx)
                arr_a = absd.reshape(arr_r.shape)
                ngh = (Nghost if nz > 1 else 0,
                       Nghost if ny > 1 else 0,
                       Nghost if nx > 1 else 0)
                interior, edge, corner = region_masks((nz, ny, nx), ngh)
                if ncomp > 1:
                    interior = np.broadcast_to(interior, arr_r.shape)
                    edge = np.broadcast_to(edge, arr_r.shape)
                    corner = np.broadcast_to(corner, arr_r.shape)
                ci = arr_r[interior].max() if interior.any() else 0.0
                ce = arr_r[edge].max() if edge.any() else 0.0
                cc = arr_r[corner].max() if corner.any() else 0.0
                a_int = max(a_int, ci); a_edge = max(a_edge, ce); a_cor = max(a_cor, cc)
                a_int_abs = max(a_int_abs, arr_a[interior].max() if interior.any() else 0.0)
                a_edge_abs = max(a_edge_abs, arr_a[edge].max() if edge.any() else 0.0)
                a_cor_abs = max(a_cor_abs, arr_a[corner].max() if corner.any() else 0.0)
                if corner.any():
                    ncells_corner += int(corner.sum())
                    nbad_corner += int((arr_r[corner] > 1e-3).sum())
                if wf is None or cc >= worst.get(field, (0,))[0]:
                    worst[field] = (cc, base, "corner")
            else:
                # fall back to flat comparison
                mx = r.max()
                a_int = max(a_int, mx)
                if wf is None or mx >= worst.get(field, (0,))[0]:
                    worst[field] = (mx, base, "flat")
            wf = True
        agg[field] = dict(interior=a_int, edge=a_edge, corner=a_cor,
                          interior_abs=a_int_abs, edge_abs=a_edge_abs, corner_abs=a_cor_abs,
                          ncorner=ncells_corner, nbadcorner=nbad_corner)
    return agg, worst


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("ref_dir")
    ap.add_argument("test_dir")
    ap.add_argument("--snap", type=int, default=2)
    ap.add_argument("--ngh", type=int, default=None)
    args = ap.parse_args()

    print("=" * 70)
    print("REF : %s" % args.ref_dir)
    print("TEST: %s" % args.test_dir)
    print("snapshot #%d" % args.snap)
    print("=" * 70)
    agg, worst = summarize(args.ref_dir, args.test_dir, args.snap, args.ngh)
    overall_corner = 0.0
    for field, d in agg.items():
        print("\n[%s]" % field)
        print("  max rel diff  interior = %.3e   edge-ghost = %.3e   CORNER-ghost = %.3e"
              % (d["interior"], d["edge"], d["corner"]))
        print("  max abs diff  interior = %.3e   edge-ghost = %.3e   CORNER-ghost = %.3e"
              % (d["interior_abs"], d["edge_abs"], d["corner_abs"]))
        if d["ncorner"]:
            print("  corner ghost cells with rel diff > 1e-3 : %d / %d"
                  % (d["nbadcorner"], d["ncorner"]))
        if field in worst:
            print("  worst file: %s (%s, rel=%.3e)" % (worst[field][1], worst[field][2], worst[field][0]))
        overall_corner = max(overall_corner, d["corner"])

    print("\n" + "=" * 70)
    THRESH = 1e-3
    if overall_corner > THRESH:
        print("RESULT: CORNER GHOST MISMATCH  (max corner rel diff = %.3e > %.1e)"
              % (overall_corner, THRESH))
        sys.exit(1)
    else:
        print("RESULT: corner ghosts agree with reference (max corner rel diff = %.3e)"
              % overall_corner)
        sys.exit(0)


if __name__ == "__main__":
    main()
