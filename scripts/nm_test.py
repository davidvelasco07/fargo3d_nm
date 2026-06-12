#!/usr/bin/env python3
"""
Test framework for the nested-mesh (JUPITER-style) machinery.

Unlike the legacy GenericTest framework (scripts/test.py, Python 2, merged
single-grid outputs), this module works with the per-patch outputs of nested
runs and compares the physical (interior) cells of every refinement level,
reassembled onto full level canvases via scripts/nm_compare.py. Ghost cells
are never written or compared: standard builds only (no WRITEGHOSTS).

A test file in test_suite/ typically does:

    import nm_test as T
    t = T.NestedMeshTest("my_test")
    cpu = t.build("nmtest", PARALLEL=1)
    ref = t.run(cpu, "setups/nmtest/nmtest.par", "np1", np=1)
    par = t.run(cpu, "setups/nmtest/nmtest.par", "np4", np=4)
    t.compare(ref, par, snap=2)          # bit-exact by default
    t.finish()

Environment variables:
    NM_MPIRUN        override the mpirun launcher (default: "mpirun")
    NM_FARGO_ARCH    FARGO_ARCH passed to the build (default: env FARGO_ARCH)
    NM_KEEP          if set, keep test output directories on success
    NM_FORCE_REBUILD if set, rebuild binaries even if cached
All artifacts live under test_out/ (gitignored).
"""

import os
import shlex
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.realpath(__file__))
ROOT = os.path.realpath(os.path.join(HERE, ".."))
sys.path.insert(0, HERE)

import nm_compare  # noqa: E402

OUTROOT = os.path.join(ROOT, "test_out")
BINDIR = os.path.join(OUTROOT, "bin")

GREEN = "\033[1;32m"
RED = "\033[1;31m"
RESET = "\033[0m"

# The build system remembers flags between invocations (sticky flags in
# std/.lastflags), so every flag that could leak from a previous build must
# be pinned explicitly.
BUILD_DEFAULTS = {
    "PARALLEL": 1,
    "GPU": 0,
    "MPICUDA": 0,
    "FULLDEBUG": 0,
    "FARGO_DISPLAY": "NONE",
}


def _flaglabel(setup, flags):
    items = sorted(flags.items())
    return setup + "_" + "_".join("%s%s" % (k, v) for k, v in items)


class TestFailure(Exception):
    pass


class NestedMeshTest:
    def __init__(self, name):
        self.name = name
        self.workdir = os.path.join(OUTROOT, name)
        self.failures = []
        self.steps = []
        os.makedirs(self.workdir, exist_ok=True)
        os.makedirs(BINDIR, exist_ok=True)
        print("=" * 70)
        print("NESTED-MESH TEST: %s" % name)
        print("=" * 70)

    # ----------------------------------------------------------- build
    def build(self, setup, **flags):
        """Build SETUP with the given make flags; cache the binary."""
        merged = dict(BUILD_DEFAULTS)
        merged.update(flags)
        label = _flaglabel(setup, merged)
        cached = os.path.join(BINDIR, "fargo3d_" + label)
        if os.path.exists(cached) and not os.environ.get("NM_FORCE_REBUILD"):
            print("[build] using cached binary %s" % os.path.relpath(cached, ROOT))
            return cached
        env = dict(os.environ)
        arch = os.environ.get("NM_FARGO_ARCH", os.environ.get("FARGO_ARCH", ""))
        if arch:
            env["FARGO_ARCH"] = arch
        cmd = ["make", "SETUP=%s" % setup] + \
              ["%s=%s" % (k, v) for k, v in sorted(merged.items())]
        print("[build] %s" % " ".join(cmd))
        log = os.path.join(self.workdir, "build_%s.log" % label)
        binary = os.path.join(ROOT, "fargo3d")
        # make.py can return success even when compilation failed, so remove
        # any stale binary first: its absence then reveals the failure.
        if os.path.exists(binary):
            os.remove(binary)
        with open(log, "w") as fh:
            rc = subprocess.call(cmd, cwd=ROOT, stdout=fh,
                                 stderr=subprocess.STDOUT, env=env)
        if rc != 0 or not os.path.exists(binary):
            with open(log) as fh:
                tail = "".join(fh.readlines()[-25:])
            raise TestFailure("build failed (rc=%d), tail of %s:\n%s"
                              % (rc, os.path.relpath(log, ROOT), tail))
        shutil.copy2(binary, cached)
        return cached

    # ------------------------------------------------------------- run
    def run(self, binary, parfile, outname, np=1, overrides=None,
            extra_args="", expect_fail=False, env=None):
        """Run a binary on parfile; output goes to <workdir>/<outname>.

        'overrides' is a dict of parameter overrides applied with -o.
        OutputDir is always overridden. Returns the output directory.
        For restarts (extra_args containing -s/-S) the existing directory
        is kept.
        """
        outdir = os.path.join(self.workdir, outname)
        restarting = any(a in ("-s", "-S") for a in shlex.split(extra_args or ""))
        if not restarting:
            shutil.rmtree(outdir, ignore_errors=True)
        os.makedirs(outdir, exist_ok=True)
        launcher = os.environ.get("NM_MPIRUN", "mpirun")
        cmd = shlex.split(launcher) + ["-np", str(np), "--bind-to", "none"]
        if hasattr(os, "geteuid") and os.geteuid() == 0:
            cmd.append("--allow-run-as-root")
        cmd += [binary]
        if extra_args:
            cmd += shlex.split(extra_args)
        for k, v in (overrides or {}).items():
            cmd += ["-o", "%s=%s" % (k, v)]
        cmd += ["-o", "OutputDir=%s" % outdir, parfile]
        print("[run  ] (np=%d) -> %s" % (np, os.path.relpath(outdir, ROOT)))
        log = os.path.join(self.workdir, "run_%s.log" % outname)
        runenv = dict(os.environ)
        runenv.update(env or {})
        with open(log, "a" if restarting else "w") as fh:
            rc = subprocess.call(cmd, cwd=ROOT, stdout=fh,
                                 stderr=subprocess.STDOUT, env=runenv)
        if expect_fail:
            if rc == 0:
                raise TestFailure("run %s succeeded but was expected to fail"
                                  % outname)
            return outdir
        if rc != 0:
            with open(log) as fh:
                tail = "".join(fh.readlines()[-25:])
            raise TestFailure("run %s failed (rc=%d), tail of %s:\n%s"
                              % (outname, rc, os.path.relpath(log, ROOT), tail))
        return outdir

    def runlog(self, outname):
        """Contents of the log of a previous run."""
        with open(os.path.join(self.workdir, "run_%s.log" % outname)) as fh:
            return fh.read()

    def clone_output(self, src_outname, dst_outname):
        """Copy a run's output directory (e.g. to restart inside the copy)."""
        src = os.path.join(self.workdir, src_outname)
        dst = os.path.join(self.workdir, dst_outname)
        shutil.rmtree(dst, ignore_errors=True)
        shutil.copytree(src, dst)
        return dst

    # --------------------------------------------------------- checks
    def compare(self, ref_dir, test_dir, snap, tol_rel=0.0, tol_abs=0.0,
                label=None):
        """Interior-cell comparison of two runs, reassembled per level."""
        label = label or "%s vs %s (snap %d)" % (os.path.basename(ref_dir),
                                                 os.path.basename(test_dir),
                                                 snap)
        try:
            ok = nm_compare.compare(ref_dir, test_dir, snap, tol_rel, tol_abs)
        except nm_compare.DescriptorError as exc:
            self._record("compare: " + label, False, "structural: %s" % exc)
            return False
        self._record("compare: " + label, ok,
                     "interior diffs exceed tol_rel=%g tol_abs=%g"
                     % (tol_rel, tol_abs))
        return ok

    def check(self, label, condition, detail=""):
        """Record an arbitrary boolean check."""
        self._record(label, bool(condition), detail)
        return bool(condition)

    def _record(self, label, ok, detail):
        self.steps.append((label, ok))
        if ok:
            print("  %sPASS%s  %s" % (GREEN, RESET, label))
        else:
            print("  %sFAIL%s  %s  (%s)" % (RED, RESET, label, detail))
            self.failures.append(label)

    # --------------------------------------------------------- finish
    def finish(self):
        print("-" * 70)
        if self.failures:
            print("%sFAILED%s  %s : %d/%d checks failed"
                  % (RED, RESET, self.name, len(self.failures),
                     len(self.steps)))
            for f in self.failures:
                print("   - %s" % f)
            sys.exit(1)
        print("%sPASSED%s  %s : %d checks"
              % (GREEN, RESET, self.name, len(self.steps)))
        if not os.environ.get("NM_KEEP"):
            for entry in os.listdir(self.workdir):
                p = os.path.join(self.workdir, entry)
                if os.path.isdir(p):
                    shutil.rmtree(p, ignore_errors=True)
        sys.exit(0)
