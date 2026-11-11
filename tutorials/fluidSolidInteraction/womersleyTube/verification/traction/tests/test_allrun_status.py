"""The tutorial Allrun exits with the status of a failing step (it used to
exit 0 after a failed solver, e.g. the serial-only diagnostic's fatal error in
a parallel run). The real solids4FoamScripts.sh is used; blockMesh, topoSet,
decomposePar, reconstructPar, solids4Foam and mpirun are stubs. Needs bash."""
import os
import shutil
import stat
import subprocess
import tempfile
import unittest
from pathlib import Path

TUTORIAL = Path(__file__).resolve().parents[3]
SCRIPTS = Path(__file__).resolve().parents[6] / "applications/scripts"

STUB = """#!/bin/bash
echo "$(basename "$0") $*" >> "$CASE_CALLS"
[ "$(basename "$0")" = "$FAIL_STEP" ] && { echo "stub failure" >&2; exit 7; }
exit 0
"""
MPIRUN = """#!/bin/bash
echo "mpirun $*" >> "$CASE_CALLS"
shift 2   # -n N
exec "$@"
"""


@unittest.skipUnless(shutil.which("bash") and SCRIPTS.is_dir(), "needs bash and the repository")
class AllrunStatus(unittest.TestCase):
    def run_allrun(self, args, fail_step=""):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            case, stubs, wm = tmp / "case", tmp / "bin", tmp / "wm"
            stubs.mkdir()
            (wm / "bin/tools").mkdir(parents=True)
            (wm / "bin/tools/RunFunctions").write_text("")
            for name in ("blockMesh", "topoSet", "decomposePar", "reconstructPar",
                         "solids4Foam"):
                (stubs / name).write_text(STUB)
            (stubs / "mpirun").write_text(MPIRUN)
            for f in stubs.iterdir():
                f.chmod(f.stat().st_mode | stat.S_IEXEC)
            shutil.copytree(TUTORIAL, case, ignore=shutil.ignore_patterns(
                "verification", "images", "work", "postProcessing"))
            for d in ("0/fluid",):
                (case / d).mkdir(parents=True, exist_ok=True)
            (case / "constant/fsiProperties.robin").write_text("")
            (case / "constant/fsiProperties.iqnils").write_text("")
            for f in ("U.robin", "p.robin", "U.dirichletNeumann", "p.dirichletNeumann",
                      "U.iqnils", "p.iqnils"):
                (case / "0/fluid" / f).write_text("")
            env = dict(os.environ, PATH=f"{stubs}:{SCRIPTS}:{os.environ['PATH']}",
                       WM_PROJECT_DIR=str(wm), WM_PROJECT="OpenFOAM",
                       WM_PROJECT_VERSION="v2512", PETSC_DIR="/stub",
                       CASE_CALLS=str(tmp / "calls"), FAIL_STEP=fail_step)
            (tmp / "calls").write_text("")
            r = subprocess.run(["bash", "./Allrun", *args], cwd=case, env=env,
                               capture_output=True, text=True)
            return r.returncode, (tmp / "calls").read_text()

    def test_success_is_zero(self):
        for args in (["robin"], ["robin", "parallel"]):
            rc, calls = self.run_allrun(args)
            self.assertEqual(rc, 0, args)
            self.assertIn("solids4Foam", calls)

    def test_failed_serial_solver(self):
        rc, calls = self.run_allrun(["robin"], "solids4Foam")
        self.assertEqual(rc, 7)

    def test_failed_parallel_solver_stops_before_reconstruct(self):
        rc, calls = self.run_allrun(["robin", "parallel"], "solids4Foam")
        self.assertEqual(rc, 7)
        self.assertIn("mpirun -n 2 solids4Foam -parallel", calls)
        self.assertNotIn("reconstructPar", calls)

    def test_failed_earlier_steps(self):
        for step in ("blockMesh", "topoSet"):
            rc, calls = self.run_allrun(["robin"], step)
            self.assertEqual(rc, 7, step)
            self.assertNotIn("solids4Foam", calls)
        rc, calls = self.run_allrun(["robin", "parallel"], "decomposePar")
        self.assertEqual(rc, 7)
        self.assertNotIn("mpirun", calls)


if __name__ == "__main__":
    unittest.main()
