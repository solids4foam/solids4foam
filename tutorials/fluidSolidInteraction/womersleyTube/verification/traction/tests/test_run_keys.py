"""Run keys with the REAL verification driver (a copy of the tutorial in a
temporary directory, with a stub Allrun instead of the solver):

- the solver-run key changes with tutorial files, settings, reference
  parameters/schedule and the solver build (by content), and does not change
  with the analysis code, the reference's exact values, a README or a touched
  binary;
- the analysis key changes with the driver and the reference JSON;
- reuse and check_runs validate settings, log, diagnostic output and build.
"""
import importlib.util
import json
import os
import shutil
import stat
import sys
import tempfile
import unittest
from pathlib import Path

sys.dont_write_bytecode = True      # edited copies must not load a stale .pyc
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "scripts"))
import check_runs  # noqa: E402
import study_runs  # noqa: E402

TUTORIAL = HERE.parents[2]          # .../womersleyTube
SPEC = {"coupling": "robin", "factor": 1, "steps": 50, "name": "robin_m1_n50"}

STUB_ALLRUN = """#!/bin/bash
echo run >> "{counter}"
mkdir -p postProcessing/interfaceTraction
echo "0 0 0" > postProcessing/interfaceTraction/fluid.dat
echo "0 0 0" > postProcessing/interfaceTraction/solid.dat
echo "{log}" > log.solids4Foam
"""


def load_driver(root: Path, tag: str):
    path = root / "verification/scripts/womersley_tube_verification.py"
    spec = importlib.util.spec_from_file_location(f"wtv_test_{tag}", path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def make_study(dest: Path, counter: Path, log="End"):
    dest.mkdir(parents=True)
    for d in ("0", "constant", "system"):
        shutil.copytree(TUTORIAL / d, dest / d, symlinks=True)
    shutil.copytree(TUTORIAL / "verification/scripts", dest / "verification/scripts",
                    ignore=shutil.ignore_patterns("__pycache__"))
    shutil.copytree(TUTORIAL / "verification/reference", dest / "verification/reference")
    (dest / "README.md").write_text("readme")
    allrun = dest / "Allrun"
    allrun.write_text(STUB_ALLRUN.format(counter=counter, log=log))
    allrun.chmod(allrun.stat().st_mode | stat.S_IXUSR)


class Keys(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)
        self.counter = self.root / "counter"
        self.bin = self.root / "bin"
        self.lib = self.root / "lib"
        self.bin.mkdir()
        self.lib.mkdir()
        (self.bin / "solids4Foam").write_text("binary-v1")
        (self.bin / "solids4Foam").chmod(0o755)
        (self.lib / "libsolids4FoamModels.so").write_text("library-v1")
        self.env = dict(os.environ)
        os.environ["PATH"] = f"{self.bin}{os.pathsep}{os.environ['PATH']}"
        os.environ["FOAM_USER_LIBBIN"] = str(self.lib)
        # a study copy as <work>/<name>/womersleyTube
        self.study = self.root / "work/mu5/womersleyTube"
        make_study(self.study, self.counter)
        self.n = 0
        self.V = self.driver()
        self.ref = json.loads(self.V.REFERENCE_FILE.read_text())

    def tearDown(self):
        os.environ.clear()
        os.environ.update(self.env)
        self.tmp.cleanup()

    def driver(self):
        self.n += 1       # a fresh module: no cached hashes or signatures
        return load_driver(self.study, str(self.n))

    def solver_key(self, V=None):
        return study_runs.solver_run_key(V or self.driver(), SPEC, self.ref)

    def analysis_key(self):
        V = self.driver()
        files = {f"s/{f.name}": f for f in (self.study / "verification/scripts").glob("*.py")}
        files["ref"] = V.REFERENCE_FILE
        return study_runs.analysis_key(files)

    # ------------------------------------------------------------ keys
    def test_analysis_and_documentation_changes_do_not_change_the_solver_key(self):
        key = self.solver_key()
        analysis = self.analysis_key()
        driver = self.study / "verification/scripts/womersley_tube_verification.py"
        driver.write_text(driver.read_text() + "\n# analysis-side edit\n")
        (self.study / "README.md").write_text("changed")
        ref_file = self.study / "verification/reference/womersleyTube_verification_references.json"
        ref = json.loads(ref_file.read_text())
        ref["reference"]["exact"]["k"] = [9.0, 9.0]       # used only by the analysis
        ref_file.write_text(json.dumps(ref))
        self.assertEqual(self.solver_key(), key)
        self.assertNotEqual(self.analysis_key(), analysis)

    def test_driver_change_changes_analysis_key_only(self):
        key, analysis = self.solver_key(), self.analysis_key()
        driver = self.study / "verification/scripts/womersley_tube_verification.py"
        driver.write_text(driver.read_text() + "\n# edit\n")
        self.assertEqual(self.solver_key(), key)
        self.assertNotEqual(self.analysis_key(), analysis)

    def test_reference_json_change_changes_analysis_key(self):
        analysis = self.analysis_key()
        ref_file = self.study / "verification/reference/womersleyTube_verification_references.json"
        ref_file.write_text(ref_file.read_text() + "\n")
        self.assertNotEqual(self.analysis_key(), analysis)

    def test_tutorial_file_change_changes_the_solver_key(self):
        key = self.solver_key()
        f = self.study / "constant/womersleyProperties"
        f.write_text(f.read_text() + "\n// changed\n")
        self.assertNotEqual(self.solver_key(), key)

    def test_physical_or_schedule_setting_changes_the_solver_key(self):
        key = self.solver_key()
        ref_file = self.study / "verification/reference/womersleyTube_verification_references.json"
        ref = json.loads(ref_file.read_text())
        ref["reference"]["parameters"]["dynamicViscosity"] *= 2
        ref_file.write_text(json.dumps(ref))
        self.ref = ref
        k2 = self.solver_key()
        self.assertNotEqual(k2, key)
        ref["study"]["periods"] += 1
        ref_file.write_text(json.dumps(ref))
        self.assertNotEqual(self.solver_key(), k2)

    def test_driver_substitution_change_changes_the_solver_key(self):
        # what configure() writes into the case is part of the key, so a
        # changed substitution (here: the sampling interval) cannot reuse a run
        key = self.solver_key()
        driver = self.study / "verification/scripts/womersley_tube_verification.py"
        text = driver.read_text()
        self.assertIn('replace_entry(control, "timePrecision", "12")', text)
        driver.write_text(text.replace('replace_entry(control, "timePrecision", "12")',
                                       'replace_entry(control, "timePrecision", "10")'))
        self.assertNotEqual(self.solver_key(), key)

    def test_build_is_keyed_on_content_not_mtime(self):
        key = self.solver_key()
        exe = self.bin / "solids4Foam"
        status = exe.stat()
        os.utime(exe, ns=(status.st_atime_ns, status.st_mtime_ns + 10**9))   # touch
        self.assertEqual(self.solver_key(), key)
        shutil.copy2(exe, self.root / "moved")                                # same bytes
        exe.write_text("binary-v2")                                           # same size
        os.utime(exe, ns=(status.st_atime_ns, status.st_mtime_ns))            # same mtime
        self.assertNotEqual(self.solver_key(), key)
        exe.write_text("binary-v1")
        self.assertEqual(self.solver_key(), key)
        (self.lib / "libsolids4FoamModels.so").write_text("library-v2")
        self.assertNotEqual(self.solver_key(), key)

    # ------------------------------------------------------------ reuse
    def run_all(self, V=None):
        V = V or self.driver()
        return V, study_runs.run_all(V, self.ref, [dict(SPEC)])

    def runs(self):
        return self.counter.read_text().count("run") if self.counter.exists() else 0

    def test_real_driver_reuses_a_valid_run_and_reruns_after_changes(self):
        V, failed = self.run_all()
        self.assertEqual(failed, [])
        self.assertEqual(self.runs(), 1)
        self.run_all()
        self.assertEqual(self.runs(), 1)                  # reused
        (self.study / "README.md").write_text("changed")
        self.run_all()
        self.assertEqual(self.runs(), 1)                  # still reused
        (self.bin / "solids4Foam").write_text("binary-v2")
        self.run_all()
        self.assertEqual(self.runs(), 2)                  # rebuilt solver
        f = self.study / "constant/womersleyProperties"
        f.write_text(f.read_text() + "\n// changed\n")
        self.run_all()
        self.assertEqual(self.runs(), 3)                  # changed tutorial file

    def check(self, *flags, case="robin_m1_n50"):
        return check_runs.main([str(self.root / "work"), *flags, f"mu5:{case}"])

    def test_check_runs_validates_settings_log_diagnostic_and_build(self):
        self.run_all()
        self.assertEqual(self.check(), 0)
        case = self.study / "verification/work/robin_m1_n50"
        keys = self.root / "keys.json"
        self.assertEqual(self.check("--keys-out", str(keys)), 0)
        record = json.loads(keys.read_text())
        self.assertIn("mu5/robin_m1_n50", record["solverRunKeys"])
        self.assertTrue(any("womersley_tube_verification.py" in k
                            for k in record["analysisKey"]))
        # killed run: no End line although fluid.dat is non-empty
        (case / "log.solids4Foam").write_text("Time = 1\n")
        self.assertEqual(self.check(), 1)
        (case / "log.solids4Foam").write_text("End\n")
        self.assertEqual(self.check(), 0)
        # empty diagnostic output
        (case / "postProcessing/interfaceTraction/fluid.dat").write_text("")
        self.assertEqual(self.check(), 1)
        (case / "postProcessing/interfaceTraction/fluid.dat").write_text("0 0 0\n")
        # stored settings of a different run (tutorial changed afterwards)
        f = self.study / "constant/womersleyProperties"
        f.write_text(f.read_text() + "\n// changed\n")
        self.assertEqual(self.check(), 1)
        f.write_text(f.read_text().replace("\n// changed\n", ""))
        self.assertEqual(self.check(), 0)
        # a missing run
        self.assertEqual(self.check(case="robin_m2_n50"), 1)

    def test_recorded_analysis_key_follows_driver_and_reference_not_solver_keys(self):
        self.run_all()
        first, second = self.root / "k1.json", self.root / "k2.json"
        self.assertEqual(self.check("--keys-out", str(first)), 0)
        ref_file = self.study / "verification/reference/womersleyTube_verification_references.json"
        ref = json.loads(ref_file.read_text())
        ref["reference"]["exact"]["k"] = [9.0, 9.0]       # analysis-only value
        ref_file.write_text(json.dumps(ref))
        driver = self.study / "verification/scripts/womersley_tube_verification.py"
        driver.write_text(driver.read_text() + "\n# edit\n")
        self.assertEqual(self.check("--keys-out", str(second)), 0)    # runs still valid
        a, b = json.loads(first.read_text()), json.loads(second.read_text())
        self.assertEqual(a["solverRunKeys"], b["solverRunKeys"])
        changed = {k for k in a["analysisKey"] if a["analysisKey"][k] != b["analysisKey"][k]}
        self.assertIn("mu5/verification/scripts/womersley_tube_verification.py", changed)
        self.assertIn("mu5/reference/womersleyTube_verification_references.json", changed)

    def test_build_mismatch_fails_unless_skipped(self):
        self.run_all()
        (self.lib / "libsolids4FoamModels.so").write_text("library-v2")
        self.assertEqual(self.check(), 1)
        self.assertEqual(self.check("--skip-build-check"), 0)
        # but never a different tutorial, even when skipping the build check
        f = self.study / "constant/womersleyProperties"
        f.write_text(f.read_text() + "\n// changed\n")
        self.assertEqual(self.check("--skip-build-check"), 1)


if __name__ == "__main__":
    unittest.main()
