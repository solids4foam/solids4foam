"""Run/skip logic of run_single.py: reuse keyed on stored settings (input
hashes, build signature, study settings), failure isolation, no deletion of
valid runs."""
import json
import shutil
import sys
import tempfile
import types
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import study_runs  # noqa: E402


class FakeDriver:
    SETTINGS_FILE = "verification_settings.json"

    def __init__(self, root, fail_names=()):
        self.WORK_ROOT = Path(root)
        self.sig = {"tutorialFiles": {"a": "h1"}, "build": {"solver": "1:2"}}
        self.fail_names = set(fail_names)
        self.runs = []
        self.analyse_calls = 0

    def settings(self, spec, reference):
        return dict(spec, **self.sig)

    def solver_log_problem(self, case):
        log = case / "log.solids4Foam"
        if not log.is_file():
            return "no log"
        return None if "End" in log.read_text() else "did not finish"

    def run_case(self, spec, reference, reuse):
        assert reuse is False
        self.runs.append(spec["name"])
        case = self.WORK_ROOT / spec["name"]
        if case.exists():
            shutil.rmtree(case)
        if spec["name"] in self.fail_names:
            raise RuntimeError("solver failed")
        diag = case / "postProcessing/interfaceTraction"
        diag.mkdir(parents=True)
        for n in study_runs.DIAGNOSTIC_FILES:
            (diag / n).write_text("1 2 3\n")
        (case / "log.solids4Foam").write_text("End\n")
        (case / self.SETTINGS_FILE).write_text(json.dumps(self.settings(spec, None)))
        return f"{spec['name']}: completed"


def spec(name):
    return {"coupling": "robin", "factor": 2, "steps": 200, "name": name}


class StudyRuns(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.V = FakeDriver(self.tmp.name)

    def tearDown(self):
        self.tmp.cleanup()

    def test_valid_run_is_reused_not_deleted(self):
        self.assertEqual(study_runs.run_all(self.V, None, [spec("a")]), [])
        marker = self.V.WORK_ROOT / "a/marker"
        marker.write_text("keep")
        self.assertEqual(study_runs.run_all(self.V, None, [spec("a")]), [])
        self.assertEqual(self.V.runs, ["a"])        # simulated once
        self.assertTrue(marker.exists())

    def test_no_analysis_is_made(self):
        self.assertFalse(hasattr(self.V, "analyse"))
        study_runs.run_all(self.V, None, [spec("a")])

    def test_failure_does_not_stop_later_runs(self):
        self.V.fail_names = {"a"}
        failed = study_runs.run_all(self.V, None, [spec("a"), spec("b")])
        self.assertEqual([n for n, _ in failed], ["a"])
        self.assertEqual(self.V.runs, ["a", "b"])
        self.assertIsNone(study_runs.reuse_problem(
            self.V, self.V.WORK_ROOT / "b", self.V.settings(spec("b"), None)))

    def test_changed_inputs_or_build_force_a_rerun(self):
        study_runs.run_all(self.V, None, [spec("a")])
        self.V.sig = {"tutorialFiles": {"a": "h2"}, "build": {"solver": "1:2"}}
        study_runs.run_all(self.V, None, [spec("a")])
        self.V.sig = {"tutorialFiles": {"a": "h2"}, "build": {"solver": "9:9"}}
        study_runs.run_all(self.V, None, [spec("a")])
        self.assertEqual(self.V.runs, ["a", "a", "a"])

    def test_incomplete_or_missing_diagnostic_is_not_reused(self):
        study_runs.run_all(self.V, None, [spec("a")])
        case = self.V.WORK_ROOT / "a"
        (case / "postProcessing/interfaceTraction/fluid.dat").write_text("")
        self.assertIn("empty", study_runs.reuse_problem(self.V, case, self.V.settings(spec("a"), None)))
        (case / "log.solids4Foam").write_text("running")
        self.assertEqual(study_runs.reuse_problem(self.V, case, self.V.settings(spec("a"), None)),
                         "did not finish")

    def test_damaged_settings_file(self):
        study_runs.run_all(self.V, None, [spec("a")])
        case = self.V.WORK_ROOT / "a"
        (case / self.V.SETTINGS_FILE).write_text("{")
        self.assertEqual(study_runs.reuse_problem(self.V, case, {}), "unreadable stored settings")


if __name__ == "__main__":
    unittest.main()
