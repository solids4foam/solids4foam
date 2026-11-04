"""prepare_study.py: counted tolerance substitutions, atomic/marked study
copies, reuse only of matching copies, no silent deletion of runs, the
diagnostic taken from the scripts directory."""
import contextlib
import io
import json
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import prepare_study as P  # noqa: E402

ROBIN = """\
    outerCorrTolerance   1e-6;
    robinPressureTolerance 1e-5;
    robinFluxTolerance     5e-3;
"""
KEYS = P.TOLERANCE_KEYS["fsiProperties.robin"]


class ScaleTolerances(unittest.TestCase):
    def test_scales_each_key_once(self):
        out = P.scale_tolerances(ROBIN, KEYS, 0.01)
        self.assertIn("outerCorrTolerance   1e-08;", out)
        self.assertIn("robinPressureTolerance 1e-07;", out)
        self.assertIn("robinFluxTolerance     5e-05;", out)

    def test_space_before_semicolon_and_exponent_forms(self):
        out = P.scale_tolerances("outerCorrTolerance 1.0e-06 ;\n", ("outerCorrTolerance",), 0.5)
        self.assertIn("5e-07 ;", out)

    def test_missing_key_is_an_error(self):
        with self.assertRaises(SystemExit):
            P.scale_tolerances(ROBIN.replace("robinFluxTolerance", "otherTolerance"), KEYS, 0.01)

    def test_key_moved_into_subdict_or_commented_is_not_matched(self):
        with self.assertRaises(SystemExit):
            P.scale_tolerances("// outerCorrTolerance 1e-6;\n", ("outerCorrTolerance",), 0.01)

    def test_duplicate_key_is_an_error(self):
        with self.assertRaises(SystemExit):
            P.scale_tolerances(ROBIN + ROBIN, KEYS, 0.01)


def make_tutorial(root: Path):
    (root / "system").mkdir(parents=True)
    (root / "constant").mkdir()
    (root / "verification/reference").mkdir(parents=True)
    (root / "system/controlDict").write_text("a 1;\nfunctions\n{\n}\n")
    (root / "constant/fsiProperties.robin").write_text(ROBIN)
    (root / "constant/fsiProperties.iqnils").write_text("outerCorrTolerance 1e-6;\n")
    (root / "verification/reference/womersleyTube_verification_references.json").write_text(
        json.dumps({"reference": {"derived": {"period": 1.0}},
                    "study": {"periods": 6, "analysisPeriods": 2}}))


def quiet(f, *a, **k):
    with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
        return f(*a, **k)


class PrepareCopies(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)
        self.src = self.root / "tut"
        make_tutorial(self.src)
        self.dest = self.root / "work/tol/womersleyTube"
        self.dest.parent.mkdir(parents=True)

    def tearDown(self):
        self.tmp.cleanup()

    def ensure(self, *args):
        return quiet(P.main, ["--ensure", str(self.src), str(self.dest), *args])

    def test_prepare_with_tolerance_and_marker(self):
        self.ensure("--coupling-tol-factor", "0.01")
        self.assertIn("1e-08", (self.dest / "constant/fsiProperties.robin").read_text())
        rec = json.loads((self.dest / P.MARKER).read_text())
        self.assertTrue(rec["complete"])
        self.assertEqual(rec["changes"], {"couplingTolFactor": 0.01})
        self.assertFalse(self.dest.with_name("womersleyTube.partial").exists())

    def test_failed_substitution_leaves_no_copy(self):
        (self.src / "constant/fsiProperties.robin").write_text("outerCorrTolerance 1e-6;\n")
        with self.assertRaises(SystemExit):
            self.ensure("--coupling-tol-factor", "0.01")
        self.assertFalse(self.dest.exists())       # half-prepared copy is never accepted
        (self.src / "constant/fsiProperties.robin").write_text(ROBIN)
        self.ensure("--coupling-tol-factor", "0.01")  # and the partial dir is cleaned up
        self.assertTrue(self.dest.is_dir())

    def test_current_copy_is_kept_untouched(self):
        self.ensure()
        keep = self.dest / "verification/work/c1"
        keep.mkdir(parents=True)
        (keep / "x").write_text("run")
        self.ensure()
        self.assertTrue((keep / "x").exists())

    def test_half_prepared_copy_without_marker_is_redone(self):
        self.dest.mkdir()
        (self.dest / "junk").write_text("x")
        self.ensure("--coupling-tol-factor", "0.01")
        self.assertFalse((self.dest / "junk").exists())
        self.assertTrue((self.dest / P.MARKER).is_file())

    def test_changed_source_or_settings_with_runs_is_refused_unless_forced(self):
        self.ensure()
        run = self.dest / "verification/work/c1"
        run.mkdir(parents=True)
        (self.src / "system/controlDict").write_text("a 2;\nfunctions\n{\n}\n")
        with self.assertRaises(SystemExit) as cm:
            self.ensure()
        self.assertEqual(cm.exception.code, 3)
        self.assertTrue(run.is_dir())
        with self.assertRaises(SystemExit) as cm:   # different settings, same source
            self.ensure("--coupling-tol-factor", "0.5")
        self.assertEqual(cm.exception.code, 3)
        self.ensure("--force")
        self.assertFalse(run.exists())
        self.assertIn("a 2;", (self.dest / "system/controlDict").read_text())

    def test_documentation_changes_do_not_make_a_copy_stale(self):
        self.ensure()
        run = self.dest / "verification/work/c1"
        run.mkdir(parents=True)
        (self.src / "README.md").write_text("docs")
        (self.src / "verification/traction").mkdir()
        (self.src / "verification/traction/README.md").write_text("docs")
        self.ensure()                      # still current: no exit, run kept
        self.assertTrue(run.is_dir())

    def test_changed_source_without_runs_is_refreshed(self):
        self.ensure()
        (self.src / "system/controlDict").write_text("a 3;\nfunctions\n{\n}\n")
        self.ensure()
        self.assertIn("a 3;", (self.dest / "system/controlDict").read_text())

    def test_diagnostic_comes_from_scripts_and_foreign_one_is_not_replaced_silently(self):
        self.ensure()
        installed = (self.dest / "system/interfaceTractionDiagnostic").read_bytes()
        self.assertEqual(installed.count(b"executeInterval 1;"), 1)
        (self.src / "system/interfaceTractionDiagnostic").write_text("other diagnostic\n")
        d2 = self.root / "work/other/womersleyTube"
        d2.parent.mkdir(parents=True)
        with self.assertRaises(SystemExit):
            quiet(P.main, [str(self.src), str(d2)])
        self.assertFalse(d2.exists())
        quiet(P.main, [str(self.src), str(d2), "--overwrite-diagnostic"])
        self.assertEqual((d2 / "system/interfaceTractionDiagnostic").read_bytes(), installed)


if __name__ == "__main__":
    unittest.main()
