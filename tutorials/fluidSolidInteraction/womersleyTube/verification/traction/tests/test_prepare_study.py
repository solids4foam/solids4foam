"""prepare_study.py: counted tolerance substitutions, atomic/marked study
copies, reuse only of matching copies, no silent deletion of runs, the
diagnostic taken from the scripts directory."""
import contextlib
import io
import json
import os
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

    def test_periods_option(self):
        diag = Path(P.DIAGNOSTIC).read_text()
        self.assertIn("timeStart", diag)
        self.ensure("--periods", "10")
        ref = json.loads((self.dest / "verification/reference/"
                          "womersleyTube_verification_references.json").read_text())
        self.assertEqual(ref["study"]["periods"], 10)
        self.assertEqual(ref["study"]["analysisPeriods"], 2)   # window: periods 9-10
        import re
        start = re.search(r"^\s*timeStart\s+([^;]+);",
                          (self.dest / "system/interfaceTractionDiagnostic").read_text(),
                          re.MULTILINE).group(1)
        self.assertEqual(float(start), 8.0)                    # (10 - 2) periods of 1.0
        rec = json.loads((self.dest / P.MARKER).read_text())
        self.assertEqual(rec["changes"], {"periods": 10})
        self.assertEqual(rec["requested"]["periods"], 10)
        # the source tutorial is untouched and a copy without --periods is another setting
        self.assertEqual(json.loads((self.src / "verification/reference/"
                                     "womersleyTube_verification_references.json"
                                     ).read_text())["study"]["periods"], 6)
        self.dest.joinpath("verification/work/c1").mkdir(parents=True)
        (self.dest / "verification/work/c1/x").write_text("run")
        with self.assertRaises(SystemExit):
            self.ensure()                                      # holds runs: refused

    def test_periods_with_negative_amplitude_on_the_real_tutorial(self):
        real = Path(P.__file__).resolve().parents[3]
        dest = self.root / "real/womersleyTube"
        quiet(P.main, [str(real), str(dest), "--amplitude=-1", "--periods", "10"])
        ref = json.loads((dest / "verification/reference/"
                          "womersleyTube_verification_references.json").read_text())
        self.assertEqual(ref["study"]["periods"], 10)
        self.assertEqual(ref["reference"]["parameters"]["pressureAmplitude"], -1.0)
        rec = json.loads((dest / P.MARKER).read_text())
        self.assertEqual(rec["changes"], {"pressureAmplitude": -1.0, "periods": 10})

    def test_periods_must_exceed_the_analysis_window(self):
        for bad in ("2", "1", "0"):
            with self.assertRaises(SystemExit):
                self.ensure("--periods", bad)
            self.assertFalse(self.dest.exists())

    def test_without_periods_the_marker_has_no_periods_key(self):
        self.ensure()
        self.assertNotIn("periods", json.loads((self.dest / P.MARKER).read_text())["requested"])

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


class AtomicPublication(unittest.TestCase):
    """Interrupted, failing and concurrent-with-runs preparations."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)
        self.src = self.root / "tut"
        make_tutorial(self.src)
        (self.src / "verification/scripts").mkdir()
        (self.src / "verification/scripts/driver.py").write_text("# analysis v1\n")
        self.dest = self.root / "work/a/womersleyTube"
        self.dest.parent.mkdir(parents=True)

    def tearDown(self):
        self.tmp.cleanup()

    def ensure(self, *args):
        return quiet(P.main, ["--ensure", str(self.src), str(self.dest), *args])

    def make_run(self):
        run = self.dest / "verification/work/c1"
        run.mkdir(parents=True)
        (run / "x").write_text("run")
        return run

    def leftovers(self):
        return sorted(p.name for p in self.dest.parent.iterdir())

    def test_no_dest_is_touched_before_the_new_copy_is_complete(self):
        self.ensure()
        (self.dest / "marker_of_old").write_text("old")
        (self.src / "system/controlDict").write_text("a 5;\nfunctions\n{\n}\n")
        real = P.validate_copy

        def broken(*a, **k):
            raise KeyboardInterrupt        # interrupted during validation
        P.validate_copy = broken
        try:
            with self.assertRaises(KeyboardInterrupt):
                self.ensure()
        finally:
            P.validate_copy = real
        self.assertTrue((self.dest / "marker_of_old").exists())   # old copy intact
        self.assertEqual(self.leftovers(), ["womersleyTube"])      # staging cleaned
        self.ensure()
        self.assertIn("a 5;", (self.dest / "system/controlDict").read_text())

    def test_interrupted_between_the_two_renames_is_recovered(self):
        self.ensure()
        run = self.make_run()
        (self.src / "system/controlDict").write_text("a 6;\nfunctions\n{\n}\n")
        real = os.replace

        def crash(a, b):
            raise KeyboardInterrupt        # dest already renamed aside
        P.os.replace = crash
        try:
            with self.assertRaises(KeyboardInterrupt):
                self.ensure("--force")
        finally:
            P.os.replace = real
        # publish() undid its own rename: the previous complete copy is back
        self.assertTrue(run.is_dir())
        self.assertEqual(self.leftovers(), ["womersleyTube"])
        # a hard kill in the same window leaves DEST absent and DEST.old-*:
        old = self.dest.with_name("womersleyTube.old-20260101T000000-1")
        os.rename(self.dest, old)
        self.assertFalse(self.dest.exists())
        with self.assertRaises(SystemExit) as cm:      # restored, then judged: stale with runs
            self.ensure()
        self.assertEqual(cm.exception.code, 3)
        self.assertTrue(run.is_dir())
        self.assertEqual(self.leftovers(), ["womersleyTube"])

    def test_leftover_staging_directory_is_discarded(self):
        part = self.dest.with_name("womersleyTube.partial")
        part.mkdir()
        (part / "junk").write_text("x")
        self.ensure()
        self.assertEqual(self.leftovers(), ["womersleyTube"])
        self.assertTrue(json.loads((self.dest / P.MARKER).read_text())["complete"])

    def test_copy_without_marker_is_rebuilt_or_refused_when_it_has_runs(self):
        self.ensure()
        (self.dest / P.MARKER).unlink()                  # an incomplete copy
        self.ensure()                                    # no runs: rebuilt
        self.assertTrue((self.dest / P.MARKER).is_file())
        (self.dest / P.MARKER).unlink()
        run = self.make_run()
        with self.assertRaises(SystemExit) as cm:
            self.ensure()
        self.assertEqual(cm.exception.code, 3)
        self.assertTrue(run.is_dir())
        self.ensure("--force")
        self.assertFalse(run.exists())
        self.assertEqual(self.leftovers(), ["womersleyTube"])   # old one deleted after success

    def test_changed_analysis_files_refresh_the_copy_and_keep_the_runs(self):
        self.ensure()
        run = self.make_run()
        (self.src / "verification/scripts/driver.py").write_text("# analysis v2\n")
        self.ensure()
        self.assertEqual((self.dest / "verification/scripts/driver.py").read_text(),
                         "# analysis v2\n")
        self.assertTrue((run / "x").exists())
        self.assertEqual(self.leftovers(), ["womersleyTube"])
        self.ensure()                                    # and is current now

    def test_missing_file_after_copy_is_detected(self):
        real = P.shutil.copytree

        def lossy(src, dst, **kw):
            P.shutil.copytree = real          # copytree recurses through this name
            real(src, dst, **kw)
            (Path(dst) / "constant/fsiProperties.iqnils").unlink()
        P.shutil.copytree = lossy
        try:
            with self.assertRaises(SystemExit):
                self.ensure()
        finally:
            P.shutil.copytree = real
        self.assertFalse(self.dest.exists())
        self.assertEqual(self.leftovers(), [])

    def test_marker_is_written_last(self):
        order = []
        real_val = P.validate_copy

        def spy(src, d, req):
            order.append(("validate", (d / P.MARKER).exists()))
            return real_val(src, d, req)
        P.validate_copy = spy
        try:
            self.ensure()
        finally:
            P.validate_copy = real_val
        self.assertEqual(order, [("validate", False)])
        self.assertTrue((self.dest / P.MARKER).exists())


if __name__ == "__main__":
    unittest.main()
