#!/usr/bin/env python3
"""Unit tests of the bendingBeamInChannel verification driver, with
synthetic data, and consistency tests of the tracked reference files and
compact results (no OpenFOAM, no SLURM, no MPI).

    python3 -m unittest discover -s scripts -p 'test_*.py'
"""

import csv
import importlib.util
import io
import math
import shutil
import sys
import tempfile
import unittest
from contextlib import redirect_stdout
from pathlib import Path
from unittest import mock

SCRIPT = Path(__file__).resolve().parent / "bending_beam_verification.py"
spec = importlib.util.spec_from_file_location("bbv", SCRIPT)
bbv = importlib.util.module_from_spec(spec)
spec.loader.exec_module(bbv)

HEADER_CUTLINK = ("# Time\tforce (x y z)\ttorque (x y z)\tinertia of the "
                  "fluid inside (x y z)\tforce from the momentum exchange "
                  "(x y z)")
HEADER_FORCING = ("# Time\tforce (x y z)\ttorque (x y z)\tinertia of the "
                  "fluid inside (x y z)")


def force_line(t, *vectors):
    return f"{t}\t" + "\t".join(" ".join(str(c) for c in v)
                                for v in vectors)


class Sandbox(unittest.TestCase):
    """Each test runs with the work directories in a temporary directory"""

    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())
        self.patches = [
            mock.patch.object(bbv, "WORK_DIR", self.tmp / "work"),
            mock.patch.object(bbv, "CASES_DIR", self.tmp / "work/cases"),
            mock.patch.object(bbv, "MANIFEST", self.tmp / "work/manifest.json"),
        ]
        for p in self.patches:
            p.start()

    def tearDown(self):
        for p in self.patches:
            p.stop()
        shutil.rmtree(self.tmp, ignore_errors=True)


class TutorialParameters(unittest.TestCase):
    def test_parameters_are_read_from_the_tutorial(self):
        p = bbv.parameters()
        self.assertAlmostEqual(p["H"], 2.0295)
        self.assertAlmostEqual(p["Ax"], 0.4)
        self.assertAlmostEqual(p["Ay"], 0.8)
        self.assertAlmostEqual(p["Lz"], 0.1)
        self.assertAlmostEqual(p["U_ref"], 1.0, places=6)
        # Cd = 9.855 Fx of the tutorial README
        self.assertAlmostEqual(p["coeff_scale"], 9.855, places=2)


class AnalyticVolumeRate(unittest.TestCase):
    def setUp(self):
        self.p = bbv.parameters()

    def test_peak_is_the_documented_value(self):
        self.assertAlmostEqual(bbv.volume_rate_amplitude(self.p), 0.028306,
                               places=6)

    def test_closed_form_agrees_with_the_outline_area(self):
        peak = bbv.volume_rate_amplitude(self.p)
        for t in (0.03, 0.125, 0.2, 0.375, 0.61, 0.875):
            self.assertAlmostEqual(
                bbv.polygon_volume_rate(self.p, t) / peak,
                bbv.volume_rate(self.p, t) / peak, places=5)

    def test_zero_at_the_extremes_of_the_bending_and_peak_between(self):
        self.assertAlmostEqual(bbv.volume_rate(self.p, 0.25), 0.0, places=12)
        self.assertAlmostEqual(bbv.volume_rate(self.p, 0.125),
                               bbv.volume_rate_amplitude(self.p))
        self.assertAlmostEqual(bbv.volume_rate(self.p, 0.375),
                               -bbv.volume_rate_amplitude(self.p))

    def test_study_b_has_no_volume_change(self):
        self.assertEqual(bbv.volume_rate_amplitude(self.p, 0.0), 0.0)
        self.assertAlmostEqual(bbv.polygon_volume_rate(self.p, 0.125, 0.0),
                               0.0, places=9)

    def test_volume_rate_is_six_point_nine_percent_of_the_throughflow(self):
        through = self.p["U_ref"] * 4.1 * self.p["Lz"]
        self.assertAlmostEqual(bbv.volume_rate_amplitude(self.p) / through,
                               0.0690, places=3)


class Plan(unittest.TestCase):
    class Args:
        runs = ["smoke", "L1", "L2", "studyB"]
        ranks = 0

    def test_plan_has_unique_named_cases(self):
        names = [c[4] for c in bbv.plan(self.Args)]
        self.assertEqual(len(names), len(set(names)))
        self.assertIn("cutLink-L1-A-t0p5", names)

    def test_incremental_keeps_its_five_outer_correctors(self):
        self.assertEqual(bbv.N_OUTER["incremental"], 5)
        self.assertEqual({bbv.N_OUTER[m] for m in ("cutLink", "penalty",
                                                   "ghostCell")}, {1})

    def test_unknown_run_is_rejected(self):
        class A(self.Args):
            runs = ["nope"]
        with self.assertRaises(SystemExit):
            bbv.plan(A)

    def test_cells_follow_the_blockmesh_formula(self):
        for level, cells in bbv.CELLS.items():
            n = 2 ** (level - 1)
            self.assertEqual(cells, (2 * 80 * n + 100 * n) * 77 * n)

    def test_default_run_is_short(self):
        self.assertEqual(bbv.DEFAULT_RUNS, ("smoke",))
        self.assertTrue(all(c[3] <= 0.5 for c in bbv.RUNS["smoke"][1]))

    def test_dry_run_prints_cost_and_writes_nothing(self):
        class A(Plan.Args):
            runs = ["L2"]
            ranks = 0
        out = io.StringIO()
        with redirect_stdout(out):
            self.assertEqual(bbv.cmd_dry_run(A), 0)
        self.assertIn("allocated rank-hours", out.getvalue())
        self.assertFalse(bbv.CASES_DIR.exists())


class CaseGeneration(Sandbox):
    def generate(self, method, study="A", ranks=1, level=1, end=0.5):
        return bbv.generate_case(method, level, study, end,
                                 bbv.case_name(method, level, study, end),
                                 ranks)

    def test_tutorial_is_not_modified_and_work_files_are_not_copied(self):
        before = (bbv.TUTORIAL_DIR / "constant/fvOptions").read_text()
        case = self.generate("penalty")
        self.assertEqual((bbv.TUTORIAL_DIR / "constant/fvOptions").read_text(),
                         before)
        self.assertTrue((case / "0/U").exists())
        self.assertFalse((case / "verification").exists())

    def test_penalty_settings_are_written_explicitly(self):
        text = (self.generate("penalty") / "constant/fvOptions").read_text()
        for entry in ("method          penalty", "penaltyCoeff    1e3",
                      "volumeFraction", "surfaceRateCoeff 3",
                      "signedDistance"):
            self.assertIn(entry.split()[0], text)
        self.assertIn("signedDistance", text)
        self.assertNotIn("forceEstimator", text)
        self.assertIn("amplitudeY  0.8", text)

    def test_cutlink_keeps_the_surface_traction_estimator(self):
        text = (self.generate("cutLink") / "constant/fvOptions").read_text()
        self.assertIn("forceEstimator", text)
        self.assertNotIn("tractionQuadrature", text)

    def test_incremental_has_five_outer_correctors_and_vertex_occupancy(self):
        case = self.generate("incremental")
        self.assertIn("vertexFraction", (case / "constant/fvOptions")
                      .read_text())
        self.assertRegex((case / "system/fvSolution").read_text(),
                         r"nOuterCorrectors\s+5;")

    def test_study_b_sets_amplitude_y_to_zero_only(self):
        a = (self.generate("cutLink", "A") / "constant/fvOptions").read_text()
        b = (self.generate("cutLink", "B") / "constant/fvOptions").read_text()
        self.assertIn("amplitudeY  0.8", a)
        self.assertIn("amplitudeY  0;", b)
        self.assertEqual(a.replace("amplitudeY  0.8;", "amplitudeY  0;"), b)

    def test_end_time_write_interval_and_decomposition(self):
        case = self.generate("cutLink", ranks=4, end=2.0)
        control = (case / "system/controlDict").read_text()
        self.assertRegex(control, r"endTime\s+2;")
        self.assertRegex(control, r"writeInterval\s+0.125;")
        self.assertIn("(4 1 1)", (case / "system/decomposeParDict")
                      .read_text())
        script = (case / "Allrun.verify").read_text()
        self.assertIn("decomposePar", script)
        self.assertIn("reconstructPar", script)
        self.assertIn("hwloc_base_binding_policy=none", script)

    def test_serial_case_has_no_mpi(self):
        script = (self.generate("cutLink", ranks=1) / "Allrun.verify") \
            .read_text()
        self.assertNotIn("runParallel", script)

    def test_existing_case_is_not_overwritten_without_force(self):
        self.generate("cutLink")
        with self.assertRaises(SystemExit):
            self.generate("cutLink")

    def test_input_hashes_change_with_the_inputs(self):
        a = bbv.input_hashes(self.generate("cutLink"))
        b = bbv.input_hashes(self.generate("penalty"))
        self.assertNotEqual(a["constant/fvOptions"], b["constant/fvOptions"])
        self.assertEqual(a["system/fvSchemes"], b["system/fvSchemes"])


class Forces(Sandbox):
    def write(self, header, lines):
        path = self.tmp / "beam.dat"
        path.write_text("# Immersed body beam\n" + header + "\n"
                        + "\n".join(lines) + "\n")
        return path

    def test_cutlink_estimators_are_kept_apart(self):
        path = self.write(HEADER_CUTLINK, [
            force_line(0.1, (1, 2, 0), (0, 0, 0), (9, 9, 9), (3, 4, 0))])
        series = bbv.hydrodynamic_forces(*bbv.read_force_file(path), 10.0)
        self.assertEqual(series["cutLink ST"], ([10.0], [20.0]))
        self.assertEqual(series["cutLink ME"], ([30.0], [40.0]))
        self.assertNotIn("total", series)

    def test_forcing_plus_inertia_is_the_hydrodynamic_force(self):
        path = self.write(HEADER_FORCING, [
            force_line(0.1, (1, 2, 0), (0, 0, 0), (5, -3, 0))])
        series = bbv.hydrodynamic_forces(*bbv.read_force_file(path), 2.0)
        self.assertEqual(series["total"], ([12.0], [-2.0]))
        self.assertEqual(series["forcing"], ([2.0], [4.0]))

    def test_unrecognised_columns_are_rejected(self):
        path = self.write("# Time\tforce (x y z)", ["0.1\t1 2 3"])
        with self.assertRaises(SystemExit):
            bbv.hydrodynamic_forces(*bbv.read_force_file(path), 1.0)


class SolverLog(Sandbox):
    def log(self, text):
        case = self.tmp / "case"
        case.mkdir(exist_ok=True)
        (case / "log.solids4Foam").write_text(text)
        return case

    def test_normal_end_with_the_startup_trapping_message_is_valid(self):
        case = self.log("sigFpe : Enabling floating point exception "
                        "trapping (FOAM_SIGFPE).\nTime = 1\n\nEnd\n")
        self.assertEqual(bbv.solver_problem(case), "")

    def test_nan_missing_end_and_missing_log_are_problems(self):
        self.assertIn("nan", bbv.solver_problem(self.log("Time = 1\nnan\n"
                                                         "\nEnd\n")))
        self.assertIn("End", bbv.solver_problem(self.log("Time = 1\n")))
        self.assertIn("no solver log",
                      bbv.solver_problem(self.tmp / "missing"))


class Metrics(unittest.TestCase):
    def series(self, fn, step=0.005, end=8.0):
        t = [i * step for i in range(int(end / step) + 1)]
        return t, [fn(x) for x in t]

    def test_nrms_of_identical_series_is_zero(self):
        t, y = self.series(lambda x: 5 * math.sin(2 * math.pi * x))
        self.assertAlmostEqual(bbv.nrms(t, y, t, y), 0.0)

    def test_nrms_of_a_scaled_series_is_the_scale_error(self):
        t, y = self.series(lambda x: 5 * math.sin(2 * math.pi * x))
        self.assertAlmostEqual(bbv.nrms(t, [1.1 * v for v in y], t, y), 10.0,
                               places=6)

    def test_harmonics_recover_amplitudes_without_phase_alignment(self):
        t, y = self.series(lambda x: 3 + 70 * math.sin(2 * math.pi * x - .6)
                           + 6 * math.sin(4 * math.pi * x + 1.1))
        h = bbv.harmonics(t, y)
        self.assertAlmostEqual(h["mean"], 3.0, delta=0.02)
        self.assertAlmostEqual(h["A1"], 70.0, delta=0.05)
        self.assertAlmostEqual(h["A2"], 6.0, delta=0.02)
        self.assertAlmostEqual(h["ph2"], math.degrees(1.1), places=1)

    def test_a_missing_second_harmonic_is_seen(self):
        t, y = self.series(lambda x: 70 * math.sin(2 * math.pi * x))
        self.assertLess(bbv.harmonics(t, y)["A2"], 1e-6)

    def test_noise_of_a_smooth_series_is_small_and_of_a_noisy_one_not(self):
        t, y = self.series(lambda x: 10 * math.sin(2 * math.pi * x), 0.002)
        smooth = bbv.noise(t, y)
        noisy = bbv.noise(t, [v + (1 if i % 2 else -1)
                              for i, v in enumerate(y)])
        self.assertLess(smooth, 0.05)
        self.assertGreater(noisy, 0.5)

    def test_interpolation_and_trapezoid(self):
        self.assertAlmostEqual(bbv.interpolate([0, 1], [0, 2], 0.25), 0.5)
        self.assertAlmostEqual(bbv.trapezoid([0, 1, 2], [0, 1, 0]), 1.0)

    def test_linear_solver(self):
        x = bbv.solve_linear([[2.0, 1.0], [1.0, 3.0]], [3.0, 5.0])
        self.assertAlmostEqual(x[0], 0.8)
        self.assertAlmostEqual(x[1], 1.4)


PHI = """\
FoamFile {{ version 2.0; format ascii; class surfaceScalarField; object phi; }}
internalField   nonuniform List<scalar> 2(1 2);
boundaryField
{{
    inlet
    {{
        type            calculated;
        value           nonuniform List<scalar>
3
(
-0.01
-0.02
-0.03
)
;
    }}
    outlet
    {{
        type            calculated;
        value           nonuniform List<scalar>
2
(
{outlet}
{outlet2}
)
;
    }}
    walls
    {{
        type            calculated;
        value           uniform 0;
    }}
}}
"""


class Flux(Sandbox):
    def make_case(self, outlets, times):
        case = self.tmp / "case"
        for t, (a, b) in zip(times, outlets):
            (case / str(t)).mkdir(parents=True)
            (case / str(t) / "phi").write_text(PHI.format(outlet=a,
                                                          outlet2=b))
        return case

    def test_patch_flux_sums_inlet_and_outlet_only(self):
        case = self.make_case([(0.04, 0.02)], [0.125])
        self.assertAlmostEqual(bbv.patch_flux(case / "0.125/phi"),
                               -0.06 + 0.06)
        self.assertAlmostEqual(bbv.patch_flux(case / "0.125/phi", ("inlet",)),
                               -0.06)

    def test_flux_report_against_the_exact_rate(self):
        p = bbv.parameters()
        peak = bbv.volume_rate_amplitude(p)
        # inlet flux -0.06 (the fixture); the outlet carries the rest
        a = [0.5 * (0.06 + q * 0.95 * peak) for q in (1, -1, 1)]
        case = self.make_case([(x, x) for x in a], [0.125, 0.375, 0.625])
        rep = bbv.flux_report(case, p, "A")
        self.assertEqual(rep["snapshots"], 3)
        self.assertAlmostEqual(rep["ratio_at_peaks"], 0.95, places=6)
        self.assertAlmostEqual(rep["peak_residual_pct"], 5.0, places=4)

    def test_study_b_reports_the_absolute_flux_only(self):
        p = bbv.parameters()
        case = self.make_case([(0.03, 0.03)], [0.125])
        rep = bbv.flux_report(case, p, "B")
        self.assertEqual(rep["peak_analytic"], 0.0)
        self.assertNotIn("ratio_at_peaks", rep)
        self.assertAlmostEqual(rep["max_abs_Q"], 0.0, places=9)

    def test_times_off_the_write_grid_are_ignored(self):
        p = bbv.parameters()
        case = self.make_case([(0.03, 0.03)], [0.13])
        self.assertEqual(bbv.flux_report(case, p, "B")["snapshots"], 0)


class CheckSemantics(unittest.TestCase):
    def row(self, **kw):
        r = {"method": "cutLink", "level": 2, "study": "A",
             "estimator": "cutLink ST", "coefficient": "Cd", "end_time": 8.0,
             "nrms_pct": 5.0, "A2": 5.7, "flux_snapshots": 4,
             "flux_peak_analytic": 0.0283, "flux_max_abs_Q": 0.0269,
             "flux_ratio_at_peaks": 0.9504}
        r.update(kw)
        return r

    TRACKED = [{"method": "cutLink", "level": "2", "study": "A",
                "estimator": "cutLink ST", "coefficient": "Cd",
                "nrms_pct": "5.04", "A2": "5.72",
                "flux_ratio_at_peaks": "0.9504"}]

    def statuses(self, row, tracked=None):
        return [s for s, _, _ in bbv.check_case(
            row, self.TRACKED if tracked is None else tracked, True)]

    def test_pass(self):
        self.assertEqual(self.statuses(self.row()), [bbv.PASS, bbv.PASS])

    def test_nrms_outside_tolerance_fails(self):
        self.assertEqual(self.statuses(self.row(nrms_pct=6.5))[0], bbv.FAIL)

    def test_second_harmonic_outside_tolerance_fails(self):
        self.assertEqual(self.statuses(self.row(A2=5.4))[0], bbv.FAIL)

    def test_flux_ratio_outside_tolerance_fails(self):
        self.assertEqual(self.statuses(self.row(flux_ratio_at_peaks=0.9))[1],
                         bbv.FAIL)

    def test_missing_tracked_result_is_skipped_not_passed(self):
        self.assertIn(bbv.SKIPPED, self.statuses(self.row(), tracked=[]))

    def test_short_run_is_skipped_for_the_force_comparison(self):
        row = self.row(end_time=0.5)
        del row["nrms_pct"], row["A2"]
        self.assertEqual(self.statuses(row)[0], bbv.SKIPPED)

    def test_penalty_must_have_zero_net_exterior_flux(self):
        row = self.row(method="penalty", estimator="total",
                       flux_max_abs_Q=1e-9)
        tracked = [dict(self.TRACKED[0], method="penalty",
                        estimator="total")]
        self.assertEqual(self.statuses(row, tracked), [bbv.PASS, bbv.PASS])
        bad = self.row(method="penalty", estimator="total",
                       flux_max_abs_Q=1e-3)
        self.assertEqual(self.statuses(bad, tracked)[1], bbv.FAIL)

    def test_missing_phi_is_skipped(self):
        row = self.row(flux_snapshots=0)
        self.assertEqual(self.statuses(row)[1], bbv.SKIPPED)

    def test_analytic_checks_pass(self):
        self.assertEqual({s for s, _, _ in bbv.check_analytic(
            bbv.parameters())}, {bbv.PASS})


class MpiSetup(unittest.TestCase):
    def test_missing_openfoam_is_reported(self):
        with mock.patch.object(shutil, "which", return_value=None):
            self.assertIn("unavailable", bbv.check_environment(1))

    def test_more_ranks_than_cores_is_reported(self):
        with mock.patch.object(shutil, "which", return_value="/bin/true"), \
                mock.patch.object(bbv.os, "sched_getaffinity",
                                  return_value={0, 1}):
            self.assertIn("share cores", bbv.check_environment(4))
            self.assertEqual(bbv.check_environment(2), "")


class TrackedData(unittest.TestCase):
    """The tracked files agree with each other"""

    def rows(self):
        with (bbv.RESULTS_DIR / "cases.csv").open() as handle:
            return list(csv.DictReader(handle))

    def test_reference_series_is_sampled_every_hundredth_of_a_second(self):
        t, cd, cl = bbv.read_series(bbv.REFERENCE_FILE)
        self.assertAlmostEqual(t[1] - t[0], 0.01, places=6)
        self.assertGreater(t[-1], 7.99)
        self.assertEqual(len(t), len(cd))

    def test_tracked_cases_have_the_documented_levels(self):
        levels = {(r["method"], int(r["level"]), r["study"])
                  for r in self.rows()}
        for key in (("cutLink", 1, "A"), ("cutLink", 3, "A"),
                    ("cutLink", 4, "A"), ("penalty", 3, "A"),
                    ("incremental", 2, "A"), ("cutLink", 3, "B")):
            self.assertIn(key, levels)

    def test_cutlink_flux_residual_halves_from_level_two(self):
        res = {int(r["level"]): float(r["flux_peak_residual_pct"])
               for r in self.rows() if r["method"] == "cutLink"
               and r["study"] == "A" and r["estimator"] == "cutLink ST"}
        self.assertAlmostEqual(res[2] / res[3], 2.0, delta=0.1)
        self.assertAlmostEqual(res[3] / res[4], 2.0, delta=0.1)
        self.assertLess(res[1], 8.0)

    def test_penalty_and_incremental_flux_is_zero(self):
        for r in self.rows():
            if r["method"] in ("penalty", "incremental", "ghostCell") \
                    and r["study"] == "A":
                self.assertLess(float(r["flux_max_abs_Q"]), 1e-7)

    def test_histories_match_the_tracked_means(self):
        for r in self.rows():
            if r["study"] != "A" or r["coefficient"] != "Cd" \
                    or not r["nrms_pct"] or int(r["level"]) > 3 \
                    or r["method"] == "ghostCell":
                continue
            name = {"cutLink ST": "Cd_cutLink_ST",
                    "total": "Cd_total"}.get(r["estimator"])
            if name is None:
                continue
            path = bbv.RESULTS_DIR / f"history-{r['method']}-L{r['level']}.csv"
            with path.open() as handle:
                data = list(csv.DictReader(handle))
            t = [float(x["t"]) for x in data]
            y = [float(x[name]) for x in data]
            inside = [v for x, v in zip(t, y) if 1 < x < 8]
            mean = sum(inside) / len(inside)
            # 0.02 s display data: the high-frequency content of the raw
            # series is undersampled, so only rough agreement is expected
            self.assertAlmostEqual(mean, float(r["mean"]), delta=0.5)

    def test_costs_are_positive_or_not_recorded(self):
        for r in self.rows():
            cost = float(r["allocated_rank_hours"])
            self.assertTrue(math.isnan(cost) or cost > 0)

    def test_nrms_to_the_reference_is_reproduced_from_the_history(self):
        t_ref, cd_ref, _ = bbv.read_series(bbv.REFERENCE_FILE)
        path = bbv.RESULTS_DIR / "history-cutLink-L2.csv"
        with path.open() as handle:
            data = list(csv.DictReader(handle))
        t = [float(x["t"]) for x in data]
        y = [float(x["Cd_cutLink_ST"]) for x in data]
        row = [r for r in self.rows() if r["method"] == "cutLink"
               and r["level"] == "2" and r["study"] == "A"
               and r["estimator"] == "cutLink ST"
               and r["coefficient"] == "Cd"][0]
        # the 0.02 s history is display data: the tracked nRMS is computed
        # from the complete series, so only rough agreement is expected
        self.assertAlmostEqual(bbv.nrms(t, y, t_ref, cd_ref),
                               float(row["nrms_pct"]), delta=0.5)


class Documentation(unittest.TestCase):
    """The numbers displayed in README.md and results.md are those of the
    tracked data, and the figures they include exist"""

    def rows(self):
        with (bbv.RESULTS_DIR / "cases.csv").open() as handle:
            return list(csv.DictReader(handle))

    def test_results_md_shows_the_tracked_nrms_and_harmonics(self):
        text = (bbv.VERIFICATION_DIR / "results.md").read_text()
        for r in self.rows():
            if r["study"] == "A" and r["nrms_pct"]:
                self.assertIn(f"{float(r['nrms_pct']):.2f}", text)
                self.assertIn(f"{float(r['A2']):.3f}", text)

    def test_results_md_shows_the_tracked_flux(self):
        text = (bbv.VERIFICATION_DIR / "results.md").read_text()
        for r in self.rows():
            if r["method"] == "cutLink" and r["study"] == "A" \
                    and r["estimator"] == "cutLink ST":
                self.assertIn(f"{float(r['flux_ratio_at_peaks']):.4f}", text)
                self.assertIn(f"{float(r['flux_peak_residual_pct']):.2f}",
                              text)

    def test_readme_table_agrees_with_the_data(self):
        text = (bbv.VERIFICATION_DIR / "README.md").read_text()
        by = {(r["method"], int(r["level"]), r["estimator"]): r
              for r in self.rows()
              if r["study"] == "A" and r["coefficient"] == "Cd"}
        line = "| cutLink, surface traction | " + " / ".join(
            f"{float(by[('cutLink', l, 'cutLink ST')]['nrms_pct']):.1f}"
            for l in (1, 2, 3)) + " % |"
        self.assertIn(line, text)
        for l in (1, 2, 3, 4):
            ratio = [r for r in self.rows() if r["method"] == "cutLink"
                     and r["study"] == "A" and int(r["level"]) == l][0]
            self.assertIn(f"{float(ratio['flux_ratio_at_peaks']):.3f}", text)

    def test_included_figures_exist(self):
        for name in ("README.md", "results.md"):
            text = (bbv.VERIFICATION_DIR / name).read_text()
            for image in set(__import__("re").findall(r"images/(\w+\.png)",
                                                       text)):
                self.assertTrue((bbv.VERIFICATION_DIR / "images" / image)
                                .exists(), image)

    def test_cost_pairs_support_the_stated_range(self):
        with (bbv.RESULTS_DIR / "cost-matched-pairs.csv").open() as handle:
            pairs = list(csv.DictReader(handle))
        ratios = [float(p["ratio"]) for p in pairs]
        text = (bbv.VERIFICATION_DIR / "results.md").read_text()
        self.assertIn(f"{min(ratios):.2f} to {max(ratios):.2f}", text)
        self.assertTrue(all(p["both_windows_after_affinity_fix"] == "True"
                            for p in pairs))

    def test_no_sensitivity_is_marked_material(self):
        with (bbv.RESULTS_DIR / "fv-sensitivity.csv").open() as handle:
            rows = list(csv.DictReader(handle))
        self.assertTrue(rows)
        self.assertEqual({r["material"] for r in rows}, {"False"})


class Wrapper(unittest.TestCase):
    def test_allverify_and_allclean_exist_and_are_executable(self):
        for name in ("Allverify", "Allclean"):
            path = bbv.VERIFICATION_DIR / name
            self.assertTrue(path.exists())
            self.assertTrue(path.stat().st_mode & 0o111)

    def test_unit_command_runs_without_openfoam(self):
        out = io.StringIO()
        with redirect_stdout(out):
            self.assertEqual(bbv.main(["unit"]), 0)
        self.assertIn("PASS", out.getvalue())

    def test_no_machine_specific_paths(self):
        text = SCRIPT.read_text()
        for needle in ("/home/", "/mnt/", "sairam"):
            self.assertNotIn(needle, text)


if __name__ == "__main__":
    unittest.main(argv=sys.argv[:1] + sys.argv[2:])
