#!/usr/bin/env python3
"""Exact tests of immersed_taylor_couette_verification.py: they need neither
OpenFOAM nor numpy. Run with ./Allverify unit."""

import contextlib
import csv
import io
import math
import os
import re
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

sys.path.insert(0, str(Path(__file__).resolve().parent))
import immersed_taylor_couette_verification as v  # noqa: E402

P = v.parameters()


def synthetic_field(level, func):
    n, h, centres = v.cell_centres(level)
    out = []
    for x, y in centres:
        r = math.hypot(x, y)
        ut = func(r)
        out.append((-ut * y / r, ut * x / r, 0.0))
    return out


def write_field(path, field):
    path.parent.mkdir(parents=True, exist_ok=True)
    body = "\n".join("(%.10g %.10g %.10g)" % u for u in field)
    path.write_text("FoamFile\n{\n}\n\ninternalField   nonuniform List<vector> "
                    f"\n{len(field)}\n(\n{body}\n)\n;\n\nboundaryField\n{{\n}}\n")


def write_forces(path, times, torque):
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = ["# Immersed body x",
             "# Time\tforce (x y z)\ttorque (x y z)\tinertia of the fluid "
             "inside (x y z)"]
    for t, tz in zip(times, torque):
        lines.append(f"{t}\t0 0 0\t0 0 {tz}\t0 0 0")
    path.write_text("\n".join(lines) + "\n")


class ParameterTests(unittest.TestCase):
    def test_parameters_are_read_from_the_tutorial(self):
        self.assertAlmostEqual(P["R1"], 0.05, places=12)
        self.assertAlmostEqual(P["R2"], 0.1, places=12)
        self.assertEqual(P["omega"], 1.0)
        self.assertEqual(P["nu"], 1e-3)
        self.assertEqual(P["rho"], 1.0)
        self.assertAlmostEqual(P["Lz"], 0.1, places=12)
        self.assertEqual(P["axis"], (0.0, 0.0, 1.0))
        self.assertEqual(P["origin"], (0.0, 0.0, 0.0))
        self.assertEqual(P["CofR"], (0.0, 0.0, 0.0))

    def test_stl_radii(self):
        r_min, r_max, z_min, z_max = v.stl_radii(
            v.TUTORIAL_DIR / "constant/triSurface/rotor.stl")
        self.assertAlmostEqual(r_max, 0.05, places=9)
        self.assertEqual(r_min, 0.0)
        self.assertLess(z_min, 0.0)
        self.assertGreater(z_max, 0.1)

    def test_dict_scalar_skips_the_dimensions(self):
        text = "nu              nu [0 2 -1 0 0 0 0] 1e-3;\nrho rho [1 -3 0 0 0 0 0] 1;"
        self.assertEqual(v.dict_scalar(text, "nu"), 1e-3)
        self.assertEqual(v.dict_scalar(text, "rho"), 1.0)
        with self.assertRaises(SystemExit), contextlib.redirect_stderr(io.StringIO()):
            v.dict_scalar(text, "mu")


class ExactSolutionTests(unittest.TestCase):
    def test_boundary_conditions(self):
        self.assertAlmostEqual(v.exact_utheta(P, P["R1"]), P["omega"] * P["R1"], 14)
        self.assertAlmostEqual(v.exact_utheta(P, P["R2"]), 0.0, 14)

    def test_satisfies_the_azimuthal_momentum_equation(self):
        # u'' + u'/r - u/r^2 = 0 (steady, no pressure gradient in theta)
        h = 1e-5
        for r in (0.055, 0.07, 0.085, 0.095):
            up = (v.exact_utheta(P, r + h) - v.exact_utheta(P, r - h)) / (2 * h)
            upp = (v.exact_utheta(P, r + h) - 2 * v.exact_utheta(P, r)
                   + v.exact_utheta(P, r - h)) / h ** 2
            res = upp + up / r - v.exact_utheta(P, r) / r ** 2
            self.assertAlmostEqual(res, 0.0, delta=1e-4)

    def test_shear_is_mu_r_d_over_dr_of_u_over_r(self):
        h = 1e-6
        for r in (0.06, 0.08, 0.099):
            d = (v.exact_utheta(P, r + h) / (r + h)
                 - v.exact_utheta(P, r - h) / (r - h)) / (2 * h)
            self.assertAlmostEqual(P["rho"] * P["nu"] * r * d,
                                   v.exact_shear(P, r), delta=1e-9)

    def test_angular_momentum_flux_r2_tau_is_constant(self):
        values = [r * r * v.exact_shear(P, r) for r in (0.051, 0.07, 0.09, 0.0999)]
        for x in values:
            self.assertAlmostEqual(x, values[0], delta=1e-15)

    def test_torque_value_sign_and_closed_form(self):
        t = v.exact_rotor_torque(P)
        self.assertAlmostEqual(t, -4.188790205e-6, delta=1e-14)
        self.assertLess(t, 0.0)
        self.assertAlmostEqual(
            t, -4 * math.pi * 1e-3 * 1.0 * 0.05 ** 2 * 0.1 ** 2
            / (0.1 ** 2 - 0.05 ** 2) * 0.1, delta=1e-18)
        self.assertAlmostEqual(v.integrated_rotor_torque(P), t, delta=1e-18)

    def test_torque_scalings_and_sign_with_the_rotation(self):
        base = v.exact_rotor_torque(P)
        for key, factor in (("rho", 2.0), ("nu", 3.0), ("omega", 0.5),
                            ("Lz", 1.5)):
            q = dict(P)
            q[key] = P[key] * factor
            self.assertAlmostEqual(v.exact_rotor_torque(q) / base, factor, 12)
        q = dict(P)
        q["omega"] = -1.0
        self.assertAlmostEqual(v.exact_rotor_torque(q), -base, 20)

    def test_stator_torque_balances_the_rotor(self):
        self.assertEqual(v.exact_stator_torque(P), -v.exact_rotor_torque(P))

    def test_exact_torque_equals_the_viscous_dissipation_balance(self):
        # power input |T| omega = viscous dissipation of the gap, the
        # integral of tau^2/mu over the volume (an independent integral)
        n, a, b = 20000, P["R1"], P["R2"]
        h = (b - a) / n
        mu = P["rho"] * P["nu"]
        power = 0.0
        for i in range(n):
            r = a + (i + 0.5) * h
            power += v.exact_shear(P, r) ** 2 / mu * 2 * math.pi * r * h
        power *= P["Lz"]
        self.assertAlmostEqual(power, abs(v.exact_rotor_torque(P)) * P["omega"],
                               delta=1e-6 * power)


class StatisticsTests(unittest.TestCase):
    def test_constant_series(self):
        t = [0.0, 3.0, 7.9, 8.5, 9.1, 10.0, 10.5]
        mean, var = v.window_stats(t, [2.0] * len(t))
        self.assertAlmostEqual(mean, 2.0, 14)
        self.assertAlmostEqual(var, 0.0, 14)

    def test_linear_series_and_window_ends(self):
        t = [0.01 * k for k in range(0, 1201)]
        y = [2 * x for x in t]
        mean, var = v.window_stats(t, y)
        self.assertAlmostEqual(mean, 18.0, 12)
        self.assertAlmostEqual(var, 4.0 * 4.0 / 12.0, 3)

    def test_adaptive_steps_are_time_weighted(self):
        # a long step with the value 1 and a short one with the value 5
        t = [8.0, 9.5, 10.0]
        y = [1.0, 1.0, 5.0]
        mean, _ = v.window_stats(t, y)
        self.assertAlmostEqual(mean, (1.5 * 1.0 + 0.5 * 3.0) / 2.0, 12)

    def test_interpolate(self):
        self.assertEqual(v.interpolate([0, 1, 2], [0, 10, 20], 1.5), 15.0)
        self.assertEqual(v.interpolate([0, 1, 2], [0, 10, 20], -1), 0)
        self.assertEqual(v.interpolate([0, 1, 2], [0, 10, 20], 3), 20)

    def test_median(self):
        self.assertEqual(v.median([3, 1, 2]), 2)
        self.assertEqual(v.median([4, 1, 2, 3]), 2.5)


class FieldTests(unittest.TestCase):
    def test_cell_centres(self):
        n, h, c = v.cell_centres(1)
        self.assertEqual(n, 32)
        self.assertAlmostEqual(h, 0.01)
        self.assertEqual(len(c), 1024)
        self.assertAlmostEqual(c[0][0], -0.155)
        self.assertAlmostEqual(c[1][0], -0.145)        # x runs fastest
        self.assertAlmostEqual(c[32][1], -0.145)

    def test_exact_field_has_no_error_and_zero_radial_velocity(self):
        for level in (1, 2):
            f = synthetic_field(level, lambda r: v.exact_utheta(P, r))
            e = v.velocity_errors(P, level, f, f)
            for key in ("ut_L1", "ut_L2", "ut_Linf", "ur_L2", "ur_Linf",
                        "du_5_10"):
                self.assertLess(e[key], 1e-12, key)

    def test_known_perturbation_norms(self):
        level = 2
        scale = P["omega"] * P["R1"]
        base = synthetic_field(level, lambda r: v.exact_utheta(P, r))
        shifted = synthetic_field(level, lambda r: v.exact_utheta(P, r) + 0.01 * scale)
        e = v.velocity_errors(P, level, shifted, base)
        for key in ("ut_L1", "ut_L2", "ut_Linf", "ut_L2_far", "ut_L2_iface"):
            self.assertAlmostEqual(e[key], 0.01, 12, key)
        # a purely radial perturbation changes only the radial norms
        n, h, centres = v.cell_centres(level)
        radial = [(b[0] + 0.02 * scale * x / math.hypot(x, y),
                   b[1] + 0.02 * scale * y / math.hypot(x, y), 0.0)
                  for b, (x, y) in zip(base, centres)]
        e = v.velocity_errors(P, level, radial, base)
        self.assertAlmostEqual(e["ur_L2"], 0.02, 12)
        self.assertLess(e["ut_L2"], 1e-12)

    def test_only_fluid_cells_are_counted(self):
        level = 2
        n, h, centres = v.cell_centres(level)
        expected = sum(1 for x, y in centres
                       if P["R1"] < math.hypot(x, y) < P["R2"])
        f = synthetic_field(level, lambda r: v.exact_utheta(P, r))
        self.assertEqual(v.velocity_errors(P, level, f, f)["fluid_cells"], expected)
        # an error confined to the fictitious interior is not counted
        bad = [(100.0, 0, 0) if math.hypot(x, y) < P["R1"] else u
               for u, (x, y) in zip(f, centres)]
        self.assertLess(v.velocity_errors(P, level, bad, f)["ut_Linf"], 1e-12)

    def test_field_torque_of_the_exact_solution(self):
        for level, tol in ((2, 5e-3), (3, 2e-3)):
            f = synthetic_field(level, lambda r: v.exact_utheta(P, r))
            for frac in (0.3, 0.5, 0.7):
                rs = P["R1"] + (P["R2"] - P["R1"]) * frac
                t = v.field_torque(P, level, f, rs)
                self.assertAlmostEqual(t / v.exact_rotor_torque(P), 1.0,
                                       delta=tol)

    def test_field_torque_of_a_rigid_rotation_is_zero(self):
        f = synthetic_field(2, lambda r: 3.0 * r)
        self.assertLess(abs(v.field_torque(P, 2, f, 0.075)
                            / v.exact_rotor_torque(P)), 1e-3)

    def test_read_vector_field(self):
        with tempfile.TemporaryDirectory() as d:
            path = Path(d) / "U"
            write_field(path, [(1, 2, 3), (4.5, -6e-3, 0)])
            self.assertEqual(v.read_vector_field(path),
                             [(1, 2, 3), (4.5, -6e-3, 0)])
            path.write_text("internalField   uniform (0 0 0);\n")
            with self.assertRaises(SystemExit), \
                    contextlib.redirect_stderr(io.StringIO()):
                v.read_vector_field(path)

    def test_read_force_file_by_header(self):
        with tempfile.TemporaryDirectory() as d:
            path = Path(d) / "x.dat"
            path.write_text(
                "# Immersed body x\n# Time\tforce (x y z)\ttorque (x y z)\t"
                "inertia of the fluid inside (x y z)\t"
                "force from the momentum exchange (x y z)\n"
                "0.5\t1 2 3\t4 5 6\t7 8 9\t10 11 12\n")
            c = v.read_force_file(path)
            self.assertEqual(c["t"], [0.5])
            self.assertEqual(c["torque"], [[4, 5, 6]])
            self.assertEqual(c["inertia"], [[7, 8, 9]])
            self.assertEqual(c["momentum"], [[10, 11, 12]])


class AnalyseCaseTests(unittest.TestCase):
    def build(self, d, scale=1.0, drift=0.0):
        case = Path(d)
        f = synthetic_field(1, lambda r: v.exact_utheta(P, r))
        write_field(case / "10/U", f)
        write_field(case / "5/U", f)
        times = [0.1 * k for k in range(1, 101)]
        exact = v.exact_rotor_torque(P)
        write_forces(case / "postProcessing/immersedBoundary/0/rotor.dat", times,
                     [exact * scale * (1 + drift * t / 10) for t in times])
        write_forces(case / "postProcessing/immersedBoundary/0/stator.dat", times,
                     [-exact * scale for _ in times])
        return case

    def test_exact_flow_is_steady_with_zero_difference(self):
        with tempfile.TemporaryDirectory() as d:
            row, (t, tz) = v.analyse_case(self.build(d), "penalty", "", 1)
        self.assertAlmostEqual(row["T_relerr"], 0.0, 12)
        self.assertEqual(row["steady"], 1)
        self.assertAlmostEqual(row["T_stator_balance"], 0.0, 12)
        self.assertLess(abs(row["T_field_rel"]), 3e-2)
        self.assertEqual(len(t), 100)

    def test_torque_error_sign_is_relative_to_the_exact_torque(self):
        with tempfile.TemporaryDirectory() as d:
            row, _ = v.analyse_case(self.build(d, scale=1.1), "cutLink", "ME", 1)
        self.assertAlmostEqual(row["T_relerr"], 0.10, 12)     # too large
        with tempfile.TemporaryDirectory() as d:
            row, _ = v.analyse_case(self.build(d, scale=0.9), "cutLink", "ME", 1)
        self.assertAlmostEqual(row["T_relerr"], -0.10, 12)    # too small

    def test_a_drifting_torque_is_not_steady(self):
        with tempfile.TemporaryDirectory() as d:
            row, _ = v.analyse_case(self.build(d, drift=0.02), "cutLink", "ME", 1)
        self.assertEqual(row["steady"], 0)
        self.assertGreater(row["T_drift_rel"], v.STEADY_DRIFT)


class CaseGenerationTests(unittest.TestCase):
    def test_fv_options_of_every_method(self):
        tut = (v.TUTORIAL_DIR / "constant/fvOptions").read_text()
        for method, est in v.VARIANTS:
            with tempfile.TemporaryDirectory() as d:
                path = Path(d) / "fvOptions"
                path.write_text(tut)
                v.write_fv_options(path, method, est)
                text = path.read_text()
            self.assertIn("solidBodyRotation", text)
            self.assertEqual(text.count("surface "), 2)
            if method == "penalty":
                self.assertRegex(text, r"method\s+penalty;")
                self.assertEqual(len(re.findall(r"occupancy\s+signedDistance;",
                                                text)), 2)
            if method == "incremental":
                self.assertRegex(text, r"couplingCoeff\s+0.2;")
                self.assertEqual(len(re.findall(r"occupancy\s+vertexFraction;",
                                                text)), 2)
            if method == "ghostCell":
                self.assertRegex(text, r"imageDistance\s+1.5;")
            if method == "cutLink" and est == "ST":
                self.assertRegex(text, r"forceEstimator\s+surfaceTraction;")
            if method == "cutLink" and est == "ME":
                self.assertNotIn("forceEstimator", text)

    def test_generated_case_changes_only_the_method_inputs(self):
        before = {p.name: p.read_bytes() for p in v.TUTORIAL_DIR.rglob("*")
                  if p.is_file() and "verification" not in p.parts}
        with tempfile.TemporaryDirectory() as d:
            case = v.generate_case("incremental", "", 2, "x", 4, cases_dir=Path(d))
            self.assertTrue((case / "Allrun.verify").exists())
            self.assertIn("MESH_LEVEL=2", (case / "Allrun.verify").read_text())
            self.assertIn("(1 4 1)",
                          (case / "system/decomposeParDict").read_text())
            self.assertRegex((case / "system/fvSolution").read_text(),
                             r"nOuterCorrectors\s+5;")
            self.assertEqual((case / "system/blockMeshDict").read_bytes(),
                             (v.TUTORIAL_DIR / "system/blockMeshDict").read_bytes())
            self.assertFalse((case / "verification").exists())
            hashes = v.input_hashes(case)
            self.assertIn("constant/fvOptions", hashes)
        after = {p.name: p.read_bytes() for p in v.TUTORIAL_DIR.rglob("*")
                 if p.is_file() and "verification" not in p.parts}
        self.assertEqual(before, after)

    def test_time_step_cap_is_written_only_when_requested(self):
        with tempfile.TemporaryDirectory() as d:
            plain = v.generate_case("penalty", "", 1, "a", 1, cases_dir=Path(d))
            capped = v.generate_case("penalty", "", 1, "b", 1, cases_dir=Path(d),
                                     dt_cap=0.01)
            self.assertNotIn("maxDeltaT",
                             (plain / "system/controlDict").read_text())
            text = (capped / "system/controlDict").read_text()
            self.assertRegex(text, r"maxCo\s+0.5;\nmaxDeltaT\s+0.01;")
            self.assertRegex(text, r"adjustTimeStep\s+yes;")

    def test_existing_case_is_not_overwritten_without_force(self):
        with tempfile.TemporaryDirectory() as d:
            v.generate_case("penalty", "", 1, "x", 1, cases_dir=Path(d))
            with self.assertRaises(SystemExit), \
                    contextlib.redirect_stderr(io.StringIO()):
                v.generate_case("penalty", "", 1, "x", 1, cases_dir=Path(d))
            v.generate_case("penalty", "", 1, "x", 1, force=True,
                            cases_dir=Path(d))

    def test_plan_selection(self):
        class A:
            runs, methods, levels, ranks = ["matrix"], [], [], 0
        self.assertEqual(len(v.plan(A)), 18)
        A.methods, A.levels = ["penalty"], [2]
        self.assertEqual([x[3] for x in v.plan(A)], ["penalty-L2"])
        A.runs, A.methods, A.levels = ["L4"], [], []
        self.assertEqual(len(v.plan(A)), 5)
        self.assertTrue(all(x[4] == 4 for x in v.plan(A)))
        A.runs = ["dt"]
        self.assertEqual(len(v.plan(A)), 40)
        self.assertIn("penalty-L2-dt0.0025", [x[3] for x in v.plan(A)])
        self.assertTrue(all(x[5] > 0 for x in v.plan(A)))


class DriverBehaviourTests(unittest.TestCase):
    def run_main(self, argv, env=None):
        out = io.StringIO()
        with contextlib.redirect_stdout(out):
            code = v.main(argv)
        return code, out.getvalue()

    def test_dry_run_starts_nothing(self):
        with mock.patch.object(v.subprocess, "call") as call:
            code, out = self.run_main(["dry-run", "--runs", "matrix"])
        self.assertEqual(code, 0)
        call.assert_not_called()
        self.assertIn("-4.188790e-06", out)

    def test_run_without_openfoam_is_skipped(self):
        with tempfile.TemporaryDirectory() as d:
            work = Path(d) / "work"
            with mock.patch.object(v, "WORK_DIR", work), \
                    mock.patch.object(v, "CASES_DIR", work / "cases"), \
                    mock.patch.object(v, "MANIFEST", work / "manifest.json"), \
                    mock.patch.dict(os.environ, {}, clear=True), \
                    mock.patch.object(v.subprocess, "call") as call:
                code, out = self.run_main(["run"])
        self.assertEqual(code, 0)
        call.assert_not_called()
        self.assertIn("SKIPPED cutLink-ME-L1", out)

    def test_check_semantics_count_skipped_apart_from_pass(self):
        with tempfile.TemporaryDirectory() as d:
            work = Path(d) / "work"
            with mock.patch.object(v, "WORK_DIR", work), \
                    mock.patch.object(v, "CASES_DIR", work / "cases"), \
                    mock.patch.object(v, "MANIFEST", work / "manifest.json"):
                code, out = self.run_main(["check"])
        self.assertEqual(code, 0)
        self.assertRegex(out, r"\d+ PASS, 0 FAIL, \d+ SKIPPED")
        self.assertIn("SKIPPED: no completed case", out)

    def test_solver_problem_ignores_the_fpe_trapping_banner(self):
        banner = ("trapFpe: Floating point exception trapping enabled "
                  "(FOAM_SIGFPE).")
        with tempfile.TemporaryDirectory() as d:
            case = Path(d)
            (case / "log.solids4Foam").write_text(
                f"{banner}\nTime = 1\n\nEnd\n")
            self.assertEqual(v.solver_problem(case), "")
            (case / "log.solids4Foam").write_text(
                f"{banner}\nTime = 1\nsigFpe: floating point exception\nEnd\n")
            self.assertIn("NaN", v.solver_problem(case))

    def test_solver_problem_detection(self):
        with tempfile.TemporaryDirectory() as d:
            case = Path(d)
            self.assertEqual(v.solver_problem(case), "no solver log")
            (case / "log.solids4Foam").write_text("Time = 1\nEnd\n")
            self.assertEqual(v.solver_problem(case), "")
            (case / "log.solids4Foam").write_text("Time = 1\n")
            self.assertIn("did not finish", v.solver_problem(case))
            (case / "log.solids4Foam").write_text("nan\nEnd\n")
            self.assertIn("NaN", v.solver_problem(case))


class TrackedResultTests(unittest.TestCase):
    def test_tracked_rows_follow_the_exact_solution_and_the_criteria(self):
        rows = v.read_results()
        if not rows:
            self.skipTest("no tracked results")
        bad = [x for x in v.check_tracked(rows, P) if x[0] == "FAIL"]
        self.assertEqual(bad, [])
        for r in rows:
            self.assertGreaterEqual(int(r["level"]), 1)
            self.assertEqual(int(r["cells"]),
                             v.MESH_CELLS_ACROSS[int(r["level"])] ** 2)

    def test_every_method_and_level_of_the_matrix_is_tracked(self):
        rows = v.read_results()
        if not rows:
            self.skipTest("no tracked results")
        keys = {v.row_key(r) for r in rows}
        for method, est in v.VARIANTS:
            for level in (1, 2, 3):
                self.assertIn((method, est, level, 0.0), keys)

    def test_results_md_contains_the_generated_table(self):
        rows = v.read_results()
        md = v.VERIFICATION_DIR / "results.md"
        if not rows or not md.exists():
            self.skipTest("no tracked results or results.md")
        table = v.master_table(rows)
        self.assertIn(table.strip(), md.read_text())


class ReviewFixTests(unittest.TestCase):
    def test_check_reports_a_failed_case_as_fail(self):
        with tempfile.TemporaryDirectory() as d:
            work = Path(d) / "work"
            work.mkdir()
            (work / "manifest.json").write_text(
                '{"cases": {"x": {"state": "failed"}}}')
            with mock.patch.object(v, "WORK_DIR", work), \
                    mock.patch.object(v, "CASES_DIR", work / "cases"), \
                    mock.patch.object(v, "MANIFEST", work / "manifest.json"):
                out = io.StringIO()
                with contextlib.redirect_stdout(out):
                    code = v.main(["check"])
        self.assertEqual(code, 1)
        self.assertIn("FAIL: x: the run failed", out.getvalue())

    def test_unsteady_cells_are_marked_in_the_time_step_table(self):
        base = {"method": "penalty", "estimator": "", "level": 1, "dt_cap": "",
                "T_relerr": "-0.1", "ut_L2": "0.01", "steady": "0"}
        table = v.time_step_table([base])
        self.assertIn("-10.00 % / 1.00e-02 *", table)

    def test_steady_flag_comparison_is_skipped_near_the_threshold(self):
        row = {"method": "penalty", "estimator": "", "level": 1,
               "T_relerr": "0", "ut_L2": "0", "steady": "1"}
        near = dict(row, steady="0", T_drift_rel="1.2e-3")
        self.assertEqual(v.check_case(row, near)[-1][0], "SKIPPED")
        far = dict(row, steady="0", T_drift_rel="1e-2")
        self.assertEqual(v.check_case(row, far)[-1][0], "FAIL")


class ReviewOfPullRequestTests(unittest.TestCase):
    """The findings of the review of the pull request"""

    def run_main(self, argv):
        out = io.StringIO()
        with contextlib.redirect_stdout(out):
            code = v.main(argv)
        return code, out.getvalue()

    # 1. the documented form ./Allverify --runs L1 L2
    def test_run_is_inserted_before_parsing_when_not_a_command(self):
        self.assertEqual(v.normalise_argv([]), ["run"])
        self.assertEqual(v.normalise_argv(["--runs", "L1", "L2"]),
                         ["run", "--runs", "L1", "L2"])
        self.assertEqual(v.normalise_argv(["--force"]), ["run", "--force"])
        for command in v.COMMANDS:
            self.assertEqual(v.normalise_argv([command, "x"]), [command, "x"])
        self.assertEqual(v.normalise_argv(["-h"]), ["-h"])
        self.assertEqual(v.normalise_argv(["--help"]), ["--help"])

    def test_main_accepts_runs_without_the_command(self):
        seen = {}

        def fake_run(args):
            seen["args"] = args
            return 0
        with mock.patch.object(v, "cmd_run", fake_run):
            self.assertEqual(v.main(["--runs", "L1", "L2"]), 0)
            self.assertEqual(seen["args"].runs, ["L1", "L2"])
            self.assertEqual(v.main(["--runs", "dt-extended", "--end-time",
                                     "20", "--levels", "3"]), 0)
            self.assertEqual(seen["args"].runs, ["dt-extended"])
            self.assertEqual(seen["args"].end_time, 20.0)
            self.assertEqual(seen["args"].levels, [3])
            self.assertEqual(v.main([]), 0)
            self.assertEqual(seen["args"].runs, ["smoke"])

    def test_main_with_the_documented_form_runs_the_cases(self):
        with tempfile.TemporaryDirectory() as d:
            work = Path(d) / "work"
            with mock.patch.object(v, "WORK_DIR", work), \
                    mock.patch.object(v, "CASES_DIR", work / "cases"), \
                    mock.patch.object(v, "MANIFEST", work / "manifest.json"), \
                    mock.patch.dict(os.environ, {}, clear=True):
                code, out = self.run_main(["--runs", "L1", "L2", "--methods",
                                           "penalty"])
        self.assertEqual(code, 0)
        self.assertIn("SKIPPED penalty-L1", out)
        self.assertIn("SKIPPED penalty-L2", out)

    def test_help_and_unknown_commands(self):
        with self.assertRaises(SystemExit) as cm, \
                contextlib.redirect_stdout(io.StringIO()):
            v.main(["--help"])
        self.assertEqual(cm.exception.code, 0)
        with self.assertRaises(SystemExit) as cm, \
                contextlib.redirect_stderr(io.StringIO()):
            v.main(["--bogus"])
        self.assertEqual(cm.exception.code, 2)

    # 3. wording
    def test_documents_do_not_claim_that_run_analyses_or_checks(self):
        readme = (v.VERIFICATION_DIR / "README.md").read_text()
        self.assertNotIn("then check", readme)
        self.assertNotIn("then analyse", v.__doc__)
        self.assertNotIn("0.05 (Levels", readme)

    # 4. observed cost in dry-run
    def test_dry_run_prints_the_observed_cost(self):
        code, out = self.run_main(["dry-run", "--runs", "matrix"])
        self.assertNotIn("observed n/a", out)
        total = float(re.search(r"total of the listed cases: ([0-9.]+)",
                                out).group(1))
        rows = {v.row_key(r): float(r["rank_h"]) for r in v.read_results()}
        expected = sum(rows[(m, e, lv, 0.0)] for lv in (1, 2, 3)
                       for m, e in v.VARIANTS)
        self.assertAlmostEqual(total, expected, places=2)
        self.assertGreater(total, 1.0)

    def test_observed_rank_hours_come_from_the_tracked_files(self):
        obs = v.load_observed_rank_hours()
        self.assertAlmostEqual(obs[("penalty", "", 3, 0.0)], 0.0617)
        self.assertAlmostEqual(obs[("cutLink", "ME", 4, 0.005)], 5.5656)
        self.assertAlmostEqual(obs[("penalty", "", 3, 0.00125)], 1.5108)
        for name in ("dt", "dt-extended"):
            _, rows = v.RUNS[name]
            self.assertTrue(all(tuple(r[:4]) in obs for r in rows), name)

    def test_dry_run_of_a_twenty_second_run_has_no_observed_cost(self):
        code, out = self.run_main(["dry-run", "--runs", "steady20"])
        self.assertEqual(out.count("observed n/a"), 15)

    # 5. reproducibility of the controls
    def test_dt_extended_covers_exactly_the_rows_of_cases_timestep(self):
        path = v.RESULTS_DIR / "controls" / "cases-timestep.csv"
        tracked = {v.row_key(r) for r in v.read_results(path)}
        class A:
            runs, methods, levels, ranks = ["dt-extended"], [], [], 0
        planned = {(m, e, lv, cap) for m, e, lv, _, _, cap, _ in v.plan(A)}
        self.assertEqual(planned, tracked)
        self.assertEqual(len(tracked), 14)
        # no key of the extended run is already in the matrix or in `dt`
        base = {v.row_key(r) for r in v.read_results()}
        self.assertFalse(planned & base)

    def test_dt_and_dt_extended_cover_every_tracked_cap_row(self):
        class A:
            methods, levels, ranks = [], [], 0
        A.runs = ["matrix", "L4", "dt", "dt-extended"]
        planned = {(m, e, lv, cap) for m, e, lv, _, _, cap, _ in v.plan(A)}
        tracked = {v.row_key(r) for r in v.read_results()}
        tracked |= {v.row_key(r) for r in v.read_results(
            v.RESULTS_DIR / "controls" / "cases-timestep.csv")}
        self.assertEqual(tracked - planned, set())

    def test_steady20_covers_exactly_the_rows_of_steady_extended(self):
        import taylor_couette_steady as steady
        path = v.RESULTS_DIR / "controls" / "steady-extended.csv"
        with path.open() as handle:
            tracked = {r["case"] for r in csv.DictReader(handle)}
        class A:
            runs, methods, levels, ranks = ["steady20"], [], [], 0
        labels = {steady.legacy_label(m, e, lv, ranks, cap, end)
                  for m, e, lv, _, ranks, cap, end in v.plan(A)}
        self.assertEqual(labels, tracked)

    def test_end_time_is_written_and_tags_the_name(self):
        class A:
            runs, methods, levels, ranks = ["smoke"], [], [], 0
            end_time = 20.0
        (item,) = v.plan(A)
        self.assertEqual(item[3], "cutLink-ME-L1-t20")
        self.assertEqual(item[6], 20.0)
        with tempfile.TemporaryDirectory() as d:
            plain = v.generate_case("penalty", "", 1, "a", 1, cases_dir=Path(d))
            long = v.generate_case("penalty", "", 1, "b", 1, cases_dir=Path(d),
                                   dt_cap=0.01, end_time=20.0)
            a = (plain / "system/controlDict").read_text()
            b = (long / "system/controlDict").read_text()
        self.assertRegex(a, r"endTime\s+10;")
        self.assertRegex(b, r"endTime\s+20;")
        self.assertEqual(a.replace("endTime         10;", "endTime         20;"),
                         b.replace("maxCo           0.5;\nmaxDeltaT      0.01;",
                                   "maxCo           0.5;"))

    def test_end_time_option_overrides_a_named_run(self):
        class A:
            runs, methods, levels, ranks = ["L1"], ["penalty"], [], 0
            end_time = 20.0
        self.assertEqual([x[3] for x in v.plan(A)], ["penalty-L1-t20"])
        self.assertEqual(v.case_name("penalty", "", 1, 0.0025, 20.0),
                         "penalty-L1-dt0.0025-t20")
        self.assertEqual(v.case_name("penalty", "", 1), "penalty-L1")

    # 6. robustness of the field analysis
    def test_guard_passes_on_the_tutorial(self):
        states = [s for s, _ in v.check_analytic(P)]
        self.assertEqual(states, ["PASS"] * len(states))
        self.assertEqual(v.mesh_guard(P)[0], "PASS")

    def test_guard_fails_when_the_rotation_is_not_centred_at_the_origin(self):
        moved = dict(P, origin=(0.01, 0.0, 0.0), CofR=(0.01, 0.0, 0.0))
        out = v.check_analytic(moved)
        self.assertIn("FAIL", [s for s, m in out if "centred at (0, 0)" in m])

    def test_guard_fails_when_the_mesh_changes(self):
        for change in (dict(mesh_x=(-0.2, 0.2)), dict(mesh_y=(-0.16, 0.17)),
                       dict(mesh_x=(-0.15, 0.17)),
                       dict(cells_formula=(16.0, 2.0, 1.0)),
                       dict(cells_formula=None), dict(square_block=False)):
            self.assertEqual(v.mesh_guard(dict(P, **change))[0], "FAIL", change)

    def test_guard_reads_a_modified_blockmeshdict(self):
        with tempfile.TemporaryDirectory() as d:
            tut = Path(d) / "t"
            import shutil
            shutil.copytree(v.TUTORIAL_DIR, tut,
                            ignore=shutil.ignore_patterns("verification",
                                                          "regressionTests"))
            bm = tut / "system/blockMeshDict"
            text = bm.read_text()
            self.assertEqual(v.mesh_guard(v.parameters(tut))[0], "PASS")
            bm.write_text(text.replace("(-0.16 -0.16 0)", "(-0.2 -0.2 0)"))
            self.assertEqual(v.mesh_guard(v.parameters(tut))[0], "FAIL")
            bm.write_text(text.replace("32*pow(2,", "40*pow(2,"))
            self.assertEqual(v.mesh_guard(v.parameters(tut))[0], "FAIL")

    # 7. OpenFOAM flavour
    def test_openfoam_com_detection(self):
        com = {"WM_PROJECT": "OpenFOAM", "WM_PROJECT_VERSION": "v2412"}
        self.assertTrue(v.is_openfoam_com(com))
        self.assertTrue(v.is_openfoam_com(dict(com, WM_PROJECT_VERSION="v2306")))
        for env in ({"WM_PROJECT": "OpenFOAM", "WM_PROJECT_VERSION": "11"},
                    {"WM_PROJECT": "OpenFOAM", "WM_PROJECT_VERSION": "dev"},
                    {"WM_PROJECT": "OpenFOAM"},
                    {"WM_PROJECT": "foam", "WM_PROJECT_VERSION": "v4.1"},
                    {"WM_PROJECT_VERSION": "v2412"}, {}):
            self.assertFalse(v.is_openfoam_com(env), env)

    def test_check_environment_reasons(self):
        base = {"WM_PROJECT_DIR": "/opt/OpenFOAM-11", "WM_PROJECT": "OpenFOAM",
                "WM_PROJECT_VERSION": "11"}
        with mock.patch.dict(os.environ, base, clear=True):
            self.assertIn("OpenFOAM.com", v.check_environment(1))
        # a path that does not contain 'openfoam.org' is still recognised
        with mock.patch.dict(os.environ, dict(base, WM_PROJECT_DIR="/opt/of"),
                             clear=True):
            self.assertIn("OpenFOAM.com", v.check_environment(1))
        with mock.patch.dict(os.environ, dict(base, WM_PROJECT_VERSION="v2412"),
                             clear=True), \
                mock.patch.object(v.shutil, "which", return_value=None):
            self.assertIn("blockMesh not found", v.check_environment(1))
        with mock.patch.dict(os.environ, {}, clear=True):
            self.assertIn("WM_PROJECT_DIR", v.check_environment(1))

    def test_run_with_openfoam_org_is_skipped_never_failed(self):
        org = {"WM_PROJECT_DIR": "/opt/OpenFOAM-11", "WM_PROJECT": "OpenFOAM",
               "WM_PROJECT_VERSION": "11"}
        with tempfile.TemporaryDirectory() as d:
            work = Path(d) / "work"
            with mock.patch.object(v, "WORK_DIR", work), \
                    mock.patch.object(v, "CASES_DIR", work / "cases"), \
                    mock.patch.object(v, "MANIFEST", work / "manifest.json"), \
                    mock.patch.dict(os.environ, org, clear=True), \
                    mock.patch.object(v.subprocess, "call") as call:
                code, out = self.run_main(["run"])
                code_check, out_check = self.run_main(["check"])
        self.assertEqual(code, 0)
        call.assert_not_called()
        self.assertIn("SKIPPED cutLink-ME-L1: OpenFOAM.com (ESI) is required",
                      out)
        self.assertNotIn("FAIL", out)
        self.assertEqual(code_check, 0)
        self.assertRegex(out_check, r"0 FAIL")


if __name__ == "__main__":
    unittest.main()
