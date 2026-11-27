#!/usr/bin/env python3
"""Exact tests of the generators of the control data, taylor_couette_steady.py
and taylor_couette_st_audit.py (no OpenFOAM, no numpy). The exact-field columns
of the surface-traction audit are compared with the tracked st-audit.csv."""

import contextlib
import io
import math
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import immersed_taylor_couette_verification as v  # noqa: E402
import taylor_couette_controls as tc  # noqa: E402
import taylor_couette_steady as steady  # noqa: E402
import taylor_couette_st_audit as st  # noqa: E402
import test_immersed_taylor_couette_verification as base  # noqa: E402

P = v.parameters()


def write_case(case: Path, level: int, torque, fields: dict) -> None:
    """A synthetic case: rotor.dat of the torque function of time and the
    velocity fields at the stated times"""
    times = [0.01 * (k + 1) for k in range(2000)]
    base.write_forces(case / "postProcessing/immersedBoundary/0/rotor.dat",
                      times, [torque(t) for t in times])
    for tk, func in fields.items():
        base.write_field(case / f"{tk}/U", base.synthetic_field(level, func))


class SteadyGeneratorTests(unittest.TestCase):
    def test_label_of_the_tracked_files(self):
        self.assertEqual(steady.legacy_label("penalty", "", 2, 1, 0.0025),
                         "TC-penalty-L2-np1-dt0025-t20")
        self.assertEqual(steady.legacy_label("incremental", "", 1, 1, 0.005),
                         "TC-incremental-L1-np1-dt005-t20")
        self.assertEqual(steady.legacy_label("incremental", "", 1, 1, 0.0),
                         "TC-incremental-L1-np1-t20")
        self.assertEqual(steady.legacy_label("ghostCell", "", 1, 1, 0.01),
                         "TC-ghostCell-L1-np1-dt01-t20")
        self.assertEqual(steady.legacy_label("cutLink", "ME", 3, 1, 0.0, 20.0),
                         "TC-cutLink-ME-L3-np1-t20")

    def test_tracked_labels_round_trip_through_the_named_run(self):
        rows = tc.read_rows(tc.CONTROLS / "steady-extended.csv")
        self.assertEqual(len(rows), 15)
        self.assertEqual([r["case"] for r in rows],
                         sorted(r["case"] for r in rows))

    def test_exact_steady_flow_has_zero_difference_and_is_steady(self):
        exact = v.exact_rotor_torque(P)

        def ut(r):
            return v.exact_utheta(P, r)
        with tempfile.TemporaryDirectory() as d:
            case = Path(d)
            write_case(case, 1, lambda t: exact, {10: ut, 15: ut, 20: ut})
            row = steady.steady_row(case, "x", 1, P)
        for key in ("T_8_10", "T_13_15", "T_18_20"):
            self.assertAlmostEqual(float(row[key]), 0.0, places=9)
        self.assertEqual(float(row["dU_10_15"]), 0.0)
        self.assertEqual(float(row["dU_15_20"]), 0.0)
        self.assertAlmostEqual(float(row["dLdt_15_20_rel"]), 0.0, places=12)
        self.assertEqual(row["steady_by_criterion"], 1)
        self.assertEqual(float(row["drift_13_15_to_18_20"]), 0.0)

    def test_windows_drift_and_velocity_change(self):
        exact = v.exact_rotor_torque(P)
        scale = P["omega"] * P["R1"]

        def ut(r):
            return v.exact_utheta(P, r)

        def ut_plus(r):                  # +1e-2 omega R1 azimuthal, uniform
            return ut(r) + 1e-2 * scale
        with tempfile.TemporaryDirectory() as d:
            case = Path(d)
            # the torque doubles in magnitude by 1 % per 5 s after 10 s
            write_case(case, 1, lambda t: exact * (1 + 0.01 * max(0, t - 10) / 5),
                       {10: ut, 15: ut, 20: ut_plus})
            row = steady.steady_row(case, "x", 1, P)
        self.assertAlmostEqual(float(row["T_8_10"]), 0.0, places=9)
        self.assertAlmostEqual(float(row["T_13_15"]), 0.01 * 4 / 5, places=3)
        self.assertAlmostEqual(float(row["T_18_20"]), 0.01 * 9 / 5, places=3)
        self.assertEqual(row["steady_by_criterion"], 0)
        self.assertAlmostEqual(float(row["dU_15_20"]), 1e-2, places=6)
        self.assertEqual(float(row["dU_10_15"]), 0.0)
        self.assertGreater(float(row["dLdt_15_20_rel"]), 0.0)

    def test_history_rows_keep_the_written_time(self):
        with tempfile.TemporaryDirectory() as d:
            case = Path(d)
            base.write_forces(case / "postProcessing/immersedBoundary/0/rotor.dat",
                              [0.0011999, 0.5], [-1.0083818e-03, 3e-6])
            rows = steady.history_rows(case, "L")
        self.assertEqual(rows[0], ["L", "0.0011999", "-1.0083818e-03"])
        self.assertEqual(rows[1][2], "3.0000000e-06")

    def test_only_completed_t20_cases_are_taken(self):
        with tempfile.TemporaryDirectory() as d:
            cases = Path(d)
            (cases / "a").mkdir()
            (cases / "b").mkdir()
            manifest = {"cases": {
                "a": {"state": "completed", "end_time": 20.0, "method": "penalty",
                      "estimator": "", "level": 1, "ranks": 1, "dt_cap": 0.01},
                "b": {"state": "completed", "end_time": 0.0, "method": "penalty",
                      "estimator": "", "level": 1, "ranks": 1},
                "c": {"state": "failed", "end_time": 20.0}}}
            found = steady.completed_t20(manifest, cases)
        self.assertEqual([x[0] for x in found], ["TC-penalty-L1-np1-dt01-t20"])


class StAuditTests(unittest.TestCase):
    def test_linear_solver(self):
        a = [[2.0, 1.0, 0.0], [1.0, 3.0, 1.0], [0.0, 1.0, 4.0]]
        x = [[1.0, -1.0], [2.0, 0.5], [-3.0, 2.0]]
        b = [[sum(a[i][k] * x[k][c] for k in range(3)) for c in range(2)]
             for i in range(3)]
        got = st.solve(a, b)
        for i in range(3):
            for c in range(2):
                self.assertAlmostEqual(got[i][c], x[i][c], places=12)
        # pivoting: a zero on the diagonal
        got = st.solve([[0.0, 1.0], [1.0, 0.0]], [[2.0], [3.0]])
        self.assertEqual(got, [[3.0], [2.0]])

    def test_stencil_selects_the_cells_within_three_widths(self):
        for level in (1, 2):
            n, h, centres = v.cell_centres(level)
            for point in ((0.05, 0.0, 0.05), (0.0, -0.0499, 0.05),
                          (0.0353, 0.0353, 0.05)):
                brute = [i for i, (x, y) in enumerate(centres)
                         if (x - point[0]) ** 2 + (y - point[1]) ** 2
                         < (3 * h) ** 2]
                self.assertEqual(st.stencil(level, point), brute)

    def test_rigid_rotation_has_no_viscous_traction(self):
        n, h, centres = v.cell_centres(2)
        rigid = [(-P["omega"] * y, P["omega"] * x) for x, y in centres]
        flat = [0.0] * len(centres)
        pts = [(P["R1"], 0.0, 0.05), (0.0, P["R1"], 0.05)]
        res = st.traction_fit(2, rigid, flat, pts, [(1.0, 0.0), (0.0, 1.0)])
        for tp, tv, pw, ncells in res:
            self.assertAlmostEqual(tv[0], 0.0, places=15)
            self.assertAlmostEqual(tv[1], 0.0, places=15)
            self.assertEqual(pw, 0.0)
            self.assertGreater(ncells, 5)

    def test_traction_of_the_exact_field_is_the_exact_shear(self):
        n, h, _ = v.cell_centres(3)
        ue, pe = st.exact_fields(3)
        pts = [(P["R1"], 0.0, 0.05)]
        (tp, tv, pw, _), = st.traction_fit(3, ue, pe, pts, [(1.0, 0.0)])
        # t_theta on the rotor is the shear -2 mu B / R1^2, y component
        self.assertAlmostEqual(tv[1] / v.exact_shear(P, P["R1"]), 1.0, delta=0.01)
        self.assertAlmostEqual(tv[0], 0.0, delta=1e-3 * abs(tv[1]))

    def test_thinning_keeps_the_area(self):
        tris = st.read_triangles(v.TUTORIAL_DIR / "constant/triSurface/rotor.stl")
        pts, fac, area, h = st.solver_points(tris, 1)
        self.assertEqual(len(pts), 42)
        self.assertAlmostEqual(sum(area), 2 * math.pi * P["R1"] * P["Lz"],
                               delta=1e-4 * 2 * math.pi * P["R1"] * P["Lz"])

    def test_exact_field_columns_reproduce_the_tracked_audit(self):
        tracked = {int(r["level"]): r
                   for r in tc.read_rows(tc.CONTROLS / "st-audit.csv")}
        tris = st.read_triangles(v.TUTORIAL_DIR / "constant/triSurface/rotor.stl")
        for level in (1, 2, 3, 4):
            row = st.audit_level(level, tris)
            ref = tracked[level]
            self.assertEqual(row["thinned_points"], int(ref["thinned_points"]))
            for key in ("recipe_exact_field_circle_wall_pct",
                        "recipe_exact_field_polygon_wall_pct"):
                self.assertAlmostEqual(row[key], float(ref[key]), delta=0.0006,
                                       msg=(level, key))

    def test_numerical_columns_with_an_exact_run(self):
        """A case whose fields are the exact ones and whose rotor torque is
        the exact torque: the numerical recipe equals the exact-field recipe
        and replacing the first layer by the exact field changes nothing"""
        tris = st.read_triangles(v.TUTORIAL_DIR / "constant/triSurface/rotor.stl")
        exact = v.exact_rotor_torque(P)
        with tempfile.TemporaryDirectory() as d:
            case = Path(d)
            write_case(case, 1, lambda t: exact,
                       {10: lambda r: v.exact_utheta(P, r)})
            _, _, centres = v.cell_centres(1)
            # the solver pressure has the gauge p(0.075) = 0, like the exact one
            body = "\n".join("%.12g" % st.p_gauge(math.hypot(x, y))
                             for x, y in centres)
            (case / "10").mkdir(exist_ok=True)
            (case / "10/p").write_text(
                "FoamFile\n{\n}\n\ninternalField   nonuniform List<scalar> "
                f"\n{len(centres)}\n(\n{body}\n)\n;\n")
            row = st.audit_level(1, tris, case)
        self.assertAlmostEqual(row["solver_surface_traction_pct"], 0.0, places=6)
        self.assertAlmostEqual(row["recipe_numerical_field_pct"],
                               row["recipe_exact_field_polygon_wall_pct"],
                               places=3)
        self.assertAlmostEqual(row["first_fluid_layer_replaced_by_exact_pct"],
                               row["recipe_exact_field_polygon_wall_pct"],
                               places=3)

    def test_main_writes_the_exact_columns_without_a_run(self):
        with tempfile.TemporaryDirectory() as d:
            out = Path(d) / "audit.csv"
            with contextlib.redirect_stderr(io.StringIO()):
                self.assertEqual(st.main(["--out", str(out), "--levels", "1",
                                          "--cases-dir", d]), 0)
            lines = out.read_text().splitlines()
        self.assertEqual(lines[0].split(","), st.COLUMNS)
        self.assertTrue(lines[1].startswith("1,42,-0.972,-0.931,,,"))


class CheckFindsControlRows(unittest.TestCase):
    def test_control_rows_are_tracked_for_check(self):
        main = {v.row_key(r) for r in v.read_results()}
        ctl = {v.row_key(r) for r in v.read_results(
            v.RESULTS_DIR / "controls" / "cases-timestep.csv")}
        key = ("penalty", "", 3, 0.00125)
        self.assertNotIn(key, main)
        self.assertIn(key, ctl)


if __name__ == "__main__":
    unittest.main()
