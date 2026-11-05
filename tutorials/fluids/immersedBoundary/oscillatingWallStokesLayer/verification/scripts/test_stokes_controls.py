#!/usr/bin/env python3
"""Exact tests of stokes_controls.py (no OpenFOAM, no numpy)."""

import math
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import oscillating_wall_stokes_layer_verification as v  # noqa: E402
import stokes_controls as c  # noqa: E402

P = v.parameters()


def geometric(limit, first, ratio, n):
    """values limit + e_k with the error e_k halving by 'ratio' per step"""
    return [limit + first / ratio ** k for k in range(n)]


class TimeStepSeries(unittest.TestCase):
    def test_first_order_series_gives_limit_and_bound(self):
        vals = geometric(1.0, 0.8, 2.0, 6)
        a = c.assess(vals)
        self.assertEqual(a["verdict"], "first-order")
        self.assertAlmostEqual(a["ratio"], 2.0, places=9)
        self.assertAlmostEqual(a["limit"], 1.0, places=9)
        self.assertAlmostEqual(a["bound"], abs(vals[-1] - 1.0), places=9)

    def test_second_order_series(self):
        a = c.assess(geometric(0.0, 1.0, 4.0, 5))
        self.assertEqual(a["verdict"], "second-order")
        self.assertAlmostEqual(a["limit"], 0.0, places=9)

    def test_two_agreeing_values_are_not_convergence(self):
        self.assertEqual(c.assess([1.0, 0.5])["verdict"], "insufficient")
        self.assertEqual(c.assess([1.0, 0.5, 0.5])["verdict"], "insufficient")

    def test_sign_change_is_non_monotone(self):
        self.assertEqual(c.assess([2.0, 1.0, 0.5, 0.9, 1.1])["verdict"],
                         "non-monotone")

    def test_ratio_outside_the_thresholds(self):
        # differences 8, 4, 1: ratios 2 and 4, not both in one band
        self.assertEqual(c.assess([0, 8, 12, 13])["verdict"],
                         "ratio outside the thresholds")

    def test_saturated_series(self):
        a = c.assess([10.0, 4.0, 1.000, 1.0005, 1.0005 + 1e-6])
        self.assertEqual(a["verdict"], "saturated")
        self.assertAlmostEqual(a["bound"], 1e-6, places=12)

    def test_halving_series_stops_at_the_first_missing_step(self):
        rows = [{"method": "m", "estimator": "", "level": "1", "dt": str(dt),
                 "q": str(i)} for i, dt in enumerate((0.001, 0.0005, 0.000125))]
        s = c.halving_series(rows, "m", "", 1, "q")
        self.assertEqual([x[0] for x in s], [0.001, 0.0005])

    def test_merged_cases_have_unique_keys_and_keep_the_baseline(self):
        rows = c.merged_cases()
        keys = [c.row_key(r) for r in rows]
        self.assertEqual(len(keys), len(set(keys)))
        base = c.read_rows(c.RESULTS / "cases.csv")
        self.assertGreater(len(rows), len(base))
        self.assertEqual(rows[:len(base)], base)

    def test_tracked_temporal_table_is_what_the_code_computes(self):
        tracked = c.read_rows(c.CONTROLS / "temporal-convergence.csv")
        fresh = c.temporal_table(c.merged_cases())
        self.assertEqual(len(tracked), len(fresh))
        for a, b in zip(tracked, fresh):
            self.assertEqual((a["method"], a["estimator"], int(a["level"]), a["quantity"], a["verdict"]),
                             (b["method"], b["estimator"], b["level"], b["quantity"], b["verdict"]))
            self.assertAlmostEqual(float(a["value_last"]), b["value_last"], delta=abs(b["value_last"]) * 1e-5 + 1e-12)


class WallOffset(unittest.TestCase):
    def test_nearest_face_offsets_of_the_tutorial_levels(self):
        got = [c.wall_offset(P, lv) for lv in (1, 2, 3, 4)]
        for g, e in zip(got, (0.064, 0.128, 0.256, -0.488)):
            self.assertAlmostEqual(g, e, places=3)

    def test_wall_on_a_face_has_no_offset(self):
        q = dict(P)
        q["yw"] = 0.125
        for lv in (1, 2, 3, 4):
            self.assertAlmostEqual(c.wall_offset(q, lv), 0.0, places=9)

    def test_acceleration_regression_recovers_the_coefficient(self):
        t = [3.0 + i * 0.001 for i in range(1001)]
        coeff = 0.37 * 0.0125
        diff = [coeff * c.slab_acceleration(P, x) for x in t]
        got, res = c.acceleration_coefficient(P, t, diff)
        self.assertAlmostEqual(got, coeff, places=12)
        self.assertLess(res, 1e-9)

    def test_acceleration_is_the_derivative_of_the_slab_velocity(self):
        t, d = 3.37, 1e-6
        num = (v.slab_velocity(P, t + d) - v.slab_velocity(P, t - d)) / (2 * d)
        self.assertAlmostEqual(c.slab_acceleration(P, t), num, places=5)

    def test_tracked_offset_audit(self):
        rows = c.read_rows(c.CONTROLS / "offset-audit.csv")
        for r in rows:
            fit, pred = float(r["fitted_ME_over_h"]), float(r["predicted_offset_over_h"])
            if r["case"] == "tutorial" and int(r["level"]) >= 2:
                self.assertAlmostEqual(fit, pred, delta=0.005)
            if r["case"] == "wall on a mesh face":
                self.assertAlmostEqual(fit, 0.0, delta=0.005)
                self.assertLess(abs(float(r["amp_err_ME"])), 1e-3)


class CellAverage(unittest.TestCase):
    modes = None

    @classmethod
    def setUpClass(cls):
        cls.modes = v.cached_modes(P)

    def test_average_of_a_fluid_cell_is_close_to_the_centre_value(self):
        hf = P["height"] - P["yw"]
        h, s = v.cell_heights(P, 3)
        j = next(i for i, x in enumerate(s) if x > 3 * h)
        a = c.cell_average(P, s[j], h, 3.5, self.modes)
        e = v.exact_u_ivp(P, s[j], 3.5, hf, self.modes)
        self.assertLess(abs(a - e), 0.01 * P["U0"])

    def test_straddling_cell_mixes_the_slab_velocity(self):
        h = 0.003125
        s = 0.002 * h               # centre just above the wall
        t = 3.3
        a = c.cell_average(P, s, h, t, self.modes)
        wall = v.slab_velocity(P, t)
        # about half of the cell is slab: the average lies near the slab velocity
        self.assertLess(abs(a - wall), 0.05 * P["U0"])

    def test_gauss_weights_integrate_one(self):
        self.assertAlmostEqual(sum(w for _, w in c.GAUSS), 2.0, places=12)
        self.assertAlmostEqual(sum(w * g * g for g, w in c.GAUSS), 2.0 / 3.0, places=12)

    def test_tracked_cell_average_table(self):
        rows = c.read_rows(c.CONTROLS / "cell-average.csv")
        self.assertTrue(rows)
        for r in rows:
            for k in ("centre", "average", "centre_whole", "average_whole"):
                self.assertGreaterEqual(float(r[k]), 0.0)
        fresh = c.cellavg_table(P)
        self.assertEqual(len(rows), len(fresh))
        for a, b in zip(rows, fresh):
            self.assertAlmostEqual(float(a["centre"]), b["centre"], delta=b["centre"] * 1e-5)


class Claims(unittest.TestCase):
    """The statements of controls.md that are numbers, against the tracked rows"""

    def test_controls_md_is_generated_from_the_tracked_tables(self):
        import make_controls_doc as m
        self.assertEqual(m.DOC.read_text(), m.build())

    def test_baseline_rows_are_not_overridden_by_controls(self):
        base = {c.row_key(r) for r in c.read_rows(c.RESULTS / "cases.csv")}
        extra = [c.row_key(r) for r in c.read_rows(c.CONTROLS / "cases-timestep.csv")]
        self.assertFalse(base & set(extra))
        self.assertEqual(len(extra), len(set(extra)))

    def test_velocity_regions_are_consistent(self):
        rows = c.read_rows(c.CONTROLS / "velocity-regions.csv")
        self.assertTrue(rows)
        for r in rows:
            self.assertEqual(int(r["rows_near"]) + int(r["rows_bulk"]), int(r["rows_all"]))
            self.assertGreaterEqual(float(r["L2_all"]), float(r["L2_bulk"]) * (int(r["rows_bulk"]) / int(r["rows_all"])) ** 0.5 - 1e-12)

    def test_wall_on_face_removes_the_me_offset(self):
        rows = [r for r in c.read_rows(c.CONTROLS / "offset-audit.csv") if r["case"] == "wall on a mesh face"]
        self.assertEqual(sorted(r["level"] for r in rows), ["2", "3", "4"])
        for r in rows:
            self.assertLess(abs(float(r["fitted_ME_over_h"])), 2e-3)
            self.assertLess(abs(float(r["amp_err_ME"])), 3e-4)

    def test_tutorial_offset_matches_prediction(self):
        for r in c.read_rows(c.CONTROLS / "offset-audit.csv"):
            if r["case"] == "tutorial":
                self.assertAlmostEqual(float(r["fitted_ME_over_h"]), float(r["predicted_offset_over_h"]), delta=3e-3)


if __name__ == "__main__":
    unittest.main()
