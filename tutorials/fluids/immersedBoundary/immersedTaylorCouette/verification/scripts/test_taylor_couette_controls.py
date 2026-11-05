#!/usr/bin/env python3
"""Exact tests of taylor_couette_controls.py (no OpenFOAM, no numpy)."""

import math
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import taylor_couette_controls as c  # noqa: E402


def row(method, est, level, cap, err, steady="1"):
    return {"method": method, "estimator": est, "level": str(level),
            "dt_cap": "" if cap == 0 else str(cap), "T_relerr": str(err),
            "steady": steady, "ut_L2": "0"}


class Series(unittest.TestCase):
    def test_cap_series_stops_at_the_first_missing_cap(self):
        rows = [row("m", "", 1, 0.01, 1.0), row("m", "", 1, 0.005, 2.0),
                row("m", "", 1, 0.00125, 3.0)]
        self.assertEqual(c.cap_series(rows, "m", "", 1), [(0.01, 1.0), (0.005, 2.0)])

    def test_tutorial_row_is_key_zero_and_not_a_cap(self):
        rows = [row("m", "", 1, 0, 5.0), row("m", "", 1, 0.01, 1.0)]
        self.assertEqual(c.row_key(rows[0])[3], 0.0)
        self.assertEqual(c.cap_series(rows, "m", "", 1), [(0.01, 1.0)])

    def test_temporal_change_small_threshold(self):
        small = c.temporal_small([(0.01, 1.00), (0.005, 1.05)])
        self.assertTrue(small["small"])
        big = c.temporal_small([(0.01, 1.00), (0.005, 1.50)])
        self.assertFalse(big["small"])
        self.assertAlmostEqual(big["relative"], 0.5 / 1.5)
        self.assertIsNone(c.temporal_small([(0.01, 1.0)])["small"])

    def test_observed_rates(self):
        e = [0.8, 0.2, 0.05, 0.0125]            # h^2
        for r in c.observed_rates(e):
            self.assertAlmostEqual(r, 2.0, places=12)
        self.assertTrue(math.isnan(c.observed_rates([0.0, 1.0])[0]))


class SpatialStatement(unittest.TestCase):
    def rows(self, errs, caps_errs=None, steady="1"):
        rows = []
        for lv, e in enumerate(errs, 1):
            rows.append(row("m", "", lv, 0.005, e, steady))
            rows.append(row("m", "", lv, 0.01, e * (caps_errs or 1.0), steady))
        return rows

    def statement(self, rows):
        out = c.spatial_table(rows)
        return next(r for r in out if r["cap"] == 0.005)["statement"]

    def test_rate_is_stated_for_clean_first_order_data(self):
        self.assertEqual(self.statement(self.rows([0.08, 0.04, 0.02, 0.01], 1.02)),
                         "observed rate 1.00")

    def test_no_rate_when_the_temporal_change_is_not_small(self):
        self.assertEqual(self.statement(self.rows([0.08, 0.04, 0.02, 0.01], 1.5)),
                         "trend only")

    def test_no_rate_when_rates_scatter(self):
        self.assertEqual(self.statement(self.rows([0.08, 0.04, 0.01, 0.0099], 1.02)),
                         "trend only")

    def test_no_rate_for_a_non_monotone_error(self):
        self.assertEqual(self.statement(self.rows([0.08, 0.02, 0.03, 0.01], 1.02)),
                         "trend only")

    def test_all_steady_flag_is_reported_not_hidden(self):
        out = c.spatial_table(self.rows([0.08, 0.04, 0.02, 0.01], 1.02, steady="0"))
        self.assertTrue(all(r["all_steady"] == "no" for r in out))


class Tracked(unittest.TestCase):
    def test_merged_cases_keep_the_baseline_and_have_unique_keys(self):
        rows = c.merged_cases()
        base = c.read_rows(c.RESULTS / "cases.csv")
        keys = [c.row_key(r) for r in rows]
        self.assertEqual(len(keys), len(set(keys)))
        self.assertEqual(rows[:len(base)], base)
        self.assertGreater(len(rows), len(base))

    def test_tracked_tables_are_what_the_code_computes(self):
        rows = c.merged_cases()
        for name, fresh in (("timestep-table.csv", c.timestep_table(rows)),
                            ("spatial-table.csv", c.spatial_table(rows))):
            tracked = c.read_rows(c.CONTROLS / name)
            self.assertEqual(len(tracked), len(fresh), name)
            for a, b in zip(tracked, fresh):
                for k in ("method", "estimator", "level", "cap", "statement", "temporal_change_small"):
                    if k in a and k in b:
                        self.assertEqual(str(a[k]).rstrip("0").rstrip(".") if k == "cap" else a[k],
                                         str(b[k]).rstrip("0").rstrip(".") if k == "cap" else str(b[k]), (name, k))

    def test_control_rows_have_the_requested_cap_and_the_expected_cells(self):
        cells = {1: 1024, 2: 4096, 3: 16384, 4: 65536}
        extra = c.read_rows(c.CONTROLS / "cases-timestep.csv")
        self.assertTrue(extra)
        for r in extra:
            self.assertEqual(int(r["cells"]), cells[int(r["level"])])
            self.assertGreater(float(r["dt_cap"] or 0), 0.0)
            self.assertEqual(r["steps"] != "", True)

    def test_surface_traction_audit_is_consistent(self):
        rows = c.read_rows(c.CONTROLS / "st-audit.csv")
        self.assertEqual([int(r["level"]) for r in rows], [1, 2, 3, 4])
        circle = [float(r["recipe_exact_field_circle_wall_pct"]) for r in rows]
        self.assertTrue(all(abs(a) > abs(b) for a, b in zip(circle, circle[1:])))
        for r in rows:
            # the re-implementation reproduces the solver's written torque
            self.assertAlmostEqual(float(r["recipe_numerical_field_pct"]),
                                   float(r["solver_surface_traction_pct"]), delta=0.1)
        num = [float(r["recipe_numerical_field_pct"]) for r in rows]
        self.assertTrue(all(abs(a) < abs(b) for a, b in zip(num, num[1:])))

    def test_extended_steady_state_flag_follows_the_unrelaxed_criterion(self):
        rows = c.read_rows(c.CONTROLS / "steady-extended.csv")
        self.assertTrue(rows)
        for r in rows:
            expect = (float(r["drift_13_15_to_18_20"]) < 1e-3 and float(r["dU_15_20"]) < 1e-3)
            self.assertEqual(r["steady_by_criterion"], "1" if expect else "0", r["case"])


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

    def test_only_cutlink_me_has_a_stated_spatial_rate(self):
        stated = [r for r in c.read_rows(c.CONTROLS / "spatial-table.csv") if r["statement"].startswith("observed")]
        self.assertTrue(stated)
        for r in stated:
            self.assertEqual((r["method"], r["estimator"]), ("cutLink", "ME"))
            self.assertEqual(r["levels"], "4")

    def test_penalty_fine_meshes_are_not_time_converged(self):
        rows = {(r["method"], int(r["level"])): r for r in c.read_rows(c.CONTROLS / "timestep-table.csv")
                if r["estimator"] == ""}
        self.assertEqual(rows[("penalty", 3)]["temporal_change_small"], "no")
        self.assertEqual(rows[("penalty", 4)]["temporal_change_small"], "no")
        self.assertEqual(rows[("penalty", 2)]["temporal_change_small"], "yes")

    def test_every_control_step_cap_is_finer_than_the_tutorial_step_on_level_4(self):
        for r in c.read_rows(c.CONTROLS / "cases-timestep.csv"):
            if r["level"] == "4":
                self.assertLess(float(r["dt_cap"]), 0.0126 + 1e-9)


if __name__ == "__main__":
    unittest.main()
