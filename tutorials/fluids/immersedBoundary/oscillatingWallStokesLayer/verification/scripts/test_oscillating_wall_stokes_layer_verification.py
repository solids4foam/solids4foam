#!/usr/bin/env python3
"""Exact tests of oscillating_wall_stokes_layer_verification.py: they need
neither OpenFOAM nor numpy. Run with ./Allverify unit."""

import contextlib
import io
import json
import math
import os
import re
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

sys.path.insert(0, str(Path(__file__).resolve().parent))
import oscillating_wall_stokes_layer_verification as v  # noqa: E402

P = v.parameters()
H = P["height"] - P["yw"]
LEVEL = 1


def write_field(path, field):
    path.parent.mkdir(parents=True, exist_ok=True)
    body = "\n".join("(%.12g %.12g %.12g)" % u for u in field)
    path.write_text("FoamFile\n{\n}\n\ninternalField   nonuniform List<vector> "
                    f"\n{len(field)}\n(\n{body}\n)\n;\n\nboundaryField\n{{\n}}\n")


def write_forces(path, times, force, inertia=None, traction=None):
    path.parent.mkdir(parents=True, exist_ok=True)
    head = "# Time\tforce (x y z)\ttorque (x y z)\tinertia of the fluid inside (x y z)"
    if traction is not None:
        head += "\tforce from the surface traction (x y z)"
    lines = ["# Immersed body slab", head]
    for i, t in enumerate(times):
        line = f"{t}\t{force[i]} 0 0\t0 0 0\t{(inertia or force)[i] * 0 if inertia is None else inertia[i]} 0 0"
        if traction is not None:
            line += f"\t{traction[i]} 0 0"
        lines.append(line)
    path.write_text("\n".join(lines) + "\n")


def field_from(level, func, t):
    ny, nx, _, _ = v.mesh_shape(P, level)
    h = P["height"] / ny
    out = []
    for j in range(ny):
        s = (j + 0.5) * h - P["yw"]
        u = func(s, t) if s > 0 else 0.0
        out += [(u, 0.0, 0.0)] * nx
    return out


def periodic_profiles(level, shift=0.0):
    return {tk: v.row_profile(P, level, field_from(
        level, lambda s, t: v.exact_u(P, s, t) + shift, tk))
        for tk in v.snapshot_times()}


class ParameterTests(unittest.TestCase):
    def test_parameters_are_read_from_the_tutorial(self):
        self.assertEqual(P["nu"], 0.00785398)
        self.assertEqual(P["rho"], 1.0)
        self.assertAlmostEqual(P["U0"], 0.1, places=8)
        self.assertAlmostEqual(P["omega"], 2 * math.pi, 12)
        self.assertAlmostEqual(P["delta"], 0.05, places=7)
        self.assertAlmostEqual(P["yw"], 0.1234, 9)
        self.assertEqual(P["height"], 0.5)
        self.assertAlmostEqual(P["depth"], 0.1, 12)
        self.assertEqual(P["direction"], (1.0, 0.0, 0.0))
        self.assertEqual(P["dt"], 0.001)

    def test_mesh_rule_and_wall_area(self):
        for level, ny in ((1, 20), (2, 40), (3, 80), (4, 160)):
            n, nx, lx, area = v.mesh_shape(P, level)
            self.assertEqual(n, ny)
            self.assertEqual(nx, ny // 10)
            self.assertAlmostEqual(lx, 0.05, 12)
            self.assertAlmostEqual(area, 0.005, 12)

    def test_slab_velocity_and_phase_convention(self):
        # position A sin(omega t + phase), phase -pi/2: velocity U0 sin(omega t)
        self.assertAlmostEqual(P["phase"], -math.pi / 2, places=9)
        for t in (0.0, 0.1, 0.25, 0.6, 0.9):
            self.assertAlmostEqual(v.slab_velocity(P, t),
                                   P["U0"] * math.sin(P["omega"] * t), 8)
        h = 1e-6
        pos = lambda t: P["amplitude"] * math.sin(P["omega"] * t + P["phase"])
        self.assertAlmostEqual((pos(0.3 + h) - pos(0.3 - h)) / (2 * h),
                               v.slab_velocity(P, 0.3), 6)

    def test_stl_top_and_dict_scalar(self):
        self.assertAlmostEqual(v.stl_top(
            v.TUTORIAL_DIR / "constant/triSurface/slab.stl"), 0.1234, 9)
        text = "nu [ 0 2 -1 0 0 0 0 ] 0.5;\nrho rho [1 -3 0 0 0 0 0] 2;"
        self.assertEqual(v.dict_scalar(text, "nu"), 0.5)
        self.assertEqual(v.dict_scalar(text, "rho"), 2.0)
        with self.assertRaises(SystemExit), contextlib.redirect_stderr(io.StringIO()):
            v.dict_scalar(text, "mu")


class ExactSolutionTests(unittest.TestCase):
    def test_boundary_values(self):
        for t in (0.0, 0.2, 0.9):
            self.assertAlmostEqual(v.exact_u(P, 0.0, t), v.slab_velocity(P, t), 8)
            self.assertAlmostEqual(v.exact_u_finite(P, 0.0, t, H),
                                   v.slab_velocity(P, t), 8)
        self.assertLess(abs(v.exact_u(P, 20 * P["delta"], 0.3)), 1e-8)

    def test_satisfies_the_diffusion_equation(self):
        h, dt = 1e-4, 1e-5
        for s in (0.01, 0.04, 0.09):
            for t in (0.1, 0.45):
                ut = (v.exact_u(P, s, t + dt) - v.exact_u(P, s, t - dt)) / (2 * dt)
                uss = (v.exact_u(P, s + h, t) - 2 * v.exact_u(P, s, t)
                       + v.exact_u(P, s - h, t)) / h ** 2
                self.assertAlmostEqual(ut, P["nu"] * uss, delta=2e-4)

    def test_finite_height_solution_equations(self):
        h, dt = 1e-4, 1e-5
        for s in (0.02, 0.2, 0.35):
            t = 0.3
            ut = (v.exact_u_finite(P, s, t + dt, H)
                  - v.exact_u_finite(P, s, t - dt, H)) / (2 * dt)
            uss = (v.exact_u_finite(P, s + h, t, H)
                   - 2 * v.exact_u_finite(P, s, t, H)
                   + v.exact_u_finite(P, s - h, t, H)) / h ** 2
            self.assertAlmostEqual(ut, P["nu"] * uss, delta=2e-4)
        # slip lid: du/ds = 0 at s = H
        d = (v.exact_u_finite(P, H, 0.3, H) - v.exact_u_finite(P, H - 1e-6, 0.3, H)) / 1e-6
        self.assertAlmostEqual(d, 0.0, delta=1e-4)

    def test_tau_is_rho_nu_du_ds_at_the_wall(self):
        h = 1e-6
        for t in (0.0, 0.3, 0.77):
            d = (v.exact_u(P, h, t) - v.exact_u(P, -h, t)) / (2 * h)
            self.assertAlmostEqual(P["rho"] * P["nu"] * d, v.exact_tau(P, t), 9)

    def test_units_and_density_factor(self):
        # tau is a stress: rho [kg/m3] nu [m2/s] U0/delta [1/s] = Pa; the
        # tutorial has rho = 1, so the value is also the kinematic one
        q = dict(P)
        q["rho"] = 3.0
        self.assertAlmostEqual(v.exact_tau(q, 0.2) / v.exact_tau(P, 0.2), 3.0, 12)
        self.assertAlmostEqual(P["nu"] * P["U0"] / P["delta"], 0.0157079633, 7)

    def test_amplitude_and_phase(self):
        a, phi = v.exact_amp_phase(P)
        self.assertAlmostEqual(a, P["rho"] * P["nu"] / P["delta"] * P["U0"]
                               * math.sqrt(2), 12)
        for t in (0.0, 0.13, 0.5, 0.81):
            self.assertAlmostEqual(a * math.sin(P["omega"] * t + phi),
                                   v.exact_tau(P, t), 10)

    def test_finite_height_effect_is_negligible(self):
        worst = max(abs(v.exact_tau_finite(P, i / 100.0, H) - v.exact_tau(P, i / 100.0))
                    for i in range(100))
        self.assertLess(worst / v.exact_amp_phase(P)[0], 1e-5)
        self.assertLess(max(abs(v.exact_u_finite(P, s, 0.4, H)
                                - v.exact_u(P, s, 0.4))
                            for s in (0.01, 0.1, 0.3)) / P["U0"], 1e-3)


class InitialValueTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.modes = v.cached_modes(P)

    def test_initial_condition_is_zero(self):
        for s in (0.01, 0.05, 0.1, 0.2, 0.3):
            self.assertAlmostEqual(v.exact_u_ivp(P, s, 0.0, H, self.modes), 0.0,
                                   delta=2e-4 * P["U0"])

    def test_is_periodic_plus_decaying_modes(self):
        for s in (0.02, 0.15):
            for t in (1.0, 3.5):
                tr = sum(a * math.exp(-P["nu"] * k * k * t) * math.sin(k * s)
                         for k, a in self.modes)
                self.assertAlmostEqual(
                    v.exact_u_ivp(P, s, t, H, self.modes),
                    v.exact_u_finite(P, s, t, H) + tr, 14)
        # the transient decays
        a1 = abs(v.exact_u_ivp(P, 0.1, 8.0, H, self.modes)
                 - v.exact_u_finite(P, 0.1, 8.0, H))
        a0 = abs(v.exact_u_ivp(P, 0.1, 1.5, H, self.modes)
                 - v.exact_u_finite(P, 0.1, 1.5, H))
        self.assertLess(a1, a0)

    def test_satisfies_the_diffusion_equation_and_boundaries(self):
        h, dt = 1e-4, 1e-5
        s, t = 0.07, 2.0
        f = lambda s, t: v.exact_u_ivp(P, s, t, H, self.modes)
        ut = (f(s, t + dt) - f(s, t - dt)) / (2 * dt)
        uss = (f(s + h, t) - 2 * f(s, t) + f(s - h, t)) / h ** 2
        self.assertAlmostEqual(ut, P["nu"] * uss, delta=5e-4)
        self.assertAlmostEqual(f(0.0, 2.3), v.slab_velocity(P, 2.3), 8)

    def test_transient_in_the_fourth_period_is_quantified(self):
        amp = v.exact_amp_phase(P)[0]
        ts = [3.0 + 0.01 * k for k in range(1, 101)]
        dtau = [v.exact_tau_ivp(P, t, H, self.modes) - v.exact_tau(P, t) for t in ts]
        rel = v.rms(dtau) / v.rms([v.exact_tau(P, t) for t in ts])
        self.assertGreater(rel, 0.003)
        self.assertLess(rel, 0.006)
        du = [v.exact_u_ivp(P, s, 3.5, H, self.modes) - v.exact_u(P, s, 3.5)
              for s in [0.01 * k for k in range(1, 37)]]
        self.assertGreater(max(abs(x) for x in du) / P["U0"], 0.005)
        self.assertLess(max(abs(x) for x in du) / P["U0"], 0.03)
        self.assertGreater(amp, 0)


class StatisticsTests(unittest.TestCase):
    def test_harmonic_fit_recovers_coefficients(self):
        w = P["omega"]
        t = [3.0 + 0.001 * k for k in range(1, 1001)]
        y = [0.5 + 2 * math.sin(w * x) - 3 * math.cos(w * x)
             + 0.25 * math.sin(2 * w * x) + 0.1 * math.cos(3 * w * x) for x in t]
        c, fit = v.harmonic_fit(t, y, w, 3)
        for got, want in zip(c, (0.5, 2, -3, 0.25, 0, 0, 0.1)):
            self.assertAlmostEqual(got, want, 9)
        self.assertLess(max(abs(a - b) for a, b in zip(fit, y)), 1e-9)

    def test_solve_linear(self):
        x = v.solve_linear([[2.0, 1.0], [1.0, 3.0]], [5.0, 10.0])
        self.assertAlmostEqual(x[0], 1.0)
        self.assertAlmostEqual(x[1], 3.0)

    def test_noise_definition_and_shear_metrics_of_exact_signal(self):
        t = [0.001 * k for k in range(1, 4001)]
        tau = [v.exact_tau(P, x) for x in t]
        m = v.shear_metrics(P, t, tau)
        self.assertLess(abs(m["shear_amp_err"]), 1e-12)
        self.assertLess(abs(m["shear_phase_err_deg"]), 1e-9)
        self.assertLess(m["shear_rms_err"], 1e-12)
        self.assertLess(m["shear_noise"], 1e-12)
        self.assertAlmostEqual(m["shear_amp"], m["shear_amp_exact"], 10)
        # the periodic signal differs from the initial-value target
        self.assertGreater(m["shear_rms_err_ivp"], 0.003)

    def test_known_amplitude_phase_and_noise(self):
        t = [0.001 * k for k in range(1, 4001)]
        a, phi = v.exact_amp_phase(P)
        dphi = math.radians(2.0)
        tau = [1.05 * a * math.sin(P["omega"] * x + phi + dphi)
               + 0.001 * math.sin(20 * P["omega"] * x) for x in t]
        m = v.shear_metrics(P, t, tau)
        self.assertAlmostEqual(m["shear_amp_err"], 0.05, 4)
        self.assertAlmostEqual(m["shear_phase_err_deg"], 2.0, 3)
        rms_ref = a / math.sqrt(2)
        self.assertAlmostEqual(m["shear_noise"], 0.001 / math.sqrt(2) / rms_ref, 4)

    def test_phase_error_sign_is_a_lead(self):
        t = [0.001 * k for k in range(1, 4001)]
        a, phi = v.exact_amp_phase(P)
        tau = [a * math.sin(P["omega"] * (x + 0.01) + phi) for x in t]
        self.assertGreater(v.shear_metrics(P, t, tau)["shear_phase_err_deg"], 3.0)

    def test_mean_offset_is_not_an_amplitude_error(self):
        t = [0.001 * k for k in range(1, 4001)]
        tau = [v.exact_tau(P, x) + 0.003 for x in t]
        m = v.shear_metrics(P, t, tau)
        self.assertLess(abs(m["shear_amp_err"]), 1e-10)
        self.assertAlmostEqual(m["shear_mean"], 0.003, 9)


class EstimatorTests(unittest.TestCase):
    def columns(self):
        with tempfile.TemporaryDirectory() as d:
            path = Path(d) / "slab.dat"
            path.write_text(
                "# Immersed body slab\n# Time\tforce (x y z)\ttorque (x y z)\t"
                "inertia of the fluid inside (x y z)\t"
                "force from the momentum exchange (x y z)\t"
                "force from the surface traction (x y z)\n"
                "0.5\t1 0 0\t0 0 0\t2 0 0\t4 0 0\t8 0 0\n"
                "1.0\t3 0 0\t0 0 0\t5 0 0\t6 0 0\t9 0 0\n")
            return v.read_force_file(path)

    def test_columns_by_header(self):
        c = self.columns()
        self.assertEqual(c["force"], [[1, 0, 0], [3, 0, 0]])
        self.assertEqual(c["inertia"], [[2, 0, 0], [5, 0, 0]])
        self.assertEqual(c["momentum"], [[4, 0, 0], [6, 0, 0]])
        self.assertEqual(c["traction"], [[8, 0, 0], [9, 0, 0]])

    def test_selection(self):
        c = self.columns()
        a = 0.5
        self.assertEqual(v.shear_series(c, "cutLink", "ME", a)[1], [2.0, 6.0])
        self.assertEqual(v.shear_series(c, "cutLink", "FO", a)[1], [2.0, 6.0])
        self.assertEqual(v.shear_series(c, "cutLink", "ST", a)[1], [16.0, 18.0])
        # forcing plus the inertia of the fluid inside the slab
        self.assertEqual(v.shear_series(c, "penalty", "", a)[1], [6.0, 16.0])
        self.assertEqual(v.shear_series(c, "ghostCell", "", a)[1], [6.0, 16.0])
        self.assertEqual(v.shear_series(c, "incremental", "", a)[0], [0.5, 1.0])


class FieldTests(unittest.TestCase):
    def test_cell_centre_sampling(self):
        h, s = v.cell_heights(P, 1)
        self.assertAlmostEqual(h, 0.025)
        self.assertAlmostEqual(s[0], 0.0125 - 0.1234)
        self.assertEqual(sum(1 for x in s if x > 0), 15)    # centres above y_w

    def test_exact_field_has_no_periodic_error(self):
        n = v.velocity_norms(P, LEVEL, periodic_profiles(LEVEL))
        for k in ("u_L1_all", "u_L2_all", "u_Linf_all", "u_L2_near",
                  "u_Linf_iface"):
            self.assertLess(n[k], 1e-12, k)
        # the same field differs from the initial-value solution
        self.assertGreater(n["u_L2_all_ivp"], 1e-3)

    def test_known_perturbation_norms(self):
        n = v.velocity_norms(P, LEVEL, periodic_profiles(LEVEL, 0.02 * P["U0"]))
        for k in ("u_L1_all", "u_L2_all", "u_Linf_all", "u_L1_near",
                  "u_L2_near", "u_Linf_iface"):
            self.assertAlmostEqual(n[k], 0.02, 12, k)

    def test_only_fluid_cells_are_counted(self):
        prof = periodic_profiles(LEVEL)
        for tk in prof:
            prof[tk][0] = 123.0      # a cell inside the slab
        n = v.velocity_norms(P, LEVEL, prof)
        self.assertLess(n["u_Linf_all"], 1e-12)

    def test_subsets(self):
        h, s = v.cell_heights(P, 2)
        prof = {tk: [0.0] * len(s) for tk in v.snapshot_times()}
        exact = v.velocity_norms(P, 2, {tk: [
            v.exact_u(P, x, tk) if x > 0 else 0.0 for x in s]
            for tk in prof})
        self.assertLess(exact["u_Linf_all"], 1e-12)
        # an error in the first fluid cell only: iface and all see it, the
        # near set contains it, the Linf is the same
        first = next(j for j, x in enumerate(s) if x > 0)
        bad = {tk: [v.exact_u(P, x, tk) if x > 0 else 0.0 for x in s]
               for tk in prof}
        for tk in bad:
            bad[tk][first] += 0.01 * P["U0"]
        n = v.velocity_norms(P, 2, bad)
        self.assertAlmostEqual(n["u_Linf_all"], 0.01, 12)
        self.assertAlmostEqual(n["u_Linf_iface"], 0.01, 12)
        self.assertAlmostEqual(n["u_Linf_near"], 0.01, 12)
        self.assertLess(n["u_L2_all"], 0.01)

    def test_row_profile_averages_the_cyclic_direction(self):
        ny, nx, _, _ = v.mesh_shape(P, LEVEL)
        f = [(float(j * 10 + i), 0.0, 0.0) for j in range(ny) for i in range(nx)]
        pr = v.row_profile(P, LEVEL, f)
        self.assertEqual(pr[3], 30.5)
        with self.assertRaises(SystemExit), contextlib.redirect_stderr(io.StringIO()):
            v.row_profile(P, LEVEL, f[:-1])

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


class AnalyseCaseTests(unittest.TestCase):
    def build(self, d, scale=1.0):
        case = Path(d)
        for tk in v.snapshot_times():
            write_field(case / f"{tk:g}" / "U",
                        field_from(LEVEL, lambda s, t: v.exact_u(P, s, t), tk))
        _, _, _, area = v.mesh_shape(P, LEVEL)
        t = [0.001 * k for k in range(1, 4001)]
        force = [scale * v.exact_tau(P, x) * area for x in t]
        write_forces(case / "postProcessing/immersedBoundary/0/slab.dat", t,
                     force, traction=[f * 1.5 for f in force])
        (case / "system").mkdir()
        (case / "system/controlDict").write_text("deltaT 0.001;\n")
        return case

    def test_cutlink_case_gives_me_and_st_rows(self):
        with tempfile.TemporaryDirectory() as d:
            rows, series = v.analyse_case(self.build(d), "cutLink", "ME", LEVEL)
        self.assertEqual([r["estimator"] for r in rows], ["ME", "ST"])
        self.assertLess(abs(rows[0]["shear_amp_err"]), 1e-9)
        self.assertAlmostEqual(rows[1]["shear_amp_err"], 0.5, 6)
        self.assertEqual(rows[0]["cells"], 40)
        self.assertEqual(rows[0]["dt"], 0.001)
        self.assertEqual(set(series), {"ME", "ST"})

    def test_forcing_run_gives_one_row(self):
        with tempfile.TemporaryDirectory() as d:
            rows, _ = v.analyse_case(self.build(d), "cutLink", "FO", LEVEL)
        self.assertEqual([r["estimator"] for r in rows], ["FO"])

    def test_other_methods_give_one_row_without_estimator(self):
        with tempfile.TemporaryDirectory() as d:
            case = self.build(d)
            t = [0.001 * k for k in range(1, 4001)]
            _, _, _, area = v.mesh_shape(P, LEVEL)
            force = [0.9 * v.exact_tau(P, x) * area for x in t]
            write_forces(case / "postProcessing/immersedBoundary/0/slab.dat",
                         t, force, inertia=[0.0] * len(t))
            rows, _ = v.analyse_case(case, "penalty", "", LEVEL)
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["estimator"], "")
        self.assertAlmostEqual(rows[0]["shear_amp_err"], -0.10, 9)

    def test_results_csv_round_trip(self):
        with tempfile.TemporaryDirectory() as d:
            rows, _ = v.analyse_case(self.build(Path(d) / "c"), "cutLink", "ME", LEVEL)
            path = Path(d) / "r.csv"
            v.write_results(path, rows)
            back = v.read_results(path)
        self.assertEqual(len(back), 2)
        self.assertEqual(set(v.FIELDS) - {"rank_h"}, set(back[0]) - {"rank_h"})
        self.assertNotIn("\r", path.name)


class CaseGenerationTests(unittest.TestCase):
    def test_fv_options_of_every_method(self):
        tut = (v.TUTORIAL_DIR / "constant/fvOptions").read_text()
        for method, est in v.VARIANTS:
            with tempfile.TemporaryDirectory() as d:
                path = Path(d) / "fvOptions"
                path.write_text(tut)
                v.write_fv_options(path, method, est)
                text = path.read_text()
            self.assertIn("sinusoidalTranslation", text)
            self.assertEqual(len(re.findall(r'surface\s+"slab.stl";', text)), 1)
            if method == "cutLink":
                self.assertRegex(text, r"writeSurfaceTraction\s+yes;")
            else:
                self.assertNotIn("writeSurfaceTraction", text)
            if method == "penalty":
                self.assertRegex(text, r"method\s+penalty;")
                self.assertRegex(text, r"occupancy\s+signedDistance;")
            if method == "incremental":
                self.assertRegex(text, r"couplingCoeff\s+0.2;")
                self.assertRegex(text, r"occupancy\s+vertexFraction;")
            if method == "ghostCell":
                self.assertRegex(text, r"imageDistance\s+1.5;")
            if method == "cutLink" and est == "FO":
                self.assertRegex(text, r"forceEstimator\s+forcing;")
            if method == "cutLink" and est == "ME":
                self.assertNotIn("forceEstimator", text)

    def test_generated_case_changes_only_the_method_inputs(self):
        before = {p.name: p.read_bytes() for p in v.TUTORIAL_DIR.rglob("*")
                  if p.is_file() and "verification" not in p.parts}
        with tempfile.TemporaryDirectory() as d:
            case = v.generate_case("incremental", "", 2, "x", 4, cases_dir=Path(d))
            self.assertIn("MESH_LEVEL=2", (case / "Allrun.verify").read_text())
            self.assertIn("(1 4 1)",
                          (case / "system/decomposeParDict").read_text())
            self.assertRegex((case / "system/fvSolution").read_text(),
                             r"nOuterCorrectors\s+5;")
            self.assertEqual((case / "system/blockMeshDict").read_bytes(),
                             (v.TUTORIAL_DIR / "system/blockMeshDict").read_bytes())
            self.assertEqual((case / "system/controlDict").read_bytes(),
                             (v.TUTORIAL_DIR / "system/controlDict").read_bytes())
            self.assertFalse((case / "verification").exists())
            self.assertIn("constant/fvOptions", v.input_hashes(case))
        after = {p.name: p.read_bytes() for p in v.TUTORIAL_DIR.rglob("*")
                 if p.is_file() and "verification" not in p.parts}
        self.assertEqual(before, after)

    def test_serial_case_has_no_mpi(self):
        with tempfile.TemporaryDirectory() as d:
            case = v.generate_case("cutLink", "ME", 1, "x", 1, cases_dir=Path(d))
            self.assertNotIn("runParallel", (case / "Allrun.verify").read_text())
            self.assertNotIn("decomposePar", (case / "Allrun.verify").read_text())

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
        self.assertEqual(len(v.plan(A)), 20)
        A.methods, A.levels = ["penalty"], [2]
        self.assertEqual([x[3] for x in v.plan(A)], ["penalty-L2"])
        A.runs, A.methods, A.levels = ["smoke"], [], []
        self.assertEqual([x[3] for x in v.plan(A)], ["cutLink-ME-L1"])


class TimeStepControlTests(unittest.TestCase):
    class A:
        runs, methods, levels, ranks = ["dt"], [], [], 0

    def test_tutorial_step_constant_matches_the_tutorial(self):
        self.assertEqual(v.TUTORIAL_DT, P["dt"])

    def test_plan_count_and_names(self):
        plan = v.plan(self.A)
        self.assertEqual(len(plan), 36)
        names = [x[3] for x in plan]
        self.assertEqual(len(set(names)), 36)
        self.assertIn("penalty-L4-dt0.00025", names)
        self.assertIn("cutLink-ME-L2-dt0.0025", names)
        self.assertTrue(all(x[5] in v.DT_STEPS for x in plan))
        self.assertEqual({x[2] for x in plan}, {2, 3, 4})
        for dt in v.DT_STEPS:     # on the phase grid of the write interval
            self.assertAlmostEqual(0.125 / dt, round(0.125 / dt), 9)

    def test_tutorial_cases_have_no_tag(self):
        self.assertEqual(v.case_name("penalty", "", 2), "penalty-L2")
        self.assertEqual(v.case_name("penalty", "", 2, 0.0005), "penalty-L2-dt0.0005")

    def test_generated_control_dict(self):
        tut = (v.TUTORIAL_DIR / "system/controlDict").read_text()
        with tempfile.TemporaryDirectory() as d:
            case = v.generate_case("penalty", "", 4, "x", 1,
                                   cases_dir=Path(d), dt=0.00025)
            text = (case / "system/controlDict").read_text()
            base = v.generate_case("penalty", "", 4, "y", 1, cases_dir=Path(d))
            self.assertEqual((base / "system/controlDict").read_text(), tut)
        self.assertEqual(v.dict_scalar(text, "deltaT"), 0.00025)
        self.assertRegex(text, r"adjustTimeStep\s+no;")
        self.assertEqual(v.dict_scalar(text, "writeInterval"),
                         v.dict_scalar(tut, "writeInterval"))
        self.assertEqual(v.dict_scalar(text, "endTime"),
                         v.dict_scalar(tut, "endTime"))
        self.assertEqual(text.replace("0.00025", "0.001"), tut)

    def test_row_key_separates_the_time_steps(self):
        a = {"method": "penalty", "estimator": "", "level": 2, "dt": "0.001"}
        b = dict(a, dt="0.0005")
        c = dict(a, dt="")
        self.assertNotEqual(v.row_key(a), v.row_key(b))
        self.assertEqual(v.row_key(a), v.row_key(c))     # blank is the tutorial

    def test_master_table_lists_the_tutorial_step_only(self):
        base = {"method": "penalty", "estimator": "", "cells": 160, "ranks": 1,
                "shear_amp_exact": 0.02, "shear_amp": 0.02,
                "shear_amp_err_ivp": 0.0, "shear_phase_err_deg_ivp": 0.0,
                "shear_rms_err_ivp": 0.0, "u_L1_all_ivp": 0.0,
                "u_L2_all_ivp": 0.0, "u_Linf_all_ivp": 0.0,
                "shear_noise": 0.0, "rank_h": "",
                "shear_amp_err": 0.01, "shear_phase_err_deg": 0.5,
                "shear_rms_err": 0.02, "u_L1_all": 0.0, "u_L2_all": 0.0,
                "u_Linf_all": 0.0}
        rows = [dict(base, level=2, dt="0.001"), dict(base, level=2, dt="")
                .copy(), dict(base, level=3, dt="0.0005")]
        table = v.master_table(rows)
        self.assertEqual(table.count("| penalty |"), 2)
        self.assertNotIn("| penalty | 3 |", table)
        # the two targets give different columns of the same rows
        self.assertIn("+1.00 %", table)
        self.assertIn("+0.00 %", v.master_table(rows, "ivp"))

    def test_results_round_trip_keeps_both_steps(self):
        rows = [{"method": "penalty", "estimator": "", "level": 2, "dt": dt}
                for dt in (0.001, 0.0005)]
        with tempfile.TemporaryDirectory() as d:
            path = Path(d) / "r.csv"
            v.write_results(path, rows)
            back = v.read_results(path)
        self.assertEqual(len({v.row_key(r) for r in back}), 2)
        self.assertEqual(v.row_dt(back[0]), v.TUTORIAL_DT)    # tutorial first


class DriverBehaviourTests(unittest.TestCase):
    def run_main(self, argv):
        out = io.StringIO()
        with contextlib.redirect_stdout(out):
            code = v.main(argv)
        return code, out.getvalue()

    def isolated(self, d):
        work = Path(d) / "work"
        return (mock.patch.object(v, "WORK_DIR", work),
                mock.patch.object(v, "CASES_DIR", work / "cases"),
                mock.patch.object(v, "MANIFEST", work / "manifest.json"))

    def test_dry_run_starts_nothing(self):
        with mock.patch.object(v.subprocess, "call") as call:
            code, out = self.run_main(["dry-run", "--runs", "matrix"])
        self.assertEqual(code, 0)
        call.assert_not_called()
        self.assertRegex(out, r"Exact wall shear amplitude: 2\.2214\d*e-02")

    def test_run_without_openfoam_is_skipped(self):
        with tempfile.TemporaryDirectory() as d:
            a, b, c = self.isolated(d)
            with a, b, c, mock.patch.dict(os.environ, {}, clear=True), \
                    mock.patch.object(v.subprocess, "call") as call:
                code, out = self.run_main(["run"])
        self.assertEqual(code, 0)
        call.assert_not_called()
        self.assertIn("SKIPPED cutLink-ME-L1", out)

    def test_check_semantics_count_skipped_apart_from_pass(self):
        with tempfile.TemporaryDirectory() as d:
            a, b, c = self.isolated(d)
            with a, b, c:
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

    def test_check_case_tolerances(self):
        row = {"method": "penalty", "estimator": "", "level": 1,
               "shear_amp_err": "-0.2", "u_L2_all_ivp": "0.05"}
        ref = {"shear_amp_err": "-0.202", "u_L2_all_ivp": "0.05005"}
        self.assertTrue(all(s == "PASS" for s, _ in v.check_case(row, ref)))
        ref = {"shear_amp_err": "-0.30", "u_L2_all_ivp": "0.06"}
        self.assertEqual([s for s, _ in v.check_case(row, ref)], ["FAIL", "FAIL"])

    def test_check_analytic_passes(self):
        self.assertTrue(all(s == "PASS" for s, _ in v.check_analytic(P)))


class TrackedResultTests(unittest.TestCase):
    def test_tracked_rows_follow_the_exact_solution(self):
        rows = v.read_results()
        if not rows:
            self.skipTest("no tracked results")
        self.assertEqual([x for x in v.check_tracked(rows, P) if x[0] == "FAIL"], [])

    def test_every_method_and_level_of_the_matrix_is_tracked(self):
        rows = v.read_results()
        if not rows:
            self.skipTest("no tracked results")
        keys = {v.row_key(r) for r in rows}
        for method, est in v.VARIANTS:
            for level in v.LEVELS:
                self.assertIn((method, est, level, v.TUTORIAL_DT), keys)

    def test_results_md_contains_the_generated_table(self):
        rows = v.read_results()
        md = v.VERIFICATION_DIR / "results.md"
        if not rows or not md.exists():
            self.skipTest("no tracked results or results.md")
        self.assertIn(v.master_table(rows).strip(), md.read_text())


class CoveredVolumeTests(unittest.TestCase):
    def test_covered_volume_ratio(self):
        p = v.parameters()
        ny, nx, _, _ = v.mesh_shape(p, 2)
        h = p["height"] / ny
        rows = 5                      # five rows of cells below y_w = 0.1234
        occ = [1.0 if j < rows else 0.0 for j in range(ny) for _ in range(nx)]
        self.assertAlmostEqual(v.covered_volume_ratio(p, 2, occ),
                               rows * h / p["yw"], 12)
        # a partially covered row adds its fraction
        occ2 = [0.5 if j == rows else x for j, x in
                zip((j for j in range(ny) for _ in range(nx)), occ)]
        self.assertAlmostEqual(v.covered_volume_ratio(p, 2, occ2),
                               (rows + 0.5) * h / p["yw"], 12)

    def test_read_scalar_field(self):
        with tempfile.TemporaryDirectory() as d:
            path = Path(d) / "lambda"
            path.write_text("internalField   nonuniform List<scalar> \n3\n"
                            "(\n1\n0.5\n0\n)\n;\n")
            self.assertEqual(v.read_scalar_field(path), [1.0, 0.5, 0.0])
            path.write_text("internalField   uniform 0;\n")
            with self.assertRaises(SystemExit), \
                    contextlib.redirect_stderr(io.StringIO()):
                v.read_scalar_field(path)


class ReviewFixTests(unittest.TestCase):
    def test_aligned_inertia_is_a_half_step_shift(self):
        # a linear ramp of the inertia: the average of two steps is the
        # ramp half a step later; a constant is unchanged
        cols = {"t": [0.1 * i for i in range(5)],
                "force": [[0.0, 0, 0]] * 5,
                "inertia": [[float(i), 0, 0] for i in range(5)]}
        t, tau = v.aligned_shear_series(cols, 1.0)
        self.assertEqual(tau[:4], [0.5, 1.5, 2.5, 3.5])
        self.assertEqual(tau[4], 4.0)
        cols["inertia"] = [[2.0, 0, 0]] * 5
        self.assertEqual(v.aligned_shear_series(cols, 2.0)[1], [1.0] * 5)

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


class ReviewFixTests(unittest.TestCase):
    """Findings of the review of the pull request"""

    def run_main(self, argv):
        out = io.StringIO()
        with contextlib.redirect_stdout(out):
            code = v.main(argv)
        return code, out.getvalue()

    # -- the documented shorthand ./Allverify --runs L1 L2
    def test_normalise_argv(self):
        n = v.normalise_argv
        self.assertEqual(n([]), ["run"])
        self.assertEqual(n(["--runs", "L1", "L2"]), ["run", "--runs", "L1", "L2"])
        self.assertEqual(n(["--levels", "2"]), ["run", "--levels", "2"])
        for c in v.COMMANDS + ("-h", "--help"):
            self.assertEqual(n([c, "x"]), [c, "x"])

    def test_main_accepts_the_shorthand(self):
        with tempfile.TemporaryDirectory() as d:
            work = Path(d) / "work"
            with mock.patch.object(v, "WORK_DIR", work), \
                    mock.patch.object(v, "CASES_DIR", work / "cases"), \
                    mock.patch.object(v, "MANIFEST", work / "manifest.json"), \
                    mock.patch.dict(os.environ, {}, clear=True), \
                    mock.patch.object(v.subprocess, "call") as call:
                code, out = self.run_main(["--runs", "L1", "L2"])
        self.assertEqual(code, 0)
        call.assert_not_called()
        for name in ("cutLink-ME-L1", "penalty-L1", "ghostCell-L2"):
            self.assertIn(f"SKIPPED {name}", out)
        self.assertNotIn("SKIPPED cutLink-ME-L3", out)

    # -- OpenFOAM.com against OpenFOAM.org
    def env(self, **kw):
        base = {"WM_PROJECT_DIR": "/opt/openfoam"}
        base.update(kw)
        return mock.patch.dict(os.environ, base, clear=True)

    def test_openfoam_org_is_detected_by_the_version(self):
        for version in ("11", "12", "10"):
            with self.env(WM_PROJECT="OpenFOAM", WM_PROJECT_VERSION=version):
                self.assertIn("OpenFOAM.com", v.check_environment(1))
        with self.env(WM_PROJECT="OpenFOAM"):          # no version at all
            self.assertIn("OpenFOAM.com", v.check_environment(1))
        with self.env(WM_PROJECT="foam", WM_PROJECT_VERSION="v2412"):
            self.assertIn("OpenFOAM.com", v.check_environment(1))

    def test_openfoam_com_passes_the_version_rule(self):
        with self.env(WM_PROJECT="OpenFOAM", WM_PROJECT_VERSION="v2412"), \
                mock.patch.object(v.shutil, "which", return_value="/bin/x"):
            self.assertEqual(v.check_environment(1), "")
            self.assertEqual(v.check_environment(4), "")

    def test_openfoam_org_run_is_skipped_not_failed(self):
        with tempfile.TemporaryDirectory() as d:
            work = Path(d) / "work"
            with mock.patch.object(v, "WORK_DIR", work), \
                    mock.patch.object(v, "CASES_DIR", work / "cases"), \
                    mock.patch.object(v, "MANIFEST", work / "manifest.json"), \
                    self.env(WM_PROJECT="OpenFOAM", WM_PROJECT_VERSION="12"), \
                    mock.patch.object(v.subprocess, "call") as call:
                code, out = self.run_main(["run"])
        self.assertEqual(code, 0)
        call.assert_not_called()
        self.assertIn("SKIPPED cutLink-ME-L1", out)
        self.assertNotIn("FAILED", out)

    def test_exit_zero_without_a_log_is_skipped_not_failed(self):
        """Allrun.verify exits 0 without a solver log through the guards of
        solids4FoamScripts.sh (caseDoesNotRunWithOpenFOAMOrg)"""
        with tempfile.TemporaryDirectory() as d:
            work = Path(d) / "work"
            with mock.patch.object(v, "WORK_DIR", work), \
                    mock.patch.object(v, "CASES_DIR", work / "cases"), \
                    mock.patch.object(v, "MANIFEST", work / "manifest.json"), \
                    mock.patch.object(v, "check_environment", return_value=""), \
                    mock.patch.object(v.subprocess, "call", return_value=0):
                code, out = self.run_main(["run"])
                state = v.load_manifest()["cases"]["cutLink-ME-L1"]["state"]
        self.assertEqual(code, 0)
        self.assertIn("SKIPPED cutLink-ME-L1", out)
        self.assertEqual(state, "skipped")

    def test_a_failing_run_still_fails(self):
        with tempfile.TemporaryDirectory() as d:
            work = Path(d) / "work"
            with mock.patch.object(v, "WORK_DIR", work), \
                    mock.patch.object(v, "CASES_DIR", work / "cases"), \
                    mock.patch.object(v, "MANIFEST", work / "manifest.json"), \
                    mock.patch.object(v, "check_environment", return_value=""), \
                    mock.patch.object(v.subprocess, "call", return_value=3):
                code, out = self.run_main(["run"])
        self.assertEqual(code, 1)
        self.assertIn("FAILED cutLink-ME-L1: Allrun status 3", out)

    # -- the observed cost
    def test_observed_rank_hours_come_from_the_tracked_rows(self):
        obs = v.observed_rank_hours()
        rows = v.read_results()
        self.assertGreater(len(obs), 60)
        r = next(x for x in rows if x["method"] == "penalty"
                 and x["level"] == 4 and v.row_dt(x) == v.TUTORIAL_DT)
        self.assertAlmostEqual(obs[("penalty", "", 4, v.TUTORIAL_DT)],
                               float(r["rank_h"]))
        self.assertIn(("incremental", "", 4, 3.125e-05), obs)

    def test_dry_run_prints_the_observed_cost(self):
        code, out = self.run_main(["dry-run", "--runs", "matrix"])
        self.assertEqual(code, 0)
        self.assertNotIn("n/a", out)
        total = float(re.search(r"Observed total of the listed cases: "
                                r"([0-9.]+) rank-h", out).group(1))
        tracked = sum(float(r["rank_h"]) for r in v.read_results()
                      if v.row_dt(r) == v.TUTORIAL_DT
                      and r["estimator"] != "ST")
        self.assertGreater(total, 0.1)
        self.assertAlmostEqual(total, tracked, places=2)

    def test_dry_run_of_a_control_prints_the_tracked_cost(self):
        code, out = self.run_main(["dry-run", "--runs", "dt-3e-5"])
        self.assertEqual(code, 0)
        self.assertRegex(out, r"incremental-L4-dt3\.125e-05\s+2560 cells, "
                         r"1 rank\(s\), observed 1\.8058 allocated")

    def test_dry_run_without_tracked_cost_says_so(self):
        code, out = self.run_main(["dry-run", "--runs", "T8"])
        self.assertIn("n/a", out)


class ControlOptionTests(unittest.TestCase):
    """The control runs of reference/results/controls have a generator"""

    class A:
        methods, levels, ranks = [], [], 0
        end_time = wall_y = 0.0

        def __init__(self, *runs, **kw):
            self.runs = list(runs)
            self.__dict__.update(kw)

    def test_every_tracked_control_row_is_planned(self):
        planned = {(m, e, lv, dt) for m, e, lv, _, _, dt, _, _ in
                   v.plan(self.A("dt-L1", "dt-fine", "dt-6e-5", "dt-3e-5"))}
        tracked = set()
        for r in v.read_results(v.CONTROLS_CSV):
            e = "ME" if (r["method"] == "cutLink" and r["estimator"] == "ST") \
                else r["estimator"]
            tracked.add((r["method"], e, r["level"], v.row_dt(r)))
        self.assertEqual(planned, tracked)
        rows = v.read_results(v.CONTROLS_CSV)       # ST rows are the ME runs
        self.assertEqual(len(planned), sum(1 for r in rows
                                           if r["estimator"] != "ST"))

    def test_dt_steps_lie_on_the_write_grid(self):
        for _, _, _, _, _, dt, _, _ in v.plan(self.A("controls")):
            if dt:
                self.assertAlmostEqual(0.125 / dt, round(0.125 / dt), 9)

    def test_named_controls(self):
        wall = v.plan(self.A("wall-face"))
        self.assertEqual([x[3] for x in wall][:2],
                         ["cutLink-ME-L2-yw0.125", "cutLink-FO-L2-yw0.125"])
        self.assertEqual({x[2] for x in wall}, {2, 3, 4})
        self.assertTrue(all(x[7] == 0.125 and x[6] == 0 for x in wall))
        t8 = v.plan(self.A("T8"))
        self.assertEqual([x[3] for x in t8],
                         ["cutLink-ME-L3-T8", "penalty-L3-T8"])

    def test_end_time_and_wall_options_tag_the_name(self):
        plan = v.plan(self.A("smoke", end_time=8.0))
        self.assertEqual(plan[0][3], "cutLink-ME-L1-T8")
        self.assertEqual(plan[0][6], 8.0)
        plan = v.plan(self.A("L2", wall_y=0.125, methods=["penalty"]))
        self.assertEqual(plan[0][3], "penalty-L2-yw0.125")
        plan = v.plan(self.A("smoke", ranks=4, end_time=1.0))
        self.assertEqual(plan[0][3], "cutLink-ME-L1-T1-np4")
        self.assertEqual(v.plan(self.A("smoke"))[0][3], "cutLink-ME-L1")

    def test_end_time_changes_only_the_end_time(self):
        tut = (v.TUTORIAL_DIR / "system/controlDict").read_text()
        with tempfile.TemporaryDirectory() as d:
            case = v.generate_case("penalty", "", 3, "x", 1, end_time=8.0,
                                   cases_dir=Path(d))
            text = (case / "system/controlDict").read_text()
        self.assertEqual(v.dict_scalar(text, "endTime"), 8.0)
        self.assertEqual(text.replace("endTime         8;", "endTime         4;"),
                         tut)

    def test_wall_option_moves_only_the_top_vertices(self):
        src = v.TUTORIAL_DIR / "constant/triSurface/slab.stl"
        with tempfile.TemporaryDirectory() as d:
            case = v.generate_case("cutLink", "ME", 2, "x", 1,
                                   cases_dir=Path(d), wall_y=0.125)
            new = (case / "constant/triSurface/slab.stl").read_text()
            self.assertEqual(v.stl_top(case / "constant/triSurface/slab.stl"),
                             0.125)
        old = src.read_text()
        ys = lambda t: [float(x) for x in re.findall(
            r"vertex\s+\S+\s+(\S+)\s+\S+", t)]
        top = v.stl_top(src)
        self.assertEqual(len(ys(old)), len(ys(new)))
        for a, b in zip(ys(old), ys(new)):
            if abs(a - top) < 1e-6:
                self.assertEqual(b, 0.125)
            else:
                self.assertEqual(a, b)
        self.assertGreater(sum(abs(a - top) < 1e-6 for a in ys(old)), 0)
        self.assertEqual(len(old.splitlines()), len(new.splitlines()))
        # the rest of the case is the tutorial's
        self.assertNotEqual(old, new)

    def test_wall_outside_the_mesh_is_refused(self):
        with tempfile.TemporaryDirectory() as d, \
                self.assertRaises(SystemExit), \
                contextlib.redirect_stderr(io.StringIO()):
            v.generate_case("cutLink", "ME", 2, "x", 1, cases_dir=Path(d),
                            wall_y=0.7)

    def test_analysis_windows_and_exact_solutions_follow_the_wall(self):
        self.assertEqual(v.analysis_windows(0.0), [(3.0, 4.0)])
        self.assertEqual(v.analysis_windows(4.0), [(3.0, 4.0)])
        self.assertEqual(v.analysis_windows(8.0), [(3.0, 4.0), (7.0, 8.0)])
        self.assertEqual(v.snapshot_times((7.0, 8.0))[0], 7.125)

    def test_check_series_against_a_tracked_history(self):
        with tempfile.TemporaryDirectory() as d:
            path = Path(d) / "shear.csv"
            path.write_text("time,cutLink-ME\n3,1.0e-02\n3.001,2.0e-02\n")
            t = [3.0, 3.001, 3.002]
            self.assertEqual(v.check_series(path, "cutLink-ME", t,
                                            [1.0e-2, 2.0e-2, 0.0], "x")[0],
                             "PASS")
            self.assertEqual(v.check_series(path, "cutLink-ME", t,
                                            [1.0e-2, 2.1e-2, 0.0], "x")[0],
                             "FAIL")
            self.assertEqual(v.check_series(path, "cutLink-ME", [3.0],
                                            [1.0e-2], "x")[0], "FAIL")
            self.assertEqual(v.check_series(path, "nope", t, [0] * 3,
                                            "x")[0], "SKIPPED")

    def test_check_long_run(self):
        rows = [{"window": "3-4", "shear_amp_err_ivp": 0.0100},
                {"window": "7-8", "shear_amp_err_ivp": 0.01005}]
        self.assertEqual(v.check_long_run(rows, "c")[0][0], "PASS")
        rows[1]["shear_amp_err_ivp"] = 0.0115
        self.assertEqual(v.check_long_run(rows, "c")[0][0], "FAIL")
        self.assertEqual(v.check_long_run([{"window": "3-4"}], "c"), [])

    def test_the_tracked_wall_histories_match_the_exact_target(self):
        """The tracked shear-yw0125 histories are analysed against the
        solution of the moved wall, not the tutorial wall"""
        q = dict(P)
        q["yw"] = 0.125
        s = v.read_series_csv(v.CONTROLS_DIR / "shear-yw0125-L3.csv")
        m = v.shear_metrics(q, s["time"], s["cutLink-FO"])
        self.assertLess(abs(m["shear_amp_err_ivp"]), 5e-3)
        m0 = v.shear_metrics(P, s["time"], s["cutLink-FO"])
        self.assertGreater(abs(m0["shear_amp_err_ivp"]),
                           abs(m["shear_amp_err_ivp"]))

    def test_work_dir_option_separates_the_runs(self):
        with tempfile.TemporaryDirectory() as d:
            saved = (v.WORK_DIR, v.CASES_DIR, v.MANIFEST)
            try:
                with mock.patch.dict(os.environ, {}, clear=True), \
                        contextlib.redirect_stdout(io.StringIO()):
                    code = v.main(["run", "--work-dir", d])
                self.assertEqual(code, 0)
                self.assertTrue((Path(d) / "manifest.json").exists())
                self.assertTrue((Path(d) / "cases/cutLink-ME-L1").exists())
            finally:
                v.set_work_dir(saved[0])
        self.assertEqual((v.WORK_DIR, v.CASES_DIR, v.MANIFEST), saved)


def synthetic_case(case, level, yw=None, end_time=4.0, scale=1.0, dt=0.001):
    """A completed case whose fields and force are the periodic solution of
    a wall at yw (the tutorial wall if None), written for the fourth period
    and, for end_time > 4, the last"""
    q = dict(P)
    if yw:
        q["yw"] = yw
    ny, nx, _, area = v.mesh_shape(q, level)
    h = q["height"] / ny
    windows = v.analysis_windows(end_time)
    for period in windows:
        for tk in v.snapshot_times(period) + [period[0] + x for x in
                                              (0.0, 0.25, 0.5, 0.75)]:
            out = []
            for j in range(ny):
                s = (j + 0.5) * h - q["yw"]
                out += [(v.exact_u(q, s, tk) if s > 0 else 0.0, 0.0, 0.0)] * nx
            write_field(case / f"{tk:g}" / "U", out)
    n = int(round(end_time / dt))
    t = [dt * k for k in range(1, n + 1)]
    force = [scale * v.exact_tau(q, x) * area for x in t]
    write_forces(case / "postProcessing/immersedBoundary/0/slab.dat", t, force,
                 traction=force)
    (case / "system").mkdir(parents=True, exist_ok=True)
    (case / "system/controlDict").write_text(f"deltaT {dt};\n")
    (case / "wallclock.txt").write_text("wall_s=36 ranks=1 allrun_status=0\n")
    return case


class GeneratorTests(unittest.TestCase):
    """analyse_case for the controls and the export scripts"""

    def test_wall_case_is_analysed_against_the_moved_wall(self):
        with tempfile.TemporaryDirectory() as d:
            case = synthetic_case(Path(d), 2, yw=0.125)
            rows, series = v.analyse_case(case, "cutLink", "ME", 2,
                                          wall_y=0.125)
            wrong, _ = v.analyse_case(case, "cutLink", "ME", 2)
        self.assertEqual([r["window"] for r in rows], ["3-4", "3-4"])
        self.assertEqual(rows[0]["wall_y"], 0.125)
        self.assertLess(rows[0]["u_L2_all"], 1e-9)
        self.assertGreater(wrong[0]["u_L2_all"], 1e-3)
        self.assertEqual(set(series), {"ME", "ST"})

    def test_end_time_case_has_a_row_per_period(self):
        with tempfile.TemporaryDirectory() as d:
            case = synthetic_case(Path(d), 1, end_time=8.0)
            rows, _ = v.analyse_case(case, "cutLink", "FO", 1, end_time=8.0)
        self.assertEqual([r["window"] for r in rows], ["3-4", "7-8"])
        for r in rows:
            self.assertLess(abs(r["shear_amp_err"]), 1e-9)
            self.assertLess(r["u_L2_all"], 1e-9)
        self.assertEqual(rows[0]["steps"], 8000)
        # the fields are the periodic solution, so the error against the
        # initial-value solution (a decaying transient) differs by period
        self.assertEqual(v.check_long_run(rows, "c")[0][0], "FAIL")

    def test_baseline_case_rows_have_no_window(self):
        with tempfile.TemporaryDirectory() as d:
            rows, _ = v.analyse_case(synthetic_case(Path(d), 1), "cutLink",
                                     "FO", 1)
        self.assertNotIn("window", rows[0])

    def test_check_variant_compares_the_wall_case_with_the_tracked_history(self):
        with tempfile.TemporaryDirectory() as d:
            case = synthetic_case(Path(d), 2, yw=0.125)
            entry = {"level": 2, "wall_y": 0.125, "end_time": 0.0}
            rows, series = v.analyse_case(case, "cutLink", "ME", 2,
                                          wall_y=0.125)
            out = v.check_variant("c", entry, rows, series, {})
        # the synthetic history is the exact periodic shear, the tracked one
        # is a solver history: they differ, which must be a FAIL
        self.assertEqual({o[0] for o in out}, {"FAIL"})
        self.assertIn("shear-yw0125-L2.csv", out[0][1])

    def test_export_writes_the_tracked_file_formats(self):
        import export_results as ex
        with tempfile.TemporaryDirectory() as d:
            work = Path(d) / "work"
            names = {"cutLink-ME-L1": ("cutLink", "ME", 0.001),
                     "cutLink-FO-L1": ("cutLink", "FO", 0.001),
                     "penalty-L1-dt0.0025": ("penalty", "", 0.0025)}
            manifest = {"cases": {}}
            for name, (m, e, dt) in names.items():
                case = synthetic_case(work / "cases" / name, 1, dt=dt)
                manifest["cases"][name] = {
                    "method": m, "estimator": e, "level": 1, "dt": dt,
                    "ranks": 1, "end_time": 0.0, "wall_y": 0.0,
                    "state": "completed"}
            (work / "manifest.json").write_text(json.dumps(manifest))
            out = Path(d) / "out"
            with mock.patch.object(v, "WORK_DIR", work):
                info = ex.export_cases(ex.completed([work]), out)
            self.assertEqual(info, {"rows": 4, "cases.csv": 3, "controls": 1})
            main = v.read_results(out / "cases.csv")
            self.assertEqual({v.row_key(r)[:3] for r in main},
                             {("cutLink", "ME", 1), ("cutLink", "ST", 1),
                              ("cutLink", "FO", 1)})
            self.assertEqual(main[0]["rank_h"], "0.0100")   # 1 rank x 36 s
            ctrl = v.read_results(out / "controls/cases-timestep.csv")
            self.assertEqual([v.row_dt(r) for r in ctrl], [0.0025])
            shear = (out / "shear-L1.csv").read_text().splitlines()
            self.assertEqual(shear[0], "time,exact_periodic,exact_ivp,"
                             "cutLink-ME,cutLink-ST,cutLink-FO")
            self.assertEqual(len(shear), 1002)        # 3.0 ... 4.0 s
            self.assertTrue(shear[1].startswith("3,"))
            prof = (out / "profiles-L1.csv").read_text().splitlines()
            self.assertEqual(prof[0], "phase,s,exact_periodic,exact_ivp,cutLink")
            self.assertEqual(len(prof), 1 + 4 * 15)

    def test_export_wall_and_scaling(self):
        import export_results as ex
        with tempfile.TemporaryDirectory() as d:
            work = Path(d) / "work"
            manifest = {"cases": {}}
            for lab, est in (("ME", "ME"), ("FO", "FO")):
                name = f"cutLink-{lab}-L2-yw0.125"
                synthetic_case(work / "cases" / name, 2, yw=0.125)
                manifest["cases"][name] = {
                    "method": "cutLink", "estimator": est, "level": 2,
                    "dt": 0.0, "ranks": 1, "end_time": 0.0, "wall_y": 0.125,
                    "state": "completed"}
            for n, clock in ((1, 14), (2, 14), (4, 13)):
                name = f"penalty-L4-T1-np{n}"
                case = work / "cases" / name
                case.mkdir(parents=True)
                steps = "".join(f"Time = {k}\n" for k in range(1, 1001))
                (case / "log.solids4Foam").write_text(
                    f"{steps}ExecutionTime = 1 s  ClockTime = {clock} s\n")
                manifest["cases"][name] = {
                    "method": "penalty", "estimator": "", "level": 4,
                    "dt": 0.0, "ranks": n, "end_time": 1.0, "wall_y": 0.0,
                    "state": "completed"}
            (work / "manifest.json").write_text(json.dumps(manifest))
            cases = ex.completed([work])
            out = Path(d) / "out"
            files = ex.export_wall(cases, out)
            self.assertEqual([p.name for p in files], ["shear-yw0125-L2.csv"])
            lines = files[0].read_text().splitlines()
            self.assertEqual(lines[0], "time,cutLink-ME,cutLink-FO")
            self.assertEqual(len(lines), 1002)
            text = ex.export_scaling(cases, out / "scaling.csv")
        self.assertEqual(text.splitlines()[0],
                         "ranks,steps,solve_wall_s,speedup,efficiency")
        self.assertEqual(text.splitlines()[1:], ["1,1000,14,1.00,1.00",
                                                 "2,1000,14,1.00,0.50",
                                                 "4,1000,13,1.08,0.27"])

    def test_completed_refuses_a_case_in_two_work_directories(self):
        import export_results as ex
        with tempfile.TemporaryDirectory() as d:
            for k in "ab":
                (Path(d) / k).mkdir()
                (Path(d) / k / "manifest.json").write_text(json.dumps(
                    {"cases": {"x": {"state": "completed"}}}))
            with self.assertRaises(SystemExit), \
                    contextlib.redirect_stderr(io.StringIO()):
                ex.completed([Path(d) / "a", Path(d) / "b"])

    def test_scripts_have_no_numpy_and_no_absolute_home_path(self):
        here = Path(__file__).resolve().parent
        for path in here.glob("*.py"):
            if path.name.startswith("test_"):
                continue
            text = path.read_text()
            self.assertNotIn("/home/", text, path.name)
            self.assertNotRegex(text, r"(?m)^\s*(import|from)\s+(numpy|scipy)", path.name)


if __name__ == "__main__":
    unittest.main()
