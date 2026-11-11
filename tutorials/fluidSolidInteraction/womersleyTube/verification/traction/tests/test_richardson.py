"""Observed order / Richardson guards, and their use in decomposition_summary.py
and traction_checks.py (synthetic data only)."""
import json
import math
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

SCRIPTS = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(SCRIPTS))
import richardson as R  # noqa: E402


def series(limit, c, p, phase=1.0):
    """errors on levels h = 1, 1/2, 1/4 of limit + c h^p (c complex)."""
    return [limit + c * phase * (0.5 ** (p * k)) for k in range(3)]


class Observed(unittest.TestCase):
    def test_clean_second_order(self):
        e = series(1 + 2j, 0.3 - 0.1j, 2.0)
        info, lim2, limp = R.richardson_limits(*e)
        self.assertTrue(info["valid"])
        self.assertAlmostEqual(info["order"], 2.0, 12)
        self.assertAlmostEqual(abs(limp - (1 + 2j)), 0, 12)
        self.assertAlmostEqual(abs(lim2 - (1 + 2j)), 0, 12)

    def test_equal_fine_levels_no_zero_division(self):
        info, lim2, limp = R.richardson_limits(1 + 0j, 0.5 + 0j, 0.5 + 0j)
        self.assertFalse(info["valid"])
        self.assertTrue(math.isnan(info["order"]))
        self.assertEqual(limp, lim2)            # falls back to the formal limit

    def test_near_zero_order_is_not_used(self):
        # |d1|/|d2| = 2^0.05 -> p = 0.05 would give 1/(2^p - 1) ~ 29
        d2 = 1.0
        e4, e2 = 0.0, d2
        e1 = e2 + d2 * 2 ** 0.05
        info, lim2, limp = R.richardson_limits(e1, e2, e4)
        self.assertFalse(info["valid"])
        self.assertEqual(limp, lim2)
        self.assertLess(abs(limp - e4), 1.0)    # bounded, not blown up

    def test_negative_order_is_not_used(self):          # d1 < d2
        info, lim2, limp = R.richardson_limits(0.6 + 0j, 0.5 + 0j, 0.0 + 0j)
        self.assertFalse(info["valid"])
        self.assertLess(info["order"], 0)
        self.assertEqual(limp, lim2)

    def test_phase_inconsistent_series_is_non_asymptotic(self):
        # increments rotate by 60 degrees per level; magnitude ratio 4 (p = 2)
        d2 = 1.0 + 0j
        d1 = 4.0 * complex(math.cos(math.radians(60)), math.sin(math.radians(60)))
        e4 = 0j
        e2 = e4 + d2
        e1 = e2 + d1
        info = R.observed_order(e1, e2, e4)
        self.assertAlmostEqual(info["order"], 2.0, 12)
        self.assertAlmostEqual(info["arg_deg"], 60.0, 9)
        self.assertFalse(info["valid"])
        self.assertIn("non-asymptotic", info["reason"])

    def test_small_phase_still_valid(self):
        d1 = 4.0 * complex(math.cos(math.radians(10)), math.sin(math.radians(10)))
        self.assertTrue(R.observed_order(d1 + 1, 1 + 0j, 0j)["valid"])

    def test_non_finite(self):
        self.assertFalse(R.observed_order(complex("nan"), 1, 0)["valid"])


def write(path, obj):
    path.write_text(json.dumps(obj))


QS = ("attenuation", "speed", "wallMid_amp", "wallMid_phase")
TQS = ("Vn_ratio", "VnEst_ratio", "Vn_err/|p|", "P_err/|p|", "Sn_err/|p|", "Vx_err/|Vx|")
SIGQ = ("P_err/|p|", "Vn_err/|p|", "Sn_err/|p|", "Vx_err/|Vx|", "solid_vs_fluid_n/|p|")


def pair(z):
    return [z.real, z.imag]


def make_inputs(work: Path, mesh_levels, slope_stored=2.5e-3, slope_values=None):
    """Synthetic inputs for decomposition_summary.py; mesh_levels = (e1, e2, e4)
    of every traction quantity; the amplitude values are linear with
    slope_values but the stored fit says slope_stored (to prove it is read)."""
    slope_values = slope_stored if slope_values is None else slope_values
    lim = {"limits": {q: {"E_m4_n200": -2e-3, "E_m4_n400": -1.5e-3,
                          "space_increment_m4_p2": 1e-4, "space_increment_m4_pobs": 1.2e-4}
                      for q in QS},
           "cases": {n: {q: -2.1e-3 for q in QS} for n in ("robin_m1_n200", "robin_m2_n200")},
           "prediction_modified_normal_viscous": {"f=0.333": dict(
               {q: 0.0 for q in QS}, finiteTube={q: 0.0 for q in QS})}}
    write(work / "attenuation_limit_mu5.json", lim)
    amp = {str(a): {q: slope_values * a for q in QS} for a in (0.5, 1.0, 2.0, 4.0)}
    write(work / "amplitude_series.json", {
        "case": "robin_m2_n200", "values": amp,
        "linearFit": {q: {"E0": 0.0, "slope_per_Pa": slope_stored, "misfit": 7e-9,
                          "m_observed": []} for q in QS}})
    e1, e2, e4 = mesh_levels
    ex = {"P": [1.0, 0.0], "Vn": [3e-4, 1e-4], "Sn": [1, 0], "Vx": [1, 0]}
    cases = {}
    for name, e in (("robin_m1_n200", e1), ("robin_m2_n200", e2), ("robin_m4_n200", e4),
                    ("robin_m4_n400", e4 * 0.99)):
        cases[name] = {"errors": {q: pair(e) for q in TQS + SIGQ}}
    for study in ("mu5", "amp2", "tol"):
        d = work / study / "womersleyTube/verification/postProcessing"
        d.mkdir(parents=True)
        c = dict(cases) if study == "mu5" else {"robin_m2_n200": cases["robin_m2_n200"]}
        write(d / "interfaceTraction_summary.json", {"exact": ex, "cases": c})


class DecompositionSummary(unittest.TestCase):
    def run_script(self, levels, **kw):
        with tempfile.TemporaryDirectory() as tmp:
            work = Path(tmp) / "w"
            work.mkdir()
            make_inputs(work, levels, **kw)
            out = Path(tmp) / "out"
            r = subprocess.run([sys.executable, "-I", str(SCRIPTS / "decomposition_summary.py"),
                                str(work), str(out)], capture_output=True, text=True)
            res = json.loads((out / "findings_summary.json").read_text()) if r.returncode == 0 else None
            return r, res

    def test_uses_stored_linear_fit_slope(self):
        r, res = self.run_script((0.4 + 0j, 0.2 + 0j, 0.1 + 0j), slope_stored=2.5e-3, slope_values=1.0e-3)
        self.assertEqual(r.returncode, 0, r.stderr)
        for q in QS:
            self.assertEqual(res["waveQoIs"][q]["amplitude_slope_per_Pa_m2"], 2.5e-3)
            self.assertEqual(res["waveQoIs"][q]["amplitude_fit_misfit"], 7e-9)

    def test_missing_stored_fit_is_a_clear_error(self):
        with tempfile.TemporaryDirectory() as tmp:
            work = Path(tmp) / "w"
            work.mkdir()
            make_inputs(work, (0.4 + 0j, 0.2 + 0j, 0.1 + 0j))
            f = work / "amplitude_series.json"
            d = json.loads(f.read_text())
            del d["linearFit"]
            f.write_text(json.dumps(d))
            r = subprocess.run([sys.executable, "-I", str(SCRIPTS / "decomposition_summary.py"),
                                str(work), str(Path(tmp) / "o")], capture_output=True, text=True)
            self.assertNotEqual(r.returncode, 0)
            self.assertIn("linearFit", r.stderr)

    def test_converged_fine_levels_no_crash(self):          # e2 == e4
        r, res = self.run_script((0.4 + 0j, 0.1 + 0j, 0.1 + 0j))
        self.assertEqual(r.returncode, 0, r.stderr)
        t = res["traction"]["Vn_err/|p|"]
        self.assertFalse(t["mesh_order_valid"])
        self.assertEqual(t["limit_a1_pobs"], t["limit_a1_p2"])

    def test_near_zero_order_is_bounded(self):
        r, res = self.run_script((0.1 + 1.034 * 0.1 + 0j, 0.1 + 0j, 0j))
        self.assertEqual(r.returncode, 0, r.stderr)
        t = res["traction"]["Vn_err/|p|"]
        self.assertFalse(t["mesh_order_valid"])
        self.assertLess(abs(complex(*t["limit_a1_pobs"])), 1.0)

    def test_valid_series_is_extrapolated_with_observed_order(self):
        r, res = self.run_script((0.4 + 0j, 0.2 + 0j, 0.1 + 0j))   # p = 1
        t = res["traction"]["Vn_err/|p|"]
        self.assertTrue(t["mesh_order_valid"])
        self.assertAlmostEqual(t["mesh_order"], 1.0, 9)


class TractionChecks(unittest.TestCase):
    def run_checks(self, vx):
        with tempfile.TemporaryDirectory() as tmp:
            d = Path(tmp) / "verification/postProcessing"
            d.mkdir(parents=True)
            zero = [0.0, 0.0]
            err = {"mom_vs_transfer_n/|Vn|": zero, "mom_vs_transfer_x/|Vx|": zero,
                   "solid_vs_fluid_n/|p|": zero, "solid_vs_fluid_x/|Vx|": zero}
            write(d / "interfaceTraction_summary.json", {"cases": {
                n: {"errors": dict(err, **{"Vx_err/|Vx|": pair(v)})}
                for n, v in zip(("robin_m1_n200", "robin_m2_n200", "robin_m4_n200"), vx)}})
            return subprocess.run([sys.executable, "-I", str(SCRIPTS / "traction_checks.py"), tmp],
                                  capture_output=True, text=True)

    def test_second_order_same_phase_passes(self):
        r = self.run_checks([4e-4 + 0j, 1e-4 + 0j, 0.25e-4 + 0j])
        self.assertEqual(r.returncode, 0, r.stdout)
        self.assertIn("arg(ratio)", r.stdout)

    def test_magnitude_only_second_order_with_rotating_phase_fails(self):
        ph = complex(math.cos(math.radians(70)), math.sin(math.radians(70)))
        r = self.run_checks([4e-4 * ph + 1e-4, 1e-4 + 0j, 0j + 0.0])   # d1/d2 = 4*ph*(...)
        self.assertEqual(r.returncode, 1)
        self.assertIn("FAIL: shear traction mesh order", r.stdout)


if __name__ == "__main__":
    unittest.main()
