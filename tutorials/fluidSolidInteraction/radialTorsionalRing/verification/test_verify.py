#!/usr/bin/env python3
"""Automated tests of the verification driver without CFD (synthetic logs in $TMPDIR, plus read-only
checks of existing logs when the paths in the environment variable RT_KNOWN_RUNS are available).

Run:  python3 -I verification/test_verify.py            (synthetic tests only)
      RT_KNOWN_RUNS=/path/to/runs python3 -I verification/test_verify.py   (also known-good/bad real logs)
"""
import sys
sys.dont_write_bytecode = True
import hashlib
import math
import json
import os
import shutil
import subprocess
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "scripts"))
import acceptance as ACC  # noqa: E402
import family as FAM  # noqa: E402
import verify  # noqa: E402

with open(os.path.join(HERE, "thresholds.json")) as _f:
    THR = json.load(_f)


def make_log(n=60, dt=0.01, iters=3, bad_step=None, truncate=False, fatal=None, nan_step=None,
             accept="r2", abs_step=None, snes_bad_step=None, p_per=3, end=True, stall_step=None, r2_accepted_above=None, per_step=None):
    """Synthetic solids4Foam log in the real format. accept: criterion of an ordinary step."""
    L = ["Create time", ""]
    for i in range(1, n + 1):
        t = round(i * dt, 10)
        L += [f"Time = {t:g}", "", "Setting traction on solid interfaces", "Solving the momentum equation for D using PETSc SNES"]
        L += ["Nonlinear solve converged due to CONVERGED_FNORM_RELATIVE iterations 3",
              "FSI relative residual1 norm for interface 0: 0.5", "FSI residual2 norm for interface 0: 0.5",
              "FSI absolute residual (max point) for interface 0: 1e-3", ""]
        k_end = iters if bad_step != i else 40
        for k in range(1, k_end + 1):
            L.append(f"Time = {t:g}, iteration: {k}")
            for _pass in range(2):
                L.append("DILUPBiCG:  Solving for Ux, Initial residual = %g, Final residual = 5e-14, No Iterations 3" % (0.5 ** k))
                L.append("DILUPBiCG:  Solving for Uy, Initial residual = 1e-3, Final residual = 5e-14, No Iterations 3")
                for _c in range(3):
                    for _p in range(p_per if not (nan_step is None and False) else 3):
                        L.append("GAMG:  Solving for p, Initial residual = 1e-3, Final residual = 5e-13, No Iterations 4")
                    L.append("time step continuity errors : sum local = 1e-18, global = 1e-19, cumulative = 1e-18")
            L.append("Nonlinear solve converged due to CONVERGED_FNORM_RELATIVE iterations 3" if snes_bad_step != i else
                     "Nonlinear solve did not converge due to DIVERGED_MAX_IT iterations 50")
            last = k == k_end
            if not last:
                r1, r2, ab = 0.5, 0.5, 1e-3
            elif stall_step == i:
                r1, r2, ab = 0.5, 0.5, 1e-3
            elif abs_step == i:
                r1, r2, ab = 3e-3, 4e-3, 5e-13
            elif r2_accepted_above == i:
                r1, r2, ab = 2e-6, 3e-6, 1e-3
            elif accept == "both":
                r1, r2, ab = 5e-7, 5e-7, 1e-3
            elif accept == "r1":
                r1, r2, ab = 5e-7, 5e-3, 1e-3
            else:
                r1, r2, ab = 5e-3, 5e-7, 1e-3
            if per_step and i in per_step and last:
                r1, r2, ab = per_step[i][:3]
            if nan_step == i and last:
                r1 = float("nan")
            L += [f"FSI relative residual1 norm for interface 0: {r1:g}", f"FSI residual2 norm for interface 0: {r2:g}",
                  f"FSI absolute residual (max point) for interface 0: {ab:g}"]
            if last and (abs_step == i or (per_step and i in per_step and len(per_step[i]) > 3 and per_step[i][3])):
                L.append("FSI interface 0 accepted by outerCorrAbsoluteTolerance: 5e-13 < 1e-12")
        if stall_step == i:
            L += ["", "--> FOAM FATAL ERROR: (openfoam-2512)", "FSI coupling did not converge after 100 outer corrections. Final residual: 0.5, tolerance: 1e-06"]
            return "\n".join(L) + "\n"
        if fatal and fatal[0] == i:
            L += ["[0] --> FOAM FATAL ERROR: (openfoam-2512)", fatal[1]]
            return "\n".join(L) + "\n"
        if truncate and i == n:
            L = L[:-8]
            return "\n".join(L) + "\n"
        L += [f"ExecutionTime = {i} s  ClockTime = {i} s", ""]
    if end:
        L += ["End", ""]
    return "\n".join(L) + "\n"


def run_a(text, **cfg):
    d = tempfile.mkdtemp(dir=os.environ.get("TMPDIR"))
    p = os.path.join(d, "log.solids4Foam")
    with open(p, "w") as fh:
        fh.write(text)
    c = dict(endTime=cfg.pop("endTime", 0.6), deltaT=0.01, nNonOrth=2, tol=1e-6, absTol=1e-12, nOuterCorr=100)
    c.update(cfg)
    return ACC.run_layer_a(d, THR, cfg_override=c)


class Synthetic(unittest.TestCase):
    def test_good_run_is_conditional_not_pass_when_single_measure(self):
        r = run_a(make_log())
        self.assertEqual(r["status"], "CONDITIONAL")
        self.assertEqual(r["checks"]["termination"]["status"], "PASS")
        self.assertEqual(r["criterion_counts"]["r2-only"], 60)

    def test_both_measures_is_pass(self):
        r = run_a(make_log(accept="both"))
        self.assertEqual(r["checks"]["fsi_step_acceptance"]["status"], "PASS")
        self.assertEqual(r["status"], "PASS" if all(v["status"] in ("PASS", "INFO", "NOT-VERIFIED") for v in r["checks"].values()) else r["status"])

    def test_r1_only_is_flagged_with_r2_above_tolerance(self):
        r = run_a(make_log(accept="r1"))
        self.assertEqual(r["criterion_counts"]["r1-only"], 60)
        self.assertEqual(r["r1_only_steps_with_r2_above_tol"]["count"], 60)
        self.assertEqual(r["checks"]["fsi_step_acceptance"]["status"], "CONDITIONAL")

    def test_step_above_iteration_limit_fails(self):
        r = run_a(make_log(bad_step=50, accept="both"))
        self.assertEqual(r["checks"]["fsi_iterations"]["status"], "FAIL")
        self.assertEqual(r["status"], "FAIL")

    def test_iteration_above_limit_in_startup_is_excluded(self):
        r = run_a(make_log(bad_step=20, accept="both"))     # t = 0.2 s < 0.25 s
        self.assertEqual(r["checks"]["fsi_iterations"]["status"], "PASS")

    def test_accepted_step_with_residuals_above_tolerance_is_flagged(self):
        r = run_a(make_log(accept="both", abs_step=30))
        self.assertEqual(r["criterion_counts"]["absolute"], 1)
        self.assertIn("both relative residuals above tolerance", next(s["flags"] for s in r["per_step"] if s["criterion"] == "absolute"))
        self.assertEqual(r["checks"]["fsi_step_acceptance"]["status"], "CONDITIONAL")
        self.assertEqual(r["abs_accepted_times"], [0.3])

    def test_step_ended_without_any_criterion_fails(self):
        # final iteration with residuals above tolerance and no absolute acceptance, but the step closes
        t = make_log(accept="both")
        t = t.replace("FSI residual2 norm for interface 0: 5e-07", "FSI residual2 norm for interface 0: 5e-03", 3)
        t = t.replace("FSI relative residual1 norm for interface 0: 5e-07", "FSI relative residual1 norm for interface 0: 5e-03", 3)
        r = run_a(t)
        self.assertGreaterEqual(r["criterion_counts"]["none"], 1)
        self.assertEqual(r["checks"]["fsi_step_acceptance"]["status"], "FAIL")
        self.assertEqual(r["status"], "FAIL")

    def test_did_not_converge_fatal_fails(self):
        r = run_a(make_log(accept="both", stall_step=30))
        self.assertEqual(r["status"], "FAIL")
        self.assertEqual(r["checks"]["termination"]["status"], "FAIL")

    def test_truncated_log_without_error_is_incomplete(self):
        r = run_a(make_log(accept="both", truncate=True, end=False))
        self.assertEqual(r["checks"]["termination"]["status"], "INCOMPLETE")
        self.assertEqual(r["status"], "INCOMPLETE")

    def test_end_but_short_of_endtime_is_incomplete(self):
        r = run_a(make_log(accept="both", n=40), endTime=0.6)
        self.assertEqual(r["checks"]["termination"]["status"], "INCOMPLETE")

    def test_nan_fails(self):
        r = run_a(make_log(accept="both", nan_step=10))
        self.assertEqual(r["status"], "FAIL")

    def test_snes_failure_fails(self):
        r = run_a(make_log(accept="both", snes_bad_step=12))
        self.assertEqual(r["status"], "FAIL")

    def test_p_solves_per_continuity(self):
        r = run_a(make_log(accept="both", p_per=2))
        self.assertEqual(r["checks"]["p_solves_per_continuity"]["status"], "FAIL")

    def test_late_linear_final_exceedance_fails(self):
        t = make_log(accept="both")
        t = t.replace("Final residual = 5e-13, No Iterations 4", "Final residual = 3e-9, No Iterations 1000")
        r = run_a(t)
        self.assertEqual(r["checks"]["linear_finals"]["status"], "FAIL")

    def test_layerB_invalid_input_when_layerA_fails(self):
        import verify
        runs = {("torsion", m): {"status": "FAIL" if m == "M3" else "PASS"} for m in ("M1", "M2", "M3")}
        rep = verify.layer_b(tempfile.mkdtemp(), runs, THR)
        self.assertEqual(rep["families"]["torsion"]["status"], "INVALID-INPUT")
        self.assertEqual(rep["claims"]["torsion_primary_orders"]["status"], "INVALID-INPUT")

    def test_not_verified_list_is_printed(self):
        import verify
        d = tempfile.mkdtemp(dir=os.environ.get("TMPDIR"))
        os.makedirs(os.path.join(d, "torsion_M1"))
        with open(os.path.join(d, "torsion_M1", "log.solids4Foam"), "w") as fh:
            fh.write(make_log(accept="both"))
        out = subprocess.run([sys.executable, "-I", os.path.join(HERE, "scripts", "verify.py"), d, "--out", os.path.join(d, "o")],
                             capture_output=True, text=True)
        self.assertIn("NOT-VERIFIED (printed in every run)", out.stdout)
        self.assertIn("three-way radial term d_u3", out.stdout)
        self.assertEqual(out.returncode, 1)     # run is short of endTime 2 s: INCOMPLETE


TUT = os.path.abspath(os.path.join(HERE, ".."))
MAKE = os.path.join(HERE, "scripts", "makeCase.py")


def make_family(cfg, root, absTol=None, mods=None):
    """Three generated cases root/<cfg>_<mesh> (dictionaries only, serial); mods: {mesh: [makeCase options]}."""
    mods = mods or {}
    for m in ("M1", "M2", "M3"):
        cmd = [sys.executable, "-B", "-I", MAKE, os.path.join(root, f"{cfg}_{m}"), "--mesh", m, "--config", cfg, "--np", "1", "--no-reference"]
        if absTol is not None:
            cmd += ["--absTol", str(absTol)]
        cmd += mods.get(m, [])
        subprocess.run(cmd, check=True, capture_output=True)


class FamilyConsistency(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.root = tempfile.mkdtemp(dir=os.environ.get("TMPDIR"))
        cls.base = os.path.join(cls.root, "base")
        make_family("torsion", cls.base, absTol=1e-12)

    def fresh(self):
        d = tempfile.mkdtemp(dir=self.root)
        shutil.copytree(self.base, d, dirs_exist_ok=True)
        return d

    def test_identical_settings_pass(self):
        r = FAM.check_family(self.base, "torsion", THR)
        self.assertEqual(r["status"], "PASS", [f for f in r["findings"] if f["status"] != "PASS"])

    def test_documented_tolerance_pattern_is_conditional(self):
        d = tempfile.mkdtemp(dir=self.root)
        make_family("torsion", d)       # default per-mesh absTol: 1e-12, 1e-12, 1e-10
        r = FAM.check_family(d, "torsion", THR)
        self.assertEqual(r["status"], "CONDITIONAL")
        self.assertEqual([f["check"] for f in r["findings"] if f["status"] == "CONDITIONAL"], ["coupling_outerCorrAbsoluteTolerance"])

    def test_undocumented_tolerance_is_invalid(self):
        d = tempfile.mkdtemp(dir=self.root)
        make_family("torsion", d, absTol=1e-12, mods={"M3": ["--absTol", "1e-9"]})
        self.assertEqual(FAM.check_family(d, "torsion", THR)["status"], "INVALID-INPUT")

    def _mutated(self, args):
        d = tempfile.mkdtemp(dir=self.root)
        make_family("torsion", d, absTol=1e-12, mods={"M2": args})
        return FAM.check_family(d, "torsion", THR)

    def test_dt_mismatch_invalid(self):
        self.assertEqual(self._mutated(["--dt", "0.005"])["status"], "INVALID-INPUT")

    def test_nnonorth_mismatch_invalid(self):
        self.assertEqual(self._mutated(["--nNonOrth", "0"])["status"], "INVALID-INPUT")

    def test_drive_mismatch_invalid(self):
        r = self._mutated(["--theta0", "0.05"])
        self.assertEqual(r["status"], "INVALID-INPUT")
        self.assertTrue(any(f["check"] == "drives_M2" and f["status"] == "INVALID-INPUT" for f in r["findings"]))

    def test_wrong_refinement_ratio_invalid(self):
        r = self._mutated(["--ntheta", "96"])
        self.assertTrue(any(f["check"] == "mesh_family" and f["status"] == "INVALID-INPUT" for f in r["findings"]))

    def test_missing_run_incomplete(self):
        d = self.fresh()
        shutil.rmtree(os.path.join(d, "torsion_M3"))
        self.assertEqual(FAM.check_family(d, "torsion", THR)["status"], "INCOMPLETE")

    def test_unreadable_dictionary_incomplete(self):
        d = self.fresh()
        os.remove(os.path.join(d, "torsion_M2", "constant", "fsiProperties"))
        self.assertEqual(FAM.check_family(d, "torsion", THR)["status"], "INCOMPLETE")

    def test_conditional_family_caps_layer_b_pass(self):
        self.assertEqual(verify.cap("PASS", "CONDITIONAL"), "CONDITIONAL")
        self.assertEqual(verify.cap("PASS", "PASS"), "PASS")
        self.assertEqual(verify.cap("FAIL", "CONDITIONAL"), "FAIL")


class Bands(unittest.TestCase):
    def test_marginal_values_are_conditional(self):
        B, N = THR["layerB_deterministic"], THR["numerical_noise"]
        nz = N["combined_max_norm_noise_of_E"]["value"]
        lo, hi = B["combined_floor_band_M3"]["value"]
        self.assertEqual(verify.band_status(5.15e-3, lo, hi, 0.2, nz, nz), "CONDITIONAL")      # sig_rt E_M3
        lim = B["mixed_term_max_error_over_signal_M3_0_1p8"]["value"]["d_sig_rt"]
        a = N["mixed_term_allowance"]["documented_values"]
        self.assertEqual(verify.band_status(0.0336, 0.0, lim, 0.2, 0.0, 0.00862), "CONDITIONAL")   # allowance of 0-1.8 s: 0.86%
        self.assertEqual(verify.band_status(0.0247, 0.0, B["mixed_term_max_error_over_signal_M3_0_2"]["value"]["d_T"], 0.2, 0.0, a["d_T"]), "CONDITIONAL")
        self.assertEqual(verify.band_status(0.0265, 0.0, B["mixed_term_max_error_over_signal_M3_0_1p8"]["value"]["d_T"], 0.2, 0.0, 0.00877), "CONDITIONAL")
        self.assertEqual(verify.band_status(0.0155, 0.0, 0.025, 0.2, 0.0, a["d_alpha"]), "PASS")

    def test_far_outside_band_fails(self):
        self.assertEqual(verify.band_status(0.02, 2.5e-3, 5.3e-3, 0.2), "FAIL")

    def test_deterministic_and_noise_sections_are_separate(self):
        self.assertIn("layerB_deterministic", THR)
        self.assertIn("numerical_noise", THR)
        self.assertNotIn("noise", " ".join(THR["layerB_deterministic"]).replace("_comment", ""))


class CorrectedCriteria(unittest.TestCase):
    def test_late_absolute_acceptance_is_conditional_not_failure(self):
        r = run_a(make_log(accept="both", abs_step=50))           # t = 0.5 s, after the start-up window
        self.assertEqual(r["checks"]["fsi_step_acceptance"]["status"], "CONDITIONAL")
        self.assertNotEqual(r["status"], "FAIL")
        self.assertEqual(r["checks"]["fsi_step_acceptance"]["abs_accepted_after_startup"], [0.5])

    def test_start_up_iteration_exception_matches_m2_combined(self):
        r = run_a(make_log(accept="both", bad_step=9))            # 40 iterations at t = 0.09 s, like the 34 at 0.0925 s
        self.assertEqual(r["checks"]["fsi_iterations"]["status"], "PASS")
        self.assertEqual(r["checks"]["fsi_iterations"]["max_overall"], 40)

    def test_readme_thresholds_driver_agree(self):
        A = THR["layerA"]
        with open(os.path.join(TUT, "README.md")) as fh:
            txt = fh.read()
        self.assertIn(f"<= {A['max_fsi_iterations']['value']}", txt)
        self.assertIn(f"{A['max_fsi_iterations_from_t']['value']} s", txt)
        self.assertIn("CONDITIONAL", txt)
        self.assertIn("34", txt)
        self.assertEqual(A["abs_accept_startup_end_t"]["value"], 0.35)


class AllrunPaths(unittest.TestCase):
    """Allrun option logic with stubbed OpenFOAM helper functions (no solver, no CFD)."""

    @classmethod
    def setUpClass(cls):
        cls.root = tempfile.mkdtemp(dir=os.environ.get("TMPDIR"))
        cls.stub = os.path.join(cls.root, "stub")
        os.makedirs(os.path.join(cls.stub, "wm", "bin", "tools"))
        with open(os.path.join(cls.stub, "wm", "bin", "tools", "RunFunctions"), "w") as fh:
            fh.write("# stub\n")
        with open(os.path.join(cls.stub, "solids4FoamScripts.sh"), "w") as fh:
            fh.write("""
solids4Foam::caseDoesNotRunWithFoamExtend() { :; }
solids4Foam::caseDoesNotRunWithOpenFOAMOrg() { :; }
solids4Foam::requirePetscOrExitSilently() { :; }
solids4Foam::convertCaseFormat() { :; }
solids4Foam::runApplication() { echo "runApplication $*" >> "$PWD/stub.calls"; }
solids4Foam::runParallel() { echo "runParallel $*" >> "$PWD/stub.calls"; }
""")

    def case(self):
        d = tempfile.mkdtemp(dir=self.root)
        for item in os.listdir(TUT):
            if item in ("regressionTests",):
                continue
            src = os.path.join(TUT, item)
            if item == "verification":
                os.makedirs(os.path.join(d, "verification"))
                shutil.copytree(os.path.join(src, "scripts"), os.path.join(d, "verification", "scripts"),
                                ignore=shutil.ignore_patterns("__pycache__"))
            else:
                (shutil.copytree if os.path.isdir(src) else shutil.copy2)(src, os.path.join(d, item))
        return d

    def run_allrun(self, d, *args):
        env = dict(os.environ, WM_PROJECT_DIR=os.path.join(self.stub, "wm"), PATH=self.stub + os.pathsep + os.environ["PATH"],
                   PYTHONDONTWRITEBYTECODE="1")
        return subprocess.run(["bash", "./Allrun", *args], cwd=d, env=env, capture_output=True, text=True)

    @staticmethod
    def sha(d, rel):
        with open(os.path.join(d, rel), "rb") as fh:
            return hashlib.sha256(fh.read()).hexdigest()

    def test_default_uses_stored_dictionaries(self):
        d = self.case()
        before = {r: self.sha(d, r) for r in ("system/controlDict", "system/fluid/fvSolution", "constant/fsiProperties", "system/fluid/blockMeshDict")}
        r = self.run_allrun(d)
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        self.assertIn("Using the stored dictionaries", r.stdout)
        self.assertNotIn("regenerating", r.stdout)
        self.assertEqual(before, {k: self.sha(d, k) for k in before})
        self.assertTrue(os.path.isfile(os.path.join(d, "reference.json")))
        self.assertFalse(os.path.exists(os.path.join(d, "caseSetup.json")))      # makeCase.py was not run

    def test_alltest_style_edits_survive_default_run(self):
        d = self.case()
        p = os.path.join(d, "system", "controlDict")
        with open(p) as fh:
            txt = fh.read()
        txt = txt.replace("stopAt          endTime;", "stopAt           nextWrite;").replace("writeInterval   0.1;", "writeInterval 1;").replace("adjustableRunTime", "timeStep")
        with open(p, "w") as fh:
            fh.write(txt)
        h = self.sha(d, "system/controlDict")
        r = self.run_allrun(d)
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        self.assertEqual(h, self.sha(d, "system/controlDict"))
        self.assertIn("nextWrite", open(p).read())

    def test_override_regenerates(self):
        d = self.case()
        r = self.run_allrun(d, "M1", "combined", "serial", "endTime=0.05")
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        self.assertIn("regenerating", r.stdout)
        self.assertRegex(open(os.path.join(d, "system", "controlDict")).read(), r"endTime\s+0.05;")

    def test_each_override_option_regenerates(self):
        for opt in ("nNonOrth=0", "absTol=1e-10", "regenerate", "M2", "torsion"):
            with self.subTest(opt=opt):
                d = self.case()
                r = self.run_allrun(d, opt)
                self.assertIn("regenerating", r.stdout, opt)

    def test_stored_set_must_be_m1(self):
        d = self.case()
        r = self.run_allrun(d, "M2", "endTime=0.01")
        self.assertEqual(r.returncode, 0)
        r2 = self.run_allrun(d)            # now the stored mesh is M2: default path must refuse
        self.assertNotEqual(r2.returncode, 0)
        self.assertIn("not the M1 set", r2.stdout + r2.stderr)


MACHINE = ("/ho" + "me/", "/m" + "nt/", "sai" + "ram")   # written in pieces so that this file passes its own check


class Hygiene(unittest.TestCase):
    def test_no_bytecode_or_machine_paths_in_tree(self):
        bad, paths = [], []
        for dp, dn, fn in os.walk(TUT):
            if "__pycache__" in dp:
                bad.append(dp)
            if "regressionTests" in dp or os.sep + "work" + os.sep in dp:
                continue
            for f in fn:
                if f.endswith(".pyc"):
                    bad.append(os.path.join(dp, f))
                if f.endswith((".py", ".sh", ".md", ".json")) or f in ("Allrun", "Allclean", "Allverify"):
                    with open(os.path.join(dp, f), errors="replace") as fh:
                        t = fh.read()
                    if f != "test_verify.py" and (MACHINE[0] in t or MACHINE[1] in t or MACHINE[2] in t.lower()):
                        paths.append(os.path.join(dp, f))
        self.assertEqual(bad, [])
        self.assertEqual(paths, [])

import random
import numpy as np
import status as ST  # noqa: E402
import convergence as CV  # noqa: E402
import make_reference  # noqa: E402
import rt_reference as rt  # noqa: E402


class NoiseNormalisation(unittest.TestCase):
    """M1: the allowance is noise x max|reference quantity| / mixed signal, from the actual magnitudes (reference.json data)."""

    @classmethod
    def setUpClass(cls):
        cls.root = tempfile.mkdtemp(dir=os.environ.get("TMPDIR"))
        for cfg, th, P0 in (("combined", 0.1, -3000.0), ("torsion", 0.1, 0.0)):
            p = rt.Params(theta0=th, Omega=3.0, P0=P0, T_ramp=1.0, ramp="Cinf")
            make_reference.write_reference(os.path.join(cls.root, cfg), p, 256, 0.0025, 2.0)
        cls.ref = {c: json.load(open(os.path.join(cls.root, c, "reference.json"))) for c in ("combined", "torsion")}

    def arr(self, cfg, k):
        r = self.ref[cfg]
        t = np.array(r["t"])
        s = (t > 2 * 0.0025 + 1e-9) & (t <= 2.0 + 1e-9)
        return np.array(r[k])[s]

    def recompute(self, key, mixed):
        N = THR["numerical_noise"]
        noise = N["noise_E_by_qoi"]["value"]["combined"]
        q = {"d_T": ("T_B", "T_B"), "d_sig_rt": ("sig_rt_b", "sig_rt"), "d_alpha": ("alpha_Rf", "alpha_B"), "d_Hf": ("H_f", "H_f")}[key]
        C, T = self.arr("combined", q[0]), self.arr("torsion", q[0])
        signal = float(np.max(np.abs(C - T)))
        mx = float(np.max(np.abs(C)))
        return noise[q[1]] * mx / signal, mx, signal

    def test_max_torque_is_about_454_not_155(self):
        mxC = float(np.max(np.abs(self.arr("combined", "T_B"))))
        mxT = float(np.max(np.abs(self.arr("torsion", "T_B"))))
        self.assertAlmostEqual(mxC, 453.98, delta=1.0)
        self.assertAlmostEqual(mxT, 443.28, delta=1.0)
        self.assertGreater(mxC, 400.0)

    def test_documented_allowances_match_recomputation(self):
        N = THR["numerical_noise"]["mixed_term_allowance"]
        for key, doc in N["documented_values"].items():
            a, mx, sig = self.recompute(key, True)
            self.assertAlmostEqual(a / doc, 1.0, delta=0.02, msg=f"{key}: recomputed {a:.4g} against documented {doc}")

    def test_stated_magnitudes_and_units_are_consistent(self):
        M = THR["numerical_noise"]["mixed_term_allowance"]["magnitudes_from_reference_json"]
        a, mx, sig = self.recompute("d_T", True)
        self.assertAlmostEqual(mx, M["max_abs_T_B_combined"], delta=0.5)
        self.assertAlmostEqual(sig, M["signal_d_T"], delta=0.1)
        a, mx, sig = self.recompute("d_sig_rt", True)
        self.assertAlmostEqual(mx, M["max_abs_sig_rt_combined"], delta=0.5)
        self.assertAlmostEqual(sig, M["signal_d_sig_rt"], delta=0.1)
        self.assertIn("per radian per metre of depth", M["units"])

    def test_original_wrong_normalisation_would_be_caught(self):
        doc = THR["numerical_noise"]["mixed_term_allowance"]["documented_values"]["d_T"]
        wrong = 1.3e-3 * 155.0 / 74.43          # the original mistake: max|T_B| of about 155
        self.assertGreater(abs(wrong - doc) / doc, 0.5)
        self.assertGreater(doc, 7.5e-3)
        self.assertLess(doc, 8.5e-3)



class NoisePropagation(unittest.TestCase):
    def test_interval_matches_audit_for_torque(self):
        iv = CV.order_interval([2.836e-2, 7.242e-3, 1.894e-3], 2.4e-4)
        self.assertAlmostEqual(iv["p23"][2], 1.75, delta=0.02)
        self.assertAlmostEqual(iv["p23"][3], 2.12, delta=0.02)
        self.assertLessEqual(iv["p23"][4], iv["p23"][2])        # worst case contains the first-order interval
        self.assertGreaterEqual(iv["p23"][5], iv["p23"][3])

    def test_noise_not_below_error_is_unbounded(self):
        self.assertIsNone(CV.order_interval([1e-2, 1e-3, 2e-4], 3e-4)["p23"])
        self.assertIsNone(CV.self_order_interval(1e-3, 1e-4, 1e-4))

    def test_allowance_never_creates_a_pass(self):
        rnd = random.Random(1)
        for _ in range(500):
            lo = rnd.uniform(0, 3e-3)
            hi = lo + rnd.uniform(1e-3, 5e-3)
            v = rnd.uniform(0, 1.5 * hi)
            n1, n2 = rnd.uniform(0, 2e-3), rnd.uniform(0, 2e-3)
            a = ST.severity(verify.band_status(v, lo, hi, 0.2))
            b = ST.severity(verify.band_status(v, lo, hi, 0.2, n1, n2))
            self.assertGreaterEqual(b, a)

    def test_delta_per_quantity(self):
        CV.NOISE["torsion"] = dict(THR["numerical_noise"]["noise_E_by_qoi"]["value"]["torsion"])
        self.assertEqual(CV.delta_for("torsion", "T_B"), 2.4e-4)
        self.assertLess(CV.delta_for("torsion", "alpha_B"), 1e-5)


class SharedStatus(unittest.TestCase):
    def test_one_ordered_definition(self):
        self.assertIs(ACC.SEV, ST.ORDER)
        self.assertIs(FAM.SEVERITY, ST.ORDER)
        o = ST.ORDER
        self.assertLess(o["PASS"], o["CONDITIONAL"])
        self.assertLess(o["CONDITIONAL"], o["UNKNOWN"])
        self.assertLess(o["UNKNOWN"], o["INVALID-INPUT"])
        self.assertLess(o["INVALID-INPUT"], o["INCOMPLETE"])
        self.assertLess(o["INCOMPLETE"], o["FAIL"])
        self.assertEqual(o["DISCLOSED"], o["CONDITIONAL"])
        self.assertEqual(o["ABOVE-FLOOR"], o["CONDITIONAL"])

    def test_unknown_status_is_safe_and_not_pass(self):
        self.assertEqual(ST.norm("weird"), "UNKNOWN")
        self.assertEqual(ST.severity("weird"), ST.ORDER["UNKNOWN"])
        self.assertEqual(ST.worst(["PASS", "weird"]), "UNKNOWN")
        self.assertEqual(ACC.worst(["PASS", "nonsense", "CONDITIONAL"]), "UNKNOWN")
        self.assertEqual(verify.worst(["PASS", None]), "UNKNOWN")
        self.assertIn("UNKNOWN", ST.BLOCKING)
        self.assertEqual(ST.worst([]), "NOT-VERIFIED")           # nothing evaluated: never PASS (behaviour documented in status.py)
        self.assertEqual(ST.worst(["PASS", "weird", "FAIL"]), "FAIL")

    def test_cap_pass(self):
        self.assertEqual(ST.cap_pass("PASS", "CONDITIONAL"), "CONDITIONAL")
        self.assertEqual(ST.cap_pass("FAIL", "CONDITIONAL"), "FAIL")
        self.assertEqual(ST.cap_pass("PASS", "PASS"), "PASS")


def overlay_case(src, dst):
    """Copy the small dictionaries of an accepted case, symlink the logs, diagnostics and reference."""
    os.makedirs(dst)
    for top in ("constant", "system", "0"):
        for dp, dn, fn in os.walk(os.path.join(src, top)):
            dn[:] = [d for d in dn if d not in ("polyMesh",)]
            for f in fn:
                if f.endswith(".withDefaultValues"):
                    continue
                rel = os.path.relpath(os.path.join(dp, f), src)
                os.makedirs(os.path.dirname(os.path.join(dst, rel)), exist_ok=True)
                shutil.copy2(os.path.join(dp, f), os.path.join(dst, rel))
    for f in ("log.solids4Foam", "rtDiagFluid.dat", "rtDiagSolid.dat", "reference.json"):
        if os.path.exists(os.path.join(src, f)):
            os.symlink(os.path.join(src, f), os.path.join(dst, f))


class FamilyInjection(unittest.TestCase):
    """Failure injection: each change of a documented-benchmark setting must give INVALID-INPUT (never PASS), in M2 only and uniformly."""

    @classmethod
    def setUpClass(cls):
        cls.root = tempfile.mkdtemp(dir=os.environ.get("TMPDIR"))
        cls.base = os.path.join(cls.root, "base")
        make_family("torsion", cls.base, absTol=1e-12)

    def mutated(self, edits, meshes=("M2",)):
        d = tempfile.mkdtemp(dir=self.root)
        shutil.copytree(self.base, d, dirs_exist_ok=True)
        for m in meshes:
            for rel, old, new in edits:
                p = os.path.join(d, f"torsion_{m}", *rel.split("/"))
                with open(p) as fh:
                    txt = fh.read()
                self.assertIn(old, txt, rel)
                with open(p, "w") as fh:
                    fh.write(txt.replace(old, new, 1))
        return FAM.check_family(d, "torsion", THR)

    CASES = [
        ("aitken", [("constant/fsiProperties", "fluidSolidInterface    IQNILS;", "fluidSolidInterface    Aitken;")]),
        ("transfer", [("constant/fsiProperties", "interfaceTransferMethod directMap;", "interfaceTransferMethod GGI;")]),
        ("fluid_predictor", [("constant/fsiProperties", "predictor            yes;", "predictor            no;")]),
        ("solid_predictor_switch", [("constant/fsiProperties", "predictSolid         yes;", "predictSolid         no;")]),
        ("relaxation", [("constant/fsiProperties", "relaxationFactor     0.05;", "relaxationFactor     0.5;")]),
        ("reuse", [("constant/fsiProperties", "couplingReuse        5;", "couplingReuse        0;")]),
        ("scaleFactor", [("constant/solid/solidProperties", "scaleFactor 0.5;", "scaleFactor 0.1;")]),
        ("stabilisation_type", [("constant/solid/solidProperties", "diffStencilLaplacian", "Laplacian")]),
        ("solid_model", [("constant/solid/solidProperties", "nonLinearGeometryTotalLagrangianTotalDisplacement;", "linearGeometryTotalDisplacement;")]),
        ("solution_algorithm", [("constant/solid/solidProperties", "solutionAlgorithm PETScSNES;", "solutionAlgorithm newtonRaphson;")]),
        ("law", [("constant/solid/mechanicalProperties", "neoHookeanElastic", "linearElastic")]),
        ("E", [("constant/solid/mechanicalProperties", "20000.0", "30000.0")]),
        ("nu", [("constant/solid/mechanicalProperties", "0.3;", "0.4;")]),
        ("momentum_predictor", [("system/fluid/fvSolution", "momentumPredictor   on;", "momentumPredictor   off;")]),
        ("flux_consistent", [("system/fluid/fvSolution", "fluxConsistentPatches (outer);", "fluxConsistentPatches ();")]),
        ("ncorrectors", [("system/fluid/fvSolution", "nCorrectors         3;", "nCorrectors         2;")]),
        ("diffusivity", [("constant/fluid/dynamicMeshDict", "quadratic inverseDistance (interface);\n}", "uniform;\n}")]),
        ("motion_solver", [("constant/fluid/dynamicMeshDict", "motionSolver    displacementLaplacian;", "motionSolver    velocityLaplacian;")]),
        ("U_bc_type", [("0/fluid/U", "type            codedMixed;", "type            fixedValue;")]),
        ("U_bc_code", [("0/fluid/U", "(Up - 2.0*nf*(nf & Up))/rf", "(Up - 1.0*nf*(nf & Up))/rf")]),
        ("p_bc_code", [("0/fluid/p", "const scalar seed = 0.0;", "const scalar seed = 5.0;")]),
        ("pointDisplacement_bc", [("0/fluid/pointDisplacement", "type            fixedValue;", "type            zeroGradient;")]),
        ("D_bc", [("0/solid/D", "type            fixedRotation;", "type            fixedDisplacement;")]),
        ("dt", [("system/controlDict", "deltaT          0.0025;", "deltaT          0.005;")]),
        ("fluid_scheme", [("system/fluid/fvSchemes", "div(phi,U)      Gauss linear;", "div(phi,U)      Gauss upwind;")]),
        ("solid_linear_solver", [("system/solid/fvSolution", "ksp_type lgmres;", "ksp_type gmres;")]),
        ("fluid_p_solver", [("system/fluid/fvSolution", "tolerance       1e-12;", "tolerance       1e-10;")]),
    ]

    def test_each_change_is_invalid_in_m2(self):
        for name, edits in self.CASES:
            with self.subTest(change=name):
                try:
                    r = self.mutated(edits)
                except AssertionError as e:
                    self.fail(f"edit target not found for {name}: {e}")
                self.assertEqual(r["status"], "INVALID-INPUT", (name, [f for f in r["findings"] if f["status"] != "PASS"]))

    def test_uniform_departure_from_benchmark_is_invalid(self):
        for name in ("aitken", "scaleFactor", "momentum_predictor", "law"):
            edits = dict(self.CASES)[name]
            with self.subTest(change=name):
                r = self.mutated(edits, meshes=("M1", "M2", "M3"))
                self.assertEqual(r["status"], "INVALID-INPUT")
                self.assertTrue(any("documented benchmark" in f["detail"] or "documented" in f["detail"] for f in r["findings"] if f["status"] == "INVALID-INPUT"))

    def test_nnonorth_absTol_mesh_theta_P0_invalid(self):
        for name, mods in (("nNonOrth", ["--nNonOrth", "0"]), ("absTol", ["--absTol", "1e-9"]), ("mesh_ratio", ["--ntheta", "96"]),
                           ("theta0", ["--theta0", "0.05"]), ("P0", ["--P0", "-1000"])):
            with self.subTest(change=name):
                d = tempfile.mkdtemp(dir=self.root)
                make_family("torsion", d, absTol=1e-12, mods={"M2": mods})
                self.assertEqual(FAM.check_family(d, "torsion", THR)["status"], "INVALID-INPUT")

    def test_documented_equivalent_differences_stay_conditional(self):
        # older seeded-pressure code with seed 0 on M2, absTol pattern: CONDITIONAL, not INVALID
        d = tempfile.mkdtemp(dir=self.root)
        shutil.copytree(self.base, d, dirs_exist_ok=True)
        p = os.path.join(d, "torsion_M2", "0", "fluid", "p")
        txt = open(p).read()
        i = txt.index("scalar pert = 0;")
        j = txt.index("pv[i] = pc + seed*bump*pert;")
        legacy = txt[:i] + "pv[i] = pc + seed*bump*Foam::cos(2*th);" + txt[j + len("pv[i] = pc + seed*bump*pert;"):]
        legacy = legacy.replace("const labelList modes({2});", "")
        with open(p, "w") as fh:
            fh.write(legacy)
        r = FAM.check_family(d, "torsion", THR)
        self.assertEqual(r["status"], "CONDITIONAL", [f for f in r["findings"] if f["status"] != "PASS"])


class MalformedInputs(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.root = tempfile.mkdtemp(dir=os.environ.get("TMPDIR"))
        cls.base = os.path.join(cls.root, "base")
        make_family("torsion", cls.base, absTol=1e-12)

    def fresh(self):
        d = tempfile.mkdtemp(dir=self.root)
        shutil.copytree(self.base, d, dirs_exist_ok=True)
        return d

    def status(self, d):
        r = FAM.check_family(d, "torsion", THR)
        self.assertNotEqual(r["status"], "PASS")
        self.assertIn(r["status"], ("INCOMPLETE", "INVALID-INPUT"), r["findings"])
        return r

    def test_missing_dictionary(self):
        d = self.fresh()
        os.remove(os.path.join(d, "torsion_M3", "constant", "solid", "mechanicalProperties"))
        self.assertEqual(self.status(d)["status"], "INCOMPLETE")

    def test_empty_file(self):
        d = self.fresh()
        open(os.path.join(d, "torsion_M1", "system", "fluid", "fvSolution"), "w").close()
        self.assertEqual(self.status(d)["status"], "INCOMPLETE")

    def test_binary_garbage(self):
        d = self.fresh()
        with open(os.path.join(d, "torsion_M2", "constant", "fsiProperties"), "wb") as fh:
            fh.write(bytes(range(256)) * 4)
        self.assertEqual(self.status(d)["status"], "INCOMPLETE")

    def test_garbled_entry(self):
        d = self.fresh()
        p = os.path.join(d, "torsion_M2", "constant", "solid", "mechanicalProperties")
        txt = open(p).read().replace("20000.0", "twenty thousand")
        open(p, "w").write(txt)
        self.assertEqual(self.status(d)["status"], "INCOMPLETE")

    def test_unbalanced_braces(self):
        d = self.fresh()
        p = os.path.join(d, "torsion_M2", "system", "fluid", "fvSolution")
        txt = open(p).read().replace("PIMPLE\n{", "PIMPLE\n")
        open(p, "w").write(txt)
        self.status(d)

    def test_truncated_drive_table(self):
        d = self.fresh()
        p = os.path.join(d, "torsion_M3", "constant", "timeVsAngle")
        lines = open(p).read().splitlines()
        open(p, "w").write("\n".join(lines[:300]) + "\n)\n")
        self.status(d)

    def test_missing_mesh_dict_and_missing_run(self):
        d = self.fresh()
        os.remove(os.path.join(d, "torsion_M1", "system", "solid", "blockMeshDict"))
        self.assertEqual(self.status(d)["status"], "INCOMPLETE")
        d2 = self.fresh()
        shutil.rmtree(os.path.join(d2, "torsion_M2"))
        self.assertEqual(self.status(d2)["status"], "INCOMPLETE")

    def test_driver_never_tracebacks_and_never_passes(self):
        d = self.fresh()
        open(os.path.join(d, "torsion_M2", "constant", "fsiProperties"), "w").write("garbled {{{ ;;; \x00")
        os.remove(os.path.join(d, "torsion_M1", "constant", "pcTable"))
        for m in ("M1", "M2", "M3"):
            with open(os.path.join(d, f"torsion_{m}", "log.solids4Foam"), "w") as fh:
                fh.write(make_log(accept="both", n=60))
            with open(os.path.join(d, f"torsion_{m}", "rtDiagSolid.dat"), "w") as fh:
                fh.write("1 2 3\n")                       # truncated, malformed diagnostics
            with open(os.path.join(d, f"torsion_{m}", "rtDiagFluid.dat"), "w") as fh:
                fh.write("")
        out = subprocess.run([sys.executable, "-B", "-I", os.path.join(HERE, "scripts", "verify.py"), d, "--out", os.path.join(d, "o")],
                             capture_output=True, text=True)
        self.assertNotIn("Traceback", out.stderr + out.stdout)
        self.assertEqual(out.returncode, 1)
        rep = json.load(open(os.path.join(d, "o", "verification_report.json")))
        self.assertNotIn("PASS", [r["status"] for r in rep["layerA"].values()])
        self.assertNotEqual(rep["layerB"]["claims"]["torsion_primary_orders"]["status"], "PASS")

    def test_truncated_diagnostics_series_is_not_pass(self):
        d = tempfile.mkdtemp(dir=self.root)
        with open(os.path.join(d, "log.solids4Foam"), "w") as fh:
            fh.write(make_log(accept="both", n=60))
        rows = "\n".join(" ".join(["%g" % (0.01 * (i + 1))] + ["0.0"] * 22) for i in range(60))
        with open(os.path.join(d, "rtDiagSolid.dat"), "w") as fh:
            fh.write(rows + "\n")
        with open(os.path.join(d, "rtDiagFluid.dat"), "w") as fh:
            fh.write("\n".join(rows.splitlines()[:30]) + "\n")
        r = ACC.run_layer_a(d, THR, cfg_override=dict(endTime=0.6, deltaT=0.01, nNonOrth=2, tol=1e-6, absTol=1e-12, nOuterCorr=100))
        self.assertNotEqual(r["status"], "PASS")
        self.assertIn(r["checks"]["diagnostics_rows"]["status"], ("FAIL", "INCOMPLETE"))

    def test_garbled_log_is_incomplete_not_traceback(self):
        d = tempfile.mkdtemp(dir=self.root)
        with open(os.path.join(d, "log.solids4Foam"), "wb") as fh:
            fh.write(bytes(range(256)) * 50)
        r = ACC.run_layer_a(d, THR)
        self.assertIn(r["status"], ("INCOMPLETE", "FAIL"))
        d2 = tempfile.mkdtemp(dir=self.root)
        open(os.path.join(d2, "log.solids4Foam"), "w").close()
        self.assertEqual(ACC.run_layer_a(d2, THR)["status"], "INCOMPLETE")


class FloorByQuantity(unittest.TestCase):
    """By-quantity floor statuses on real accepted data (needs RT_KNOWN_RUNS)."""

    @unittest.skipUnless(os.environ.get("RT_KNOWN_RUNS"), "RT_KNOWN_RUNS not set")
    def test_floor_statuses(self):
        root = os.environ["RT_KNOWN_RUNS"]
        W = tempfile.mkdtemp(dir=os.environ.get("TMPDIR"))
        names = {"torsion_M1": "F2_M1T_nNonOrth2", "torsion_M2": "F2_M2T_nNonOrth2", "torsion_M3": "F1_M3T_nNonOrth2_clean",
                 "combined_M1": "F3_M1C_nNonOrth2", "combined_M2": "F3_M2C_nNonOrth2", "combined_M3": "F3_M3C_nNonOrth2"}
        for n, c in names.items():
            overlay_case(os.path.join(root, c), os.path.join(W, n))
        runs = {tuple(n.split("_")): ACC.run_layer_a(os.path.join(W, n), THR) for n in names}
        rep = verify.layer_b(W, runs, THR)
        fq = rep["floor_by_quantity"]
        for q in ("alpha_B", "T_B", "sig_rt", "T_A", "H_f"):
            self.assertTrue(fq[q]["in_floor_statement"])
            self.assertNotEqual(fq[q]["status"], "PASS")
        for q, e3 in (("H_s", 6.9e-3), ("w_wall", 6.3e-3), ("KE_s", 1.06e-2), ("P_wall", 4.8e-2)):
            self.assertFalse(fq[q]["in_floor_statement"])
            self.assertEqual(fq[q]["status"], "ABOVE-FLOOR")
            self.assertAlmostEqual(fq[q]["E_M3"] / e3, 1.0, delta=0.02)
        self.assertAlmostEqual(min(fq[q]["E_M3"] for q in ("alpha_B", "T_B", "sig_rt", "T_A", "H_f")), 3.06e-3, delta=2e-5)
        self.assertNotEqual(rep["claims"]["combined_reference_error_floor"]["status"], "PASS")
        self.assertNotEqual(rep["claims"]["mixed_terms_M3"]["status"], "PASS")
        text = rep["claims"]["mixed_terms_M3"]["detail"]
        self.assertIn("d_T 2.47% (limit 3.0%, allowance 0.793%", text)
        self.assertIn("CONDITIONAL", text.split("0-2 d_T")[1].split(";")[0])
        # accepted copies with an injected setting: INVALID-INPUT and layer B stops
        p = os.path.join(W, "combined_M2", "constant", "fsiProperties")
        with open(p) as fh:
            txt = fh.read()
        with open(p, "w") as fh:
            fh.write(txt.replace("fluidSolidInterface    IQNILS;", "fluidSolidInterface    Aitken;"))
        rep2 = verify.layer_b(W, runs, THR)
        self.assertEqual(rep2["families"]["combined"]["status"], "INVALID-INPUT")
        self.assertEqual(rep2["claims"]["mixed_terms_M3"]["status"], "INVALID-INPUT")
        # accepted torsion family itself is CONDITIONAL, never INVALID
        self.assertEqual(FAM.check_family(W, "torsion", THR)["status"], "CONDITIONAL")


class ReadmeConsistency(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        with open(os.path.join(TUT, "README.md")) as fh:
            cls.readme = fh.read()
        with open(os.path.join(TUT, "verification", "README.md")) as fh:
            cls.vreadme = fh.read()

    def test_what_is_checked_is_stated(self):
        for key in ("coupling algorithm", "interfaceTransferMethod", "momentumPredictor", "fluxConsistentPatches", "dynamic mesh motion",
                    "boundary conditions", "fingerprint", "Not checked", "never a crash"):
            self.assertIn(key, self.readme, key)
        doc = FAM.__doc__
        for key in ("WHAT IS CHECKED", "NOT checked here", "coupling", "boundary", "safety net"):
            self.assertIn(key, doc)

    def test_stale_statements_removed(self):
        for text in (self.readme, self.vreadme):
            self.assertNotIn("85-93", text)
            self.assertNotIn("Local preparation", text)
            self.assertNotIn("IQN-ILS coupling diverges", text)
            self.assertNotIn("A stricter test was not", text)
        self.assertNotIn("99 s (98.56", self.readme)
        self.assertIn("98.6 s", self.readme)

    def test_floor_wording_by_quantity(self):
        fb = THR["layerB_deterministic"]["combined_floor_band_M3"]
        self.assertEqual(fb["qoi"], ["alpha_B", "T_B", "sig_rt", "T_A", "H_f"])
        for key in ("6.9e-3", "6.3e-3", "1.06e-2", "4.8e-2", "3.06e-3", "lower edge"):
            self.assertIn(key, self.readme, key)
        self.assertIn("ONLY", " ".join(open(os.path.join(HERE, "scripts", "verify.py")).read().split()))

    def test_noise_statement_consistent(self):
        N = THR["numerical_noise"]
        self.assertEqual(N["torsion_max_norm_noise_of_E"]["value"], 2.4e-4)
        self.assertEqual(N["combined_max_norm_noise_of_E"]["value"], 1.3e-3)
        self.assertIn("2.4e-4", self.readme)
        self.assertIn("1.3e-3", self.readme)
        self.assertIn("combined-run value, not the torsion one", " ".join(self.readme.split()))
        self.assertIn("2.4e-4", N["_comment"])
        self.assertIn("1.3e-3", N["_comment"])

    def test_readme_tables_match_computation(self):
        iv = CV.order_interval([2.836e-2, 7.242e-3, 1.894e-3], 2.4e-4)["p23"]
        self.assertIn("1.94 (%.2f-%.2f noise-limited)" % (iv[2], iv[3]), self.readme)
        self.assertIn("97.5%", self.readme)
        self.assertIn("93.0-97.5%", self.readme)

    def test_no_work_paths_outside_tutorial(self):
        self.assertNotRegex(self.readme, r"work/(F[0-9]|G[A-Z]|P[0-9]|runs)")
        self.assertNotRegex(self.vreadme, r"work/(F[0-9]|G[A-Z]|P[0-9]|runs)")
        self.assertNotRegex(json.dumps(THR), r"work/|parallel_vs_serial|[A-Z][A-Z0-9]?_audit|PROVENANCE|Phase G|F[123]_[A-Z0-9]")


class NoiseDriverPath(unittest.TestCase):
    """The mixed-term allowance is exercised through the DRIVER code (verify.mixed_claim, convergence.mixed_terms) on a small
    generated fixture (reference solutions, no CFD, no RT_KNOWN_RUNS): a wrong normalisation inside the driver is caught here."""

    @classmethod
    def setUpClass(cls):
        root = tempfile.mkdtemp(dir=os.environ.get("TMPDIR"))
        cls.series = {}
        for cfg, th, P0 in (("combined", 0.1, -3000.0), ("torsion", 0.1, 0.0)):
            p = rt.Params(theta0=th, Omega=3.0, P0=P0, T_ramp=1.0, ramp="Cinf")
            make_reference.write_reference(os.path.join(root, cfg), p, 256, 0.0025, 2.0)
            r = json.load(open(os.path.join(root, cfg, "reference.json")))
            t = np.array(r["t"])
            ref = {"alpha_B": np.array(r["alpha_Rf"]), "T_B": np.array(r["T_B"]), "sig_rt": np.array(r["sig_rt_b"]), "H_f": np.array(r["H_f"])}
            num = {k: v * (1.0 + 0.01 * np.sin(7.0 * t)) for k, v in ref.items()}      # a deterministic 1% perturbation
            cls.series[cfg] = dict(t=t, num=num, ref=ref, dt=0.0025)

    def mixed(self):
        CV.OVERRIDE.clear()
        for cfg, S in self.series.items():
            CV.OVERRIDE[(cfg, "M3")] = S
        try:
            return {"M3": CV.mixed_terms("/nonexistent", "M3")}
        finally:
            CV.OVERRIDE.clear()

    def test_claim_allowances_follow_the_formula(self):
        st, detail = verify.mixed_claim(self.mixed(), THR)
        self.assertIn("0-2 d_T", detail)
        self.assertRegex(detail, r"d_T [0-9.]+% \(limit 3.0%, allowance 0\.793% = 0\.0013 x 454 / 74\.4")
        self.assertRegex(detail, r"d_sig_rt [0-9.]+% \(limit 3.0%, allowance 0\.56[01]% = 0\.00121 x 418\.4 / 90\.3")

    def test_wrong_normalisation_inside_the_driver_is_caught(self):
        good = verify.mixed_claim(self.mixed(), THR)[1]
        orig = verify.mixed_allowance
        try:
            verify.mixed_allowance = lambda noise, max_ref, signal: noise * 155.0 / signal      # the original mistake
            bad = verify.mixed_claim(self.mixed(), THR)[1]
        finally:
            verify.mixed_allowance = orig
        self.assertNotEqual(good, bad)
        self.assertNotIn("allowance 0.793%", bad)
        self.assertIn("allowance 0.793%", good)
        self.assertAlmostEqual(verify.mixed_allowance(1.3e-3, 453.98, 74.427), 7.93e-3, delta=5e-5)

    def test_allowance_never_upgrades_a_claim(self):
        m = self.mixed()
        st, _ = verify.mixed_claim(m, THR)
        for a_ in (0.0, 1e-3, 1e-2):
            for key in ("0-2", "0-1.8"):
                for term in m["M3"][key]:
                    v = m["M3"][key][term]["error_over_signal"]
                    lim = {"d_alpha": 0.025, "d_T": 0.03, "d_sig_rt": 0.03}.get(term)
                    if lim and v is not None:
                        self.assertGreaterEqual(ST.severity(verify.band_status(v, 0.0, lim, 0.2, 0.0, a_)), ST.severity(verify.band_status(v, 0.0, lim, 0.2, 0.0, 0.0)))


class IntervalLabels(unittest.TestCase):
    def test_unresolved_interval_is_a_downgrade_with_a_label(self):
        self.assertEqual(verify.interval_status(None, 1.65, 2.10), "CONDITIONAL")
        self.assertIn("unresolved (noise >= error)", verify.interval_text(None, 0.15))
        self.assertTrue(verify.noise_limited(None, 0.15))
        self.assertIsNone(CV.order_interval([1e-2, 1e-3, 2e-4], 5e-4)["p23"])

    def test_noise_limited_and_sharp_labels(self):
        self.assertIn("noise-limited", verify.interval_text([2.0, 0.6, 1.39, 2.60], 0.15))
        self.assertTrue(verify.noise_limited([2.0, 0.6, 1.39, 2.60], 0.15))
        self.assertFalse(verify.noise_limited([1.97, 0.01, 1.96, 1.97], 0.15))
        self.assertIn("sharp", verify.interval_text([1.97, 0.01, 1.96, 1.98], 0.15))
        self.assertEqual(verify.interval_status([1.9, 0.1, 1.8, 2.0], 1.65, 2.10), "PASS")      # neutral: never an upgrade by itself
        self.assertEqual(verify.interval_status([2.0, 0.6, 1.39, 2.60], 1.65, 2.10), "CONDITIONAL")

    def test_documented_bound(self):
        self.assertEqual(THR["numerical_noise"]["order_uncertainty"]["noise_limited_halfwidth"]["value"], 0.15)

    def test_combined_self_orders_are_labelled_from_real_noise(self):
        # the documented self-order intervals of the combined run: alpha_B and H_f sharp, T_B, sig_rt, T_A noise-limited
        d = {"alpha_B": (7.61e-3, 1.95e-3, 5.7e-6), "H_f": (1.41e-2, 3.68e-3, 1.11e-5), "T_B": (1.81e-2, 4.53e-3, 1.3e-3),
             "sig_rt": (1.74e-2, 4.46e-3, 1.21e-3), "T_A": (1.71e-2, 4.03e-3, 8.9e-4)}
        lab = {q: verify.noise_limited(CV.self_order_interval(a, b, n), 0.15) for q, (a, b, n) in d.items()}
        self.assertEqual(lab, {"alpha_B": False, "H_f": False, "T_B": True, "sig_rt": True, "T_A": True})


def run_driver(args, cwd=None):
    return subprocess.run([sys.executable, "-B", "-I", os.path.join(HERE, "scripts", "verify.py")] + args, capture_output=True, text=True, cwd=cwd)


class LogInputs(unittest.TestCase):
    """--log with case directories, relative and absolute log file paths, missing and malformed logs (P2)."""

    @classmethod
    def setUpClass(cls):
        cls.root = tempfile.mkdtemp(dir=os.environ.get("TMPDIR"))
        cls.work = os.path.join(cls.root, "work")
        os.makedirs(cls.work)
        cls.case = os.path.join(cls.root, "caseA")
        os.makedirs(cls.case)
        with open(os.path.join(cls.case, "log.solids4Foam"), "w") as fh:
            fh.write(make_log(accept="both", n=60))

    def go(self, log_arg, cwd=None, out=None):
        out = out or os.path.join(self.root, "o_" + str(abs(hash(log_arg)) % 10 ** 6))
        r = run_driver([self.work, "--log", log_arg, "--out", out], cwd=cwd)
        self.assertNotIn("Traceback", r.stdout + r.stderr, r.stderr)
        return r, out

    def files_under_out(self, out):
        return sorted(os.listdir(out))

    def test_case_directory(self):
        r, out = self.go(self.case)
        self.assertEqual(r.returncode, 1)          # 0.6 s of a 2 s run: INCOMPLETE by design of the synthetic log
        self.assertTrue(any(f.startswith("layerA_steps_") for f in self.files_under_out(out)))

    def test_absolute_log_path(self):
        r, out = self.go(os.path.join(self.case, "log.solids4Foam"))
        rep = json.load(open(os.path.join(out, "verification_report.json")))
        self.assertEqual(list(rep["layerA_extra_logs"]), ["caseA_log.solids4Foam"])
        self.assertFalse(any("/" in f for f in self.files_under_out(out)))
        self.assertTrue(os.path.isfile(os.path.join(out, "layerA_steps_caseA_log.solids4Foam.csv")))

    def test_relative_log_path(self):
        r, out = self.go("caseA/log.solids4Foam", cwd=self.root, out=os.path.join(self.root, "o_rel"))
        rep = json.load(open(os.path.join(out, "verification_report.json")))
        self.assertIn("caseA_log.solids4Foam", rep["layerA_extra_logs"])
        r2, out2 = self.go("../caseA/log.solids4Foam", cwd=self.work, out=os.path.join(self.root, "o_rel2"))
        self.assertTrue(os.path.isfile(os.path.join(out2, "layerA_steps_caseA_log.solids4Foam.csv")))

    def test_output_stays_under_the_output_directory(self):
        before = set(os.listdir(self.root))
        r, out = self.go(os.path.join(self.case, "log.solids4Foam"), out=os.path.join(self.root, "o_clean"))
        after = set(os.listdir(self.root))
        self.assertEqual(after - before, {"o_clean"})

    def test_missing_log_is_a_clean_incomplete(self):
        r, out = self.go(os.path.join(self.root, "does", "not", "exist.log"))
        self.assertEqual(r.returncode, 1)
        rep = json.load(open(os.path.join(out, "verification_report.json")))
        st = [v["status"] for v in rep["layerA_extra_logs"].values()]
        self.assertEqual(st, ["INCOMPLETE"])
        r2, out2 = self.go("no_such_dir/log.solids4Foam", cwd=self.root)
        self.assertEqual(r2.returncode, 1)

    def test_malformed_log(self):
        bad = os.path.join(self.root, "bad.log")
        with open(bad, "wb") as fh:
            fh.write(bytes(range(256)) * 80)
        r, out = self.go(bad)
        rep = json.load(open(os.path.join(out, "verification_report.json")))
        self.assertIn(list(rep["layerA_extra_logs"].values())[0]["status"], ("INCOMPLETE", "FAIL"))
        self.assertEqual(r.returncode, 1)
        empty = os.path.join(self.root, "empty.log")
        open(empty, "w").close()
        r2, out2 = self.go(empty)
        self.assertEqual(r2.returncode, 1)

    def test_run_label_is_sanitised_and_unique(self):
        lab = verify.run_label("/a/b c/log.solids4Foam")
        self.assertNotIn("/", lab)
        self.assertEqual(lab, "b_c_log.solids4Foam")
        self.assertNotEqual(verify.run_label("/x/y/log.a", taken=["y_log"]), "y_log")
        self.assertNotIn("/", verify.run_label("../../weird name/::log"))


class StatusWorst(unittest.TestCase):
    def test_not_verified_is_never_pass(self):
        self.assertEqual(ST.worst(["NOT-VERIFIED"]), "NOT-VERIFIED")
        self.assertNotEqual(ST.worst(["NOT-VERIFIED"]), "PASS")
        self.assertEqual(ST.worst([]), "NOT-VERIFIED")
        self.assertEqual(ST.decide([]), "NOT-VERIFIED")
        self.assertEqual(ST.decide(["INFO", "NOT-VERIFIED"]), "NOT-VERIFIED")

    def test_every_combination_with_pass(self):
        for s_ in ST.ORDER:
            w = ST.worst(["PASS", s_])
            self.assertEqual(w == "PASS", s_ == "PASS", (s_, w))
            self.assertEqual(ST.worst([s_, "PASS"]), w)
        for s_ in ("weird", None, 3, ""):
            self.assertEqual(ST.worst(["PASS", s_]), "UNKNOWN")
        for a_ in ST.ORDER:
            for b_ in ST.ORDER:
                self.assertEqual(ST.severity(ST.worst([a_, b_])), max(ST.severity(a_), ST.severity(b_)))

    def test_decide_ignores_neutral_items_but_not_failures(self):
        self.assertEqual(ST.decide(["PASS", "INFO", "NOT-VERIFIED"]), "PASS")
        self.assertEqual(ST.decide(["PASS", "NOT-VERIFIED", "CONDITIONAL"]), "CONDITIONAL")
        self.assertEqual(ST.decide(["PASS", "weird"]), "UNKNOWN")


class RejectionFractions(unittest.TestCase):
    """(A) fail the both-measures test, (B) strict-mode reject: recomputed from a synthetic log with known counts."""

    def test_known_counts(self):
        per = {}
        for i in range(1, 61):
            if i <= 10:
                continue                                   # 10 steps: both (default accept='both')
            elif i <= 20:
                per[i] = (5e-7, 5e-3, 1e-3)                # 10 r1-only, absolute residual high
            elif i <= 30:
                per[i] = (5e-3, 5e-7, 1e-3)                # 10 r2-only, absolute high
            elif i <= 35:
                per[i] = (5e-3, 5e-7, 5e-13)               # 5 r2-only rescued by the absolute residual
            elif i <= 45:
                per[i] = (3e-3, 4e-3, 5e-13, True)         # 10 absolute-tolerance steps
            else:
                per[i] = (5e-7, 5e-7, 1e-3)                # 15 more 'both'
        r = run_a(make_log(accept="both", per_step=per))
        rj = r["rejection"]
        self.assertEqual(r["criterion_counts"], {"both": 25, "r1-only": 10, "r2-only": 15, "absolute": 10, "none": 0})
        self.assertEqual(rj["steps"], 60)
        self.assertEqual(rj["A_fail_both_measures"], 35)                       # everything not 'both'
        self.assertEqual(rj["B_strict_mode_reject"], 20)                       # r1-only and r2-only with a high absolute residual
        self.assertEqual(rj["single_measure"], 25)
        self.assertEqual(rj["single_measure_rescued_by_absolute"], 5)
        self.assertAlmostEqual(rj["A_pct"], 100 * 35 / 60)
        self.assertAlmostEqual(rj["B_pct"], 100 * 20 / 60)
        self.assertNotAlmostEqual(rj["A_pct"], rj["B_pct"])


class FamilyGeometryAndControlDict(unittest.TestCase):
    """P3: mesh geometry, grading, patches and every controlDict entry (before and after the functions block) are semantic checks."""

    @classmethod
    def setUpClass(cls):
        cls.root = tempfile.mkdtemp(dir=os.environ.get("TMPDIR"))
        cls.base = os.path.join(cls.root, "base")
        make_family("torsion", cls.base, absTol=1e-12)

    def edit(self, mesh_edits=(), meshes=("M2",)):
        d = tempfile.mkdtemp(dir=self.root)
        shutil.copytree(self.base, d, dirs_exist_ok=True)
        for m in meshes:
            for rel, fn in mesh_edits:
                p = os.path.join(d, f"torsion_{m}", *rel.split("/"))
                with open(p) as fh:
                    txt = fh.read()
                new = fn(txt)
                self.assertNotEqual(new, txt, rel)
                with open(p, "w") as fh:
                    fh.write(new)
        return d

    @staticmethod
    def scale_radius(target, factor):
        import re as _re

        def fn(txt):
            def sub(m):
                x, y, z = (float(m.group(i)) for i in (1, 2, 3))
                if abs(math.hypot(x, y) - target) < 1e-9:
                    return "(%.15g %.15g %.15g)" % (x * factor, y * factor, z)
                return m.group(0)
            return _re.sub(r"\(([-+0-9.e]+) ([-+0-9.e]+) ([-+0-9.e]+)\)", sub, txt)
        return fn

    def status(self, d):
        return FAM.check_family(d, "torsion", THR)

    def test_outer_radius_change_is_invalid(self):
        for meshes in (("M2",), ("M1", "M2", "M3")):
            r = self.status(self.edit([("system/fluid/blockMeshDict", self.scale_radius(2.0, 1.25))], meshes))
            self.assertEqual(r["status"], "INVALID-INPUT")
            self.assertTrue(any(f["check"] in ("geometry_fluid", "mesh_geometry_fluid") and f["status"] == "INVALID-INPUT" for f in r["findings"]))

    def test_inner_radius_change_is_invalid(self):
        r = self.status(self.edit([("system/solid/blockMeshDict", self.scale_radius(0.5, 1.2))]))
        self.assertEqual(r["status"], "INVALID-INPUT")
        self.assertTrue(any(f["check"] == "geometry_solid" and f["status"] == "INVALID-INPUT" for f in r["findings"]))

    def test_interface_radius_change_is_invalid(self):
        r = self.status(self.edit([("system/fluid/blockMeshDict", self.scale_radius(1.0, 1.1)), ("system/solid/blockMeshDict", self.scale_radius(1.0, 1.1))]))
        self.assertEqual(r["status"], "INVALID-INPUT")

    def test_grading_and_patch_changes_are_invalid(self):
        r = self.status(self.edit([("system/fluid/blockMeshDict", lambda t: t.replace("simpleGrading (1 1 1)", "simpleGrading (1 2 1)", 1))]))
        self.assertTrue(any(f["check"] == "grading_fluid" and f["status"] == "INVALID-INPUT" for f in r["findings"]))
        r = self.status(self.edit([("system/solid/blockMeshDict", lambda t: t.replace("type wall;", "type patch;", 1))]))
        self.assertEqual(r["status"], "INVALID-INPUT")
        r = self.status(self.edit([("system/fluid/blockMeshDict", lambda t: t.replace("simpleGrading (1 1 1)", "simpleGrading (1 1 1)").replace("(0 1 3 2 8 9 11 10)", "(0 1 3 2 8 9 10 11)", 1))]))
        self.assertEqual(r["status"], "INVALID-INPUT")

    def test_controlDict_before_and_after_functions(self):
        before = [("deltaT          0.0025;", "deltaT          0.0025;\nadjustTimeStep  yes;"), ("writeControl    adjustableRunTime;", "writeControl    timeStep;"),
                  ("startFrom       startTime;", "startFrom       latestTime;"), ("runTimeModifiable no;", "runTimeModifiable yes;"),
                  ("writePrecision  14;", "writePrecision  6;")]
        for old, new in before:
            with self.subTest(before=old):
                r = self.status(self.edit([("system/controlDict", lambda t, o=old, n=new: t.replace(o, n, 1))]))
                self.assertEqual(r["status"], "INVALID-INPUT", [f for f in r["findings"] if f["status"] != "PASS"])
        after = [lambda t: t + "\nadjustTimeStep yes;\n", lambda t: t + "\ndeltaT 0.005;\n", lambda t: t + "\nmaxDeltaT 0.01;\n", lambda t: t + "\nwriteControl timeStep;\n"]
        for k, fn in enumerate(after):
            with self.subTest(after=k):
                r = self.status(self.edit([("system/controlDict", fn)]))
                self.assertEqual(r["status"], "INVALID-INPUT", [f for f in r["findings"] if f["status"] != "PASS"])

    def test_function_object_parameters_and_code(self):
        r = self.status(self.edit([("system/controlDict", lambda t: t.replace("executeControl  timeStep;", "executeControl  writeTime;", 1))]))
        self.assertEqual(r["status"], "INVALID-INPUT")
        r = self.status(self.edit([("system/controlDict", lambda t: t.replace("const scalar rho = 1000.0;", "const scalar rho = 999.0;", 1))]))
        self.assertEqual(r["status"], "INVALID-INPUT")
        self.assertTrue(any(f["check"] == "diagnostics_function_objects" and f["status"] == "INVALID-INPUT" for f in r["findings"]))
        r = self.status(self.edit([("system/controlDict", lambda t: t.replace("region          fluid;", "region          solid;", 1))]))
        self.assertEqual(r["status"], "INVALID-INPUT")

    def test_endtime_and_writeinterval_are_allowed_to_differ(self):
        r = self.status(self.edit([("system/controlDict", lambda t: t.replace("writeInterval   0.1;", "writeInterval   0.0125;", 1))]))
        self.assertEqual(r["status"], "PASS", [f for f in r["findings"] if f["status"] != "PASS"])

    def test_harmless_reformatting_stays_compatible(self):
        import re as _re

        def ctl(t):
            t = t.replace("deltaT          0.0025;", "")                       # move an entry after the functions block
            t = t.replace("application     solids4Foam;", "/* a block\n comment */ application   solids4Foam; // trailing")
            t = _re.sub(r"[ ]{2,}", "     ", t)
            return t + "\n// end\ndeltaT   2.5e-3 ;\n"

        def mesh(t):
            t = t.replace("convertToMeters 1;", "// units\nconvertToMeters   1.0 ;")
            t = t.replace("simpleGrading (1 1 1)", "simpleGrading ( 1.0   1.0 1 )")
            t = _re.sub(r"\n", "\n  ", t, count=40)
            return t

        def fvs(t):
            return t.replace("momentumPredictor   on;", "momentumPredictor on ; // c").replace("nCorrectors         3;", "nCorrectors 3.0;")

        d = self.edit([("system/controlDict", ctl), ("system/fluid/blockMeshDict", mesh), ("system/solid/blockMeshDict", mesh), ("system/fluid/fvSolution", fvs)], meshes=("M1", "M2", "M3"))
        r = self.status(d)
        self.assertEqual(r["status"], "PASS", [f for f in r["findings"] if f["status"] != "PASS"])

    def test_geometry_change_never_gives_a_layer_b_pass(self):
        d = self.edit([("system/fluid/blockMeshDict", self.scale_radius(2.0, 1.25))])
        r = FAM.check_family(d, "torsion", THR)
        self.assertNotIn(r["status"], ("PASS", "CONDITIONAL"))

    def test_canonical_form_ignores_order_and_comments(self):
        a = "FoamFile { version 2.0; }\nsolvers { p { tol 1e-12; solver GAMG; } }\nx 1.0;\ny (1 2 3);"
        b = "// c\ny ( 1.0  2 3 ) ;\n/* blk */ x 1 ;\nsolvers { p { solver GAMG; tol 1e-12 ; } }"
        self.assertEqual(FAM.canon(FAM.parse_dict(a)), FAM.canon(FAM.parse_dict(b)))
        self.assertNotEqual(FAM.canon(FAM.parse_dict(a)), FAM.canon(FAM.parse_dict(b.replace("GAMG", "PCG"))))
        for bad in ("x 1;;", "x { y 1; ", "x 1", "}", "x ( 1 2;"):
            with self.assertRaises(ValueError):
                FAM.parse_dict(bad)

    def test_accepted_geometry_is_documented_geometry(self):
        c = FAM.read_case(TUT)
        G = THR["family_consistency"]["documented_geometry"]["value"]
        self.assertEqual(c["mesh_solid"]["radii"], [G["A"], G["B"]])
        self.assertEqual(c["mesh_fluid"]["radii"], [G["B"], G["c"]])
        self.assertEqual(c["mesh_fluid"]["z"], [0.0, G["depth"]])


KNOWN = {   # path relative to RT_KNOWN_RUNS -> expected layer-A status
    "F2_M1T_nNonOrth2": "CONDITIONAL", "F2_M2T_nNonOrth2": "CONDITIONAL", "F1_M3T_nNonOrth2_clean": "CONDITIONAL",
    "F3_M1C_nNonOrth2": "CONDITIONAL", "F3_M2C_nNonOrth2": "CONDITIONAL", "F3_M3C_nNonOrth2": "CONDITIONAL",
    "F3B_M2C_nNonOrth2_absTol1e-10": "CONDITIONAL",
    "F3_M2R_nNonOrth2": "FAIL",
    "F1_M3T_nNonOrth2_clean_killed_wsl_t0p185": "INCOMPLETE",
    "F1_M3T_nNonOrth2_clean/log.killed_attempt2_dynamicCode": "INCOMPLETE",
    "F1_M3T_nNonOrth2_clean/log.failed_attempt3_stale_o": "FAIL",
    "E2a_correctPhi_yes": "FAIL",
    "M3_T_attempt2_IQNILS_divergence_t1.81": "FAIL",
}


@unittest.skipUnless(os.environ.get("RT_KNOWN_RUNS"), "RT_KNOWN_RUNS not set")
class KnownRuns(unittest.TestCase):
    def test_known_statuses(self):
        root = os.environ["RT_KNOWN_RUNS"]
        for rel, want in KNOWN.items():
            p = os.path.join(root, rel)
            if not os.path.exists(p):
                self.skipTest(f"missing {rel}")
            with self.subTest(run=rel):
                got = ACC.run_layer_a(p, THR)["status"]
                self.assertEqual(got, want, rel)


class NoiseAllowanceSynthetic(unittest.TestCase):
    """The documented mixed-term allowance, tested with hand-built inputs (no CFD, no accepted-run data, no RT_KNOWN_RUNS).
    d_T: deterministic limit 3.0% of the signal; allowance 0.0013 x 453.98 / 74.427 = 0.793%, so the unconditional PASS region is
    error <= 3.0% - 0.793% = 2.207%; 2.207-3.0% is inside the limit but marginal (CONDITIONAL); beyond the limit by at most the
    margin CONDITIONAL, beyond that FAIL."""
    MAXREF, SIGNAL = 453.98, 74.427
    LIM = 0.03

    def claim(self, v, term="d_T"):
        return verify.mixed_claim({"M3": {"0-2": {term: dict(error_over_signal=v, max_abs_ref_combined=self.MAXREF, signal_max_abs=self.SIGNAL)}}}, THR)[0]

    def expected(self, v, a=None):
        a = verify.mixed_allowance(1.3e-3, self.MAXREF, self.SIGNAL) if a is None else a
        if v <= self.LIM - a:
            return "PASS"
        if v <= self.LIM * (1 + THR["layerB_deterministic"]["band_conditional_margin"]["value"]):
            return "CONDITIONAL"
        return "FAIL"

    GRID = (0.0, 0.01, 0.0220, 0.0221, 0.0225, 0.025, 0.0299, 0.03, 0.0301, 0.035, 0.0359, 0.0361, 0.05)

    def test_allowance_value(self):
        self.assertAlmostEqual(verify.mixed_allowance(1.3e-3, self.MAXREF, self.SIGNAL), 7.93e-3, delta=5e-5)

    def test_both_sides_of_both_thresholds(self):
        a = verify.mixed_allowance(1.3e-3, self.MAXREF, self.SIGNAL)
        self.assertEqual(self.claim(self.LIM - a - 1e-4), "PASS")            # inside the noise-adjusted threshold
        self.assertEqual(self.claim(self.LIM - a + 1e-4), "CONDITIONAL")     # inside the deterministic limit, inside the allowance
        self.assertEqual(self.claim(self.LIM - 1e-4), "CONDITIONAL")         # marginal results stay CONDITIONAL
        self.assertEqual(self.claim(self.LIM + 1e-4), "CONDITIONAL")         # just beyond the limit (margin)
        self.assertEqual(self.claim(0.05), "FAIL")
        for v in self.GRID:
            self.assertEqual(self.claim(v), self.expected(v), v)

    def test_removing_the_allowance_fails_a_test(self):
        orig = verify.mixed_allowance
        verify.mixed_allowance = lambda noise, max_ref, signal: 0.0
        try:
            got = [self.claim(v) for v in self.GRID]
        finally:
            verify.mixed_allowance = orig
        self.assertNotEqual(got, [self.expected(v) for v in self.GRID])
        self.assertEqual(self.claim(0.025), "CONDITIONAL")                    # would be PASS without the allowance

    def test_adding_the_allowance_to_the_limit_fails_a_test(self):
        orig = verify.band_status
        verify.band_status = lambda v, lo, hi, margin, noise_lo=0.0, noise_hi=0.0: orig(v, lo, hi + noise_hi, margin, noise_lo, 0.0)
        try:
            got = [self.claim(v) for v in self.GRID]
        finally:
            verify.band_status = orig
        self.assertNotEqual(got, [self.expected(v) for v in self.GRID])

    def test_allowance_cannot_manufacture_pass(self):
        for a in (0.0, 1e-3, 1e-2, 0.1, 10.0):
            for v in self.GRID:
                with_a = verify.band_status(v, 0.0, self.LIM, 0.2, 0.0, a)
                without = verify.band_status(v, 0.0, self.LIM, 0.2, 0.0, 0.0)
                self.assertGreaterEqual(ST.severity(with_a), ST.severity(without), (a, v))
                if v > self.LIM:
                    self.assertNotEqual(with_a, "PASS")
        # an allowance larger than the limit leaves no PASS region at all
        self.assertEqual(verify.band_status(0.0, 0.0, self.LIM, 0.2, 0.0, 10.0), "CONDITIONAL")
        # the lower-edge allowance likewise only shrinks the band
        self.assertEqual(verify.band_status(1.0, 1.0, 2.0, 0.2, 0.05, 0.0), "CONDITIONAL")
        self.assertEqual(verify.band_status(1.06, 1.0, 2.0, 0.2, 0.05, 0.0), "PASS")

    def test_other_terms_use_their_own_allowance(self):
        a_s = verify.mixed_allowance(THR["numerical_noise"]["noise_E_by_qoi_combined"]["value"]["sig_rt"], 418.4, 90.3)
        lim = THR["layerB_deterministic"]["mixed_term_max_error_over_signal_M3_0_2"]["value"]["d_sig_rt"]
        mk = lambda v: verify.mixed_claim({"M3": {"0-2": {"d_sig_rt": dict(error_over_signal=v, max_abs_ref_combined=418.4, signal_max_abs=90.3)}}}, THR)[0]
        self.assertEqual(mk(lim - a_s - 1e-4), "PASS")
        self.assertEqual(mk(lim - a_s + 1e-4), "CONDITIONAL")

    def test_empty_mixed_input_is_not_pass(self):
        self.assertEqual(verify.mixed_claim({}, THR)[0], "NOT-VERIFIED")


class CommandLineRobustness(unittest.TestCase):
    """verify.py option errors give a one-line message and a nonzero status, never a traceback (L-gg-4)."""

    @classmethod
    def setUpClass(cls):
        cls.root = tempfile.mkdtemp(dir=os.environ.get("TMPDIR"))
        cls.work = os.path.join(cls.root, "w")
        os.makedirs(cls.work)

    def check(self, args, option):
        r = run_driver(args)
        self.assertEqual(r.returncode, 2, r.stdout + r.stderr)
        self.assertNotIn("Traceback", r.stdout + r.stderr)
        self.assertNotIn("IndexError", r.stdout + r.stderr)
        self.assertIn(option, r.stderr)
        self.assertIn("requires a value", r.stderr)

    def test_log_without_value(self):
        self.check([self.work, "--log"], "--log")

    def test_out_without_value(self):
        self.check([self.work, "--out"], "--out")

    def test_thresholds_without_value(self):
        self.check([self.work, "--thresholds"], "--thresholds")

    def test_option_followed_by_option(self):
        self.check([self.work, "--log", "--out", os.path.join(self.root, "o")], "--log")
        self.check([self.work, "--out", "--log", "x"], "--out")
        self.assertFalse(os.path.exists(os.path.join(self.root, "o")))

    def test_empty_out_and_bad_thresholds_file(self):
        for args in ([self.work, "--out", ""], [self.work, "--thresholds", os.path.join(self.root, "nope.json")]):
            r = run_driver(args)
            self.assertEqual(r.returncode, 2, r.stdout + r.stderr)
            self.assertNotIn("Traceback", r.stdout + r.stderr)
            self.assertIn("error:", r.stderr)
        bad = os.path.join(self.root, "bad.json")
        with open(bad, "w") as fh:
            fh.write("{ not json")
        r = run_driver([self.work, "--thresholds", bad])
        self.assertEqual(r.returncode, 2)
        self.assertNotIn("Traceback", r.stdout + r.stderr)

    def test_unknown_option_and_no_arguments(self):
        r = run_driver([self.work, "--bogus"])
        self.assertEqual(r.returncode, 2)
        self.assertNotIn("Traceback", r.stdout + r.stderr)
        r = run_driver([])
        self.assertEqual(r.returncode, 2)
        self.assertNotIn("Traceback", r.stdout + r.stderr)


class InvalidResidualEvidence(unittest.TestCase):
    """An unreadable FSI residual in a completed step is INCOMPLETE evidence, never PASS or CONDITIONAL (L-gg-5)."""

    @staticmethod
    def mutate(text, marker, occurrence, new):
        lines, seen = text.split("\n"), 0
        for i, ln in enumerate(lines):
            if ln.startswith(marker):
                seen += 1
                if seen == occurrence:
                    lines[i] = new
                    return "\n".join(lines)
        raise AssertionError("marker not found")

    R1 = "FSI relative residual1 norm for interface 0: "
    R2 = "FSI residual2 norm for interface 0: "
    AB = "FSI absolute residual (max point) for interface 0: "
    LAST = 40          # last iteration of the 10th step (4 residual prints per step)

    def base(self):
        return make_log(accept="both", n=60)

    def test_baseline_passes(self):
        r = run_a(self.base())
        self.assertEqual(r["status"], "PASS")
        self.assertNotIn("fsi_residual_evidence", r["checks"])

    def assert_incomplete(self, text):
        r = run_a(text)
        self.assertEqual(r["status"], "INCOMPLETE", r["status"])
        self.assertEqual(r["checks"]["fsi_residual_evidence"]["status"], "INCOMPLETE")
        self.assertNotIn(r["status"], ("PASS", "CONDITIONAL"))

    def test_infinity(self):
        self.assert_incomplete(self.mutate(self.base(), self.R1, self.LAST, self.R1 + "inf"))
        self.assert_incomplete(self.mutate(self.base(), self.R2, self.LAST, self.R2 + "-inf"))
        self.assert_incomplete(self.mutate(self.base(), self.AB, self.LAST, self.AB + "inf"))

    def test_nan_is_never_pass(self):
        r = run_a(self.mutate(self.base(), self.R2, self.LAST, self.R2 + "nan"))
        self.assertIn(r["status"], ("FAIL", "INCOMPLETE"))

    def test_missing_residual(self):
        for marker in (self.R2, self.AB):
            lines = self.base().split("\n")
            seen, out = 0, []
            for ln in lines:
                if ln.startswith(marker):
                    seen += 1
                    if seen == self.LAST:
                        continue
                out.append(ln)
            self.assert_incomplete("\n".join(out))

    def test_missing_residual1_line_in_the_deciding_iteration(self):
        lines, seen, out = self.base().split("\n"), 0, []
        for ln in lines:
            if ln.startswith(self.R1):
                seen += 1
                if seen == self.LAST:
                    continue
            out.append(ln)
        self.assert_incomplete("\n".join(out))
        # a previous iteration at 4e-7 must not be inherited by the step
        t = self.mutate("\n".join(out), self.R1, self.LAST - 1, self.R1 + "4e-7")
        self.assertNotIn(run_a(t)["status"], ("PASS", "CONDITIONAL"))

    def test_malformed_residual(self):
        for bad in ("1e-", "-", ".", "abc", "1.2.3e", "--5"):
            self.assert_incomplete(self.mutate(self.base(), self.R1, self.LAST, self.R1 + bad))

    def test_truncated_numerical_line(self):
        self.assert_incomplete(self.mutate(self.base(), self.R1, self.LAST, "FSI relative residual1 norm for interface 0:"))
        self.assert_incomplete(self.mutate(self.base(), self.R2, self.LAST, "FSI residual2 norm for interf"))
        self.assert_incomplete(self.mutate(self.base(), self.AB, self.LAST, "FSI absolute residual (max point) for interface 0: "))

    def test_unreadable_in_a_non_deciding_iteration_is_still_flagged_if_malformed(self):
        self.assert_incomplete(self.mutate(self.base(), self.R1, self.LAST - 1, self.R1 + "1e-"))

    def test_final_incomplete_step_of_a_truncated_log_is_not_judged_twice(self):
        r = run_a(make_log(accept="both", truncate=True, end=False))
        self.assertEqual(r["status"], "INCOMPLETE")
        self.assertEqual(r["checks"]["termination"]["status"], "INCOMPLETE")

    def test_never_pass_with_bad_evidence_in_any_acceptance_mode(self):
        for acc in ("both", "r1", "r2"):
            text = self.mutate(make_log(accept=acc, n=60), self.R2, self.LAST, self.R2 + "1e-")
            self.assertNotIn(run_a(text)["status"], ("PASS", "CONDITIONAL"))


class FamilyDuplicateAndDeadCode(unittest.TestCase):
    def test_family_normalise_removed_and_unreferenced(self):
        self.assertFalse(hasattr(FAM, "normalise"))
        for name in sorted(os.listdir(os.path.join(HERE, "scripts"))):
            if name.endswith(".py"):
                txt = open(os.path.join(HERE, "scripts", name)).read()
                self.assertNotIn("normalise(", txt.replace("_normalise(", ""), name)


class ThresholdsProvenanceText(unittest.TestCase):
    def test_json_valid_and_no_garbled_provenance(self):
        txt = open(os.path.join(HERE, "thresholds.json")).read()
        json.loads(txt)
        for bad in ("independent. ", "of the independent.", "1e-16.1e-13", "cross_check_the absolute", "  ."):
            self.assertNotIn(bad, txt)
        for k_, v_ in THR["layerB_deterministic"].items():
            if isinstance(v_, dict) and "provenance" in v_:
                self.assertNotRegex(v_["provenance"], r"\b(\w+) \1\b", k_)

    def test_self_convergence_provenance_is_qualified_and_values_unchanged(self):
        B = THR["layerB_deterministic"]
        self.assertEqual(B["combined_self_order_band"]["value"], [1.8, 2.15])
        prov = B["combined_self_order_band"]["provenance"]
        for need in ("sharp for alpha_B", "H_f", "noise-limited", "1.39-2.60", "1.39-2.54", "1.63-2.55", "not a proof"):
            self.assertIn(need, prov)
        self.assertNotIn("1.97-2.09", prov)

    def test_readme_wording_matches_provenance(self):
        rd = open(os.path.join(HERE, "..", "README.md")).read()
        self.assertNotIn("self-convergence about second order for five quantities", rd)
        self.assertNotIn("Near-second-order self-convergence belongs to the five quantities", rd)
        for need in ("1.39-2.60", "1.39-2.54", "1.63-2.55", "1.42-2.92"):
            self.assertIn(need, rd)
        src = open(os.path.join(HERE, "scripts", "verify.py")).read()
        self.assertIn("noise-limited, interval 1.42-2.92", src)
        self.assertNotIn("KE_s self-order 2.17 but", src)


if __name__ == "__main__":
    unittest.main(verbosity=2)
