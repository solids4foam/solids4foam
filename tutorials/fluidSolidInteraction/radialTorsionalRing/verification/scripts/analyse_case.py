#!/usr/bin/env python3
"""QoI extraction, reference errors, conservation diagnostics and coupling statistics
for one radialTorsionalRing case (portable: uses only the case directory).

Usage: analyse_case.py CASE_DIR        (writes CASE_DIR/analysis.json and CASE_DIR/series.csv)

Inputs (case dir): rtDiagFluid.dat, rtDiagSolid.dat, reference.json, log.solids4Foam.
All extensive quantities are converted to per radian, per unit depth (divide by 2 pi dz).
Errors: E = max_t |num - ref| / max_t |ref| over t > 2 dt (max-norm, same convention as the
study). Conservation-balance residuals are also reported per time window.

Derived balances (all from the CFD's own data):
  fluid AM:     T_bal = -dH_f/dt + T_c - F_c            (torque of the fluid on the wall)
  solid AM:     T_A,num = dH_s/dt - T_wall              (drive torque applied at R = A)
  fluid energy: res_Ef = dKE_f/dt + P_wall - P_c + D_f  (P_wall = power fluid -> solid)
  solid energy: P_A,num = d(KE_s + W_s)/dt - P_wall      (drive power at R = A)
"""
import sys as _sys
_sys.dont_write_bytecode = True
import json
import math
import os
import re
import sys

import numpy as np


WINDOWS = [(0, 0.5), (0.5, 1.0), (1.0, 1.5), (1.5, 1.8), (1.8, 2.0), (0, 2.0)]


def load(C, dz=0.1):
    F = np.loadtxt(os.path.join(C, "rtDiagFluid.dat"), ndmin=2)
    S = np.loadtxt(os.path.join(C, "rtDiagSolid.dat"), ndmin=2)
    n = min(len(F), len(S))
    F, S = F[:n], S[:n]
    ref = json.load(open(os.path.join(C, "reference.json")))
    tr = np.array(ref["t"])
    t = F[:, 0]
    assert np.allclose(t, S[:, 0])
    k2 = 2 * math.pi * dz
    R = {k: np.interp(t, tr, np.array(ref[k])) for k in ref if isinstance(ref[k], list) and len(ref[k]) == len(tr)}
    N = ref["args"]["ntheta"]
    cpol = math.cos(math.pi / N)
    num = dict(
        u_B=S[:, 2], alpha_B=S[:, 3], alpha_A=S[:, 9], u_A=S[:, 8], n2_amp=np.hypot(S[:, 4], S[:, 5]),
        bdot=F[:, 17] / cpol, w_wall=F[:, 16] / cpol,
        p_b=F[:, 3], sig_nv=F[:, 5], sig_rr=F[:, 4], sig_rt=F[:, 6], T_B=F[:, 7] / k2, Q=F[:, 14] / k2,
        H_f=F[:, 9] / k2, KE_f=F[:, 10] / k2, D_f=F[:, 11] / k2, P_wall=F[:, 8] / k2, P_c=F[:, 15] / k2,
        F_c=F[:, 12] / k2, T_c=F[:, 13] / k2,
        H_s=S[:, 13] / k2, KE_s=S[:, 14] / k2, W_s=S[:, 15] / k2, lmin=S[:, 16], Jmin=S[:, 17], Jmax=S[:, 18])
    d = lambda x: np.gradient(x, t, edge_order=2)  # noqa: E731
    num["T_bal"] = -d(num["H_f"]) + num["T_c"] - num["F_c"]
    num["T_A"] = d(num["H_s"]) - num["T_B"]
    num["P_A"] = d(num["KE_s"] + num["W_s"]) - num["P_wall"]
    num["res_Ef"] = d(num["KE_f"]) + num["P_wall"] - num["P_c"] + num["D_f"]
    refq = dict(
        u_B=R["u_Rf"], alpha_B=R["alpha_Rf"], alpha_A=R["alpha_A"], u_A=0 * t, n2_amp=0 * t,
        bdot=R["bdot"], w_wall=R["b"] * R["omega_B"],
        p_b=R["p_b"], sig_nv=R["sig_nv_b"], sig_rr=R["sig_rr_b"], sig_rt=R["sig_rt_b"], T_B=R["T_B"], Q=R["Q"],
        H_f=R["H_f"], KE_f=R["KE_f"], D_f=R["D_f"], P_wall=R["sig_rr_b"] * R["b"] * R["bdot"] + R["T_B"] * R["omega_B"],
        P_c=R["P_c_in"], F_c=R["H_flux_out"], T_c=0 * t, H_s=R["H_s"], KE_s=R["KE_s"], W_s=R["V_s"],
        T_bal=R["T_B"], T_A=R["T_A"], P_A=R["P_A"], res_Ef=0 * t)
    return t, num, refq, ref


def coupling_stats(C):
    """Outer coupling iterations and the criterion that accepted each step (raw log, same definition as
    acceptance.py / verify.py: the criterion of the LAST iteration of the step).

    The solver accepts a step when min(residual1, residual2) <= outerCorrTolerance (legacy minimum test,
    requireAllResidualMeasures unset) or when the absolute point mismatch is below outerCorrAbsoluteTolerance. With
    allowUnconvergedCoupling unset an unconverged step aborts the run (FOAM FATAL ERROR). The audits' printed-line
    counts of absolute acceptances (76 and 83 for M1 and M2 torsion) include the solid-predictor evaluation before
    iteration 1 and therefore differ from the counts here (68 and 70)."""
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import acceptance as ACC
    thr = dict(layerA={"p_final_residual_max": {"value": 1e-12}, "U_final_residual_max": {"value": 1e-13},
                       "final_residual_slack": {"value": 1e-6}, "outerCorrTolerance_default": {"value": 1e-6},
                       "nOuterCorr_default": {"value": 100}})
    cfg = ACC.case_settings(C, thr)
    rec = ACC.parse_log(os.path.join(C, "log.solids4Foam"), thr, cfg)
    cnt = {"both": 0, "r1-only": 0, "r2-only": 0, "absolute": 0, "none": 0}
    its = []
    for st in rec["steps"]:
        c, _ = ACC.classify_step(st, cfg["tol"], cfg["absTol"])
        cnt[c] += 1
        its.append(len(st["trip"]) - 1 if st["trip"] else 0)
    its = np.array(its)
    return dict(steps=int(its.size), outer_iters_mean=float(its.mean()), outer_iters_max=int(its.max()),
                criterion_counts=cnt, steps_accepted_by_absolute_tol=cnt["absolute"],
                steps_passed_by_single_measure_only=cnt["r1-only"] + cnt["r2-only"])


def main():
    C = sys.argv[1]
    t, num, refq, ref = load(C)
    out = {"case": os.path.basename(os.path.abspath(C)), "args": ref["args"], "t_end_reached": float(t[-1])}
    sel = t > 2 * (t[1] - t[0])
    for k in num:
        if k not in refq:
            continue
        sc = float(np.max(np.abs(refq[k][sel])))
        e = float(np.max(np.abs(num[k][sel] - refq[k][sel])))
        out[k] = dict(max_abs_err=e, max_ref=sc, rel=e / sc if sc > 0 else None, max_num=float(np.max(np.abs(num[k][sel]))))
    # normal traction decomposition relative to max |sigma_rr|
    srr = out["sig_rr"]["max_ref"]
    out["normal_decomposition_rel_srr"] = dict(
        pressure=out["p_b"]["max_abs_err"] / srr, viscous=out["sig_nv"]["max_abs_err"] / srr,
        total=out["sig_rr"]["max_abs_err"] / srr)
    nvr = np.abs(refq["sig_nv"]) > 0.2 * np.max(np.abs(refq["sig_nv"])) if np.max(np.abs(refq["sig_nv"])) > 0 else None
    if nvr is not None and nvr.any():
        out["nv_fraction_median"] = float(np.median(num["sig_nv"][nvr] / refq["sig_nv"][nvr]))
    Tsc = out["T_B"]["max_ref"]
    if Tsc > 0:
        out["wall_torque_vs_AM_balance_rel"] = float(np.max(np.abs(num["T_B"][sel] - num["T_bal"][sel])) / Tsc)
    E = np.max(np.abs(refq["KE_f"])) + np.max(np.abs(refq["KE_s"])) + np.max(np.abs(refq["W_s"]))
    out["fluid_energy_residual_max_rel_to_power"] = float(np.max(np.abs(num["res_Ef"][sel])) / max(np.max(np.abs(refq["P_wall"])), 1e-300))
    out["notes"] = {"sig_nv": "sig_nv (viscous normal wall stress) differs from the reference by O(1) (about 0.7 of the reference): the CFD transmits about 0.29-0.32 of the reference value "
                              "(ratio 0.286, 0.307, 0.317 on M1, M2, M3 combined), without mesh convergence. It is a documented discrepancy, not an error of this script; "
                              "it is excluded from every order claim (see README).",
                    "torque_units": "T_B, T_A, T_c, H_f, H_s are per radian and per metre of depth (rtDiag value / (2 pi dz)); multiply by 2 pi for the whole annulus per metre of depth"}
    out["coupling"] = coupling_stats(C)
    # momentum / energy balance windows (relative to the reference torque / interface power in the window)
    dt0 = float(t[1] - t[0])
    win = {}
    for a, b in WINDOWS:
        s_ = (t > max(a, 2 * dt0) + 1e-9) & (t <= b + 1e-9)
        if not s_.any():
            continue
        Ts = float(np.max(np.abs(refq["T_B"][s_])))
        Ps = float(np.max(np.abs(refq["P_wall"][s_])))
        win[f"{a:g}-{b:g}"] = dict(
            AM_balance_rel_T=float(np.max(np.abs(num["T_B"][s_] - num["T_bal"][s_])) / Ts) if Ts > 0 else None,
            fluid_energy_residual_rel_P=float(np.max(np.abs(num["res_Ef"][s_])) / Ps) if Ps > 0 else None,
            solid_energy_residual_rel_P=float(np.max(np.abs(num["P_A"][s_] - refq["P_A"][s_])) / Ps) if Ps > 0 else None)
    out["balance_windows"] = win
    Psc = float(np.max(np.abs(refq["P_wall"][sel])))
    out["solid_energy_residual_rel_to_max_interface_power"] = float(np.max(np.abs(num["P_A"][sel] - refq["P_A"][sel])) / Psc)
    out["fluid_energy_residual_rel_to_max_interface_power"] = float(np.max(np.abs(num["res_Ef"][sel])) / Psc)
    out["min_principal_stretch"] = float(np.min(num["lmin"]))
    out["J_range"] = [float(np.min(num["Jmin"])), float(np.max(num["Jmax"]))]
    out["n2_amp_max"] = float(np.max(num["n2_amp"]))
    out["n2_amp_series_0p1s"] = [[float(x), float(y)] for x, y in zip(t[::40], num["n2_amp"][::40])]
    with open(os.path.join(C, "analysis.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    keys = list(num)
    with open(os.path.join(C, "series.csv"), "w") as fh:
        fh.write("t," + ",".join(f"num_{k},ref_{k}" for k in keys) + "\n")
        for i in range(t.size):
            fh.write(f"{t[i]:.6g}," + ",".join(f"{num[k][i]:.9e},{(refq[k][i] if k in refq else float('nan')):.9e}" for k in keys) + "\n")
    short = {k: (round(v["rel"], 6) if isinstance(v, dict) and v.get("rel") is not None else v)
             for k, v in out.items() if k not in ("args", "n2_amp_series_0p1s")}
    print(json.dumps(short, indent=1, default=str))


if __name__ == "__main__":
    main()
