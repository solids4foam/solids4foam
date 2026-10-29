#!/usr/bin/env python3
"""Convergence analysis and mixed-term analysis of radialTorsionalRing runs (portable).

Reads the case directories WORK/<config>_<mesh> (config = torsion|radial|combined,
mesh = M1|M2|M3) produced by verification/Allverify (or by ./Allrun, copied or linked
there), using analyse_case.load for the QoIs and the reference of each case.

Conventions (identical to the study):
  E_m    = max_{t in W, t > 2 dt} |q_m - q_ref| / max_{t in W, t > 2 dt} |q_ref|   (max-norm)
  order  p12 = log2(E1/E2), p23 = log2(E2/E3), reported only if E1 > E2 > E3 > 0
  self   d12 = max|q1 - q2|, d23 = max|q2 - q3| on the common grid, normalised with the
         M2 reference scale; p_self = log2(d12/d23)   (reference independent)
  mixed  d_alpha   = alpha_C - alpha_T            d_T       = T_B,C - T_B,T
         d_sig_rt  = sigma_rt,C - sigma_rt,T      d_u3      = u_C - u_R - u_T
         d_sig_rr3 = sigma_rr,C - sigma_rr,R - sigma_rr,T
         error/signal = max|mixed_num - mixed_ref| / max|mixed_ref|
Windows W: 0-2 s and 0-1.8 s.
Primary rotational QoIs (order claims): alpha_B, w_wall, T_B, sig_rt, T_A, H_s, H_f,
P_wall, KE_s. Secondary (reported, EXCLUDED from order claims because of the
normal-traction discrepancy): sig_rr, p_b, u_B.

Usage: convergence.py [WORK_DIR]   (default: ./work next to this tutorial's verification/)
Writes WORK_DIR/convergence.json and WORK_DIR/convergence.csv.
"""
import sys as _sys
_sys.dont_write_bytecode = True
import csv
import json
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import analyse_case as AC  # noqa: E402

PRIMARY = ["alpha_B", "w_wall", "T_B", "sig_rt", "T_A", "H_s", "H_f", "P_wall", "KE_s"]
SECONDARY = ["sig_rr", "p_b", "u_B"]
QOI = {"torsion": PRIMARY + SECONDARY,
       "combined": ["u_B", "alpha_B", "T_B", "sig_rt", "sig_rr", "T_A", "T_bal", "p_b", "H_s", "H_f", "w_wall", "P_wall", "KE_s"],
       "radial": ["u_B", "sig_rr", "p_b", "bdot", "P_wall", "KE_s", "W_s"]}
WINDOWS = {"0-2": (0.0, 2.0), "0-1.8": (0.0, 1.8)}
MESHES = ("M1", "M2", "M3")


OVERRIDE = {}   # (config, mesh) -> already loaded case (used by verify.py)


def load_case(work, config, mesh):
    if (config, mesh) in OVERRIDE:
        return OVERRIDE[(config, mesh)]
    c = os.path.join(work, f"{config}_{mesh}")
    if not os.path.isfile(os.path.join(c, "rtDiagFluid.dat")):
        return None
    t, num, ref, meta = AC.load(c)
    return dict(t=t, num=num, ref=ref, dt=float(t[1] - t[0]))


def select(S, a, b):
    t = S["t"]
    return (t > max(a, 2 * S["dt"]) + 1e-9) & (t <= b + 1e-9)


def common_grid(cases, a, b):
    """Boolean masks selecting the same time levels in every case (rounded to 1e-9 s)."""
    keys = [set(np.round(S["t"][select(S, a, b)], 9)) for S in cases]
    common = set.intersection(*keys)
    return [np.array([round(x, 9) in common for x in S["t"]]) & select(S, a, b) for S in cases]


NOISE = {"torsion": {"default": 2.4e-4}, "combined": {"default": 1.3e-3}}   # max-norm noise of E in units of max|ref| per QoI (thresholds.json, numerical_noise)


def delta_for(config, q):
    n = NOISE.get(config, {"default": 0.0})
    return n.get(q, n["default"])


def order_interval(E, delta):
    """First-order propagation of the noise delta (units of max|ref|) into the orders log2(E1/E2), log2(E2/E3).
    Returns {'p12': [p, sigma, lo, hi, worst_lo, worst_hi], 'p23': [...]} (None entries if delta is not below the error)."""
    out = {}
    for name, (a, b) in (("p12", (E[0], E[1])), ("p23", (E[1], E[2]))):
        if not (a > b > 0) or delta >= b:
            out[name] = None
            continue
        p = math.log2(a / b)
        sg = math.sqrt((delta / a) ** 2 + (delta / b) ** 2) / math.log(2)
        wlo = math.log2((a - delta) / (b + delta))
        whi = math.log2((a + delta) / (b - delta))
        out[name] = [p, sg, p - sg, p + sg, wlo, whi]
    return out


def self_order_interval(d12, d23, delta):
    if not (d12 > d23 > 0) or math.sqrt(2) * delta >= d23:
        return None
    dd = math.sqrt(2) * delta
    p = math.log2(d12 / d23)
    sg = math.sqrt((dd / d12) ** 2 + (dd / d23) ** 2) / math.log(2)
    return [p, sg, p - sg, p + sg]


def porder(e):
    if e[0] > e[1] > e[2] > 0:
        return [math.log2(e[0] / e[1]), math.log2(e[1] / e[2])]
    return "non-monotone"


def convergence(work, config, noise=None):
    S = [load_case(work, config, m) for m in MESHES]
    if any(s is None for s in S):
        return {"status": "incomplete: missing " + ",".join(m for m, s in zip(MESHES, S) if s is None)}
    out = {}
    for wn, (a, b) in WINDOWS.items():
        masks = common_grid(S, a, b)
        res = {}
        for q in QOI[config]:
            if q not in S[0]["num"] or q not in S[0]["ref"]:
                continue
            sc = [np.max(np.abs(s["ref"][q][m])) for s, m in zip(S, masks)]
            if sc[1] == 0:
                continue
            E = [float(np.max(np.abs(s["num"][q][m] - s["ref"][q][m])) / x) for s, m, x in zip(S, masks, sc)]
            n = [s["num"][q][m] for s, m in zip(S, masks)]
            L2 = [float(np.sqrt(np.mean((s["num"][q][m] - s["ref"][q][m]) ** 2)) / x) for s, m, x in zip(S, masks, sc)]
            d12 = float(np.max(np.abs(n[0] - n[1])) / sc[1])
            d23 = float(np.max(np.abs(n[1] - n[2])) / sc[1])
            delta = (noise if noise is not None else delta_for(config, q))
            res[q] = dict(E=E, E_L2=L2, observed_order=porder(E), observed_order_L2=porder(L2),
                          order_interval_max_norm=order_interval(E, delta), noise_delta=delta,
                          self_order_interval=self_order_interval(d12, d23, delta), d12=d12, d23=d23,
                          self_order=(math.log2(d12 / d23) if d12 > 0 and d23 > 0 else None),
                          order_claim_allowed=q in PRIMARY and config == "torsion")
        out[wn] = res
    return out


def mixed_terms(work, mesh):
    C, T, R = (load_case(work, c, mesh) for c in ("combined", "torsion", "radial"))
    if C is None or T is None:
        return {"status": "needs combined and torsion"}
    out = {}
    for wn, (a, b) in WINDOWS.items():
        cases = [x for x in (C, T, R) if x is not None]
        masks = common_grid(cases, a, b)
        mc, mt = masks[0], masks[1]
        g = lambda S, m, q, w: S[w][q][m]  # noqa: E731
        defs = {"d_alpha": ("alpha_B", [(C, mc, 1), (T, mt, -1)]),
                "d_Hf": ("H_f", [(C, mc, 1), (T, mt, -1)]),
                "d_T": ("T_B", [(C, mc, 1), (T, mt, -1)]),
                "d_sig_rt": ("sig_rt", [(C, mc, 1), (T, mt, -1)])}
        if R is not None:
            mr = masks[2]
            defs["d_u3"] = ("u_B", [(C, mc, 1), (R, mr, -1), (T, mt, -1)])
            defs["d_sig_rr3"] = ("sig_rr", [(C, mc, 1), (R, mr, -1), (T, mt, -1)])
        res = {}
        for k, (q, terms) in defs.items():
            ref = sum(sg * g(S, m, q, "ref") for S, m, sg in terms)
            num = sum(sg * g(S, m, q, "num") for S, m, sg in terms)
            sig = float(np.max(np.abs(ref)))
            err = float(np.max(np.abs(num - ref)))
            refmax = float(np.max(np.abs(g(C, mc, q, "ref"))))
            res[k] = dict(signal_max_abs=sig, error_max_abs=err, error_over_signal=(err / sig if sig > 0 else None),
                          qoi=q, max_abs_ref_combined=refmax)
        out[wn] = res
    return out


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    work = sys.argv[1] if len(sys.argv) > 1 else os.path.join(here, "..", "work")
    out = {"convergence": {c: convergence(work, c) for c in ("torsion", "combined", "radial")},
           "mixed_terms": {m: mixed_terms(work, m) for m in MESHES}}
    with open(os.path.join(work, "convergence.json"), "w") as fh:
        json.dump(out, fh, indent=1, default=lambda x: x.item() if hasattr(x, "item") else str(x))
    with open(os.path.join(work, "convergence.csv"), "w", newline="") as fh:
        wr = csv.writer(fh)
        wr.writerow(["config", "window", "qoi", "E_M1", "E_M2", "E_M3", "p12", "p23", "d12", "d23", "p_self", "order_claim_allowed"])
        for c, byw in out["convergence"].items():
            if "status" in byw:
                continue
            for w, res in byw.items():
                for q, r in res.items():
                    o = r["observed_order"]
                    wr.writerow([c, w, q] + [f"{x:.4e}" for x in r["E"]] + ([f"{x:.3f}" for x in o] if isinstance(o, list) else ["n/m", "n/m"])
                                + [f"{r['d12']:.4e}", f"{r['d23']:.4e}", "" if r["self_order"] is None else f"{r['self_order']:.3f}", r["order_claim_allowed"]])
    for c, byw in out["convergence"].items():
        print(f"== {c}: {byw.get('status', 'complete')}")
        if "status" in byw:
            continue
        for w, res in byw.items():
            for q, r in res.items():
                o = r["observed_order"]
                ps = "n/a" if r["self_order"] is None else f"{r['self_order']:.2f}"
                print(f"  {w:6s} {q:8s} E = " + " ".join(f"{x:.2e}" for x in r["E"])
                      + "  p = " + (" ".join(f"{x:.2f}" for x in o) if isinstance(o, list) else o)
                      + f"  p_self = {ps}")
    for m, byw in out["mixed_terms"].items():
        if "status" in byw:
            continue
        for w, res in byw.items():
            print(f"== mixed terms {m} {w}: " + " ".join(f"{k}: {v['error_over_signal']:.3f}" for k, v in res.items()))


if __name__ == "__main__":
    main()
