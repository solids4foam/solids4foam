#!/usr/bin/env python3
"""Convergence analysis of the smooth-mesh replay family (ht_mesh_smooth.py).

Uses second_order_analysis.analyse (five periods of the motion ending at
--t-end, the protocol of the second-order study) on the cases
<work>/sm_1x, sm_2x, sm_3x, sm_4x and the check runs, and reports:

  * the values per level and the differences;
  * the 1x/2x/4x triplet order p = log2(d21/d32) where both differences have
    the same sign and shrink (no Richardson value otherwise);
  * the 2x/3x/4x three-point fit Q(h) = Q_inf + c h^p (h = 1/2, 1/3, 1/4),
    where monotone, and the ratio (Q2 - Q3)/(Q3 - Q4) against 20/7 = 2.857
    (the value for p = 2);
  * for the complex lift fundamental: the successive difference phasors,
    their magnitude ratio and the angle between them (a dominant
    single-order error gives aligned differences);
  * window-end sensitivity (window ending one period earlier), drift per
    second (second_order_analysis), the time-step check (2x at dt/2), the
    outlet-tail check and the converged-solve check, each as a fraction of the
    finest relevant difference |Q4 - Q3| and |Q4 - Q2|.

    python3 smooth_family_analysis.py --work <results> --json smooth_family.json
"""
from __future__ import annotations

import argparse
import cmath
import json
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import second_order_analysis as so  # noqa: E402

QOIS = ["Q_in", "Q_quad", "lift1_amp", "drag_mean", "lift1_re", "lift1_im"]
LEVELS = {"1x": 1, "2x": 2, "3x": 3, "4x": 4}


def fit3(q, h):
    """p and Q_inf of Q = Q_inf + c h^p through three points, if monotone."""
    d1, d2 = q[1] - q[0], q[2] - q[1]
    if d1 * d2 <= 0 or abs(d2) >= abs(d1) * (h[2] - h[1]) / (h[1] - h[0]) * 10:
        return None
    r = d1 / d2

    def g(p):
        return (h[0] ** p - h[1] ** p) / (h[1] ** p - h[2] ** p) - r
    lo, hi = 0.05, 12.0
    if g(lo) * g(hi) > 0:
        return None
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if g(lo) * g(mid) <= 0:
            hi = mid
        else:
            lo = mid
    p = 0.5 * (lo + hi)
    c = d2 / (h[2] ** p - h[1] ** p)
    return {"p": p, "q_inf": q[2] - c * h[2] ** p}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--work", type=Path, required=True)
    ap.add_argument("--t-end", type=float, default=30.2)
    ap.add_argument("--prefix", default="sm_")
    ap.add_argument("--json", type=Path)
    a = ap.parse_args()
    runs = {}
    for tag in ["1x", "2x", "3x", "4x", "2x_dt2.5e-4", "2x_tail90", "1x_conv"]:
        c = a.work / f"{a.prefix}{tag}"
        if (c / "energy.dat").is_file():
            r = so.analyse(c, a.t_end, 1.0)
            if "error" not in r:
                # window-end sensitivity: the same window one period earlier
                col, omega, a_tip, _ = so.load(c)
                e = so.window(col, omega, a_tip, a.t_end - 2 * math.pi / omega)
                r["window_shift"] = {q: r[q] - e[q] for q in QOIS} if e else None
            runs[tag] = r
    out = {"runs": runs}
    lv = [t for t in LEVELS if t in runs]
    print("level  " + "".join(f"{q:>12s}" for q in QOIS) + "      dt   dQin/dt")
    for t in lv:
        r = runs[t]
        print(f"{t:6s} " + "".join(f"{r[q]:12.3f}" for q in QOIS) + f"  {r['dt']:.2e}  {r.get('dQ_in_dt', float('nan')):7.3f}")
    res = {}
    for q in QOIS:
        v = {t: runs[t][q] for t in lv}
        e = {"values": v}
        if all(t in v for t in ("1x", "2x", "4x")):
            e["triplet_124"] = so.order(v["1x"], v["2x"], v["4x"])
        if all(t in v for t in ("2x", "3x", "4x")):
            d23, d34 = v["3x"] - v["2x"], v["4x"] - v["3x"]
            e["d23"], e["d34"] = d23, d34
            e["ratio_23_34"] = d23 / d34 if d34 else None
            e["fit_234"] = fit3([v["2x"], v["3x"], v["4x"]], [1 / 2, 1 / 3, 1 / 4])
        if all(t in v for t in ("1x", "2x", "3x")):
            e["fit_123"] = fit3([v["1x"], v["2x"], v["3x"]], [1.0, 1 / 2, 1 / 3])
        if all(t in v for t in ("1x", "3x", "4x")):
            e["fit_134"] = fit3([v["1x"], v["3x"], v["4x"]], [1.0, 1 / 3, 1 / 4])
        if all(t in v for t in ("1x", "2x", "3x", "4x")):
            seq = [v[t] for t in ("1x", "2x", "3x", "4x")]
            dif = [seq[i + 1] - seq[i] for i in range(3)]
            e["monotone_1234"] = all(d * dif[0] > 0 for d in dif) and abs(dif[0]) > abs(dif[1]) > abs(dif[2])
        res[q] = e
    out["qoi"] = res
    print("\nthree-point fits Q = Q_inf + c h^p (p, Q_inf) on 1x/2x/3x, 1x/3x/4x, 2x/3x/4x")
    for q in QOIS:
        e = res[q]
        print(f"  {q:10s} " + "   ".join(
            (f"{k[4:]}: p {e[k]['p']:5.2f} Q_inf {e[k]['q_inf']:9.3f}" if e.get(k) else f"{k[4:]}: -")
            for k in ("fit_123", "fit_134", "fit_234")))
    print("\nqoi        d21(1x-2x)  d42(2x-4x)  p(1,2,4)   d23      d34     ratio(20/7=2.857)  fit p(2,3,4)  Q_inf")
    for q in QOIS:
        e = res[q]
        tr = e.get("triplet_124", {})
        f = e.get("fit_234")
        print(f"{q:10s} {tr.get('d21', float('nan')):10.3f} {tr.get('d32', float('nan')):10.3f} "
              f"{tr.get('order', float('nan')):8.2f} {e.get('d23', float('nan')):8.3f} {e.get('d34', float('nan')):8.3f} "
              f"{(e.get('ratio_23_34') or float('nan')):10.3f}         "
              + (f"{f['p']:6.2f}   {f['q_inf']:9.3f}" if f else "   -") + f"   monotone {e.get('monotone_1234')}")
    # complex lift phasor differences
    ph = {t: complex(runs[t]["lift1_re"], runs[t]["lift1_im"]) for t in lv}
    seq = [t for t in ("1x", "2x", "3x", "4x") if t in ph]
    diffs = []
    print("\nlift fundamental phasor differences (N/m, deg)")
    for i in range(len(seq) - 1):
        d = ph[seq[i + 1]] - ph[seq[i]]
        diffs.append({"from": seq[i], "to": seq[i + 1], "re": d.real, "im": d.imag, "abs": abs(d),
                      "arg_deg": math.degrees(cmath.phase(d))})
        print(f"  {seq[i]}->{seq[i + 1]}: ({d.real:7.3f}, {d.imag:7.3f})  |d| {abs(d):6.3f}  arg {math.degrees(cmath.phase(d)):7.1f}")
    for t in seq:
        print(f"  {t}: phase of the lift fundamental {math.degrees(cmath.phase(ph[t])):7.3f} deg, amplitude {abs(ph[t]):8.3f}")
    out["lift_phasor_differences"] = diffs
    # checks relative to the finest differences
    checks = {}
    pairs = (("time step (2x dt/2 - 2x)", "2x_dt2.5e-4", "2x"), ("outlet tail (2x tail90 - 2x)", "2x_tail90", "2x"),
             ("converged solve (1x conv - 1x)", "1x_conv", "1x"))
    print("\nchecks: change, and as a fraction of |4x-3x| and |4x-2x|")
    for name, a1, b1 in pairs:
        if a1 in runs and b1 in runs:
            ch = {}
            for q in QOIS:
                dq = runs[a1][q] - runs[b1][q]
                f34 = abs(dq) / abs(res[q]["d34"]) if res[q].get("d34") else None
                f24 = abs(dq) / abs(runs["4x"][q] - runs["2x"][q]) if "4x" in runs and "2x" in runs else None
                ch[q] = {"change": dq, "frac_d34": f34, "frac_d42": f24}
            checks[name] = ch
            print(f"  {name}")
            for q, c in ch.items():
                print(f"     {q:10s} {c['change']:9.4f}   "
                      + (f"{c['frac_d34']:6.2f}" if c["frac_d34"] is not None else "   -") + "  "
                      + (f"{c['frac_d42']:6.2f}" if c["frac_d42"] is not None else "   -"))
    print("  window-end sensitivity (window ending one period earlier)")
    for t in lv:
        ws = runs[t].get("window_shift")
        if ws:
            print(f"     {t:4s} " + "  ".join(f"{q} {ws[q]:+.3f}" for q in QOIS))
    out["checks"] = checks
    # sliding-window spread (windows shifted by one period over the last second)
    import qin_slide
    spread = {}
    print("\nsliding five-period windows ending in the last second: mean, std, range (max-min)")
    for t in lv:
        w = [r for r in qin_slide.slide(a.work / f"{a.prefix}{t}") if r["t_end"] > a.t_end - 1.0 - 1e-6]
        spread[t] = {}
        for q, k in (("Q_in", "Q_in"), ("Q_quad", "Q_quad"), ("lift1_amp", "lift_amp1"), ("drag_mean", "drag_mean")):
            v = [r[k] for r in w]
            spread[t][q] = {"n": len(v), "mean": sum(v) / len(v), "std": float(__import__("numpy").std(v)),
                            "range": max(v) - min(v)}
        print(f"  {t:3s} n={len(w)} " + "  ".join(f"{q} {spread[t][q]['std']:.3f}/{spread[t][q]['range']:.3f}"
                                               for q in spread[t]))
    out["sliding_spread"] = spread
    # time-step share along the constant-Courant path (dt proportional to 1/L): with BDF2 and a
    # level-independent coefficient the temporal error at level L is e(L) = (4/3) D (2/L)^2, D the
    # measured 2x change Q(dt/2) - Q(dt), sign such that Q - e removes it. An assumption, not a
    # measurement at 3x/4x.
    if "2x_dt2.5e-4" in runs and "2x" in runs:
        tc = {}
        print("\ntime-corrected values Q - e_t(L), e_t(L) = -(4/3) [Q(2x,dt/2) - Q(2x,dt)] (2/L)^2")
        for q in QOIS:
            D = runs["2x_dt2.5e-4"][q] - runs["2x"][q]
            e = {t: -(4.0 / 3.0) * D * (2.0 / LEVELS[t]) ** 2 for t in lv}
            v = {t: runs[t][q] - e[t] for t in lv}
            ent = {"e_t": e, "values": v}
            if all(t in v for t in ("2x", "3x", "4x")):
                ent["d23"], ent["d34"] = v["3x"] - v["2x"], v["4x"] - v["3x"]
                ent["fit_234"] = fit3([v["2x"], v["3x"], v["4x"]], [1 / 2, 1 / 3, 1 / 4])
            if all(t in v for t in ("1x", "2x", "4x")):
                ent["triplet_124"] = so.order(v["1x"], v["2x"], v["4x"])
            tc[q] = ent
            f = ent.get("fit_234")
            tr = ent.get("triplet_124", {})
            print(f"  {q:10s} " + " ".join(f"{t} {v[t]:9.3f}" for t in lv)
                  + f"  d23 {ent.get('d23', float('nan')):7.3f} d34 {ent.get('d34', float('nan')):7.3f}"
                  + (f"  p(2,3,4) {f['p']:.2f}" if f else "  p(2,3,4) -")
                  + (f"  p(1,2,4) {tr['order']:.2f}" if "order" in tr else "  p(1,2,4) -"))
        out["time_corrected"] = tc
    if a.json:
        a.json.write_text(json.dumps(out, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
