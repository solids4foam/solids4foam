#!/usr/bin/env python3
"""Force harmonics of replay cases in fixed windows, drift, and observed orders.

For every case written by hron_turek_second_order.py / hron_turek_energy.py
(energy.dat), evaluates over five periods of the prescribed motion ending at
--t-end (the same window on every level, so that a slow drift cannot bias the
differences), the quantities of energy_analysis.py:

  Q_in, Q_quad      first-harmonic generalised force on the replayed y shape,
                    normalised by the tip amplitude (N/m); _p: pressure only
  lift1_re/im       complex first harmonic of the plate lift, phase referenced
                    to the motion (N/m), and its amplitude
  drag_mean         mean plate drag (N/m)

and the drift of Q_in and drag_mean, as the change per second between the
window ending at --t-end and one ending --drift-span s earlier.

Families: --family NAME=case1x,case2x,case4x prints the differences, the ratio
and the observed order p = log2(d21/d32) where both differences have the same
sign, and writes everything to --json/--csv.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import energy_analysis as ea  # noqa: E402

QOIS = ["Q_in", "Q_quad", "Q_in_p", "Q_quad_p", "lift1_re", "lift1_im", "lift1_amp", "drag_mean"]


def load(case: Path):
    s = ea.setup(case)
    k = float(s["amplitude_factor"])
    data = np.loadtxt(case / "energy.dat")
    _, idx = np.unique(data[:, 0], return_index=True)
    data = data[np.sort(idx)]
    col = {n: data[:, i] for i, n in enumerate(ea.COLS)}
    omega, d_tip, *_ = ea.plate_kinematics(case, k)
    return col, omega, abs(d_tip) * k, s


def window(col, omega, a_tip, t_end, periods=5):
    t = col["t"]
    period = 2 * math.pi / omega
    if t_end > t[-1] + 1e-9 or t_end - periods * period < t[0] + 0.2 - 1e-9:
        return None
    sel = (t > t_end - periods * period + 1e-9) & (t <= t_end + 1e-9)
    n = int(sel.sum())
    ph = np.exp(-1j * omega * col["T"][sel])
    r = {"t_end": float(t_end), "n": n}
    for tag, a, b in (("", "Q1ya", "Q1yb"), ("_p", "Qp1ya", "Qp1yb")):
        g = 2.0 / n * np.sum((col[a][sel] + 1j * col[b][sel]) * ph)
        r["Q_in" + tag] = float(g.real / a_tip)
        r["Q_quad" + tag] = float(g.imag / a_tip)
    fy = col["Fpy"][sel] + col["Fvy"][sel]
    fx = col["Fpx"][sel] + col["Fvx"][sel]
    lift = 2.0 / n * np.sum(fy * ph)
    r["lift1_re"], r["lift1_im"], r["lift1_amp"] = float(lift.real), float(lift.imag), float(abs(lift))
    r["drag_mean"] = float(fx.mean())
    return r


def analyse(case: Path, t_end: float | None, drift_span: float):
    col, omega, a_tip, s = load(case)
    te = t_end if t_end is not None else float(col["t"][-1])
    res = window(col, omega, a_tip, te)
    if res is None:
        return {"case": case.name, "error": f"window ending {te} not available (t {col['t'][0]}..{col['t'][-1]})"}
    res["case"] = case.name
    res["variants"] = s.get("variants", "base")
    res["dt"] = float(np.median(np.diff(col["t"])))
    early = window(col, omega, a_tip, te - drift_span)
    if early:
        res["dQ_in_dt"] = (res["Q_in"] - early["Q_in"]) / drift_span
        res["ddrag_dt"] = (res["drag_mean"] - early["drag_mean"]) / drift_span
    return res


def order(v1, v2, v3):
    d21, d32 = v2 - v1, v3 - v2
    out = {"v1": v1, "v2": v2, "v3": v3, "d21": d21, "d32": d32}
    if d21 != 0 and d32 != 0 and d21 * d32 > 0 and abs(d21) > abs(d32):
        out["ratio"] = d21 / d32
        out["order"] = math.log2(d21 / d32)
        out["richardson_limit"] = v3 + d32 / (d21 / d32 - 1)
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("cases", nargs="*", type=Path)
    ap.add_argument("--work", type=Path, default=Path.home() / "ht_2nd/work")
    ap.add_argument("--t-end", type=float)
    ap.add_argument("--drift-span", type=float, default=1.0)
    ap.add_argument("--family", action="append", default=[], help="NAME=c1,c2,c3")
    ap.add_argument("--json", type=Path)
    ap.add_argument("--csv", type=Path)
    a = ap.parse_args()
    names = [str(c) for c in a.cases]
    fams = []
    for f in a.family:
        name, cs = f.split("=", 1)
        fams.append((name, cs.split(",")))
        names += [c for c in cs.split(",") if c not in names]
    runs = {}
    for c in names:
        p = Path(c) if Path(c).is_absolute() else a.work / c
        runs[c] = analyse(p, a.t_end, a.drift_span)
        r = runs[c]
        if "error" in r:
            print(f"{c}: {r['error']}")
            continue
        print(f"{c:34s} t_end {r['t_end']:.3f}  Q_in {r['Q_in']:8.3f}  Q_quad {r['Q_quad']:8.3f}  "
              f"lift1 {r['lift1_amp']:8.2f}  drag {r['drag_mean']:8.3f}  "
              f"dQin/dt {r.get('dQ_in_dt', float('nan')):7.3f}  ddrag/dt {r.get('ddrag_dt', float('nan')):7.3f}")
    families = {}
    for name, cs in fams:
        if any("error" in runs[c] for c in cs):
            continue
        families[name] = {"cases": cs, "qoi": {q: order(*[runs[c][q] for c in cs]) for q in QOIS}}
        print(f"-- {name}: {cs}")
        for q, o in families[name]["qoi"].items():
            print(f"   {q:10s} {o['v1']:9.3f} {o['v2']:9.3f} {o['v3']:9.3f}  d21 {o['d21']:8.3f}  d32 {o['d32']:8.3f}"
                  + (f"  p {o['order']:.2f}  lim {o['richardson_limit']:.2f}" if "order" in o else "  p -"))
    if a.json:
        a.json.write_text(json.dumps({"runs": runs, "families": families}, indent=2))
    if a.csv:
        with a.csv.open("w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["family", "qoi", "v1", "v2", "v3", "d21", "d32", "order", "richardson_limit"])
            for name, f in families.items():
                for q, o in f["qoi"].items():
                    w.writerow([name, q] + [f"{o[k]:.6g}" if k in o else "" for k in
                                            ("v1", "v2", "v3", "d21", "d32", "order", "richardson_limit")])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
