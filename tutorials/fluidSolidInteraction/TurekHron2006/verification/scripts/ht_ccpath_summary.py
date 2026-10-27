#!/usr/bin/env python3
"""Constant-Courant path summary of the coupled FSI3 study (2x' dt 1e-3, 4x' dt 5e-4, 8x' dt 2.5e-4).

Reads the JSON written by ht_coupled_analysis.py (cases named as below) and prints, per
quantity, the three levels, the differences d1 = q(4x') - q(2x') and d2 = q(8x') - q(4x'), the
ratio d2/d1, and the observed order p = ln(d1/d2)/ln 2 when the sequence is monotone with
shrinking differences. Each difference is compared with the time sensitivity at 2x' (the dt
1e-3 -> 5e-4 change, the step the path takes between levels) and with the coupling sensitivity
(1x' tol 1e-4 against 1e-5). Optional: the fixed-dt sequence (1x', 2x', 4x' at dt 5e-4) for
reference.

    python3 ht_ccpath_summary.py ccpath.json --csv ccpath_summary.csv
"""
from __future__ import annotations

import argparse
import csv
import json
import math

PATH = ("2x_dt1e-3", "4x_dt5e-4", "8x_dt2.5e-4")
FIXED = ("1x_dt5e-4", "2x_dt5e-4", "4x_dt5e-4")
QUANTITIES = (  # key, label, scale
    ("uy_amplitude", "u_y amplitude (mm)", 1e3),
    ("ux_midrange", "u_x (max+min)/2 (mm)", 1e3),
    ("ux_amplitude", "u_x amplitude (mm)", 1e3),
    ("ux_timemean", "u_x time mean (mm)", 1e3),
    ("frequency", "frequency (Hz)", 1.0),
    ("drag_midrange", "drag (max+min)/2 (N/m)", 1.0),
    ("drag_amplitude_smoothed", "drag amplitude (N/m)", 1.0),
    ("lift_midrange", "lift (max+min)/2 (N/m)", 1.0),
    ("lift_amplitude_smoothed", "lift amplitude (N/m)", 1.0),
)


def load(path: str) -> dict:
    d = json.load(open(path))
    runs = d if isinstance(d, list) else d.get("runs", list(d.values()))
    return {r["name"]: r for r in runs if isinstance(r, dict) and "name" in r}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("json")
    ap.add_argument("--csv")
    a = ap.parse_args()
    r = load(a.json)
    rows = []
    for key, label, s in QUANTITIES:
        if not all(n in r for n in PATH):
            continue
        q = [r[n][key] * s for n in PATH]
        d1, d2 = q[1] - q[0], q[2] - q[1]
        mono = d1 * d2 > 0 and abs(d2) < abs(d1)
        p = math.log(abs(d1 / d2)) / math.log(2) if mono else float("nan")
        t = (r["2x_dt5e-4"][key] - r["2x_dt1e-3"][key]) * s if {"2x_dt5e-4", "2x_dt1e-3"} <= r.keys() else float("nan")
        c = (r["1x_dt5e-4_tol1e-4"][key] - r["1x_dt5e-4"][key]) * s \
            if {"1x_dt5e-4_tol1e-4", "1x_dt5e-4"} <= r.keys() else float("nan")
        fx = [r[n][key] * s for n in FIXED] if all(n in r for n in FIXED) else []
        row = {"quantity": label, "2x'_dt1e-3": q[0], "4x'_dt5e-4": q[1], "8x'_dt2.5e-4": q[2],
               "d1": d1, "d2": d2, "d2/d1": d2 / d1 if d1 else float("nan"), "p": p,
               "time_2x_1e-3_to_5e-4": t, "coupling_1x_tol": c,
               "|d2|/|time|": abs(d2 / t) if t else float("nan"),
               "|d2|/|coupling|": abs(d2 / c) if c else float("nan")}
        if fx:
            row.update({"fixed_1x": fx[0], "fixed_2x": fx[1], "fixed_4x": fx[2]})
        rows.append(row)
    w = max(len(x["quantity"]) for x in rows)
    print(f"{'quantity':{w}s} {'2x dt1e-3':>11s} {'4x dt5e-4':>11s} {'8x dt2.5e-4':>11s} {'d1':>9s} {'d2':>9s}"
          f" {'d2/d1':>7s} {'p':>5s} {'time':>9s} {'coupl':>9s}")
    k2, k4, k8 = "2x'_dt1e-3", "4x'_dt5e-4", "8x'_dt2.5e-4"
    for x in rows:
        print(f"{x['quantity']:{w}s} {x[k2]:11.4f} {x[k4]:11.4f} {x[k8]:11.4f} {x['d1']:+9.4f} {x['d2']:+9.4f}"
              f" {x['d2/d1']:7.2f} {x['p']:5.2f} {x['time_2x_1e-3_to_5e-4']:+9.4f} {x['coupling_1x_tol']:+9.4f}")
    if a.csv:
        with open(a.csv, "w", newline="") as fh:
            wr = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
            wr.writeheader()
            wr.writerows(rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
