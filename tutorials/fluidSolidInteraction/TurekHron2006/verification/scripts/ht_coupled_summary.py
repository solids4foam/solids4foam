#!/usr/bin/env python3
"""Spatial, temporal and coupling table of the production coupled FSI3 runs.

Reads the JSON of ht_coupled_analysis.py (all runs in the common window of
five u_y periods ending at the same time) and writes, per quantity:

- the spatial sequence 1x'/2x'/4x' at dt 5e-4 (solid 2x/4x/8x), the
  differences d21, d32 and the observed order p = log2(d21/d32) where both
  differences have one sign and shrink;
- the temporal sequence at 2x' (dt 1e-3/5e-4/2.5e-4), its observed order and
  the temporal-error estimate at the production dt: 4/3 of the 5e-4 -> 2.5e-4
  difference (Richardson for BDF2, the factor Astra asked for);
- the 1x' dt-halving difference and the coupling-tolerance difference
  (tol 1e-4 minus tol 1e-5 at 1x', an upper bound for the iterative error at
  tol 1e-5), and the 1x' solid-refinement difference (solid 4x minus 2x);
- the ratios to D = |q(4x') - q(2x')| (Astra criterion: time < 0.2 D,
  iteration < 0.05 D).

    python3 ht_coupled_summary.py coupled.json coupled_summary.csv
"""
import csv
import json
import math
import sys

QOIS = [("ux_timemean", 1e3, "mm"), ("ux_midrange", 1e3, "mm"), ("ux_amplitude", 1e3, "mm"),
        ("uy_amplitude", 1e3, "mm"), ("frequency", 1, "Hz"), ("drag_timemean", 1, "N/m"),
        ("drag_amplitude", 1, "N/m"), ("lift_timemean", 1, "N/m"), ("lift_amplitude", 1, "N/m"),
        ("drag_amplitude_smoothed", 1, "N/m"), ("lift_amplitude_smoothed", 1, "N/m")]


def main() -> int:
    runs = {r["name"]: r for r in json.load(open(sys.argv[1]))}

    def g(name, q, s):
        return runs[name][q] * s if name in runs and q in runs[name] else math.nan

    rows = []
    for q, s, u in QOIS:
        v1, v2, v4 = g("1x_dt5e-4", q, s), g("2x_dt5e-4", q, s), g("4x_dt5e-4", q, s)
        d21, d32 = v2 - v1, v4 - v2
        p = math.log2(d21 / d32) if d21 * d32 > 0 and abs(d21) > abs(d32) else math.nan
        t1, t3 = g("2x_dt1e-3", q, s), g("2x_dt2.5e-4", q, s)
        terr2 = (v2 - t3) * 4 / 3
        pt = math.log2((v2 - t1) / (t3 - v2)) if (v2 - t1) * (t3 - v2) > 0 else math.nan
        terr1 = (g("1x_dt1e-3", q, s) - v1) / 3
        tol = g("1x_dt5e-4_tol1e-4", q, s) - v1
        sol = g("1x_s4x_dt5e-4", q, s) - v1
        D = abs(d32)
        rows.append([q, u, v1, v2, v4, d21, d32, p, t1, t3, pt, terr2, terr1, tol, sol, D,
                     abs(terr2) / D if D else math.nan, abs(tol) / D if D else math.nan])
    with open(sys.argv[2], "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["qoi", "unit", "1x", "2x", "4x", "d21", "d32", "p_space", "2x_dt1e-3", "2x_dt2.5e-4",
                    "p_time_2x", "time_err_2x", "time_err_1x", "tol1e-4_minus_tol1e-5_1x",
                    "solid4x_minus_solid2x_1x", "D", "time_err_over_D", "tol_over_D"])
        for r in rows:
            w.writerow(r[:2] + [f"{x:.6g}" for x in r[2:]])
    for r in rows:
        print(f"{r[0]:15s} " + " ".join(f"{x:9.4f}" for x in r[2:8]) + " | t " + " ".join(f"{x:8.4f}" for x in r[10:13])
              + f" | tol {r[13]:8.4f} sol {r[14]:8.4f} | tE/D {r[16]:.2f} tolE/D {r[17]:.2f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
