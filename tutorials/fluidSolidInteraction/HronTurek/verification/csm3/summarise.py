#!/usr/bin/env python3
"""Collect the CSM3 / CSM2 structural convergence results into CSV + JSON.

Reads work/<case>/ (written by run_case.sh) and writes results/*.csv and
results/csm3_structural_convergence.json. Observed order of a monotone
three-level sequence with ratio r = 2: p = log2((Q2-Q1)/(Q3-Q2)).
"""
import csv, glob, json, math, os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from analyse import analyse

here = os.path.dirname(os.path.abspath(__file__))
work = os.path.join(here, "work"); out = os.path.join(here, "results")
os.makedirs(out, exist_ok=True)
L, NX1, NY1, TH = 0.35101, 105, 6, 0.02          # flag length, 1x cells, thickness

# Featflow (Turek & Hron) references, mm. CSM3: dt = 0.005 s, level 4+0
REF3 = {"ux_mean": -14.305, "ux_amp": 14.305, "uy_mean": -63.607,
        "uy_amp": 65.160, "freq": 1.0995}
REF2 = {"ux": -0.469000, "uy": -16.9739}          # CSM2 static, level 5+1

def order(q):
    """q = [Q1, Q2, Q3] on h, h/2, h/4: (ratio, p) or (nan, 'non-monotone')."""
    d1, d2 = q[1] - q[0], q[2] - q[1]
    if d1 * d2 <= 0:
        return float("nan"), float("nan")
    r = d1 / d2
    return r, math.log2(r)

def level_rows(prefix, qois):
    rows = []
    for f in (1, 2, 4, 8):
        c = os.path.join(work, prefix.format(f=f))
        if os.path.exists(os.path.join(c, "runtime.s")):
            rows.append((f, c))
    return rows

# ---------------------------------------------------------------- CSM3
csm3 = []
for f, c in level_rows("L{f}x_dt1e-3", None):
    a = analyse(c)
    csm3.append({"level": f, "nx": NX1*f, "ny": NY1*f, "cells": NX1*NY1*f*f,
                 "h_mm": 1e3*L/(NX1*f), "dt": 1e-3, "snes_rtol": 1e-6,
                 "runtime_s": a.get("runtime_s"),
                 "ux_mean": 1e3*a["ux"]["mean"], "ux_amp": 1e3*a["ux"]["amp"],
                 "uy_mean": 1e3*a["uy"]["mean"], "uy_amp": 1e3*a["uy"]["amp"],
                 "freq": a["uy"]["freq"]})
keys = list(REF3)
with open(os.path.join(out, "csm3_levels.csv"), "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["level", "nx", "ny", "cells", "h_mm", "dt_s", "snes_rtol",
                "runtime_s"] + keys + [k + "_err_vs_featflow_pct" for k in keys])
    for r in csm3:
        w.writerow([r["level"], r["nx"], r["ny"], r["cells"], f"{r['h_mm']:.4f}",
                    r["dt"], r["snes_rtol"], r["runtime_s"]]
                   + [f"{r[k]:.5f}" for k in keys]
                   + [f"{100*(r[k]-REF3[k])/abs(REF3[k]):+.2f}" for k in keys])
orders3 = []
for k in keys:
    for i in range(len(csm3) - 2):
        q = [csm3[i+j][k] for j in range(3)]
        r, p = order(q)
        orders3.append({"qoi": k, "levels": f"{csm3[i]['level']}x-{csm3[i+2]['level']}x",
                        "d12": q[1]-q[0], "d23": q[2]-q[1], "ratio": r, "p": p,
                        "rel_change_last_pct": 100*(q[2]-q[1])/abs(q[2])})
with open(os.path.join(out, "csm3_observed_order.csv"), "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["qoi", "levels", "diff_12", "diff_23", "diff_ratio", "observed_p",
                "rel_change_last_pct"])
    for o in orders3:
        w.writerow([o["qoi"], o["levels"], f"{o['d12']:.5f}", f"{o['d23']:.5f}",
                    f"{o['ratio']:.3f}", f"{o['p']:.2f}", f"{o['rel_change_last_pct']:+.2f}"])

# time-step check
dtrows = []
for name, nx, dt in (("L2x_dt1e-3", 210, 1e-3), ("L2x_dt5e-4", 210, 5e-4),
                     ("L2x_dt2.5e-4", 210, 2.5e-4), ("L4x_dt1e-3", 420, 1e-3),
                     ("L4x_dt5e-4", 420, 5e-4), ("I_L2x_tight", 210, 1e-3)):
    c = os.path.join(work, name)
    if os.path.exists(os.path.join(c, "runtime.s")):
        a = analyse(c)
        dtrows.append({"case": name, "dt": dt,
                       **{k: 1e3*a[k.split("_")[0]][k.split("_")[1]] for k in keys[:4]},
                       "freq": a["uy"]["freq"]})
with open(os.path.join(out, "csm3_dt_check.csv"), "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["case", "dt_s"] + keys)
    for r in dtrows:
        w.writerow([r["case"], r["dt"]] + [f"{r[k]:.5f}" for k in keys])

# ---------------------------------------------------------------- CSM2 static
def static(name):
    c = os.path.join(work, name)
    if not os.path.exists(os.path.join(c, "runtime.s")):
        return None
    d = np.loadtxt(glob.glob(c + "/postProcessing/*/*.dat")[0], comments="#")
    d = np.atleast_2d(d)[-1]
    if abs(d[2]) == 0:
        return None
    return {"ux": 1e3*d[1], "uy": 1e3*d[2], "runtime_s": int(open(c + "/runtime.s").read())}

variants = {"baseline": "S2_csm2_{f}x", "stab0.25": "D_stab025_{f}x",
            "stab1.0": "D_stab1_{f}x", "highOrder_p2": "D_ho2_{f}x",
            "highOrder_p3": "D_ho3_{f}x"}
static_rows, static_orders = [], []
for vname, pat in variants.items():
    series = []
    for f in (1, 2, 4, 8):
        s = static(pat.format(f=f))
        if s:
            series.append((f, s))
            static_rows.append({"variant": vname, "level": f, "nx": NX1*f, "ny": NY1*f,
                                "cells": NX1*NY1*f*f, **s,
                                "uy_err_vs_featflow_pct": 100*(s["uy"]-REF2["uy"])/abs(REF2["uy"]),
                                "ux_err_vs_featflow_pct": 100*(s["ux"]-REF2["ux"])/abs(REF2["ux"])})
    for k in ("ux", "uy"):
        for i in range(len(series) - 2):
            q = [series[i+j][1][k] for j in range(3)]
            r, p = order(q)
            static_orders.append({"variant": vname, "qoi": k,
                                  "levels": f"{series[i][0]}x-{series[i+2][0]}x",
                                  "ratio": r, "p": p})
with open(os.path.join(out, "csm2_static_levels.csv"), "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["variant", "level", "nx", "ny", "cells", "ux_mm", "uy_mm",
                "ux_err_vs_featflow_pct", "uy_err_vs_featflow_pct", "runtime_s"])
    for r in static_rows:
        w.writerow([r["variant"], r["level"], r["nx"], r["ny"], r["cells"],
                    f"{r['ux']:.5f}", f"{r['uy']:.5f}", f"{r['ux_err_vs_featflow_pct']:+.2f}",
                    f"{r['uy_err_vs_featflow_pct']:+.2f}", r["runtime_s"]])
with open(os.path.join(out, "csm2_static_observed_order.csv"), "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["variant", "qoi", "levels", "diff_ratio", "observed_p"])
    for o in static_orders:
        w.writerow([o["variant"], o["qoi"], o["levels"], f"{o['ratio']:.3f}", f"{o['p']:.2f}"])

# anisotropic resolution (CSM2 static, baseline discretisation)
aniso = []
for name in sorted(glob.glob(work + "/T_nx*")):
    s = static(os.path.basename(name))
    if s:
        nx, ny = [int(x) for x in os.path.basename(name)[4:].replace("ny", "").split("_")]
        aniso.append({"nx": nx, "ny": ny, "uy_mm": s["uy"], "ux_mm": s["ux"]})
base = {r["level"]: r for r in static_rows if r["variant"] == "baseline"}
if 1 in base:
    aniso.append({"nx": NX1, "ny": NY1, "uy_mm": base[1]["uy"], "ux_mm": base[1]["ux"]})
aniso.sort(key=lambda r: (r["nx"], r["ny"]))
with open(os.path.join(out, "csm2_static_anisotropic.csv"), "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["nx", "ny", "ux_mm", "uy_mm", "uy_err_vs_featflow_pct"])
    for r in aniso:
        w.writerow([r["nx"], r["ny"], f"{r['ux_mm']:.5f}", f"{r['uy_mm']:.5f}",
                    f"{100*(r['uy_mm']-REF2['uy'])/abs(REF2['uy']):+.2f}"])

json.dump({"reference_csm3_featflow_dt0.005_L4": REF3, "reference_csm2_featflow_L5+1": REF2,
           "csm3_levels": csm3, "csm3_observed_order": orders3, "csm3_dt_check": dtrows,
           "csm2_static": static_rows, "csm2_static_observed_order": static_orders,
           "csm2_static_anisotropic": aniso},
          open(os.path.join(out, "csm3_structural_convergence.json"), "w"), indent=1,
          default=float)
print("written", out)
