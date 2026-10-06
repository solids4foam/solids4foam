#!/usr/bin/env python3
"""CSM3 statistics for tip point A from solidPointDisplacement output.

usage: analyse.py <caseDir>... ; prints mean, amplitude, frequency of ux, uy
over the closing window [T-1, T] (contains >= 1 full period, 0.91 s) (extrema refined by a parabola).
"""
import sys, glob, json
import numpy as np

def extremum(t, y, kind):
    i = np.argmax(y) if kind == "max" else np.argmin(y)
    if 0 < i < len(y) - 1:
        a, b, c = np.polyfit(t[i-1:i+2] - t[i], y[i-1:i+2], 2)
        if a != 0:
            return c - b*b/(4*a)
    return y[i]

def stats(t, y, tfull, yfull):
    mx, mn = extremum(t, y, "max"), extremum(t, y, "min")
    # frequency from upward crossings of the window mean over the whole run
    m = 0.5*(mx+mn)
    up = [tfull[i] + (m-yfull[i])*(tfull[i+1]-tfull[i])/(yfull[i+1]-yfull[i])
          for i in range(len(yfull)-1) if yfull[i] < m <= yfull[i+1]]
    f = (len(up)-1)/(up[-1]-up[0]) if len(up) > 1 else float("nan")
    return 0.5*(mx+mn), 0.5*(mx-mn), f

def analyse(case, window=1.0):
    f = glob.glob(f"{case}/postProcessing/*/solidPointDisplacement_pointDisp.dat")[0]
    d = np.loadtxt(f, comments="#")
    t = d[:, 0]; T = t[-1]
    sel = t >= T - window - 1e-12
    out = {"case": case.rstrip("/").split("/")[-1], "T": T}
    for name, col in (("ux", 1), ("uy", 2)):
        m, a, fr = stats(t[sel], d[sel, col], t, d[:, col])
        out[name] = {"mean": m, "amp": a, "freq": fr}
    try:
        out["runtime_s"] = int(open(f"{case}/runtime.s").read())
    except OSError:
        pass
    return out

if __name__ == "__main__":
    res = [analyse(c) for c in sys.argv[1:]]
    print(json.dumps(res, indent=1))
