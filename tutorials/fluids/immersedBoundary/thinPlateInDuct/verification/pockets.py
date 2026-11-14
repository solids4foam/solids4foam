#!/usr/bin/env python3
"""Pressure metrics of pocketCase runs. usage: pockets.py T_EVENT CASEDIR...
Windows: before (1 s .. T_EVENT-0.3), early (0.1-0.8 s), event (T_EVENT +- 0.15),
after (T_EVENT+0.5 .. +2.5). p_max: largest |p| of the field; dpAB: largest
|pA - pB| of the area-averaged pressure at x = 1.5 m and 1.7 m (a difference
within the downstream region, independent of the gauge of a pocket). The
estimate L dU/dt of an incompressible column of length L = 1 m (B to outlet)
with a one-step jump of the outlet velocity is printed for the event window."""
import sys, numpy as np
te = float(sys.argv[1])
for d in sys.argv[2:]:
    t, lo, hi = [], [], []
    for l in open(d + "/postProcessing/pMM/0/fieldMinMax.dat"):
        if l.startswith("#"):
            continue
        f = l.replace("(", "").replace(")", "").split()
        t.append(float(f[0])); lo.append(float(f[2])); hi.append(float(f[6]))
    t, lo, hi = map(np.array, (t, lo, hi))
    pa = np.loadtxt(d + "/postProcessing/p_fz_a/0/surfaceFieldValue.dat", comments="#")
    pb = np.loadtxt(d + "/postProcessing/p_fz_b/0/surfaceFieldValue.dat", comments="#")
    n = min(len(pa), len(pb)); tt = pa[:n, 0]; dp = abs(pa[:n, 1] - pb[:n, 1])
    out = "%s steps %d " % (d, len(t))
    for lab, (a, b) in (("before", (1, te - 0.3)), ("early", (0.1, 0.8)), ("event", (te - 0.15, te + 0.15)), ("after", (te + 0.5, te + 2.5))):
        m = (t >= a) & (t <= b); m2 = (tt >= a) & (tt <= b)
        if m.sum() and m2.sum():
            out += "| %s p_max %.3g dpAB %.3g " % (lab, max(hi[m].max(), -lo[m].min()), dp[m2].max())
    print(out)
