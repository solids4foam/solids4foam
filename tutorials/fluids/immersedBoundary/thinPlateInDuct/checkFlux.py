#!/usr/bin/env python3
"""Q/(w*A) of the outlet flux, time-average over 2 <= t <= 8 s (mean of the
per-time-step values, which are not equally spaced when the time step is
adjusted); usage: checkFlux.py CASE [MIN MAX]"""
import sys, numpy as np
a = np.loadtxt(sys.argv[1] + "/postProcessing/fluxOut/0/surfaceFieldValue.dat", comments="#")
q = a[(a[:, 0] >= 2) & (a[:, 0] <= 8), 1]/(0.1*0.2*0.01)
print("%s: Q/(wA) mean over 2-8 s = %.4f (min %.4f, max %.4f)" % (sys.argv[1], q.mean(), q.min(), q.max()))
if len(sys.argv) > 3:
    ok = float(sys.argv[2]) <= q.mean() <= float(sys.argv[3])
    print("  expected in [%s, %s]: %s" % (sys.argv[2], sys.argv[3], "ok" if ok else "FAILED"))
    sys.exit(0 if ok else 1)
