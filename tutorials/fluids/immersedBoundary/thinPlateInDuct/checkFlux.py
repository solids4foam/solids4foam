#!/usr/bin/env python3
"""Q/(w*A) of the outlet flux, mean over 2 <= t <= 6 s; usage: checkFlux.py CASE MIN MAX"""
import sys, numpy as np
a = np.loadtxt(sys.argv[1] + "/postProcessing/fluxOut/0/surfaceFieldValue.dat", comments="#")
q = a[(a[:, 0] >= 2) & (a[:, 0] <= 6), 1].mean()/(0.1*0.2*0.01)
ok = float(sys.argv[2]) <= q <= float(sys.argv[3])
print("%s: Q/(wA) = %.4f, expected in [%s, %s]: %s" % (sys.argv[1], q, sys.argv[2], sys.argv[3], "ok" if ok else "FAILED"))
sys.exit(0 if ok else 1)
