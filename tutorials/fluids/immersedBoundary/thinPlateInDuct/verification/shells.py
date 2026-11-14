#!/usr/bin/env python3
"""Q/(wA) of the shell cases: mean over 2-8 s and the values at 1, 4, 10 s."""
import glob, numpy as np
for d in sorted(glob.glob("../run_*_0.3_crossedFaces_*_*")):
    a = np.loadtxt(d + "/postProcessing/fluxOut/0/surfaceFieldValue.dat", comments="#")
    q = a[:, 1]/2e-4; t = a[:, 0]
    print("%-50s mean %.4f | %s" % (d[3:], q[(t >= 2) & (t <= 8)].mean(), " ".join("%.3f" % q[abs(t - x).argmin()] for x in (1, 4, 10))))
