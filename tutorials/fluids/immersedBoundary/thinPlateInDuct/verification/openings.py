#!/usr/bin/env python3
"""Q(6 s)/Q_free of the opening cases. usage: openings.py APERTURE REFINE
Q_free = 6.0e-5 m3/s: the flux of the duct without a plate at 6 s (inlet
p = 0.01, 2 m long, slip walls: dU/dt = 0.005 m/s2, A = 2e-3 m2). There is no
exact reference for the flux through an opening narrower than a cell."""
import sys, glob, numpy as np
ap, ref = sys.argv[1], sys.argv[2]
for d in sorted(glob.glob("../run_slip_0.43_crossedFaces_*_%s_h%s" % (ap, ref))):
    a = np.loadtxt(d + "/postProcessing/fluxOut/0/surfaceFieldValue.dat", comments="#")
    print("%-60s Q(6)/Qfree %.4f" % (d[3:], a[abs(a[:, 0] - 6).argmin(), 1]/6.0e-5))
