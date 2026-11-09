#!/usr/bin/env python3
"""Metrics of gateCase runs. usage: gate.py CASEDIR...
flux error: max |Q_out/Q_in - 1| for t >= 0.5 s (Q_in = 1e-4 m3/s, the flux is
prescribed, so this is the volume conservation of the run).
dp(t) = area-averaged p at x = 0.15 m minus that at x = 1.5 m.
ripple: dp(t) is interpolated to 400 uniform times in [1, 8] s, log|dp| is
fitted with a polynomial of degree 6 and the ripple is
max |dp - sign(dp) exp(fit)| / max |dp| (it is large when dp changes sign or
is not smooth in time; smooth dp gives a value near 0). dp_end: dp at 8 s."""
import sys, numpy as np
Qin = 0.05*0.2*0.01
for d in sys.argv[1:]:
    ld = lambda n: np.loadtxt(d + "/postProcessing/%s/0/surfaceFieldValue.dat" % n, comments="#")
    up, dn, fo = ld("pUp"), ld("pDn"), ld("fluxOut")
    n = min(len(up), len(dn)); t = up[:n, 0]; dp = up[:n, 1] - dn[:n, 1]
    tg = np.linspace(1, 8, 400); dpg = np.interp(tg, t, dp)
    c = np.polyfit(tg, np.log(np.maximum(abs(dpg), 1e-12)), 6)
    rip = np.abs(dpg - np.sign(dpg)*np.exp(np.polyval(c, tg))).max()/abs(dpg).max()
    m = fo[:, 0] >= 0.5
    print("%-24s steps %5d flux error %.1e ripple %.3f max|dp| %.3g dp_end %.3g" % (d, len(t), np.abs(abs(fo[m, 1])/Qin - 1).max(), rip, abs(dpg).max(), dp[-1]))
