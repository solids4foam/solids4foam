#!/usr/bin/env python3
"""Q_in and plate force harmonics in sliding 5-period windows (drift check)."""
import sys, math, json
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
import energy_analysis as ea

def slide(case, periods=5, step_periods=1.0, tstart=None):
    s = ea.setup(case); k = float(s["amplitude_factor"])
    data = np.loadtxt(case / "energy.dat")
    _, idx = np.unique(data[:, 0], return_index=True); data = data[np.sort(idx)]
    col = {n: data[:, i] for i, n in enumerate(ea.COLS)}
    t = col["t"]
    omega, d_tip, *_ = ea.plate_kinematics(case, k)
    P = 2*math.pi/omega
    out = []
    t1 = t[-1]
    t0min = (tstart if tstart is not None else t[0]+0.2)
    while t1 - periods*P >= t0min:
        sel = (t > t1-periods*P) & (t <= t1); n = sel.sum()
        ph = np.exp(-1j*omega*col["T"][sel])
        g = 2.0/n*np.sum((col["Q1ya"][sel]+1j*col["Q1yb"][sel])*ph)
        fy = col["Fpy"][sel]+col["Fvy"][sel]; fx = col["Fpx"][sel]+col["Fvx"][sel]
        L = 2.0/n*np.sum(fy*ph)
        out.append(dict(t_end=float(t1), Q_in=float(g.real/abs(d_tip)), Q_quad=float(g.imag/abs(d_tip)),
                        lift_amp1=float(abs(L)), drag_mean=float(fx.mean())))
        t1 -= step_periods*P
    return out[::-1]

if __name__ == "__main__":
    for c in sys.argv[1:]:
        print("==", c)
        for r in slide(Path(c)):
            print("  t_end %.3f  Q_in %.3f  Q_quad %.3f  lift1 %.2f  drag %.3f" % (r["t_end"], r["Q_in"], r["Q_quad"], r["lift_amp1"], r["drag_mean"]))
