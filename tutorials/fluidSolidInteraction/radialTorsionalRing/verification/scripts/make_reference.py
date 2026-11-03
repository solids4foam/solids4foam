#!/usr/bin/env python3
"""Semi-analytical reference solution for the radial-torsional ring (portable).

Solves the finite-amplitude reduced model of rt_reference.py (neo-Hookean solid on
Legendre-Gauss-Lobatto nodes, Chebyshev-collocated fluid, Radau time integration with
rtol 1e-11, 20 solid and 32 fluid nodes) and stores, at every solver time step, the
quantities of interest compared by analyse_case.py in CASE_DIR/reference.json.

The reference is only numerically resolved (spectral refinement and tolerance
sensitivity); it is not validated by a second, independent nonlinear implementation.

Tested with Python 3.12.3, numpy 1.26.4, scipy 1.11.4 (the versions of the accepted
study; other versions change the reference only at the integrator tolerance level).
Run time: about 2 s.

Usage: make_reference.py CASE_DIR [--config torsion|radial|combined] [--mesh M1|M2|M3]
                         [--dt 0.0025] [--tend 2]
"""
import sys as _sys
_sys.dont_write_bytecode = True
import argparse
import json
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import rt_reference as rt  # noqa: E402

NTHETA = {"M1": 64, "M2": 128, "M3": 256}
CONFIGS = {"torsion": (0.1, 0.0), "radial": (0.0, -3000.0), "combined": (0.1, -3000.0)}


def write_reference(C, p, ntheta, dt, tend, args=None):
    """Write C/reference.json for rt.Params p on the time grid k*dt, k = 0 .. round(tend/dt)+2."""
    n = int(round(tend / dt))
    tt = np.arange(n + 3) * dt
    # reference (fine resolution) at every time step of the run
    m, sol = rt.solve(p, tt[-1], 20, 32, rtol=1e-11, t_eval=tt)
    H = rt.history(m, sol)
    keys = ("t", "b", "u_B", "bdot", "alpha_B", "alpha_A", "omega_B", "Q", "p_b", "sig_rr_b", "sig_nv_b", "sig_rt_b",
            "T_B", "T_A", "pc", "H_s", "H_f", "H_flux_out", "KE_s", "V_s", "KE_f", "D_f", "P_A", "P_c_in")
    ref = {k: H[k].tolist() for k in keys}
    # reference solid field at the solid interface face-centre radius (polygon) for u, alpha
    Rf = p.B * math.cos(math.pi / ntheta)
    from numpy.polynomial import legendre as L
    xg = 2 * (m.R - p.A) / (p.B - p.A) - 1
    xf = 2 * (Rf - p.A) / (p.B - p.A) - 1
    uf, af = [], []
    for k in range(tt.size):
        u, al, _, _, _ = m.unpack(sol.y[:, k], tt[k])
        uf.append(float(L.legval(xf, L.legfit(xg, u, 20))))
        af.append(float(L.legval(xf, L.legfit(xg, al, 20))))
    ref["u_Rf"], ref["alpha_Rf"], ref["R_face_solid"] = uf, af, Rf
    os.makedirs(C, exist_ok=True)
    meta = dict(ntheta=ntheta, dt=dt, tend=tend, theta0=p.theta0, P0=p.P0, Omega=p.Omega, Tramp=p.T_ramp)
    with open(os.path.join(C, "reference.json"), "w") as fh:
        json.dump(dict(params=repr(p), args=meta, **ref), fh)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("case")
    ap.add_argument("--config", choices=sorted(CONFIGS), default="combined")
    ap.add_argument("--mesh", choices=sorted(NTHETA), default="M1", help="sets ntheta (polygon face radius of u_Rf, alpha_Rf)")
    ap.add_argument("--dt", type=float, default=0.0025)
    ap.add_argument("--tend", type=float, default=2.0)
    a = ap.parse_args()
    th, P0 = CONFIGS[a.config]
    p = rt.Params(theta0=th, Omega=3.0, P0=P0, T_ramp=1.0, ramp="Cinf")
    write_reference(a.case, p, NTHETA[a.mesh], a.dt, a.tend)
    print("written", os.path.join(a.case, "reference.json"))


if __name__ == "__main__":
    main()
