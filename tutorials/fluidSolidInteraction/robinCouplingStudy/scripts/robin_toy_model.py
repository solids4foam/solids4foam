#!/usr/bin/env python3
"""
Toy Fourier model of the Robin-Neumann FSI fixed-point iteration used by
solids4foam's elasticWallPressure boundary condition.

Per mode k (acceleration form, all impedances in kg/m^2):

    rho(alpha, k) = (S_f/S_s) |S_s - alpha| / (S_f + alpha)

S_f: fluid added-mass operator, S_s: solid time-discrete interface impedance,
alpha = rho_s*hs: the Robin parameter.

Solid model: through-thickness elastodynamic slab (implicit time step gives a
modified Helmholtz operator with decay length l = c_p*dt_eff), free back face:
    S_slab = rho_s * l * tanh(H/l)
plus a thin-shell stiffness (hoop + bending) contribution dt_eff^2*K(k).

Fluid model:
    tube (axisymmetric, inviscid): S_f = rho_f I0(kR)/(k I1(kR))
    flat wall facing a fluid layer of depth Lf: S_f = rho_f coth(k Lf)/k
"""

import numpy as np


def i0(x):
    return np.i0(x)


def i1(x):
    # Modified Bessel I1 via I1 = dI0/dx (numerical, adequate for a toy model)
    h = 1e-6*np.maximum(np.abs(x), 1.0)
    return (np.i0(x + h) - np.i0(x - h))/(2*h)


def pwave_modulus(E, nu):
    return E*(1 - nu)/((1 + nu)*(1 - 2*nu))


def case_params():
    return {
        "3dTube": dict(
            E=3e5, nu=0.3, rho_s=1200, rho_f=1000, H=1e-3, R=5e-3,
            L=0.05, dt=2.5e-5, ddt="Euler", fluid="tube", two_sided=False,
            dz=0.05/100,
        ),
        "3dTube_dtx16": dict(
            E=3e5, nu=0.3, rho_s=1200, rho_f=1000, H=1e-3, R=5e-3,
            L=0.05, dt=16*2.5e-5, ddt="Euler", fluid="tube", two_sided=False,
            dz=0.05/100,
        ),
        "cerebralAneurysm": dict(
            E=640e3, nu=0.45, rho_s=1200, rho_f=1050, H=3e-4, R=2e-3,
            L=0.03, dt=5e-5, ddt="Euler", fluid="tube", two_sided=False,
            dz=2e-4,
        ),
        "fillingElasticContainer": dict(
            E=21e6, nu=0.3, rho_s=20, rho_f=1000, H=0.2, R=None,
            L=5.0, dt=1e-3, ddt="Euler", fluid="layer", Lf=2.0,
            two_sided=False, dz=0.1,
        ),
        "beamInCrossFlow_modified": dict(
            E=1e4, nu=0.4, rho_s=1000, rho_f=1000, H=0.1, R=None,
            L=1.0, dt=0.05, ddt="backward", fluid="layer", Lf=0.5,
            two_sided=True, dz=0.02,
        ),
        "beamInCrossFlow_original": dict(
            E=1.4e6, nu=0.4, rho_s=1000, rho_f=1000, H=0.1, R=None,
            L=1.0, dt=0.05, ddt="backward", fluid="layer", Lf=0.5,
            two_sided=True, dz=0.02,
        ),
    }


def analyse(name, p, verbose=True):
    M = pwave_modulus(p["E"], p["nu"])
    c = np.sqrt(M/p["rho_s"])
    dt_eff = p["dt"]/1.5 if p["ddt"] == "backward" else p["dt"]
    ell = c*dt_eff
    H_eff = 0.5*p["H"] if p["two_sided"] else p["H"]
    hs_default = ell
    hs_tanh = ell*np.tanh(H_eff/ell)

    # Modes from domain length to grid scale
    k = np.logspace(np.log10(np.pi/p["L"]), np.log10(np.pi/p["dz"]), 400)

    if p["fluid"] == "tube":
        R = p["R"]
        Sf = p["rho_f"]*i0(k*R)/(k*i1(k*R))
        Kh = p["E"]*p["H"]/(R**2*(1 - p["nu"]**2))
    else:
        Sf = p["rho_f"]/(k*np.tanh(k*p["Lf"]))
        Kh = 0.0
        if p["two_sided"]:
            Sf = 0.5*Sf  # per-face share of a two-sided load on one plate

    D = p["E"]*p["H"]**3/(12*(1 - p["nu"]**2))
    Ss = p["rho_s"]*hs_tanh + dt_eff**2*(Kh + D*k**4)

    def rate(alpha):
        return np.max((Sf/Ss)*np.abs(Ss - alpha)/(Sf + alpha))

    hs_grid = np.logspace(np.log10(hs_tanh) - 2, np.log10(hs_default) + 1, 2000)
    rates = np.array([rate(p["rho_s"]*h) for h in hs_grid])
    i_opt = np.argmin(rates)
    hs_opt = hs_grid[i_opt]
    # Largest hs that still contracts
    ok = hs_grid[rates < 1.0]
    hs_max = ok.max() if ok.size else np.nan

    res = dict(
        name=name, c=c, dt_eff=dt_eff, ell=ell, H=p["H"], H_eff=H_eff,
        hs_default=hs_default, hs_tanh=hs_tanh, hs_opt=hs_opt,
        rate_default=rate(p["rho_s"]*hs_default),
        rate_tanh=rate(p["rho_s"]*hs_tanh), rate_opt=rates[i_opt],
        hs_divergence=hs_max,
        addedMassRatio=np.max(Sf/Ss),
    )

    if verbose:
        print(
            f"{name:26s} c={c:8.3g} l={ell:8.3g} H_eff={H_eff:8.3g} "
            f"| hs: default={hs_default:8.3g} tanh={hs_tanh:8.3g} "
            f"opt={hs_opt:8.3g} diverge>{hs_max:8.3g} "
            f"| rate: default={res['rate_default']:6.3f} "
            f"tanh={res['rate_tanh']:6.3f} opt={res['rate_opt']:6.3f} "
            f"| max Sf/Ss={res['addedMassRatio']:8.3g}"
        )
    return res


if __name__ == "__main__":
    for name, p in case_params().items():
        analyse(name, p)
