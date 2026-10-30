#!/usr/bin/env python3
"""Finite-amplitude reduced reference for the breathing-and-twisting ring.

Problem (plane strain, rotationally symmetric, per unit depth, per radian):

  solid  A <= R <= B, motion x = phi(R,t) e_r(Theta + alpha(R,t)),
         compressible neo-Hookean (solids4foam form) or StVK, rho_s;
         R = A: clamped radially (u = 0), prescribed rotation alpha_A(t)
                (solids4foam: fixedRotation);
         R = B: fluid-solid interface.
  fluid  b(t) = phi(B,t) <= r <= c, incompressible Newtonian rho_f, mu_f,
         v = Q/r e_r + (L/r) e_theta, Q = b bdot (exact),
         L_t + (Q/r) L_r = nu (L_rr - L_r/r)  on the moving gap,
         L(b) = b^2 alpha_t(B)  (no slip),
         r = c: p = p_c(t), d(u_r)/dr = -u_r/c (exact), sigma_rt = 0  (open,
                tangentially traction-free);  =>  L_r(c) = 2 L(c)/c.
  interface:  B P_rR = b sigma_rr(b),   B phi P_tR = b^2 sigma_rt(b)
         sigma_rr(b) = -p(b) - 2 mu Q/b^2,  sigma_rt = mu (L_r/r - 2L/r^2),
         p(b) = p_c + rho[Qdot ln(c/b) + Q^2/2 (1/c^2 - 1/b^2)] - rho int_b^c L^2/r^3 dr.

Discretisation (independent of the earlier radial-only reference):
  solid: discrete Lagrangian on Legendre-Gauss-Lobatto nodes (Galerkin with GLL
         quadrature: diagonal mass, exact discrete rotational invariance);
  fluid: Chebyshev collocation in xi = (r - b)/(c - b);
  time:  scipy solve_ivp (Radau / DOP853 / BDF), dense output.

Units SI.  All energies, torques, angular momenta are PER RADIAN of the full
circle and per unit depth (multiply by 2 pi for the whole annulus).
"""
from __future__ import annotations

import dataclasses
import math
import time as _time

import numpy as np
from numpy.polynomial import legendre as npleg
from scipy.integrate import solve_ivp


# ----------------------------------------------------------------------------
# parameters and drives
# ----------------------------------------------------------------------------
@dataclasses.dataclass
class Params:
    A: float = 0.5
    B: float = 1.0
    c: float = 2.0
    rho_s: float = 1000.0
    E: float = 2.0e4
    nu_s: float = 0.3
    rho_f: float = 1000.0
    mu_f: float = 300.0
    law: str = "neoHookean"           # "neoHookean" | "StVK" | "linear"
    # drives
    theta0: float = 0.0               # rotation amplitude of the inner surface [rad]
    Omega: float = 8.0                # rotation frequency [rad/s]
    T_ramp: float = 0.5               # C3 ramp duration [s]
    P0: float = 0.0                   # outer-boundary pressure amplitude [Pa]
    torsion_mode: str = "ramped_sine"  # "ramped_sine" | "ramp_hold"
    ramp: str = "C3"                  # "C3" (smoothstep7) | "Cinf" (exp-based, all derivatives vanish)
    pressure_mode: str = "ramp_hold"   # "ramp_hold" | "ramped_sine"
    Omega_p: float = 8.0
    # model switches (for sensitivity studies only)
    f_normal_viscous: float = 1.0     # fraction of 2 mu d(u_r)/dr transmitted
    fluid_convection: bool = True
    # fluid=False: solid-only variant (Test A): outer surface R = B carries the
    # follower pressure p_c(t) on the current area and zero shear; no fluid
    fluid: bool = True
    # emulated implementation defects (sensitivity studies only; "" = exact model):
    #  "no_L_convection"   fluid: radial convection of angular momentum (Q/r) L_r omitted
    #  "no_centrifugal_p"  fluid: centrifugal pressure  -rho int w^2/r dr  omitted
    #  "reference_area"    interface tractions applied per REFERENCE area (B instead of b)
    #  "unrotated_traction" interface traction applied on the reference basis (E_R, E_T)
    #                       instead of the current basis (e_r, e_t) (missing normal rotation)
    #  "no_solid_coriolis" solid: Coriolis 2 phi_t alpha_t and centripetal phi alpha_t^2 omitted
    #  "normal_viscous_1_3" only 1/3 of 2 mu d(u_r)/dr transmitted (code-path reading, cf. Womersley)
    defect: str = ""

    @property
    def G(self):
        return self.E / (2 * (1 + self.nu_s))

    @property
    def lam(self):                     # plane strain uses the 3-D lambda
        return self.E * self.nu_s / ((1 + self.nu_s) * (1 - 2 * self.nu_s))

    @property
    def kappa(self):
        return self.lam + 2.0 / 3.0 * self.G

    @property
    def nu_f(self):
        return self.mu_f / self.rho_f


def smoothstep7(tau):
    """C3 step: 0 for tau<=0, 1 for tau>=1; S, S', S'', S''' vanish at both ends."""
    tau = np.clip(tau, 0.0, 1.0)
    return tau**4 * (35 - 84 * tau + 70 * tau**2 - 20 * tau**3)


def smoothstep7_d(tau, k):
    """k-th derivative (k = 1, 2) of smoothstep7 with respect to tau (0 outside)."""
    inside = (tau > 0) & (tau < 1) if np.ndim(tau) else (0 < tau < 1)
    t = np.clip(tau, 0.0, 1.0)
    if k == 1:
        d = 140 * t**3 - 420 * t**4 + 420 * t**5 - 140 * t**6
    elif k == 2:
        d = 420 * t**2 - 1680 * t**3 + 2100 * t**4 - 840 * t**5
    else:
        raise ValueError
    return np.where(inside, d, 0.0) if np.ndim(tau) else (d if inside else 0.0)


def _psi(x):
    x = np.asarray(x, dtype=float)
    with np.errstate(divide="ignore", over="ignore", invalid="ignore"):
        return np.where(x > 0, np.exp(-1.0 / np.where(x > 0, x, 1.0)), 0.0)


def ramp_cinf(tau):
    """C-infinity step S and its first two tau-derivatives (all derivatives vanish at 0 and 1)."""
    tau = float(tau)
    if tau <= 0.0:
        return 0.0, 0.0, 0.0
    if tau >= 1.0:
        return 1.0, 0.0, 0.0
    x, y = tau, 1.0 - tau
    f, g = math.exp(-1 / x), math.exp(-1 / y)
    f1, g1 = f / x**2, -g / y**2                      # d/dtau
    f2, g2 = f * (1 / x**4 - 2 / x**3), g * (1 / y**4 - 2 / y**3)
    D = f + g
    S = f / D
    N = f1 * g - f * g1
    S1 = N / D**2
    N1 = f2 * g - f * g2
    D1 = f1 + g1
    S2 = (N1 * D - 2 * N * D1) / D**3
    return S, S1, S2


def ramp(p, t):
    """S(t/T_ramp) and its first two time derivatives."""
    Tr = p.T_ramp
    tau = t / Tr
    if p.ramp == "Cinf":
        S, S1, S2 = ramp_cinf(tau)
    else:
        S, S1, S2 = float(smoothstep7(tau)), float(smoothstep7_d(tau, 1)), float(smoothstep7_d(tau, 2))
    return S, S1 / Tr, S2 / Tr**2


def rotation_drive(p: Params, t):
    """alpha_A(t), its first and second time derivatives."""
    if p.theta0 == 0.0:
        return 0.0, 0.0, 0.0
    S, S1, S2 = ramp(p, t)
    if p.torsion_mode == "ramp_hold":
        return p.theta0 * S, p.theta0 * S1, p.theta0 * S2
    W = p.Omega
    s, co = math.sin(W * t), math.cos(W * t)
    a = p.theta0 * S * s
    ad = p.theta0 * (S1 * s + S * W * co)
    add = p.theta0 * (S2 * s + 2 * S1 * W * co - S * W * W * s)
    return a, ad, add


def pressure_drive(p: Params, t):
    if p.P0 == 0.0:
        return 0.0
    S = ramp(p, t)[0]
    if p.pressure_mode == "ramp_hold":
        return p.P0 * S
    return p.P0 * S * math.sin(p.Omega_p * t)


# ----------------------------------------------------------------------------
# constitutive laws: cancellation-free in the small quantities
#   a = phi' = 1 + e1,  h = phi/R = 1 + e2,  g = phi alpha'
# returns (W, P_rR = dW/da, P_tR = dW/dg, P_tT = dW/dh)
# ----------------------------------------------------------------------------
def law_forces(law, G, lam, kappa, e1, e2, g, want_energy=False):
    a = 1.0 + e1
    h = 1.0 + e2
    a2m1 = e1 * (2.0 + e1)
    h2m1 = e2 * (2.0 + e2)
    g2 = g * g
    if law == "neoHookean":
        Jm1 = e1 + e2 + e1 * e2
        J = 1.0 + Jm1
        J23 = np.exp(-2.0 / 3.0 * np.log1p(Jm1))         # J^(-2/3)
        hyd = 0.5 * kappa * Jm1 * (J + 1.0)              # K/2 (J^2 - 1)
        t_rr = G * J23 * (2 * a2m1 - h2m1 - g2) / 3.0 + hyd
        t_tt = G * J23 * (2 * h2m1 + 2 * g2 - a2m1) / 3.0 + hyd
        t_rt = G * J23 * a * g
        PrR = t_rr / a
        PtR = t_rt / a
        PtT = (t_tt - t_rt * g / a) / h
        if want_energy:
            # mu/2 (J^-2/3 tr C - 3) + K/4 (J^2 - 1 - 2 ln J)
            # Ibar1 - 3 without O(eps) cancellation: with log-stretches x1 = ln a,
            # x2 = ln h, x3 = 0 and y_i = x_i - (x1 + x2)/3 (sum y_i = 0),
            # Ibar1 - 3 = sum_i E(2 y_i) + J^(-2/3) g^2,  E(z) = e^z - 1 - z
            x1, x2 = np.log1p(e1), np.log1p(e2)
            mbar = (x1 + x2) / 3.0
            Ibar1m3 = _E(2 * (x1 - mbar)) + _E(2 * (x2 - mbar)) + _E(-2 * mbar) + J23 * g2
            Wiso = 0.5 * G * Ibar1m3
            lnJ = np.log1p(Jm1)
            # J^2 - 1 - 2 ln J = 2(J-1) + (J-1)^2 - 2 ln J ; 2[(J-1) - ln J] via series-safe form
            x = Jm1
            d = np.where(np.abs(x) < 1e-3,
                         x**2 / 2 - x**3 / 3 + x**4 / 4 - x**5 / 5 + x**6 / 6,
                         x - lnJ)
            Wvol = 0.25 * kappa * (2.0 * d + x * x)
            return Wiso + Wvol, PrR, PtR, PtT
        return None, PrR, PtR, PtT
    if law == "StVK":
        ERR = 0.5 * (a2m1 + g2)
        ERT = 0.5 * g * h
        ETT = 0.5 * h2m1
        trE = ERR + ETT
        SRR = lam * trE + 2 * G * ERR
        SRT = 2 * G * ERT
        STT = lam * trE + 2 * G * ETT
        PrR = a * SRR
        PtR = g * SRR + h * SRT
        PtT = g * SRT + h * STT
        if want_energy:
            Wv = 0.5 * lam * trE**2 + G * (ERR**2 + 2 * ERT**2 + ETT**2)
            return Wv, PrR, PtR, PtT
        return None, PrR, PtR, PtT
    if law == "linear":
        # small-strain linear elasticity written on the same generalised
        # coordinates (u_r = phi - R, u_t = R alpha): only for linear checks
        raise NotImplementedError
    raise ValueError(law)


def _E(z):
    """e^z - 1 - z, accurate for small |z| (series below 1e-2)."""
    z = np.asarray(z, dtype=float)
    small = np.abs(z) < 1e-2
    ser = z * z * (0.5 + z * (1 / 6 + z * (1 / 24 + z * (1 / 120 + z * (1 / 720 + z / 5040)))))
    return np.where(small, ser, np.expm1(z) - z)


# ----------------------------------------------------------------------------
# spectral building blocks
# ----------------------------------------------------------------------------
def gll(N):
    """Legendre-Gauss-Lobatto nodes (ascending), weights, differentiation matrix on [-1,1]."""
    PN = npleg.Legendre.basis(N)
    x = np.concatenate(([-1.0], np.sort(PN.deriv().roots().real), [1.0]))
    Pv = PN(x)
    w = 2.0 / (N * (N + 1) * Pv**2)
    D = np.zeros((N + 1, N + 1))
    for i in range(N + 1):
        for j in range(N + 1):
            if i != j:
                D[i, j] = Pv[i] / (Pv[j] * (x[i] - x[j]))
    D[0, 0] = -N * (N + 1) / 4.0
    D[N, N] = N * (N + 1) / 4.0
    return x, w, D


def cheb(N):
    """Chebyshev-Gauss-Lobatto nodes x_k = -cos(pi k/N) (ascending), D, Clenshaw-Curtis weights."""
    k = np.arange(N + 1)
    x = -np.cos(np.pi * k / N)
    cc = np.ones(N + 1)
    cc[0] = cc[-1] = 2.0
    cc = cc * (-1.0) ** k
    X = np.tile(x, (N + 1, 1)).T
    dX = X - X.T
    D = np.outer(cc, 1.0 / cc) / (dX + np.eye(N + 1))
    D = D - np.diag(D.sum(axis=1))
    # Clenshaw-Curtis weights on [-1, 1]
    theta = np.pi * k / N
    w = np.zeros(N + 1)
    v = np.ones(N - 1)
    if N % 2 == 0:
        w[0] = w[N] = 1.0 / (N**2 - 1)
        for kk in range(1, N // 2):
            v -= 2 * np.cos(2 * kk * theta[1:-1]) / (4 * kk**2 - 1)
        v -= np.cos(N * theta[1:-1]) / (N**2 - 1)
    else:
        w[0] = w[N] = 1.0 / N**2
        for kk in range(1, (N - 1) // 2 + 1):
            v -= 2 * np.cos(2 * kk * theta[1:-1]) / (4 * kk**2 - 1)
    w[1:-1] = 2 * v / N
    return x, D, w


# ----------------------------------------------------------------------------
# the semi-discrete system
# ----------------------------------------------------------------------------
class Model:
    def __init__(self, p: Params, Ns=16, Nf=24):
        self.p = p
        self.Ns, self.Nf = Ns, Nf
        x, w, D = gll(Ns)
        hR = 0.5 * (p.B - p.A)
        self.R = p.A + hR * (x + 1.0)
        self.wq = w * hR
        self.Ds = D / hR
        self.Ms = p.rho_s * self.wq * self.R          # nodal masses per radian
        xi, Dc, wc = cheb(Nf)
        self.xi = 0.5 * (xi + 1.0)                    # 0 at interface, 1 at r = c
        self.Dxi = 2.0 * Dc
        self.wxi = 0.5 * wc
        self.n = 4 * Ns + (Nf - 1)
        self.nfev = 0

    # state layout: u[1..Ns], al[1..Ns], v[1..Ns], om[1..Ns], L[1..Nf-1]
    def unpack(self, y, t):
        Ns, Nf = self.Ns, self.Nf
        u = np.empty(Ns + 1); al = np.empty(Ns + 1); v = np.empty(Ns + 1); om = np.empty(Ns + 1)
        u[0] = 0.0; v[0] = 0.0
        aA, adA, _ = rotation_drive(self.p, t)
        al[0] = aA; om[0] = adA
        u[1:] = y[0:Ns]; al[1:] = y[Ns:2 * Ns]; v[1:] = y[2 * Ns:3 * Ns]; om[1:] = y[3 * Ns:4 * Ns]
        Li = y[4 * Ns:]
        return u, al, v, om, Li

    def fluid_L(self, b, om_B, Li):
        """Full nodal L including the boundary values (Dirichlet at b, Robin at c)."""
        p = self.p
        Nf = self.Nf
        ell = p.c - b
        L = np.empty(Nf + 1)
        L[0] = b * b * om_B
        L[1:Nf] = Li
        row = self.Dxi[Nf] / ell
        # row . L = 2 L_N / c  ->  L_N (row_N - 2/c) = - row[:N] . L[:N]
        L[Nf] = -np.dot(row[:Nf], L[:Nf]) / (row[Nf] - 2.0 / p.c)
        return L, ell

    def solid_forces(self, u, al, want_energy=False):
        p = self.p
        R, Ds, wq = self.R, self.Ds, self.wq
        phi = R + u
        e1 = Ds @ u
        e2 = u / R
        alp = Ds @ al
        g = phi * alp
        Wd, PrR, PtR, PtT = law_forces(p.law, p.G, p.lam, p.kappa, e1, e2, g, want_energy)
        wR = wq * R
        # dV/dphi_j = sum_q wR_q [PrR_q D_qj + (PtT_q/R_q + PtR_q alpha'_q) delta_qj]
        Gphi = Ds.T @ (wR * PrR) + wR * (PtT / R + PtR * alp)
        # dV/dalpha_j = sum_q wR_q PtR_q phi_q D_qj
        Gal = Ds.T @ (wR * PtR * phi)
        V = np.sum(wR * Wd) if want_energy else None
        return phi, Gphi, Gal, V, PrR, PtR

    def interface(self, t, b, vB, L, ell, vBdot=None):
        """Interface loads.  Returns F_r0 (radial force per radian excluding the added-mass
        term), m_add (added mass per radian), torque T_B = b^2 sigma_rt, and extras."""
        p = self.p
        if not p.fluid:
            pc = pressure_drive(p, t)
            z = np.zeros_like(L)
            return -b * pc, 0.0, 0.0, dict(Q=b * vB, Lr=z, r=b + self.xi * ell, Iw=0.0, pc=pc,
                                           p_b_wo=pc, sig_nv=0.0, sig_rt=0.0)
        rho, mu = p.rho_f, p.mu_f
        c = p.c
        Q = b * vB
        r = b + self.xi * ell
        Lr = (self.Dxi @ L) / ell
        # int_b^c L^2/r^3 dr
        Iw = ell * np.dot(self.wxi, L**2 / r**3)
        pc = pressure_drive(p, t)
        # np.log: an invalid explicit trial stage (b >= c) must give NaN (step rejected), not raise
        with np.errstate(invalid="ignore", divide="ignore"):
            lnr = float(np.log(c / b)) if c / b > 0 else float("nan")
        conv = 1.0 if p.fluid_convection else 0.0
        # p(b) = pc + rho[ Qdot ln(c/b) + Q^2/2 (1/c^2 - 1/b^2) ] - rho Iw ; Qdot = vB^2 + b vBdot
        cen = 0.0 if p.defect == "no_centrifugal_p" else 1.0
        p_b_wo = pc + rho * (conv * vB * vB * lnr + conv * 0.5 * Q * Q * (1.0 / c**2 - 1.0 / b**2)) - rho * conv * cen * Iw
        fnv = p.f_normal_viscous * (1.0 / 3.0 if p.defect == "normal_viscous_1_3" else 1.0)
        sig_nv = -2.0 * mu * Q / b**2 * fnv
        F_r0 = b * (-p_b_wo + sig_nv)
        m_add = rho * b * b * lnr
        sig_rt = mu * (Lr[0] / b - 2.0 * L[0] / b**2)
        T_B = b * b * sig_rt
        if p.defect == "reference_area":
            F_r0 *= p.B / b
            m_add *= p.B / b
            T_B *= p.B / b
        return F_r0, m_add, T_B, dict(Q=Q, Lr=Lr, r=r, Iw=Iw, pc=pc, p_b_wo=p_b_wo,
                                      sig_nv=sig_nv, sig_rt=sig_rt)

    def rhs(self, t, y):
        self.nfev += 1
        p = self.p
        Ns = self.Ns
        u, al, v, om, Li = self.unpack(y, t)
        phi, Gphi, Gal, _, _, _ = self.solid_forces(u, al)
        b, vB, omB = phi[-1], v[-1], om[-1]
        L, ell = self.fluid_L(b, omB, Li)
        F_r0, m_add, T_B, ex = self.interface(t, b, vB, L, ell)
        M = self.Ms
        kin = 0.0 if p.defect == "no_solid_coriolis" else 1.0
        if p.defect == "unrotated_traction":
            # traction b(s_rr E_R + s_rt E_T) with E_R = cos(a) e_r - sin(a) e_t,
            # E_T = sin(a) e_r + cos(a) e_t, a = alpha_B (normal rotation ignored)
            aB = al[-1]
            s_rr_force, s_rt_force = F_r0, T_B / b          # per radian, excl. added mass
            F_r0 = s_rr_force * math.cos(aB) + s_rt_force * math.sin(aB)
            T_B = b * (-s_rr_force * math.sin(aB) + s_rt_force * math.cos(aB))
        # radial: M (vdot - phi om^2) = -Gphi + delta_N F_r ; F_r = F_r0 - m_add vdot_N
        vdot = (kin * M[1:] * phi[1:] * om[1:] ** 2 - Gphi[1:]) / M[1:]
        vdot[-1] = (kin * M[-1] * phi[-1] * om[-1] ** 2 - Gphi[-1] + F_r0) / (M[-1] + m_add)
        # azimuthal: M (phi^2 omdot + 2 phi v om) = -Gal + delta_N T_B
        rhs_t = -Gal[1:] - kin * 2.0 * M[1:] * phi[1:] * v[1:] * om[1:]
        rhs_t[-1] += T_B
        omdot = rhs_t / (M[1:] * phi[1:] ** 2)
        # fluid interior
        r = ex["r"]
        Lr = ex["Lr"]
        Lrr = (self.Dxi @ Lr) / ell
        Q = ex["Q"]
        conv = 1.0 if p.fluid_convection else 0.0
        mesh_vel = vB * (1.0 - self.xi)
        lconv = 0.0 if p.defect == "no_L_convection" else conv
        Ldot = (-(lconv * Q / r) + mesh_vel) * Lr + p.nu_f * (Lrr - Lr / r)
        if not p.fluid:
            Ldot = np.zeros_like(Ldot)
        return np.concatenate((v[1:], om[1:], vdot, omdot, Ldot[1:self.Nf]))

    # ------------------------------------------------------------------ outputs
    def diagnostics(self, t, y, ydot=None):
        """Physical quantities at time t (per radian, per unit depth)."""
        p = self.p
        u, al, v, om, Li = self.unpack(y, t)
        phi, Gphi, Gal, V, PrR, PtR = self.solid_forces(u, al, want_energy=True)
        b, vB, omB = phi[-1], v[-1], om[-1]
        L, ell = self.fluid_L(b, omB, Li)
        F_r0, m_add, T_B, ex = self.interface(t, b, vB, L, ell)
        if ydot is None:
            ydot = self.rhs(t, y)
        Ns = self.Ns
        vdotB = ydot[3 * Ns - 1]
        omdot = np.concatenate(([rotation_drive(p, t)[2]], ydot[3 * Ns:4 * Ns]))
        vdot = np.concatenate(([0.0], ydot[2 * Ns:3 * Ns]))
        M = self.Ms
        # reactions at the clamped inner node (what the drive must supply)
        T_A = M[0] * (phi[0] ** 2 * omdot[0] + 2 * phi[0] * v[0] * om[0]) + Gal[0]
        F_A = M[0] * (vdot[0] - phi[0] * om[0] ** 2) + Gphi[0]
        rho, mu, c = p.rho_f, p.mu_f, p.c
        Q = ex["Q"]
        r = ex["r"]
        Lr = ex["Lr"]
        F_r = F_r0 - m_add * vdotB
        sig_rr_b = F_r / b
        p_b = ex["p_b_wo"] + rho * b * vdotB * math.log(c / b)
        KE_s = 0.5 * np.sum(M * (v**2 + phi**2 * om**2))
        KE_f = 0.5 * rho * (Q * Q * math.log(c / b) + ell * np.dot(self.wxi, L**2 / r))
        D_f = mu * (2 * Q * Q * (1 / b**2 - 1 / c**2) * p.f_normal_viscous
                    + ell * np.dot(self.wxi, (Lr / r - 2 * L / r**2) ** 2 * r))
        pc = ex["pc"]
        sig_rr_c = -pc - 2 * mu * Q / c**2
        w_c = L[-1] / c
        P_c_in = Q * sig_rr_c - rho * Q * (Q * Q / c**2 + w_c**2) / 2
        P_A = T_A * om[0] + F_A * v[0]
        H_s = np.sum(M * phi**2 * om)
        H_f = rho * ell * np.dot(self.wxi, L * r)
        H_flux_out = rho * Q * L[-1]
        if not p.fluid:
            KE_f = D_f = H_f = H_flux_out = w_c = 0.0
            P_c_in = -pc * b * vB          # follower pressure power on the solid
            p_b = pc
        return dict(t=t, b=b, u_B=b - p.B, bdot=vB, alpha_B=al[-1], alpha_A=al[0], omega_B=omB,
                    Q=Q, p_b=p_b, sig_rr_b=sig_rr_b, sig_nv_b=ex["sig_nv"], sig_rt_b=ex["sig_rt"],
                    T_B=T_B, T_A=T_A, F_A=F_A, w_c=w_c, L_c=L[-1], pc=pc,
                    KE_s=KE_s, V_s=V, KE_f=KE_f, D_f=D_f, P_A=P_A, P_c_in=P_c_in,
                    H_s=H_s, H_f=H_f, H_flux_out=H_flux_out,
                    ang_mom_rate_in=T_A - H_flux_out)


def solve(p: Params, t_end, Ns=16, Nf=24, method="Radau", rtol=1e-10, atol_scale=1e-12,
          t_eval=None, max_wall=600.0, dense=False, first_step=None):
    """Integrate from rest.  Returns (model, solution).  atol is scaled per component."""
    m = Model(p, Ns, Nf)
    y0 = np.zeros(m.n)
    # absolute tolerances scaled by the expected amplitudes
    th = max(abs(p.theta0), 1e-30)
    Pscale = max(abs(p.P0) / (p.lam + 2 * p.G), 1e-30) * (p.B - p.A)
    ua = max(Pscale, th * 1e-2 * p.B, 1e-30)
    Ls = max(th * p.Omega * p.B**2, ua * 10, 1e-30)
    atol = np.concatenate((np.full(m.Ns, ua), np.full(m.Ns, max(th, ua)), np.full(m.Ns, ua * 10),
                           np.full(m.Ns, max(th, ua) * 10), np.full(m.Nf - 1, Ls))) * atol_scale / rtol * rtol
    t0 = _time.time()

    def guarded(t, y):
        if _time.time() - t0 > max_wall:
            raise RuntimeError(f"wall-clock limit {max_wall}s exceeded at t={t}")
        return m.rhs(t, y)

    kw = dict(rtol=rtol, atol=atol, dense_output=dense, t_eval=t_eval)
    if first_step:
        kw["first_step"] = first_step
    sol = solve_ivp(guarded, (0.0, t_end), y0, method=method, **kw)
    sol.wall = _time.time() - t0
    sol.nfev_model = m.nfev
    if not sol.success:
        raise RuntimeError(sol.message)
    return m, sol


def history(m: Model, sol, keys=None):
    """Diagnostics at every stored time; returns dict of arrays."""
    out = None
    for k in range(sol.t.size):
        d = m.diagnostics(sol.t[k], sol.y[:, k])
        if out is None:
            out = {kk: [] for kk in d}
        for kk, vv in d.items():
            out[kk].append(vv)
    return {kk: np.array(vv) for kk, vv in out.items()}


def solve_with_budgets(p: Params, t_end, Ns=16, Nf=24, method="Radau", rtol=1e-11, t_eval=None,
                       max_wall=1200.0):
    """Integrate the system together with the time integrals of the energy and
    angular-momentum sources, so the balance residuals are integrator-controlled:
        E(t) - E(0) - int (P_A + P_c_in - D_f) dt         (energy)
        H(t) - H(0) - int (T_A - rho Q L(c)) dt            (angular momentum)
    Returns (model, sol, history) with history['res_E'], history['res_H']."""
    m = Model(p, Ns, Nf)
    n = m.n
    t0 = _time.time()

    def f(t, z):
        if _time.time() - t0 > max_wall:
            raise RuntimeError("wall-clock limit")
        y = z[:n]
        yd = m.rhs(t, y)
        d = m.diagnostics(t, y, yd)
        return np.concatenate((yd, [d["P_A"] + d["P_c_in"] - d["D_f"], d["ang_mom_rate_in"], d["D_f"]]))

    th = max(abs(p.theta0), 1e-30)
    us = max(abs(p.P0) / (p.lam + 2 * p.G) * (p.B - p.A), th * 1e-2 * p.B, 1e-30)
    atol = np.concatenate((np.full(Ns, us), np.full(Ns, max(th, us)), np.full(Ns, us * 10),
                           np.full(Ns, max(th, us) * 10), np.full(Nf - 1, max(th, us) * 10),
                           [1e-30, 1e-30, 1e-30])) * 1e-2 * rtol
    # budget states: absolute tolerance relative to the work scale (avoids a
    # degenerate numerical-Jacobian column scaling)
    Escale = max(abs(p.P0) * us, (p.G * th * th) * (p.B**2 - p.A**2), 1e-30)
    atol[-3:] = Escale * rtol
    sol = solve_ivp(f, (0.0, t_end), np.zeros(n + 3), method=method, rtol=rtol, atol=atol, t_eval=t_eval)
    if not sol.success:
        raise RuntimeError(sol.message)
    sol.wall = _time.time() - t0
    sol.y_full = sol.y
    sol.y = sol.y[:n]
    H = history(m, sol)
    E = H["KE_s"] + H["V_s"] + H["KE_f"]
    Hm = H["H_s"] + H["H_f"]
    H["E_tot"] = E
    H["H_tot"] = Hm
    H["res_E"] = E - E[0] - sol.y_full[n]
    H["res_H"] = Hm - Hm[0] - sol.y_full[n + 1]
    H["diss_int"] = sol.y_full[n + 2]
    return m, sol, H
