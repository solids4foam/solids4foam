#!/usr/bin/env python3
"""Usage: python3 stokes_rigid_ring_check.py [E in Pa, default 2e7]   (about 1 minute at 2e7, 6 minutes at 2e8; no CFD)
No-CFD limiting-case check of the nonlinear semi-analytical reference: RIGID RING, small amplitude, swirl only = Stokes layer in the annulus,
compared in the TIME DOMAIN over the same 0-2 s as the benchmark (exact Laplace-domain solution, numerically inverted; no long-time steady state needed).
Assumptions that must hold (arranged and checked below):
  * very stiff solid (neo-Hookean, small strain = linear; E given on the command line), so the ring rotates with the inner rotation;
  * theta0 = 1e-3 rad, P0 = 0, fluid_convection = False: radial displacement O(theta0^2), the swirl problem is linear;
  * fluid: w(r,t) on B <= r <= c, no slip at B, sigma_r-theta = mu (w' - w/r) = 0 at r = c (the benchmark's open outer boundary).
Closed form: w^(r,s) = A I1(k r) + B K1(k r), k = sqrt(s/nu); w^(B,s) = B*omega^_B(s); sigma^_r-theta(B,s) = G(s) omega^_B(s).
The wall angular velocity omega_B(t) is the one computed by the reference itself (its interface rotation rate), whose agreement with the prescribed rigid
rotation is reported; the Laplace transform of that smooth function is obtained from the analytic drive (step response of the transfer function, convolved with d omega/dt of the analytic drive).
Compared: sigma_r-theta(B,t) (reference 'sig_rt_b'), T_B = B^2 sigma_r-theta ('T_B'), and H_f = rho int w r^2 dr ('H_f')."""
import sys, numpy as np
sys.dont_write_bytecode = True
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import rt_reference as rt
from scipy.special import iv, kv
from scipy.integrate import quad
E = float(sys.argv[1]) if len(sys.argv) > 1 else 2.0e7
th0, Om, Tr = 1e-3, 3.0, 1.0
p = rt.Params(E=E, law="neoHookean", theta0=th0, Omega=Om, T_ramp=Tr, ramp="Cinf", P0=0.0, fluid_convection=False)
te = np.linspace(0.05, 2.0, 40)
m, sol = rt.solve(p, 2.0, Ns=16, Nf=32, rtol=1e-10, t_eval=te, max_wall=1500)
h = rt.history(m, sol); t = h["t"]
nu, B, c, mu, rho = p.nu_f, p.B, p.c, p.mu_f, p.rho_f
# drive: alpha_A = theta0 S(t) sin(Om t); omega_A = d/dt
wA = lambda tt: rt.rotation_drive(p, tt)[1]
from numpy.polynomial.legendre import leggauss
_gx, _gw = leggauss(48)
aA = lambda tt: rt.rotation_drive(p, tt)[2]          # d omega_A / dt
def trans(s, r):
    """overflow-safe (exponentially scaled Bessel functions): I(z) = ive(z) e^{Re z}, K(z) = kve(z) e^{-z}  (complex phase kept): W(x) = At ive1(kx) e^{a(x-c)} + Bt kve1(kx) e^{-k(x-B)}, a = Re k"""
    from scipy.special import ive, kve
    k = np.sqrt(s / nu); a = k.real
    eAB, eBc = np.exp(-a * (c - B)), np.exp(-a * (c - B))
    Mx = np.array([[ive(1, k * B) * np.exp(a * (B - c)), kve(1, k * B)],
                   [k * ive(0, k * c) - 2 * ive(1, k * c) / c, (-k * kve(0, k * c) - 2 * kve(1, k * c) / c) * np.exp(-k * (c - B))]], dtype=complex)
    At, Bt = np.linalg.solve(Mx, np.array([B, 0], dtype=complex))
    W = lambda x: At * ive(1, k * x) * np.exp(a * (x - c)) + Bt * kve(1, k * x) * np.exp(-k * (x - B))
    dW = lambda x: At * k * (ive(0, k * x) - ive(1, k * x) / (k * x)) * np.exp(a * (x - c)) + Bt * k * (-kve(0, k * x) - kve(1, k * x) / (k * x)) * np.exp(-k * (x - B))
    return mu * (dW(B) - W(B) / B), W
def talbot(F, tau, M=32):
    """Abate-Valko fixed Talbot inverse Laplace transform of F at time tau > 0"""
    r = 2.0 * M / (5.0 * tau)
    acc = 0.5 * F(r + 0j) * np.exp(r * tau)
    for k in range(1, M):
        th = k * np.pi / M; sg = r * th * (1 / np.tan(th) + 1j); tau_k = th + (th / np.tan(th) - 1) / np.tan(th)
        acc += (np.exp(tau * sg) * F(sg) * (1 + 1j * tau_k)).real
    return r / M * acc
gx2, gw2 = leggauss(60); xx = 0.5 * (c - B) * gx2 + 0.5 * (c + B); ww = 0.5 * (c - B) * gw2
GH = lambda s: rho * np.sum(ww * np.array([trans(s, B)[1](x) * x * x for x in xx]))
def response(Gfun, tt):
    """quantity(t) = int_0^t g(t-x) omega'(x) dx with g = L^-1[G(s)/s] (step response; omega(0)=0); substitution x = t - u^2 removes the integrable u^-1 singularity of g(u^2)"""
    umax = np.sqrt(tt); u = 0.5 * umax * (_gx + 1); w = 0.5 * umax * _gw
    return float(np.sum(w * np.array([2.0 * ui * talbot(lambda s: Gfun(s) / s, ui * ui) * aA(tt - ui * ui) for ui in u])))
sig = np.array([response(lambda s: trans(s, B)[0], tt) for tt in t])
Hf = np.array([response(GH, tt) for tt in t])
TB = B**2 * sig
print("E = %.1e Pa: checks: |b - B| max %.2e m; wall omega vs prescribed rigid rotation: max|omega_B - omega_A|/max|omega_A| = %.2e" % (E, np.max(np.abs(h["b"] - B)), np.max(np.abs(h["omega_B"] - np.array([wA(x) for x in t]))) / np.max(np.abs(np.array([wA(x) for x in t])))))
print("convention: median reference T_B / (B^2 sig_rt_b) = %.9f" % np.median(h["T_B"][np.abs(h["sig_rt_b"]) > 1e-3 * np.max(np.abs(h["sig_rt_b"]))] / (B**2 * h["sig_rt_b"][np.abs(h["sig_rt_b"]) > 1e-3 * np.max(np.abs(h["sig_rt_b"]))])))
for name, ref, cf in (("sigma_r-theta(B)", h["sig_rt_b"], sig), ("T_B", h["T_B"], TB), ("H_f", h["H_f"], Hf)):
    print("%-16s peak(closed form) %.6e   max|reference - closed form| / peak = %.3e" % (name, np.max(np.abs(cf)), np.max(np.abs(ref - cf)) / np.max(np.abs(cf))))
print("solve wall %.1f s, %d output times in 0.05-2 s" % (sol.wall, t.size))
