#!/usr/bin/env python3
"""
Manufactured solution for 2-D transient ALE FSI: a channel with a compliant,
finite-strain top wall.

Fluid  : 0 < x < L, 0 < y < h(x,t)            (Eulerian, moving domain)
Solid  : 0 < X < L, H < Y < H + T              (reference configuration)
Interface (reference): Y = H

Construction (see the tutorial README.md and the reference given there):
  1. Choose the solid displacement u(X,Y,t) in closed form.  Its trace on Y = H
     fixes the interface motion x_G(X,t) = X + xi(X,t), y_G(X,t) = H + eta(X,t).
     The tangential part xi is a Moebius map of [0,L] onto itself, so the
     interface parametrisation can be inverted in closed form: X = Xi(x,t).
  2. The fluid velocity comes from a stream function psi(x,y,t) = Psi(x,s,t)
     with s = y/h(x,t), so div v = 0 exactly.  Hermite polynomials in s impose
     no-slip at y = 0 and v = du/dt on the deformed interface (s = 1).
  3. A term s^2 (1-s)^2 Phi(x,t) changes only the second normal derivative of
     psi at s = 1, i.e. the fluid shear traction there.  Phi is solved
     (linearly) so that the fluid shear traction equals the solid shear
     traction.  The pressure on the interface is set so that the normal
     traction also matches.  Hence sigma_f n = sigma_s n holds exactly and no
     interface source is needed.
  4. Fluid and solid momentum sources follow from the residuals.

Run:  ./run_checks.sh (see the tutorial README.md)
"""
import argparse
import os
import random
import sys
import time as _time

import sympy as sp

parser = argparse.ArgumentParser()
parser.add_argument("--law", default="neoHookean", choices=["neoHookean", "StVK"])
parser.add_argument("--quick", action="store_true", help="skip slow symbolic checks")
parser.add_argument("--no-codegen", action="store_true")
parser.add_argument("--rho-f", default="1", help="fluid density (rational string)")
parser.add_argument("--rho-s", default="1", help="solid reference density (rational string)")
parser.add_argument("--outdir", default="generated")
parser.add_argument("--mu-f", default="1/100", help="fluid dynamic viscosity (rational string)")
parser.add_argument("--variant", default="full", choices=["full", "flat"],
                    help="full: tilting, travelling, stretching interface; "
                         "flat: the interface stays flat and moves only vertically (orthogonal ALE mesh)")
args = parser.parse_args()

R = sp.Rational
_T0 = _time.time()
def tick(msg):
    print(f'[{_time.time()-_T0:7.1f}s] {msg}', flush=True)
pi = sp.pi

# ---------------------------------------------------------------------------
# Parameters (exact rationals so that the checks are exact)
# ---------------------------------------------------------------------------
P = dict(
    L=R(2), H=R(1), T=R(1, 2),        # channel length, height, wall thickness
    omega=2 * pi,                     # angular frequency (period 1)
    a_eta=R(1, 10),                   # interface normal amplitude (0.1 H, x1.5 with the wave)
    b_eta=R(1, 2),                    # travelling-wave part of eta
    c_kap=R(1, 10),                   # tangential (Moebius) stretch amplitude
    gam=R(1, 4),                      # through-thickness shape of solid field
    e_sh=R(1, 50),                    # small solid interface shear (0: zero shear traction)
    Q0=R(1),                          # mean through-flow (per unit depth)
    A_Om=R(8),                        # vortical stream-function amplitude
    pi0=R(1, 2),                      # free pressure amplitude
    rho_f=sp.Rational(args.rho_f), mu_f=sp.Rational(args.mu_f),  # fluid density, dynamic viscosity
    rho_s=sp.Rational(args.rho_s),    # solid reference density (ratio 1: added mass)
    E_s=R(10), nu_s=R(3, 10),         # solid Young's modulus, Poisson ratio
)

L, H, T = P["L"], P["H"], P["T"]
om = P["omega"]
rho_f, mu_f, rho_s = P["rho_f"], P["mu_f"], P["rho_s"]
E_s, nu_s = P["E_s"], P["nu_s"]
mu_s = E_s / (2 * (1 + nu_s))
lam_s = nu_s * E_s / ((1 + nu_s) * (1 - 2 * nu_s))
K_s = lam_s + R(2, 3) * mu_s          # 3-D bulk modulus

x, y, t = sp.symbols("x y t", real=True)
X, Y = sp.symbols("X Y", real=True)
s_ = sp.symbols("s", real=True)

# Smooth start: r(0) = r'(0) = 0, so u = du/dt = v = 0 and p = 0 at t = 0.
r1 = (1 - sp.cos(om * t)) / 2
r2 = sp.sin(om * t) ** 2                     # tangential stretch: twice the frequency

# ---------------------------------------------------------------------------
# 1. Solid displacement and interface motion
# ---------------------------------------------------------------------------
Sx = X / L
kap = 1 + (P["c_kap"] if args.variant == "full" else 0) * r2   # kappa(t) > 0
# Moebius map of [0,1] onto itself with fixed ends: tangential interface motion
xG_over_L = Sx / (Sx + kap * (1 - Sx))
xG = L * xG_over_L                                      # deformed abscissa of interface point X
xi = xG - X                                             # tangential interface displacement

# Closed-form inverse of the tangential map: X = Xi(x,t)
xL = x / L
Xi = L * kap * xL / (1 - xL + kap * xL)

def at_x(expr):
    """Interface quantity (function of X,t) expressed at Eulerian abscissa x."""
    return expr.subs(X, Xi)

# Normal interface motion is specified in EULERIAN form, as the channel height
# h(x,t) = H + etahat(x,t).  The solid displacement then follows by composition,
# eta(X,t) = etahat(xG(X,t),t).  This keeps the swept area int h dx elementary.
# sin^2 gives zero wall slope at the clamped ends (needed for the shear design below)
if args.variant == "full":
    etahat = P["a_eta"] * sp.sin(pi * x / L) ** 2 * r1 * (1 + P["b_eta"] * sp.sin(2 * pi * x / L - om * t))
else:
    # flat variant: uniform vertical motion, no tangential stretch; the solid
    # sides are then not clamped (exact Dirichlet data instead)
    etahat = P["a_eta"] * r1 * (1 + P["b_eta"] * sp.sin(om * t))
eta = etahat.subs(x, xG)
h = H + etahat

yh = (Y - H) / T                                        # 0 at interface, 1 at the clamped top
g = (1 - yh) * (1 + P["gam"] * yh)                      # g(0)=1, g(1)=0
# Interface shear design.  For the solids4foam neo-Hookean law the tangential
# traction on the interface is mu J^(-2/3) tau.(F N)/(...) because the
# volumetric and I1 terms of P are parallel to F^-T N, i.e. to the normal.
# It therefore vanishes iff the deformed material fibre dx/dY stays orthogonal
# to the deformed interface tangent dx/dX.  m0 enforces this; e_sh then adds a
# small, smooth, exactly known shear so that the fluid wall shear stays O(mu U/H).
u_y = eta * g
gp0 = sp.diff(g, Y).subs(Y, H)                          # g'(Y=H)
m0 = -sp.diff(eta, X) * (1 + eta * gp0) / (1 + sp.diff(xi, X)) - xi * gp0
esh = sp.Symbol("e_sh")
m_extra = esh * sp.sin(2 * pi * Sx) * r1
u_x = xi * g + (Y - H) * (1 - yh) ** 2 * (m0 + m_extra)
U = sp.Matrix([u_x, u_y])

# Interface traces (functions of X, t)
yG = H + eta
VG = sp.Matrix([sp.diff(xi, t), sp.diff(eta, t)])       # interface velocity
Vx = at_x(VG[0])
Vy = at_x(VG[1])

# ---------------------------------------------------------------------------
tick('kinematics built')
# 2. Solid stress (plane strain, 3-D law with F_zz = 1)
# ---------------------------------------------------------------------------
F2 = sp.eye(2) + U.jacobian([X, Y])
F = sp.Matrix([[F2[0, 0], F2[0, 1], 0], [F2[1, 0], F2[1, 1], 0], [0, 0, 1]])
J = F2[0, 0] * F2[1, 1] - F2[0, 1] * F2[1, 0]

def cauchy_neoHookean(F, J):
    # As neoHookeanElasticMechanicalConstitutiveLaw::evaluate:
    # (mu/J)*dev(J^(-2/3)*b) + 0.5*K*(J^2 - 1)/J*I, with the 3-D K
    b = F * F.T
    bbar = J ** R(-2, 3) * b
    devbbar = bbar - bbar.trace() / 3 * sp.eye(3)
    return mu_s / J * devbbar + K_s / 2 * (J ** 2 - 1) / J * sp.eye(3)

def cauchy_StVK(F, J):
    Ee = (F.T * F - sp.eye(3)) / 2
    S = 2 * mu_s * Ee + lam_s * Ee.trace() * sp.eye(3)
    return F * S * F.T / J

sigma_s = cauchy_neoHookean(F, J) if args.law == "neoHookean" else cauchy_StVK(F, J)
# explicit inverse (avoid sympy's simplifying Matrix.inv on large expressions)
Finv = sp.Matrix([[F2[1, 1] / J, -F2[0, 1] / J, 0], [-F2[1, 0] / J, F2[0, 0] / J, 0], [0, 0, 1]])
Pk1 = J * sigma_s * Finv.T                              # first Piola-Kirchhoff

# Outward normal of the fluid on the interface (points into the solid, ~ +y),
# from Nanson: n da = J F^{-T} N dA with the solid reference normal N = (0,-1,0).
Finv_T = Finv.T
nS = Finv_T * sp.Matrix([0, -1, 0])                     # solid outward (into the fluid)
nF_ref = -nS / sp.sqrt((nS.T * nS)[0])                  # fluid outward unit normal
# sympy writes sqrt(a**2) as Abs(a) when a component vanishes identically (flat
# variant); the arguments are positive (F11 > 0, J > 0), and symengine cannot
# differentiate Abs, so drop it. The C++ check of the normal against the FD
# tangent of the interface would catch a sign error.
nF_ref = nF_ref.replace(sp.Abs, lambda a: a)
tS_on_nF = sigma_s * nF_ref                             # solid Cauchy stress . n_f
# Evaluate on the interface Y = H as functions of (X,t), then of (x,t)
tS_G = tS_on_nF.subs(Y, H)
tS_x = at_x(tS_G[0])
tS_y = at_x(tS_G[1])

# ---------------------------------------------------------------------------
tick('solid stress built')
# 3. Fluid stream function
# ---------------------------------------------------------------------------
# Q(x,t) = psi on the interface, Q_x = -h_t (volume change).  The OUTLET flux
# Q(L,t) = Q_out(t) = Q0 r(t) >= 0 is prescribed and the inlet absorbs the
# volume change of the channel, so that the outlet (a Neumann velocity
# boundary in the solver) never sees backflow, which makes the discrete
# problem unstable there (classical open-boundary backflow instability).
xd = sp.Dummy('xd', real=True)
A_xL = sp.integrate(etahat.subs(x, xd), (xd, x, L))     # swept area from x to the outlet
Q_out = P["Q0"] * r1
Q = Q_out + sp.diff(A_xL, t)

tick('area/Q built')
Phi = sp.Function("Phi")(x, t)
B = 3 * s_ ** 2 - 2 * s_ ** 3          # B(0)=B'(0)=0, B(1)=1, B'(1)=0
C = s_ ** 3 - s_ ** 2                  # C(0)=C'(0)=C(1)=0, C'(1)=1
W2 = s_ ** 2 * (1 - s_) ** 2           # zero value/slope at both ends, W2''(1)=2
W3 = s_ ** 2 * (1 - s_) ** 3           # zero value/slope/curvature at s=1
# travelling vortex train, faded out at the inlet and outlet by sin^2(pi x/L)
# so that it does not reverse the flow at the outlet
Om = P["A_Om"] * r1 * sp.sin(pi * x / L) ** 2 * sp.sin(2 * pi * x / L - om * t) * (1 + R(1, 2) * sp.cos(pi * s_))
Psi_s = Q * B + h * Vx * C + W2 * Phi + W3 * Om
psi_Phi = Psi_s.subs(s_, y / h)

def vel(psi):
    return sp.Matrix([sp.diff(psi, y), -sp.diff(psi, x)])

def grad2(v):
    return v.jacobian([x, y])

# Interface geometry in Eulerian form: unit normal of the curve y = h(x,t)
hx = sp.diff(h, x)
nrm = sp.sqrt(1 + hx ** 2)
nG = sp.Matrix([-hx, 1]) / nrm              # fluid outward unit normal
tG = sp.Matrix([1, hx]) / nrm               # unit tangent

# Tangential traction balance -> solve for Phi on the interface (linear)
v_Phi = vel(psi_Phi)
Gv_Phi = grad2(v_Phi)
D_Phi = (Gv_Phi + Gv_Phi.T) / 2
D_Phi_G = D_Phi.subs(y, h)
tau_f = (tG.T * (2 * mu_f * D_Phi_G) * nG)[0]
tau_s = tG[0] * tS_x + tG[1] * tS_y
phi_sym = sp.Symbol("phi_sym")
# Derivatives of Phi multiply W2(1) = W2'(1) = 0 and must be absent BEFORE Phi
# is replaced by a symbol (checked on the unevaluated expression)
assert not any(d.has(Phi) for d in tau_f.atoms(sp.Derivative)), \
    "tangential traction depends on derivatives of Phi"
tau_f_lin = tau_f.subs(Phi, phi_sym).doit()
c1 = sp.diff(tau_f_lin, phi_sym)
c0 = tau_f_lin.subs(phi_sym, 0)
assert sp.diff(c1, phi_sym) == 0
Phi_trac = (tau_s - c0) / c1
c1_expected = 2 * mu_f * (1 + hx ** 2) / h ** 2       # derived by hand (see README)

tick('Phi_trac solved')
# ---------------------------------------------------------------------------
# From here on: symengine (fast differentiation of large shared expressions)
# ---------------------------------------------------------------------------
import symengine as se
xs, ys, ts, Xs, Ys = [se.Symbol(n) for n in ("x", "y", "t", "X", "Y")]
def S(e, esh_val=None):
    e = e.subs(esh, P["e_sh"] if esh_val is None else esh_val)
    return se.sympify(e).subs({se.sympify(x): xs, se.sympify(y): ys, se.sympify(t): ts,
                               se.sympify(X): Xs, se.sympify(Y): Ys})
def SM(M):
    return se.Matrix([[S(M[i, j]) for j in range(M.shape[1])] for i in range(M.shape[0])])

psi = S(psi_Phi.subs(Phi, Phi_trac))
h_e, Pi_e = S(h), S(P["pi0"] * r1 * sp.cos(pi * x / L))
nG_e = [S(nG[0]), S(nG[1])]
tS_e = [S(tS_x), S(tS_y)]
v = [se.diff(psi, ys), -se.diff(psi, xs)]
Gv = [[se.diff(v[i], xs), se.diff(v[i], ys)] for i in range(2)]
Dm = [[(Gv[i][j] + Gv[j][i]) / 2 for j in range(2)] for i in range(2)]
tick('velocity built')

# Pressure: p = p_G(x,t) + (1 - s) Pi(x,t), with p_G from the normal balance
DG = [[Dm[i][j].subs({ys: h_e}) for j in range(2)] for i in range(2)]
p_G = 2 * se.sympify(mu_f) * sum(nG_e[i] * DG[i][j] * nG_e[j] for i in range(2) for j in range(2)) \
    - (nG_e[0] * tS_e[0] + nG_e[1] * tS_e[1])
p = p_G + (1 - ys / h_e) * Pi_e
sigma_f = [[-p * (1 if i == j else 0) + 2 * se.sympify(mu_f) * Dm[i][j] for j in range(2)] for i in range(2)]
tick('pressure built')

# ---------------------------------------------------------------------------
# 4. Sources
# ---------------------------------------------------------------------------
# Fluid (per unit volume): rho (dv/dt + (v.grad) v) - div sigma_f = f
rf = se.sympify(rho_f)
f_fluid = [rf * (se.diff(v[i], ts) + Gv[i][0] * v[0] + Gv[i][1] * v[1])
           - se.diff(sigma_f[i][0], xs) - se.diff(sigma_f[i][1], ys) for i in range(2)]
tick('fluid source built')

# Solid (per unit reference volume): rho0 d2u/dt2 - Div P = b0
Ue = [S(U[0]), S(U[1])]
Ute = [se.diff(Ue[i], ts) for i in range(2)]
Pk = [[S(Pk1[i, j]) for j in range(2)] for i in range(2)]
b0_solid = [se.sympify(rho_s) * se.diff(Ute[i], ts) - se.diff(Pk[i][0], Xs) - se.diff(Pk[i][1], Ys)
            for i in range(2)]
sig_s = [[S(sigma_s[i, j]) for j in range(3)] for i in range(3)]
tick('solid source built')

# ---------------------------------------------------------------------------
# 5. Verification (symengine with MPFR, 133 bits ~ 40 digits)
# ---------------------------------------------------------------------------
PREC = 133
results = []

def report(name, ok, detail=""):
    results.append((name, ok, detail))
    print(f"[{'PASS' if ok else 'FAIL'}] {name} {detail}", flush=True)

rng = random.Random(1234)
def rnd(a, b):
    return se.sympify(a) + (se.sympify(b) - se.sympify(a)) * se.Rational(rng.randint(1, 10**6), 10**6)

def ev(e, d):
    # substitute floats (exact rational substitution can leave symbolic
    # constants that symengine's MPFR evalf does not implement)
    dd = {k: se.RealMPFR(str(se.sympify(v).evalf(PREC, real=True)), PREC) for k, v in d.items()}
    return se.sympify(e).subs(dd).evalf(PREC, real=True)

# (a) inverse map: Xi(xG(X,t),t) == X exactly (rational identity)
inv_err = sp.cancel(sp.together(Xi.subs(x, xG) - X))
report("closed-form inverse Xi(xG(X,t),t) = X (symbolic)", inv_err == 0, str(inv_err))

# (b) continuity: v = (psi_y, -psi_x) -> div v = 0 identically (symbolic, generic psi)
fgen = sp.Function("f")(x, y)
report("div curl(psi) = 0 (symbolic, generic psi)",
       sp.diff(sp.diff(fgen, y), x) - sp.diff(sp.diff(fgen, x), y) == 0)

xG_e, yG_e = S(xG), S(yG)
nS_e = [S(nF_ref[i].subs(Y, H)) for i in range(2)]
tSs_e = [S(tS_on_nF[i].subs(Y, H)) for i in range(2)]
Ut_G = [Ute[i].subs({Ys: S(H)}) for i in range(2)]
divv = se.diff(v[0], xs) + se.diff(v[1], ys)

maxdiv = maxkin = maxwall = maxtrac = maxh = se.Integer(0)
NPT = 8
for k in range(NPT):
    Xv, tv = rnd(0, L), rnd(0, 1)
    dX = {Xs: Xv, ts: tv}
    xg = ev(xG_e, dX); yg = ev(yG_e, dX)
    dx = {xs: xg, ys: yg, ts: tv}
    vf = [ev(v[i], dx) for i in range(2)]
    vs = [ev(Ut_G[i], dX) for i in range(2)]
    maxkin = max(maxkin, abs(vf[0] - vs[0]), abs(vf[1] - vs[1]))
    maxh = max(maxh, abs(ev(h_e, {xs: xg, ts: tv}) - yg))
    vw = [ev(v[i], {xs: xg, ys: 0, ts: tv}) for i in range(2)]
    maxwall = max(maxwall, abs(vw[0]), abs(vw[1]))
    nf = [ev(nS_e[i], dX) for i in range(2)]
    tsol = [ev(tSs_e[i], dX) for i in range(2)]
    sf = [[ev(sigma_f[i][j], dx) for j in range(2)] for i in range(2)]
    tf = [sf[i][0] * nf[0] + sf[i][1] * nf[1] for i in range(2)]
    maxtrac = max(maxtrac, abs(tf[0] - tsol[0]), abs(tf[1] - tsol[1]))
    yv = rnd(0, 1) * yg
    maxdiv = max(maxdiv, abs(ev(divv, {xs: xg, ys: yv, ts: tv})))
    print(f"  point {k}: kin {float(maxkin):.1e} trac {float(maxtrac):.1e} |t_s| {float(abs(tsol[0])+abs(tsol[1])):.3e}", flush=True)
report("fluid domain top = deformed interface, h(xG) = yG", maxh < 1e-30, f"max {float(maxh):.1e}")
report(f"div v = 0 at 40 digits ({NPT} random points)", maxdiv < 1e-28, f"max {float(maxdiv):.1e}")
report("kinematic v_f(xG,yG) = du/dt(X,H) (40 digits)", maxkin < 1e-28, f"max {float(maxkin):.1e}")
report("no-slip at y = 0 (40 digits)", maxwall < 1e-28, f"max {float(maxwall):.1e}")
report("traction sigma_f n = sigma_s n (40 digits)", maxtrac < 1e-25, f"max {float(maxtrac):.1e}")

# (b2) shear design: with e_sh = 0 the solid tangential traction on the
#      interface vanishes identically (neo-Hookean); with e_sh > 0 it is small.
tG_s = [S(tS_G[i], 0) for i in range(2)]
tvec = [se.diff(xG_e, Xs), se.diff(yG_e, Xs)]           # deformed tangent dx/dX
maxsh0 = maxsh = maxtn = se.Integer(0)
for _ in range(6):
    dX = {Xs: rnd(0, L), Ys: S(H), ts: rnd(0, 1)}
    tv_ = [ev(tvec[i], dX) for i in range(2)]
    nt = se.sqrt(tv_[0] ** 2 + tv_[1] ** 2)
    maxsh0 = max(maxsh0, abs(sum(ev(tG_s[i], dX) * tv_[i] for i in range(2)) / nt))
    maxsh = max(maxsh, abs(sum(ev(tSs_e[i], dX) * tv_[i] for i in range(2)) / nt))
    maxtn = max(maxtn, abs(sum(ev(tSs_e[i], dX) * ev(nS_e[i], dX) for i in range(2))))
if args.law == "neoHookean":
    report("solid interface shear traction = 0 when e_sh = 0 (40 digits)", maxsh0 < 1e-28, f"max {float(maxsh0):.1e}")
print(f"  interface |shear traction| <= {float(maxsh):.3e}, |normal traction| <= {float(maxtn):.3e} (sampled)")

# (b3) coefficient of Phi in the fluid shear traction: c1 = 2 mu (1 + h_x^2)/h^2
c1_err = se.Integer(0)
c1_e, c1x_e = S(c1), S(c1_expected)
for _ in range(5):
    dd = {xs: rnd(0, L), ts: rnd(0, 1)}
    c1_err = max(c1_err, abs(ev(c1_e, dd) - ev(c1x_e, dd)))
report("Phi coefficient c1 = 2 mu (1 + h_x^2)/h^2 (40 digits)", c1_err < 1e-30, f"max {float(c1_err):.1e}")

# (c) solid clamped boundaries: u = 0 at X = 0, X = L, Y = H + T
maxclamp = se.Integer(0)
for _ in range(5):
    tv = rnd(0, 1); Xv = rnd(0, L); Yv = rnd(H, H + T)
    for pt in ([(0, Yv), (L, Yv), (Xv, H + T)] if args.variant == "full" else [(Xv, H + T)]):
        for i in range(2):
            maxclamp = max(maxclamp, abs(ev(Ue[i], {Xs: pt[0], Ys: pt[1], ts: tv})))
report("solid clamped on X=0, X=L, Y=H+T" if args.variant == "full" else "solid clamped on Y=H+T", maxclamp < 1e-30, f"max {float(maxclamp):.1e}")

# (d) zero initial state
mi = se.Integer(0)
d0 = {Xs: rnd(0, L), Ys: rnd(H, H + T), xs: rnd(0, L), ys: rnd(0, 1), ts: 0}
for e in [*Ue, *Ute, *v, p]:
    mi = max(mi, abs(ev(e, d0)))
report("zero initial state u = u_t = v = p = 0 at t = 0", mi < 1e-30, f"max {float(mi):.1e}")

# (e) Ue, sigma_s consistency: the solid law in the derivation is the
#     solids4foam neoHookeanElastic / StVenantKirchhoffElastic formula
print("checks done", flush=True)

# (f) reference values at fixed points (40 digits) for the direct comparison
#     with the generated C++ header (check_generated.C reads this file)
refpts = [(R(37, 100), R(41, 100), R(23, 100)), (R(151, 100), R(83, 100), R(71, 100)),
          (R(1, 1), R(1, 2), R(1, 2)), (R(9, 50), R(97, 100), R(9, 10))]
os.makedirs(os.path.join(os.path.dirname(os.path.abspath(__file__)), args.outdir), exist_ok=True)
with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), args.outdir, "reference_values.txt"), "w") as fh:
    for (a_, b_, c_) in refpts:
        dF = {xs: a_ * 2, ys: b_, ts: c_}
        dS = {Xs: a_ * 2, Ys: S(H) + b_ * S(T), ts: c_}
        vals = [ev(e, dF) for e in (*v, p / se.sympify(rho_f), *[fi / se.sympify(rho_f) for fi in f_fluid])]
        vals += [ev(e, dS) for e in (*Ue, *Ute, *b0_solid)]
        fh.write(" ".join(str(q) for q in (float(a_ * 2), float(b_), float(c_), float(S(H) + b_ * S(T)))))
        fh.write(" " + " ".join(str(se.Float(q, PREC)) for q in vals) + "\n")
tick("reference values written")

# ---------------------------------------------------------------------------
# 6. Code generation: C++ header with plain double functions
# ---------------------------------------------------------------------------
if not args.no_codegen:
    import codegen
    outdir = os.path.join(os.path.dirname(os.path.abspath(__file__)), args.outdir)
    os.makedirs(outdir, exist_ok=True)
    rho_f_e = se.sympify(rho_f)
    funcs = [
        # name, args, outputs (list of (cname, expr))
        ("fluid", ["x", "y", "t"], [("U", v), ("p", [p / rho_f_e]), ("f", [fi / rho_f_e for fi in f_fluid]),
                                     ("gradU", [Gv[0][0], Gv[0][1], Gv[1][0], Gv[1][1]]),
                                     ("gradp", [se.diff(p, xs) / rho_f_e, se.diff(p, ys) / rho_f_e])]),
        ("channelHeight", ["x", "t"], [("h", [h_e])]),
        ("solid", ["X", "Y", "t"], [("D", Ue), ("U", Ute), ("b0", b0_solid),
                                     ("sigma", [sig_s[0][0], sig_s[0][1], sig_s[0][2], sig_s[1][1], sig_s[1][2], sig_s[2][2]])]),
        ("interface", ["X", "t"], [("x", [xG_e, yG_e]), ("traction", tSs_e), ("n", nS_e)]),
    ]
    codegen.write_header(os.path.join(outdir, "mmsFsiExact.H"), funcs, P, args.law, tick, args.variant)
    tick("code generated")

nfail = sum(1 for r in results if not r[1])
print(f"{len(results)-nfail}/{len(results)} checks passed")
sys.exit(1 if nfail else 0)
