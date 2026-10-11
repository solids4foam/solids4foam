#!/usr/bin/env python3
"""Exact linear n-mode frequencies of an elastic ring in a fluid annulus.

The ring is a plane-strain linear elastic annulus Ri < r < Ro, traction free
on its inner surface, surrounded by incompressible fluid in Ro < r < b with a
rigid outer wall at r = b. Small-amplitude motion is sought in the form

    u_r = U(r) cos(n theta) exp(s t),  u_theta = V(r) sin(n theta) exp(s t)

and the radial eigenproblem of the solid is solved by Chebyshev collocation,
so the ring is an exact two-dimensional continuum (no thin-ring
approximation). The fluid enters through the interface tractions:

- Inviscid potential flow. The closed-form annulus solution gives the
  interface pressure p = m_a(Ro) s^2 U(Ro) cos(n theta), with the added mass
  per unit area

      m_a(R) = rho_f R (b^2n + R^2n) / (n (b^2n - R^2n)),

  and the tangential ring motion does not couple.
- Linearised incompressible Navier-Stokes flow with no slip on both walls.
  The stream function Psi(r) sin(n theta) is
      Psi = A r^n + B r^-n + C I_n(k r) + D K_n(k r),  k = sqrt(s / nu_f),
  whose modified Bessel functions are evaluated with their large-argument
  expansions (the Stokes layer is thin). The interface tractions then depend
  on s, and the eigenvalue is found by iterating on s. It is complex, and
  gives the viscous frequency shift and damping.

The thin-ring (inextensional) reference, with the fluid acting on the
interface radius Ro, is also given:

    omega_dry^2 = D n^2 (n^2 - 1)^2 / (rho_s h a^4 (n^2 + 1)),
    omega_wet / omega_dry = [1 + m_a(Ro) Ro n^2 / (rho_s h a (n^2 + 1))]^-1/2,

with D = E h^3 / (12 (1 - nu^2)), mean radius a and thickness h.

This script needs numpy. The verification driver does not: it reads the
values this script prints with --json, which are stored in
reference/ringAddedMass_verification_references.json.

usage: ring_exact_frequencies.py [--json] [--check]
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np

SCRIPT = Path(__file__).resolve()
REFERENCE_FILE = (
    SCRIPT.parents[1] / "reference" / "ringAddedMass_verification_references.json"
)


def chebyshev(n_intervals: int, r_min: float, r_max: float):
    """Chebyshev-Gauss-Lobatto points on [r_min, r_max] and the derivative
    matrix; the first point is r_max and the last r_min."""
    x = np.cos(np.pi * np.arange(n_intervals + 1) / n_intervals)
    c = np.hstack([2.0, np.ones(n_intervals - 1), 2.0]) \
        * (-1.0) ** np.arange(n_intervals + 1)
    dx = np.subtract.outer(x, x)
    d = np.outer(c, 1.0 / c) / (dx + np.eye(n_intervals + 1))
    d -= np.diag(d.sum(axis=1))
    r = r_min + 0.5 * (x + 1.0) * (r_max - r_min)
    return r, d * 2.0 / (r_max - r_min)


def added_mass(n: int, radius: float, outer: float, rho_f: float) -> float:
    """Potential-flow added mass per unit area of mode n on a cylinder of
    radius `radius` inside a rigid cylinder of radius `outer`."""
    return rho_f * radius / n * (outer ** (2 * n) + radius ** (2 * n)) \
        / (outer ** (2 * n) - radius ** (2 * n))


def solid_matrices(n, ri, ro, young, poisson, rho_s, n_solid):
    """Collocation matrices A, B of s x = ... for x = (U, V, W_U, W_V), with
    traction-free rows at both surfaces; returns the row and column of the
    interface traction equations and velocities."""
    lam = young * poisson / ((1 + poisson) * (1 - 2 * poisson))
    mu = young / (2 * (1 + poisson))
    r, d = chebyshev(n_solid, ri, ro)
    ns = n_solid + 1
    eye = np.eye(ns)
    rinv = np.diag(1.0 / r)
    # Dilatation U' + U/r + nV/r and rotation V' + V/r + nU/r
    del_u, del_v = d + rinv, n * rinv
    rot_u, rot_v = n * rinv, d + rinv
    # Navier equations
    lru = (lam + 2 * mu) * d @ del_u - mu * n * rinv @ rot_u
    lrv = (lam + 2 * mu) * d @ del_v - mu * n * rinv @ rot_v
    ltu = -(lam + 2 * mu) * n * rinv @ del_u + mu * d @ rot_u
    ltv = -(lam + 2 * mu) * n * rinv @ del_v + mu * d @ rot_v
    # Tractions sigma_rr, sigma_rtheta
    srr_u, srr_v = lam * del_u + 2 * mu * d, lam * del_v
    srt_u, srt_v = -mu * n * rinv, mu * (d - rinv)
    k_outer, k_inner = 0, n_solid

    a = np.zeros((4 * ns, 4 * ns), dtype=complex)
    b = np.zeros((4 * ns, 4 * ns), dtype=complex)
    iu, iv, iwu, iwv = 0, ns, 2 * ns, 3 * ns
    # s u = w, rho_s s w = L u
    a[iu:iu + ns, iwu:iwu + ns] = eye
    b[iu:iu + ns, iu:iu + ns] = eye
    a[iv:iv + ns, iwv:iwv + ns] = eye
    b[iv:iv + ns, iv:iv + ns] = eye
    a[iwu:iwu + ns, iu:iu + ns] = lru
    a[iwu:iwu + ns, iv:iv + ns] = lrv
    b[iwu:iwu + ns, iwu:iwu + ns] = rho_s * eye
    a[iwv:iwv + ns, iu:iu + ns] = ltu
    a[iwv:iwv + ns, iv:iv + ns] = ltv
    b[iwv:iwv + ns, iwv:iwv + ns] = rho_s * eye
    # Traction conditions replace the momentum rows at the surfaces
    for row, su, sv in ((iwu, srr_u, srr_v), (iwv, srt_u, srt_v)):
        for k in (k_inner, k_outer):
            a[row + k, :] = 0
            b[row + k, :] = 0
            a[row + k, iu:iu + ns] = su[k]
            a[row + k, iv:iv + ns] = sv[k]
    rows = (iwu + k_outer, iwv + k_outer)
    return a, b, rows, rows


def bessel_series(z: complex, order: int, sign: int) -> complex:
    """sum_k sign^k a_k(order) / z^k of the large-argument expansions
    I_n(z) ~ e^z S(-1) / sqrt(2 pi z) and K_n(z) ~ sqrt(pi/(2z)) e^-z S(+1)."""
    total, term, k = 1.0 + 0j, 1.0 + 0j, 0
    while True:
        k += 1
        term *= sign * (4 * order ** 2 - (2 * k - 1) ** 2) / (k * 8 * z)
        if abs(term) < 1e-17 * abs(total) or k > 60:
            break
        total += term
    return total


def viscous_impedance(n: int, ro: float, outer: float, rho_f: float,
                      nu_f: float, s: complex) -> np.ndarray:
    """2x2 matrix Z with (sigma_rr, sigma_rtheta) = Z (w_r, w_theta) on the
    ring for fluid velocities w on the ring surface, at eigenvalue s."""
    kappa = np.sqrt(s / nu_f)
    if abs(kappa) * ro < 30:
        raise ValueError("Stokes layer too thick for the Bessel expansions")
    mu_f = rho_f * nu_f

    def f_i(r):
        """I_n(kappa r) / I_n(kappa b) and its r derivative."""
        z, zb = kappa * r, kappa * outer
        value = np.exp(z - zb) * np.sqrt(zb / z) \
            * bessel_series(z, n, -1) / bessel_series(zb, n, -1)
        ratio = bessel_series(z, n + 1, -1) / bessel_series(z, n, -1)
        return value, kappa * value * (ratio + n / z)

    def f_k(r):
        """K_n(kappa r) / K_n(kappa Ro) and its r derivative."""
        z, zo = kappa * r, kappa * ro
        value = np.exp(zo - z) * np.sqrt(zo / z) \
            * bessel_series(z, n, 1) / bessel_series(zo, n, 1)
        ratio = bessel_series(z, n + 1, 1) / bessel_series(z, n, 1)
        return value, kappa * value * (-ratio + n / z)

    def basis(r):
        """Psi, Psi' and D_n Psi of the four solutions at r."""
        p1 = (r / ro) ** n
        p2 = (ro / r) ** n
        fi, fi1 = f_i(r)
        fk, fk1 = f_k(r)
        psi = np.array([p1, p2, fi, fk])
        dpsi = np.array([n * p1 / r, -n * p2 / r, fi1, fk1])
        dn_psi = np.array([0, 0, kappa ** 2 * fi, kappa ** 2 * fk])
        dn_dpsi = np.array([0, 0, kappa ** 2 * fi1, kappa ** 2 * fk1])
        return psi, dpsi, dn_psi, dn_dpsi

    psi_b, dpsi_b, _, _ = basis(outer)
    psi_o, dpsi_o, dn_o, dn_d_o = basis(ro)
    # No slip on the rigid wall; u_r = n Psi / r and u_theta = -Psi' on the
    # ring
    m = np.array([psi_b, dpsi_b, n * psi_o / ro, -dpsi_o])
    z = np.zeros((2, 2), dtype=complex)
    for j in range(2):
        rhs = np.zeros(4, dtype=complex)
        rhs[2 + j] = 1.0
        coeff = np.linalg.solve(m, rhs)
        psi, dpsi = psi_o @ coeff, dpsi_o @ coeff
        omega = -(dn_o @ coeff)
        domega = -(dn_d_o @ coeff)
        d2psi = -omega - dpsi / ro + n * n * psi / ro ** 2
        f, g = n * psi / ro, -dpsi
        df = n * (dpsi / ro - psi / ro ** 2)
        dg = -d2psi
        # theta momentum: rho_f s G = n P / r + mu_f Omega'
        p = (ro / n) * (rho_f * s * g - mu_f * domega)
        z[0, j] = -p + 2 * mu_f * df
        z[1, j] = mu_f * (-n * f / ro + dg - g / ro)
    return z


def nearest_eigenvalue(a, b, guess: complex) -> complex:
    """Eigenvalue of A x = s B x closest to guess, by shift-invert."""
    mu_ev = np.linalg.eigvals(np.linalg.solve(a - guess * b, b))
    mu_ev = mu_ev[np.abs(mu_ev) > 1e-14]
    s = guess + 1.0 / mu_ev
    return complex(s[np.argmin(np.abs(s - guess))])


def eigenvalue(n: int, ri: float, ro: float, young: float, poisson: float,
               rho_s: float, rho_f: float = 0.0, outer: float | None = None,
               nu_f: float = 0.0, n_solid: int = 24,
               guess: float = 1.0) -> complex:
    """Eigenvalue s (exp(s t)) of the n-mode closest to s = i*guess."""
    a, b, rows, cols = solid_matrices(n, ri, ro, young, poisson, rho_s,
                                      n_solid)
    s = 1j * guess
    if rho_f > 0 and nu_f == 0:
        # sigma_rr(Ro) = -p = -m_a s^2 U(Ro) = -m_a s w_r(Ro)
        b[rows[0], cols[0]] = -added_mass(n, ro, outer, rho_f)
        return nearest_eigenvalue(a, b, s)
    if rho_f == 0:
        return nearest_eigenvalue(a, b, s)
    # The tractions are Z(s) w = s (Z(s)/s) w, where Z/s (the added mass in
    # the inviscid limit) depends only weakly on s, through the Stokes layer:
    # keep the factor s in the eigenproblem and iterate on Z/s
    for _ in range(100):
        z = viscous_impedance(n, ro, outer, rho_f, nu_f, s) / s
        b_s = b.copy()
        for i in range(2):
            for j in range(2):
                b_s[rows[i], cols[j]] = z[i, j]
        s_new = nearest_eigenvalue(a, b_s, s)
        if abs(s_new - s) < 1e-11 * abs(s):
            return s_new
        s = s_new
    raise RuntimeError("Viscous eigenvalue iteration did not converge")


def thin_ring(n: int, a: float, h: float, young: float, poisson: float,
              rho_s: float, rho_f: float = 0.0, radius: float = 1.0,
              outer: float = 2.0) -> tuple[float, float, float]:
    """Thin-ring dry frequency, wet/dry frequency ratio and added-mass ratio,
    with the fluid acting at `radius`."""
    rigidity = young * h ** 3 / (12 * (1 - poisson ** 2))
    omega_dry = math.sqrt(rigidity * n ** 2 * (n ** 2 - 1) ** 2
                          / (rho_s * h * a ** 4 * (n ** 2 + 1)))
    mass_ratio = 0.0
    if rho_f > 0:
        mass_ratio = added_mass(n, radius, outer, rho_f) * radius * n ** 2 \
            / (rho_s * h * a * (n ** 2 + 1))
    return omega_dry, 1.0 / math.sqrt(1.0 + mass_ratio), mass_ratio


def compute(spec: dict) -> dict:
    geo, mat = spec["geometry"], spec["material"]
    n = spec["mode"]
    ri, ro, outer = geo["innerRadius"], geo["outerRadius"], geo["wallRadius"]
    a, h = 0.5 * (ri + ro), ro - ri
    young, poisson, rho_s = mat["E"], mat["nu"], mat["rhoSolid"]
    nu_f = spec["fluid"]["nu"]

    def solve(**kwargs):
        # Two collocation resolutions: report the finer, check the change
        coarse = eigenvalue(n, ri, ro, young, poisson, rho_s, n_solid=16,
                            **kwargs)
        fine = eigenvalue(n, ri, ro, young, poisson, rho_s, n_solid=24,
                          **kwargs)
        return fine, abs(fine / coarse - 1)

    thin_dry, _, _ = thin_ring(n, a, h, young, poisson, rho_s)
    dry, change = solve(guess=thin_dry)
    out = {"dry": {"omega": dry.imag, "omegaThinRing": thin_dry,
                   "collocationChange": change},
           "levels": {}}
    for name, level in spec["levels"].items():
        rho_f = level["rhoFluid"]
        _, thin_ratio, mass_ratio = thin_ring(
            n, a, h, young, poisson, rho_s, rho_f, ro, outer
        )
        wet, change = solve(rho_f=rho_f, outer=outer,
                            guess=dry.imag * thin_ratio)
        visc, visc_change = solve(rho_f=rho_f, outer=outer, nu_f=nu_f,
                                  guess=wet.imag)
        out["levels"][name] = {
            "massRatio": mass_ratio,
            "omega": wet.imag,
            "ratio": wet.imag / dry.imag,
            "ratioThinRing": thin_ratio,
            "collocationChange": change,
            "viscousOmega": visc.imag,
            "viscousFrequencyShift": visc.imag / wet.imag - 1,
            "viscousDampingRatio": -visc.real / abs(visc),
            "viscousCollocationChange": visc_change,
            "stokesLayerOverGap": math.sqrt(2 * nu_f / wet.imag) / (outer - ro),
        }
    return out


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--json", action="store_true",
                        help="print the 'exact' block of the reference JSON")
    parser.add_argument("--check", action="store_true",
                        help="check the stored values against a fresh solve")
    args = parser.parse_args()
    spec = json.loads(REFERENCE_FILE.read_text())
    exact = compute(spec)
    if args.json:
        print(json.dumps(exact, indent=4))
        return 0
    print(f"dry: omega = {exact['dry']['omega']:.8f} rad/s "
          f"(thin ring {exact['dry']['omegaThinRing']:.8f})")
    for name, level in exact["levels"].items():
        print(f"{name}: m_a/m_s = {level['massRatio']:.4f}, "
              f"omega = {level['omega']:.8f}, ratio = {level['ratio']:.8f} "
              f"(thin ring {level['ratioThinRing']:.8f}); viscous: frequency "
              f"shift {level['viscousFrequencyShift']:.2e}, damping ratio "
              f"{level['viscousDampingRatio']:.2e}")
    if args.check:
        stored = spec["exact"]
        bad = []
        if abs(stored["dry"]["omega"] / exact["dry"]["omega"] - 1) > 1e-6:
            bad.append("dry.omega")
        for name, level in exact["levels"].items():
            for key in ("omega", "ratio", "ratioThinRing", "viscousOmega"):
                if abs(stored["levels"][name][key] / level[key] - 1) > 1e-6:
                    bad.append(f"{name}.{key}")
        if bad:
            print("MISMATCH: " + ", ".join(bad))
            return 1
        print("The stored exact values agree with a fresh solve")
    return 0


if __name__ == "__main__":
    sys.exit(main())
