#!/usr/bin/env python3
"""Exact linear solution for pulsatile flow in an elastic tube.

A harmonic pressure wave travels along an infinitely long, straight elastic
tube, with every field proportional to exp(i (omega t - k x)), x being the
axial coordinate and k the complex wave number. Two solutions are computed.

Womersley's thin-walled tube (the classical reference). The fluid is the
long-wave solution in a tube of radius R,

    u_x(r) = k P / (rho omega) [1 - J0(Lambda r/R) / J0(Lambda)]
             + i omega xi J0(Lambda r/R) / J0(Lambda),
    Lambda = i^(3/2) alpha,  alpha = R sqrt(omega / nu),
    g = F10 = 2 J1(Lambda) / (Lambda J0(Lambda)),

xi being the axial wall displacement, and the wall is a thin elastic membrane.
For a free (untethered) wall the frequency equation is Womersley's quadratic
(Filonova et al. 2020, Eq. 41)

    (1 - g)(1 - nu^2) v^2 - [2 + m (1 - g) + g (1/2 - 2 nu)] v + g + 2 m = 0,
    v = E h / ((1 - nu^2) rho R c^2),  m = rho_s h / (rho R),

whose root nearer the tethered value is the pressure wave, c = omega/k. For a
longitudinally tethered wall (xi = 0) it reduces to

    k^2 = 2 rho omega^2 / (R (1 - g) (K - rho_s h omega^2)),
    K = E h / ((1 - nu^2) R^2).

The exact continuum solution (the reference for the simulations). The fluid
is the linearised incompressible Navier-Stokes solution in the whole tube,
not the long-wave one,

    p   = A I0(k r),
    u_x = k A / (rho omega) I0(k r) + B I0(s r),
    u_r = i k A / (rho omega) I1(k r) + i k B I1(s r) / s,
    s^2 = k^2 + i omega / nu,

and the wall is a linear elastic cylinder Ri < r < Ro of finite thickness,
whose radial problem is solved by Chebyshev collocation. The fluid and the
wall share the velocity and the traction at r = Ri (linearised about the
undeformed position). The outer surface is traction free, or free of normal
traction with u_x = 0 when tethered. The wave number is the root of the
resulting frequency equation, found by the secant method from the thin-wall
value. As h/R -> 0 at fixed E h and rho_s h, the exact wave number tends to
Womersley's, with an error proportional to h/R (see --thickness).

The Bessel functions of complex argument are evaluated by their power series,
which is accurate for the moderate arguments used here. This script needs
numpy; the verification driver does not, as it reads the values this script
stores in reference/womersleyTube_verification_references.json.

usage: womersley_exact.py [--json] [--check] [--write-case DIR]
"""

from __future__ import annotations

import argparse
import cmath
import json
import math
import sys
from pathlib import Path

import numpy as np

SCRIPT = Path(__file__).resolve()
REFERENCE_FILE = (
    SCRIPT.parents[1] / "reference" / "womersleyTube_verification_references.json"
)

# Tutorial parameters (SI units)
PARAMETERS = {
    "innerRadius": 1.0,
    "thickness": 0.1,
    "length": 15.0,
    "youngsModulus": 2.0e4,
    "poissonsRatio": 0.3,
    "solidDensity": 1000.0,
    "fluidDensity": 1000.0,
    "dynamicViscosity": 5.0,
    "frequency": 0.02,
    "pressureAmplitude": 1.0,
    "tethered": False,
}

# Chebyshev collocation intervals across the wall; the wave number is
# converged to round-off with 8
COLLOCATION = 16

FOAM_HEADER = r"""/*--------------------------------*- C++ -*----------------------------------*\
| solids4foam: solid mechanics and fluid-solid interaction simulations        |
| Version:     v2.4                                                           |
| Web:         https://solids4foam.github.io                                  |
| Disclaimer:  This offering is not approved or endorsed by OpenCFD Limited,  |
|              producer and distributor of the OpenFOAM software via          |
|              www.openfoam.com, and owner of the OPENFOAM® and OpenCFD®      |
|              trade marks.                                                   |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      womersleyProperties;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
"""


def bessel(n: int, z: complex, modified: bool) -> complex:
    """J_n(z) or I_n(z), n = 0 or 1, by the power series."""
    z = complex(z)
    term = (z / 2) ** n / math.factorial(n)
    total = term
    q = (z / 2) ** 2 * (1 if modified else -1)
    for m in range(1, 500):
        term *= q / (m * (m + n))
        total += term
        if abs(term) <= 1e-17 * abs(total):
            return total
    raise RuntimeError(f"Bessel series did not converge for z = {z}")


def chebyshev(n: int, a: float, b: float):
    """Chebyshev-Gauss-Lobatto points on [a, b], first point b, and the
    derivative matrix."""
    x = np.cos(np.pi * np.arange(n + 1) / n)
    c = np.hstack([2.0, np.ones(n - 1), 2.0]) * (-1.0) ** np.arange(n + 1)
    dx = np.subtract.outer(x, x)
    d = np.outer(c, 1 / c) / (dx + np.eye(n + 1))
    d -= np.diag(d.sum(axis=1))
    return a + 0.5 * (x + 1) * (b - a), d * 2 / (b - a)


def cvec(z: complex) -> list[float]:
    return [z.real, z.imag]


class Tube:
    def __init__(self, prm: dict, collocation: int = COLLOCATION):
        self.prm = dict(prm)
        self.Ri = prm["innerRadius"]
        self.h = prm["thickness"]
        self.Ro = self.Ri + self.h
        self.E = prm["youngsModulus"]
        self.nu = prm["poissonsRatio"]
        self.rhos = prm["solidDensity"]
        self.rho = prm["fluidDensity"]
        self.mu = prm["dynamicViscosity"]
        self.nuf = self.mu / self.rho
        self.omega = 2 * math.pi * prm["frequency"]
        self.P = prm["pressureAmplitude"]
        self.tethered = bool(prm["tethered"])
        self.N = collocation
        self.G = self.E / (2 * (1 + self.nu))
        self.lam = self.E * self.nu / ((1 + self.nu) * (1 - 2 * self.nu))
        self.alpha = self.Ri * math.sqrt(self.omega / self.nuf)
        self.k = None

    # -- Womersley's thin-walled tube ---------------------------------------

    def F10(self) -> complex:
        lam = self.alpha * cmath.exp(3j * math.pi / 4)
        return 2 * bessel(1, lam, False) / (lam * bessel(0, lam, False))

    def womersley_tethered_k(self) -> complex:
        K = self.E * self.h / ((1 - self.nu**2) * self.Ri**2)
        k = cmath.sqrt(2 * self.rho * self.omega**2 / (
            self.Ri * (1 - self.F10()) * (K - self.rhos * self.h * self.omega**2)))
        return k if k.real > 0 else -k

    def womersley_free_k(self) -> complex:
        nu, g = self.nu, self.F10()
        m = self.rhos * self.h / (self.rho * self.Ri)
        roots = np.roots([(1 - g) * (1 - nu**2),
                          -(2 + m * (1 - g) + g * (0.5 - 2 * nu)),
                          g + 2 * m])
        ks = []
        for v in roots:
            c = cmath.sqrt(self.E * self.h / ((1 - nu**2) * self.rho * self.Ri * v))
            k = self.omega / c
            ks.append(k if k.real > 0 else -k)
        return min(ks, key=lambda k: abs(k - self.womersley_tethered_k()))

    def womersley_k(self) -> complex:
        return self.womersley_tethered_k() if self.tethered \
            else self.womersley_free_k()

    # -- Exact continuum solution -------------------------------------------

    def s(self, k: complex) -> complex:
        s = cmath.sqrt(k * k + 1j * self.omega / self.nuf)
        return s if s.real > 0 else -s

    def fluid_rows(self, k: complex, r: float) -> dict:
        """(p, u_x, u_r, sigma_rr, sigma_rx) at r > 0 per unit (A, B)."""
        w, rho, mu, s = self.omega, self.rho, self.mu, self.s(k)
        I0k, I1k = bessel(0, k * r, True), bessel(1, k * r, True)
        I0s, I1s = bessel(0, s * r, True), bessel(1, s * r, True)
        dI1k, dI1s = I0k - I1k / (k * r), I0s - I1s / (s * r)
        ux = np.array([k / (rho * w) * I0k, I0s])
        ur = np.array([1j * k / (rho * w) * I1k, 1j * k / s * I1s])
        dur = np.array([1j * k * k / (rho * w) * dI1k, 1j * k * dI1s])
        dux = np.array([k * k / (rho * w) * I1k, s * I1s])
        p = np.array([I0k, 0])
        return {"p": p, "ux": ux, "ur": ur,
                "srr": -p + 2 * mu * dur, "srx": mu * (dux - 1j * k * ur)}

    def matrix(self, k: complex) -> np.ndarray:
        """Collocation equations for (U, W, A, B), U and W being the radial
        and axial wall displacements at the Chebyshev points (Ro first)."""
        N, n = self.N, self.N + 1
        _, d = chebyshev(N, self.Ri, self.Ro)
        r = chebyshev(N, self.Ri, self.Ro)[0]
        rinv, eye = np.diag(1 / r), np.eye(n)
        lam2G, G, rw2 = self.lam + 2 * self.G, self.G, self.rhos * self.omega**2
        # Dilatation U' + U/r - i k W and rotation -i k U - W'
        divU, divW = d + rinv, -1j * k * eye
        rotU, rotW = -1j * k * eye, -d
        rows = []

        def row(u=None, w=None, fluid=None, u_at=None, w_at=None, value=1):
            out = np.zeros(2 * n + 2, dtype=complex)
            if u is not None:
                out[:n] = u
            if w is not None:
                out[n:2 * n] = w
            if u_at is not None:
                out[u_at] = value
            if w_at is not None:
                out[n + w_at] = value
            if fluid is not None:
                out[2 * n:] = -fluid
            rows.append(out)

        # Navier equations, rho_s omega^2 u + (lambda + 2G) grad div u
        # - G curl curl u = 0
        radialU = rw2 * eye + lam2G * d @ divU - G * 1j * k * rotU
        radialW = lam2G * d @ divW - G * 1j * k * rotW
        axialU = -lam2G * 1j * k * divU - G * (d + rinv) @ rotU
        axialW = rw2 * eye - lam2G * 1j * k * divW - G * (d + rinv) @ rotW
        srrU, srrW = self.lam * divU + 2 * G * d, self.lam * divW
        srxU, srxW = -1j * k * G * eye, G * d
        for i in range(1, N):
            row(radialU[i], radialW[i])
            row(axialU[i], axialW[i])
        # Outer surface: no normal traction, and tethered or free of shear
        row(srrU[0], srrW[0])
        if self.tethered:
            row(w_at=0)
        else:
            row(srxU[0], srxW[0])
        # Interface: traction and velocity continuity
        f = self.fluid_rows(k, self.Ri)
        row(srrU[N], srrW[N], f["srr"])
        row(srxU[N], srxW[N], f["srx"])
        row(w_at=N, value=1j * self.omega, fluid=f["ux"])
        row(u_at=N, value=1j * self.omega, fluid=f["ur"])
        return np.array(rows)

    def residual(self, k: complex) -> tuple[complex, np.ndarray]:
        """Solve with A = 1 without the last equation and return that
        equation's residual, which vanishes at a root."""
        m = self.matrix(k)
        # Equilibrate: the rows mix stiffness, traction and velocity scales
        m /= np.abs(m).max(axis=1, keepdims=True)
        iA = 2 * (self.N + 1)
        cols = [j for j in range(m.shape[1]) if j != iA]
        x = np.linalg.solve(m[:-1][:, cols], -m[:-1, iA])
        x = np.insert(x, iA, 1.0)
        return m[-1] @ x, x

    def solve(self) -> complex:
        k0 = self.womersley_k()
        k1 = k0 * (1 + 1e-4)
        f0, f1 = self.residual(k0)[0], self.residual(k1)[0]
        # The secant iterates settle to round-off, about 1e-10 relative
        for _ in range(100):
            k0, k1 = k1, k1 - f1 * (k1 - k0) / (f1 - f0)
            f0, f1 = f1, self.residual(k1)[0]
            if abs(k1 - k0) <= 1e-10 * abs(k1) or f1 == f0:
                break
        else:
            raise RuntimeError("Secant iteration did not converge")
        self.k = k1
        _, x = self.residual(k1)
        n = self.N + 1
        # Amplitudes for a pressure amplitude P on the axis at x = 0
        self.U = self.P * x[:n]
        self.W = self.P * x[n:2 * n]
        self.A, self.B = self.P * x[2 * n], self.P * x[2 * n + 1]
        return k1

    def flow_rate(self) -> complex:
        """Q = 2 pi int_0^Ri u_x r dr at x = 0."""
        k, s, R = self.k, self.s(self.k), self.Ri
        return 2 * math.pi * R * (
            self.A / (self.rho * self.omega) * bessel(1, k * R, True)
            + self.B * bessel(1, s * R, True) / s)

    # -- Output ----------------------------------------------------------------

    def properties(self) -> str:
        """constant/womersleyProperties, read by the coded boundary
        conditions and initial fields. The pressure coefficient is
        kinematic (p/rho), as the pimpleFluid pressure is."""
        def c(z: complex) -> str:
            return f"({z.real:.16e} {z.imag:.16e})"
        lines = [
            FOAM_HEADER.rstrip("\n"),
            "",
            "// Generated by verification/scripts/womersley_exact.py from the",
            "// parameters in that script: do not edit by hand",
            "",
            f"omega           {self.omega:.16e};",
            f"k               {c(self.k)};",
            f"s               {c(self.s(self.k))};",
            f"A               {c(self.A / self.rho)};",
            f"B               {c(self.B)};",
            f"innerRadius     {self.Ri:.16e};",
            f"outerRadius     {self.Ro:.16e};",
            "",
            "// Radial and axial wall displacement at the Chebyshev points",
            "// r = Ri + (1 + cos(pi j/N)) (Ro - Ri)/2, j = 0..N",
            "solidRadialDisplacement",
            "(",
        ]
        lines += [f"    {c(complex(v))}" for v in self.U]
        lines += [");", "", "solidAxialDisplacement", "("]
        lines += [f"    {c(complex(v))}" for v in self.W]
        lines += [");", "", "// " + "*" * 73 + " //"]
        return "\n".join(lines) + "\n"


def wave(tube: Tube, k: complex) -> dict:
    return {"k": cvec(k), "waveSpeed": tube.omega / k.real,
            "wavelength": 2 * math.pi / k.real,
            "amplitudeRatioPerWavelength":
                math.exp(2 * math.pi * k.imag / k.real)}


def thickness_study(base: dict) -> list[dict]:
    """h/R sequence at fixed E h and rho_s h, for the free and the
    tethered tube."""
    rows = []
    for ratio in (0.2, 0.1, 0.05, 0.02, 0.01, 0.005):
        entry = {"hOverR": ratio}
        for tethered in (False, True):
            prm = dict(base, tethered=tethered)
            h = ratio * base["innerRadius"]
            prm["thickness"] = h
            prm["youngsModulus"] = base["youngsModulus"] * base["thickness"] / h
            prm["solidDensity"] = base["solidDensity"] * base["thickness"] / h
            tube = Tube(prm)
            k = tube.solve()
            key = "tethered" if tethered else "free"
            entry[key] = {"k": cvec(k),
                          "kOverWomersley": cvec(k / tube.womersley_k())}
        rows.append(entry)
    return rows


def reference() -> dict:
    tube = Tube(PARAMETERS)
    k = tube.solve()
    tethered = Tube(dict(PARAMETERS, tethered=True))
    kt = tethered.solve()
    c0 = math.sqrt(tube.E * tube.h / (2 * tube.rho * tube.Ri))
    R = tube.Ri
    return {
        "parameters": PARAMETERS,
        "derived": {
            "omega": tube.omega,
            "period": 1 / PARAMETERS["frequency"],
            "womersleyNumber": tube.alpha,
            "F10": cvec(tube.F10()),
            "moensKortewegSpeed": c0,
            "massRatio": tube.rhos * tube.h / (tube.rho * R),
        },
        "exact": {
            **wave(tube, k),
            "s": cvec(tube.s(k)),
            # Pressure p = A I0(k r) and velocity coefficients, in Pa and m/s
            "A": cvec(tube.A),
            "B": cvec(tube.B),
            "wallRadialDisplacement": cvec(complex(tube.U[-1])),
            "wallAxialDisplacement": cvec(complex(tube.W[-1])),
            "flowRate": cvec(tube.flow_rate()),
        },
        "womersleyThinWall": wave(tube, tube.womersley_k()),
        "tetheredExact": wave(tethered, kt),
        "tetheredWomersleyThinWall": wave(tethered, tethered.womersley_k()),
        "collocationCheck": {
            "k8": cvec(Tube(PARAMETERS, 8).solve()),
            "k24": cvec(Tube(PARAMETERS, 24).solve()),
        },
        "thickness": thickness_study(PARAMETERS),
    }


def compare(a, b, path: str = "") -> list[str]:
    """Differences between two reference trees beyond 1e-8 relative."""
    if isinstance(a, dict) and isinstance(b, dict):
        out = []
        for key in sorted(set(a) | set(b)):
            if key not in a or key not in b:
                out.append(f"{path}/{key}: missing")
            else:
                out += compare(a[key], b[key], f"{path}/{key}")
        return out
    if isinstance(a, list) and isinstance(b, list):
        if len(a) != len(b):
            return [f"{path}: length {len(a)} != {len(b)}"]
        return [e for i, (x, y) in enumerate(zip(a, b))
                for e in compare(x, y, f"{path}[{i}]")]
    if isinstance(a, bool) or not isinstance(a, (int, float)):
        return [] if a == b else [f"{path}: {a!r} != {b!r}"]
    if abs(a - b) > 1e-8 * max(abs(a), abs(b), 1e-300):
        return [f"{path}: {a!r} != {b!r}"]
    return []


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--json", action="store_true",
                        help="print the reference values as JSON")
    parser.add_argument("--check", action="store_true",
                        help="compare with the stored reference file")
    parser.add_argument("--write-case", metavar="DIR",
                        help="write DIR/constant/womersleyProperties")
    args = parser.parse_args()
    if args.write_case:
        tube = Tube(PARAMETERS)
        tube.solve()
        path = Path(args.write_case) / "constant" / "womersleyProperties"
        path.write_text(tube.properties())
        print(f"Wrote {path}")
    data = reference() if (args.json or args.check or not args.write_case) \
        else None
    if args.check:
        stored = json.loads(REFERENCE_FILE.read_text())["reference"]
        differences = compare(data, stored)
        for line in differences:
            print(line)
        print("Reference file up to date" if not differences
              else f"{len(differences)} differences")
        return 1 if differences else 0
    if data is not None:
        print(json.dumps(data, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
