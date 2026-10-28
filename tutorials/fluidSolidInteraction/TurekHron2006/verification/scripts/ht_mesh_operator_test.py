#!/usr/bin/env python3
"""Discrete-operator test of 2-D meshes with a smooth quadratic field.

For phi = sin(a x) cos(b y) style smooth field (here a quadratic plus a
trigonometric part, so that no operator is exact) evaluated at the cell
centroids, computes on the internal faces / cells:

  e_int   face value of OpenFOAM's linear interpolation (weights from the
          face-normal projections, surfaceInterpolation::makeWeights)
          minus phi at the face centre;
  e_grad  cell gradient of OpenFOAM's leastSquares scheme (weights 1/|d|^2,
          boundary faces included with the exact face value) minus the
          exact gradient;
  e_div   truncation error of Gauss linear divergence of the exact vector
          field grad(phi) (cell values, linear face interpolation) against the
          exact Laplacian of phi, over internal cells (all faces internal);
  e_lap   truncation error of the Gauss linear corrected Laplacian of phi
          (orthogonal part (phi_N - phi_O)/(n.d) plus the non-orthogonal
          correction with the linearly interpolated leastSquares gradient);

the error sets: "corner" (within --radius of the four valence-5 ring corners
and the two flag-root cylinder vertices), "roi" (x 0.1-1.0 m, y 0.08-0.32 m,
not corner: cylinder ring, flag and near wake) and "rest" (elsewhere). For the
divergence and Laplacian only cells without boundary faces are used.

reported as max and rms, separately for the "corner" set (faces/cells within
--radius of the four valence-5 ring corners and the two flag-root cylinder
vertices) and for the rest, together with the observed orders between the
cases given (which must be one family at successive levels L, passed with
--levels).

    python3 ht_mesh_operator_test.py --levels 1,2,4 sm_1 sm_2 sm_4
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np

import ht_mesh_audit as au

CORNERS = np.array([[0.12, 0.12], [0.28, 0.12], [0.12, 0.28], [0.28, 0.28],
                    [0.2353553, 0.1646447], [0.2353553, 0.2353553]])


def field(x, y):
    k = 2 * math.pi / 0.2
    f = np.sin(k * x) * np.cos(k * y) + 4.0 * (x - 0.3) ** 2
    gx = k * np.cos(k * x) * np.cos(k * y) + 8.0 * (x - 0.3)
    gy = -k * np.sin(k * x) * np.sin(k * y)
    return f, gx, gy


def laplacian(x, y):
    k = 2 * math.pi / 0.2
    return -2 * k * k * np.sin(k * x) * np.cos(k * y) + 8.0


def geometry(case: Path):
    d = case / "constant/polyMesh"
    P, F = au.read_points(d / "points"), au.read_faces(d / "faces")
    O, N = au.read_labels(d / "owner"), au.read_labels(d / "neighbour")
    nc = int(O.max()) + 1
    C = np.zeros((nc, 2))
    V = np.zeros(nc)
    zmax = P[:, 2].max()
    side = []  # in-plane faces: (index, centre, normal*length)
    for fi, f in enumerate(F):
        q = P[f]
        if np.allclose(q[:, 2], zmax):
            x, y = q[:, 0], q[:, 1]
            cr = x * np.roll(y, -1) - np.roll(x, -1) * y
            a6 = 3.0 * cr.sum()
            C[O[fi]] = [((x + np.roll(x, -1)) * cr).sum() / a6, ((y + np.roll(y, -1)) * cr).sum() / a6]
            V[O[fi]] = abs(a6) / 6.0
        elif np.ptp(q[:, 2]) > 0:
            e = q[np.argsort(q[:, 2])][:2, :2]
            t = e[1] - e[0]
            side.append((fi, 0.5 * (e[0] + e[1]), np.array([t[1], -t[0]])))
    return C, V, O, N, side


def test(case: Path, radius: float):
    C, V, O, N, side = geometry(case)
    ni = len(N)
    phiC = field(C[:, 0], C[:, 1])[0]
    idx = np.array([s[0] for s in side])
    Cf = np.array([s[1] for s in side])
    Sf = np.array([s[2] for s in side])
    o = O[idx]
    # orient Sf from owner to the outside
    sgn = np.sign(np.einsum("ij,ij->i", Sf, Cf - C[o]))
    Sf = Sf * sgn[:, None]
    internal = idx < ni
    # interpolation error on internal faces
    fi = idx[internal]
    on, nn = O[fi], N[fi]
    sf, cf = Sf[internal], Cf[internal]
    dO = np.einsum("ij,ij->i", sf, cf - C[on])
    dN = np.einsum("ij,ij->i", sf, C[nn] - cf)
    w = dN / (dO + dN)
    e_int = w * phiC[on] + (1 - w) * phiC[nn] - field(cf[:, 0], cf[:, 1])[0]
    # least-squares gradient
    nc = len(C)
    G = np.zeros((nc, 2, 2))
    b = np.zeros((nc, 2))
    d = C[nn] - C[on]
    wt = 1.0 / np.einsum("ij,ij->i", d, d)
    dd = np.einsum("i,ij,ik->ijk", wt, d, d)
    dphi = phiC[nn] - phiC[on]
    np.add.at(G, on, dd)
    np.add.at(G, nn, dd)
    np.add.at(b, on, (wt * dphi)[:, None] * d)
    np.add.at(b, nn, (wt * dphi)[:, None] * d)
    bf = ~internal
    ob = o[bf]
    db = Cf[bf] - C[ob]
    wb = 1.0 / np.einsum("ij,ij->i", db, db)
    np.add.at(G, ob, np.einsum("i,ij,ik->ijk", wb, db, db))
    np.add.at(b, ob, (wb * (field(Cf[bf, 0], Cf[bf, 1])[0] - phiC[ob]))[:, None] * db)
    g = np.linalg.solve(G, b[..., None])[..., 0]
    _, gx, gy = field(C[:, 0], C[:, 1])
    e_grad = np.hypot(g[:, 0] - gx, g[:, 1] - gy)
    # divergence of the exact gradient field with linear face interpolation
    ux, uy = gx, gy
    uf = np.stack([w * ux[on] + (1 - w) * ux[nn], w * uy[on] + (1 - w) * uy[nn]], axis=1)
    flux = np.einsum("ij,ij->i", uf, sf)
    div = np.zeros(nc)
    np.add.at(div, on, flux)
    np.add.at(div, nn, -flux)
    # Gauss linear corrected Laplacian
    mag = np.linalg.norm(sf, axis=1)
    nh = sf / mag[:, None]
    nd = np.maximum(np.einsum("ij,ij->i", nh, d), 0.05 * np.linalg.norm(d, axis=1))
    kcorr = nh - d / nd[:, None]
    gf = w[:, None] * g[on] + (1 - w)[:, None] * g[nn]
    sn = (phiC[nn] - phiC[on]) / nd + np.einsum("ij,ij->i", kcorr, gf)
    lflux = sn * mag
    lap = np.zeros(nc)
    np.add.at(lap, on, lflux)
    np.add.at(lap, nn, -lflux)
    inner = np.ones(nc, bool)
    inner[ob] = False
    exact = laplacian(C[:, 0], C[:, 1])
    e_div = np.abs(div / V - exact)
    e_lap = np.abs(lap / V - exact)

    def near(xy):
        return np.min(np.linalg.norm(xy[:, None, :] - CORNERS[None], axis=2), axis=1) < radius
    def roi(xy):
        return (xy[:, 0] > 0.1) & (xy[:, 0] < 1.0) & (xy[:, 1] > 0.08) & (xy[:, 1] < 0.32)
    fn, cn = near(cf), near(C)
    fr, cr = roi(cf) & ~fn, roi(C) & ~cn
    stats = {}
    for name, e, m, r, valid in (("int", np.abs(e_int), fn, fr, np.ones(len(fn), bool)),
                                 ("grad", e_grad, cn, cr, np.ones(nc, bool)),
                                 ("div", e_div, cn, cr, inner), ("lap", e_lap, cn, cr, inner)):
        for tag, sel in (("corner", m & valid), ("roi", r & valid), ("rest", ~m & ~r & valid)):
            stats[f"{name}_{tag}_max"] = float(e[sel].max())
            stats[f"{name}_{tag}_rms"] = float(np.sqrt(np.mean(e[sel] ** 2)))
    # where the roi maximum of the divergence and Laplacian errors sits
    for name, e in (("div", e_div), ("lap", e_lap)):
        sel = np.where(cr & inner)[0]
        i = sel[int(np.argmax(e[sel]))]
        stats[f"{name}_roi_max_at"] = [float(C[i, 0]), float(C[i, 1])]
    stats["cells"] = nc
    return stats


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("cases", nargs="+", type=Path)
    ap.add_argument("--levels", required=True)
    ap.add_argument("--radius", type=float, default=0.02)
    ap.add_argument("--json", type=Path)
    a = ap.parse_args()
    lv = [float(x) for x in a.levels.split(",")]
    res = [test(c, a.radius) for c in a.cases]
    keys = [k for k in res[0] if k != "cells" and not k.endswith("_at")]
    print(f"{'quantity':22s}" + "".join(f"{str(c.name):>12s}" for c in a.cases) + "   orders")
    for k in keys:
        v = [r[k] for r in res]
        ords = [math.log(v[i] / v[i + 1]) / math.log(lv[i + 1] / lv[i]) for i in range(len(v) - 1)]
        print(f"{k:22s}" + "".join(f"{x:12.3e}" for x in v) + "   " + " ".join(f"{p:5.2f}" for p in ords))
    for c, r in zip(a.cases, res):
        print(f"{c.name}: div roi max at {r['div_roi_max_at']}, lap roi max at {r['lap_roi_max_at']}")
    if a.json:
        a.json.write_text(json.dumps({str(c): r for c, r in zip(a.cases, res)}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
