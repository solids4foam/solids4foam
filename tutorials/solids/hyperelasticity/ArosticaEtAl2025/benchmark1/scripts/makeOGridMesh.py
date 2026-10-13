#!/usr/bin/env python3
"""
Write an all-hexahedral OpenFOAM mesh (constant/polyMesh) of the truncated
ellipsoid of Arostica et al. (2025), benchmark 1, with an O-grid at the apex.

Why not the revolved slice of LandEtAl2015/problem3: revolving a meridian
slice about the long axis puts a fan of wedge cells, all meeting on the axis,
at the apex. Those cells are badly shaped for a cell-centred finite-volume
gradient, and in this benchmark (strong apical compression, rapid relaxation)
they are where the solution first loses accuracy and finally inverts.

Construction, on the (s, theta) parameter plane of each surface, where s is
the normalised arc length of the meridian from the apex (s = 0) to the base
(s = 1) and theta the angle about the long (x) axis:
  * the plane (u, v) = s (cos theta, sin theta) is the unit disk, apex at
    the centre and base on the rim;
  * a central square of nCap x nCap cells (no singular point at the apex);
  * a ring of 4 nCap x nRing cells whose rows blend from the square's
    boundary to the circle s = 1 (the base);
  * nWall layers through the wall, by straight-line interpolation between the
    endocardial and epicardial points with the same (s, theta).
Both surfaces are exact ellipsoids:
    x = -r_long cos(phi),  r = r_short sin(phi),  phi in [0, phi_base]
and the base, joining the two rims, is the straight conical surface of the
reference gmsh geometry.

Patches: endo, epi, base.

Usage:
    python3 scripts/makeOGridMesh.py [level]     # 1..4, default 2
Needs numpy.
"""
import math
import os
import sys

import numpy as np

R_SHORT_ENDO, R_LONG_ENDO = 0.025, 0.090
R_SHORT_EPI, R_LONG_EPI = 0.035, 0.097
PHI_BASE_ENDO = math.pi - math.acos(5.0 / 17.0)
PHI_BASE_EPI = math.pi - math.acos(5.0 / 20.0)

# nCap (cells per side of the apical square; 4 nCap around), nRing (cells
# from the square to the base), nWall (through the wall)
LEVELS = {
    1: (8, 18, 3),
    2: (16, 36, 5),
    3: (32, 72, 10),
    4: (64, 144, 20),
}
C_SQ = 0.18    # half-width of the apical square, in units of the meridian
SMOOTH_ITERS = 0      # optional Laplacian smoothing sweeps of the apical grid
SMOOTH_ROWS = 0.35    # fraction of the ring rows that may be smoothed
TRANSITION = 0.45     # fraction of the ring over which square -> circle
SHAPE_POWER = 1.0     # < 1 rounds the rows off the square sooner (less robust)


class Meridian:
    """phi(s) for an ellipse meridian parametrised by normalised arc length."""

    def __init__(self, rs, rl, phib, n=20001):
        phi = np.linspace(0.0, phib, n)
        dx = rl * np.sin(phi)
        dr = rs * np.cos(phi)
        ds = np.sqrt(dx**2 + dr**2)
        s = np.concatenate([[0.0], np.cumsum(0.5 * (ds[1:] + ds[:-1]) * np.diff(phi))])
        self.s = s / s[-1]
        self.phi = phi
        self.rs, self.rl = rs, rl

    def point(self, s, theta):
        phi = np.interp(s, self.s, self.phi)
        x = -self.rl * np.cos(phi)
        r = self.rs * np.sin(phi)
        return np.stack([x, r * np.cos(theta), r * np.sin(theta)], axis=-1)


def surface_params(nCap, nRing):
    """(s, theta) of the surface nodes, and the quads joining them.

    The (s, theta) plane is treated as a unit disk, (u, v) = s (cos, sin)
    theta, with the apex at its centre and the base on its rim. It is meshed
    as a classic O-grid: a square of half-width C_SQ in the middle, and a ring
    whose rows blend from the square's boundary to the circle s = 1.
    """
    n = nCap
    nAround = 4 * n
    st = []
    index = {}

    def add(key, u, v):
        index[key] = len(st)
        s = math.hypot(u, v)
        st.append((s, math.atan2(v, u) if s > 0 else 0.0))

    X = np.linspace(-C_SQ, C_SQ, n + 1)
    for i in range(n + 1):
        for j in range(n + 1):
            add(("c", i, j), X[i], X[j])

    # Square boundary, counter-clockwise starting at (C, -C)
    ring0 = [("c", n, j) for j in range(0, n)] \
        + [("c", i, n) for i in range(n, 0, -1)] \
        + [("c", 0, j) for j in range(n, 0, -1)] \
        + [("c", i, 0) for i in range(0, n)]
    Q = np.array([[X[k[1]], X[k[2]]] for k in ring0])

    rows = [ring0]
    thq = np.unwrap(np.arctan2(Q[:, 1], Q[:, 0]))
    thu = thq[0] + 2.0 * math.pi * np.arange(nAround) / nAround
    for k in range(1, nRing + 1):
        w = k / nRing
        # Smooth transition over the first TRANSITION of the ring: the rows
        # are superellipses |u|^p + |v|^p = r^p, from a square (p = inf) to a
        # circle (p = 2), and the node angles go from the square's to uniform
        x = min(1.0, w / TRANSITION)
        bShape = math.sin(0.5 * math.pi * x) ** SHAPE_POWER   # rounds early
        b = x * x * (3.0 - 2.0 * x)                  # angles: smoothstep
        p = 2.0 / max(bShape, 1e-12)                 # inf -> 2
        r = C_SQ + (1.0 - C_SQ) * w
        keys = []
        for m in range(nAround):
            th = thq[m] + b * (thu[m] - thq[m])
            d = np.array([math.cos(th), math.sin(th)])
            if p > 200:
                norm = np.abs(d).max()
            else:
                norm = (np.abs(d) ** p).sum() ** (1.0 / p)
            P = r * d / norm
            key = ("r", k, m)
            add(key, P[0], P[1])
            keys.append(key)
        rows.append(keys)

    quads = []
    for i in range(n):
        for j in range(n):
            quads.append([index[("c", i, j)], index[("c", i + 1, j)],
                          index[("c", i + 1, j + 1)], index[("c", i, j + 1)]])
    baseEdges = []
    for k in range(nRing):
        for m in range(nAround):
            m1 = (m + 1) % nAround
            a, b = rows[k][m], rows[k][m1]
            c, d = rows[k + 1][m1], rows[k + 1][m]
            quads.append([index[a], index[b], index[c], index[d]])
            if k == nRing - 1:
                baseEdges.append((index[d], index[c]))
    quads = np.array(quads)

    # Smooth the algebraic grid: Laplacian smoothing of every node but those
    # on the rim, in the (u, v) plane. This removes the kinks of the blend and
    # leaves a smooth, nearly orthogonal O-grid
    st = np.array(st)
    uv = np.c_[st[:, 0] * np.cos(st[:, 1]), st[:, 0] * np.sin(st[:, 1])]
    fixed = np.zeros(len(uv), bool)
    # Only the apical region is smoothed: rows beyond SMOOTH_ROWS of the ring
    # are already circles with uniform angles and are held fixed, which also
    # keeps the meridional cell size uniform
    for k in range(int(round(SMOOTH_ROWS * nRing)), nRing + 1):
        fixed[[index[key] for key in rows[k]]] = True
    edges = np.concatenate([quads[:, [0, 1]], quads[:, [1, 2]],
                            quads[:, [2, 3]], quads[:, [3, 0]]])
    edges = np.unique(np.sort(edges, axis=1), axis=0)
    deg = np.bincount(edges.ravel(), minlength=len(uv)).astype(float)
    for _ in range(SMOOTH_ITERS):
        acc = np.zeros_like(uv)
        np.add.at(acc, edges[:, 0], uv[edges[:, 1]])
        np.add.at(acc, edges[:, 1], uv[edges[:, 0]])
        new = acc / deg[:, None]
        uv[~fixed] = uv[~fixed] + 0.8 * (new[~fixed] - uv[~fixed])

    st = np.c_[np.hypot(uv[:, 0], uv[:, 1]), np.arctan2(uv[:, 1], uv[:, 0])]
    return st, quads, baseEdges


HEADER = r"""/*--------------------------------*- C++ -*----------------------------------*\
| solids4foam: solid mechanics and fluid-solid interaction simulations        |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       %s;
    location    "constant/polyMesh";
    object      %s;
}
// Generated by scripts/makeOGridMesh.py
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

"""


def build(level):
    nCap, nRing, nWall = LEVELS[level]
    st, quads, baseEdges = surface_params(nCap, nRing)
    nS = len(st)

    endo = Meridian(R_SHORT_ENDO, R_LONG_ENDO, PHI_BASE_ENDO)
    epi = Meridian(R_SHORT_EPI, R_LONG_EPI, PHI_BASE_EPI)
    pEndo = endo.point(st[:, 0], st[:, 1])
    pEpi = epi.point(st[:, 0], st[:, 1])

    tw = np.linspace(0.0, 1.0, nWall + 1)
    pts = np.concatenate([(1 - t) * pEndo + t * pEpi for t in tw])

    def pid(layer, node):
        return layer * nS + node

    # Hexes: bottom quad on layer l, top quad on layer l+1
    cells = []
    for l in range(nWall):
        for q in quads:
            cells.append([pid(l, v) for v in q] + [pid(l + 1, v) for v in q])
    cells = np.array(cells)

    # Orient every hex right-handed (positive volume)
    # Simple, robust test: signed volume from the bottom-face normal and the
    # wall direction
    p = pts[cells]
    nb = np.cross(p[:, 2] - p[:, 0], p[:, 3] - p[:, 1])
    up = p[:, 4:].mean(1) - p[:, :4].mean(1)
    flip = np.einsum("ij,ij->i", nb, up) < 0
    cells[flip] = cells[flip][:, [0, 3, 2, 1, 4, 7, 6, 5]]

    # Hex faces, outward for a right-handed hex (OpenFOAM's hex model)
    HF = [(0, 4, 7, 3), (1, 2, 6, 5), (0, 1, 5, 4),
          (3, 7, 6, 2), (0, 3, 2, 1), (4, 5, 6, 7)]
    faceMap = {}
    for ci, c in enumerate(cells):
        for f in HF:
            verts = [int(c[k]) for k in f]
            key = tuple(sorted(verts))
            if key in faceMap:
                faceMap[key][1] = ci
            else:
                faceMap[key] = [ci, -1, verts]

    baseKeys = set()
    for (a, b) in baseEdges:
        for l in range(nWall):
            baseKeys.add(tuple(sorted([pid(l, a), pid(l, b),
                                       pid(l + 1, a), pid(l + 1, b)])))

    internal, patches = [], {"base": [], "endo": [], "epi": []}
    for key, (o, nbr, verts) in faceMap.items():
        if nbr >= 0:
            internal.append((o, nbr, verts))
        else:
            if all(v < nS for v in verts):
                patches["endo"].append((o, verts))
            elif all(v >= nWall * nS for v in verts):
                patches["epi"].append((o, verts))
            elif key in baseKeys:
                patches["base"].append((o, verts))
            else:
                raise RuntimeError("Unclassified boundary face")
    internal.sort(key=lambda f: (f[0], f[1]))
    for k in patches:
        patches[k].sort(key=lambda f: f[0])
    return pts, cells, internal, patches


def check(pts, cells, internal, patches):
    """Closedness, positive volumes and non-orthogonality, like checkMesh."""
    nC = len(cells)
    Sf_sum = np.zeros((nC, 3))
    cc = pts[cells].mean(1)
    allf = [(o, n, v) for o, n, v in internal] + \
        [(o, -1, v) for k in patches for o, v in patches[k]]
    vol = np.zeros(nC)
    nonOrth = []
    for o, n, v in allf:
        P = pts[v]
        fc = P.mean(0)
        Sf = 0.5 * np.cross(P[2] - P[0], P[3] - P[1])
        Sf_sum[o] += Sf
        vol[o] += np.dot(Sf, fc - cc[o]) / 3.0
        if n >= 0:
            Sf_sum[n] -= Sf
            vol[n] -= np.dot(Sf, fc - cc[n]) / 3.0
            d = cc[n] - cc[o]
            c = np.dot(d, Sf) / (np.linalg.norm(d) * np.linalg.norm(Sf))
            nonOrth.append(math.degrees(math.acos(max(-1, min(1, c)))))
        else:
            if np.dot(Sf, fc - cc[o]) <= 0:
                raise RuntimeError("Inward boundary face")
    nonOrth = np.array(nonOrth)
    return dict(cells=nC, minVol=vol.min(), maxVol=vol.max(), totVol=vol.sum(),
                maxNonOrth=nonOrth.max(), meanNonOrth=nonOrth.mean(),
                maxOpen=np.abs(Sf_sum).max(), nNeg=int((nonOrth > 90).sum()))


def write(pts, cells, internal, patches, root="constant/polyMesh"):
    os.makedirs(root, exist_ok=True)
    order = ["base", "endo", "epi"]
    faces = [v for _, _, v in internal] + [v for k in order for _, v in patches[k]]
    owner = [o for o, _, _ in internal] + [o for k in order for o, _ in patches[k]]
    neigh = [n for _, n, _ in internal]
    note = f'note "nPoints:{len(pts)} nCells:{len(cells)} nFaces:{len(faces)} nInternalFaces:{len(internal)}";\n'

    with open(f"{root}/points", "w") as fh:
        fh.write(HEADER % ("vectorField", "points"))
        fh.write(f"{len(pts)}\n(\n")
        fh.write("".join(f"({p[0]:.10g} {p[1]:.10g} {p[2]:.10g})\n" for p in pts))
        fh.write(")\n")
    with open(f"{root}/faces", "w") as fh:
        fh.write(HEADER % ("faceList", "faces"))
        fh.write(f"{len(faces)}\n(\n")
        fh.write("".join(f"4({f[0]} {f[1]} {f[2]} {f[3]})\n" for f in faces))
        fh.write(")\n")
    for name, data in (("owner", owner), ("neighbour", neigh)):
        with open(f"{root}/{name}", "w") as fh:
            fh.write((HEADER % ("labelList", name)).replace(
                "FoamFile\n{\n", "FoamFile\n{\n    " + note))
            fh.write(f"{len(data)}\n(\n")
            fh.write("".join(f"{x}\n" for x in data))
            fh.write(")\n")
    with open(f"{root}/boundary", "w") as fh:
        fh.write(HEADER % ("polyBoundaryMesh", "boundary"))
        fh.write(f"{len(order)}\n(\n")
        start = len(internal)
        for k in order:
            fh.write(f"    {k}\n    {{\n        type            patch;\n"
                     f"        nFaces          {len(patches[k])};\n"
                     f"        startFace       {start};\n    }}\n")
            start += len(patches[k])
        fh.write(")\n")


def main():
    level = int(sys.argv[1]) if len(sys.argv) > 1 else 2
    if level not in LEVELS:
        sys.exit(f"Unknown level {level}; choose from {sorted(LEVELS)}")
    pts, cells, internal, patches = build(level)
    q = check(pts, cells, internal, patches)
    print(f"Mesh level {level} (nCap, nRing, nWall) = {LEVELS[level]}: "
          f"{q['cells']} hexahedra, volume {q['totVol']*1e6:.3f} ml")
    print(f"    cell volume min/max {q['minVol']:.3e} / {q['maxVol']:.3e} m^3, "
          f"non-orthogonality max {q['maxNonOrth']:.1f}, mean {q['meanNonOrth']:.1f} deg")
    print(f"    patch faces: " + ", ".join(f"{k} {len(v)}" for k, v in patches.items()))
    if q["minVol"] <= 0:
        sys.exit("Negative or zero cell volume")
    write(pts, cells, internal, patches)


if __name__ == "__main__":
    main()
