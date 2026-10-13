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
Needs only the Python 3 standard library.
"""
import bisect
import math
import os
import sys

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
TRANSITION = 0.45     # fraction of the ring over which square -> circle
SHAPE_POWER = 1.0     # < 1 rounds the rows off the square sooner (less robust)


def linspace(a, b, n):
    return [a + (b - a) * i / (n - 1) for i in range(n)]


def sub(a, b):
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def cross(a, b):
    return (a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0])


def dot(a, b):
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def mean(P):
    n = len(P)
    return (sum(p[0] for p in P) / n, sum(p[1] for p in P) / n,
            sum(p[2] for p in P) / n)


class Meridian:
    """phi(s) for an ellipse meridian parametrised by normalised arc length."""

    def __init__(self, rs, rl, phib, n=20001):
        phi = linspace(0.0, phib, n)
        ds = [math.hypot(rl * math.sin(p), rs * math.cos(p)) for p in phi]
        s = [0.0]
        for i in range(1, n):
            s.append(s[-1] + 0.5 * (ds[i] + ds[i - 1]) * (phi[i] - phi[i - 1]))
        self.s = [x / s[-1] for x in s]
        self.phi = phi
        self.rs, self.rl = rs, rl

    def phi_of(self, s):
        """Linear interpolation of phi(s), as numpy.interp."""
        if s <= self.s[0]:
            return self.phi[0]
        if s >= self.s[-1]:
            return self.phi[-1]
        i = bisect.bisect_right(self.s, s) - 1
        w = (s - self.s[i]) / (self.s[i + 1] - self.s[i])
        return self.phi[i] + w * (self.phi[i + 1] - self.phi[i])

    def point(self, s, theta):
        phi = self.phi_of(s)
        r = self.rs * math.sin(phi)
        return (-self.rl * math.cos(phi), r * math.cos(theta), r * math.sin(theta))


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

    X = linspace(-C_SQ, C_SQ, n + 1)
    for i in range(n + 1):
        for j in range(n + 1):
            add(("c", i, j), X[i], X[j])

    # Square boundary, counter-clockwise starting at (C, -C)
    ring0 = [("c", n, j) for j in range(0, n)] \
        + [("c", i, n) for i in range(n, 0, -1)] \
        + [("c", 0, j) for j in range(n, 0, -1)] \
        + [("c", i, 0) for i in range(0, n)]

    # Angles of the square's boundary nodes, unwrapped to increase
    # continuously, as numpy.unwrap
    thq = []
    for k in ring0:
        th = math.atan2(X[k[2]], X[k[1]])
        if thq:
            while th - thq[-1] > math.pi:
                th -= 2.0 * math.pi
            while th - thq[-1] < -math.pi:
                th += 2.0 * math.pi
        thq.append(th)
    thu = [thq[0] + 2.0 * math.pi * m / nAround for m in range(nAround)]

    rows = [ring0]
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
            d = (math.cos(th), math.sin(th))
            if p > 200:
                norm = max(abs(d[0]), abs(d[1]))
            else:
                norm = (abs(d[0]) ** p + abs(d[1]) ** p) ** (1.0 / p)
            key = ("r", k, m)
            add(key, r * d[0] / norm, r * d[1] / norm)
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
    pEndo = [endo.point(s, th) for s, th in st]
    pEpi = [epi.point(s, th) for s, th in st]

    pts = []
    for t in linspace(0.0, 1.0, nWall + 1):
        pts.extend(tuple((1 - t) * a[i] + t * b[i] for i in range(3))
                   for a, b in zip(pEndo, pEpi))

    def pid(layer, node):
        return layer * nS + node

    # Hexes: bottom quad on layer l, top quad on layer l+1, oriented
    # right-handed (positive volume): the bottom-face normal must point
    # through the wall, towards the top face
    cells = []
    for l in range(nWall):
        for q in quads:
            c = [pid(l, v) for v in q] + [pid(l + 1, v) for v in q]
            p = [pts[v] for v in c]
            nb = cross(sub(p[2], p[0]), sub(p[3], p[1]))
            up = sub(mean(p[4:]), mean(p[:4]))
            if dot(nb, up) < 0:
                c = [c[k] for k in (0, 3, 2, 1, 4, 7, 6, 5)]
            cells.append(c)

    # Hex faces, outward for a right-handed hex (OpenFOAM's hex model)
    HF = [(0, 4, 7, 3), (1, 2, 6, 5), (0, 1, 5, 4),
          (3, 7, 6, 2), (0, 3, 2, 1), (4, 5, 6, 7)]
    faceMap = {}
    for ci, c in enumerate(cells):
        for f in HF:
            verts = [c[k] for k in f]
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
    """Positive volumes, outward boundary faces and non-orthogonality.

    checkMesh, which Allrun runs next, checks the mesh fully; this only stops
    a bad mesh from being written.
    """
    nC = len(cells)
    cc = [mean([pts[v] for v in c]) for c in cells]
    vol = [0.0] * nC
    maxNonOrth = 0.0
    sumNonOrth = 0.0
    allf = [(o, n, v) for o, n, v in internal] + \
        [(o, -1, v) for k in patches for o, v in patches[k]]
    for o, n, v in allf:
        P = [pts[i] for i in v]
        fc = mean(P)
        Sf = tuple(0.5 * x for x in cross(sub(P[2], P[0]), sub(P[3], P[1])))
        vol[o] += dot(Sf, sub(fc, cc[o])) / 3.0
        if n >= 0:
            vol[n] -= dot(Sf, sub(fc, cc[n])) / 3.0
            d = sub(cc[n], cc[o])
            c = dot(d, Sf) / math.sqrt(dot(d, d) * dot(Sf, Sf))
            a = math.degrees(math.acos(max(-1.0, min(1.0, c))))
            maxNonOrth = max(maxNonOrth, a)
            sumNonOrth += a
        elif dot(Sf, sub(fc, cc[o])) <= 0:
            raise RuntimeError("Inward boundary face")
    return dict(cells=nC, minVol=min(vol), maxVol=max(vol), totVol=sum(vol),
                maxNonOrth=maxNonOrth,
                meanNonOrth=sumNonOrth / max(len(internal), 1))


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
