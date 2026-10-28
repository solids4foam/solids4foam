#!/usr/bin/env python3
"""Smoothly graded fluid blockMeshDict for the HronTurek FSI3 replay/coupled cases.

The tutorial's fluid blockMeshDict (24 blocks, cylinder O-grid ring, flag
rows, wake columns) is kept: vertices, arcs, block vertex order and patches
are unchanged, so the flag surface, the patch names and the plate faces are
those of the tutorial. Only the divisions and gradings are replaced.

Design (all lengths in mm at level 1):

* Every block vertex carries a target cell size for each grid-line family
  through it ('x', 'y', 'r' radial, 'c' circumferential on the cylinder);
  the default is isotropic, the walls and the outlet are anisotropic
  (SIZES below, from the few parameters of --help).
* Every block edge gets a distribution whose first and last cell sizes are
  the targets at its two end vertices: two geometric sections whose sizes
  meet at an interior point (blockMesh multi-grading, written with
  edgeGrading so that opposite edges of a block may differ). A section with
  n cells, end sizes h0 and h1 and total expansion e = h1/h0 has length
  n (h1 - h0)/ln(h1/h0) in the limit n -> inf (the logarithmic mean), which
  is the mapping blockMesh refines to when the divisions are doubled with
  the expansion kept. Matching the end sizes of these limit mappings makes
  the cell-size ratio across every block interface tend to 1 with the
  level; at level 1 it is 1 + O(1/n).
* The division count of each topological direction (chain of block edges
  that must share it) is the largest continuum count of its edges, rounded;
  the other edges of the chain absorb the difference in the interior
  sizes (their two sections no longer share one expansion).
* Level L multiplies every division by L and keeps every grading, so the
  levels are one fixed mapping refined uniformly by L in-plane.

    python3 ht_mesh_smooth.py --base system/fluid/blockMeshDict --level 2 --out blockMeshDict
    python3 ht_mesh_smooth.py --base ... --level 1 --report   # chains, sizes, cells
"""
from __future__ import annotations

import argparse
import math
import re
from collections import defaultdict
from pathlib import Path

# Mirror of a vertex about y = 200 mm (the block topology is symmetric)
MIRROR = {0: 27, 1: 28, 2: 29, 3: 30, 64: 69, 4: 31, 5: 22, 6: 23, 7: 24, 8: 25, 65: 68, 9: 26,
          10: 16, 11: 17, 12: 18, 13: 19, 14: 20, 66: 67, 15: 21}
MIRROR.update({v: k for k, v in list(MIRROR.items())})

# Parameters (mm, level 1); see --help
DEFAULTS = dict(
    h_band=5.5,      # ring outer corners (6, 7, 23, 24), band edges y = 120/280
    h_flag=4.5,      # along the flag at the end of the junction blocks (13, 19)
    h_tip=3.8,       # along x at the flag tip (x = 600)
    h_normal=3.4,    # normal to the flag faces and in the wake row y 190-210
    h_wake=8.0,      # along x at x = 1000
    h_out=150.0,     # along x at the outlet
    h_inlet=15.0,    # along x at the inlet
    h_wall=15.0,     # normal to the channel walls
    h_cyl_r=1.2,     # radial first cell on the cylinder (and along the flag root)
    h_cyl_c=2.9,     # circumferential on the cylinder at 10/16 (135 and 225 deg)
    h_root_c=2.0,    # circumferential on the cylinder at 11/17 (45 and 315 deg)
    h_corner=1.6,    # along the cylinder at the flag-cylinder corners 12/18
    scale=1.1,       # multiplies every size; 1.1: 12 568 cells at level 1
    equal_corners=0,  # 1: reshape the ring-corner arcs (CORNER_ARCS), 0: tutorial arcs
)


def vertex_sizes(p: dict) -> dict:
    """Target sizes per vertex (lower half; the upper half is mirrored)."""
    s = p["scale"]
    hb, hf, ht, hn = p["h_band"] * s, p["h_flag"] * s, p["h_tip"] * s, p["h_normal"] * s
    hw, ho, hi, hwall = p["h_wake"] * s, p["h_out"] * s, p["h_inlet"] * s, p["h_wall"] * s
    hr, hc, hrc, hcn = p["h_cyl_r"] * s, p["h_cyl_c"] * s, p["h_root_c"] * s, p["h_corner"] * s
    S = {
        # channel wall y = 0: x as on the line above, y = wall-normal size
        0: {"x": hi, "y": hwall}, 1: {"x": hb, "y": hwall}, 2: {"x": hb, "y": hwall},
        3: {"x": ht, "y": hwall}, 64: {"x": hw, "y": hwall}, 4: {"x": ho, "y": hwall},
        # band edge y = 120
        5: {"x": hi, "y": hb}, 6: hb, 7: hb, 8: {"x": ht, "y": hb}, 65: {"x": hw, "y": hb},
        9: {"x": ho, "y": hb},
        # cylinder
        10: {"c": hc, "r": hr}, 11: {"c": hrc, "r": hr}, 12: {"c": hcn, "x": hr},
        # flag lower face and the wake row
        13: {"x": hf, "y": hn}, 14: {"x": ht, "y": hn}, 66: {"x": hw, "y": hn}, 15: {"x": ho, "y": hn},
    }
    for v in list(S):
        S[MIRROR[v]] = S[v]
    return S


def edge_class(a: int, b: int, pts: dict, arcs: dict) -> str:
    cyl = {10, 11, 12, 16, 17, 18}
    if (a, b) in arcs or (b, a) in arcs:
        if a in cyl and b in cyl:
            return "c"
        pa, pb = pts[a], pts[b]
        return "x" if abs(pb[0] - pa[0]) > abs(pb[1] - pa[1]) else "y"
    pa, pb = pts[a], pts[b]
    dx, dy = abs(pb[0] - pa[0]), abs(pb[1] - pa[1])
    if (a in cyl) != (b in cyl) and dx > 1e-6 and dy > 1e-6:
        return "r"
    return "x" if dx > dy else "y"


def size_at(S: dict, v: int, cls: str) -> float:
    s = S[v]
    if not isinstance(s, dict):
        return s
    if cls in s:
        return s[cls]
    # radial edges at the ring outer corners and the flag face edges of the
    # junction blocks continue the 'x' lines
    if cls == "r" and "x" in s:
        return s["x"]
    raise KeyError(f"no size for vertex {v} class {cls}")


def logmean(a: float, b: float) -> float:
    return a if abs(a - b) < 1e-12 * a else (b - a) / math.log(b / a)


def arc_length(p0, pm, p1) -> float:
    """Length of the circular arc through p0, pm, p1 (2-D)."""
    ax, ay = p0
    bx, by = pm
    cx, cy = p1
    d = 2 * (ax * (by - cy) + bx * (cy - ay) + cx * (ay - by))
    ux = ((ax**2 + ay**2) * (by - cy) + (bx**2 + by**2) * (cy - ay) + (cx**2 + cy**2) * (ay - by)) / d
    uy = ((ax**2 + ay**2) * (cx - bx) + (bx**2 + by**2) * (ax - cx) + (cx**2 + cy**2) * (bx - ax)) / d
    r = math.hypot(ax - ux, ay - uy)
    t0, tm, t1 = (math.atan2(q[1] - uy, q[0] - ux) for q in (p0, pm, p1))

    def ccw(a, b):
        return (b - a) % (2 * math.pi)
    span = ccw(t0, t1)
    if ccw(t0, tm) > span:  # the mid point is on the other way round
        span = 2 * math.pi - span
    return r * span


class Dist:
    """Two-section distribution of n cells on an edge a->b, end sizes ha, hb."""

    def __init__(self, length: float, n: int, ha: float, hb: float):
        self.L, self.n, self.ha, self.hb = length, n, ha, hb
        n1 = n // 2
        n2 = n - n1
        self.n1, self.n2 = n1, n2

        def f(hm):
            return n1 * logmean(ha, hm) + n2 * logmean(hm, hb) - length
        lo, hi = 1e-6, 1e6
        for _ in range(200):
            mid = math.sqrt(lo * hi)
            if f(mid) > 0:
                hi = mid
            else:
                lo = mid
        self.hm = math.sqrt(lo * hi)
        self.L1 = n1 * logmean(ha, self.hm)
        self.e1, self.e2 = self.hm / ha, hb / self.hm

    def per_cell(self) -> float:
        """Largest cell-to-cell expansion (limit mapping, this level)."""
        q1 = self.e1 ** (1 / max(self.n1, 1))
        q2 = self.e2 ** (1 / max(self.n2, 1))
        return max(q1, 1 / q1, q2, 1 / q2)

    def grading(self, reverse: bool = False) -> str:
        secs = [(self.L1 / self.L, self.n1 / self.n, self.e1), (1 - self.L1 / self.L, self.n2 / self.n, self.e2)]
        if reverse:
            secs = [(fl, fn, 1 / e) for fl, fn, e in reversed(secs)]
        return "(" + " ".join(f"({fl:.12g} {fn:.12g} {e:.12g})" for fl, fn, e in secs) + ")"


def arc_mid(p0, p1, t0_deg: float):
    """Mid point (on the arc) of the circular arc p0 -> p1 whose tangent at p0 has direction t0_deg."""
    t = (math.cos(math.radians(t0_deg)), math.sin(math.radians(t0_deg)))
    d = (p1[0] - p0[0], p1[1] - p0[1])
    n = (-t[1], t[0])  # left normal
    nd = n[0] * d[0] + n[1] * d[1]
    r = (d[0] ** 2 + d[1] ** 2) / (2 * nd)  # signed radius, centre at p0 + r n
    c = (p0[0] + r * n[0], p0[1] + r * n[1])
    m = (0.5 * (p0[0] + p1[0]), 0.5 * (p0[1] + p1[1]))
    v = (m[0] - c[0], m[1] - c[1])
    lv = math.hypot(*v)
    return (c[0] + abs(r) * v[0] / lv, c[1] + abs(r) * v[1] / lv)


# Option (not the default): reshape the ring outer arcs and the junction arcs
# so that the four non-Cartesian blocks at each valence-5 ring corner (6, 7,
# 23, 24) meet at 67.5 deg (tutorial arcs: 61.9 to 73.1 deg next to the 90 deg
# of the Cartesian block). Tested: it lowers the corner jump (1.13 -> 1.08)
# but flattens the arcs, so the transfinite blending of the ring and of the
# outer blocks no longer matches along the arcs (1.21 at the arc mid points,
# level-independent); the tutorial arcs balance the two and are kept.
CORNER_ARCS = {  # (a, b): tangent direction at a (deg)
    (6, 7): -22.5, (6, 23): 112.5, (7, 13): 67.5,
    (23, 24): 22.5, (24, 19): -67.5,
}


def equalise_corners(base: str) -> str:
    pts, arcs, _ = parse(base)
    out = base
    for (a, b), t in CORNER_ARCS.items():
        m = arc_mid(pts[a], pts[b], t)
        for za, zb, z in ((a, b, 0), (a + 32, b + 32, 15)):
            pat = re.compile(rf"arc\s+({za}\s+{zb}|{zb}\s+{za})\s*\([^)]*\)")
            new = f"arc {za} {zb} ({m[0]:.6f} {m[1]:.6f} {z})"
            out, k = pat.subn(new, out)
            if k != 1:
                raise SystemExit(f"arc {za} {zb} not found")
    return out


def parse(base: str):
    vtxt = base[base.index("vertices"):base.index("blocks")]
    pts = [tuple(float(x) for x in m.groups()) for m in
           re.finditer(r"\(\s*([-\d.eE+]+)\s+([-\d.eE+]+)\s+([-\d.eE+]+)\s*\)", vtxt)]
    etxt = base[base.index("edges"):base.index("patches")]
    arcs = {}
    for m in re.finditer(r"arc\s+(\d+)\s+(\d+)\s*\(\s*([-\d.eE+]+)\s+([-\d.eE+]+)\s+([-\d.eE+]+)\s*\)", etxt):
        a, b = int(m.group(1)), int(m.group(2))
        arcs[(a, b)] = (float(m.group(3)), float(m.group(4)))
    btxt = base[base.index("blocks"):base.index("edges")]
    blocks = [list(map(int, m.group(1).split())) for m in re.finditer(r"hex\s+\(([\d ]+)\)", btxt)]
    if len(blocks) != 24:
        raise SystemExit(f"expected 24 blocks, found {len(blocks)}")
    return {i: (p[0], p[1]) for i, p in enumerate(pts)}, arcs, blocks


def design(base: str, params: dict):
    if params.get("equal_corners", 0):
        base = equalise_corners(base)
    pts, arcs, blocks = parse(base)
    S = vertex_sizes(params)

    def length(a, b):
        if (a, b) in arcs:
            return arc_length(pts[a], arcs[(a, b)], pts[b])
        if (b, a) in arcs:
            return arc_length(pts[b], arcs[(b, a)], pts[a])
        return math.dist(pts[a], pts[b])

    # chains of edges that share a division count (union-find on z = 0 edges)
    parent = {}

    def key(a, b):
        return (min(a, b), max(a, b))

    def find(e):
        parent.setdefault(e, e)
        while parent[e] != e:
            parent[e] = parent[parent[e]]
            e = parent[e]
        return e

    def union(e, f):
        parent[find(e)] = find(f)
    for v in blocks:
        union(key(v[0], v[1]), key(v[3], v[2]))
        union(key(v[0], v[3]), key(v[1], v[2]))
    chains = defaultdict(list)
    for e in list(parent):
        chains[find(e)].append(e)
    # division count per chain: largest continuum count of its edges
    info = {}
    nchain = {}
    for c, edges in chains.items():
        best, best_c = 0.0, 0.0
        for a, b in edges:
            cls = edge_class(a, b, pts, arcs)
            L = length(a, b)
            ha, hb = size_at(S, a, cls), size_at(S, b, cls)
            ncont = L / logmean(ha, hb)
            info[(a, b)] = (cls, L, ha, hb, ncont)
            if cls == "c":  # cylinder arcs follow the ring outer edges
                best_c = max(best_c, ncont)
            else:
                best = max(best, ncont)
        nchain[c] = max(1, round(best or best_c))
    dists = {}
    for c, edges in chains.items():
        for a, b in edges:
            cls, L, ha, hb, _ = info[(a, b)]
            dists[(a, b)] = Dist(L, nchain[c], ha, hb)
    return pts, arcs, blocks, chains, nchain, dists, find, key


def build(base: str, level: int, params: dict | None = None) -> str:
    p = dict(DEFAULTS, **(params or {}))
    pts, arcs, blocks, chains, nchain, dists, find, key = design(base, p)

    def g(a, b):  # a, b: vertices of the z = 0 face (the first four of each hex)
        k = key(a, b)
        return dists[k].grading(reverse=(k != (a, b)))

    lines = []
    for v in blocks:
        nx = nchain[find(key(v[0], v[1]))] * level
        ny = nchain[find(key(v[0], v[3]))] * level
        eg = [g(v[0], v[1]), g(v[3], v[2]), g(v[3], v[2]), g(v[0], v[1]),
              g(v[0], v[3]), g(v[1], v[2]), g(v[1], v[2]), g(v[0], v[3]),
              "1", "1", "1", "1"]
        lines.append(f"    hex ({' '.join(map(str, v))}) ({nx} {ny} 1)\n        edgeGrading\n        (\n"
                     + "".join(f"            {x}\n" for x in eg) + "        )")
    if p.get("equal_corners", 0):
        base = equalise_corners(base)
    i0 = base.index("blocks")
    i1 = base.index("edges", i0)
    head = base[:i0]
    tail = base[i1:]
    return head + "blocks\n(\n" + "\n".join(lines) + "\n);\n\n" + tail


def report(base: str, params: dict) -> None:
    p = dict(DEFAULTS, **params)
    pts, arcs, blocks, chains, nchain, dists, find, key = design(base, p)
    cells = 0
    for v in blocks:
        cells += nchain[find(key(v[0], v[1]))] * nchain[find(key(v[0], v[3]))]
    print(f"cells at level 1: {cells}  (x4 per level: {cells * 4}, {cells * 9}, {cells * 16})")
    for c, edges in sorted(chains.items(), key=lambda kv: kv[0]):
        print(f"chain {c}: n = {nchain[c]}")
        for e in sorted(edges):
            d = dists[e]
            print(f"   {e[0]:2d}-{e[1]:2d}  L {d.L:7.2f}  h {d.ha:6.2f} -> {d.hm:6.2f} -> {d.hb:6.2f}  "
                  f"max cell ratio {d.per_cell():.3f}")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--base", type=Path, required=True, help="tutorial fluid blockMeshDict")
    ap.add_argument("--level", type=int, default=1)
    ap.add_argument("--out", type=Path)
    ap.add_argument("--report", action="store_true")
    for k, v in DEFAULTS.items():
        ap.add_argument("--" + k.replace("_", "-"), type=float, default=v)
    a = ap.parse_args()
    params = {k: getattr(a, k) for k in DEFAULTS}
    base = a.base.read_text()
    if a.report:
        report(base, params)
    if a.out:
        a.out.write_text(build(base, a.level, params))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
