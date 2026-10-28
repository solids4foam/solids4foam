#!/usr/bin/env python3
"""Instantaneous flag pressure load per region at one time, compared across levels.

For each case (reconstructed ascii p at --time, moving mesh points in
<time>/polyMesh/points), the plate face pressure is the fixedGradient value
of movingWallPressure, p_f = p_c + g (n . (C_f - C_c)), and the pressure
force rho p_f S_f/depth (N/m) is summed over the regions of
traction_regions.py (x bins on the upper and lower faces of the undeformed
flag, face x-extent overlap weights, plus the end face). Printed: the y force
per region and level, and the change between successive levels per region
with its share of the total change. One instant only (the same motion phase
on every level): a localisation aid, not a harmonic decomposition.

    python3 snapshot_plate.py --time 30.2 sm_2x sm_3x sm_4x
"""
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import numpy as np

import ht_mesh_audit as au
from traction_regions import EDGES, regions  # noqa: E402

RHO, DEPTH = 1000.0, 0.015


def patch(boundary: Path, name: str):
    t = boundary.read_text()
    m = re.search(rf"\n\s*{name}\s*\{{[^}}]*nFaces\s+(\d+);[^}}]*startFace\s+(\d+);", t)
    return int(m.group(2)), int(m.group(1))


def plate_load(case: Path, time: str):
    d = case / "constant/polyMesh"
    P0 = au.read_points(d / "points")
    P = au.read_points(case / time / "polyMesh/points")
    F = au.read_faces(d / "faces")
    O = au.read_labels(d / "owner")
    start, n = patch(d / "boundary", "plate")
    nc = int(O.max()) + 1
    # cell centres (front polygon centroids) of the plate owner cells only
    zmax = P[:, 2].max()
    owners = set(O[start:start + n].tolist())
    Cc = {}
    for fi, f in enumerate(F):
        if O[fi] in owners and fi >= 0:
            q = P[f]
            if np.allclose(q[:, 2], zmax):
                x, y = q[:, 0], q[:, 1]
                cr = x * np.roll(y, -1) - np.roll(x, -1) * y
                a6 = 3.0 * cr.sum()
                Cc[O[fi]] = np.array([((x + np.roll(x, -1)) * cr).sum() / a6, ((y + np.roll(y, -1)) * cr).sum() / a6])
    t = (case / time / "p").read_text()
    m = re.search(r"internalField\s+nonuniform\s+List<scalar>\s*(\d+)\s*\(", t)
    pc = np.array(t[m.end():].split(")", 1)[0].split()[: int(m.group(1))], float)
    assert len(pc) == nc
    b = t[t.index("\n    plate"):]
    m = re.search(r"gradient\s+nonuniform\s+List<scalar>\s*(\d+)\s*\(", b)
    g = np.array(b[m.end():].split(")", 1)[0].split()[: int(m.group(1))], float)
    rows = []
    for i in range(n):
        f = F[start + i]
        q, q0 = P[f], P0[f]
        cf = q[:, :2].mean(0)
        # face area vector (2-D: edge of length L times depth), outward from the fluid cell
        e = q[np.argsort(q[:, 2])][:2, :2]
        tvec = e[1] - e[0]
        nvec = np.array([tvec[1], -tvec[0]])
        if np.dot(nvec, cf - Cc[O[start + i]]) < 0:
            nvec = -nvec
        nh = nvec / np.linalg.norm(nvec)
        pf = pc[O[start + i]] + g[i] * np.dot(nh, cf - Cc[O[start + i]])
        force = RHO * pf * nvec * DEPTH / DEPTH  # |nvec| = edge length; force per unit depth
        rows.append((q0[:, 0].min(), q0[:, 0].max(), q0[:, 1].min(), q0[:, 1].max(), force[0], force[1]))
    return np.array(rows)


def sums(rows):
    x0, x1, y0, y1, fx, fy = rows.T
    end = (x1 - x0) < 1e-9
    kind = np.where(end, "end", np.where(0.5 * (y0 + y1) > 0.2, "upper", "lower"))
    out = {}
    for name, k, iv in regions():
        if k == "end":
            w = (kind == "end").astype(float)
        else:
            ov = np.clip(np.minimum(x1, iv[1]) - np.maximum(x0, iv[0]), 0, None)
            w = np.where(kind == k, ov / np.maximum(x1 - x0, 1e-12), 0.0)
        out[name] = {"fy": float(np.sum(w * fy)), "fx": float(np.sum(w * fx))}
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("cases", nargs="+", type=Path)
    ap.add_argument("--time", default="30.2")
    ap.add_argument("--json", type=Path)
    a = ap.parse_args()
    data = {c.name: sums(plate_load(c, a.time)) for c in a.cases}
    names = [r[0] for r in regions()]
    lv = [c.name for c in a.cases]
    tot = {l: sum(data[l][n]["fy"] for n in names) for l in lv}
    print("pressure lift (N/m) at t =", a.time, " ".join(f"{l} {tot[l]:.3f}" for l in lv))
    pairs = list(zip(lv[:-1], lv[1:]))
    print("region".ljust(18) + "".join(f"{l:>10s}" for l in lv) + "".join(f"   d {p[0]}->{p[1]} (share)" for p in pairs))
    out = {"time": a.time, "totals": tot, "regions": {}}
    for nme in names:
        v = [data[l][nme]["fy"] for l in lv]
        s = nme.ljust(18) + "".join(f"{x:10.3f}" for x in v)
        for p in pairs:
            d = data[p[1]][nme]["fy"] - data[p[0]][nme]["fy"]
            s += f"   {d:8.3f} ({d / (tot[p[1]] - tot[p[0]]):6.2f})"
        out["regions"][nme] = v
        print(s)
    # net per x bin (upper + lower: a uniform pressure shift cancels)
    print("net per x bin (upper + lower), change and share of the total change")
    out["net"] = {}
    for lo, hi in zip(EDGES[:-1], EDGES[1:]):
        k = f"{lo:.2f}-{hi:.2f}"
        v = [data[l][f"upper {k}"]["fy"] + data[l][f"lower {k}"]["fy"] for l in lv]
        out["net"][k] = v
        print(f"  x {k}" + "".join(f"{x:10.3f}" for x in v) + "".join(
            f"   {v[i + 1] - v[i]:8.3f} ({(v[i + 1] - v[i]) / (tot[lv[i + 1]] - tot[lv[i]]):6.2f})" for i in range(len(lv) - 1)))
    if a.json:
        a.json.write_text(json.dumps(out, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
