#!/usr/bin/env python3
"""Neighbour cell-size audit of a 2-D (one cell thick) OpenFOAM mesh.

Reads constant/polyMesh (ascii) of each case and reports, over the internal
faces, the neighbour cell-area ratio r = max(A_o/A_n, A_n/A_o) (the measure of
~/.cache/ht_meshaudit/ratio.py used in the coordinator's audit): its maximum
and location, the number of faces above 1.2, 1.5 and 1.9, the 99.9th
percentile, the maximum in the region of interest (x 0.1-1.0 m, y 0.08-0.32
m: cylinder ring, flag and near wake) and the lines along which faces with
r > --line-threshold cluster (block-interface jumps).

    python3 ht_mesh_audit.py <case> [<case> ...] [--json audit.json]
"""
from __future__ import annotations

import argparse
import json
import re
from collections import defaultdict
from pathlib import Path

import numpy as np


def _body(path: Path) -> str:
    t = path.read_text()
    t = re.sub(r"/\*.*?\*/", "", t, flags=re.S)
    i = t.index("FoamFile")
    return t[t.index("}", i) + 1:]


def _count_open(t: str):
    m = re.search(r"(\d+)\s*\(", t)
    return int(m.group(1)), m.end()


def read_points(path: Path) -> np.ndarray:
    t = _body(path)
    n, e = _count_open(t)
    v = np.array(re.findall(r"\(([-\d.eE+]+) ([-\d.eE+]+) ([-\d.eE+]+)\)", t[e - 1:])[:n], float)
    return v


def read_faces(path: Path) -> list[list[int]]:
    t = _body(path)
    n, e = _count_open(t)
    return [list(map(int, f.split())) for f in re.findall(r"\d+\(([\d ]+)\)", t[e:])][:n]


def read_labels(path: Path) -> np.ndarray:
    t = _body(path)
    n, e = _count_open(t)
    return np.array(t[e:].split(")")[0].split(), int)[:n]


def audit(case: Path, threshold: float = 1.5) -> dict:
    d = case / "constant/polyMesh"
    P, F = read_points(d / "points"), read_faces(d / "faces")
    O, N = read_labels(d / "owner"), read_labels(d / "neighbour")
    nc = int(O.max()) + 1
    A = np.zeros(nc)
    C = np.zeros((nc, 2))
    zmax = P[:, 2].max()
    for fi, f in enumerate(F):
        q = P[f]
        if np.allclose(q[:, 2], zmax):
            x, y = q[:, 0], q[:, 1]
            A[O[fi]] = 0.5 * abs(np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1)))
            C[O[fi]] = q[:, :2].mean(0)
    # area centroids of the cells (front polygons)
    for fi, f in enumerate(F):
        q = P[f]
        if np.allclose(q[:, 2], zmax):
            x, y = q[:, 0], q[:, 1]
            cr = x * np.roll(y, -1) - np.roll(x, -1) * y
            a6 = 3.0 * cr.sum()
            C[O[fi]] = [((x + np.roll(x, -1)) * cr).sum() / a6, ((y + np.roll(y, -1)) * cr).sum() / a6]
    ni = len(N)
    r = np.maximum(A[O[:ni]] / A[N], A[N] / A[O[:ni]])
    # non-orthogonality of the in-plane internal faces: angle between the face
    # normal and the owner-neighbour centroid vector
    nonorth = np.zeros(ni)
    for fi in range(ni):
        q = P[F[fi]]
        if np.ptp(q[:, 2]) > 0:
            e = q[np.argsort(q[:, 2])][:2, :2]  # the two points of the z = min edge
            t = e[1] - e[0]
            n = np.array([t[1], -t[0]])
            d = C[N[fi]] - C[O[fi]]
            c = abs(n @ d) / (np.linalg.norm(n) * np.linalg.norm(d))
            nonorth[fi] = np.degrees(np.arccos(min(1.0, c)))
    mid = 0.5 * (C[O[:ni]] + C[N])
    k = int(np.argmax(r))
    roi = (mid[:, 0] > 0.1) & (mid[:, 0] < 1.0) & (mid[:, 1] > 0.08) & (mid[:, 1] < 0.32)
    out = {
        "case": str(case), "cells": nc, "internal_faces": ni,
        "max_ratio": float(r.max()), "max_at": [float(mid[k, 0]), float(mid[k, 1])],
        "p999_ratio": float(np.quantile(r, 0.999)),
        "max_ratio_roi": float(r[roi].max()),
        "n_gt_1.2": int((r > 1.2).sum()), "n_gt_1.5": int((r > 1.5).sum()), "n_gt_1.9": int((r > 1.9).sum()),
        "min_area": float(A.min()), "max_area": float(A.max()),
        "max_nonorth": float(nonorth.max()),
        "max_nonorth_at": [float(mid[int(np.argmax(nonorth)), 0]), float(mid[int(np.argmax(nonorth)), 1])],
    }
    sel = np.where(r > threshold)[0]
    lines = defaultdict(list)
    for i in sel:
        dvec = C[N[i]] - C[O[i]]
        if abs(dvec[1]) > abs(dvec[0]):
            lines[("y", round(float(mid[i, 1]), 3))].append(float(mid[i, 0]))
        else:
            lines[("x", round(float(mid[i, 0]), 3))].append(float(mid[i, 1]))
    out["jump_lines"] = [{"line": f"{a}={b:.3f}", "n": len(v), "span": [min(v), max(v)]}
                         for (a, b), v in sorted(lines.items(), key=lambda kv: -len(kv[1]))[:10]]
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("cases", nargs="+", type=Path)
    ap.add_argument("--line-threshold", type=float, default=1.5)
    ap.add_argument("--json", type=Path)
    a = ap.parse_args()
    res = []
    for c in a.cases:
        r = audit(c, a.line_threshold)
        res.append(r)
        print(f"{c}: cells {r['cells']}  max ratio {r['max_ratio']:.3f} at ({r['max_at'][0]:.3f}, {r['max_at'][1]:.3f})"
              f"  roi max {r['max_ratio_roi']:.3f}  p99.9 {r['p999_ratio']:.3f}"
              f"  faces >1.2 {r['n_gt_1.2']}  >1.5 {r['n_gt_1.5']}  >1.9 {r['n_gt_1.9']}"
              f"  non-orth max {r['max_nonorth']:.1f} at ({r['max_nonorth_at'][0]:.3f}, {r['max_nonorth_at'][1]:.3f})")
        for ln in r["jump_lines"]:
            print(f"   line {ln['line']}  n={ln['n']}  span {ln['span'][0]:.3f}..{ln['span'][1]:.3f}")
    if a.json:
        a.json.write_text(json.dumps(res, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
