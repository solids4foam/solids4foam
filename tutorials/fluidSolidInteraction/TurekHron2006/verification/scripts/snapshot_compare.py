#!/usr/bin/env python3
"""Where do two levels of a replay family differ? (common-grid snapshot comparison)

Reads reconstructed ascii U and p at one time (same motion phase on every
level) on the moving meshes (<time>/polyMesh/points), carries each coarse
field to the cell centroids of the fine mesh by a linear reconstruction from
the nearest coarse cell (value + leastSquares gradient . offset), and splits
the squared difference integrals int |dU|^2 dA and int dp^2 dA over regions:

  corners   within 0.03 m of the four valence-5 ring corners (0.12/0.28, 0.12/0.28)
  ring      within 0.115 m of the cylinder centre, not corners
  flag      x 0.25-0.6 m, y 0.12-0.28 m (flag rows), not ring/corners
  nearwake  x 0.6-1.0 m, y 0.1-0.3 m
  farwake   x > 1.0 m, y 0.1-0.3 m
  outer     the rest
  (and, overlapping, "flagwall": within 0.01 m of the flag surface)

The pressure difference is taken after removing its area-weighted mean over
the near field (x < 1 m, 0.08 < y < 0.32 m): a uniform pressure offset does
not change the flag loads.

    python3 snapshot_compare.py --time 30.2 --fine sm_4x --coarse sm_2x sm_3x --png diff.png
"""
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import numpy as np

import ht_mesh_audit as au

CORNERS = np.array([[0.12, 0.12], [0.28, 0.12], [0.12, 0.28], [0.28, 0.28]])


def read_field(path: Path, nc: int) -> np.ndarray:
    t = path.read_text()
    m = re.search(r"internalField\s+nonuniform\s+List<(\w+)>\s*(\d+)\s*\(", t)
    n = int(m.group(2))
    body = t[m.end():]
    if m.group(1) == "scalar":
        v = np.array(body.split(")", 1)[0].split()[:n], float)
    else:
        v = np.array(re.findall(r"\(([-\d.eE+]+) ([-\d.eE+]+) ([-\d.eE+]+)\)", body[: body.index(";")])[:n], float)
    assert len(v) == nc, (path, len(v), nc)
    return v


def mesh(case: Path, time: str):
    d = case / "constant/polyMesh"
    P = au.read_points(case / time / "polyMesh/points")
    F = au.read_faces(d / "faces")
    O, N = au.read_labels(d / "owner"), au.read_labels(d / "neighbour")
    nc = int(O.max()) + 1
    C = np.zeros((nc, 2))
    A = np.zeros(nc)
    zmax = P[:, 2].max()
    polys = [None] * nc
    for fi, f in enumerate(F):
        q = P[f]
        if np.allclose(q[:, 2], zmax):
            x, y = q[:, 0], q[:, 1]
            cr = x * np.roll(y, -1) - np.roll(x, -1) * y
            a6 = 3.0 * cr.sum()
            C[O[fi]] = [((x + np.roll(x, -1)) * cr).sum() / a6, ((y + np.roll(y, -1)) * cr).sum() / a6]
            A[O[fi]] = abs(a6) / 6.0
            polys[O[fi]] = q[:, :2]
    return C, A, O[: len(N)], N, polys


def lsq_grad(C, O, N, phi):
    nc = len(C)
    d = C[N] - C[O]
    w = 1.0 / np.einsum("ij,ij->i", d, d)
    G = np.zeros((nc, 2, 2))
    b = np.zeros((nc, 2) + phi.shape[1:])
    dd = np.einsum("i,ij,ik->ijk", w, d, d)
    np.add.at(G, O, dd)
    np.add.at(G, N, dd)
    dphi = phi[N] - phi[O]
    wd = w[:, None] * d
    contrib = np.einsum("ij,i...->ij...", wd, dphi)
    np.add.at(b, O, contrib)
    np.add.at(b, N, contrib)
    G += 1e-30 * np.eye(2)
    Ginv = np.linalg.inv(G)
    return np.einsum("ijk,ik...->ij...", Ginv, b)


def nearest(src, tgt, cell=0.004):
    """Index of the nearest src point for every tgt point (uniform-bin search)."""
    key = np.floor(src / cell).astype(np.int64)
    order = np.lexsort((key[:, 1], key[:, 0]))
    bins = {}
    for i in order:
        bins.setdefault((key[i, 0], key[i, 1]), []).append(i)
    bins = {k: np.array(v) for k, v in bins.items()}
    tk = np.floor(tgt / cell).astype(np.int64)
    out = np.empty(len(tgt), np.int64)
    for j, (p, k) in enumerate(zip(tgt, tk)):
        best, bi, r, stop = np.inf, -1, 0, None
        while stop is None or r <= stop:
            ring = [(k[0] + dx, k[1] + dy) for dx in range(-r, r + 1) for dy in range(-r, r + 1)
                    if max(abs(dx), abs(dy)) == r]
            cand = [bins[q] for q in ring if q in bins]
            if cand:
                c = np.concatenate(cand)
                dd = np.sum((src[c] - p) ** 2, axis=1)
                m = int(np.argmin(dd))
                if dd[m] < best:
                    best, bi = dd[m], int(c[m])
            if bi >= 0 and stop is None:
                # a closer point can lie at most sqrt(2) times further out in bins
                stop = int(np.ceil((r + 1) * 1.415))
            r += 1
        out[j] = bi
    return out


def region(C):
    x, y = C[:, 0], C[:, 1]
    corner = np.min(np.linalg.norm(C[:, None, :] - CORNERS[None], axis=2), axis=1) < 0.03
    ring = (np.hypot(x - 0.2, y - 0.2) < 0.115) & ~corner
    flag = (x > 0.25) & (x < 0.6) & (y > 0.12) & (y < 0.28) & ~ring & ~corner
    near = (x >= 0.6) & (x < 1.0) & (y > 0.1) & (y < 0.3)
    far = (x >= 1.0) & (y > 0.1) & (y < 0.3)
    outer = ~(corner | ring | flag | near | far)
    fx = np.clip(x, 0.249, 0.6)
    fy = np.clip(y, 0.19, 0.21)
    wall = (np.hypot(x - fx, y - fy) < 0.01) & (x > 0.24)
    return {"corners": corner, "ring": ring, "flag": flag, "nearwake": near, "farwake": far, "outer": outer,
            "flagwall": wall}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--time", default="30.2")
    ap.add_argument("--fine", type=Path, required=True)
    ap.add_argument("--coarse", type=Path, nargs="+", required=True)
    ap.add_argument("--json", type=Path)
    ap.add_argument("--png", type=Path)
    a = ap.parse_args()
    Cf, Af, Of, Nf, polys = mesh(a.fine, a.time)
    Uf = read_field(a.fine / a.time / "U", len(Cf))[:, :2]
    pf = read_field(a.fine / a.time / "p", len(Cf))
    regs = region(Cf)
    res, maps = {}, {}
    for cc in a.coarse:
        Cc, Ac, Oc, Nc, _ = mesh(cc, a.time)
        Uc = read_field(cc / a.time / "U", len(Cc))[:, :2]
        pc = read_field(cc / a.time / "p", len(Cc))
        gU = lsq_grad(Cc, Oc, Nc, Uc)          # (nc, 2, 2): d/dx_j of U_k -> [cell, j, k]
        gp = lsq_grad(Cc, Oc, Nc, pc)          # (nc, 2)
        idx = nearest(Cc, Cf)
        off = Cf - Cc[idx]
        Ui = Uc[idx] + np.einsum("ij,ijk->ik", off, gU[idx])
        pi = pc[idx] + np.einsum("ij,ij->i", off, gp[idx])
        dU2 = np.sum((Uf - Ui) ** 2, axis=1)
        dp = pf - pi
        nf = (Cf[:, 0] < 1.0) & (Cf[:, 1] > 0.08) & (Cf[:, 1] < 0.32)
        offset = float((dp * Af)[nf].sum() / Af[nf].sum())
        dp2 = (dp - offset) ** 2
        tot_U, tot_p = float((dU2 * Af).sum()), float((dp2 * Af).sum())
        r = {"L2_dU": tot_U ** 0.5, "L2_dp": tot_p ** 0.5, "p_offset": offset, "regions": {}}
        for name, m in regs.items():
            r["regions"][name] = {"area": float(Af[m].sum()), "share_dU2": float((dU2 * Af)[m].sum() / tot_U),
                                  "share_dp2": float((dp2 * Af)[m].sum() / tot_p),
                                  "rms_dU": float(np.sqrt((dU2 * Af)[m].sum() / Af[m].sum())),
                                  "rms_dp": float(np.sqrt((dp2 * Af)[m].sum() / Af[m].sum()))}
        res[f"{a.fine.name}-{cc.name}"] = r
        maps[cc.name] = np.sqrt(dU2)
        print(f"{a.fine.name} - {cc.name}: L2 |dU| {r['L2_dU']:.4e} m^2/s, L2 dp {r['L2_dp']:.4e} "
              f"(near-field mean dp {offset:+.4f} removed)")
        for name, q in r["regions"].items():
            print(f"   {name:9s} area {q['area']:.4f}  share dU^2 {q['share_dU2']:6.3f}  dp^2 {q['share_dp2']:6.3f}"
                  f"  rms |dU| {q['rms_dU']:.4f}  rms dp {q['rms_dp']:.4f}")
    if a.json:
        a.json.write_text(json.dumps(res, indent=2))
    if a.png:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.collections import PolyCollection
        fig, ax = plt.subplots(len(maps), 1, figsize=(11, 3.2 * len(maps)), squeeze=False)
        for axi, (name, v) in zip(ax[:, 0], maps.items()):
            pcol = PolyCollection(polys, array=v, cmap="magma", edgecolor="face", linewidth=0.1, antialiased=False)
            pcol.set_clim(0, np.quantile(v, 0.995))
            axi.add_collection(pcol)
            axi.set_xlim(0, 1.6)
            axi.set_ylim(0, 0.41)
            axi.set_aspect("equal")
            axi.plot(CORNERS[:, 0], CORNERS[:, 1], "c+", ms=10)
            axi.set_title(f"|U({a.fine.name}) - U({name})| at t = {a.time} s (m/s); + ring corners")
            fig.colorbar(pcol, ax=axi, fraction=0.02)
        fig.tight_layout()
        fig.savefig(a.png, dpi=130)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
