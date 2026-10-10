#!/usr/bin/env python3
"""Compare the #546 A/B restarts (old vs fixed library) and the original run.

FSI3 2x: flag interface point positions at every step (plateMotion.dat of each
processor), total drag/lift (fluid/forces) and the point-A displacement
(solidPointDisplacement). CFD3 2x: total drag/lift. Prints maximum absolute
and relative differences over the common time steps.
"""
import sys
from pathlib import Path

import numpy as np


def plate_motion(case: Path):
    out = {}
    for proc in sorted(case.glob("processor*")):
        f = proc / "postProcessing/plateMotion.dat"
        if not f.is_file():
            continue
        lines = f.read_text().splitlines()
        data = np.array([[float(v) for v in line.split()] for line in lines[1:] if line.strip()])
        out[proc.name] = data
    return out


def forces(case: Path, sub: str):
    fs = sorted(case.glob(f"postProcessing/{sub}/*/force.dat"))
    rows = []
    for f in fs:
        for line in f.read_text().splitlines():
            if line.startswith("#") or not line.strip():
                continue
            v = line.replace("(", " ").replace(")", " ").split()
            rows.append([float(x) for x in v[:3]])
    a = np.array(rows)
    _, i = np.unique(a[:, 0], return_index=True)
    return a[i] / np.array([1, 0.015, 0.015])


def point_disp(case: Path):
    fs = sorted(case.glob("postProcessing/*/solidPointDisplacement*.dat"))
    rows = []
    for f in fs:
        for line in f.read_text().splitlines():
            if line.startswith("#") or not line.strip():
                continue
            v = line.replace("(", " ").replace(")", " ").split()
            rows.append([float(x) for x in v[:3]])
    a = np.array(rows)
    _, i = np.unique(a[:, 0], return_index=True)
    return a[i]


def common(a, b):
    ta = np.round(a[:, 0], 7)
    tb = np.round(b[:, 0], 7)
    t = np.intersect1d(ta, tb)
    return a[np.isin(ta, t)], b[np.isin(tb, t)]


def report(name, a, b, labels):
    a, b = common(a, b)
    for j, lab in enumerate(labels, start=1):
        d = np.abs(a[:, j] - b[:, j])
        scale = np.max(np.abs(b[:, j]))
        print(f"  {name:10s} {lab:6s} n={len(a):5d}  max|diff| = {d.max():.3e}  "
              f"(max|value| {scale:.4g}; rel {d.max() / scale:.2e})  t-range {a[0, 0]:.4f}-{a[-1, 0]:.4f}")


def main():
    kind, A, B = sys.argv[1], Path(sys.argv[2]), Path(sys.argv[3])
    print(f"{kind}: A = {A}  B = {B}")
    if kind == "fsi":
        report("forces", forces(A, "fluid/forces"), forces(B, "fluid/forces"), ["drag", "lift"])
        try:
            report("pointA", point_disp(A), point_disp(B), ["ux", "uy"])
        except Exception as exc:  # noqa: BLE001
            print("  pointA unavailable:", exc)
        pa, pb = plate_motion(A), plate_motion(B)
        worst = 0.0
        n = 0
        for p in pa:
            if p in pb:
                if pa[p].ndim != 2 or pb[p].ndim != 2 or pa[p].shape[1] < 2:
                    continue
                x, y = common(pa[p], pb[p])
                if len(x):
                    worst = max(worst, float(np.max(np.abs(x[:, 1:] - y[:, 1:]))))
                    n = len(x)
        print(f"  interface points: max|position diff| = {worst:.3e} m over {n} steps")
    else:
        report("forces", forces(A, "forces"), forces(B, "forces"), ["drag", "lift"])


if __name__ == "__main__":
    main()
