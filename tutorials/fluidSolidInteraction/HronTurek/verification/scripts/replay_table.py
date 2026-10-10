#!/usr/bin/env python3
"""Write constant/plateMotion.tab of a replay case for any fluid mesh.

The source is the replay table of the coupled 2x motion
(reference/energy_balance/plateMotion_2x.tab: one row per 2x plate point,
keyed by its reference xy in 10 um, with the mean and four harmonics in x and
y). Its rows are ordered along the flag surface (lower face, end face, upper
face), and the coefficients are interpolated linearly along that polyline to
the plate points of the target mesh, as hron_turek_replay.py does for the
standard levels. The replayed boundary is therefore the same piecewise-linear
shape (the 2x polyline) on every mesh whose plate points lie on the flag.

    python3 replay_table.py --source plateMotion_2x.tab --mesh <case>/constant/polyMesh \
        --out <case>/constant/plateMotion.tab
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from hron_turek_replay import plate_polyline, interpolate, write_table  # noqa: E402

KEY = 1e5
Y_LO, Y_HI, X_TIP = 0.19, 0.21, 0.6


def read_source(path: Path):
    lines = path.read_text().split("\n")
    omega, n, ncoef = lines[0].split()
    keys, coef = [], []
    for line in lines[1:]:
        if line.strip():
            v = line.split()
            keys.append((int(v[0]), int(v[1])))
            coef.append([float(x) for x in v[2:]])
    xy = np.array(keys, float) / KEY
    coef = np.array(coef).reshape(len(keys), 2, int(ncoef))
    lo = [i for i, p in enumerate(xy) if abs(p[1] - Y_LO) < 1e-7]
    hi = [i for i, p in enumerate(xy) if abs(p[1] - Y_HI) < 1e-7]
    end = [i for i, p in enumerate(xy) if abs(p[0] - X_TIP) < 1e-7 and Y_LO + 1e-7 < p[1] < Y_HI - 1e-7]
    order = (sorted(lo, key=lambda i: xy[i][0]) + sorted(end, key=lambda i: xy[i][1])
             + sorted(hi, key=lambda i: -xy[i][0]))
    if len(order) != len(xy):
        raise SystemExit(f"source rows not all on the flag surface ({len(order)} of {len(xy)})")
    return float(omega), xy[order], coef[order]


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--source", type=Path, required=True)
    ap.add_argument("--mesh", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args()
    omega, path, coef = read_source(a.source)
    target, _ = plate_polyline(a.mesh)
    write_table(a.out, target, interpolate(path, coef, target), omega)
    print(f"{a.out}: {len(target)} plate points")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
