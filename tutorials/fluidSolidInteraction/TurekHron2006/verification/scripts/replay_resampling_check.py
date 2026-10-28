#!/usr/bin/env python3
"""How far the replayed wall of a mesh departs from the source 2x polyline.

replay_table.py samples the coefficients of the 2x replay table (linear along
the 2x polyline) at the plate points of the target mesh; the wall imposed on
the target is then the piecewise-linear interpolant of those samples, which
differs from the displaced 2x polyline wherever a 2x knot falls inside a
target face. Over 32 phases of the period this reports the largest distance
from a displaced 2x knot to the displaced target polyline (the largest
deviation between the two piecewise-linear walls), its rms over the knots,
both relative to the tip amplitude, and the change in the first-harmonic tip
y phasor at x = 0.6 m.

    python3 replay_resampling_check.py --source plateMotion_2x.tab <case>/constant/plateMotion.tab ...
"""
from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from replay_table import read_source  # noqa: E402


def shape(xy, coef, omega, t):
    nc = coef.shape[2]
    out = xy.copy()
    for comp in (0, 1):
        c = coef[:, comp, :]
        v = c[:, 0].copy()
        for h in range(1, nc // 2 + 1):
            v += c[:, 2 * h - 1] * math.cos(h * omega * t) + c[:, 2 * h] * math.sin(h * omega * t)
        out[:, comp] += v
    return out


def dist_to_polyline(p, path):
    seg = path[1:] - path[:-1]
    l2 = (seg ** 2).sum(1)
    w = np.clip(((p[:, None, :] - path[None, :-1]) * seg[None]).sum(2) / l2[None], 0, 1)
    foot = path[None, :-1] + w[..., None] * seg[None]
    return np.min(np.linalg.norm(foot - p[:, None, :], axis=2), axis=1)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--source", type=Path, required=True)
    ap.add_argument("tables", nargs="+", type=Path)
    a = ap.parse_args()
    omega, sxy, scoef = read_source(a.source)
    tip = np.argmin(np.hypot(sxy[:, 0] - 0.6, sxy[:, 1] - 0.19))
    amp = math.hypot(scoef[tip, 1, 1], scoef[tip, 1, 2])
    print(f"tip y first-harmonic amplitude {amp * 1e3:.3f} mm")
    for tab in a.tables:
        _, txy, tcoef = read_source(tab)
        worst, sq, n = 0.0, 0.0, 0
        for k in range(32):
            t = 2 * math.pi / omega * k / 32
            src = shape(sxy, scoef, omega, t)
            tgt = shape(txy, tcoef, omega, t)
            d = dist_to_polyline(src, tgt)
            worst = max(worst, float(d.max()))
            sq += float((d ** 2).sum())
            n += len(d)
        ttip = np.argmin(np.hypot(txy[:, 0] - 0.6, txy[:, 1] - 0.19))
        dph = complex(tcoef[ttip, 1, 1], tcoef[ttip, 1, 2]) - complex(scoef[tip, 1, 1], scoef[tip, 1, 2])
        print(f"{tab}: {len(txy)} points  max deviation {worst * 1e6:.2f} um ({worst / amp:.2e} of the tip "
              f"amplitude)  rms {math.sqrt(sq / n) * 1e6:.2f} um  tip phasor change {abs(dph) / amp:.1e}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
