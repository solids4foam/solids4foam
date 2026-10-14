#!/usr/bin/env python3
"""Compare restarted coupled FSI3 runs with the continuous run (issue #555 / PR #562 check).

For each restart (case, restart time) against the continuous case:
  - the first steps after the restart time: solid SNES norms of the predictor and of the
    first coupling iterations, FSI residual1 per iteration, iterations per step;
  - the tip displacement (solidPointDisplacement) over the overlap: max |difference| of
    u_x and u_y, relative to the u_y amplitude, and the first step with a difference.

    python3 ht_restart_check.py --continuous ctl_1x --restart rst_r4.5=4.5 --restart rst_r4.75=4.75
"""
from __future__ import annotations

import argparse
import re
from pathlib import Path

import numpy as np


def disp(case: Path) -> np.ndarray:
    rows = []
    for f in sorted(case.glob("postProcessing/**/solidPointDisplacement*.dat")):
        for line in f.read_text().splitlines():
            if line.strip() and not line.startswith("#"):
                rows.append([float(v) for v in line.split()[:3]])
    a = np.array(rows)
    _, idx = np.unique(a[::-1, 0], return_index=True)
    return a[::-1][idx]


def steps(case: Path, t0: float, nsteps: int) -> list[tuple[str, list[str]]]:
    """Log lines (SNES norms, FSI residuals) of the first nsteps steps after t0."""
    out, cur = [], None
    pat = re.compile(r"^Time = ([0-9.eE+-]+)(, iteration: \d+)?\s*$")
    for log in sorted(case.glob("log.solids4Foam*")):
        with log.open(errors="replace") as fh:
            for line in fh:
                m = pat.match(line)
                if m:
                    t = float(m.group(1))
                    if t <= t0 + 1e-12:
                        cur = None
                        continue
                    if not m.group(2):
                        if len(out) >= nsteps:
                            return out
                        cur = (m.group(1), [])
                        out.append(cur)
                    elif cur is not None:
                        cur[1].append(line.strip())
                elif cur is not None and ("SNES Function norm" in line or "FSI relative residual1" in line
                                           or "Linear predictor" in line):
                    cur[1].append(line.strip())
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--continuous", type=Path, required=True)
    ap.add_argument("--restart", action="append", required=True, help="case=restart time")
    ap.add_argument("--steps", type=int, default=2)
    a = ap.parse_args()
    c = disp(a.continuous)
    amp = 0.5 * (c[c[:, 0] > c[-1, 0] - 1.0, 2].max() - c[c[:, 0] > c[-1, 0] - 1.0, 2].min())
    print(f"continuous: {a.continuous}, tip u_y amplitude over its last second {amp * 1e3:.3f} mm")
    for spec in a.restart:
        name, t0 = spec.split("=")
        t0 = float(t0)
        r = disp(Path(name))
        r = r[r[:, 0] > t0 + 1e-12]
        cm = {round(t, 7): (x, y) for t, x, y in c}
        common = [(t, x, y, *cm[round(t, 7)]) for t, x, y in r if round(t, 7) in cm]
        d = np.array([[t, x - cx, y - cy] for t, x, y, cx, cy in common])
        print(f"\n== {name}: restart at {t0}, {len(d)} common steps ({d[0, 0]:g} to {d[-1, 0]:g} s)")
        nz = np.nonzero(np.abs(d[:, 1:]).max(axis=1) > 0)[0]
        print(f"  bit-identical tip displacement on {len(d) if not len(nz) else nz[0]} steps;"
              f" max |du_x| {np.abs(d[:, 1]).max():.3e} m, max |du_y| {np.abs(d[:, 2]).max():.3e} m"
              f" ({np.abs(d[:, 2]).max() / amp:.2e} of the u_y amplitude)")
        for k in (1, 10, 100, len(d)):
            if k <= len(d):
                print(f"  after {k:5d} steps (t = {d[k - 1, 0]:.4f}): max |du_y| so far {np.abs(d[:k, 2]).max():.3e} m")
        for label, case in (("continuous", a.continuous), ("restart", Path(name))):
            for t, lines in steps(case, t0, a.steps):
                print(f"  {label:10s} Time = {t}")
                for line in lines:
                    print(f"      {line}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
