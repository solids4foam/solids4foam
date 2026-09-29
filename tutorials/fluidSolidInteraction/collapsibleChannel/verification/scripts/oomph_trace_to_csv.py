#!/usr/bin/env python3
"""Write the reference CSV from two runs of the oomph-lib reference driver.

usage: oomph_trace_to_csv.py <coarse/trace.dat> [<fine/trace.dat>] <output.csv>

With two runs, the fine run must use half the time step of the coarse one.
Both use BDF2, so the wall positions are Richardson-extrapolated in time,
f = f_fine + (f_fine - f_coarse)/3, at the times of the coarse run. With one
run, its wall positions are written as they are.

trace.dat columns: time, x2 at the wall midpoint, inflow and outflow
centreline velocities, p_ext, then x1 and x2 at 25% and 75% of the elastic
wall and x1 at its midpoint. The CSV holds the vertical displacements,
x2 - 1, of the three wall points.
"""

import math
import sys


def read(path):
    rows = {}
    previous = -math.inf
    for number, line in enumerate(open(path), 1):
        fields = line.split()
        if not fields:
            continue
        if len(fields) < 9:
            sys.exit(f"ERROR: {path}:{number} has {len(fields)} fields, not 9")
        try:
            values = [float(field) for field in fields]
        except ValueError:
            sys.exit(f"ERROR: {path}:{number} is not numeric")
        if not all(math.isfinite(v) for v in values):
            sys.exit(f"ERROR: non-finite values in {path}:{number}")
        if values[0] <= previous:
            sys.exit(f"ERROR: {path}:{number}: times are not increasing")
        previous = values[0]
        rows[round(values[0], 9)] = (
            values[6] - 1.0, values[1] - 1.0, values[8] - 1.0
        )
    if len(rows) < 3:
        sys.exit(f"ERROR: fewer than three samples in {path}")
    return rows


def step(rows, path):
    """The uniform time step of a run"""
    times = sorted(rows)
    steps = [b - a for a, b in zip(times, times[1:])]
    if max(steps) - min(steps) > 1e-6*max(steps):
        sys.exit(f"ERROR: {path} does not have a uniform time step")
    return steps[0]


def main() -> int:
    if len(sys.argv) == 3:
        run = read(sys.argv[1])
        step(run, sys.argv[1])
        with open(sys.argv[2], "w") as handle:
            handle.write(
                "# oomph-lib solution: vertical wall displacement (m) at 25%,"
                " 50% and 75% of the elastic wall\n"
                f"# from {sys.argv[1]}\n"
            )
            handle.write("time,wallQuarter,wallMid,wallThreeQuarter\n")
            for t in sorted(run):
                handle.write(
                    f"{t:.10g}," + ",".join(f"{v:.8f}" for v in run[t]) + "\n"
                )
        return 0
    coarse, fine = read(sys.argv[1]), read(sys.argv[2])
    target = sys.argv[3]
    dt_coarse, dt_fine = step(coarse, sys.argv[1]), step(fine, sys.argv[2])
    if abs(dt_coarse/dt_fine - 2.0) > 1e-6:
        sys.exit(
            f"ERROR: the time steps {dt_coarse:g} and {dt_fine:g} are not in the"
            " ratio 2:1 that the Richardson extrapolation assumes"
        )
    if abs(max(coarse) - max(fine)) > 1e-9:
        sys.exit(
            f"ERROR: the runs end at different times, {max(coarse):g} and"
            f" {max(fine):g}"
        )
    missing = [t for t in coarse if t not in fine]
    if missing:
        sys.exit(
            f"ERROR: {len(missing)} coarse samples, from t = {min(missing):g},"
            f" are missing in {sys.argv[2]}"
        )
    times = sorted(coarse)
    correction = max(
        abs(f - c)/3.0
        for t in times for f, c in zip(fine[t], coarse[t])
    )
    with open(target, "w") as handle:
        handle.write(
            "# oomph-lib reference: vertical wall displacement (m) at 25%, 50%"
            " and 75% of the elastic wall\n"
            f"# Richardson extrapolation in time of {sys.argv[1]} and"
            f" {sys.argv[2]}\n"
            f"# largest Richardson correction: {correction:.3e} m\n"
        )
        handle.write("time,wallQuarter,wallMid,wallThreeQuarter\n")
        for t in times:
            values = [f + (f - c)/3.0 for f, c in zip(fine[t], coarse[t])]
            handle.write(
                f"{t:.10g}," + ",".join(f"{v:.8f}" for v in values) + "\n"
            )
    print(f"{len(times)} times, largest Richardson correction {correction:.3e} m")
    return 0


if __name__ == "__main__":
    sys.exit(main())
