#!/usr/bin/env python3
"""Write a reference CSV from runs of the oomph-lib reference driver.

usage: oomph_trace_to_csv.py <coarse/trace.dat> [<fine/trace.dat>] <output.csv>

With two runs, the fine run must use half the time step of the coarse one.
Both use BDF2, so the tip displacements are Richardson-extrapolated in time,
f = f_fine + (f_fine - f_coarse)/3, at the times of the coarse run. With one
run, its tip displacements are written as they are.

trace.dat columns: time, flux, x and y of the leaflet tip, their time
derivatives, the number of fluid elements and the number of unknowns. The
CSV holds the tip displacements, x - 1 and y - 0.5.
"""

import math
import sys

X0 = 1.0
HEIGHT = 0.5


def read(path):
    rows = {}
    previous = -math.inf
    with open(path) as handle:
        for number, line in enumerate(handle, 1):
            fields = line.split()
            if not fields:
                continue
            if len(fields) != 8:
                sys.exit(f"ERROR: {path}:{number} has {len(fields)} fields, not 8")
            try:
                values = [float(field) for field in fields]
            except ValueError:
                sys.exit(f"ERROR: {path}:{number} is not numeric")
            if not all(math.isfinite(v) for v in values):
                sys.exit(f"ERROR: non-finite values in {path}:{number}")
            if values[0] <= previous:
                sys.exit(f"ERROR: {path}:{number}: times are not increasing")
            previous = values[0]
            rows[round(values[0], 9)] = (values[2] - X0, values[3] - HEIGHT)
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


def write(target, header, rows):
    with open(target, "w") as handle:
        for line in header:
            handle.write(f"# {line}\n")
        handle.write("time,tipDx,tipDy\n")
        for t in sorted(rows):
            handle.write(
                f"{t:.10g}," + ",".join(f"{v:.9f}" for v in rows[t]) + "\n"
            )


def main() -> int:
    if len(sys.argv) == 3:
        run = read(sys.argv[1])
        step(run, sys.argv[1])
        write(
            sys.argv[2],
            ["oomph-lib solution: leaflet tip displacement (m)",
             f"from {sys.argv[1]}"],
            run,
        )
        return 0
    if len(sys.argv) != 4:
        sys.exit(__doc__)
    coarse, fine = read(sys.argv[1]), read(sys.argv[2])
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
    extrapolated = {
        t: tuple(f + (f - c)/3.0 for f, c in zip(fine[t], coarse[t]))
        for t in coarse
    }
    correction = max(
        abs(f - c)/3.0 for t in coarse for f, c in zip(fine[t], coarse[t])
    )
    write(
        sys.argv[3],
        ["oomph-lib reference: leaflet tip displacement (m)",
         f"Richardson extrapolation in time of {sys.argv[1]} and {sys.argv[2]}",
         f"largest Richardson correction: {correction:.3e} m"],
        extrapolated,
    )
    print(f"{len(coarse)} times, largest Richardson correction {correction:.3e} m")
    return 0


if __name__ == "__main__":
    sys.exit(main())
