#!/usr/bin/env python3
"""Write a reference CSV from runs of the oomph-lib reference driver.

usage: oomph_trace_to_csv.py [--space-ratio R] RUN [RUN] <output.csv>
       oomph_trace_to_csv.py --static [--space-ratio R] RUN

RUN is a trace.dat, or FINE.dat:COARSE.dat for two runs at the same time
step on meshes refined in the ratio R (default 1.5, i.e. --refine 3 and
--refine 2). The oomph-lib solution converges at first order in space, so a
pair is extrapolated in space as f = f_fine + (f_fine - f_coarse)/(R - 1),
i.e. 3 f_3 - 2 f_2 for --refine 3 and 2.

With two RUNs, the second must use half the time step of the first. Both use
BDF2, so the wall positions are Richardson-extrapolated in time,
f = f_fine + (f_fine - f_coarse)/3, at the times of the coarse run. With one
RUN, its wall positions are written as they are.

--static prints the wall displacements of the last sample of RUN, i.e. the
static deflection of a --q 0 --steady run, for the static study's reference
in collapsibleChannel_verification_references.json.

trace.dat columns: time, x2 at the wall midpoint, inflow and outflow
centreline velocities, p_ext, then x1 and x2 at 25% and 75% of the elastic
wall and x1 at its midpoint. The CSV holds the vertical displacements,
x2 - 1, of the three wall points.
"""

import argparse
import math
import sys


def read(path, steady=False):
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
        # A steady run writes t = 0 for every load step: key it by line
        key = number if steady else round(values[0], 9)
        if not steady:
            if values[0] <= previous:
                sys.exit(f"ERROR: {path}:{number}: times are not increasing")
            previous = values[0]
        rows[key] = (values[6] - 1.0, values[1] - 1.0, values[8] - 1.0)
    if len(rows) < (1 if steady else 3):
        sys.exit(f"ERROR: too few samples in {path}")
    return rows


def load(run, ratio, steady=False):
    """A run, or the spatial extrapolation of a FINE:COARSE pair, with the
    largest spatial correction (0 for a single run)"""
    if ":" not in run:
        rows = read(run, steady)
        return ({0: rows[max(rows)]} if steady else rows), 0.0
    fine_path, coarse_path = run.split(":", 1)
    fine, coarse = read(fine_path, steady), read(coarse_path, steady)
    if steady:
        fine, coarse = {0: fine[max(fine)]}, {0: coarse[max(coarse)]}
    if sorted(fine) != sorted(coarse):
        sys.exit(
            f"ERROR: {fine_path} and {coarse_path} do not have the same"
            " sample times"
        )
    factor = 1.0/(ratio - 1.0)
    rows = {
        t: tuple(f + factor*(f - c) for f, c in zip(fine[t], coarse[t]))
        for t in fine
    }
    correction = max(
        factor*abs(f - c) for t in fine for f, c in zip(fine[t], coarse[t])
    )
    return rows, correction


def step(rows, path):
    """The uniform time step of a run"""
    times = sorted(rows)
    steps = [b - a for a, b in zip(times, times[1:])]
    if max(steps) - min(steps) > 1e-6*max(steps):
        sys.exit(f"ERROR: {path} does not have a uniform time step")
    return steps[0]


HEADER = "wall displacement (m) at 25%, 50% and 75% of the elastic wall\n"
COLUMNS = "time,wallQuarter,wallMid,wallThreeQuarter\n"


def source_note(run, ratio, correction):
    if ":" not in run:
        return f"# from {run}\n"
    return (
        f"# first-order extrapolation in space (ratio {ratio:g}) of"
        f" {run.replace(':', ' and ')}; largest spatial correction:"
        f" {correction:.3e} m\n"
    )


def write_rows(handle, rows, times):
    for t in times:
        handle.write(
            f"{t:.10g}," + ",".join(f"{v:.8f}" for v in rows[t]) + "\n"
        )


def main() -> int:
    parser = argparse.ArgumentParser(
        description=__doc__.split("\n\n")[0],
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__.split("\n\n", 1)[1],
    )
    parser.add_argument("--space-ratio", type=float, default=1.5)
    parser.add_argument("--static", action="store_true")
    parser.add_argument("paths", nargs="+")
    args = parser.parse_args()
    ratio = args.space_ratio
    if ratio <= 1.0:
        sys.exit("ERROR: --space-ratio must be greater than 1")

    if args.static:
        if len(args.paths) != 1:
            sys.exit("ERROR: --static takes one RUN")
        rows, correction = load(args.paths[0], ratio, steady=True)
        names = ("wallQuarter", "wallMid", "wallThreeQuarter")
        for name, value in zip(names, rows[0]):
            print(f'"{name}": {value:.8f}')
        if correction:
            print(f"largest spatial correction {correction:.3e} m")
        return 0

    if len(args.paths) == 2:
        run, target = args.paths
        rows, correction = load(run, ratio)
        step(rows, run)
        with open(target, "w") as handle:
            handle.write("# oomph-lib solution: vertical " + HEADER)
            handle.write(source_note(run, ratio, correction))
            handle.write(COLUMNS)
            write_rows(handle, rows, sorted(rows))
        return 0

    if len(args.paths) != 3:
        parser.error("give one or two RUNs and the output CSV")
    coarse_run, fine_run, target = args.paths
    coarse, coarse_space = load(coarse_run, ratio)
    fine, fine_space = load(fine_run, ratio)
    dt_coarse, dt_fine = step(coarse, coarse_run), step(fine, fine_run)
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
            f" are missing in {fine_run}"
        )
    times = sorted(coarse)
    correction = max(
        abs(f - c)/3.0
        for t in times for f, c in zip(fine[t], coarse[t])
    )
    rows = {
        t: tuple(f + (f - c)/3.0 for f, c in zip(fine[t], coarse[t]))
        for t in times
    }
    with open(target, "w") as handle:
        handle.write("# oomph-lib reference: vertical " + HEADER)
        handle.write(
            "# Richardson extrapolation in time of the two runs below\n"
        )
        handle.write(source_note(coarse_run, ratio, coarse_space))
        handle.write(source_note(fine_run, ratio, fine_space))
        handle.write(f"# largest Richardson correction: {correction:.3e} m\n")
        handle.write(COLUMNS)
        write_rows(handle, rows, times)
    print(
        f"{len(times)} times, largest Richardson correction {correction:.3e} m,"
        f" largest spatial correction {max(coarse_space, fine_space):.3e} m"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
