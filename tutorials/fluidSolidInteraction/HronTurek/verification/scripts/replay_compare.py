#!/usr/bin/env python3
"""Compare replay runs with different plate wall-pressure conditions.

Prints, for each available level at the matched time steps (and 2x at 2.5e-4),
the complex first harmonic of the force on the flag and on cylinder + flag
(in phase with the tip displacement, in phase with the tip velocity, phase),
the plate mean drag and the pressure work per cycle, for the tutorial's
zeroGradient (tag `_replay`) and the given variant tag, with the successive
differences. Usage: python3 scripts/replay_compare.py [--tag _replay_movingWallPressure]
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import hron_turek_verification as driver  # noqa: E402
import replay_analysis as ra  # noqa: E402

LEVELS = [(1, 1e-3), (2, 5e-4), (4, 2.5e-4)]
KEYS = ["plate_lift_inphase", "plate_lift_quad", "plate_lift_amp1", "plate_drag_mean",
        "total_lift_inphase", "total_lift_quad", "total_lift_amp1", "total_drag_mean",
        "work_per_cycle"]


def phase(r: dict, group: str) -> float:
    return math.degrees(math.atan2(r[f"{group}_lift_quad"], r[f"{group}_lift_inphase"]))


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--tag", default="_replay_movingWallPressure")
    parser.add_argument("--out", default=None)
    args = parser.parse_args()
    data: dict = {}
    for label, tag in (("zeroGradient", "_replay"), (args.tag, args.tag)):
        runs = {}
        for level, dt in LEVELS + [(2, 2.5e-4)]:
            r = ra.analyse_run(level, dt, tag)
            if r:
                r["plate_phase"] = phase(r, "plate")
                r["total_phase"] = phase(r, "total")
                runs[f"{level}x_dt{dt:g}"] = r
        data[label] = runs
    for label, runs in data.items():
        print(f"\n== {label}")
        names = list(runs)
        print(f"{'':22s}" + "".join(f"{n:>16s}" for n in names))
        for k in KEYS + ["plate_phase", "total_phase"]:
            print(f"{k:22s}" + "".join(f"{runs[n][k]:16.4g}" for n in names))
    if args.out:
        Path(args.out).write_text(json.dumps(data, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
