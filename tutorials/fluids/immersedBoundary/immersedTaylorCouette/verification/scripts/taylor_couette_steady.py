#!/usr/bin/env python3
"""Generator of reference/results/controls/steady-extended.csv and of
history-extended.csv: the steady-state diagnostics of the 20 s runs.

The cases are those of the named run 'steady20' of the driver, which writes
endTime 20 into controlDict:

    ./Allverify run --runs steady20      # about 2 serial hours (15 cases)
    python3 scripts/taylor_couette_steady.py --out steady.csv
    python3 scripts/taylor_couette_steady.py --history-out history.csv

For each completed case of work/cases with an end time of 20 s it reports

 - the mean rotor torque (the force estimator of the method; for incremental,
   penalty and ghostCell the forcing torque as written) in the windows 8-10,
   13-15 and 18-20 s, as a difference from the exact torque, (T - T_exact)
   / T_exact (the plain mean of the written samples);
 - the endpoint drift between the windows 13-15 and 18-20 s,
   |T_18_20 - T_13_15| / |T_18_20|, in the torque itself;
 - the relative L2 change of the azimuthal-plane velocity over the fluid cells
   between the snapshots 10 and 15 s and 15 and 20 s, as a fraction of
   omega R1;
 - the rate of change of the fluid angular momentum between the snapshots as a
   fraction of |T_exact|: a steady interior needs it to vanish before the
   forcing torque can be equated with the torque on the rotor;
 - steady_by_criterion: the drift and the change of the velocity between 15 and
   20 s are both below the predefined 1e-3 of the driver.

Only the standard library and the analysis functions of the driver are used;
nothing is written into the tutorial. The case column is the label of the
tracked CSV, TC-<method>[-<estimator>]-L<level>-np<ranks>[-dt<cap>]-t20, where
<cap> is the digits after the decimal point of the step cap (0.0025 -> 0025).
"""

import argparse
import csv
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import immersed_taylor_couette_verification as v  # noqa: E402

END_TIME = 20.0
WINDOWS = ((8.0, 10.0), (13.0, 15.0), (18.0, 20.0))
SNAPSHOTS = (10, 15, 20)
COLUMNS = ["case", "T_8_10", "T_13_15", "T_18_20", "drift_13_15_to_18_20",
           "dU_10_15", "dU_15_20", "dLdt_10_15_rel", "dLdt_15_20_rel",
           "steady_by_criterion"]
# the cases whose torque histories are tracked in history-extended.csv
HISTORY_CASES = ("TC-incremental-L1-np1-dt0025-t20",
                 "TC-incremental-L1-np1-dt005-t20",
                 "TC-incremental-L1-np1-dt01-t20",
                 "TC-penalty-L2-np1-dt005-t20")


def legacy_label(method: str, est: str, level: int, ranks: int,
                 dt_cap: float, end_time: float = END_TIME) -> str:
    """The case label of the tracked CSV files"""
    cap = f"-dt{v.fmt(dt_cap).split('.')[1]}" if dt_cap else ""
    return (f"TC-{method}{'-' + est if est else ''}-L{level}-np{ranks}{cap}"
            f"-t{v.fmt(end_time)}")


def window_mean(t: list, y: list, a: float, b: float) -> float:
    x = [q for tt, q in zip(t, y) if a <= tt <= b]
    return sum(x) / len(x)


def steady_row(case: Path, label: str, level: int, p: dict) -> dict:
    """The diagnostics of one completed 20 s case"""
    rotor = v.read_force_file(case / "postProcessing/immersedBoundary/0/rotor.dat")
    t, tz = rotor["t"], [q[2] for q in rotor["torque"]]
    exact = v.exact_rotor_torque(p)
    means = [window_mean(t, tz, a, b) for a, b in WINDOWS]
    n, h, centres = v.cell_centres(level)
    fluid = [i for i, (x, y) in enumerate(centres)
             if p["R1"] < math.hypot(x, y) < p["R2"]]
    cx = [centres[i] for i in fluid]
    snaps = {}
    for tk in SNAPSHOTS:
        u = v.read_vector_field(case / f"{tk}/U")
        snaps[tk] = [u[i] for i in fluid]

    def d_u(a, b):
        s = sum((q[0] - r[0]) ** 2 + (q[1] - r[1]) ** 2
                for q, r in zip(snaps[a], snaps[b]))
        return math.sqrt(s / len(fluid)) / (p["omega"] * p["R1"])

    def angular_momentum(tk):
        return sum(x * u[1] - y * u[0] for (x, y), u in zip(cx, snaps[tk])) \
            * h * h * p["Lz"] * p["rho"]
    dl1 = (angular_momentum(15) - angular_momentum(10)) / 5 / abs(exact)
    dl2 = (angular_momentum(20) - angular_momentum(15)) / 5 / abs(exact)
    drift = abs(means[2] - means[1]) / abs(means[2])
    steady = drift < v.STEADY_DRIFT and d_u(15, 20) < v.STEADY_DU
    rel = [(m - exact) / exact for m in means]
    return {"case": label, "T_8_10": f"{rel[0]:+.4e}",
            "T_13_15": f"{rel[1]:+.4e}", "T_18_20": f"{rel[2]:+.4e}",
            "drift_13_15_to_18_20": f"{drift:.3e}",
            "dU_10_15": f"{d_u(10, 15):.3e}", "dU_15_20": f"{d_u(15, 20):.3e}",
            "dLdt_10_15_rel": f"{dl1:+.3e}", "dLdt_15_20_rel": f"{dl2:+.3e}",
            "steady_by_criterion": int(steady)}


def history_rows(case: Path, label: str) -> list:
    """The torque series of a case, as the raw time and %.7e of the torque"""
    rotor = v.read_force_file(case / "postProcessing/immersedBoundary/0/rotor.dat")
    text = (case / "postProcessing/immersedBoundary/0/rotor.dat").read_text()
    times = [line.split()[0] for line in text.splitlines()
             if line and not line.startswith("#")]
    return [[label, tt, f"{q[2]:.7e}"] for tt, q in zip(times, rotor["torque"])]


def completed_t20(manifest: dict, cases_dir: Path) -> list:
    """(label, case path, entry) of the completed 20 s cases"""
    out = []
    for name, e in sorted(manifest["cases"].items()):
        if e.get("state") == "completed" and e.get("end_time") == END_TIME \
                and (cases_dir / name).exists():
            out.append((legacy_label(e["method"], e["estimator"], e["level"],
                                     e["ranks"], e.get("dt_cap", 0.0)),
                        cases_dir / name, e))
    return sorted(out, key=lambda x: x[0])


def write_csv(path, header: list, rows: list) -> None:
    out = sys.stdout if path is None else open(path, "w", newline="")
    try:
        w = csv.writer(out, lineterminator="\n")
        w.writerow(header)
        w.writerows(rows)
    finally:
        if path is not None:
            out.close()


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", help="steady-extended.csv (default: stdout)")
    ap.add_argument("--history-out",
                    help="write the torque histories of HISTORY_CASES here")
    ap.add_argument("--pattern", default="", help="only labels containing this")
    ap.add_argument("--cases-dir", type=Path, default=v.CASES_DIR)
    args = ap.parse_args(argv)
    manifest = v.load_manifest()
    p = v.parameters()
    found = [x for x in completed_t20(manifest, args.cases_dir)
             if args.pattern in x[0]]
    if not found:
        print("No completed 20 s case in work/: run ./Allverify run --runs "
              "steady20", file=sys.stderr)
        return 1
    if args.history_out:
        rows = [r for label, case, _ in found if label in HISTORY_CASES
                for r in history_rows(case, label)]
        write_csv(args.history_out, ["case", "time", "rotor_torque_z"], rows)
    else:
        rows = [steady_row(case, label, e["level"], p) for label, case, e in found]
        write_csv(args.out, COLUMNS, [[r[c] for c in COLUMNS] for r in rows])
    return 0


if __name__ == "__main__":
    sys.exit(main())
