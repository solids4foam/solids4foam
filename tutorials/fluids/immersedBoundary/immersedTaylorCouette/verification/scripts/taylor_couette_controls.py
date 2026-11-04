#!/usr/bin/env python3
"""Controls of the immersedTaylorCouette verification: the separation of the
mesh and the time step, the steady-state extension and the surface-traction
audit. Everything is computed from the tracked data in reference/results/ (no
CFD, no numpy).

    taylor_couette_controls.py timestep   torque and velocity against the step cap
    taylor_couette_controls.py spatial    error against the mesh at each step cap,
                                          with the rate only where it is supported
    taylor_couette_controls.py all        both, written to reference/results/controls

The rows come from the runs of the driver (./Allverify run --runs dt dt-extended
steady20 st-audit; controls.md lists the generator of every file).

The tutorial step is limited by the Courant number and halves with each mesh
level (mean 0.095, 0.049, 0.025, 0.013 s on Levels 1-4); a maxDeltaT cap is a
step that is constant after the start-up ramp (checked on the written times).

Thresholds, the constants SMALL_FRACTION and RATE_SPREAD below, set when the
step-cap runs were planned and before their rows were interpreted (they are not
tuned: the statements of controls.md are checked against them by
test_taylor_couette_controls.py, and a threshold that is changed changes those
statements):
 - a temporal change is SMALL if the last halving of the step changes the
   torque error by less than 10 % of that error;
 - a rate between meshes is stated only if the temporal change is small on
   every mesh involved, the errors decrease monotonically and the successive
   rates agree within 0.3; otherwise the table says 'trend only'.
"""

import csv
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

RESULTS = Path(__file__).resolve().parent.parent / "reference" / "results"
CONTROLS = RESULTS / "controls"
CAPS = (0.01, 0.005, 0.0025, 0.00125)
SMALL_FRACTION = 0.10
RATE_SPREAD = 0.3


def read_rows(path: Path) -> list:
    with path.open() as f:
        return list(csv.DictReader(f))


def row_key(row: dict) -> tuple:
    cap = float(row["dt_cap"]) if row["dt_cap"] not in ("", None) else 0.0
    return (row["method"], row["estimator"], int(row["level"]), round(cap, 12))


def merged_cases(results: Path = RESULTS) -> list:
    """The tutorial-baseline rows and the control rows, each key once"""
    rows = read_rows(results / "cases.csv")
    extra = results / "controls" / "cases-timestep.csv"
    if extra.exists():
        seen = {row_key(r) for r in rows}
        rows += [r for r in read_rows(extra) if row_key(r) not in seen]
    return rows


def find(rows: list, method: str, estimator: str, level: int, cap: float):
    for r in rows:
        if row_key(r) == (method, estimator, level, round(cap, 12)):
            return r
    return None


def cap_series(rows: list, method: str, estimator: str, level: int,
               quantity: str = "T_relerr") -> list:
    """[(cap, value)] for the unbroken halving sequence 0.01, 0.005, ..."""
    out = []
    for cap in CAPS:
        r = find(rows, method, estimator, level, cap)
        if r is None or r[quantity] == "":
            break
        out.append((cap, float(r[quantity])))
    return out


def temporal_small(series: list) -> dict:
    """The last halving of the step: change and its size against the error"""
    if len(series) < 2:
        return {"change": float("nan"), "relative": float("nan"), "small": None}
    change = series[-1][1] - series[-2][1]
    err = max(abs(series[-1][1]), 1e-30)
    rel = abs(change) / err
    return {"change": change, "relative": rel, "small": rel < SMALL_FRACTION}


def timestep_table(rows: list) -> list:
    out = []
    keys = sorted({(r["method"], r["estimator"], int(r["level"])) for r in rows
                   if r["estimator"] != "FO"})
    for method, est, level in keys:
        s = cap_series(rows, method, est, level)
        if len(s) < 2:
            continue
        t = temporal_small(s)
        tut = find(rows, method, est, level, 0.0)
        out.append({"method": method, "estimator": est, "level": level,
                    "tutorial_error": float(tut["T_relerr"]) if tut else "",
                    **{f"error_cap_{c:g}": next((x[1] for x in s if x[0] == c), "")
                       for c in CAPS},
                    "last_cap": s[-1][0], "last_change": t["change"],
                    "last_change_over_error": t["relative"],
                    "temporal_change_small": "yes" if t["small"] else "no"})
    return out


def observed_rates(errors: list) -> list:
    """log2 of the error ratio between successive levels (h halves)"""
    return [math.log(abs(a) / abs(b), 2.0) if a != 0 and b != 0 else float("nan")
            for a, b in zip(errors, errors[1:])]


def spatial_table(rows: list, levels=(1, 2, 3, 4)) -> list:
    out = []
    methods = sorted({(r["method"], r["estimator"]) for r in rows if r["estimator"] != "FO"})
    for method, est in methods:
        for cap in CAPS[:3]:
            vals, small, steady = [], [], []
            for lv in levels:
                r = find(rows, method, est, lv, cap)
                if r is None:
                    break
                vals.append(float(r["T_relerr"]))
                t = temporal_small(cap_series(rows, method, est, lv))
                small.append(t["small"])
                steady.append(r["steady"] == "1")
            if len(vals) < 3:
                continue
            mags = [abs(x) for x in vals]
            monotone = all(a > b for a, b in zip(mags, mags[1:]))
            rates = observed_rates(vals)
            ok = (all(small) and monotone and all(math.isfinite(x) for x in rates)
                  and (max(rates) - min(rates)) <= RATE_SPREAD)
            out.append({"method": method, "estimator": est, "cap": cap,
                        "levels": len(vals), **{f"error_L{i + 1}": x for i, x in enumerate(vals)},
                        "rates": " ".join(f"{x:.2f}" for x in rates),
                        "temporal_change_small_on_all": "yes" if all(small) else "no",
                        "all_steady": "yes" if all(steady) else "no",
                        "monotone": "yes" if monotone else "no",
                        "statement": "observed rate %.2f" % (sum(rates) / len(rates)) if ok else "trend only"})
    return out


def write_csv(path: Path, rows: list) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    keys = []
    for r in rows:
        for k in r:
            if k not in keys:
                keys.append(k)
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys, restval="", lineterminator="\n")
        w.writeheader()
        for r in rows:
            w.writerow({k: (f"{x:.6g}" if isinstance(x, float) else x) for k, x in r.items()})


def main(argv=None) -> int:
    cmd = (argv or sys.argv[1:] or ["all"])[0]
    if cmd not in ("timestep", "spatial", "all"):
        print(__doc__)
        return 2
    rows = merged_cases()
    if cmd in ("timestep", "all"):
        write_csv(CONTROLS / "timestep-table.csv", timestep_table(rows))
    if cmd in ("spatial", "all"):
        write_csv(CONTROLS / "spatial-table.csv", spatial_table(rows))
    print(f"written to {CONTROLS}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
