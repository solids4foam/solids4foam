#!/usr/bin/env python3
"""Analyse the CFD3 (rigid flag) runs of hron_turek_cfd3.py.

For every run found in verification/work the script extracts, over the closing
window and with the benchmark's definitions (mean = (max+min)/2, amplitude =
(max-min)/2 over the last full period of the lift), the drag and lift, their
pressure and viscous parts, the lift frequency and the Strouhal number. The
period of the lift delimits the periods of all quantities, because the drag
oscillates at twice the lift frequency with alternating troughs.

Three paths are then analysed, each with successive differences, their ratio
and the observed order where it is defined:

- matched: 1x dt=1e-3, 2x dt=5e-4, 4x dt=2.5e-4 (the FSI3 space-time path);
- fixed dt: 1x, 2x, 4x all at dt = 2.5e-4 (the spatial path);
- temporal control on 2x: dt = 5e-4, 2.5e-4, 1.25e-4.

Usage, from the verification directory:  python3 scripts/cfd3_analysis.py
"""

from __future__ import annotations

import csv
import json
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import hron_turek_verification as driver  # noqa: E402
from hron_turek_cfd3 import case_name  # noqa: E402

THICKNESS = 0.015
DIAMETER = 0.1
MEAN_INFLOW = 2.0
WINDOW = 1.0
SMOOTHING_WINDOW = 0.004
# Featflow level 4+0, dt = 0.005 (the finest published CFD3 values)
REFERENCE = {
    "drag_mean": 439.45, "drag_amplitude": 5.6183, "drag_frequency": 4.3956,
    "lift_mean": -11.893, "lift_amplitude": 437.81, "lift_frequency": 4.3956,
}
PATHS = {
    "matched": [(1, 1e-3), (2, 5e-4), (4, 2.5e-4)],
    "fixed_dt": [(1, 2.5e-4), (2, 2.5e-4), (4, 2.5e-4)],
    "temporal_2x": [(2, 5e-4), (2, 2.5e-4), (2, 1.25e-4)],
}
# Quantities whose mean is near zero against the amplitude: no relative change
# or order is meaningful
NEAR_ZERO = {"lift_mean", "lift_p_mean", "lift_v_mean"}


def force_columns(case: Path) -> dict[str, list[float]]:
    path = driver.find_file(case, ("postProcessing/**/force.dat",
                                   "postProcessing/**/forces.dat"), "force history")
    cols: dict[str, list[float]] = {k: [] for k in
        ("t", "drag", "lift", "drag_p", "lift_p", "drag_v", "lift_v")}
    for f in driver.numeric_rows(path):
        v = [float(x) for x in f]
        if len(v) < 10:
            continue
        cols["t"].append(v[0])
        cols["drag"].append(v[1] / THICKNESS)
        cols["lift"].append(v[2] / THICKNESS)
        cols["drag_p"].append(v[4] / THICKNESS)
        cols["lift_p"].append(v[5] / THICKNESS)
        cols["drag_v"].append(v[7] / THICKNESS)
        cols["lift_v"].append(v[8] / THICKNESS)
    return cols


def moving_average(values: list[float], width: int) -> list[float]:
    half, n = width // 2, len(values)
    return [sum(values[max(0, i - half):min(n, i + half + 1)]) /
            (min(n, i + half + 1) - max(0, i - half)) for i in range(n)]


def analyse_run(level: int, delta_t: float, end_time: float | None = None) -> dict | None:
    case = driver.WORK_ROOT / case_name(level, delta_t)
    log = case / "log.solids4Foam"
    if not log.is_file() or not driver.re.search(r"^End\s*$", log.read_text(errors="replace"),
                                                 driver.re.MULTILINE):
        return None
    c = force_columns(case)
    t = c["t"]
    end = t[-1]
    stat = driver.periodic_statistics(t, c["lift"], WINDOW)
    bounds = stat["crossings"]
    out: dict = {"level": level, "delta_t": delta_t, "end_time": end,
                 "periods_in_window": len(bounds) - 1}
    text = (case / "system/blockMeshDict").read_text()
    out["fluid_cells"] = sum(math.prod(int(v) for v in m.group(1).split())
                             for m in driver.re.finditer(r"hex\s+\([^)]*\)\s+\(([^()]+)\)", text))
    out["cores"] = len(list(case.glob("processor*"))) or 1
    out["runtime_s"] = driver.execution_time(case)
    out["core_hours"] = out["cores"] * out["runtime_s"] / 3600
    start = next(i for i, x in enumerate(t) if x >= end - WINDOW)
    width = max(3, int(round(SMOOTHING_WINDOW / (t[1] - t[0]))) | 1)
    for name in ("drag", "lift", "drag_p", "lift_p", "drag_v", "lift_v"):
        label = name.replace("_p", "_p").replace("_v", "_v")
        y = c[name][start:]
        tw = t[start:]

        def seg(values, a=bounds[-2], b=bounds[-1]):
            return [v for x, v in zip(tw, values) if a <= x <= b]

        s = seg(y)
        out[f"{label}_mean"] = 0.5 * (max(s) + min(s))
        out[f"{label}_amplitude"] = 0.5 * (max(s) - min(s))
        # all full periods of the window
        per = [(0.5 * (max(q) + min(q)), 0.5 * (max(q) - min(q)))
               for q in (seg(y, a, b) for a, b in zip(bounds[:-1], bounds[1:]))]
        out[f"{label}_mean_allperiods"] = sum(p[0] for p in per) / len(per)
        out[f"{label}_amplitude_allperiods"] = sum(p[1] for p in per) / len(per)
        out[f"{label}_cycle_spread_mean"] = max(p[0] for p in per) - min(p[0] for p in per)
        out[f"{label}_cycle_spread_amplitude"] = max(p[1] for p in per) - min(p[1] for p in per)
        sm = moving_average(y, width)
        ss = seg(sm)
        out[f"{label}_amplitude_smoothed"] = 0.5 * (max(ss) - min(ss))
        out[f"{label}_mean_smoothed"] = 0.5 * (max(ss) + min(ss))
    f = stat["frequency"]
    out["lift_frequency"] = f
    out["drag_frequency"] = 2 * f
    out["strouhal"] = f * DIAMETER / MEAN_INFLOW
    periods = [b - a for a, b in zip(bounds[:-1], bounds[1:])]
    out["frequency_cycle_spread"] = (max(periods) - min(periods)) / (sum(periods) / len(periods)) ** 2
    # Settling check: lift amplitude of the preceding 2 s of the run against
    # the closing window, from the per-period amplitudes over the last 3 s
    late = driver.periodic_statistics(t, c["lift"], 3.0)["crossings"]
    amps = []
    for a, b in zip(late[:-1], late[1:]):
        q = [v for x, v in zip(t, c["lift"]) if a <= x <= b]
        amps.append(0.5 * (max(q) - min(q)))
    out["lift_amplitude_late_variability"] = max(amps) - min(amps)
    out["lift_amplitude_late_drift"] = amps[-1] - amps[0]
    return out


def successive(values: list[float], noise: float, column: str) -> dict:
    f1, f2, f3 = values
    d21, d32 = f2 - f1, f3 - f2
    out = {"d21": d21, "d32": d32, "ratio": None, "order": None, "status": ""}
    if column in NEAR_ZERO:
        out["status"] = "undefined: mean near zero against the amplitude"
        return out
    if d21 == 0:
        out["status"] = "undefined: no change"
        return out
    ratio = d32 / d21
    out["ratio"] = ratio
    if abs(d32) <= noise:
        out["status"] = "unreliable: the last change is within the cycle-to-cycle spread"
    elif abs(d21) <= noise:
        out["status"] = "unreliable: a change is within the cycle-to-cycle spread"
    elif ratio <= 0:
        out["status"] = "undefined: differences change sign (non-monotone)"
    elif ratio >= 1:
        out["status"] = "undefined: differences do not decrease"
    else:
        out["order"] = math.log(1 / ratio) / math.log(2)
        if out["order"] < 1.5:
            out["status"] = "monotone; sub-nominal (observed order below 1.5)"
        elif out["order"] <= 2.5:
            out["status"] = "monotone; consistent with the formal second order"
        else:
            out["status"] = ("monotone; order above 2, i.e. pre-asymptotic error "
                             "reduction, not an asymptotic order")
    return out


COLUMNS = ["drag_mean", "drag_amplitude", "drag_amplitude_smoothed",
           "lift_mean", "lift_amplitude", "lift_amplitude_smoothed",
           "lift_frequency", "strouhal",
           "drag_p_mean", "drag_v_mean", "drag_p_amplitude", "drag_v_amplitude",
           "lift_p_amplitude", "lift_v_amplitude"]


def noise_of(runs: list[dict], column: str) -> float:
    base = column.replace("_smoothed", "")
    key = base.replace("_mean", "_cycle_spread_mean").replace("_amplitude", "_cycle_spread_amplitude")
    if column.endswith("frequency") or column == "strouhal":
        return max(r["frequency_cycle_spread"] for r in runs) * (
            DIAMETER / MEAN_INFLOW if column == "strouhal" else (2 if column.startswith("drag") else 1))
    return max(r.get(key, 0.0) for r in runs)


def main() -> int:
    results: dict = {"reference": REFERENCE, "runs": {}, "paths": {}}
    for path_name, spec in PATHS.items():
        runs = []
        for level, dt in spec:
            key = case_name(level, dt)
            run = results["runs"].get(key) or analyse_run(level, dt)
            if run is None:
                break
            results["runs"][key] = run
            runs.append(run)
        if len(runs) < 3:
            continue
        table = {}
        for column in COLUMNS:
            values = [r[column] for r in runs]
            entry = successive(values, noise_of(runs, column), column)
            entry["values"] = values
            if column in REFERENCE:
                entry["reference_error_pct"] = [100 * (v - REFERENCE[column]) / abs(REFERENCE[column])
                                                for v in values]
            entry["relative_change_pct"] = (
                None if column in NEAR_ZERO else
                [100 * entry["d21"] / abs(values[0]), 100 * entry["d32"] / abs(values[1])])
            table[column] = entry
        results["paths"][path_name] = {"cases": [case_name(*s) for s in spec], "quantities": table}
    out = driver.OUTPUT_ROOT
    out.mkdir(exist_ok=True)
    (out / "cfd3_study.json").write_text(json.dumps(results, indent=2))
    with (out / "cfd3_study_runs.csv").open("w", newline="") as handle:
        rows = list(results["runs"].values())
        if rows:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
    with (out / "cfd3_study_paths.csv").open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["path", "quantity", "v1", "v2", "v3", "d21", "d32", "ratio",
                         "order", "status", "ref_err1_pct", "ref_err2_pct", "ref_err3_pct"])
        for p, body in results["paths"].items():
            for q, e in body["quantities"].items():
                err = e.get("reference_error_pct", [None] * 3)
                writer.writerow([p, q, *e["values"], e["d21"], e["d32"], e["ratio"],
                                 e["order"], e["status"], *err])
    for p, body in results["paths"].items():
        print(f"\n== {p}: {', '.join(body['cases'])}")
        for q, e in body["quantities"].items():
            v = " ".join(f"{x:11.5g}" for x in e["values"])
            o = f"{e['order']:.2f}" if e["order"] is not None else "  - "
            r = f"{e['ratio']:.3f}" if e["ratio"] is not None else "  -  "
            print(f"{q:26s} {v}  ratio {r} order {o}  {e['status']}")
    print()
    for k, r in results["runs"].items():
        print(k, f"cells {r['fluid_cells']} periods {r['periods_in_window']} "
              f"core-h {r['core_hours']:.2f} late-lift-amp-var {r['lift_amplitude_late_variability']:.2f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
