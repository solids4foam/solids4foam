#!/usr/bin/env python3
"""Three-level analysis of the FSI3 IQN-ILS mesh study.

Reads the completed runs of `Allverify --levels 1,2,4` (and, optionally, a run
on the 2x mesh with a smaller time step, `--delta-t`) from verification/work,
and writes, for every benchmark quantity:

- the value on each level, using the driver's extraction (last full u_y
  period of the closing window), and its cycle-to-cycle spread over all the
  full periods of the window;
- the successive differences 1x -> 2x and 2x -> 4x, their ratio, and the
  observed order along the refinement path where it is defined;
- the error against the Featflow level-4 reference (dt = 0.00025 s) and
  whether it decreases monotonically;
- the coupling iterations, final interface residuals and cost of each run.

The time step is halved with the mesh spacing, so the observed order is that
of the combined space-time refinement path, not a purely spatial order.

Usage, from the verification directory:

    python3 scripts/fsi3_refinement_analysis.py \
        [--diagnostic iqnils_mesh_2x_dt0.00025] [--diagnostic iqnils_mesh_2x_tol1e-06]
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import statistics
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import hron_turek_verification as driver  # noqa: E402

QUANTITIES = driver.QUANTITIES
STATISTICS = driver.STATISTICS
# Means that are close to zero against the amplitude: their relative errors
# and orders are not meaningful
NEAR_ZERO = {"uy_mean", "lift_mean"}
LEVELS = (1, 2, 4)
# Observed orders within this band are taken as consistent with the formal
# second order of the discretisation
FORMAL_ORDER_BAND = (1.5, 2.5)


def per_cycle(time: list[float], values: list[float], bounds: list[float]) -> list[tuple[float, float]]:
    """(mean, amplitude) over each full period delimited by bounds."""
    cycles = []
    for start, end in zip(bounds[:-1], bounds[1:]):
        segment = [v for t, v in zip(time, values) if start <= t <= end]
        cycles.append((0.5 * (max(segment) + min(segment)),
                       0.5 * (max(segment) - min(segment))))
    return cycles


# Moving-average window for the noise diagnostic. The interface tolerance
# leaves step-to-step noise in the forces, which grows as the time step falls
# and inflates extrema-based amplitudes. A 4 ms centred average attenuates the
# u_y fundamental (about 5.5 Hz) by less than 0.1%.
SMOOTHING_WINDOW = 0.004


def moving_average(values: list[float], width: int) -> list[float]:
    half = width // 2
    out = []
    for i in range(len(values)):
        a, b = max(0, i - half), min(len(values), i + half + 1)
        out.append(sum(values[a:b]) / (b - a))
    return out


def residual_summary(case: Path, window_start: float, tolerance: float,
                     end_time: float) -> dict:
    """Final interface residual of each coupled step, overall and in the window."""
    path = case / "postProcessing/fsiResiduals.dat"
    final: dict[float, tuple[float, float]] = {}
    for fields in driver.numeric_rows(path):
        final[float(fields[0])] = (float(fields[1]), float(fields[2]))
    start = driver.dictionary_scalar(case / "constant/fsiProperties.iqnils", "couplingStartTime")
    # A run that stopped early also records the step it failed in, which is
    # not part of the analysed history
    coupled = {t: v for t, v in final.items() if start < t <= end_time + 1e-9}
    window = {t: v for t, v in coupled.items() if t >= window_start}
    def stats(steps: dict) -> dict:
        iterations = [n for n, _ in steps.values()]
        residuals = [r for _, r in steps.values()]
        return {
            "steps": len(steps),
            "mean_iterations": sum(iterations) / len(iterations),
            "max_iterations": max(iterations),
            "max_final_residual": max(residuals),
            "steps_above_tolerance": sum(r > tolerance for r in residuals),
        }
    return {"coupled": stats(coupled), "window": stats(window)}


# Start of the late-time variability measure. After the coupling starts at
# t = 2 s, the amplitudes saturate by t = 4 s, but the means, amplitudes and
# frequency keep wandering slowly, by about 1% over periods of a second or
# more, which the 1 s analysis window does not see. The range of the per-period
# values from this time to the end of the run is used as the variability of a
# run; differences between runs smaller than it are not resolved.
VARIABILITY_START = 4.5


def analyse_case(name: str, spec: dict, window: float, end_time: float,
                 variability_start: float = VARIABILITY_START) -> dict:
    case = driver.WORK_ROOT / name
    row = driver.extract(case, spec, window, end_time)
    driver.periodicity_failures(row, spec["periodicity"]["amplitudeTolerance"])
    time_d, ux, uy = driver.displacement_history(case)
    time_f, drag, lift = driver.force_history(case, spec["thickness_m"])
    fundamental = driver.periodic_statistics(time_d, uy, window)["crossings"]
    histories = {"ux": (time_d, ux), "uy": (time_d, uy),
                 "drag": (time_f, drag), "lift": (time_f, lift)}
    result: dict = {"case": name, "values": {}, "cycle_spread": {}, "cycles": {},
                    "amplitude_change": {}}
    for quantity, (time, values) in histories.items():
        cycles = per_cycle(time, values, fundamental)
        result["cycles"][quantity] = cycles
        for index, statistic in enumerate(("mean", "amplitude")):
            series = [c[index] for c in cycles]
            # Range over the full periods of the closing window
            result["cycle_spread"][f"{quantity}_{statistic}"] = max(series) - min(series)
        result["amplitude_change"][quantity] = row[f"{quantity}_amplitude_change"]
        for statistic in STATISTICS:
            result["values"][f"{quantity}_{statistic}"] = row[f"{quantity}_{statistic}"]
    # Noise diagnostic: statistics of the smoothed histories over the same last
    # full u_y period, and the rms of the removed step-to-step component
    result["smoothed"], result["noise_rms"], result["noise_max"] = {}, {}, {}
    for quantity, (time, values) in histories.items():
        start = next(i for i, t in enumerate(time) if t >= end_time - window)
        t, y = time[start:], values[start:]
        width = max(3, int(round(SMOOTHING_WINDOW / (t[1] - t[0]))) | 1)
        smooth = moving_average(y, width)
        segment = [v for tv, v in zip(t, smooth) if fundamental[-2] <= tv <= fundamental[-1]]
        result["smoothed"][f"{quantity}_mean"] = 0.5 * (max(segment) + min(segment))
        result["smoothed"][f"{quantity}_amplitude"] = 0.5 * (max(segment) - min(segment))
        noise = [a - b for a, b in zip(y, smooth)]
        result["noise_rms"][quantity] = math.sqrt(sum(v * v for v in noise) / len(noise))
        result["noise_max"][quantity] = max(abs(v) for v in noise)
    periods = [b - a for a, b in zip(fundamental[:-1], fundamental[1:])]
    spread = (max(periods) - min(periods)) / statistics.mean(periods) ** 2
    for quantity in QUANTITIES:
        factor = 2 if quantity in ("ux", "drag") else 1
        result["cycle_spread"][f"{quantity}_frequency"] = factor * spread
    # Late-time variability: the range of the per-period values over all the
    # full u_y periods from variability_start to the end of the run
    late = driver.periodic_statistics(time_d, uy, end_time - variability_start)["crossings"]
    result["late_periods"] = len(late) - 1
    result["late_variability"] = {}
    result["late_variability_smoothed"] = {}
    for quantity, (time, values) in histories.items():
        cycles = per_cycle(time, values, late)
        start = next(i for i, t in enumerate(time) if t >= late[0] - 0.01)
        width = max(3, int(round(SMOOTHING_WINDOW / (time[1] - time[0]))) | 1)
        smooth_cycles = per_cycle(time[start:], moving_average(values[start:], width), late)
        for index, statistic in enumerate(("mean", "amplitude")):
            series = [c[index] for c in cycles]
            result["late_variability"][f"{quantity}_{statistic}"] = max(series) - min(series)
            series = [c[index] for c in smooth_cycles]
            result["late_variability_smoothed"][f"{quantity}_{statistic}"] = max(series) - min(series)
        factor = 2 if quantity in ("ux", "drag") else 1
        frequencies = [factor / (b - a) for a, b in zip(late[:-1], late[1:])]
        result["late_variability"][f"{quantity}_frequency"] = max(frequencies) - min(frequencies)
    result["periods_in_window"] = len(fundamental) - 1
    tolerance = driver.dictionary_scalar(case / "constant/fsiProperties.iqnils", "outerCorrTolerance")
    result["outer_corr_tolerance"] = tolerance
    result["residuals"] = residual_summary(case, end_time - window, tolerance, end_time)
    result["cores"] = len(list(case.glob("processor*"))) or 1
    result["execution_time_s"] = driver.execution_time(case)
    result["core_hours"] = result["cores"] * result["execution_time_s"] / 3600
    control = case / "system/controlDict"
    result["delta_t"] = driver.dictionary_scalar(control, "deltaT")
    counts = {}
    for region in ("fluid", "solid"):
        text = (case / "system" / region / "blockMeshDict").read_text()
        counts[region] = sum(
            math.prod(int(v) for v in m.group(1).split())
            for m in driver.re.finditer(r"hex\s+\([^)]*\)\s+\(([^()]+)\)", text)
        )
    result["fluid_cells"] = counts["fluid"]
    result["solid_cells"] = counts["solid"]
    log = (case / "log.solids4Foam").read_text(errors="replace")
    result["status"] = "completed" if driver.re.search(r"^End\s*$", log, driver.re.MULTILINE) else "incomplete"
    return result


def write_history(case_name: str, spec: dict, start: float, path: Path,
                  sampling: float = 0.001) -> None:
    """Write the histories from start to the end of a run, subsampled.

    The stored histories make the extraction reproducible without the run
    directories; 1 ms sampling keeps the files small. The noise diagnostic
    needs the full sampling, so it is computed from the run directories.
    """
    case = driver.WORK_ROOT / case_name
    time_d, ux, uy = driver.displacement_history(case)
    time_f, drag, lift = driver.force_history(case, spec["thickness_m"])
    forces = {round(t, 9): (d, l) for t, d, l in zip(time_f, drag, lift)}
    with path.open("w", newline="") as handle:
        handle.write(f"# {case_name}: point-A displacement (m) and drag and lift on the "
                     f"cylinder and plate (N/m), subsampled to {sampling:g} s\n")
        writer = csv.writer(handle)
        writer.writerow(["time", "drag", "lift", "ux", "uy"])
        for t, x, y in zip(time_d, ux, uy):
            if t < start - 1e-9 or abs(t / sampling - round(t / sampling)) > 1e-6:
                continue
            force = forces.get(round(t, 9))
            if force is not None:
                writer.writerow([f"{t:.5f}", f"{force[0]:.6e}", f"{force[1]:.6e}",
                                 f"{x:.6e}", f"{y:.6e}"])


def classify(column: str, f1: float, f2: float, f4: float, noise: float) -> dict:
    """Successive differences, their ratio and the observed order."""
    d21, d42 = f2 - f1, f4 - f2
    out = {"d21": d21, "d42": d42, "ratio": None, "order": None,
           "extrapolated": None, "extrapolation_meaningful": False, "status": ""}
    if column in NEAR_ZERO:
        out["status"] = "undefined: quantity near zero"
        return out
    if d21 == 0:
        out["status"] = "undefined: no change from 1x to 2x"
        return out
    ratio = d42 / d21
    out["ratio"] = ratio
    if abs(d42) <= noise:
        out["status"] = "unreliable: 2x->4x change within the late-time variability"
        return out
    if ratio <= 0:
        out["status"] = "undefined: differences change sign (non-monotone)"
        return out
    if ratio >= 1:
        out["status"] = "undefined: differences do not decrease"
        return out
    order = math.log(1 / ratio) / math.log(2)
    out["order"] = order
    out["extrapolated"] = f4 + d42 / (2 ** order - 1)
    # The schemes are formally second order in space and time. An observed
    # order far from that means the sequence is not in the asymptotic range,
    # and a Richardson extrapolation with it is not meaningful
    out["extrapolation_meaningful"] = FORMAL_ORDER_BAND[0] <= order <= FORMAL_ORDER_BAND[1]
    if abs(d21) <= noise:
        out["status"] = "unreliable: 1x->2x change within the late-time variability"
    elif order < FORMAL_ORDER_BAND[0]:
        out["status"] = "monotone; order well below the formal 2: not asymptotic"
    elif order > FORMAL_ORDER_BAND[1]:
        out["status"] = ("monotone; order above the formal 2: possible error "
                         "cancellation, not demonstrably asymptotic")
    else:
        out["status"] = "monotone; order consistent with the formal 2"
    return out


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--coupling", default="iqnils")
    parser.add_argument("--diagnostic", action="append", default=[],
                        help="additional run on the 2x mesh, e.g. with a smaller time step "
                        "(--delta-t) or interface tolerance (--outer-corr-tolerance), as name "
                        "or name:end_time for a run that stopped early; repeatable")
    parser.add_argument("--output", type=Path, default=driver.OUTPUT_ROOT)
    parser.add_argument("--history-start", type=float, default=VARIABILITY_START,
                        help="also write each run's histories from this time, at 1 ms")
    args = parser.parse_args()
    spec = json.loads(driver.REFERENCE_FILE.read_text())
    window = spec["mesh"]["analysisWindow"]
    end_time = spec["mesh"]["endTime"]
    references = {column: entry["value"] for column, entry in spec["references"].items()}
    featflow = {row["level"]: row for row in spec["featflowTables"]["0.00025"]}

    runs = {level: analyse_case(f"{args.coupling}_mesh_{level}x", spec, window, end_time)
            for level in LEVELS}
    # A diagnostic is given as name or name:end_time, the latter for a run that
    # stopped before the study's end time; its window then closes there
    diagnostics = {}
    for item in args.diagnostic:
        name, _, stop = item.partition(":")
        diagnostics[name] = analyse_case(name, spec, window, float(stop) if stop else end_time)
        diagnostics[name]["end_time"] = float(stop) if stop else end_time

    table = []
    for quantity in QUANTITIES:
        for statistic in STATISTICS:
            column = f"{quantity}_{statistic}"
            f1, f2, f4 = (runs[level]["values"][column] for level in LEVELS)
            # Differences smaller than the late-time variability of the two
            # finer runs are not resolved
            noise = max(runs[2]["late_variability"][column], runs[4]["late_variability"][column])
            entry = {"quantity": column, "near_zero": column in NEAR_ZERO,
                     "value_1x": f1, "value_2x": f2, "value_4x": f4,
                     "cycle_spread_1x": runs[1]["cycle_spread"][column],
                     "cycle_spread_2x": runs[2]["cycle_spread"][column],
                     "cycle_spread_4x": runs[4]["cycle_spread"][column],
                     **{f"late_variability_{level}x": runs[level]["late_variability"][column]
                        for level in LEVELS}}
            entry.update(classify(column, f1, f2, f4, noise))
            reference = references[column]
            entry["featflow_l4"] = reference
            entry["featflow_l3"] = featflow["3+0"][column]
            for level, value in zip(LEVELS, (f1, f2, f4)):
                entry[f"error_{level}x"] = (value - reference) / abs(reference)
            errors = [abs(entry[f"error_{level}x"]) for level in LEVELS]
            entry["error_monotone"] = errors[0] > errors[1] > errors[2]
            if entry["extrapolated"] is not None:
                entry["extrapolated_vs_featflow_l4"] = (entry["extrapolated"] - reference) / abs(reference)
            # Featflow level 3 -> 4 change at dt = 0.00025 s: a measure of the
            # reference's own discretisation error
            entry["featflow_l3_to_l4"] = (reference - entry["featflow_l3"]) / abs(reference)
            if column in runs[1]["smoothed"]:
                s1, s2, s4 = (runs[level]["smoothed"][column] for level in LEVELS)
                for level, value in zip(LEVELS, (s1, s2, s4)):
                    entry[f"smoothed_{level}x"] = value
                    entry[f"smoothed_error_{level}x"] = (value - reference) / abs(reference)
                # The smoothed values remove the step-to-step coupling noise,
                # which at dt = 0.00025 s inflates the extrema of the forces;
                # their threshold is the late-time variability of the smoothed
                # per-period values
                smoothed_noise = max(runs[2]["late_variability_smoothed"][column],
                                     runs[4]["late_variability_smoothed"][column])
                entry["smoothed_late_variability_4x"] = runs[4]["late_variability_smoothed"][column]
                smoothed = classify(column, s1, s2, s4, smoothed_noise)
                for key in ("d21", "d42", "ratio", "order", "extrapolated",
                            "extrapolation_meaningful", "status"):
                    entry[f"smoothed_{key}"] = smoothed[key]
            for name, run in diagnostics.items():
                entry[f"value_{name}"] = run["values"][column]
                entry[f"change_{name}_vs_2x"] = run["values"][column] - f2
                if column in run["smoothed"]:
                    entry[f"smoothed_{name}"] = run["smoothed"][column]
            table.append(entry)

    args.output.mkdir(parents=True, exist_ok=True)
    stem = f"fsi3_{args.coupling}_refinement_study"
    for run in list(runs.values()) + list(diagnostics.values()):
        write_history(run["case"], spec, args.history_start,
                      args.output / f"fsi3_{run['case']}_history.csv")
    level_rows = []
    for level in LEVELS + tuple(diagnostics):
        run = runs[level] if level in runs else diagnostics[level]
        level_rows.append({
            "level": f"{level}x" if level in runs else level,
            "end_time": run.get("end_time", end_time),
            "case": run["case"],
            "fluid_cells": run["fluid_cells"], "solid_cells": run["solid_cells"],
            "delta_t": run["delta_t"], "cores": run["cores"], "status": run["status"],
            "periods_in_window": run["periods_in_window"],
            "late_variability_periods": run["late_periods"],
            **{f"periodicity_{q}": run["amplitude_change"][q] for q in QUANTITIES},
            "coupled_steps": run["residuals"]["coupled"]["steps"],
            "mean_fsi_iterations": run["residuals"]["coupled"]["mean_iterations"],
            "max_fsi_iterations": run["residuals"]["coupled"]["max_iterations"],
            "max_final_residual": run["residuals"]["coupled"]["max_final_residual"],
            "window_mean_fsi_iterations": run["residuals"]["window"]["mean_iterations"],
            "window_max_final_residual": run["residuals"]["window"]["max_final_residual"],
            "steps_above_tolerance": run["residuals"]["coupled"]["steps_above_tolerance"],
            **{f"{q}_noise_rms": run["noise_rms"][q] for q in QUANTITIES},
            **{f"{q}_noise_max": run["noise_max"][q] for q in QUANTITIES},
            **{f"{column}_smoothed": value for column, value in run["smoothed"].items()},
            "outer_corr_tolerance": run["outer_corr_tolerance"],
            "execution_time_s": run["execution_time_s"], "core_hours": run["core_hours"],
            **{column: run["values"][column] for column in references},
            **{f"{column}_cycle_spread": run["cycle_spread"][column] for column in references},
            **{f"{column}_late_variability": run["late_variability"][column] for column in references},
        })

    with (args.output / f"{stem}_levels.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(level_rows[0]))
        writer.writeheader()
        writer.writerows(level_rows)
    with (args.output / f"{stem}_quantities.csv").open("w", newline="") as handle:
        columns = []
        for entry in table:
            columns += [key for key in entry if key not in columns]
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        writer.writerows(table)
    (args.output / f"{stem}.json").write_text(json.dumps({
        "description": (
            "FSI3 IQN-ILS study on three meshes with dt halved with h (a combined "
            "space-time refinement path). Values from the last full u_y period of "
            "the closing 1 s of a run to t = 7 s; cycle_spread is the range over all "
            "full periods of that window. Errors are signed, relative to the Featflow "
            "level-4 table at dt = 0.00025 s."),
        "levels": level_rows,
        "quantities": table,
        "cycles": {str(level): runs[level]["cycles"] for level in LEVELS},
    }, indent=2))

    def fmt(column: str, value: float | None) -> str:
        if value is None:
            return "-"
        if column.endswith("frequency"):
            return f"{value:.3f}"
        if column.startswith(("ux", "uy")):
            return f"{1000 * value:.3f}"
        return f"{value:.2f}"

    lines = ["| Quantity | 1x | 2x | 4x | Featflow L4 | Error 1x | Error 2x | Error 4x | "
             "2x-1x | 4x-2x | Ratio | Order | Status |",
             "|---|" + "---:|" * 11 + "---|"]
    for e in table:
        q = e["quantity"]
        errors = " | ".join("-" if e["near_zero"] else f"{100 * e[f'error_{level}x']:+.2f}%"
                            for level in LEVELS)
        ratio = "-" if e["ratio"] is None else f"{e['ratio']:.3f}"
        order = "-" if e["order"] is None else f"{e['order']:.2f}"
        lines.append(
            f"| {q} | {fmt(q, e['value_1x'])} | {fmt(q, e['value_2x'])} | {fmt(q, e['value_4x'])} | "
            f"{fmt(q, e['featflow_l4'])} | {errors} | {fmt(q, e['d21'])} | {fmt(q, e['d42'])} | "
            f"{ratio} | {order} | {e['status']} |")
    (args.output / f"{stem}.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))
    for row in level_rows:
        print({k: row[k] for k in ("level", "status", "cores", "execution_time_s",
                                   "mean_fsi_iterations", "max_fsi_iterations",
                                   "max_final_residual", "steps_above_tolerance")})
    return 0


if __name__ == "__main__":
    sys.exit(main())
