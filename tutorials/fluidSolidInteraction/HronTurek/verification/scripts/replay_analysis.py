#!/usr/bin/env python3
"""Analyse the replay runs of hron_turek_replay.py.

The replayed coupled trajectory is periodic at the fitted frequency f. Over the
last five periods of each run the script evaluates mean, first and second
harmonics of the force on the flag ("plate") and on cylinder + flag ("total"),
per unit depth; the first harmonic is split into the part in phase with the
tip displacement (`inphase`) and in phase with the tip velocity (`quad`); and
the work done by the pressure on the flag per cycle (`POWR` log lines,
integrated over a period). The 2x run is also compared with the coupled
source run over its own fit window (total force only), which validates the
replay. Paths as in cfd3_analysis.py. Usage: python3 scripts/replay_analysis.py
"""

from __future__ import annotations

import csv
import json
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import hron_turek_verification as driver  # noqa: E402
import hron_turek_ale as ale  # noqa: E402
import ale_analysis as aa  # noqa: E402
import cfd3_analysis as cfd3  # noqa: E402

PERIODS = 5
PATHS = {
    "matched": [(1, 1e-3), (2, 5e-4), (4, 2.5e-4)],
    "fixed_dt": [(1, 2.5e-4), (2, 2.5e-4), (4, 2.5e-4)],
}


def tip_fundamental(case: Path, omega: float) -> tuple[float, float]:
    """(sin, cos) coefficients of the replayed tip y-displacement."""
    lines = (case / "constant/plateMotion.tab").read_text().splitlines()
    best = None
    for line in lines[1:]:
        v = line.split()
        x, y = int(v[0]) / 1e5, int(v[1]) / 1e5
        d = math.hypot(x - 0.6, y - 0.19)
        if best is None or d < best[0]:
            best = (d, [float(c) for c in v[2:]])
    coef = best[1]
    ncoef = len(coef) // 2
    return coef[ncoef + 2], coef[ncoef + 1]  # sin, cos of dy


def harmonics(t: list[float], y: list[float], end: float, f: float,
              tip: tuple[float, float]) -> dict:
    aa.F = f
    aa.T0 = ale.START
    start = end - PERIODS / f
    w = aa.window_stats(t, y, start, PERIODS)
    fs, fc = w["a1"], w["b1"]          # sin, cos components of the force
    ds, dc = tip
    norm = math.hypot(ds, dc)
    return {"mean": w["mean"], "amp1": w["amp1"], "amp2": w["amp2"],
            "inphase": (fs * ds + fc * dc) / norm, "quad": (fc * ds - fs * dc) / norm,
            "unlocked": w["unlocked_fraction"]}


def analyse_run(level: int, delta_t: float, tag: str = "_replay") -> dict | None:
    case = driver.WORK_ROOT / ale.case_name(level, delta_t, tag)
    log = case / "log.solids4Foam"
    if not log.is_file() or not driver.re.search(r"^End\s*$", log.read_text(errors="replace"),
                                                 driver.re.MULTILINE):
        return None
    fit = json.loads((case / "replay_fit.json").read_text())
    f = fit["frequency"]
    tip = tip_fundamental(case, fit["omega"])
    out: dict = {"level": level, "delta_t": delta_t, "frequency_hz": f,
                 "runtime_s": driver.execution_time(case)}
    out["cores"] = len(list(case.glob("processor*"))) or 1
    out["core_hours"] = out["cores"] * out["runtime_s"] / 3600
    for group, function in (("total", "forces"), ("plate", "forcesPlate")):
        c = aa.read_forces(case, function)
        end = c["t"][-1]
        for name in ("drag", "lift", "drag_p", "lift_p"):
            h = harmonics(c["t"], c[name], end, f, tip)
            for q, v in h.items():
                out[f"{group}_{name}_{q}"] = v
            # window-to-window change of the fundamental
            prev = aa.window_stats(c["t"], c[name], end - 2 * PERIODS / f, PERIODS)
            out[f"{group}_{name}_amp1_change_prev"] = h["amp1"] - prev["amp1"]
            out[f"{group}_{name}_mean_change_prev"] = h["mean"] - prev["mean"]
    power: dict[float, float] = {}
    for line in log.read_text(errors="replace").splitlines():
        if line.startswith("POWR "):
            _, tt, v = line.split()
            power[round(float(tt), 6)] = float(v)
    times = sorted(power)
    end = times[-1]
    y = [power[t] for t in times]
    start = end - PERIODS / f
    idx = [i for i, t in enumerate(times) if start - 1e-12 <= t <= end + 1e-12]
    mean_power = sum(0.5 * (y[i] + y[j]) * (times[j] - times[i])
                     for i, j in zip(idx[:-1], idx[1:])) / (times[idx[-1]] - times[idx[0]])
    out["work_per_cycle"] = mean_power / f
    prev_idx = [i for i, t in enumerate(times) if start - PERIODS / f - 1e-12 <= t <= start + 1e-12]
    prev_power = sum(0.5 * (y[i] + y[j]) * (times[j] - times[i])
                     for i, j in zip(prev_idx[:-1], prev_idx[1:])) / (times[prev_idx[-1]] - times[prev_idx[0]])
    out["work_per_cycle_change_prev"] = (mean_power - prev_power) / f
    return out


def source_reference() -> dict:
    """Coupled source run: total force over its fit window (5 periods)."""
    case = driver.WORK_ROOT / "replaysrc_iqnils_mesh_2x"
    fit = json.loads((driver.WORK_ROOT / "ale_2x_dt0.0005_replay/replay_fit.json").read_text())
    f = fit["frequency"]
    ta = float(driver.re.findall(r"[-+]?\d+\.\d+(?:[eE][-+]?\d+)?", str(fit["window"]))[0])
    c = aa.read_forces(case, "fluid/forces")
    out = {}
    for name in ("drag", "lift"):
        aa.F = f
        aa.T0 = ale.START
        w = aa.window_stats(c["t"], c[name], ta, PERIODS)
        out[f"total_{name}_mean"] = w["mean"]
        out[f"total_{name}_amp1"] = w["amp1"]
        out[f"total_{name}_amp2"] = w["amp2"]
    return out


COLUMNS = ["work_per_cycle", "plate_lift_amp1", "plate_lift_inphase", "plate_lift_quad",
           "plate_drag_mean", "plate_drag_amp1", "plate_drag_amp2", "plate_lift_mean",
           "total_lift_amp1", "total_lift_inphase", "total_lift_quad", "total_lift_mean",
           "total_drag_mean", "total_drag_amp2", "total_drag_amp1",
           "plate_lift_p_amp1", "plate_drag_p_mean"]


def noise_of(runs: list[dict], column: str) -> float:
    if column == "work_per_cycle":
        return max(abs(r["work_per_cycle_change_prev"]) for r in runs)
    tokens = column.split("_")
    base = "_".join(tokens[:3] if tokens[2] == "p" else tokens[:2])
    kind = "mean" if column.endswith("_mean") else "amp1"
    return max(abs(r[f"{base}_{kind}_change_prev"]) for r in runs)


def main() -> int:
    results: dict = {"runs": {}, "paths": {}, "source_coupled_2x": None}
    for path_name, spec in PATHS.items():
        runs = []
        for level, dt in spec:
            key = ale.case_name(level, dt, "_replay")
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
            entry = cfd3.successive(values, noise_of(runs, column), column)
            entry["values"] = values
            entry["noise"] = noise_of(runs, column)
            entry["relative_change_pct"] = [100 * entry["d21"] / abs(values[0]) if values[0] else None,
                                            100 * entry["d32"] / abs(values[1]) if values[1] else None]
            table[column] = entry
        results["paths"][path_name] = {"cases": [ale.case_name(*s, "_replay") for s in spec],
                                       "quantities": table}
    try:
        results["source_coupled_2x"] = source_reference()
    except Exception as exc:  # noqa: BLE001
        results["source_coupled_2x"] = f"unavailable: {exc}"
    out = driver.OUTPUT_ROOT
    out.mkdir(exist_ok=True)
    (out / "replay_study.json").write_text(json.dumps(results, indent=2))
    rows = list(results["runs"].values())
    if rows:
        with (out / "replay_study_runs.csv").open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
    with (out / "replay_study_paths.csv").open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["path", "quantity", "v1", "v2", "v3", "d21", "d32", "ratio", "order",
                         "noise", "status"])
        for p, body in results["paths"].items():
            for q, e in body["quantities"].items():
                writer.writerow([p, q, *e["values"], e["d21"], e["d32"], e["ratio"],
                                 e["order"], e["noise"], e["status"]])
    print("coupled source (2x):", results["source_coupled_2x"])
    for p, body in results["paths"].items():
        print(f"\n== {p}: {', '.join(body['cases'])}")
        for q, e in body["quantities"].items():
            v = " ".join(f"{x:11.5g}" for x in e["values"])
            o = f"{e['order']:.2f}" if e["order"] is not None else "  - "
            r = f"{e['ratio']:.3f}" if e["ratio"] is not None else "  -  "
            print(f"{q:24s} {v}  ratio {r} order {o} noise {e['noise']:.3g} {e['status']}")
    for k, r in results["runs"].items():
        print(k, f"core-h {r['core_hours']:.1f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
