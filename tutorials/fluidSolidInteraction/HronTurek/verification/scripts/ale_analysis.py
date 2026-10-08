#!/usr/bin/env python3
"""Analyse the prescribed-motion ALE runs of hron_turek_ale.py.

The flag bends with the prescribed first-mode motion at f = 5.5 Hz, so the
forces are analysed as periodic signals at the forcing frequency over the last
`PERIODS` forcing periods of the run: mean, extrema-based amplitude over the
last period, and the Fourier components at f (and 2f for the drag) split into
the part in phase with the prescribed displacement (a) and in quadrature (b,
in phase with the velocity). The quadrature part is the one that transfers
energy between fluid and flag. Forces are on the cylinder and flag together
("total") and on the flag alone ("plate"), per unit depth. A periodicity
check compares the fundamental of the last window with the preceding one, and
the fraction of the variance not captured by the first six harmonics shows
whether the flow has locked to the forcing.

Paths (as in cfd3_analysis.py): matched 1x/2x/4x with dt = 1e-3/5e-4/2.5e-4 and
fixed dt = 2.5e-4. Usage: python3 scripts/ale_analysis.py
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
import cfd3_analysis as cfd3  # noqa: E402

THICKNESS = 0.015
F = ale.FREQUENCY
T0 = ale.START
PERIODS = 5
# Generalised (modal) structural inertia force of the FSI3 flag, per unit depth,
# for the prescribed motion: rho_s t L int(phi^2) w^2 A, with int(phi^2) = 1/4 for
# a clamped-free mode normalised to unit tip value. It sets the scale against
# which differences of the fluid generalised force matter to the flag motion.
RHO_S, THICKNESS_FLAG = 1000.0, 0.02
STRUCTURAL_SCALE = (RHO_S * THICKNESS_FLAG * (ale.X_TIP - ale.X_ROOT) * 0.25
                    * (2 * math.pi * F) ** 2 * ale.AMPLITUDE)
PATHS = {
    "matched": [(1, 1e-3), (2, 5e-4), (4, 2.5e-4)],
    "fixed_dt": [(1, 2.5e-4), (2, 2.5e-4), (4, 2.5e-4)],
}


def read_forces(case: Path, function: str) -> dict[str, list[float]]:
    path = sorted(case.glob(f"postProcessing/{function}/*/force.dat"))
    # Restarts write further directories; concatenate by time
    cols = {k: [] for k in ("t", "drag", "lift", "drag_p", "lift_p", "drag_v", "lift_v")}
    seen = set()
    for p in path:
        for f in driver.numeric_rows(p):
            v = [float(x) for x in f]
            if len(v) < 10 or v[0] in seen:
                continue
            seen.add(v[0])
            cols["t"].append(v[0])
            for name, i in (("drag", 1), ("lift", 2), ("drag_p", 4), ("lift_p", 5),
                            ("drag_v", 7), ("lift_v", 8)):
                cols[name].append(v[i] / THICKNESS)
    return cols


def fourier(t: list[float], y: list[float], k: int, t_start: float, n: int) -> tuple[float, float]:
    """Components a, b of y ~ a sin(k w tau) + b cos(k w tau), tau = t - T0,
    over n forcing periods starting at t_start (trapezoidal rule)."""
    w = 2 * math.pi * F * k
    length = n / F
    a = b = 0.0
    idx = [i for i, x in enumerate(t) if t_start - 1e-12 <= x <= t_start + length + 1e-12]
    for i, j in zip(idx[:-1], idx[1:]):
        dt = t[j] - t[i]
        for (ti, yi), c in (((t[i], y[i]), 0.5), ((t[j], y[j]), 0.5)):
            a += c * dt * yi * math.sin(w * (ti - T0))
            b += c * dt * yi * math.cos(w * (ti - T0))
    return 2 * a / length, 2 * b / length


def window_stats(t: list[float], y: list[float], t_start: float, n: int) -> dict:
    length = n / F
    idx = [i for i, x in enumerate(t) if t_start - 1e-12 <= x <= t_start + length + 1e-12]
    seg = [y[i] for i in idx]
    mean = 0.0
    for i, j in zip(idx[:-1], idx[1:]):
        mean += 0.5 * (y[i] + y[j]) * (t[j] - t[i])
    mean /= length
    a1, b1 = fourier(t, y, 1, t_start, n)
    a2, b2 = fourier(t, y, 2, t_start, n)
    # variance outside the first six harmonics
    rec_var = sum(0.5 * (math.hypot(*fourier(t, y, k, t_start, n)) ** 2) for k in range(1, 7))
    var = sum((yi - mean) ** 2 for yi in seg) / len(seg)
    return {"mean": mean, "a1": a1, "b1": b1, "amp1": math.hypot(a1, b1),
            "a2": a2, "b2": b2, "amp2": math.hypot(a2, b2),
            "unlocked_fraction": max(0.0, 1 - rec_var / var) if var > 0 else 0.0}


def analyse_run(level: int, delta_t: float) -> dict | None:
    case = driver.WORK_ROOT / ale.case_name(level, delta_t)
    log = case / "log.solids4Foam"
    if not log.is_file() or not driver.re.search(r"^End\s*$", log.read_text(errors="replace"),
                                                 driver.re.MULTILINE):
        return None
    out: dict = {"level": level, "delta_t": delta_t}
    out["runtime_s"] = driver.execution_time(case)
    out["cores"] = len(list(case.glob("processor*"))) or 1
    out["core_hours"] = out["cores"] * out["runtime_s"] / 3600
    for group, function in (("total", "forces"), ("plate", "forcesPlate")):
        c = read_forces(case, function)
        t = c["t"]
        end = t[-1]
        s_last = end - PERIODS / F
        s_prev = s_last - PERIODS / F
        out["end_time"] = end
        for name in ("drag", "lift", "drag_p", "lift_p", "drag_v", "lift_v"):
            last = window_stats(t, c[name], s_last, PERIODS)
            prev = window_stats(t, c[name], s_prev, PERIODS)
            key = f"{group}_{name}"
            for q in ("mean", "a1", "b1", "amp1", "amp2", "unlocked_fraction"):
                out[f"{key}_{q}"] = last[q]
            out[f"{key}_amp1_change_prev"] = last["amp1"] - prev["amp1"]
            out[f"{key}_mean_change_prev"] = last["mean"] - prev["mean"]
            # extrema-based over the last forcing period
            seg = [y for x, y in zip(t, c[name]) if end - 1 / F <= x <= end]
            out[f"{key}_mean_extrema"] = 0.5 * (max(seg) + min(seg))
            out[f"{key}_amplitude_extrema"] = 0.5 * (max(seg) - min(seg))
    # Pressure force on the flag projected on the prescribed mode shape
    # (generalised force) and its uniform sum, from the GENF log lines
    gen: dict[float, tuple[float, float]] = {}
    for line in log.read_text(errors="replace").splitlines():
        if line.startswith("GENF "):
            _, tt, q, u = line.split()
            gen[float(tt)] = (float(q), float(u))
    times = sorted(gen)
    end = times[-1]
    s_last, s_prev = end - PERIODS / F, end - 2 * PERIODS / F
    for index, name in enumerate(("q", "u")):
        y = [gen[x][index] for x in times]
        last = window_stats(times, y, s_last, PERIODS)
        prev = window_stats(times, y, s_prev, PERIODS)
        key = f"gen_{name}"
        for q in ("mean", "a1", "b1", "amp1", "amp2", "unlocked_fraction"):
            out[f"{key}_{q}"] = last[q]
        out[f"{key}_amp1_change_prev"] = last["amp1"] - prev["amp1"]
        out[f"{key}_mean_change_prev"] = last["mean"] - prev["mean"]
    # Net work per cycle done by the fluid on the prescribed motion (J/m):
    # integral of Q dq/dt over a period, q = A sin(w tau): pi A b1
    out["gen_q_work_per_cycle"] = math.pi * ale.AMPLITUDE * out["gen_q_b1"]
    return out


COLUMNS = ["gen_q_a1", "gen_q_b1", "gen_q_amp1", "gen_q_mean", "gen_q_work_per_cycle",
           "total_drag_mean", "total_drag_amp2", "total_lift_amp1", "total_lift_a1",
           "total_lift_b1", "plate_drag_mean", "plate_lift_amp1", "plate_lift_a1",
           "plate_lift_b1", "plate_lift_amplitude_extrema", "total_lift_amplitude_extrema",
           "plate_lift_p_amp1", "plate_lift_v_amp1", "plate_drag_p_mean", "plate_drag_v_mean"]


def main() -> int:
    results: dict = {"structural_scale_N_per_m": STRUCTURAL_SCALE, "forcing": {"amplitude_m": ale.AMPLITUDE, "frequency_hz": F,
                                 "ramp_s": ale.RAMP, "start_s": T0}, "runs": {}, "paths": {}}
    for path_name, spec in PATHS.items():
        runs = []
        for level, dt in spec:
            key = ale.case_name(level, dt)
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
            # Noise: change of the fundamental amplitude (or of the mean for
            # mean quantities) between the last two windows
            tokens = column.split("_")
            base = "_".join(tokens[:3] if tokens[2] in ("p", "v") else tokens[:2])
            kind = "mean" if column.endswith("_mean") else "amp1"
            noise = max(abs(r[f"{base}_{kind}_change_prev"]) for r in runs)
            entry = cfd3.successive(values, noise, column)
            entry["values"] = values
            entry["noise"] = noise
            entry["relative_change_pct"] = [100 * entry["d21"] / abs(values[0]) if values[0] else None,
                                            100 * entry["d32"] / abs(values[1]) if values[1] else None]
            table[column] = entry
        for column in ("gen_q_a1", "gen_q_b1"):
            e = table[column]
            e["d21_over_structural_scale"] = e["d21"] / STRUCTURAL_SCALE
            e["d32_over_structural_scale"] = e["d32"] / STRUCTURAL_SCALE
        results["paths"][path_name] = {"cases": [ale.case_name(*s) for s in spec],
                                       "quantities": table}
    out = driver.OUTPUT_ROOT
    out.mkdir(exist_ok=True)
    (out / "ale_study.json").write_text(json.dumps(results, indent=2))
    rows = list(results["runs"].values())
    if rows:
        with (out / "ale_study_runs.csv").open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
    with (out / "ale_study_paths.csv").open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["path", "quantity", "v1", "v2", "v3", "d21", "d32", "ratio", "order",
                         "noise", "status"])
        for p, body in results["paths"].items():
            for q, e in body["quantities"].items():
                writer.writerow([p, q, *e["values"], e["d21"], e["d32"], e["ratio"],
                                 e["order"], e["noise"], e["status"]])
    for p, body in results["paths"].items():
        print(f"\n== {p}: {', '.join(body['cases'])}")
        for q, e in body["quantities"].items():
            v = " ".join(f"{x:11.5g}" for x in e["values"])
            o = f"{e['order']:.2f}" if e["order"] is not None else "  - "
            r = f"{e['ratio']:.3f}" if e["ratio"] is not None else "  -  "
            print(f"{q:30s} {v}  ratio {r} order {o}  {e['status']}")
    for k, r in results["runs"].items():
        print(k, f"core-h {r['core_hours']:.1f} unlocked(total lift) {r['total_lift_unlocked_fraction']:.4f} "
              f"amp1 change prev {r['total_lift_amp1_change_prev']:.3f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
