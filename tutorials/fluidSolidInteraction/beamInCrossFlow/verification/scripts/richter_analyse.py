#!/usr/bin/env python3
"""Summarise campaign cases: quantities of interest, steadiness, coupling.

For every case directory given, the script reads

- the point-A displacement history (solidPointDisplacement),
- the interface force history (forces function object),
- the coupling iterations per step (postProcessing/fsiResiduals.dat), and
- optionally the solver logs, for the two residual measures r1 and r2 at
  the accepted iteration of every step,

and writes one JSON summary per case (``campaign_summary.json`` in the case)
plus one CSV row per case to ``--csv``.

Steadiness is reported as the relative change of each quantity over the last
``--window`` seconds, and as the change over the last 1 s extrapolated with
an exponential fit (``u(t) = u_inf + a exp(-t/tau)``) over the last
``--fit-window`` seconds, where the fit is well conditioned.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
from pathlib import Path


def numeric_rows(path: Path) -> list[list[float]]:
    rows = []
    for line in path.read_text(errors="replace").splitlines():
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        fields = line.replace("(", " ").replace(")", " ").split()
        try:
            rows.append([float(value) for value in fields])
        except ValueError:
            continue
    return rows


def find_one(case: Path, pattern: str) -> Path | None:
    hits = sorted(case.glob(pattern))
    return hits[0] if hits else None


def histories(case: Path) -> dict[str, list[tuple[float, float]]]:
    result: dict[str, list[tuple[float, float]]] = {}
    paths = sorted(case.glob("postProcessing/**/solidPointDisplacement*_displacement.dat"))
    rows: dict[float, list[float]] = {}
    for path in paths:  # restarts write new time directories
        for row in numeric_rows(path):
            rows[row[0]] = row
    if rows:
        times = sorted(rows)
        for index, name in enumerate(("ux", "uy", "uz"), start=1):
            result[name] = [(t, rows[t][index]) for t in times]
    rows = {}
    for path in sorted(case.glob("postProcessing/**/force.dat")):
        for row in numeric_rows(path):
            rows[row[0]] = row
    if rows:
        times = sorted(rows)
        for index, name in enumerate(("fx", "fy", "fz"), start=1):
            result[name] = [(t, rows[t][index]) for t in times]
    return result


def value_at(series: list[tuple[float, float]], time: float) -> float:
    best = None
    for t, value in series:
        if t <= time + 1e-9:
            best = value
        else:
            break
    if best is None:
        raise ValueError(f"no value at or before t={time}")
    return best


def exponential_fit(series: list[tuple[float, float]], start: float):
    """Fit u = u_inf + a exp(-t/tau) by a scan over tau (robust, no scipy)."""
    data = [(t, u) for t, u in series if t >= start]
    if len(data) < 20:
        return None
    best = None
    for k in range(1, 400):
        tau = 0.02 * 1.03 ** k
        xs = [math.exp(-(t - data[-1][0]) / tau) for t, _ in data]
        ys = [u for _, u in data]
        n = len(xs)
        sx, sy = sum(xs), sum(ys)
        sxx = sum(x * x for x in xs)
        sxy = sum(x * y for x, y in zip(xs, ys))
        det = n * sxx - sx * sx
        if abs(det) < 1e-30:
            continue
        a = (n * sxy - sx * sy) / det
        u_inf = (sy - a * sx) / n
        sse = sum((u_inf + a * x - y) ** 2 for x, y in zip(xs, ys))
        if best is None or sse < best[0]:
            best = (sse, tau, u_inf, a)
    if best is None:
        return None
    sse, tau, u_inf, a = best
    return {"tau": tau, "u_inf": u_inf, "offset_at_end": a, "rms": math.sqrt(sse / len(data))}


def coupling_from_residual_file(case: Path) -> dict:
    path = case / "postProcessing/fsiResiduals.dat"
    if not path.is_file():
        return {}
    last: dict[float, tuple[int, float]] = {}
    for row in numeric_rows(path):
        last[row[0]] = (int(row[1]), row[2])
    return {t: v for t, v in last.items()}


def coupling_from_logs(case: Path) -> dict[float, dict]:
    """Return {time: {iterations, r1, r2}} at the accepted iteration."""
    steps: dict[float, dict] = {}
    logs = sorted(case.glob("log.solids4Foam*"), key=lambda p: (p.name != "log.solids4Foam", p.name))
    time_re = re.compile(r"^Time = ([-+0-9.eE]+)$")
    for log in logs:
        current = None
        r1 = r2 = None
        count = 0
        with log.open(errors="replace") as handle:
            for line in handle:
                if line.startswith("Time = "):
                    match = time_re.match(line.strip())
                    if match:
                        if current is not None and r1 is not None:
                            steps[current] = {"iterations": count, "r1": r1, "r2": r2}
                        current = float(match.group(1))
                        r1 = r2 = None
                        count = 0
                    continue
                if line.startswith("FSI relative residual1 norm"):
                    r1 = float(line.split()[-1])
                    count += 1
                elif line.startswith("FSI residual2 norm"):
                    r2 = float(line.split()[-1])
        if current is not None and r1 is not None:
            steps[current] = {"iterations": count, "r1": r1, "r2": r2}
    return steps


def summarise(case: Path, args: argparse.Namespace) -> dict:
    meta_path = case / "campaign_case.json"
    meta = json.loads(meta_path.read_text()) if meta_path.is_file() else {}
    series = histories(case)
    summary: dict = {"case": case.name}
    for key in ("factor", "thickness_multiplier", "dt", "coupling_tol", "strict",
                "n_outer_corr", "pimple_outer", "p_tol", "u_tol", "snes_rtol",
                "fluid_cells", "solid_cells", "rigid", "fluid_ddt", "solid_d2dt2"):
        if key in meta:
            summary[key] = meta[key]
    if not series:
        return summary
    any_series = series.get("fx") or series.get("ux")
    end = any_series[-1][0]
    summary["end_time"] = end
    for name, data in series.items():
        summary[name] = data[-1][1]
        for window in (1.0, 2.0):
            if end - window >= data[0][0]:
                start_value = value_at(data, end - window)
                scale = abs(data[-1][1]) if abs(data[-1][1]) > 0 else 1.0
                summary[f"{name}_change_last_{window:g}s"] = (data[-1][1] - start_value) / scale
        if name in ("ux", "fx") and end - args.fit_window > 1.5:
            fit = exponential_fit(data, end - args.fit_window)
            if fit:
                scale = abs(data[-1][1])
                summary[f"{name}_fit_tau"] = fit["tau"]
                summary[f"{name}_fit_inf"] = fit["u_inf"]
                summary[f"{name}_fit_remaining_rel"] = -fit["offset_at_end"] / scale
                summary[f"{name}_fit_rms_rel"] = fit["rms"] / scale
    residual_steps = coupling_from_residual_file(case)
    if residual_steps:
        iterations = [v[0] for v in residual_steps.values()]
        summary["steps"] = len(iterations)
        summary["coupling_iterations_total"] = sum(iterations)
        summary["coupling_iterations_mean"] = sum(iterations) / len(iterations)
        summary["coupling_iterations_max"] = max(iterations)
        cap = meta.get("n_outer_corr")
        if cap:
            summary["capped_steps"] = sum(1 for v in iterations if v >= cap)
        late = [v[0] for t, v in residual_steps.items() if t > end - 1.0]
        if late:
            summary["coupling_iterations_mean_last_1s"] = sum(late) / len(late)
    if args.logs:
        steps = coupling_from_logs(case)
        if steps:
            tol = meta.get("coupling_tol")
            r1 = [s["r1"] for s in steps.values()]
            r2 = [s["r2"] for s in steps.values() if s["r2"] is not None]
            summary["accepted_r1_max"] = max(r1)
            summary["accepted_r2_max"] = max(r2) if r2 else None
            late = [s for t, s in steps.items() if t > end - 1.0]
            if late:
                summary["accepted_r1_max_last_1s"] = max(s["r1"] for s in late)
                summary["accepted_r2_max_last_1s"] = max(s["r2"] for s in late)
            if tol:
                summary["steps_accepted_r1_above_tol"] = sum(1 for v in r1 if v > tol)
                summary["steps_accepted_r2_above_tol"] = sum(1 for v in r2 if v > tol)
            (case / "campaign_steps.csv").write_text(
                "time,iterations,r1,r2\n" + "".join(
                    f"{t:.10g},{s['iterations']},{s['r1']:.6g},{s['r2']:.6g}\n"
                    for t, s in sorted(steps.items())
                )
            )
    (case / "campaign_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    return summary


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("cases", nargs="+")
    parser.add_argument("--csv")
    parser.add_argument("--logs", action="store_true",
                        help="parse solver logs for r1/r2 at acceptance")
    parser.add_argument("--fit-window", type=float, default=3.0)
    args = parser.parse_args()
    rows = []
    for case in args.cases:
        rows.append(summarise(Path(case), args))
    keys: list[str] = []
    for row in rows:
        for key in row:
            if key not in keys:
                keys.append(key)
    if args.csv:
        with open(args.csv, "w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=keys)
            writer.writeheader()
            writer.writerows(rows)
    for row in rows:
        print(json.dumps(row))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
