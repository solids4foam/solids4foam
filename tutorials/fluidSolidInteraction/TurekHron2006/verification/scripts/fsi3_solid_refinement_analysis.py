#!/usr/bin/env python3
"""Error decomposition of the FSI3 mesh dependence: solid refined alone.

Reads the completed runs of the fixed-fluid study (fluid 2x; solid 2x, 4x, 8x),

    ./Allverify --levels 2 --cores N                       # solid 2x (baseline)
    ./Allverify --levels 2 --solid-refinement 4 --cores N
    ./Allverify --levels 2 --solid-refinement 8 --cores N

and writes, to postProcessing/fsi3_solid_refinement_study.{json,csv,md}, the
benchmark quantities of each run, the successive changes with the solid
refinement, the interface face counts, the force balance across the interface
and the coupling iterations. It also compares the changes with those of the
matched 1x/2x/4x sequence (runs iqnils_mesh_{1,2,4}x, if present).

The fluid mesh is fixed, so this is an error decomposition, not a spatial
order study: no convergence order is computed.

Usage, from the verification directory:

    python3 scripts/fsi3_solid_refinement_analysis.py
"""

from __future__ import annotations

import csv
import json
import math
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import hron_turek_verification as driver  # noqa: E402
import fsi3_refinement_analysis as base  # noqa: E402

RUNS = [(2, 2, "iqnils_mesh_2x"),
        (2, 4, "iqnils_mesh_2x_solid4x"),
        (2, 8, "iqnils_mesh_2x_solid8x")]
MATCHED = [(1, "iqnils_mesh_1x"), (2, "iqnils_mesh_2x"), (4, "iqnils_mesh_4x")]
QUANTITIES = ("ux_mean", "ux_amplitude", "ux_frequency", "uy_amplitude",
              "uy_frequency", "drag_mean", "drag_amplitude", "lift_amplitude")
# Smoothed (4 ms moving average) amplitudes where the extrema are noise-inflated
SMOOTHED = ("drag_amplitude", "lift_amplitude")


def patch_faces(case: Path, region: str, patch: str) -> int | None:
    path = case / "constant" / region / "polyMesh" / "boundary"
    if not path.is_file():
        return None
    match = re.search(rf"\b{patch}\s*\{{[^}}]*?nFaces\s+(\d+);", path.read_text())
    return int(match.group(1)) if match else None


def force_balance(case: Path, start: float) -> dict:
    """rms of (fluid + solid interface force) over rms fluid force, from start.

    The AMI transfers the traction fluid -> solid; the two interface totals
    are printed at every coupling iteration, so the ratio measures the
    conservation of the force transfer (not the iteration residual).
    """
    time = 0.0
    fluid = None
    num = den = 0.0
    n = 0
    vec = r"\(([-+0-9.eE]+) ([-+0-9.eE]+) ([-+0-9.eE]+)\)"
    with (case / "log.solids4Foam").open(errors="replace") as handle:
        for line in handle:
            if line.startswith("Time = "):
                time = float(line.split()[2].rstrip(","))
            elif time >= start and line.startswith("Total force on fluid interface"):
                fluid = [float(x) for x in re.search(vec, line).groups()]
            elif time >= start and fluid and line.startswith("Total force on solid interface"):
                solid = [float(x) for x in re.search(vec, line).groups()]
                num += sum((a + b) ** 2 for a, b in zip(fluid, solid))
                den += sum(a * a for a in fluid)
                n += 1
                fluid = None
    return {"force_balance_rms_ratio": math.sqrt(num / den) if den else float("nan"),
            "force_balance_samples": n}


def analyse(fluid: int, solid: int, name: str, spec: dict, window: float,
            end_time: float) -> dict:
    case = driver.WORK_ROOT / name
    result = base.analyse_case(name, spec, window, end_time)
    values = dict(result["values"])
    for quantity in SMOOTHED:
        values[f"{quantity}_smoothed"] = result["smoothed"][quantity]
    row = {"case": name, "fluid_level": fluid, "solid_level": solid,
           "fluid_cells": result["fluid_cells"], "solid_cells": result["solid_cells"],
           "dt": result["delta_t"], "coupling_tolerance": result["outer_corr_tolerance"],
           "fluid_interface_faces": patch_faces(case, "fluid", "plate"),
           "solid_interface_faces": patch_faces(case, "solid", "plate"),
           "cores": result["cores"], "runtime_s": result["execution_time_s"],
           "core_hours": result["core_hours"],
           "mean_outer_iterations": result["residuals"]["coupled"]["mean_iterations"],
           "max_outer_iterations": result["residuals"]["coupled"]["max_iterations"],
           "max_final_residual": result["residuals"]["coupled"]["max_final_residual"],
           "steps_above_tolerance": result["residuals"]["coupled"]["steps_above_tolerance"],
           "status": result["status"]}
    row.update(force_balance(case, end_time - window))
    row.update(values)
    row["_spread"] = result["cycle_spread"]
    row["_late_variability"] = result["late_variability"]
    return row


def change(a: float, b: float) -> tuple[float, float]:
    return b - a, (b - a) / a if a else float("nan")


def main() -> int:
    spec = json.loads(driver.REFERENCE_FILE.read_text())
    window = spec["mesh"]["analysisWindow"]
    end_time = spec["mesh"]["endTime"]
    rows = []
    for fluid, solid, name in RUNS:
        log = driver.WORK_ROOT / name / "log.solids4Foam"
        if log.is_file() and re.search(r"^End\s*$", log.read_text(errors="replace"), re.M):
            rows.append(analyse(fluid, solid, name, spec, window, end_time))
    columns = QUANTITIES + tuple(f"{q}_smoothed" for q in SMOOTHED)
    for previous, row in zip([None] + rows[:-1], rows):
        for column in columns:
            if previous is None:
                continue
            row[f"{column}_change"], row[f"{column}_change_rel"] = change(
                previous[column], row[column])
    # Matched 1x/2x/4x sequence and the structural CSM2 scaling, for comparison
    matched = {}
    for factor, name in MATCHED:
        if (driver.WORK_ROOT / name / "log.solids4Foam").is_file():
            matched[factor] = base.analyse_case(name, spec, window, end_time)["values"]
    out = driver.OUTPUT_ROOT
    flat = [{k: v for k, v in r.items() if not k.startswith("_")} for r in rows]
    driver.write_csv(out / "fsi3_solid_refinement_study.csv", flat)
    refs = {q: spec["references"][q]["value"] for q in QUANTITIES if q in spec["references"]}
    payload = {"fluid_level": 2, "runs": flat,
               "matched_sequence": {str(k): v for k, v in matched.items()},
               "cycle_to_cycle_spread": {r["case"]: r["_spread"] for r in rows},
               "late_variability": {r["case"]: r["_late_variability"] for r in rows},
               "featflow_L4": refs}
    (out / "fsi3_solid_refinement_study.json").write_text(json.dumps(payload, indent=1))

    unit = {"ux": 1e3, "uy": 1e3}
    def fmt(column: str, value: float) -> str:
        q = column.split("_")[0]
        scale = 1 if "frequency" in column else unit.get(q, 1)
        return f"{value * scale:.5g}"
    lines = ["| Quantity | " + " | ".join(f"f2/s{x}x" for x in (2, 4, 8))
             + " | s2→s4 | s4→s8 | matched 1x→2x | matched 2x→4x |",
             "|---|" + "---:|" * 7]
    for column in columns:
        cells = [fmt(column, r[column]) for r in rows] + [""] * (3 - len(rows))
        for i in (1, 2):
            if len(rows) > i:
                r = rows[i]
                cells.append(f"{r[column + '_change_rel'] * 100:+.2f}%")
            else:
                cells.append("")
        for a, b in ((1, 2), (2, 4)):
            key = column.removesuffix("_smoothed")
            if a in matched and b in matched and column in matched[a]:
                cells.append(f"{change(matched[a][column], matched[b][column])[1] * 100:+.2f}%")
            else:
                cells.append("")
        lines.append(f"| {column} | " + " | ".join(cells) + " |")
    (out / "fsi3_solid_refinement_study.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))
    for r in flat:
        print({k: r[k] for k in ("case", "fluid_cells", "solid_cells", "fluid_interface_faces",
                                 "solid_interface_faces", "mean_outer_iterations",
                                 "max_outer_iterations", "force_balance_rms_ratio",
                                 "runtime_s", "cores", "status")})
    return 0


if __name__ == "__main__":
    sys.exit(main())
