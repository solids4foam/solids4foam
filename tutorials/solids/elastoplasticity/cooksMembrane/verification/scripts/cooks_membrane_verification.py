#!/usr/bin/env python3
"""Run the Cook's membrane mesh-convergence verification study."""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
import shutil
import subprocess
import sys
from pathlib import Path

MESH_FAMILIES = ("structured", "triangular")
MEMBRANE_AREA_MM2 = 1440.0


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run the Cook's membrane mesh-convergence study"
    )
    parser.add_argument(
        "--verification-dir",
        type=Path,
        required=True,
        help=argparse.SUPPRESS,
    )
    parser.add_argument(
        "--mesh",
        choices=MESH_FAMILIES,
        default="structured",
        help="mesh family to verify (default: structured)",
    )
    parser.add_argument(
        "--levels",
        help=(
            "comma-separated cell counts per side for the structured family, "
            "or spacing levels for the triangular family (default: reference "
            "file)"
        ),
    )
    parser.add_argument(
        "--quick",
        action="store_true",
        help="run only the first three levels as a smoke test",
    )
    parser.add_argument(
        "--keep-going",
        action="store_true",
        help="continue the sweep after an individual case fails",
    )
    parser.add_argument(
        "--reuse",
        action="store_true",
        help="reuse completed cases under verification/work",
    )
    return parser.parse_args()


def ignored(_directory: str, names: list[str]) -> set[str]:
    ignored_names = {
        "verification",
        "regressionTests",
        "postProcessing",
    }
    ignored_names.update(name for name in names if name.startswith("processor"))
    ignored_names.update(name for name in names if name.startswith("log."))
    return ignored_names.intersection(names)


def set_mesh_level(block_mesh_dict: Path, cells_per_side: int) -> None:
    text = block_mesh_dict.read_text()
    pattern = re.compile(
        r"(hex\s*\([^\n]+\)\s*\()\s*\d+\s+\d+\s+1\s*(\))"
    )
    replacement = rf"\g<1>{cells_per_side} {cells_per_side} 1\g<2>"
    updated, count = pattern.subn(replacement, text, count=1)
    if count != 1:
        raise RuntimeError(f"could not update cell counts in {block_mesh_dict}")
    block_mesh_dict.write_text(updated)


def set_triangular_spacing(geometry: Path, spacing: float) -> None:
    text = geometry.read_text()
    updated, count = re.subn(
        r"(?m)^dx\s*=\s*[0-9.eE+-]+\s*;",
        f"dx = {spacing:.10g};",
        text,
        count=1,
    )
    if count != 1:
        raise RuntimeError(f"could not set mesh spacing in {geometry}")
    geometry.write_text(updated)


def read_cell_count(run_dir: Path) -> int:
    owner = run_dir / "constant" / "polyMesh" / "owner"
    match = re.search(r"nCells:\s*(\d+)", owner.read_text(errors="ignore"))
    if not match:
        raise RuntimeError(f"could not read the cell count from {owner}")
    return int(match.group(1))


def run_level(
    level: int,
    mesh: str,
    reference: dict,
    case_dir: Path,
    work_dir: Path,
    reuse: bool,
) -> dict[str, float | int]:
    if mesh == "structured":
        name = f"{level} x {level} mesh"
        run_dir = work_dir / f"mesh_{level}x{level}"
        allrun_arguments = reference["allrun_arguments"]
    else:
        spacings = reference["triangular"]["spacings_m"]
        if level > len(spacings):
            raise RuntimeError(
                f"triangular level {level} is unavailable; the finest "
                f"supplied level is {len(spacings)}"
            )
        name = f"triangular level {level} mesh"
        run_dir = work_dir / f"triangular_level_{level}"
        allrun_arguments = reference["triangular"]["allrun_arguments"]
    solver_log = run_dir / "log.solids4Foam"
    completed_case = (
        reuse
        and solver_log.exists()
        and re.search(r"^End\s*$", solver_log.read_text(), re.MULTILINE)
    )
    if completed_case:
        print(f"Reusing completed {name} in {run_dir}")
    else:
        if run_dir.exists():
            shutil.rmtree(run_dir)
        shutil.copytree(case_dir, run_dir, ignore=ignored, symlinks=True)
        if mesh == "structured":
            set_mesh_level(run_dir / "system" / "blockMeshDict", level)
        else:
            set_triangular_spacing(
                run_dir / "cooksMembrane.geo", float(spacings[level - 1])
            )

        print(f"Running {name} in {run_dir}", flush=True)
        command = ["./Allrun", *allrun_arguments]
        with (run_dir / "log.Allverify").open("w") as log:
            completed = subprocess.run(
                command,
                cwd=run_dir,
                stdout=log,
                stderr=subprocess.STDOUT,
                check=False,
            )
        if completed.returncode:
            raise RuntimeError(
                f"{name} run failed; see {run_dir / 'log.Allverify'}"
            )

    candidates = sorted(
        run_dir.glob(
            "postProcessing/*/solidPointDisplacement_pointDisp.dat"
        )
    )
    if not candidates:
        raise RuntimeError(f"point displacement output not found in {run_dir}")

    rows = [
        line.split()
        for line in candidates[-1].read_text().splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    ]
    if not rows or len(rows[-1]) < 4:
        raise RuntimeError(f"invalid point displacement output in {candidates[-1]}")
    displacement_y = float(rows[-1][2])
    if not math.isfinite(displacement_y):
        raise RuntimeError(f"non-finite displacement in {candidates[-1]}")

    clock_time = math.nan
    if solver_log.exists():
        matches = re.findall(r"ClockTime\s*=\s*([0-9.eE+-]+)", solver_log.read_text())
        if matches:
            clock_time = float(matches[-1])

    cells = level * level if mesh == "structured" else read_cell_count(run_dir)
    return {
        "cells_per_side" if mesh == "structured" else "level": level,
        "cells": cells,
        "spacing_mm": math.sqrt(MEMBRANE_AREA_MM2 / cells),
        "displacement_y_m": displacement_y,
        "clock_time_s": clock_time,
    }


def read_published_series(reference: dict, reference_dir: Path) -> list[dict]:
    """Read the published mesh-convergence series listed in the reference."""
    series = []
    for entry in reference.get("published_mesh_series", []):
        path = reference_dir / entry["file"]
        lines = [
            line
            for line in path.read_text().splitlines()
            if line.strip() and not line.lstrip().startswith("#")
        ]
        rows = list(csv.DictReader(lines))
        series.append(
            {
                "label": entry["label"],
                "path": path,
                "spacing_column": list(rows[0].keys()).index("spacing_mm") + 1,
                "value_column": list(rows[0].keys()).index(entry["column"]) + 1,
                "cells_per_side": [int(row["cells_per_side"]) for row in rows],
                "displacement_m": [float(row[entry["column"]]) for row in rows],
            }
        )
    return series


def create_plot(
    csv_path: Path,
    published: list[dict],
    target: float,
    output: Path,
) -> None:
    """Plot the sweep against the published series when gnuplot is found."""
    if shutil.which("gnuplot") is None:
        print("gnuplot not found; skipping the plot (the CSV results are valid)")
        return
    with csv_path.open() as stream:
        columns = next(csv.reader(stream))
    spacing = columns.index("spacing_mm") + 1
    value = columns.index("displacement_y_m") + 1
    plots = [
        f"'{csv_path}' using {spacing}:(1000*${value}) with linespoints "
        "lw 2 pt 7 lc rgb 'red' title 'solids4foam'"
    ]
    for index, series in enumerate(published):
        plots.append(
            f"'{series['path']}' using {series['spacing_column']}:"
            f"(1000*${series['value_column']}) with linespoints "
            f"dt {2 + index % 4} pt {index + 4} lc {index + 1} "
            f"title \"{series['label']}\""
        )
    plots.append(
        f"{1000 * target} with lines lc rgb 'black' dt 3 title 'Reference'"
    )
    script = (
        "set terminal pdfcairo dashed noenhanced size 7in,3.5in\n"
        f"set output '{output}'\n"
        "set datafile separator ','\n"
        "set grid\n"
        "set logscale x\n"
        "set xlabel 'Average cell spacing (mm)'\n"
        "set ylabel 'Top-right vertical displacement (mm)'\n"
        "set key outside right top box\n"
        "plot " + ", \\\n    ".join(plots) + "\n"
    )
    completed = subprocess.run(
        ["gnuplot"],
        input=script,
        text=True,
        check=False,
    )
    if completed.returncode:
        print("WARNING: gnuplot failed; the CSV results are still valid.")
    else:
        print(f"Plot: {output}")


def write_results(
    results: list[dict[str, float | int]],
    reference: dict,
    published: list[dict],
    mesh: str,
    quick: bool,
    post_dir: Path,
) -> bool:
    post_dir.mkdir(parents=True, exist_ok=True)
    prefix = "" if mesh == "structured" else f"{mesh}_"
    level_key = "cells_per_side" if mesh == "structured" else "level"
    target = float(reference["primary_reference"]["displacement_m"])
    tolerance = float(reference["primary_reference"]["relative_tolerance"])

    for row in results:
        value = float(row["displacement_y_m"])
        row["relative_error"] = abs(value - target) / abs(target)

    csv_path = post_dir / f"{prefix}mesh_convergence.csv"
    with csv_path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=results[0].keys())
        writer.writeheader()
        writer.writerows(results)

    first_error = float(results[0]["relative_error"])
    final_error = float(results[-1]["relative_error"])
    refinement_ratio = (
        float(results[0]["spacing_mm"]) / float(results[-1]["spacing_mm"])
    )
    net_order = (
        math.log(first_error / final_error) / math.log(refinement_ratio)
        if first_error > 0 and final_error > 0 and refinement_ratio > 1
        else math.nan
    )
    increasing = all(
        float(right["displacement_y_m"]) > float(left["displacement_y_m"])
        for left, right in zip(results, results[1:])
    )
    passed = (
        all(math.isfinite(float(row["displacement_y_m"])) for row in results)
        if quick
        else final_error <= tolerance and final_error < first_error and net_order > 0
    )

    final_value = float(results[-1]["displacement_y_m"])
    published_lines = ""
    if published:
        published_lines = (
            "\n## Published mesh-convergence series\n\n"
            "These series are reported for comparison only and do not enter "
            "the pass criterion.\n\n"
            "| Series | Finest cells per side | Finest displacement (mm) "
            "| Finest solids4foam difference |\n"
            "|---|---:|---:|---:|\n"
        )
        for series in published:
            value = series["displacement_m"][-1]
            published_lines += (
                f"| {series['label']} | {series['cells_per_side'][-1]} "
                f"| {1000 * value:.4g} | {(final_value - value) / value:+.2%} |\n"
            )

    summary_path = post_dir / f"{prefix}verification_summary.md"
    summary_path.write_text(
        "# Cook's membrane verification summary\n\n"
        f"- Case: {reference['case_name']}\n"
        f"- Mesh family: {mesh}\n"
        f"- Levels: {', '.join(str(row[level_key]) for row in results)}\n"
        f"- Finest displacement: {float(results[-1]['displacement_y_m']):.8g} m\n"
        f"- Primary reference: {target:.8g} m\n"
        f"- Finest relative error: {final_error:.3%}\n"
        f"- Net convergence order: {net_order:.3f}\n"
        f"- Displacement increases under refinement: {increasing}\n"
        f"- Mode: {'quick smoke test' if quick else 'full verification'}\n"
        f"- Result: {'PASS' if passed else 'FAIL'}\n"
        + published_lines
    )
    print(summary_path.read_text())
    print(f"CSV: {csv_path}")
    create_plot(
        csv_path,
        published,
        target,
        post_dir / f"{prefix}cooksMembrane_mesh_convergence.pdf",
    )
    return passed


def main() -> int:
    args = parse_args()
    verify_dir = args.verification_dir.resolve()
    case_dir = verify_dir.parent
    work_dir = verify_dir / "work"
    post_dir = verify_dir / "postProcessing"
    reference_dir = verify_dir / "reference"
    reference_file = reference_dir / "cooksMembrane_verification_references.json"
    reference = json.loads(reference_file.read_text())
    if args.mesh != "structured" and args.mesh not in reference:
        raise SystemExit(f"the {args.mesh} mesh family is not available here")
    if args.levels:
        levels = [int(value) for value in args.levels.split(",")]
    elif args.mesh == "structured":
        levels = [int(value) for value in reference["default_levels"]]
    else:
        levels = [
            int(value) for value in reference[args.mesh]["default_levels"]
        ]
    if args.quick:
        levels = levels[:3]
    if len(levels) < 2 or any(level <= 0 for level in levels):
        raise SystemExit("at least two positive mesh levels are required")
    if levels != sorted(set(levels)):
        raise SystemExit("mesh levels must be unique and increasing")

    if args.mesh == "structured":
        required = ["blockMesh", "solids4Foam"]
    else:
        required = ["gmsh", "gmshToFoam", "changeDictionary", "solids4Foam"]
    missing = [command for command in required if shutil.which(command) is None]
    if missing:
        raise SystemExit(f"required command(s) not found: {', '.join(missing)}")

    work_dir.mkdir(parents=True, exist_ok=True)
    results: list[dict[str, float | int]] = []
    failures = 0
    for level in levels:
        try:
            results.append(
                run_level(
                    level,
                    args.mesh,
                    reference,
                    case_dir,
                    work_dir,
                    args.reuse,
                )
            )
        except RuntimeError as error:
            print(f"ERROR: {error}", file=sys.stderr)
            failures += 1
            if not args.keep_going:
                return 1
    if len(results) < 2:
        print("ERROR: fewer than two mesh levels completed", file=sys.stderr)
        return 1
    published = read_published_series(reference, reference_dir)
    passed = write_results(
        results, reference, published, args.mesh, args.quick, post_dir
    )
    return 0 if passed and failures == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
