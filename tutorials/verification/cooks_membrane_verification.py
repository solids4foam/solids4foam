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
        "--levels",
        help="comma-separated cell counts per side (default: reference file)",
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


def run_level(
    level: int,
    allrun_arguments: list[str],
    case_dir: Path,
    work_dir: Path,
    reuse: bool,
) -> dict[str, float | int]:
    run_dir = work_dir / f"mesh_{level}x{level}"
    solver_log = run_dir / "log.solids4Foam"
    completed_case = (
        reuse
        and solver_log.exists()
        and re.search(r"^End\s*$", solver_log.read_text(), re.MULTILINE)
    )
    if completed_case:
        print(f"Reusing completed {level} x {level} mesh in {run_dir}")
    else:
        if run_dir.exists():
            shutil.rmtree(run_dir)
        shutil.copytree(case_dir, run_dir, ignore=ignored, symlinks=True)
        set_mesh_level(run_dir / "system" / "blockMeshDict", level)

        print(f"Running {level} x {level} mesh in {run_dir}", flush=True)
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
                f"{level} x {level} run failed; see {run_dir / 'log.Allverify'}"
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

    return {
        "cells_per_side": level,
        "cells": level * level,
        "spacing_mm": math.sqrt(1440.0 / (level * level)),
        "displacement_y_m": displacement_y,
        "clock_time_s": clock_time,
    }


def write_results(
    results: list[dict[str, float | int]],
    reference: dict,
    quick: bool,
    post_dir: Path,
) -> bool:
    post_dir.mkdir(parents=True, exist_ok=True)
    target = float(reference["primary_reference"]["displacement_m"])
    tolerance = float(reference["primary_reference"]["relative_tolerance"])

    for row in results:
        value = float(row["displacement_y_m"])
        row["relative_error"] = abs(value - target) / abs(target)

    csv_path = post_dir / "mesh_convergence.csv"
    with csv_path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=results[0].keys())
        writer.writeheader()
        writer.writerows(results)

    first_error = float(results[0]["relative_error"])
    final_error = float(results[-1]["relative_error"])
    refinement_ratio = (
        float(results[-1]["cells_per_side"])
        / float(results[0]["cells_per_side"])
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

    summary_path = post_dir / "verification_summary.md"
    summary_path.write_text(
        "# Cook's membrane verification summary\n\n"
        f"- Case: {reference['case_name']}\n"
        f"- Levels: {', '.join(str(row['cells_per_side']) for row in results)}\n"
        f"- Finest displacement: {float(results[-1]['displacement_y_m']):.8g} m\n"
        f"- Primary reference: {target:.8g} m\n"
        f"- Finest relative error: {final_error:.3%}\n"
        f"- Net convergence order: {net_order:.3f}\n"
        f"- Displacement increases under refinement: {increasing}\n"
        f"- Mode: {'quick smoke test' if quick else 'full verification'}\n"
        f"- Result: {'PASS' if passed else 'FAIL'}\n"
    )
    print(summary_path.read_text())
    print(f"CSV: {csv_path}")
    return passed


def main() -> int:
    args = parse_args()
    verify_dir = args.verification_dir.resolve()
    case_dir = verify_dir.parent
    work_dir = verify_dir / "work"
    post_dir = verify_dir / "postProcessing"
    reference_file = (
        verify_dir / "reference" / "cooksMembrane_verification_references.json"
    )
    reference = json.loads(reference_file.read_text())
    if args.levels:
        levels = [int(value) for value in args.levels.split(",")]
    else:
        levels = [int(value) for value in reference["default_levels"]]
    if args.quick:
        levels = levels[:3]
    if len(levels) < 2 or any(level <= 0 for level in levels):
        raise SystemExit("at least two positive mesh levels are required")
    if levels != sorted(set(levels)):
        raise SystemExit("mesh levels must be unique and increasing")

    required = ["blockMesh", "solids4Foam"]
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
                    reference["allrun_arguments"],
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
    passed = write_results(results, reference, args.quick, post_dir)
    return 0 if passed and failures == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
