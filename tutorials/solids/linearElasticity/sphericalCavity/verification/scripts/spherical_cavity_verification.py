#!/usr/bin/env python3
"""Run the spherical-cavity analytical mesh-convergence study."""

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


SCRIPT_DIR = Path(__file__).resolve().parent
VERIFY_DIR = SCRIPT_DIR.parent
CASE_DIR = VERIFY_DIR.parent
WORK_DIR = VERIFY_DIR / "work"
POST_DIR = VERIFY_DIR / "postProcessing"
REFERENCE_FILE = (
    VERIFY_DIR / "reference" / "sphericalCavity_verification_references.json"
)
DOMAIN_VOLUME = 0.5**3 - math.pi * 0.2**3 / 6.0


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run the sphericalCavity polyhedral mesh study"
    )
    parser.add_argument(
        "--levels",
        help="comma-separated minimum Gmsh spacings in metres",
    )
    parser.add_argument(
        "--quick",
        action="store_true",
        help="run only the first two levels as a smoke test",
    )
    parser.add_argument(
        "--reuse",
        action="store_true",
        help="reuse completed cases under verification/work",
    )
    parser.add_argument(
        "--keep-going",
        action="store_true",
        help="continue the sweep after an individual case fails",
    )
    return parser.parse_args()


def ignored(_directory: str, names: list[str]) -> set[str]:
    ignored_names = {"verification", "regressionTests", "postProcessing"}
    ignored_names.update(name for name in names if name.startswith("processor"))
    ignored_names.update(name for name in names if name.startswith("log."))
    return ignored_names.intersection(names)


def level_name(spacing: float) -> str:
    return f"mesh_{format(spacing, '.8g').replace('.', 'p')}"


def set_mesh_spacing(geometry_file: Path, minimum: float, ratio: float) -> None:
    text = geometry_file.read_text()
    substitutions = {
        "minDeltaX": minimum,
        "maxDeltaX": ratio * minimum,
    }
    for name, value in substitutions.items():
        pattern = re.compile(rf"({name}\s*=\s*)[0-9.eE+-]+(\s*;)")
        text, count = pattern.subn(rf"\g<1>{value:.10g}\g<2>", text, count=1)
        if count != 1:
            raise RuntimeError(f"could not set {name} in {geometry_file}")
    geometry_file.write_text(text)


def extract_norms(log_text: str, marker: str, component: int | None = None) -> tuple[float, float]:
    sections = log_text.split(marker)
    if len(sections) < 2:
        raise RuntimeError(f"could not find '{marker}' in solver log")
    section = sections[-1]
    if component is not None:
        match = re.search(
            rf"Component:\s*{component}\b(.*?)(?=Component:|$)",
            section,
            re.DOTALL,
        )
        if not match:
            raise RuntimeError(f"could not find component {component} after '{marker}'")
        section = match.group(1)
    match = re.search(
        r"Norms: mean L1, mean L2, LInf:\s*"
        r"([0-9.eE+-]+)\s+([0-9.eE+-]+)\s+([0-9.eE+-]+)",
        section,
    )
    if not match:
        raise RuntimeError(f"could not extract norms after '{marker}'")
    return float(match.group(2)), float(match.group(3))


def count_cells(run_dir: Path) -> int:
    completed = subprocess.run(
        ["checkMesh", "-constant"],
        cwd=run_dir,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )
    if completed.returncode:
        raise RuntimeError(f"checkMesh failed in {run_dir}")
    match = re.search(r"^\s*cells:\s*(\d+)\s*$", completed.stdout, re.MULTILINE)
    if not match:
        raise RuntimeError(f"could not extract cell count in {run_dir}")
    return int(match.group(1))


def run_level(spacing: float, reference: dict, reuse: bool) -> dict[str, float | int]:
    run_dir = WORK_DIR / level_name(spacing)
    solver_log = run_dir / "log.solids4Foam"
    completed_case = (
        reuse
        and solver_log.exists()
        and re.search(r"^End\s*$", solver_log.read_text(), re.MULTILINE)
    )
    if completed_case:
        print(f"Reusing completed spacing {spacing:g} m in {run_dir}")
    else:
        if run_dir.exists():
            shutil.rmtree(run_dir)
        shutil.copytree(CASE_DIR, run_dir, ignore=ignored, symlinks=True)
        set_mesh_spacing(
            run_dir / "sphericalCavity.geo",
            spacing,
            float(reference["max_to_min_spacing_ratio"]),
        )
        command = [
            "./Allrun",
            reference["approach"],
            reference["mesh_type"],
        ]
        print(f"Running spacing {spacing:g} m in {run_dir}", flush=True)
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
                f"spacing {spacing:g} m failed; see {run_dir / 'log.Allverify'}"
            )
        if not solver_log.exists():
            raise RuntimeError(f"solver log was not created in {run_dir}")

    log_text = solver_log.read_text()
    displacement_l2, displacement_linf = extract_norms(
        log_text, "Writing DDifference field"
    )
    stress_l2, stress_linf = extract_norms(
        log_text, "Writing cellStressDifference field", component=5
    )
    cell_count = count_cells(run_dir)
    clock_matches = re.findall(r"ClockTime\s*=\s*([0-9.eE+-]+)", log_text)

    result: dict[str, float | int] = {
        "min_spacing_m": spacing,
        "cells": cell_count,
        "effective_spacing_m": (DOMAIN_VOLUME / cell_count) ** (1.0 / 3.0),
        "displacement_l2_m": displacement_l2,
        "displacement_linf_m": displacement_linf,
        "stress_zz_l2_pa": stress_l2,
        "stress_zz_linf_pa": stress_linf,
        "clock_time_s": float(clock_matches[-1]) if clock_matches else math.nan,
    }
    if not all(math.isfinite(float(value)) for value in result.values()):
        raise RuntimeError(f"non-finite result extracted from {solver_log}")
    return result


def net_order(results: list[dict[str, float | int]], metric: str) -> float:
    coarse = float(results[0][metric])
    fine = float(results[-1][metric])
    spacing_ratio = float(results[0]["effective_spacing_m"]) / float(
        results[-1]["effective_spacing_m"]
    )
    if coarse <= 0 or fine <= 0 or spacing_ratio <= 1:
        return math.nan
    return math.log(coarse / fine) / math.log(spacing_ratio)


def write_results(results: list[dict[str, float | int]], reference: dict, quick: bool) -> bool:
    POST_DIR.mkdir(parents=True, exist_ok=True)
    metrics = reference["acceptance"]["required_metrics"]
    orders = {metric: net_order(results, metric) for metric in metrics}

    with (POST_DIR / "mesh_convergence.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=results[0].keys())
        writer.writeheader()
        writer.writerows(results)

    finite_positive = all(
        math.isfinite(float(row[metric])) and float(row[metric]) > 0
        for row in results
        for metric in metrics
    )
    if quick:
        passed = finite_positive
    else:
        minimum_order = float(reference["acceptance"]["minimum_net_order"])
        passed = finite_positive and all(
            float(results[-1][metric]) < float(results[0][metric])
            and orders[metric] > minimum_order
            for metric in metrics
        )

    level_text = ", ".join(
        f"{float(row['min_spacing_m']):g}" for row in results
    )
    lines = [
        "# Spherical-cavity verification summary",
        "",
        f"- Levels: {level_text} m",
        f"- Finest mesh: {int(results[-1]['cells'])} cells",
        f"- Mode: {'quick smoke test' if quick else 'full verification'}",
    ]
    for metric in metrics:
        lines.append(
            f"- {metric}: finest {float(results[-1][metric]):.8g}, "
            f"net order {orders[metric]:.3f}"
        )
    lines.append(f"- Result: {'PASS' if passed else 'FAIL'}")
    summary = "\n".join(lines) + "\n"
    (POST_DIR / "verification_summary.md").write_text(summary)
    print(summary)
    print(f"CSV: {POST_DIR / 'mesh_convergence.csv'}")
    return passed


def main() -> int:
    args = parse_args()
    reference = json.loads(REFERENCE_FILE.read_text())
    if args.levels:
        levels = [float(value) for value in args.levels.split(",")]
    else:
        levels = [float(value) for value in reference["default_min_spacings_m"]]
    if args.quick:
        levels = levels[:2]
    if len(levels) < 2 or any(level <= 0 for level in levels):
        raise SystemExit("at least two positive mesh spacings are required")
    if any(right >= left for left, right in zip(levels, levels[1:])):
        raise SystemExit("mesh spacings must be strictly decreasing")

    required = ["gmsh", "polyDualMesh", "checkMesh", "solids4Foam"]
    missing = [command for command in required if shutil.which(command) is None]
    if missing:
        raise SystemExit(f"required command(s) not found: {', '.join(missing)}")

    WORK_DIR.mkdir(parents=True, exist_ok=True)
    results: list[dict[str, float | int]] = []
    failures = 0
    for spacing in levels:
        try:
            results.append(run_level(spacing, reference, args.reuse))
        except RuntimeError as error:
            print(f"ERROR: {error}", file=sys.stderr)
            failures += 1
            if not args.keep_going:
                return 1
    if len(results) < 2:
        print("ERROR: fewer than two mesh levels completed", file=sys.stderr)
        return 1
    passed = write_results(results, reference, args.quick)
    return 0 if passed and failures == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
