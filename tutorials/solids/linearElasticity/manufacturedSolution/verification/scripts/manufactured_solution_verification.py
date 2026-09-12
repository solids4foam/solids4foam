#!/usr/bin/env python3
"""Run the manufactured-solution mesh-convergence variants."""

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
    VERIFY_DIR
    / "reference"
    / "manufacturedSolution_verification_references.json"
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run the manufacturedSolution verification variants"
    )
    parser.add_argument(
        "--variants",
        help="comma-separated variant names",
    )
    parser.add_argument(
        "--levels",
        help="comma-separated cells per coordinate direction",
    )
    parser.add_argument(
        "--quick",
        action="store_true",
        help="run the first two levels as a smoke test",
    )
    parser.add_argument(
        "--reuse",
        action="store_true",
        help="reuse completed cases under verification/work",
    )
    parser.add_argument(
        "--keep-going",
        action="store_true",
        help="continue after an individual case fails",
    )
    return parser.parse_args()


def ignored(directory: str, names: list[str]) -> set[str]:
    ignored_names = {"verification", "regressionTests", "postProcessing"}
    ignored_names.update(name for name in names if name.startswith("processor"))
    ignored_names.update(name for name in names if name.startswith("log."))

    directory_path = Path(directory)
    if directory_path == CASE_DIR:
        ignored_names.update(
            name
            for name in names
            if name != "0" and re.fullmatch(r"[0-9]+(?:\.[0-9]+)?", name)
        )
    if directory_path.name == "constant":
        ignored_names.add("polyMesh")
    if directory_path.name == "src":
        ignored_names.update({"lnInclude", "log.wmake"})
    if directory_path.name == "Make":
        ignored_names.update(name for name in names if name not in {"files", "options"})

    return ignored_names.intersection(names)


def set_resolution(run_dir: Path, divisions: int, domain_length: float) -> None:
    block_mesh_dict = run_dir / "system" / "blockMeshDict"
    text = block_mesh_dict.read_text()
    pattern = re.compile(
        r"(hex\s*\([^)]*\)\s*)\(\s*\d+\s+\d+\s+\d+\s*\)"
    )
    text, count = pattern.subn(
        rf"\g<1>({divisions} {divisions} {divisions})",
        text,
        count=1,
    )
    if count != 1:
        raise RuntimeError(f"could not set divisions in {block_mesh_dict}")
    block_mesh_dict.write_text(text)

    spacing = domain_length / divisions
    (run_dir / "gmsh" / "meshSpacing.geo").write_text(
        "// Mesh spacing for unstructured mesh\n" f"dx = {spacing:.10g};\n"
    )


def extract_norms(log_text: str, marker: str) -> tuple[float, float]:
    sections = log_text.split(marker)
    if len(sections) < 2:
        raise RuntimeError(f"could not find '{marker}' in solver log")
    match = re.search(
        r"Magnitude:\s*([0-9.eE+-]+)\s+([0-9.eE+-]+)\s+([0-9.eE+-]+)",
        sections[-1],
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


def run_level(
    variant_name: str,
    variant: dict,
    divisions: int,
    domain_length: float,
    reuse: bool,
) -> dict[str, float | int | str]:
    run_dir = WORK_DIR / variant_name / f"n{divisions}"
    solver_log = run_dir / "log.solids4Foam"
    completed_case = (
        reuse
        and solver_log.exists()
        and re.search(r"^End\s*$", solver_log.read_text(), re.MULTILINE)
    )

    if completed_case:
        print(f"Reusing {variant_name} n={divisions} in {run_dir}")
    else:
        if run_dir.exists():
            shutil.rmtree(run_dir)
        shutil.copytree(CASE_DIR, run_dir, ignore=ignored, symlinks=True)
        set_resolution(run_dir, divisions, domain_length)
        command = ["./Allrun", variant["approach"], variant["mesh"]]
        print(f"Running {variant_name} n={divisions} in {run_dir}", flush=True)
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
                f"{variant_name} n={divisions} failed; "
                f"see {run_dir / 'log.Allverify'}"
            )
        if not solver_log.exists():
            raise RuntimeError(f"solver log was not created in {run_dir}")

    log_text = solver_log.read_text()
    displacement_l2, displacement_linf = extract_norms(
        log_text, "Writing DDifference field"
    )
    stress_l2, stress_linf = extract_norms(
        log_text, "Writing sigmaDifference field"
    )
    cell_count = count_cells(run_dir)
    clock_matches = re.findall(r"ClockTime\s*=\s*([0-9.eE+-]+)", log_text)
    result: dict[str, float | int | str] = {
        "variant": variant_name,
        "approach": variant["approach"],
        "mesh": variant["mesh"],
        "divisions": divisions,
        "cells": cell_count,
        "effective_spacing_m": (domain_length**3 / cell_count) ** (1.0 / 3.0),
        "displacement_l2_m": displacement_l2,
        "displacement_linf_m": displacement_linf,
        "stress_l2_pa": stress_l2,
        "stress_linf_pa": stress_linf,
        "clock_time_s": float(clock_matches[-1]) if clock_matches else math.nan,
    }
    numeric_values = (
        value
        for key, value in result.items()
        if key not in {"variant", "approach", "mesh"}
    )
    if not all(math.isfinite(float(value)) for value in numeric_values):
        raise RuntimeError(f"non-finite result extracted from {solver_log}")
    return result


def net_order(
    results: list[dict[str, float | int | str]], metric: str
) -> float:
    coarse = float(results[0][metric])
    fine = float(results[-1][metric])
    spacing_ratio = float(results[0]["effective_spacing_m"]) / float(
        results[-1]["effective_spacing_m"]
    )
    if coarse <= 0 or fine <= 0 or spacing_ratio <= 1:
        return math.nan
    return math.log(coarse / fine) / math.log(spacing_ratio)


def write_results(
    results: list[dict[str, float | int | str]],
    variant_names: list[str],
    reference: dict,
    quick: bool,
) -> bool:
    POST_DIR.mkdir(parents=True, exist_ok=True)
    metrics = reference["acceptance"]["required_metrics"]
    grouped = {
        name: [row for row in results if row["variant"] == name]
        for name in variant_names
    }
    orders = {
        name: {metric: net_order(grouped[name], metric) for metric in metrics}
        for name in variant_names
    }

    with (POST_DIR / "mesh_convergence.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=results[0].keys())
        writer.writeheader()
        writer.writerows(results)

    finite_positive = all(
        math.isfinite(float(row[metric])) and float(row[metric]) > 0
        for row in results
        for metric in metrics
    )
    passed = finite_positive
    if not quick:
        minimum_order = float(reference["acceptance"]["minimum_net_order"])
        passed = passed and all(
            float(grouped[name][-1][metric]) < float(grouped[name][0][metric])
            and orders[name][metric] > minimum_order
            for name in variant_names
            for metric in metrics
        )

    lines = [
        "# Manufactured-solution verification summary",
        "",
        f"- Mode: {'quick smoke test' if quick else 'full verification'}",
        f"- Variants: {', '.join(variant_names)}",
    ]
    for name in variant_names:
        variant_results = grouped[name]
        levels = ", ".join(str(row["divisions"]) for row in variant_results)
        lines.extend(
            [
                "",
                f"## {name}",
                "",
                f"- Divisions: {levels}",
                f"- Finest mesh: {int(variant_results[-1]['cells'])} cells",
            ]
        )
        for metric in metrics:
            lines.append(
                f"- {metric}: finest {float(variant_results[-1][metric]):.8g}, "
                f"net order {orders[name][metric]:.3f}"
            )
    lines.extend(["", f"- Result: {'PASS' if passed else 'FAIL'}"])
    summary = "\n".join(lines) + "\n"
    (POST_DIR / "verification_summary.md").write_text(summary)
    print(summary)
    print(f"CSV: {POST_DIR / 'mesh_convergence.csv'}")
    return passed


def main() -> int:
    args = parse_args()
    reference = json.loads(REFERENCE_FILE.read_text())
    all_variants = reference["variants"]
    if args.variants:
        variant_names = [value.strip() for value in args.variants.split(",")]
    else:
        variant_names = list(reference["default_variants"])
    invalid = [name for name in variant_names if name not in all_variants]
    if not variant_names or invalid or len(set(variant_names)) != len(variant_names):
        available = ", ".join(all_variants)
        raise SystemExit(f"variants must be a unique selection from: {available}")

    if args.levels:
        levels = [int(value) for value in args.levels.split(",")]
    else:
        levels = [int(value) for value in reference["default_levels"]]
    if args.quick:
        levels = levels[:2]
    if len(levels) < 2 or any(level <= 0 for level in levels):
        raise SystemExit("at least two positive levels are required")
    if any(right <= left for left, right in zip(levels, levels[1:])):
        raise SystemExit("levels must be strictly increasing")

    selected = [all_variants[name] for name in variant_names]
    required = ["checkMesh", "solids4Foam"]
    if any(item["mesh"] in {"tet", "poly"} for item in selected):
        required.extend(["gmsh", "gmshToFoam", "createPatch"])
    if any(item["mesh"] == "poly" for item in selected):
        required.append("polyDualMesh")
    if any(item["mesh"] == "distHex" for item in selected):
        required.append("perturbMeshPoints")
    missing = [command for command in required if shutil.which(command) is None]
    if missing:
        raise SystemExit(f"required command(s) not found: {', '.join(missing)}")

    WORK_DIR.mkdir(parents=True, exist_ok=True)
    results: list[dict[str, float | int | str]] = []
    failures = 0
    for name in variant_names:
        for divisions in levels:
            try:
                results.append(
                    run_level(
                        name,
                        all_variants[name],
                        divisions,
                        float(reference["domain_length_m"]),
                        args.reuse,
                    )
                )
            except RuntimeError as error:
                print(f"ERROR: {error}", file=sys.stderr)
                failures += 1
                if not args.keep_going:
                    return 1
    if any(
        sum(row["variant"] == name for row in results) < 2
        for name in variant_names
    ):
        print("ERROR: fewer than two levels completed for a variant", file=sys.stderr)
        return 1
    passed = write_results(results, variant_names, reference, args.quick)
    return 0 if passed and failures == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
