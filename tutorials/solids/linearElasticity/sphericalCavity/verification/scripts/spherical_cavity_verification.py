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
        description="Run the sphericalCavity mesh-convergence study"
    )
    parser.add_argument(
        "--meshes",
        help="comma-separated mesh types (tet and/or poly)",
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


def level_name(mesh_type: str, spacing: float) -> str:
    prefix = "mesh" if mesh_type == "poly" else f"{mesh_type}_mesh"
    return f"{prefix}_{format(spacing, '.8g').replace('.', 'p')}"


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


def run_level(
    mesh_type: str, spacing: float, reference: dict, reuse: bool
) -> dict[str, float | int | str]:
    run_dir = WORK_DIR / level_name(mesh_type, spacing)
    solver_log = run_dir / "log.solids4Foam"
    completed_case = (
        reuse
        and solver_log.exists()
        and re.search(r"^End\s*$", solver_log.read_text(), re.MULTILINE)
    )
    if completed_case:
        print(
            f"Reusing completed {mesh_type} spacing {spacing:g} m in {run_dir}"
        )
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
            mesh_type,
        ]
        print(
            f"Running {mesh_type} spacing {spacing:g} m in {run_dir}",
            flush=True,
        )
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
                f"{mesh_type} spacing {spacing:g} m failed; "
                f"see {run_dir / 'log.Allverify'}"
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

    result: dict[str, float | int | str] = {
        "mesh_type": mesh_type,
        "min_spacing_m": spacing,
        "cells": cell_count,
        "effective_spacing_m": (DOMAIN_VOLUME / cell_count) ** (1.0 / 3.0),
        "displacement_l2_m": displacement_l2,
        "displacement_linf_m": displacement_linf,
        "stress_zz_l2_pa": stress_l2,
        "stress_zz_linf_pa": stress_linf,
        "clock_time_s": float(clock_matches[-1]) if clock_matches else math.nan,
    }
    numeric_values = (
        value for key, value in result.items() if key != "mesh_type"
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
    mesh_types: list[str],
    reference: dict,
    quick: bool,
) -> bool:
    POST_DIR.mkdir(parents=True, exist_ok=True)
    metrics = reference["acceptance"]["required_metrics"]
    grouped = {
        mesh_type: [row for row in results if row["mesh_type"] == mesh_type]
        for mesh_type in mesh_types
    }
    orders = {
        mesh_type: {
            metric: net_order(grouped[mesh_type], metric) for metric in metrics
        }
        for mesh_type in mesh_types
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
            float(grouped[mesh_type][-1][metric])
            < float(grouped[mesh_type][0][metric])
            and orders[mesh_type][metric] > minimum_order
            for mesh_type in mesh_types
            for metric in metrics
        )

    lines = [
        "# Spherical-cavity verification summary",
        "",
        f"- Mesh types: {', '.join(mesh_types)}",
        f"- Mode: {'quick smoke test' if quick else 'full verification'}",
    ]
    for mesh_type in mesh_types:
        mesh_results = grouped[mesh_type]
        level_text = ", ".join(
            f"{float(row['min_spacing_m']):g}" for row in mesh_results
        )
        lines.extend(
            [
                "",
                f"## {mesh_type} mesh",
                "",
                f"- Levels: {level_text} m",
                f"- Finest mesh: {int(mesh_results[-1]['cells'])} cells",
            ]
        )
        for metric in metrics:
            lines.append(
                f"- {metric}: finest {float(mesh_results[-1][metric]):.8g}, "
                f"net order {orders[mesh_type][metric]:.3f}"
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
    if args.meshes:
        mesh_types = [value.strip() for value in args.meshes.split(",")]
    else:
        mesh_types = list(reference["mesh_types"])
    invalid_meshes = [
        value for value in mesh_types if value not in {"tet", "poly"}
    ]
    if not mesh_types or invalid_meshes or len(set(mesh_types)) != len(mesh_types):
        raise SystemExit("mesh types must be a unique selection from: tet, poly")

    if args.levels:
        selected_levels = [float(value) for value in args.levels.split(",")]
        levels_by_mesh = {
            mesh_type: list(selected_levels) for mesh_type in mesh_types
        }
    else:
        levels_by_mesh = {
            mesh_type: [
                float(value)
                for value in reference["default_min_spacings_m_by_mesh"][mesh_type]
            ]
            for mesh_type in mesh_types
        }
    if args.quick:
        levels_by_mesh = {
            mesh_type: levels[:2]
            for mesh_type, levels in levels_by_mesh.items()
        }
    for levels in levels_by_mesh.values():
        if len(levels) < 2 or any(level <= 0 for level in levels):
            raise SystemExit("at least two positive mesh spacings are required")
        if any(right >= left for left, right in zip(levels, levels[1:])):
            raise SystemExit("mesh spacings must be strictly decreasing")

    required = ["gmsh", "checkMesh", "solids4Foam"]
    if "poly" in mesh_types:
        required.append("polyDualMesh")
    missing = [command for command in required if shutil.which(command) is None]
    if missing:
        raise SystemExit(f"required command(s) not found: {', '.join(missing)}")

    WORK_DIR.mkdir(parents=True, exist_ok=True)
    results: list[dict[str, float | int | str]] = []
    failures = 0
    for mesh_type in mesh_types:
        for spacing in levels_by_mesh[mesh_type]:
            try:
                results.append(
                    run_level(mesh_type, spacing, reference, args.reuse)
                )
            except RuntimeError as error:
                print(f"ERROR: {error}", file=sys.stderr)
                failures += 1
                if not args.keep_going:
                    return 1
    if any(
        sum(row["mesh_type"] == mesh_type for row in results) < 2
        for mesh_type in mesh_types
    ):
        print(
            "ERROR: fewer than two levels completed for a mesh type",
            file=sys.stderr,
        )
        return 1
    passed = write_results(results, mesh_types, reference, args.quick)
    return 0 if passed and failures == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
