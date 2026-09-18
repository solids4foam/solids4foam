#!/usr/bin/env python3
"""Run the curvedCantilever mesh-convergence verification study.

The tutorial already ships the Timoshenko curved-beam analytical solution as a
function object and samples both the computed and the analytical stress along
the radial line at theta = 45 degrees. This study refines the mesh, measures
the error against that analytical solution on the sampled line, and checks the
observed order of accuracy.
"""

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
REFERENCE_DIR = VERIFY_DIR / "reference"
REFERENCE_FILE = REFERENCE_DIR / "curvedCantilever_verification_references.json"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run the curvedCantilever verification variants"
    )
    parser.add_argument("--variants", help="comma-separated variant names")
    parser.add_argument(
        "--levels",
        help="comma-separated mesh levels (1 is the coarsest)",
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
    """Copy only the inputs a fresh run needs, never a previous result."""
    ignored_names = {"verification", "regressionTests", "postProcessing", "images"}
    ignored_names.update(name for name in names if name.startswith("processor"))
    ignored_names.update(name for name in names if name.startswith("log."))

    directory_path = Path(directory)
    if directory_path == CASE_DIR:
        ignored_names.update(
            name
            for name in names
            if name != "0" and re.fullmatch(r"[0-9]+(?:\.[0-9]+)?", name)
        )
        ignored_names.add("sigmaAtTheta45deg.png")
    if directory_path.name == "constant":
        ignored_names.add("polyMesh")
    if directory_path.name == "src":
        ignored_names.update({"lnInclude", "log.wmake"})
    if directory_path.name == "Make":
        ignored_names.update(
            name for name in names if name not in {"files", "options"}
        )

    return ignored_names.intersection(names)


def set_resolution(run_dir: Path, level: int, tutorial_level: int) -> tuple[int, int]:
    """Scale the single block of the tutorial mesh by a power of two.

    The case is two dimensional, so only the circumferential and radial
    divisions are scaled; the single spanwise cell is left alone.
    """
    block_mesh_dict = run_dir / "system" / "blockMeshDict"
    pattern = re.compile(
        r"(^\s*hex\s*\([^)]*\)\s*)\(\s*(\d+)\s+(\d+)\s+(\d+)\s*\)", re.MULTILINE
    )
    shift = level - tutorial_level
    numerator = 2**shift if shift > 0 else 1
    denominator = 2 ** (-shift) if shift < 0 else 1
    divisions: list[int] = []

    def scale(match: re.Match[str]) -> str:
        scaled = []
        for index in (2, 3):
            value = int(match.group(index)) * numerator
            if value % denominator:
                raise RuntimeError(
                    f"block division {match.group(index)} in {block_mesh_dict} "
                    f"is not divisible by {denominator}; level {level} is not "
                    "reachable from the supplied mesh"
                )
            scaled.append(value // denominator)
        divisions.extend(scaled)
        return f"{match.group(1)}({scaled[0]} {scaled[1]} {match.group(4)})"

    text, count = pattern.subn(scale, block_mesh_dict.read_text())
    if count != 1:
        raise RuntimeError(
            f"expected one block in {block_mesh_dict}, found {count}"
        )
    block_mesh_dict.write_text(text)
    return divisions[0], divisions[1]


def read_divisions(run_dir: Path) -> tuple[int, int]:
    """Read the block divisions already written into a case."""
    block_mesh_dict = run_dir / "system" / "blockMeshDict"
    match = re.search(
        r"^\s*hex\s*\([^)]*\)\s*\(\s*(\d+)\s+(\d+)\s+\d+\s*\)",
        block_mesh_dict.read_text(),
        re.MULTILINE,
    )
    if not match:
        raise RuntimeError(f"no block divisions found in {block_mesh_dict}")
    return int(match.group(1)), int(match.group(2))


def check_solver_log(solver_log: Path, label: str) -> None:
    """Reject a run that crashed, was killed, or never reached the end."""
    text = solver_log.read_text(errors="replace")
    if re.search(
        r"FOAM FATAL|FOAM aborting|^ERROR$|\[stack trace\]", text, re.MULTILINE
    ):
        raise RuntimeError(f"{label} did not complete; see {solver_log}")
    if not re.search(r"^End\s*$", text, re.MULTILINE):
        raise RuntimeError(
            f"{label} did not reach the end of the run; see {solver_log}"
        )


def read_sampled_stress(run_dir: Path, sample: dict) -> list[list[float]]:
    candidates = sorted(run_dir.glob(sample["file"]))
    if not candidates:
        raise RuntimeError(
            f"no sampled stress data found in {run_dir}; the tutorial Allrun "
            "samples the theta = 45 degree line through its own sample "
            "dictionary"
        )
    rows = []
    for line in candidates[-1].read_text().splitlines():
        fields = line.split()
        if len(fields) >= 15 and not line.lstrip().startswith("#"):
            rows.append([float(field) for field in fields])
    if not rows:
        raise RuntimeError(f"no numeric rows in {candidates[-1]}")
    return rows


def stress_errors(rows: list[list[float]], sample: dict) -> dict[str, float]:
    """Relative L2 and Linf errors of the sampled stress components.

    The components are normalised together by the L2 norm of the analytical
    stress over the line, so a single figure describes the whole profile
    rather than one component that happens to be small.
    """
    squared_error = 0.0
    squared_reference = 0.0
    max_error = 0.0
    max_reference = 0.0
    for row in rows:
        for component in sample["components"].values():
            analytical = row[int(component["analytical_column"])]
            computed = row[int(component["computed_column"])]
            squared_error += (computed - analytical) ** 2
            squared_reference += analytical**2
            max_error = max(max_error, abs(computed - analytical))
            max_reference = max(max_reference, abs(analytical))
    if squared_reference <= 0.0 or max_reference <= 0.0:
        raise RuntimeError("the analytical stress on the sample line is zero")
    return {
        "relative_l2_error": math.sqrt(squared_error / squared_reference),
        "relative_linf_error": max_error / max_reference,
        "samples": len(rows),
    }


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
    level: int,
    reference: dict,
    reuse: bool,
) -> dict:
    run_dir = WORK_DIR / variant_name / f"mesh{level}"
    solver_log = run_dir / "log.solids4Foam"
    completed_case = (
        reuse
        and solver_log.exists()
        and re.search(r"^End\s*$", solver_log.read_text(), re.MULTILINE)
    )

    if completed_case:
        # The scaling is relative to the shipped mesh, so a reused case is
        # read back rather than rescaled.
        print(f"Reusing {variant_name} mesh{level} in {run_dir}")
        circumferential, radial = read_divisions(run_dir)
    else:
        if run_dir.exists():
            shutil.rmtree(run_dir)
        shutil.copytree(CASE_DIR, run_dir, ignore=ignored, symlinks=True)
        circumferential, radial = set_resolution(
            run_dir, level, int(reference["tutorial_level"])
        )
        command = ["./Allrun"]
        if variant["approach"]:
            command.append(variant["approach"])
        print(
            f"Running {variant_name} mesh{level} "
            f"({circumferential} x {radial}) in {run_dir}",
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
                f"{variant_name} mesh{level} failed; "
                f"see {run_dir / 'log.Allverify'}"
            )
        if not solver_log.exists():
            raise RuntimeError(f"solver log was not created in {run_dir}")

    check_solver_log(solver_log, f"{variant_name} mesh{level}")

    rows = read_sampled_stress(run_dir, reference["sample"])
    cell_count = count_cells(run_dir)
    clock_matches = re.findall(
        r"ClockTime\s*=\s*([0-9.eE+-]+)", solver_log.read_text()
    )
    result: dict = {
        "variant": variant_name,
        "level": level,
        "circumferential_cells": circumferential,
        "radial_cells": radial,
        "cells": cell_count,
        "effective_spacing_m": cell_count ** (-1.0 / 2.0),
        "clock_time_s": float(clock_matches[-1]) if clock_matches else math.nan,
    }
    result.update(stress_errors(rows, reference["sample"]))
    if not all(
        math.isfinite(float(value))
        for key, value in result.items()
        if key != "variant"
    ):
        raise RuntimeError(f"non-finite result extracted from {run_dir}")
    return result


def net_order(results: list[dict], metric: str) -> float:
    coarse = float(results[0][metric])
    fine = float(results[-1][metric])
    spacing_ratio = float(results[0]["effective_spacing_m"]) / float(
        results[-1]["effective_spacing_m"]
    )
    if coarse <= 0 or fine <= 0 or spacing_ratio <= 1:
        return math.nan
    return math.log(coarse / fine) / math.log(spacing_ratio)


def write_results(
    results: list[dict],
    variant_names: list[str],
    reference: dict,
    quick: bool,
) -> bool:
    POST_DIR.mkdir(parents=True, exist_ok=True)
    acceptance = reference["acceptance"]
    metrics = ("relative_l2_error", "relative_linf_error")
    grouped = {
        name: [row for row in results if row["variant"] == name]
        for name in variant_names
    }

    with (POST_DIR / "mesh_convergence.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=results[0].keys())
        writer.writeheader()
        writer.writerows(results)

    (POST_DIR / "profiles").mkdir(parents=True, exist_ok=True)
    for name in variant_names:
        (POST_DIR / "profiles" / f"{name}_errors.txt").write_text(
            "# spacing_m relative_l2 relative_linf\n"
            + "".join(
                f"{float(row['effective_spacing_m']):.10g} "
                f"{float(row['relative_l2_error']):.10g} "
                f"{float(row['relative_linf_error']):.10g}\n"
                for row in grouped[name]
            )
        )

    orders = {
        name: {metric: net_order(grouped[name], metric) for metric in metrics}
        for name in variant_names
    }
    # The reference here is an exact analytical solution rather than a
    # digitised curve, so the error itself is the convergence measure.
    passed = all(
        math.isfinite(float(row[metric])) and float(row[metric]) > 0.0
        for row in results
        for metric in metrics
    )
    if not quick:
        for name in variant_names:
            errors = [float(row["relative_l2_error"]) for row in grouped[name]]
            passed = (
                passed
                and all(right < left for left, right in zip(errors, errors[1:]))
                and errors[-1] <= float(acceptance["relative_l2_tolerance"])
                and orders[name]["relative_l2_error"]
                > float(acceptance["minimum_net_order"])
            )

    lines = [
        "# curvedCantilever verification summary",
        "",
        f"- Mode: {'quick smoke test' if quick else 'full verification'}",
        f"- Variants: {', '.join(variant_names)}",
        "- Reference: Timoshenko curved-beam analytical stress on the"
        " theta = 45 degree line",
    ]
    for name in variant_names:
        variant_results = grouped[name]
        lines.extend(
            [
                "",
                f"## {name}",
                "",
                "- Levels: "
                + ", ".join(str(row["level"]) for row in variant_results),
                f"- Finest mesh: {int(variant_results[-1]['cells'])} cells "
                f"({int(variant_results[-1]['circumferential_cells'])} x "
                f"{int(variant_results[-1]['radial_cells'])})",
                "- Relative L2 error per level: "
                + ", ".join(
                    f"{float(row['relative_l2_error']):.4g}"
                    for row in variant_results
                ),
                "- Relative Linf error per level: "
                + ", ".join(
                    f"{float(row['relative_linf_error']):.4g}"
                    for row in variant_results
                ),
                f"- Net order (L2): {orders[name]['relative_l2_error']:.3f}",
                f"- Net order (Linf): {orders[name]['relative_linf_error']:.3f}",
            ]
        )
    lines.extend(["", f"- Result: {'PASS' if passed else 'FAIL'}"])
    summary = "\n".join(lines) + "\n"
    (POST_DIR / "verification_summary.md").write_text(summary)
    print(summary)
    print(f"CSV: {POST_DIR / 'mesh_convergence.csv'}")
    return passed


def create_plot(variant_names: list[str]) -> None:
    plot_script = SCRIPT_DIR / "plotErrorNorms.gnuplot"
    if shutil.which("gnuplot") is None or not plot_script.is_file():
        return
    output = POST_DIR / "curvedCantilever_errorNorms.pdf"
    completed = subprocess.run(
        [
            "gnuplot",
            "-e",
            f"profiles='{POST_DIR / 'profiles'}'; "
            f"variant='{variant_names[0]}'; "
            f"output='{output}'",
            str(plot_script),
        ],
        cwd=VERIFY_DIR,
        text=True,
        check=False,
    )
    if completed.returncode:
        print("WARNING: gnuplot failed; the CSV results are still valid.")
    else:
        print(f"Plot: {output}")


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

    missing = [
        command
        for command in ("blockMesh", "checkMesh", "postProcess", "solids4Foam")
        if shutil.which(command) is None
    ]
    if missing:
        raise SystemExit(f"required command(s) not found: {', '.join(missing)}")

    WORK_DIR.mkdir(parents=True, exist_ok=True)
    POST_DIR.mkdir(parents=True, exist_ok=True)
    results: list[dict] = []
    failures = 0
    for name in variant_names:
        for level in levels:
            try:
                results.append(
                    run_level(name, all_variants[name], level, reference, args.reuse)
                )
            except RuntimeError as error:
                print(f"ERROR: {error}", file=sys.stderr)
                failures += 1
                if not args.keep_going:
                    return 1
    if any(
        sum(row["variant"] == name for row in results) < 2 for name in variant_names
    ):
        print("ERROR: fewer than two levels completed for a variant", file=sys.stderr)
        return 1
    passed = write_results(results, variant_names, reference, args.quick)
    create_plot(variant_names)
    return 0 if passed and failures == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
