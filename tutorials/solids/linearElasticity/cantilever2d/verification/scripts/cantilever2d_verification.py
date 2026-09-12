#!/usr/bin/env python3
"""Run the cantilever2d mesh-convergence verification study.

The tutorial already ships the Timoshenko slender-cantilever analytical
solution as the `cantileverAnalyticalSolution` function object, which writes
the L1, L2 and LInf norms of the displacement and stress error to the solver
log at the end of the run. This study refines the mesh, reads those norms back
and checks the observed order of accuracy.

It ports the order-of-accuracy sweep from the solid-benchmarks repository
(linearElasticity/cantilever) onto the block-structured tutorial mesh.
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
REFERENCE_FILE = REFERENCE_DIR / "cantilever2d_verification_references.json"

# The stress components the function object reports: xx, xy and yy.
STRESS_COMPONENTS = ("xx", "xy", "yy")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run the cantilever2d verification variants"
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

    The case is two dimensional, so only the axial and through-thickness
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


def read_error_norms(solver_log: Path) -> dict[str, float]:
    """Read the error norms written by the cantileverAnalyticalSolution object.

    The function object prints, at the end of the run,

        Writing cellStressDifference field
            Component: 0
            Norms: mean L1, mean L2, LInf:
            <L1> <L2> <LInf>
        ...
        Writing DDifference field
            Norms: mean L1, mean L2, LInf:
            <L1> <L2> <LInf>

    where the stress components reported are xx (0), xy (1) and yy (3). The
    vertex-centred solid model reports the displacement error on the points
    rather than on the cells, so the point block is accepted as a fallback.
    """
    text = solver_log.read_text(errors="replace")
    norms = re.compile(
        r"Norms: mean L1, mean L2, LInf:\s*\n\s*"
        r"([-+0-9.eE]+)\s+([-+0-9.eE]+)\s+([-+0-9.eE]+)"
    )

    def block(header: str) -> str | None:
        starts = [match.end() for match in re.finditer(re.escape(header), text)]
        return text[starts[-1]:] if starts else None

    values: dict[str, float] = {}

    stress_block = block("Writing cellStressDifference field")
    if stress_block is None:
        raise RuntimeError(
            f"no cell stress error norms in {solver_log}; the function object "
            "must be run with cellStress yes"
        )
    stress_matches = norms.findall(stress_block[: stress_block.find("Writing D")])
    if len(stress_matches) < len(STRESS_COMPONENTS):
        raise RuntimeError(
            f"expected {len(STRESS_COMPONENTS)} stress components in "
            f"{solver_log}, found {len(stress_matches)}"
        )
    for component, match in zip(STRESS_COMPONENTS, stress_matches):
        for name, field in zip(("l1", "l2", "linf"), match):
            values[f"sigma_{component}_{name}_Pa"] = float(field)

    for header in ("Writing DDifference field", "Writing pointDDifference field"):
        displacement_block = block(header)
        if displacement_block is None:
            continue
        match = norms.search(displacement_block)
        if match:
            for name, field in zip(("l1", "l2", "linf"), match.groups()):
                values[f"D_{name}_m"] = float(field)
            break
    else:
        raise RuntimeError(f"no displacement error norms in {solver_log}")

    return values


def relative_errors(norms: dict[str, float], scales: dict) -> dict[str, float]:
    """Normalise the absolute norms by the peak analytical response.

    The displacement is normalised by the analytical tip deflection and the
    stress components are combined, then normalised by the peak analytical
    bending stress, so that a single figure describes the whole stress field
    rather than the one component that happens to be small.
    """
    displacement_scale = float(scales["displacement_m"])
    stress_scale = float(scales["stress_Pa"])
    combined_l2 = math.sqrt(
        sum(norms[f"sigma_{c}_l2_Pa"] ** 2 for c in STRESS_COMPONENTS)
    )
    combined_linf = max(norms[f"sigma_{c}_linf_Pa"] for c in STRESS_COMPONENTS)
    return {
        "relative_displacement_l2_error": norms["D_l2_m"] / displacement_scale,
        "relative_displacement_linf_error": norms["D_linf_m"] / displacement_scale,
        "relative_stress_l2_error": combined_l2 / stress_scale,
        "relative_stress_linf_error": combined_linf / stress_scale,
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
        axial, thickness = read_divisions(run_dir)
    else:
        if run_dir.exists():
            shutil.rmtree(run_dir)
        shutil.copytree(CASE_DIR, run_dir, ignore=ignored, symlinks=True)
        axial, thickness = set_resolution(
            run_dir, level, int(reference["tutorial_level"])
        )
        command = ["./Allrun"]
        if variant["approach"]:
            command.append(variant["approach"])
        print(
            f"Running {variant_name} mesh{level} "
            f"({axial} x {thickness}) in {run_dir}",
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

    norms = read_error_norms(solver_log)
    cell_count = count_cells(run_dir)
    clock_matches = re.findall(
        r"ClockTime\s*=\s*([0-9.eE+-]+)", solver_log.read_text()
    )
    result: dict = {
        "variant": variant_name,
        "level": level,
        "axial_cells": axial,
        "thickness_cells": thickness,
        "cells": cell_count,
        # The cells are square, so the axial division sets the spacing.
        "effective_spacing_m": 2.0 / axial,
        "clock_time_s": float(clock_matches[-1]) if clock_matches else math.nan,
    }
    result.update(relative_errors(norms, reference["scales"]))
    result.update(norms)
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
    metrics = (
        "relative_displacement_l2_error",
        "relative_displacement_linf_error",
        "relative_stress_l2_error",
        "relative_stress_linf_error",
    )
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
            "# spacing_m displacement_l2 displacement_linf stress_l2 stress_linf\n"
            + "".join(
                f"{float(row['effective_spacing_m']):.10g} "
                + " ".join(f"{float(row[metric]):.10g}" for metric in metrics)
                + "\n"
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
            for metric, tolerance in (
                (
                    "relative_displacement_l2_error",
                    acceptance["relative_displacement_l2_tolerance"],
                ),
                (
                    "relative_stress_l2_error",
                    acceptance["relative_stress_l2_tolerance"],
                ),
            ):
                errors = [float(row[metric]) for row in grouped[name]]
                if min(errors) <= float(acceptance["roundoff_floor"]):
                    # The highOrder variant reproduces the cubic Timoshenko
                    # solution to machine precision, so from that level on the
                    # error is roundoff and neither monotone nor ordered.
                    passed = passed and errors[-1] <= float(tolerance)
                    continue
                passed = (
                    passed
                    and all(right < left for left, right in zip(errors, errors[1:]))
                    and errors[-1] <= float(tolerance)
                    and orders[name][metric] > float(acceptance["minimum_net_order"])
                )

    lines = [
        "# cantilever2d verification summary",
        "",
        f"- Mode: {'quick smoke test' if quick else 'full verification'}",
        f"- Variants: {', '.join(variant_names)}",
        "- Reference: Timoshenko slender-cantilever analytical displacement"
        " and stress over the whole domain",
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
                f"({int(variant_results[-1]['axial_cells'])} x "
                f"{int(variant_results[-1]['thickness_cells'])})",
            ]
        )
        for metric, label in (
            ("relative_displacement_l2_error", "displacement L2"),
            ("relative_displacement_linf_error", "displacement Linf"),
            ("relative_stress_l2_error", "stress L2"),
            ("relative_stress_linf_error", "stress Linf"),
        ):
            lines.append(
                f"- Relative {label} error per level: "
                + ", ".join(
                    f"{float(row[metric]):.4g}" for row in variant_results
                )
                + f" (net order {orders[name][metric]:.3f})"
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
    output = POST_DIR / "cantilever2d_errorNorms.pdf"
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
        for command in ("blockMesh", "checkMesh", "solids4Foam")
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
