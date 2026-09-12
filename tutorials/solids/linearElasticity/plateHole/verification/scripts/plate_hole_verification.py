#!/usr/bin/env python3
"""Run the plateHole mesh-convergence verification study.

The study refines the tutorial blockMesh through the mesh family of
solid-benchmarks/linearElasticity/plateHole, reads the error norms that the
plateHoleAnalyticalSolution function object prints against the Kirsch
analytical solution, and checks the observed order of accuracy.
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
REFERENCE_FILE = REFERENCE_DIR / "plateHole_verification_references.json"

# Components of the symmetric stress tensor reported by the function object.
STRESS_XX_COMPONENT = 0


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run the plateHole verification variants"
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
            name for name in names if not name.startswith("files")
            and name != "options"
        )

    return ignored_names.intersection(names)


def set_resolution(run_dir: Path, level: int, tutorial_level: int) -> None:
    """Rescale the in-plane block divisions from the mesh the tutorial ships.

    The tutorial dictionary holds one member of the benchmark mesh family, so
    a level is reached by scaling it by a power of two. The case is two
    dimensional, so the third division, which spans the single cell between
    the empty frontAndBack patches, is left untouched, as is the grading.
    """
    block_mesh_dict = run_dir / "system" / "blockMeshDict"
    pattern = re.compile(
        r"(^\s*hex\s*\([^)]*\)\s*\S*\s*)\(\s*(\d+)\s+(\d+)\s+(\d+)\s*\)",
        re.MULTILINE,
    )
    shift = level - tutorial_level
    numerator = 2 ** shift if shift > 0 else 1
    denominator = 2 ** (-shift) if shift < 0 else 1

    def scale(match: re.Match[str]) -> str:
        divisions = []
        for index in (2, 3):
            value = int(match.group(index)) * numerator
            if value % denominator:
                raise RuntimeError(
                    f"block division {match.group(index)} in {block_mesh_dict} "
                    f"is not divisible by {denominator}; level {level} is not "
                    "reachable from the supplied mesh"
                )
            divisions.append(value // denominator)
        return (
            f"{match.group(1)}({divisions[0]} {divisions[1]} {match.group(4)})"
        )

    text, count = pattern.subn(scale, block_mesh_dict.read_text())
    if not count:
        raise RuntimeError(f"no block divisions found in {block_mesh_dict}")
    block_mesh_dict.write_text(text)


def set_convergence_tolerance(run_dir: Path, approach: str, rtol: float) -> None:
    """Tighten the outer solution tolerance of the copied case.

    The tutorial ships the solidModel default `rTol` of 1e-6, which leaves an
    iterative error of the same size as the discretisation error on the finer
    meshes and flattens the convergence curve. A verification study needs the
    iterative error to sit well below the discretisation error, so the run
    copy asks for a much tighter tolerance.
    """
    properties = run_dir / "constant" / f"solidProperties.{approach}"
    if not properties.is_file():
        raise RuntimeError(f"{properties} not found")
    text = properties.read_text()
    if re.search(r"^\s*rTol\s+", text, re.MULTILINE):
        text = re.sub(
            r"^(\s*rTol\s+)[0-9.eE+-]+;",
            rf"\g<1>{rtol:g};",
            text,
            count=1,
            flags=re.MULTILINE,
        )
    else:
        pattern = re.compile(r"(Coeffs\"?\s*\n\{\n)")
        text, count = pattern.subn(
            rf"\g<1>    // Tightened by the verification driver\n"
            rf"    rTol            {rtol:g};\n\n",
            text,
            count=1,
        )
        if not count:
            raise RuntimeError(
                f"could not find a solidModel coefficients dictionary in "
                f"{properties}"
            )
    properties.write_text(text)


def extract_norms(log_text: str, marker: str) -> tuple[float, float, float]:
    """Return the mean L1, mean L2 and LInf norms printed after a marker.

    The plateHoleAnalyticalSolution function object writes, for each
    difference field,

        Writing DDifference field
            Norms: mean L1, mean L2, LInf:
            1.23e-08 1.56e-08 6.54e-08

    and the last occurrence in the log belongs to the final time.
    """
    sections = log_text.split(marker)
    if len(sections) < 2:
        raise RuntimeError(f"could not find '{marker}' in the solver log")
    match = re.search(
        r"Norms: mean L1, mean L2, LInf:\s*\n\s*"
        r"([0-9.eE+-]+)\s+([0-9.eE+-]+)\s+([0-9.eE+-]+)",
        sections[-1],
    )
    if not match:
        raise RuntimeError(f"could not extract the norms after '{marker}'")
    return (
        float(match.group(1)),
        float(match.group(2)),
        float(match.group(3)),
    )


def extract_stress_norms(
    log_text: str, component: int
) -> tuple[float, float, float]:
    """Return the norms of one component of the cell stress difference.

    The function object reports the XX, XY and YY components in turn as
    components 0, 1 and 3 of the symmetric tensor.
    """
    sections = log_text.split("Writing cellStressDifference field")
    if len(sections) < 2:
        raise RuntimeError(
            "could not find 'Writing cellStressDifference field' in the "
            "solver log"
        )
    return extract_norms(sections[-1], f"Component: {component}")


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


def check_solver_log(solver_log: Path, label: str) -> None:
    """Reject a run that crashed, was killed, or never reached the end.

    A tutorial Allrun can return zero even when the solver aborts or is
    killed by the operating system, so the log has to be inspected directly.
    """
    text = solver_log.read_text(errors="replace")
    if re.search(r"FOAM FATAL|FOAM aborting|^ERROR$|\[stack trace\]", text,
                 re.MULTILINE):
        raise RuntimeError(
            f"{label} did not complete; see {solver_log}. A level that is "
            "too large for the available memory is the usual cause."
        )
    if not re.search(r"^End\s*$", text, re.MULTILINE):
        raise RuntimeError(
            f"{label} did not reach the end of the run; see {solver_log}"
        )


def run_level(
    variant_name: str,
    variant: dict,
    level: int,
    reference: dict,
    reuse: bool,
) -> dict[str, float | int | str]:
    run_dir = WORK_DIR / variant_name / f"mesh{level}"
    solver_log = run_dir / "log.solids4Foam"
    completed_case = (
        reuse
        and solver_log.exists()
        and re.search(r"^End\s*$", solver_log.read_text(), re.MULTILINE)
    )

    if completed_case:
        print(f"Reusing {variant_name} mesh{level} in {run_dir}")
    else:
        if run_dir.exists():
            shutil.rmtree(run_dir)
        shutil.copytree(CASE_DIR, run_dir, ignore=ignored, symlinks=True)
        set_resolution(run_dir, level, int(reference["tutorial_level"]))
        if "rtol" in variant:
            set_convergence_tolerance(
                run_dir, variant["approach"], float(variant["rtol"])
            )
        print(f"Running {variant_name} mesh{level} in {run_dir}", flush=True)
        with (run_dir / "log.Allverify").open("w") as log:
            completed = subprocess.run(
                ["./Allrun", variant["approach"]],
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
            raise RuntimeError(
                f"solver log was not created in {run_dir}; the "
                f"{variant['approach']} approach may be unavailable in this "
                "installation"
            )

    check_solver_log(solver_log, f"{variant_name} mesh{level}")

    log_text = solver_log.read_text()
    displacement = extract_norms(log_text, "Writing DDifference field")
    point_displacement = extract_norms(
        log_text, "Writing pointDDifference field"
    )
    stress_xx = extract_stress_norms(log_text, STRESS_XX_COMPONENT)
    cell_count = count_cells(run_dir)
    clock_matches = re.findall(r"ClockTime\s*=\s*([0-9.eE+-]+)", log_text)
    result: dict[str, float | int | str] = {
        "variant": variant_name,
        "approach": variant["approach"],
        "level": level,
        "cells": cell_count,
        "effective_spacing_m": math.sqrt(
            float(reference["domain_area_m2"]) / cell_count
        ),
        "displacement_l1_m": displacement[0],
        "displacement_l2_m": displacement[1],
        "displacement_linf_m": displacement[2],
        "point_displacement_l2_m": point_displacement[1],
        "stress_xx_l1_pa": stress_xx[0],
        "stress_xx_l2_pa": stress_xx[1],
        "stress_xx_linf_pa": stress_xx[2],
        "clock_time_s": float(clock_matches[-1]) if clock_matches else math.nan,
    }
    numeric_values = (
        value
        for key, value in result.items()
        if key not in {"variant", "approach"}
    )
    if not all(math.isfinite(float(value)) for value in numeric_values):
        raise RuntimeError(f"non-finite result extracted from {solver_log}")
    return result


def net_order(results: list[dict[str, float | int | str]], metric: str) -> float:
    coarse = float(results[0][metric])
    fine = float(results[-1][metric])
    spacing_ratio = float(results[0]["effective_spacing_m"]) / float(
        results[-1]["effective_spacing_m"]
    )
    if coarse <= 0 or fine <= 0 or spacing_ratio <= 1:
        return math.nan
    return math.log(coarse / fine) / math.log(spacing_ratio)


def write_profiles(
    results: list[dict[str, float | int | str]], variant_name: str
) -> None:
    """Write the data the convergence plot needs, free of column indices."""
    profiles = POST_DIR / "profiles"
    profiles.mkdir(parents=True, exist_ok=True)
    lines = [
        "# spacing_m displacement_l2_m displacement_linf_m "
        "stress_xx_l2_pa stress_xx_linf_pa"
    ]
    for row in results:
        lines.append(
            " ".join(
                f"{float(row[key]):.10g}"
                for key in (
                    "effective_spacing_m",
                    "displacement_l2_m",
                    "displacement_linf_m",
                    "stress_xx_l2_pa",
                    "stress_xx_linf_pa",
                )
            )
        )
    (profiles / f"{variant_name}_convergence.txt").write_text(
        "\n".join(lines) + "\n"
    )


def write_results(
    results: list[dict[str, float | int | str]],
    variant_names: list[str],
    reference: dict,
    quick: bool,
) -> bool:
    POST_DIR.mkdir(parents=True, exist_ok=True)
    acceptance = reference["acceptance"]
    metrics = acceptance["required_metrics"]
    minimum_order = acceptance["minimum_net_order"]
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

    # The reference is the exact analytical solution rather than a digitised
    # curve, so the error itself is meaningful: it must fall monotonically
    # and at the expected order.
    passed = all(
        math.isfinite(float(row[metric])) and float(row[metric]) > 0.0
        for row in results
        for metric in metrics
    )
    if not quick:
        for name in variant_names:
            for metric in metrics:
                errors = [float(row[metric]) for row in grouped[name]]
                passed = (
                    passed
                    and all(
                        right < left
                        for left, right in zip(errors, errors[1:])
                    )
                    and orders[name][metric] > float(minimum_order[metric])
                )

    lines = [
        "# plateHole verification summary",
        "",
        f"- Mode: {'quick smoke test' if quick else 'full verification'}",
        f"- Variants: {', '.join(variant_names)}",
        "- Reference: the analytical plate-with-hole solution evaluated by the"
        " plateHoleAnalyticalSolution function object",
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
                "- Cells: "
                + ", ".join(str(int(row["cells"])) for row in variant_results),
                f"- Finest mesh: {int(variant_results[-1]['cells'])} cells",
            ]
        )
        for metric in metrics:
            lines.append(
                f"- {metric}: coarsest {float(variant_results[0][metric]):.6g},"
                f" finest {float(variant_results[-1][metric]):.6g},"
                f" net order {orders[name][metric]:.3f}"
                f" (minimum {float(minimum_order[metric]):.2f})"
            )
    lines.extend(["", f"- Result: {'PASS' if passed else 'FAIL'}"])
    summary = "\n".join(lines) + "\n"
    (POST_DIR / "verification_summary.md").write_text(summary)
    print(summary)
    print(f"CSV: {POST_DIR / 'mesh_convergence.csv'}")
    return passed


def create_plot(variant_names: list[str]) -> None:
    plot_script = SCRIPT_DIR / "plotConvergence.gnuplot"
    if shutil.which("gnuplot") is None or not plot_script.is_file():
        return
    output = POST_DIR / "plateHole_convergence.pdf"
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
    results: list[dict[str, float | int | str]] = []
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
        variant_results = [row for row in results if row["variant"] == name]
        if variant_results:
            write_profiles(variant_results, name)
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
