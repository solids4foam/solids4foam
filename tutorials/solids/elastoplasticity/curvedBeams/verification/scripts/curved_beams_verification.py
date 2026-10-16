#!/usr/bin/env python3
"""Run the curvedBeams mesh-convergence verification study.

The study refines the tutorial blockMesh of the two curved beams, extracts the
total reaction force history on the fixed patch of the lower beam, and checks
convergence of that history towards the published solution of Neto et al.
(2016) for a given Coulomb friction coefficient.
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
REFERENCE_FILE = REFERENCE_DIR / "curvedBeams_verification_references.json"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run the curvedBeams verification variants"
    )
    parser.add_argument(
        "--variants",
        help="comma-separated variant names",
    )
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

    return ignored_names.intersection(names)


def divisions_for_level(reference: dict, level: int) -> tuple[int, int, int]:
    """Beam divisions of a mesh level, taken from the tabulated mesh family.

    The circumferential count doubles between levels, so the table is a
    factor-of-two family; it is tabulated rather than computed because the
    tutorial mesh has an odd radial count that cannot simply be halved.
    """
    divisions = reference["level_divisions"].get(str(level))
    if divisions is None:
        available = ", ".join(sorted(reference["level_divisions"]))
        raise SystemExit(f"level {level} is not defined; available levels: {available}")
    return tuple(int(value) for value in divisions)


def set_resolution(run_dir: Path, divisions: tuple[int, int, int]) -> None:
    """Scale the divisions of both beams in the m4 blockMesh template.

    The tutorial builds its blockMeshDict from system/blockMeshDict.m4, where
    the divisions of the upper and lower beam are held in two m4 macros.
    """
    block_mesh_dict = run_dir / "system" / "blockMeshDict.m4"
    if not block_mesh_dict.is_file():
        raise RuntimeError(f"no blockMeshDict.m4 found in {run_dir}/system")
    pattern = re.compile(
        r"(^\s*m4_define\(cylinder(?:Upper|Lower),\s*)\d+\s+\d+\s+\d+(\s*\))",
        re.MULTILINE,
    )
    text, count = pattern.subn(
        rf"\g<1>{divisions[0]} {divisions[1]} {divisions[2]}\g<2>",
        block_mesh_dict.read_text(),
    )
    if count != 2:
        raise RuntimeError(
            f"could not set both beam divisions in {block_mesh_dict}"
        )
    block_mesh_dict.write_text(text)


def set_friction(run_dir: Path, friction_coefficient: float) -> None:
    """Set the Coulomb friction coefficient of the contact boundary condition.

    A zero coefficient is applied as a frictionless contact model, matching the
    way the benchmark study in solid-benchmarks sets up its variants.
    """
    displacement_field = run_dir / "0" / "DD"
    text = displacement_field.read_text()
    if friction_coefficient > 0.0:
        text, count = re.subn(
            r"frictionCoeff\s+[0-9.eE+-]+\s*;",
            f"frictionCoeff   {friction_coefficient};",
            text,
        )
    else:
        text, count = re.subn(
            r"frictionContactModel\s+\w+\s*;",
            "frictionContactModel  frictionless;",
            text,
        )
    if count != 1:
        raise RuntimeError(f"could not set the friction model in {displacement_field}")
    displacement_field.write_text(text)


def read_reference_curve(path: Path) -> list[tuple[float, float]]:
    """Return the digitised curve as (displacement in mm, force in N)."""
    points = []
    for line in path.read_text().splitlines():
        if line.strip().startswith("#") or not line.strip():
            continue
        displacement, force = (float(value) for value in line.split()[:2])
        points.append((displacement, force))
    if len(points) < 2:
        raise RuntimeError(f"fewer than two reference points in {path}")
    return sorted(points)


def read_force_history(
    run_dir: Path, patch: str
) -> tuple[list[float], list[float], list[float]]:
    """Return (displacement, force_x, force_y) from the solidForces history."""
    candidates = sorted(
        (run_dir / "postProcessing").glob(f"**/solidForces{patch}.dat")
    )
    if not candidates:
        raise RuntimeError(f"no solidForces{patch}.dat found in {run_dir}")
    displacement, force_x, force_y = [], [], []
    for line in candidates[-1].read_text().splitlines():
        fields = line.split()
        if len(fields) < 3 or line.strip().startswith("#"):
            continue
        displacement.append(float(fields[0]))
        force_x.append(float(fields[1]))
        force_y.append(float(fields[2]))
    if len(displacement) < 2:
        raise RuntimeError(f"fewer than two rows in {candidates[-1]}")
    return displacement, force_x, force_y


def interpolate(points: list[float], values: list[float], point: float) -> float:
    if point <= points[0]:
        return values[0]
    if point >= points[-1]:
        return values[-1]
    for index in range(1, len(points)):
        if point <= points[index]:
            span = points[index] - points[index - 1]
            weight = (point - points[index - 1]) / span
            return values[index - 1] + weight * (values[index] - values[index - 1])
    return values[-1]


def sample_grid(displacement_mm: float, samples: int) -> list[float]:
    """Uniform displacement stations used to compare successive meshes."""
    return [displacement_mm * index / samples for index in range(1, samples + 1)]


def resample(
    displacement: list[float],
    force_x: list[float],
    force_y: list[float],
    grid: list[float],
) -> list[float]:
    """Both force components on the common grid, as one flat profile."""
    return [interpolate(displacement, force_x, station) for station in grid] + [
        interpolate(displacement, force_y, station) for station in grid
    ]


def compare_with_reference(
    displacement: list[float],
    force: list[float],
    reference: list[tuple[float, float]],
) -> tuple[float, float, float, float]:
    """Normalised RMS error, peak error, computed peak and reference peak.

    The reaction force passes through zero at the start and the end of the
    sliding, so the errors are normalised by the peak reference force rather
    than pointwise.
    """
    reference_peak = max(abs(value) for _, value in reference)
    errors = [
        interpolate(displacement, force, station) - value
        for station, value in reference
    ]
    computed_peak = max(abs(value) for value in force)
    return (
        math.sqrt(sum(error**2 for error in errors) / len(errors)) / reference_peak,
        abs(computed_peak - reference_peak) / reference_peak,
        computed_peak,
        reference_peak,
    )


def count_cells(run_dir: Path) -> int:
    match = re.search(
        r"^\s*nCells:\s*(\d+)\s*$",
        (run_dir / "log.blockMesh").read_text(errors="replace"),
        re.MULTILINE,
    )
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
            f"{label} did not complete; see {solver_log}. A diverged contact "
            "iteration is the usual cause."
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
    reference_x: list[tuple[float, float]],
    reference_y: list[tuple[float, float]],
    reuse: bool,
) -> dict:
    run_dir = WORK_DIR / variant_name / f"mesh{level}"
    solver_log = run_dir / "log.solids4Foam"
    completed_case = (
        reuse
        and solver_log.exists()
        and re.search(r"^End\s*$", solver_log.read_text(), re.MULTILINE)
    )

    divisions = divisions_for_level(reference, level)
    if completed_case:
        print(f"Reusing {variant_name} mesh{level} in {run_dir}")
    else:
        if run_dir.exists():
            shutil.rmtree(run_dir)
        shutil.copytree(CASE_DIR, run_dir, ignore=ignored, symlinks=True)
        set_resolution(run_dir, divisions)
        set_friction(run_dir, float(variant["friction_coefficient"]))
        print(
            f"Running {variant_name} mesh{level} "
            f"({divisions[0]} {divisions[1]} {divisions[2]} per beam) in {run_dir}",
            flush=True,
        )
        allverify_log = run_dir / "log.Allverify"
        with allverify_log.open("w") as log:
            completed = subprocess.run(
                ["./Allrun"],
                cwd=run_dir,
                stdout=log,
                stderr=subprocess.STDOUT,
                check=False,
            )
        if completed.returncode:
            raise RuntimeError(
                f"{variant_name} mesh{level} failed; see {allverify_log}"
            )
        if "only runs in foam-extend" in allverify_log.read_text(errors="replace"):
            raise RuntimeError(
                "the curvedBeams tutorial only runs in foam-extend; source a "
                "foam-extend environment before running this study"
            )
        if not solver_log.exists():
            raise RuntimeError(f"solver log was not created in {run_dir}")

    check_solver_log(solver_log, f"{variant_name} mesh{level}")

    displacement, force_x, force_y = read_force_history(
        run_dir, reference["force_patch"]
    )
    grid = sample_grid(
        float(reference["displacement_mm"]), int(reference["samples"])
    )
    profile = resample(displacement, force_x, force_y, grid)

    (POST_DIR / "profiles").mkdir(parents=True, exist_ok=True)
    history = POST_DIR / "profiles" / f"{variant_name}_mesh{level}_reaction.txt"
    history.write_text(
        "# displacement_mm force_X_N force_Y_N\n"
        + "".join(
            f"{station:.10g} {value_x:.10g} {value_y:.10g}\n"
            for station, value_x, value_y in zip(displacement, force_x, force_y)
        )
    )

    rms_x, peak_error_x, peak_x, reference_peak_x = compare_with_reference(
        displacement, force_x, reference_x
    )
    rms_y, peak_error_y, peak_y, reference_peak_y = compare_with_reference(
        displacement, force_y, reference_y
    )

    cell_count = count_cells(run_dir)
    clock_matches = re.findall(
        r"ClockTime\s*=\s*([0-9.eE+-]+)", solver_log.read_text(errors="replace")
    )
    row: dict = {
        "variant": variant_name,
        "friction_coefficient": float(variant["friction_coefficient"]),
        "level": level,
        "divisions": "x".join(str(value) for value in divisions[:2]),
        "cells": cell_count,
        "effective_spacing_mm": cell_count ** (-1.0 / 2.0),
        "clock_time_s": float(clock_matches[-1]) if clock_matches else math.nan,
        "peak_force_x_n": peak_x,
        "peak_force_y_n": peak_y,
        "reference_peak_force_x_n": reference_peak_x,
        "reference_peak_force_y_n": reference_peak_y,
        "reference_rms_relative_error": max(rms_x, rms_y),
        "reference_peak_relative_error": max(peak_error_x, peak_error_y),
    }
    numeric_values = (
        value for key, value in row.items() if key not in {"variant", "divisions"}
    )
    if not all(math.isfinite(float(value)) for value in numeric_values):
        raise RuntimeError(f"non-finite result extracted from {run_dir}")
    return {"row": row, "profile": profile}


def rms_change(first: list[float], second: list[float]) -> float:
    """RMS difference between two force histories on the same stations."""
    return math.sqrt(
        sum((left - right) ** 2 for left, right in zip(first, second)) / len(first)
    )


def net_order(results: list[dict], metric: str) -> float:
    finite = [row for row in results if math.isfinite(float(row[metric]))]
    if len(finite) < 2:
        return math.nan
    coarse = float(finite[0][metric])
    fine = float(finite[-1][metric])
    spacing_ratio = float(finite[0]["effective_spacing_mm"]) / float(
        finite[-1]["effective_spacing_mm"]
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
    grouped = {
        name: [row for row in results if row["variant"] == name]
        for name in variant_names
    }

    with (POST_DIR / "mesh_convergence.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=results[0].keys())
        writer.writeheader()
        writer.writerows(results)

    rms_tolerance = float(acceptance["reference_rms_relative_error_tolerance"])
    peak_tolerance = float(acceptance["reference_peak_relative_error_tolerance"])
    minimum_order = float(acceptance["minimum_net_order"])

    # The published curves are digitised from figures, so the difference from
    # them stops falling once the discretisation error drops below the
    # digitisation uncertainty. Convergence is therefore measured by the change
    # in the computed force history between successive meshes, and the
    # published curves are used to check the accuracy of the finest mesh.
    orders = {
        name: net_order(grouped[name][1:], "force_rms_change_n")
        for name in variant_names
    }
    passed = all(
        math.isfinite(float(row[metric])) and float(row[metric]) > 0.0
        for row in results
        for metric in ("peak_force_x_n", "peak_force_y_n",
                       "reference_rms_relative_error")
    )
    if not quick:
        for name in variant_names:
            changes = [
                float(row["force_rms_change_n"])
                for row in grouped[name][1:]
                if math.isfinite(float(row["force_rms_change_n"]))
            ]
            passed = (
                passed
                and float(grouped[name][-1]["reference_rms_relative_error"])
                <= rms_tolerance
                and float(grouped[name][-1]["reference_peak_relative_error"])
                <= peak_tolerance
                and len(changes) >= 2
                and all(right < left for left, right in zip(changes, changes[1:]))
                and orders[name] > minimum_order
            )

    lines = [
        "# curvedBeams verification summary",
        "",
        f"- Mode: {'quick smoke test' if quick else 'full verification'}",
        f"- Variants: {', '.join(variant_names)}",
        "- Reference: Neto et al. (2016), reaction force on the "
        f"{reference['force_patch']} patch",
    ]
    for name in variant_names:
        variant_results = grouped[name]
        levels = ", ".join(str(row["level"]) for row in variant_results)
        finest = variant_results[-1]
        lines.extend(
            [
                "",
                f"## {name}",
                "",
                f"- Levels: {levels}",
                f"- Finest mesh: {int(finest['cells'])} cells",
                "- Finest reference RMS relative error: "
                f"{float(finest['reference_rms_relative_error']):.4f}",
                "- Finest reference peak relative error: "
                f"{float(finest['reference_peak_relative_error']):.4f}",
                "- Finest peak reaction force: "
                f"{float(finest['peak_force_x_n']):.6g} N in x "
                f"(reference {float(finest['reference_peak_force_x_n']):.6g} N), "
                f"{float(finest['peak_force_y_n']):.6g} N in y "
                f"(reference {float(finest['reference_peak_force_y_n']):.6g} N)",
                "- Force history RMS change between successive levels (N): "
                + ", ".join(
                    f"{float(row['force_rms_change_n']):.4g}"
                    for row in variant_results[1:]
                ),
                f"- Net order of the force change: {orders[name]:.3f}",
            ]
        )
    lines.extend(["", f"- Result: {'PASS' if passed else 'FAIL'}"])
    summary = "\n".join(lines) + "\n"
    (POST_DIR / "verification_summary.md").write_text(summary)
    print(summary)
    print(f"CSV: {POST_DIR / 'mesh_convergence.csv'}")
    return passed


def create_plot(
    variant_names: list[str], levels: list[int], reference: dict
) -> None:
    plot_script = SCRIPT_DIR / "plotReactionForces.gnuplot"
    if shutil.which("gnuplot") is None or not plot_script.is_file():
        return
    variant = reference["variants"][variant_names[0]]
    output = POST_DIR / "curvedBeams_reaction.pdf"
    completed = subprocess.run(
        [
            "gnuplot",
            "-e",
            f"profiles='{POST_DIR / 'profiles'}'; "
            f"referenceX='{REFERENCE_DIR / variant['reference_x']}'; "
            f"referenceY='{REFERENCE_DIR / variant['reference_y']}'; "
            f"variant='{variant_names[0]}'; "
            f"levels='{' '.join(str(level) for level in levels)}'; "
            f"displacement={float(reference['displacement_mm'])}; "
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
    unknown = [level for level in levels if str(level) not in reference["level_divisions"]]
    if unknown:
        available = ", ".join(sorted(reference["level_divisions"]))
        raise SystemExit(f"undefined level(s): {unknown}; available levels: {available}")

    missing = [
        command
        for command in ("m4", "blockMesh", "solids4Foam")
        if shutil.which(command) is None
    ]
    if missing:
        raise SystemExit(f"required command(s) not found: {', '.join(missing)}")

    WORK_DIR.mkdir(parents=True, exist_ok=True)
    POST_DIR.mkdir(parents=True, exist_ok=True)
    results: list[dict] = []
    failures = 0
    for name in variant_names:
        variant = all_variants[name]
        reference_x = read_reference_curve(REFERENCE_DIR / variant["reference_x"])
        reference_y = read_reference_curve(REFERENCE_DIR / variant["reference_y"])
        previous_profile: list[float] | None = None
        for level in levels:
            try:
                outcome = run_level(
                    name,
                    variant,
                    level,
                    reference,
                    reference_x,
                    reference_y,
                    args.reuse,
                )
            except RuntimeError as error:
                print(f"ERROR: {error}", file=sys.stderr)
                failures += 1
                previous_profile = None
                if not args.keep_going:
                    return 1
                continue
            row = outcome["row"]
            row["force_rms_change_n"] = (
                rms_change(previous_profile, outcome["profile"])
                if previous_profile is not None
                else math.nan
            )
            previous_profile = outcome["profile"]
            results.append(row)
    if any(
        sum(row["variant"] == name for row in results) < 2 for name in variant_names
    ):
        print("ERROR: fewer than two levels completed for a variant", file=sys.stderr)
        return 1
    passed = write_results(results, variant_names, reference, args.quick)
    create_plot(variant_names, levels, reference)
    return 0 if passed and failures == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
