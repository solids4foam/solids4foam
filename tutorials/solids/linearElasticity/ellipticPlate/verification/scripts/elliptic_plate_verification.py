#!/usr/bin/env python3
"""Run the ellipticPlate mesh-convergence verification study.

The study refines the tutorial blockMesh through the Demirdzic et al. (1997)
mesh family, samples the equivalent (von Mises) stress along the line
r = 2.1 m, z = 0.3 m, and checks convergence towards the published solution.
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
REFERENCE_FILE = (
    REFERENCE_DIR / "ellipticPlate_verification_references.json"
)

SAMPLE_DICT_NAME = "sampleVerificationLine"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run the ellipticPlate verification variants"
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


def divisions_for_level(base_divisions: list[int], level: int) -> tuple[int, int, int]:
    factor = 2 ** (level - 1)
    return tuple(value * factor for value in base_divisions)


def set_resolution(run_dir: Path, divisions: tuple[int, int, int]) -> None:
    """Scale the block divisions, retaining the tutorial edge grading.

    The tutorial keeps a separate mesh dictionary for OpenFOAM and for
    foam-extend and its Allrun links whichever suits the loaded environment,
    so both are refined here.
    """
    pattern = re.compile(
        r"(^\s*hex\s*\([^)]*\)\s*)\(\s*\d+\s+\d+\s+\d+\s*\)", re.MULTILINE
    )
    updated = 0
    for name in ("blockMeshDict.openfoam", "blockMeshDict.foamextend"):
        block_mesh_dict = run_dir / "system" / name
        if not block_mesh_dict.is_file():
            continue
        text, count = pattern.subn(
            rf"\g<1>({divisions[0]} {divisions[1]} {divisions[2]})",
            block_mesh_dict.read_text(),
            count=1,
        )
        if count != 1:
            raise RuntimeError(
                f"could not set block divisions in {block_mesh_dict}"
            )
        block_mesh_dict.write_text(text)
        updated += 1
    if not updated:
        raise RuntimeError(f"no blockMeshDict variant found in {run_dir}/system")


def write_sample_dict(run_dir: Path, sample_line: dict) -> None:
    """Write a sets function object sampling sigmaEq along the benchmark line."""
    radius = float(sample_line["radius_m"])
    z_coordinate = float(sample_line["z_m"])
    samples = int(sample_line["samples"])
    limit = 0.5 * math.pi
    points = []
    for index in range(samples):
        theta = limit * index / (samples - 1)
        # Keep the end points marginally inside the symmetry planes so the
        # mesh search always succeeds.
        theta = min(max(theta, 1.0e-9), limit - 1.0e-9)
        points.append(
            (radius * math.cos(theta), radius * math.sin(theta), z_coordinate)
        )

    entries = "\n".join(
        f"            ({point[0]:.12g} {point[1]:.12g} {point[2]:.12g})"
        for point in points
    )
    (run_dir / "system" / SAMPLE_DICT_NAME).write_text(
        "// Sampling of sigmaEq along the line r = "
        f"{radius} m, z = {z_coordinate} m, written by the verification driver.\n"
        "type            sets;\n"
        "libs            (sampling);\n"
        "interpolationScheme cellPoint;\n"
        "setFormat       raw;\n"
        "fields          (sigmaEq);\n"
        "\n"
        "sets\n"
        "{\n"
        "    arc\n"
        "    {\n"
        "        type    cloud;\n"
        "        ordered yes;\n"
        "        axis    xyz;\n"
        "        points\n"
        "        (\n"
        f"{entries}\n"
        "        );\n"
        "    }\n"
        "}\n"
    )


def sample_angles(sample_line: dict) -> list[float]:
    samples = int(sample_line["samples"])
    return [90.0 * index / (samples - 1) for index in range(samples)]


def read_sampled_stress(run_dir: Path) -> list[float]:
    candidates = sorted(
        (run_dir / "postProcessing" / SAMPLE_DICT_NAME).glob("*/arc_sigmaEq.xy")
    )
    if not candidates:
        raise RuntimeError(f"no sampled sigmaEq data found in {run_dir}")
    values = []
    for line in candidates[-1].read_text().splitlines():
        fields = line.split()
        if len(fields) == 4:
            values.append(float(fields[3]))
    if not values:
        raise RuntimeError(f"no numeric rows in {candidates[-1]}")
    return values


def read_reference_curve(path: Path) -> list[tuple[float, float]]:
    """Return the digitised curve as (degrees, Pa), sorted by angle."""
    points = []
    for line in path.read_text().splitlines():
        if line.strip().startswith("#") or not line.strip():
            continue
        angle, stress_mpa = (float(value) for value in line.split()[:2])
        points.append((angle, 1.0e6 * stress_mpa))
    if len(points) < 2:
        raise RuntimeError(f"fewer than two reference points in {path}")
    return sorted(points)


def interpolate(angles: list[float], values: list[float], angle: float) -> float:
    if angle <= angles[0]:
        return values[0]
    if angle >= angles[-1]:
        return values[-1]
    for index in range(1, len(angles)):
        if angle <= angles[index]:
            span = angles[index] - angles[index - 1]
            weight = (angle - angles[index - 1]) / span
            return values[index - 1] + weight * (values[index] - values[index - 1])
    return values[-1]


def compare_with_reference(
    angles: list[float],
    values: list[float],
    reference: list[tuple[float, float]],
) -> dict[str, float]:
    relative_errors = []
    for angle, reference_value in reference:
        computed = interpolate(angles, values, angle)
        relative_errors.append(abs(computed - reference_value) / abs(reference_value))
    peak_reference = max(reference, key=lambda point: point[0])
    peak_computed = interpolate(angles, values, peak_reference[0])
    return {
        "reference_rms_relative_error": math.sqrt(
            sum(error**2 for error in relative_errors) / len(relative_errors)
        ),
        "reference_max_relative_error": max(relative_errors),
        "sigma_eq_peak_pa": peak_computed,
        "reference_peak_pa": peak_reference[1],
        "reference_peak_relative_error": abs(peak_computed - peak_reference[1])
        / abs(peak_reference[1]),
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
    reference_curve: list[tuple[float, float]],
    reuse: bool,
) -> dict[str, float | int | str]:
    run_dir = WORK_DIR / variant_name / f"mesh{level}"
    solver_log = run_dir / "log.solids4Foam"
    completed_case = (
        reuse
        and solver_log.exists()
        and re.search(r"^End\s*$", solver_log.read_text(), re.MULTILINE)
    )

    divisions = divisions_for_level(reference["base_divisions"], level)
    if completed_case:
        print(f"Reusing {variant_name} mesh{level} in {run_dir}")
    else:
        if run_dir.exists():
            shutil.rmtree(run_dir)
        shutil.copytree(CASE_DIR, run_dir, ignore=ignored, symlinks=True)
        set_resolution(run_dir, divisions)
        print(
            f"Running {variant_name} mesh{level} "
            f"({divisions[0]} {divisions[1]} {divisions[2]}) in {run_dir}",
            flush=True,
        )
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
            raise RuntimeError(f"solver log was not created in {run_dir}")

    check_solver_log(solver_log, f"{variant_name} mesh{level}")

    write_sample_dict(run_dir, reference["sample_line"])
    with (run_dir / "log.postProcess").open("w") as log:
        completed = subprocess.run(
            ["postProcess", "-func", SAMPLE_DICT_NAME, "-latestTime"],
            cwd=run_dir,
            stdout=log,
            stderr=subprocess.STDOUT,
            check=False,
        )
    if completed.returncode:
        raise RuntimeError(
            f"sampling failed for {variant_name} mesh{level}; "
            f"see {run_dir / 'log.postProcess'}"
        )

    values = read_sampled_stress(run_dir)
    angles = sample_angles(reference["sample_line"])
    if len(values) != len(angles):
        raise RuntimeError(
            f"expected {len(angles)} sampled values, found {len(values)} "
            f"in {run_dir}"
        )
    (POST_DIR / "profiles").mkdir(parents=True, exist_ok=True)
    profile = POST_DIR / "profiles" / f"{variant_name}_mesh{level}_sigmaEq.txt"
    profile.write_text(
        "# angle_deg sigmaEq_Pa\n"
        + "".join(
            f"{angle:.10g} {value:.10g}\n" for angle, value in zip(angles, values)
        )
    )

    cell_count = count_cells(run_dir)
    log_text = solver_log.read_text()
    clock_matches = re.findall(r"ClockTime\s*=\s*([0-9.eE+-]+)", log_text)
    result: dict[str, float | int | str] = {
        "variant": variant_name,
        "approach": variant["approach"],
        "level": level,
        "divisions": "x".join(str(value) for value in divisions),
        "cells": cell_count,
        "effective_spacing_m": cell_count ** (-1.0 / 3.0),
        "clock_time_s": float(clock_matches[-1]) if clock_matches else math.nan,
    }
    result.update(compare_with_reference(angles, values, reference_curve))
    numeric_values = (
        value
        for key, value in result.items()
        if key not in {"variant", "approach", "divisions"}
    )
    if not all(math.isfinite(float(value)) for value in numeric_values):
        raise RuntimeError(f"non-finite result extracted from {run_dir}")
    return {"row": result, "profile": values}


def rms_change(first: list[float], second: list[float]) -> float:
    """RMS difference between two sampled profiles on the same sample points."""
    return math.sqrt(
        sum((left - right) ** 2 for left, right in zip(first, second)) / len(first)
    )


def net_order(results: list[dict[str, float | int | str]], metric: str) -> float:
    finite = [row for row in results if math.isfinite(float(row[metric]))]
    if len(finite) < 2:
        return math.nan
    coarse = float(finite[0][metric])
    fine = float(finite[-1][metric])
    spacing_ratio = float(finite[0]["effective_spacing_m"]) / float(
        finite[-1]["effective_spacing_m"]
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

    # The published curve is digitised from a figure, so the difference from
    # it stops falling once the discretisation error drops below the
    # digitisation uncertainty. Convergence is therefore measured by the
    # change in the computed profile between successive meshes, and the
    # published curve is used to check the accuracy of the finest mesh.
    orders = {
        name: net_order(grouped[name][1:], "profile_rms_change_pa")
        for name in variant_names
    }
    passed = all(
        math.isfinite(float(row[metric])) and float(row[metric]) > 0.0
        for row in results
        for metric in ("reference_rms_relative_error", "sigma_eq_peak_pa")
    )
    if not quick:
        for name in variant_names:
            changes = [
                float(row["profile_rms_change_pa"])
                for row in grouped[name][1:]
                if math.isfinite(float(row["profile_rms_change_pa"]))
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
        "# ellipticPlate verification summary",
        "",
        f"- Mode: {'quick smoke test' if quick else 'full verification'}",
        f"- Variants: {', '.join(variant_names)}",
        "- Reference: Demirdzic et al. (1997), sigmaEq on r = 2.1 m, z = 0.3 m",
    ]
    for name in variant_names:
        variant_results = grouped[name]
        levels = ", ".join(str(row["level"]) for row in variant_results)
        lines.extend(
            [
                "",
                f"## {name}",
                "",
                f"- Levels: {levels}",
                f"- Finest mesh: {int(variant_results[-1]['cells'])} cells",
                "- Finest reference RMS relative error: "
                f"{float(variant_results[-1]['reference_rms_relative_error']):.4f}",
                "- Finest reference peak relative error: "
                f"{float(variant_results[-1]['reference_peak_relative_error']):.4f}",
                "- Finest peak sigmaEq: "
                f"{float(variant_results[-1]['sigma_eq_peak_pa']):.6g} Pa "
                f"(reference {float(variant_results[-1]['reference_peak_pa']):.6g} Pa)",
                "- Profile RMS change between successive levels (Pa): "
                + ", ".join(
                    f"{float(row['profile_rms_change_pa']):.4g}"
                    for row in variant_results[1:]
                ),
                f"- Net order of the profile change: {orders[name]:.3f}",
            ]
        )
    lines.extend(["", f"- Result: {'PASS' if passed else 'FAIL'}"])
    summary = "\n".join(lines) + "\n"
    (POST_DIR / "verification_summary.md").write_text(summary)
    print(summary)
    print(f"CSV: {POST_DIR / 'mesh_convergence.csv'}")
    return passed


def create_plot(variant_names: list[str], levels: list[int]) -> None:
    plot_script = SCRIPT_DIR / "plotSigmaEq.gnuplot"
    if shutil.which("gnuplot") is None or not plot_script.is_file():
        return
    output = POST_DIR / "ellipticPlate_sigmaEq.pdf"
    completed = subprocess.run(
        [
            "gnuplot",
            "-e",
            f"profiles='{POST_DIR / 'profiles'}'; "
            f"reference='{REFERENCE_DIR / 'demirdzic_sigmaEq.txt'}'; "
            f"variant='{variant_names[0]}'; "
            f"levels='{' '.join(str(level) for level in levels)}'; "
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
    reference_curve = read_reference_curve(
        REFERENCE_DIR / reference["reference_curve"]
    )
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
    results: list[dict[str, float | int | str]] = []
    failures = 0
    for name in variant_names:
        previous_profile: list[float] | None = None
        for level in levels:
            try:
                outcome = run_level(
                    name,
                    all_variants[name],
                    level,
                    reference,
                    reference_curve,
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
            row["profile_rms_change_pa"] = (
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
    create_plot(variant_names, levels)
    return 0 if passed and failures == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
