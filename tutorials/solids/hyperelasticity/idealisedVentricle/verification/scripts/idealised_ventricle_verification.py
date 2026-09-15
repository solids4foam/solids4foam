#!/usr/bin/env python3
"""Run the idealisedVentricle mesh-convergence verification study.

The study refines the tutorial blockMesh and rotational extrusion together,
samples the deformed mid-wall line of Land et al. (2015) Problem 2, and checks
that the solution converges under uniform refinement.  Land et al. publish the
benchmark as a comparison between codes rather than as a closed-form solution,
so the acceptance criterion here is self-convergence of the mid-wall line and
of the apex position, with an independent solids4foam solution retained as an
informational diagnostic.
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
    REFERENCE_DIR / "idealisedVentricle_verification_references.json"
)

SAMPLE_DICT_NAME = "sampleVerificationMidLine"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run the idealisedVentricle verification variants"
    )
    parser.add_argument(
        "--variants",
        help="comma-separated variant names",
    )
    parser.add_argument(
        "--levels",
        help="comma-separated mesh levels (1 is the tutorial mesh)",
    )
    parser.add_argument(
        "--quick",
        action="store_true",
        help="run the first two levels as a smoke test",
    )
    parser.add_argument(
        "--cores",
        default="auto",
        help="MPI ranks: auto (default), one positive integer for every "
             "level, or one comma-separated value per requested level",
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


def parse_cores(requested: str, levels: list[int]) -> dict[int, int]:
    """Map each requested level to its MPI rank count.

    `auto` scales the ranks with the cell count on a workstation-sized
    machine. A single value applies to every level. A comma-separated list is
    matched to the requested levels in order, which is what a scheduler
    allocation with one fixed rank count per mesh needs.
    """
    if requested == "auto":
        default_cores = {1: 1, 2: 8, 3: 16, 4: 16}
        return {level: default_cores.get(level, 16) for level in levels}

    values = [value.strip() for value in requested.split(",")]
    if not all(value.isdecimal() and int(value) >= 1 for value in values):
        raise SystemExit(
            "--cores must be auto, a positive integer, or one positive "
            "integer per requested level"
        )
    if len(values) == 1:
        return {level: int(values[0]) for level in levels}
    if len(values) != len(levels):
        raise SystemExit(
            f"--cores lists {len(values)} values for {len(levels)} levels"
        )
    return {level: int(value) for level, value in zip(levels, values)}


def set_resolution(run_dir: Path, level: int, reference: dict) -> tuple[int, int, int, int]:
    """Refine the block divisions and the extrusion layers together."""
    factor = 2 ** (level - 1)
    base_divisions = reference["base_divisions"]
    divisions = (
        base_divisions[0] * factor,
        base_divisions[1] * factor,
        base_divisions[2],
    )
    layers = int(reference["base_extrude_layers"]) * factor

    block_mesh_dict = run_dir / "caseOptions" / "petsc" / "system" / "blockMeshDict"
    text = block_mesh_dict.read_text()
    text, count = re.subn(
        r"(^\s*hex\s*\([^)]*\)\s*)\(\s*\d+\s+\d+\s+\d+\s*\)",
        rf"\g<1>({divisions[0]} {divisions[1]} {divisions[2]})",
        text,
        count=1,
        flags=re.MULTILINE,
    )
    if count != 1:
        raise RuntimeError(f"could not set block divisions in {block_mesh_dict}")
    block_mesh_dict.write_text(text)

    extrude_mesh_dict = run_dir / "caseOptions" / "petsc" / "system" / "extrudeMeshDict"
    text = extrude_mesh_dict.read_text()
    text, count = re.subn(
        r"(^\s*nLayers\s+)\d+\s*;",
        rf"\g<1>{layers};",
        text,
        count=1,
        flags=re.MULTILINE,
    )
    if count != 1:
        raise RuntimeError(f"could not set nLayers in {extrude_mesh_dict}")
    extrude_mesh_dict.write_text(text)

    return divisions[0], divisions[1], divisions[2], layers


def set_cores(run_dir: Path, cores: int) -> None:
    decompose_par_dict = (
        run_dir / "caseOptions" / "petsc" / "system" / "decomposeParDict"
    )
    text = decompose_par_dict.read_text()
    text, count = re.subn(
        r"(^\s*numberOfSubdomains\s+)\d+\s*;",
        rf"\g<1>{cores};",
        text,
        count=1,
        flags=re.MULTILINE,
    )
    if count != 1:
        raise RuntimeError(
            f"could not set numberOfSubdomains in {decompose_par_dict}"
        )
    decompose_par_dict.write_text(text)


def mid_line_points(mid_line: dict) -> list[tuple[float, float, float]]:
    """Return the undeformed mid-wall line of Land et al. (2015) Problem 2."""
    short_axis = mid_line["short_axis_radii_m"]
    long_axis = mid_line["long_axis_radii_m"]
    samples = int(mid_line["samples"])
    r_short = 0.5 * (short_axis[0] + short_axis[1])
    r_long = 0.5 * (long_axis[0] + long_axis[1])
    u_max = -math.acos(float(mid_line["truncation_ratio"]))
    u_min = -math.pi
    step = (u_max - u_min) / (samples - 1)
    return [
        (
            r_short * math.sin(u_min + index * step),
            0.0,
            r_long * math.cos(u_min + index * step),
        )
        for index in range(samples)
    ]


def write_sample_dict(run_dir: Path, points: list[tuple[float, float, float]]) -> None:
    entries = "\n".join(
        f"            ({point[0]:.12g} {point[1]:.12g} {point[2]:.12g})"
        for point in points
    )
    (run_dir / "system" / SAMPLE_DICT_NAME).write_text(
        "// Sampling of D along the undeformed mid-wall line, written by the\n"
        "// verification driver.\n"
        "type            sets;\n"
        "libs            (sampling);\n"
        "interpolationScheme cellPoint;\n"
        "setFormat       raw;\n"
        "fields          (D);\n"
        "\n"
        "sets\n"
        "{\n"
        "    midLine\n"
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


def read_deformed_mid_line(
    run_dir: Path, points: list[tuple[float, float, float]]
) -> list[tuple[float, float, float]]:
    candidates = sorted(
        (run_dir / "postProcessing" / SAMPLE_DICT_NAME).glob("*/midLine_D.xy")
    )
    if not candidates:
        raise RuntimeError(f"no sampled displacement data found in {run_dir}")
    deformed = []
    for line in candidates[-1].read_text().splitlines():
        fields = line.split()
        if len(fields) != 6:
            continue
        values = [float(field) for field in fields]
        deformed.append((values[0] + values[3], values[1] + values[4],
                         values[2] + values[5]))
    if len(deformed) != len(points):
        raise RuntimeError(
            f"expected {len(points)} sampled points, found {len(deformed)} "
            f"in {candidates[-1]}"
        )
    return deformed


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
    cores: int,
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
        print(f"Reusing {variant_name} mesh{level} in {run_dir}")
        divisions = set_resolution(run_dir, level, reference)
    else:
        if run_dir.exists():
            shutil.rmtree(run_dir)
        shutil.copytree(CASE_DIR, run_dir, ignore=ignored, symlinks=True)
        divisions = set_resolution(run_dir, level, reference)
        set_cores(run_dir, cores)
        command = ["./Allrun", variant["approach"]]
        if cores > 1:
            command.append("parallel")
        print(
            f"Running {variant_name} mesh{level} "
            f"({divisions[0]} {divisions[1]} {divisions[2]}, "
            f"{divisions[3]} layers) on {cores} core(s) in {run_dir}",
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

    points = mid_line_points(reference["mid_line"])
    write_sample_dict(run_dir, points)
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

    deformed = read_deformed_mid_line(run_dir, points)
    (POST_DIR / "profiles").mkdir(parents=True, exist_ok=True)
    profile = POST_DIR / "profiles" / f"{variant_name}_mesh{level}_midLine.txt"
    profile.write_text(
        "# x y z (deformed mid-wall line, m)\n"
        + "".join(f"{x:.10g} {y:.10g} {z:.10g}\n" for x, y, z in deformed)
    )

    cell_count = count_cells(run_dir)
    log_text = solver_log.read_text()
    clock_matches = re.findall(r"ClockTime\s*=\s*([0-9.eE+-]+)", log_text)
    return {
        "row": {
            "variant": variant_name,
            "approach": variant["approach"],
            "level": level,
            "divisions": "x".join(str(value) for value in divisions[:3]),
            "extrude_layers": divisions[3],
            "cells": cell_count,
            "cores": cores,
            "effective_spacing_m": cell_count ** (-1.0 / 3.0),
            "apex_z_m": min(point[2] for point in deformed),
            "clock_time_s": float(clock_matches[-1]) if clock_matches else math.nan,
        },
        "mid_line": deformed,
    }


def rms_distance(
    first: list[tuple[float, float, float]],
    second: list[tuple[float, float, float]],
) -> float:
    total = sum(
        sum((a - b) ** 2 for a, b in zip(left, right))
        for left, right in zip(first, second)
    )
    return math.sqrt(total / len(first))


def net_order(
    rows: list[dict], metric: str
) -> float:
    finite = [row for row in rows if math.isfinite(float(row[metric]))]
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
    rows: list[dict],
    variant_names: list[str],
    reference: dict,
    quick: bool,
) -> bool:
    POST_DIR.mkdir(parents=True, exist_ok=True)
    diagnostic = reference["diagnostics"]["apex_z_m"]
    minimum_order = float(reference["acceptance"]["minimum_net_order"])

    for row in rows:
        row["diagnostic_apex_z_m"] = diagnostic["value"]
        row["diagnostic_apex_relative_difference"] = abs(
            float(row["apex_z_m"]) - float(diagnostic["value"])
        ) / abs(float(diagnostic["value"]))

    grouped = {
        name: [row for row in rows if row["variant"] == name]
        for name in variant_names
    }

    with (POST_DIR / "mesh_convergence.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)

    passed = all(
        math.isfinite(float(row["apex_z_m"])) and float(row["apex_z_m"]) < 0.0
        for row in rows
    )
    orders = {}
    for name in variant_names:
        variant_rows = grouped[name]
        orders[name] = net_order(variant_rows[1:], "mid_line_rms_change_m")
        if not quick:
            changes = [
                float(row["mid_line_rms_change_m"])
                for row in variant_rows[1:]
                if math.isfinite(float(row["mid_line_rms_change_m"]))
            ]
            passed = (
                passed
                and len(changes) >= 2
                and all(right < left for left, right in zip(changes, changes[1:]))
                and orders[name] > minimum_order
            )

    lines = [
        "# idealisedVentricle verification summary",
        "",
        f"- Mode: {'quick smoke test' if quick else 'full verification'}",
        f"- Variants: {', '.join(variant_names)}",
        "- Study: self-convergence of the Land et al. (2015) Problem 2"
        " mid-wall line under uniform refinement",
    ]
    for name in variant_names:
        variant_rows = grouped[name]
        levels = ", ".join(str(row["level"]) for row in variant_rows)
        lines.extend(
            [
                "",
                f"## {name}",
                "",
                f"- Levels: {levels}",
                f"- Finest mesh: {int(variant_rows[-1]['cells'])} cells",
                "- Finest apex position: "
                f"{1.0e3 * float(variant_rows[-1]['apex_z_m']):.4f} mm",
                "- Independent solution diagnostic: "
                f"{1.0e3 * float(diagnostic['value']):.4f} mm "
                "(difference "
                f"{100.0 * float(variant_rows[-1]['diagnostic_apex_relative_difference']):.2f}%)",
                "- Mid-wall line RMS change between successive levels (m): "
                + ", ".join(
                    f"{float(row['mid_line_rms_change_m']):.4g}"
                    for row in variant_rows[1:]
                ),
                f"- Net order of the mid-wall line change: {orders[name]:.3f}",
            ]
        )
    lines.extend(["", f"- Result: {'PASS' if passed else 'FAIL'}"])
    summary = "\n".join(lines) + "\n"
    (POST_DIR / "verification_summary.md").write_text(summary)
    print(summary)
    print(f"CSV: {POST_DIR / 'mesh_convergence.csv'}")
    return passed


def create_plot(variant_names: list[str], levels: list[int]) -> None:
    plot_script = SCRIPT_DIR / "plotMidLine.gnuplot"
    if shutil.which("gnuplot") is None or not plot_script.is_file():
        return
    output = POST_DIR / "idealisedVentricle_midLine.pdf"
    completed = subprocess.run(
        [
            "gnuplot",
            "-e",
            f"profiles='{POST_DIR / 'profiles'}'; "
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
    cores = parse_cores(args.cores, levels)

    required = ["blockMesh", "checkMesh", "createPatch", "extrudeMesh",
                "postProcess", "solids4Foam"]
    if any(value > 1 for value in cores.values()):
        required.extend(["decomposePar", "reconstructPar"])
    missing = [command for command in required if shutil.which(command) is None]
    if missing:
        raise SystemExit(f"required command(s) not found: {', '.join(missing)}")

    WORK_DIR.mkdir(parents=True, exist_ok=True)
    POST_DIR.mkdir(parents=True, exist_ok=True)
    rows: list[dict] = []
    failures = 0
    for name in variant_names:
        previous_mid_line: list[tuple[float, float, float]] | None = None
        for level in levels:
            try:
                outcome = run_level(
                    name,
                    all_variants[name],
                    level,
                    reference,
                    cores[level],
                    args.reuse,
                )
            except RuntimeError as error:
                print(f"ERROR: {error}", file=sys.stderr)
                failures += 1
                previous_mid_line = None
                if not args.keep_going:
                    return 1
                continue
            row = outcome["row"]
            row["mid_line_rms_change_m"] = (
                rms_distance(previous_mid_line, outcome["mid_line"])
                if previous_mid_line is not None
                else math.nan
            )
            previous_mid_line = outcome["mid_line"]
            rows.append(row)
    if any(
        sum(row["variant"] == name for row in rows) < 2 for name in variant_names
    ):
        print("ERROR: fewer than two levels completed for a variant", file=sys.stderr)
        return 1
    passed = write_results(rows, variant_names, reference, args.quick)
    create_plot(variant_names, levels)
    return 0 if passed and failures == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
