#!/usr/bin/env python3
"""Run the cavityFlexibleBottom mesh-convergence verification study.

The study uniformly refines the tutorial fluid and solid meshes, runs each
level to a steady FSI response, and compares the steady vertical displacement
at the monitor point and the steady vertical interface force with the mesh
study of Tukovic et al. (2018).
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
    REFERENCE_DIR / "cavityFlexibleBottom_verification_references.json"
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run the cavityFlexibleBottom verification variants"
    )
    parser.add_argument("--variants", help="comma-separated variant names")
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
        "--end-time",
        type=float,
        help="override the pseudo-transient end time",
    )
    parser.add_argument(
        "--delta-t",
        type=float,
        help="override the time step, which the finest mesh needs",
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
    if directory_path.name in {"fluid", "solid"} and directory_path.parent.name == "constant":
        ignored_names.add("polyMesh")

    return ignored_names.intersection(names)


def refine_mesh(path: Path, factor: int) -> int:
    """Scale the in-plane block divisions, leaving the one-cell span alone."""
    pattern = re.compile(r"(hex\s*\([^)]*\)\s*)\(\s*(\d+)\s+(\d+)\s+(\d+)\s*\)")

    def scale(match: re.Match[str]) -> str:
        first, second, third = (int(match.group(index)) for index in (2, 3, 4))
        return (
            f"{match.group(1)}({first * factor} {second * factor} {third})"
        )

    text, count = pattern.subn(scale, path.read_text())
    if not count:
        raise RuntimeError(f"no block divisions found in {path}")
    path.write_text(text)
    return count


def set_control(run_dir: Path, key: str, value: float) -> None:
    """Set one controlDict entry in every region that declares it."""
    updated = 0
    for control_dict in (
        run_dir / "system" / "controlDict",
        run_dir / "system" / "solid" / "controlDict",
        run_dir / "system" / "fluid" / "controlDict",
    ):
        if not control_dict.is_file():
            continue
        text, count = re.subn(
            rf"(^\s*{re.escape(key)}\s+)[^;]+;",
            rf"\g<1>{value:.10g};",
            control_dict.read_text(),
            count=1,
            flags=re.MULTILINE,
        )
        if count == 1:
            control_dict.write_text(text)
            updated += 1
    if not updated:
        raise RuntimeError(f"no controlDict in {run_dir} declares {key}")


def set_coupling(run_dir: Path, coupling: str) -> None:
    fsi_properties = run_dir / "constant" / "fsiProperties"
    text, count = re.subn(
        r"(^fluidSolidInterface\s+)\w+\s*;",
        rf"\g<1>{coupling};",
        fsi_properties.read_text(),
        count=1,
        flags=re.MULTILINE,
    )
    if count != 1:
        raise RuntimeError(
            f"could not set fluidSolidInterface in {fsi_properties}"
        )
    fsi_properties.write_text(text)


def read_reference_curve(path: Path) -> list[tuple[float, float]]:
    """Return the digitised curve as (DeltaX, value), sorted by DeltaX."""
    points = []
    for line in path.read_text().splitlines():
        if line.strip().startswith("#") or not line.strip():
            continue
        delta_x, value = (float(field) for field in line.split(",")[:2])
        points.append((delta_x, value))
    if len(points) < 2:
        raise RuntimeError(f"fewer than two reference points in {path}")
    return sorted(points)


def reference_at(curve: list[tuple[float, float]], delta_x: float) -> float:
    """Nearest published mesh spacing, which the levels are chosen to match."""
    return min(curve, key=lambda point: abs(point[0] - delta_x))[1]


def numeric_rows(path: Path) -> list[list[float]]:
    rows = []
    for line in path.read_text(errors="replace").splitlines():
        fields = line.replace("(", " ").replace(")", " ").split()
        if not fields or fields[0].startswith("#"):
            continue
        try:
            rows.append([float(field) for field in fields])
        except ValueError:
            continue
    return rows


def read_displacement_history(run_dir: Path) -> list[tuple[float, float]]:
    candidates = sorted(
        run_dir.glob("postProcessing/**/solidPointDisplacement_*.dat")
    )
    if not candidates:
        raise RuntimeError(f"no solidPointDisplacement data in {run_dir}")
    rows = numeric_rows(candidates[-1])
    if not rows:
        raise RuntimeError(f"no numeric rows in {candidates[-1]}")
    # Columns are time, then the displacement vector components.
    return [(row[0], row[2]) for row in rows if len(row) >= 3]


def force_columns(path: Path) -> list[int]:
    """Return the column indices whose sum is the vertical force.

    The forces function object labels its columns in the last comment line.
    Recent OpenFOAM versions write a `total_y` column; older ones write only
    the pressure and viscous contributions, which then have to be added.
    """
    header = ""
    for line in path.read_text(errors="replace").splitlines():
        if line.startswith("#"):
            header = line
        else:
            break
    names = header.lstrip("#").split()
    if "total_y" in names:
        return [names.index("total_y")]
    contributions = [
        names.index(name)
        for name in ("pressure_y", "viscous_y")
        if name in names
    ]
    if contributions:
        return contributions
    # Without a usable header, assume the leading vector is the total force.
    return [2]


def read_force_history(run_dir: Path) -> list[tuple[float, float]]:
    candidates = sorted(run_dir.glob("postProcessing/**/force*.dat"))
    if not candidates:
        raise RuntimeError(f"no force data in {run_dir}")
    columns = force_columns(candidates[-1])
    rows = numeric_rows(candidates[-1])
    if not rows:
        raise RuntimeError(f"no numeric rows in {candidates[-1]}")
    history = [
        (row[0], sum(row[column] for column in columns))
        for row in rows
        if len(row) > max(columns)
    ]
    if not history:
        raise RuntimeError(f"could not parse force vectors in {candidates[-1]}")
    return history


def steady_value(
    history: list[tuple[float, float]], window_fraction: float
) -> tuple[float, float]:
    """Return the final value and its relative spread over the closing window."""
    if not history:
        raise RuntimeError("empty history")
    end_time = history[-1][0]
    start_time = end_time * (1.0 - window_fraction)
    window = [value for time, value in history if time >= start_time]
    if len(window) < 2:
        window = [value for _, value in history[-2:]]
    final = history[-1][1]
    scale = max(abs(final), 1.0e-30)
    return final, (max(window) - min(window)) / scale


def count_cells(run_dir: Path, region: str) -> int:
    completed = subprocess.run(
        ["checkMesh", "-constant", "-region", region],
        cwd=run_dir,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )
    if completed.returncode:
        raise RuntimeError(f"checkMesh failed for region {region} in {run_dir}")
    match = re.search(r"^\s*cells:\s*(\d+)\s*$", completed.stdout, re.MULTILINE)
    if not match:
        raise RuntimeError(f"could not extract {region} cell count in {run_dir}")
    return int(match.group(1))


def check_solver_log(solver_log: Path, end_time: float, label: str) -> None:
    """Reject a run that crashed or stopped short of the end time.

    The tutorial Allrun returns zero even when the solver aborts, so a
    diverged case would otherwise be read as a converged steady result.
    """
    text = solver_log.read_text(errors="replace")
    if re.search(r"FOAM FATAL|FOAM aborting|^ERROR$|\[stack trace\]", text,
                 re.MULTILINE):
        raise RuntimeError(
            f"{label} diverged or aborted; see {solver_log}. Partitioned "
            "coupling is unstable on the finest mesh at the tutorial time "
            "step: rerun that level with a smaller --delta-t."
        )
    times = [float(match) for match in
             re.findall(r"^Time = ([0-9.eE+-]+)", text, re.MULTILINE)]
    if not times or times[-1] < end_time - 1.0e-8:
        reached = times[-1] if times else 0.0
        raise RuntimeError(
            f"{label} stopped at t = {reached:g} before the end time "
            f"{end_time:g}; see {solver_log}"
        )


def run_level(
    variant_name: str,
    variant: dict,
    level: int,
    reference: dict,
    end_time: float,
    delta_t: float,
    reuse: bool,
) -> dict:
    run_dir = WORK_DIR / variant_name / f"mesh{level}"
    solver_log = run_dir / "log.solids4Foam"
    completed_case = (
        reuse
        and solver_log.exists()
        and re.search(r"^End\s*$", solver_log.read_text(), re.MULTILINE)
    )

    factor = 2 ** (level - 1)
    delta_x = float(reference["base_delta_x_m"]) / factor

    if completed_case:
        print(f"Reusing {variant_name} mesh{level} in {run_dir}")
    else:
        if run_dir.exists():
            shutil.rmtree(run_dir)
        shutil.copytree(CASE_DIR, run_dir, ignore=ignored, symlinks=True)
        refine_mesh(run_dir / "system" / "fluid" / "blockMeshDict", factor)
        refine_mesh(run_dir / "system" / "solid" / "blockMeshDict", factor)
        set_control(run_dir, "endTime", end_time)
        set_control(run_dir, "deltaT", delta_t)
        set_coupling(run_dir, variant["coupling"])
        print(
            f"Running {variant_name} mesh{level} "
            f"(DeltaX = {delta_x:g} m) in {run_dir}",
            flush=True,
        )
        with (run_dir / "log.Allverify").open("w") as log:
            completed = subprocess.run(
                ["./Allrun"],
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

    check_solver_log(solver_log, end_time, f"{variant_name} mesh{level}")

    window_fraction = float(reference["steady_state"]["window_fraction"])
    displacement, displacement_spread = steady_value(
        read_displacement_history(run_dir), window_fraction
    )
    force, force_spread = steady_value(
        read_force_history(run_dir), window_fraction
    )

    log_text = solver_log.read_text()
    clock_matches = re.findall(r"ClockTime\s*=\s*([0-9.eE+-]+)", log_text)
    return {
        "variant": variant_name,
        "coupling": variant["coupling"],
        "level": level,
        "delta_x_m": delta_x,
        "cells_fluid": count_cells(run_dir, "fluid"),
        "cells_solid": count_cells(run_dir, "solid"),
        "end_time_s": end_time,
        "delta_t_s": delta_t,
        "uy_m": displacement,
        "uy_steady_spread": displacement_spread,
        "fy_n": force,
        "fy_steady_spread": force_spread,
        "clock_time_s": float(clock_matches[-1]) if clock_matches else math.nan,
    }


def net_order(rows: list[dict], metric: str) -> float:
    finite = [row for row in rows if math.isfinite(float(row[metric]))]
    if len(finite) < 2:
        return math.nan
    coarse = float(finite[0][metric])
    fine = float(finite[-1][metric])
    spacing_ratio = float(finite[0]["delta_x_m"]) / float(finite[-1]["delta_x_m"])
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
    acceptance = reference["acceptance"]
    steady_tolerance = float(reference["steady_state"]["relative_tolerance"])
    grouped = {
        name: [row for row in rows if row["variant"] == name]
        for name in variant_names
    }

    with (POST_DIR / "mesh_convergence.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)

    # A separate, stable three-column file for the plot, so the gnuplot script
    # does not depend on the column order of the CSV.
    (POST_DIR / "profiles").mkdir(parents=True, exist_ok=True)
    for name in variant_names:
        (POST_DIR / "profiles" / f"{name}_meshConvergence.txt").write_text(
            # Comma separated, matching the reference curves, so one gnuplot
            # separator setting serves every file in a plot command.
            "# deltaX_m, uy_m, fy_N\n"
            + "".join(
                f"{float(row['delta_x_m']):.10g}, {float(row['uy_m']):.10g}, "
                f"{float(row['fy_n']):.10g}\n"
                for row in grouped[name]
            )
        )

    # Every level must actually have reached a steady response, otherwise the
    # comparison with the published steady values is meaningless.
    passed = all(
        float(row["uy_steady_spread"]) <= steady_tolerance
        and float(row["fy_steady_spread"]) <= steady_tolerance
        for row in rows
    )
    orders = {
        name: {
            metric: net_order(grouped[name][1:], metric)
            for metric in ("uy_change_m", "fy_change_n")
        }
        for name in variant_names
    }
    if not quick:
        for name in variant_names:
            finest = grouped[name][-1]
            changes_uy = [
                float(row["uy_change_m"])
                for row in grouped[name][1:]
                if math.isfinite(float(row["uy_change_m"]))
            ]
            passed = (
                passed
                and float(finest["uy_relative_error"])
                <= float(acceptance["displacement_relative_error_tolerance"])
                and float(finest["fy_relative_error"])
                <= float(acceptance["force_relative_error_tolerance"])
                and len(changes_uy) >= 2
                and all(
                    right < left
                    for left, right in zip(changes_uy, changes_uy[1:])
                )
                and orders[name]["uy_change_m"]
                > float(acceptance["minimum_net_order"])
            )

    lines = [
        "# cavityFlexibleBottom verification summary",
        "",
        f"- Mode: {'quick smoke test' if quick else 'full verification'}",
        f"- Variants: {', '.join(variant_names)}",
        "- Reference: Tukovic et al. (2018) steady mesh study",
    ]
    for name in variant_names:
        variant_rows = grouped[name]
        finest = variant_rows[-1]
        lines.extend(
            [
                "",
                f"## {name}",
                "",
                "- Levels: "
                + ", ".join(str(row["level"]) for row in variant_rows),
                f"- Finest mesh: DeltaX = {float(finest['delta_x_m']):g} m, "
                f"{int(finest['cells_fluid'])} fluid and "
                f"{int(finest['cells_solid'])} solid cells",
                f"- Finest Uy: {float(finest['uy_m']):.6g} m "
                f"(reference {float(finest['reference_uy_m']):.6g} m, "
                f"{100.0 * float(finest['uy_relative_error']):.2f}%)",
                f"- Finest Fy: {float(finest['fy_n']):.6g} N "
                f"(reference {float(finest['reference_fy_n']):.6g} N, "
                f"{100.0 * float(finest['fy_relative_error']):.2f}%)",
                "- Uy change between successive levels (m): "
                + ", ".join(
                    f"{float(row['uy_change_m']):.4g}" for row in variant_rows[1:]
                ),
                "- Net order of the Uy change: "
                f"{orders[name]['uy_change_m']:.3f}",
            ]
        )
    lines.extend(["", f"- Result: {'PASS' if passed else 'FAIL'}"])
    summary = "\n".join(lines) + "\n"
    (POST_DIR / "verification_summary.md").write_text(summary)
    print(summary)
    print(f"CSV: {POST_DIR / 'mesh_convergence.csv'}")
    return passed


def create_plot(variant_names: list[str]) -> None:
    plot_script = SCRIPT_DIR / "plotMeshConvergence.gnuplot"
    if shutil.which("gnuplot") is None or not plot_script.is_file():
        return
    output = POST_DIR / "cavityFlexibleBottom_meshConvergence.pdf"
    completed = subprocess.run(
        [
            "gnuplot",
            "-e",
            f"profiles='{POST_DIR / 'profiles'}'; "
            f"variant='{variant_names[0]}'; "
            f"displacements='{REFERENCE_DIR / 'TukovicDisplacements.csv'}'; "
            f"forces='{REFERENCE_DIR / 'TukovicForces.csv'}'; "
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
    displacement_curve = read_reference_curve(
        REFERENCE_DIR / reference["displacement_curve"]
    )
    force_scale = float(reference["tutorial_thickness_m"]) / float(
        reference["force_reference_thickness_m"]
    )
    force_curve = [
        (delta_x, force_scale * value)
        for delta_x, value in read_reference_curve(
            REFERENCE_DIR / reference["force_curve"]
        )
    ]

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

    end_time = (
        args.end_time
        if args.end_time is not None
        else float(reference["steady_state"]["end_time"])
    )
    if end_time <= 0:
        raise SystemExit("--end-time must be positive")
    delta_t = (
        args.delta_t
        if args.delta_t is not None
        else float(reference["steady_state"]["delta_t"])
    )
    if delta_t <= 0 or delta_t > end_time:
        raise SystemExit("--delta-t must be positive and below the end time")

    missing = [
        command
        for command in ("blockMesh", "checkMesh", "solids4Foam")
        if shutil.which(command) is None
    ]
    if missing:
        raise SystemExit(f"required command(s) not found: {', '.join(missing)}")

    WORK_DIR.mkdir(parents=True, exist_ok=True)
    POST_DIR.mkdir(parents=True, exist_ok=True)
    rows: list[dict] = []
    failures = 0
    for name in variant_names:
        previous: dict | None = None
        for level in levels:
            try:
                row = run_level(
                    name,
                    all_variants[name],
                    level,
                    reference,
                    end_time,
                    delta_t,
                    args.reuse,
                )
            except RuntimeError as error:
                print(f"ERROR: {error}", file=sys.stderr)
                failures += 1
                previous = None
                if not args.keep_going:
                    return 1
                continue
            row["reference_uy_m"] = reference_at(
                displacement_curve, row["delta_x_m"]
            )
            row["reference_fy_n"] = reference_at(force_curve, row["delta_x_m"])
            row["uy_relative_error"] = abs(
                row["uy_m"] - row["reference_uy_m"]
            ) / abs(row["reference_uy_m"])
            row["fy_relative_error"] = abs(
                row["fy_n"] - row["reference_fy_n"]
            ) / abs(row["reference_fy_n"])
            row["uy_change_m"] = (
                abs(row["uy_m"] - previous["uy_m"]) if previous else math.nan
            )
            row["fy_change_n"] = (
                abs(row["fy_n"] - previous["fy_n"]) if previous else math.nan
            )
            previous = row
            rows.append(row)
    if any(
        sum(row["variant"] == name for row in rows) < 2 for name in variant_names
    ):
        print("ERROR: fewer than two levels completed for a variant", file=sys.stderr)
        return 1
    passed = write_results(rows, variant_names, reference, args.quick)
    create_plot(variant_names)
    return 0 if passed and failures == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
