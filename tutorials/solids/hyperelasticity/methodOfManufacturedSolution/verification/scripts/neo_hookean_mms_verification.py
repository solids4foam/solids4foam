#!/usr/bin/env python3
"""Run the neo-Hookean manufactured-solution convergence studies.

The spatial study runs the quasi-static ``steady`` case on a sequence of
meshes and measures the order of the displacement and stress errors with
respect to the exact solution. The temporal study runs the ``transient`` case
on one mesh with a sequence of time-step sizes; because the spatial error is
fixed, the temporal order is measured from the differences between successive
time-step sizes (Richardson), and the errors with respect to the exact
solution are reported alongside.
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
TUTORIAL_DIR = VERIFY_DIR.parent
WORK_DIR = VERIFY_DIR / "work"
POST_DIR = VERIFY_DIR / "postProcessing"
REFERENCE_FILE = (
    VERIFY_DIR / "reference" / "neoHookeanMMS_verification_references.json"
)
SPATIAL_METRICS = [
    "displacement_l2_m",
    "displacement_linf_m",
    "stress_l2_pa",
    "stress_linf_pa",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run the methodOfManufacturedSolution verification studies"
    )
    parser.add_argument(
        "--study",
        choices=["spatial", "temporal", "all"],
        default="all",
        help="which study to run (default: all)",
    )
    parser.add_argument(
        "--variants",
        help="comma-separated variant names of the selected study",
    )
    parser.add_argument(
        "--levels",
        help=(
            "comma-separated cells per direction (spatial) or numbers of "
            "time steps (temporal)"
        ),
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
    """Files not copied from the tutorial directory into a work case."""
    ignored_names = {
        "verification",
        "regressionTests",
        "postProcessing",
        "0.tmp",
        "lnInclude",
        "log.wmake",
    }
    ignored_names.update(name for name in names if name.startswith("processor"))
    ignored_names.update(name for name in names if name.startswith("log."))
    ignored_names.update(name for name in names if name.endswith(".msh"))

    directory_path = Path(directory)
    if directory_path.name in {"steady", "transient"}:
        ignored_names.update(
            name
            for name in names
            if name != "0" and re.fullmatch(r"[0-9]+(?:\.[0-9]+)?(?:e-?[0-9]+)?", name)
        )
    if directory_path.name == "constant":
        ignored_names.add("polyMesh")
    if directory_path.name == "Make":
        ignored_names.update(name for name in names if name not in {"files", "options"})

    return ignored_names.intersection(names)


def copy_tutorial(run_dir: Path) -> None:
    if run_dir.exists():
        shutil.rmtree(run_dir)
    shutil.copytree(TUTORIAL_DIR, run_dir, ignore=ignored, symlinks=True)


def set_divisions(case_dir: Path, divisions: int, domain_length: float) -> None:
    block_mesh_dict = case_dir / "system" / "blockMeshDict"
    text = block_mesh_dict.read_text()
    pattern = re.compile(r"(hex\s*\([^)]*\)\s*)\(\s*\d+\s+\d+\s+\d+\s*\)")
    text, count = pattern.subn(
        rf"\g<1>({divisions} {divisions} {divisions})", text, count=1
    )
    if count != 1:
        raise RuntimeError(f"could not set divisions in {block_mesh_dict}")
    block_mesh_dict.write_text(text)

    spacing = domain_length / divisions
    (case_dir / "gmsh" / "meshSpacing.geo").write_text(
        "// Mesh spacing for the Gmsh meshes\n" f"dx = {spacing:.10g};\n"
    )


def set_time_step(case_dir: Path, end_time: float, time_steps: int) -> None:
    control_dict = case_dir / "system" / "controlDict"
    text = control_dict.read_text()
    delta_t = end_time / time_steps
    text, count = re.subn(
        r"(?m)^(deltaT\s+)[^;]+;", rf"\g<1>{delta_t:.15g};", text
    )
    if count != 1:
        raise RuntimeError(f"could not set deltaT in {control_dict}")
    text, count = re.subn(
        r"(?m)^(endTime\s+)[^;]+;", rf"\g<1>{end_time:.15g};", text
    )
    if count != 1:
        raise RuntimeError(f"could not set endTime in {control_dict}")
    control_dict.write_text(text)


def set_time_function(case_dir: Path, time_function: str) -> None:
    properties = case_dir / "constant" / "neoHookeanMMSProperties"
    text = properties.read_text()
    text, count = re.subn(
        r"(?m)^(timeFunction\s+)\w+;", rf"\g<1>{time_function};", text
    )
    if count != 1:
        raise RuntimeError(f"could not set timeFunction in {properties}")
    properties.write_text(text)


def run_case(case_dir: Path, command: list[str]) -> Path:
    solver_log = case_dir / "log.solids4Foam"
    print(f"Running {' '.join(command)} in {case_dir}", flush=True)
    with (case_dir / "log.Allverify").open("w") as log:
        completed = subprocess.run(
            command,
            cwd=case_dir,
            stdout=log,
            stderr=subprocess.STDOUT,
            check=False,
        )
    if completed.returncode:
        raise RuntimeError(
            f"{' '.join(command)} failed in {case_dir}; "
            f"see {case_dir / 'log.Allverify'}"
        )
    if not solver_log.exists():
        raise RuntimeError(f"solver log was not created in {case_dir}")
    return solver_log


def check_log(solver_log: Path) -> str:
    log_text = solver_log.read_text()
    if not re.search(r"^End\s*$", log_text, re.MULTILINE):
        raise RuntimeError(f"incomplete solver log: {solver_log}")
    if "DIVERGED_" in log_text:
        raise RuntimeError(f"PETSc failed to converge: {solver_log}")
    if re.search(r"did not converge", log_text):
        raise RuntimeError(f"momentum equation did not converge: {solver_log}")
    return log_text


def extract_norms(log_text: str, marker: str) -> tuple[float, float]:
    """Return the final mean-L2 and L-infinity error norms after a marker."""
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


def extract_exact_errors(log_text: str) -> dict[str, float]:
    displacement_l2, displacement_linf = extract_norms(
        log_text, "Writing DDifference field"
    )
    stress_l2, stress_linf = extract_norms(log_text, "Writing sigmaDifference field")
    return {
        "displacement_l2_m": displacement_l2,
        "displacement_linf_m": displacement_linf,
        "stress_l2_pa": stress_l2,
        "stress_linf_pa": stress_linf,
    }


def count_cells(case_dir: Path) -> int:
    completed = subprocess.run(
        ["checkMesh", "-constant"],
        cwd=case_dir,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )
    if completed.returncode:
        raise RuntimeError(f"checkMesh failed in {case_dir}")
    match = re.search(r"^\s*cells:\s*(\d+)\s*$", completed.stdout, re.MULTILINE)
    if not match:
        raise RuntimeError(f"could not extract cell count in {case_dir}")
    return int(match.group(1))


def clock_time(log_text: str) -> float:
    matches = re.findall(r"ClockTime\s*=\s*([0-9.eE+-]+)", log_text)
    return float(matches[-1]) if matches else math.nan


def read_internal_vectors(field_file: Path) -> list[tuple[float, float, float]]:
    """Read the internal values of an ascii volVectorField."""
    text = field_file.read_text()
    match = re.search(
        r"internalField\s+nonuniform\s+List<vector>\s*(\d+)\s*\(", text
    )
    if not match:
        raise RuntimeError(f"could not find a nonuniform internalField in {field_file}")
    count = int(match.group(1))
    start = match.end()
    values: list[tuple[float, float, float]] = []
    pattern = re.compile(r"\(\s*([^()\s]+)\s+([^()\s]+)\s+([^()\s]+)\s*\)")
    for item in pattern.finditer(text, start):
        values.append((float(item.group(1)), float(item.group(2)), float(item.group(3))))
        if len(values) == count:
            break
    if len(values) != count:
        raise RuntimeError(f"expected {count} values in {field_file}, read {len(values)}")
    return values


def field_difference_norms(
    coarse: list[tuple[float, float, float]],
    fine: list[tuple[float, float, float]],
) -> tuple[float, float]:
    if len(coarse) != len(fine):
        raise RuntimeError("fields of different size cannot be compared")
    sum_sqr = 0.0
    max_mag = 0.0
    for a, b in zip(coarse, fine):
        mag_sqr = sum((x - y) ** 2 for x, y in zip(a, b))
        sum_sqr += mag_sqr
        max_mag = max(max_mag, math.sqrt(mag_sqr))
    return math.sqrt(sum_sqr / len(coarse)), max_mag


def final_time_directory(case_dir: Path, end_time: float) -> Path:
    candidates = []
    for item in case_dir.iterdir():
        if not item.is_dir():
            continue
        try:
            value = float(item.name)
        except ValueError:
            continue
        candidates.append((abs(value - end_time), item))
    if not candidates:
        raise RuntimeError(f"no time directories in {case_dir}")
    candidates.sort(key=lambda pair: pair[0])
    if candidates[0][0] > 1e-6 * max(end_time, 1.0):
        raise RuntimeError(f"no time directory at t = {end_time} in {case_dir}")
    return candidates[0][1]


def finite_positive(value: float) -> bool:
    return math.isfinite(value) and value > 0


def fmt(value: float, spec: str) -> str:
    """Format a number, printing '-' for values that are not finite."""
    return format(value, spec) if math.isfinite(value) else "-"


# -----------------------------------------------------------------------------
# Spatial study
# -----------------------------------------------------------------------------

def run_spatial_level(
    name: str, variant: dict, divisions: int, reference: dict, reuse: bool
) -> dict:
    run_dir = WORK_DIR / "spatial" / name / f"n{divisions}"
    case_dir = run_dir / "steady"
    solver_log = case_dir / "log.solids4Foam"
    completed = (
        reuse
        and solver_log.exists()
        and re.search(r"^End\s*$", solver_log.read_text(), re.MULTILINE)
    )
    if completed:
        print(f"Reusing {name} n={divisions} in {run_dir}")
    else:
        copy_tutorial(run_dir)
        set_divisions(case_dir, divisions, float(reference["domain_length_m"]))
        run_case(case_dir, ["./Allrun", variant["approach"], variant["mesh"]])

    log_text = check_log(solver_log)
    result: dict = {
        "variant": name,
        "approach": variant["approach"],
        "mesh": variant["mesh"],
        "divisions": divisions,
        "cells": count_cells(case_dir),
    }
    result["effective_spacing_m"] = (
        float(reference["domain_length_m"]) ** 3 / result["cells"]
    ) ** (1.0 / 3.0)
    result.update(extract_exact_errors(log_text))
    result["clock_time_s"] = clock_time(log_text)
    if not all(finite_positive(result[metric]) for metric in SPATIAL_METRICS):
        raise RuntimeError(f"non-finite or zero error extracted from {solver_log}")
    return result


def net_order(results: list[dict], metric: str) -> float:
    coarse = float(results[0][metric])
    fine = float(results[-1][metric])
    ratio = float(results[0]["effective_spacing_m"]) / float(
        results[-1]["effective_spacing_m"]
    )
    if coarse <= 0 or fine <= 0 or ratio <= 1:
        return math.nan
    return math.log(coarse / fine) / math.log(ratio)


def summarise_spatial(
    results: list[dict], variant_names: list[str], reference: dict, quick: bool
) -> tuple[list[str], bool]:
    grouped = {
        name: [row for row in results if row["variant"] == name]
        for name in variant_names
    }
    passed = True
    lines = ["## Spatial convergence (steady case)", ""]
    lines.append(
        "Errors with respect to the exact solution at the end of the load ramp."
    )
    for name in variant_names:
        rows = grouped[name]
        variant = reference["variants"][name]
        lines.extend(
            [
                "",
                f"### {name}",
                "",
                f"- Divisions: {', '.join(str(row['divisions']) for row in rows)}",
                f"- Finest mesh: {int(rows[-1]['cells'])} cells",
            ]
        )
        for metric in SPATIAL_METRICS:
            minimum = variant["minimum_net_order"][metric.split("_")[0]]
            order = net_order(rows, metric)
            metric_passed = all(finite_positive(float(row[metric])) for row in rows)
            if not quick:
                metric_passed = (
                    metric_passed
                    and float(rows[-1][metric]) < float(rows[0][metric])
                    and math.isfinite(order)
                    and order >= minimum
                )
            passed = passed and metric_passed
            status = "PASS" if metric_passed else "FAIL"
            if quick:
                status += " (order not checked)"
            lines.append(
                f"- {metric}: finest {float(rows[-1][metric]):.8g}, "
                f"net order {fmt(order, '.3f')}, minimum {minimum:g}: {status}"
            )
        lines.append("")
        lines.append("| divisions | cells | disp L2 [m] | disp Linf [m] | stress L2 [Pa] | stress Linf [Pa] | clock [s] |")
        lines.append("| --- | --- | --- | --- | --- | --- | --- |")
        for row in rows:
            lines.append(
                f"| {row['divisions']} | {row['cells']} | "
                f"{row['displacement_l2_m']:.4e} | {row['displacement_linf_m']:.4e} | "
                f"{row['stress_l2_pa']:.4e} | {row['stress_linf_pa']:.4e} | "
                f"{fmt(row['clock_time_s'], '.0f')} |"
            )
    return lines, passed


# -----------------------------------------------------------------------------
# Temporal study
# -----------------------------------------------------------------------------

def run_temporal_level(
    name: str, variant: dict, time_steps: int, reference: dict, reuse: bool
) -> dict:
    run_dir = WORK_DIR / "temporal" / name / f"nt{time_steps}"
    case_dir = run_dir / "transient"
    solver_log = case_dir / "log.solids4Foam"
    end_time = float(reference["end_time_s"])
    completed = (
        reuse
        and solver_log.exists()
        and re.search(r"^End\s*$", solver_log.read_text(), re.MULTILINE)
    )
    if completed:
        print(f"Reusing {name} nt={time_steps} in {run_dir}")
    else:
        copy_tutorial(run_dir)
        set_divisions(case_dir, int(reference["mesh_divisions"]), 1.0)
        set_time_step(case_dir, end_time, time_steps)
        set_time_function(case_dir, variant["time_function"])
        run_case(
            case_dir,
            ["./Allrun", variant["approach"], variant["scheme"], variant["mesh"]],
        )

    log_text = check_log(solver_log)
    result: dict = {
        "variant": name,
        "approach": variant["approach"],
        "scheme": variant["scheme"],
        "time_function": variant["time_function"],
        "mesh": variant["mesh"],
        "divisions": int(reference["mesh_divisions"]),
        "time_steps": time_steps,
        "delta_t_s": end_time / time_steps,
    }
    result.update(extract_exact_errors(log_text))
    result["clock_time_s"] = clock_time(log_text)
    result["_D"] = read_internal_vectors(
        final_time_directory(case_dir, end_time) / "D"
    )
    if not all(finite_positive(result[metric]) for metric in SPATIAL_METRICS):
        raise RuntimeError(f"non-finite or zero error extracted from {solver_log}")
    return result


def add_richardson_columns(rows: list[dict]) -> None:
    """Differences between successive time-step sizes and the observed order.

    With e_k = ||D(dt_k) - D(dt_{k+1})||, the observed order between levels
    k and k+1 is log(e_k/e_{k+1})/log(dt_k/dt_{k+1}).
    """
    for k, row in enumerate(rows):
        row["difference_l2_m"] = math.nan
        row["difference_linf_m"] = math.nan
        row["observed_order_l2"] = math.nan
        row["observed_order_linf"] = math.nan
        if k + 1 < len(rows):
            l2, linf = field_difference_norms(row["_D"], rows[k + 1]["_D"])
            row["difference_l2_m"] = l2
            row["difference_linf_m"] = linf
    for k in range(len(rows) - 2):
        ratio = rows[k]["delta_t_s"] / rows[k + 1]["delta_t_s"]
        for norm in ("l2", "linf"):
            coarse = rows[k][f"difference_{norm}_m"]
            fine = rows[k + 1][f"difference_{norm}_m"]
            if coarse > 0 and fine > 0 and ratio > 1:
                rows[k]["observed_order_" + norm] = math.log(coarse / fine) / math.log(ratio)


def summarise_temporal(
    results: list[dict], variant_names: list[str], reference: dict, quick: bool
) -> tuple[list[str], bool]:
    grouped = {
        name: [row for row in results if row["variant"] == name]
        for name in variant_names
    }
    passed = True
    lines = ["## Temporal convergence (transient case)", ""]
    lines.append(
        f"Fixed {reference['mesh_divisions']}^3 hexahedral mesh, end time "
        f"{reference['end_time_s']} s. The observed order is measured from the "
        "differences between the displacement fields of successive time-step "
        "sizes; the errors with respect to the exact solution include the "
        "fixed spatial error."
    )
    for name in variant_names:
        rows = grouped[name]
        variant = reference["variants"][name]
        add_richardson_columns(rows)
        orders = [row["observed_order_l2"] for row in rows if math.isfinite(row["observed_order_l2"])]
        final_order = orders[-1] if orders else math.nan
        lines.extend(
            [
                "",
                f"### {name}",
                "",
                f"- Time steps: {', '.join(str(row['time_steps']) for row in rows)}",
                f"- Time function: {variant['time_function']}",
                f"- Nominal order: {variant['nominal_order']}",
            ]
        )
        variant_passed = all(
            finite_positive(float(row[metric])) for row in rows for metric in SPATIAL_METRICS
        )
        if quick or len(rows) < 3:
            status = "PASS (order not checked)" if variant_passed else "FAIL"
        else:
            variant_passed = (
                variant_passed
                and math.isfinite(final_order)
                and final_order >= variant["minimum_order"]
            )
            status = "PASS" if variant_passed else "FAIL"
        passed = passed and variant_passed
        lines.append(
            f"- Observed order (finest three levels, L2): {fmt(final_order, '.3f')}, "
            f"minimum {variant['minimum_order']:g}: {status}"
        )
        lines.append("")
        lines.append("| steps | dt [s] | disp L2 vs exact [m] | stress L2 vs exact [Pa] | diff to next dt, L2 [m] | diff to next dt, Linf [m] | order L2 | order Linf | clock [s] |")
        lines.append("| --- | --- | --- | --- | --- | --- | --- | --- | --- |")
        for row in rows:
            lines.append(
                f"| {row['time_steps']} | {row['delta_t_s']:.6g} | "
                f"{row['displacement_l2_m']:.4e} | {row['stress_l2_pa']:.4e} | "
                f"{fmt(row['difference_l2_m'], '.4e')} | "
                f"{fmt(row['difference_linf_m'], '.4e')} | "
                f"{fmt(row['observed_order_l2'], '.3f')} | "
                f"{fmt(row['observed_order_linf'], '.3f')} | "
                f"{fmt(row['clock_time_s'], '.0f')} |"
            )
    return lines, passed


# -----------------------------------------------------------------------------
# Driver
# -----------------------------------------------------------------------------

def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        return
    fieldnames = [key for key in rows[0].keys() if not key.startswith("_")]
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def select_variants(args: argparse.Namespace, reference: dict) -> list[str]:
    all_variants = reference["variants"]
    if args.variants:
        names = [value.strip() for value in args.variants.split(",")]
    else:
        names = list(reference["default_variants"])
    invalid = [name for name in names if name not in all_variants]
    if not names or invalid or len(set(names)) != len(names):
        raise SystemExit(
            f"variants must be a unique selection from: {', '.join(all_variants)}"
        )
    return names


def select_levels(args: argparse.Namespace, default: list[int]) -> list[int]:
    if args.levels:
        levels = [int(value) for value in args.levels.split(",")]
    else:
        levels = [int(value) for value in default]
    if args.quick:
        levels = levels[:2]
    if len(levels) < 2 or any(level <= 0 for level in levels):
        raise SystemExit("at least two positive levels are required")
    if any(right <= left for left, right in zip(levels, levels[1:])):
        raise SystemExit("levels must be strictly increasing")
    return levels


def main() -> int:
    args = parse_args()
    reference = json.loads(REFERENCE_FILE.read_text())
    studies = ["spatial", "temporal"] if args.study == "all" else [args.study]
    if args.study == "all" and (args.variants or args.levels):
        raise SystemExit("--variants and --levels require --study spatial or temporal")

    required = ["checkMesh", "solids4Foam"]
    missing = [command for command in required if shutil.which(command) is None]
    if missing:
        raise SystemExit(f"required command(s) not found: {', '.join(missing)}")

    WORK_DIR.mkdir(parents=True, exist_ok=True)
    POST_DIR.mkdir(parents=True, exist_ok=True)
    mode = "quick smoke test" if args.quick else "full verification"
    overall_passed = True
    failures = 0

    for study in studies:
        study_reference = reference[study]
        variant_names = select_variants(args, study_reference)
        if study == "spatial":
            levels = select_levels(args, study_reference["default_levels"])
            runner = run_spatial_level
        else:
            levels = select_levels(args, study_reference["default_time_steps"])
            runner = run_temporal_level
        tet_or_perturbed = any(
            study_reference["variants"][name]["mesh"] in {"tet", "distHex"}
            for name in variant_names
        )
        if tet_or_perturbed:
            extra = ["perturbMeshPoints"]
            if any(study_reference["variants"][name]["mesh"] == "tet" for name in variant_names):
                extra.extend(["gmsh", "gmshToFoam", "createPatch"])
            missing = [command for command in extra if shutil.which(command) is None]
            if missing:
                raise SystemExit(f"required command(s) not found: {', '.join(missing)}")

        results: list[dict] = []
        for name in variant_names:
            for level in levels:
                try:
                    results.append(
                        runner(name, study_reference["variants"][name], level, study_reference, args.reuse)
                    )
                except RuntimeError as error:
                    print(f"ERROR: {error}", file=sys.stderr)
                    failures += 1
                    if not args.keep_going:
                        return 1
        completed_names = [
            name for name in variant_names
            if sum(row["variant"] == name for row in results) >= 2
        ]
        if len(completed_names) < len(variant_names):
            print("ERROR: fewer than two levels completed for a variant", file=sys.stderr)
            overall_passed = False
        if not completed_names:
            continue
        if study == "spatial":
            lines, passed = summarise_spatial(results, completed_names, study_reference, args.quick)
            write_csv(POST_DIR / "spatial_convergence.csv", results)
        else:
            lines, passed = summarise_temporal(results, completed_names, study_reference, args.quick)
            write_csv(POST_DIR / "temporal_convergence.csv", results)
        overall_passed = overall_passed and passed
        lines.extend(["", f"- Mode: {mode}", f"- Result: {'PASS' if passed else 'FAIL'}"])
        (POST_DIR / f"{study}_summary.md").write_text("\n".join(lines) + "\n")

    # Combine the latest summary of each study, so that running one study
    # does not discard the other's results
    summary = ["# Manufactured-solution verification summary"]
    for study in ("spatial", "temporal"):
        study_summary = POST_DIR / f"{study}_summary.md"
        if study_summary.exists():
            summary.extend(["", study_summary.read_text().rstrip("\n")])
    text = "\n".join(summary) + "\n"
    (POST_DIR / "verification_summary.md").write_text(text)
    print(text)
    return 0 if overall_passed and not failures else 1


if __name__ == "__main__":
    raise SystemExit(main())
