#!/usr/bin/env python3
"""Run the multiMaterialPlate verification studies.

The thickness study runs the tutorial mesh for each mid-layer thickness,
compares the order of the free-edge singularity with the analytical order of
Bogy (1970), and plots the normal stress along the interfaces against the
boundary element results of Ioka et al. (2007). Those curves are digitised
from a figure, so they are drawn for a visual comparison only.

The mesh study coarsens and refines the tutorial mesh for one mid-layer
thickness and checks that the stress along the interfaces converges.
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
PROFILE_DIR = POST_DIR / "profiles"
REFERENCE_DIR = VERIFY_DIR / "reference"
REFERENCE_FILE = REFERENCE_DIR / "multiMaterialPlate_verification_references.json"

STUDIES = ("thickness", "mesh")

Profile = list[tuple[float, float]]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run the multiMaterialPlate verification studies"
    )
    parser.add_argument(
        "--studies",
        help="comma-separated studies: thickness, mesh (default: both)",
    )
    parser.add_argument("--variants", help="comma-separated variant names")
    parser.add_argument(
        "--thicknesses",
        help="comma-separated h/W values for the thickness study",
    )
    parser.add_argument(
        "--levels",
        help="comma-separated mesh levels for the mesh study (1 is the "
        "coarsest)",
    )
    parser.add_argument(
        "--quick",
        action="store_true",
        help="run the thinnest and thickest mid-layers and the two coarsest "
        "levels as a smoke test",
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
    ignored_names = {
        "verification",
        "regressionTests",
        "postProcessing",
        "images",
        "case.foam",
    }
    ignored_names.update(name for name in names if name.startswith("processor"))
    ignored_names.update(name for name in names if name.startswith("log."))
    ignored_names.update(name for name in names if name.endswith(".png"))

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


def scaled(value: int, shift: int, required: bool, what: str) -> int:
    """Scale a block division by 2^shift."""
    if shift >= 0:
        return value * 2**shift
    denominator = 2 ** (-shift)
    if value % denominator:
        if required:
            raise RuntimeError(
                f"{what} = {value} is not divisible by {denominator}: the "
                "level is not reachable from the tutorial mesh"
            )
        return max(1, round(value / denominator))
    return value // denominator


def set_resolution(
    run_dir: Path, level: int, tutorial_level: int, thickness: str
) -> None:
    """Rescale the block divisions in the plane from the tutorial mesh.

    The tutorial dictionaries hold one member of a mesh family, so a level is
    reached by scaling them by a power of two. The grading is left untouched,
    which keeps the ratio between the largest and the smallest cell. The
    single cell through the thickness is not scaled.
    """
    shift = level - tutorial_level
    if shift == 0:
        return

    hex_pattern = re.compile(
        r"(^\s*hex\s*\([^)]*\)\s*\w+\s*)\(\s*(\w+)\s+(\d+)\s+(\d+)\s*\)",
        re.MULTILINE,
    )
    mid_pattern = re.compile(
        r"(m4_define\(nxMid,\s*m4_ifelse\(hOverW,\s*)([0-9.]+)(,\s*)(\d+)"
        r"(,\s*)(\d+)(\)\))"
    )

    def scale_hex(match: re.Match[str]) -> str:
        nx = match.group(2)
        if nx.isdigit():
            nx = str(scaled(int(nx), shift, True, "block division"))
        ny = scaled(int(match.group(3)), shift, True, "block division")
        return f"{match.group(1)}({nx} {ny} {match.group(4)})"

    def scale_mid(match: re.Match[str]) -> str:
        # The first count is used for the thickness named in the template,
        # the second for the others
        first_used = thickness == match.group(2)
        first = scaled(int(match.group(4)), shift, first_used, "nxMid")
        second = scaled(int(match.group(6)), shift, not first_used, "nxMid")
        return (
            f"{match.group(1)}{match.group(2)}{match.group(3)}{first}"
            f"{match.group(5)}{second}{match.group(7)}"
        )

    for name, has_mid in (
        ("blockMeshDict.twoMaterials", False),
        ("blockMeshDict.threeMaterials.m4", True),
    ):
        block_mesh_dict = run_dir / "system" / name
        text, count = hex_pattern.subn(scale_hex, block_mesh_dict.read_text())
        if not count:
            raise RuntimeError(f"no block divisions found in {block_mesh_dict}")
        if has_mid:
            text, count = mid_pattern.subn(scale_mid, text)
            if count != 1:
                raise RuntimeError(
                    f"no mid-layer division found in {block_mesh_dict}"
                )
        block_mesh_dict.write_text(text)


def check_solver_log(solver_log: Path, label: str) -> None:
    """Reject a run that crashed, was killed, or never reached the end.

    A tutorial Allrun can return zero even when the solver aborts or is
    killed by the operating system, so the log has to be inspected directly.
    """
    text = solver_log.read_text(errors="replace")
    if re.search(
        r"FOAM FATAL|FOAM aborting|^ERROR$|\[stack trace\]", text, re.MULTILINE
    ):
        raise RuntimeError(
            f"{label} did not complete; see {solver_log}. A level that is "
            "too large for the available memory is the usual cause."
        )
    if not re.search(r"^End\s*$", text, re.MULTILINE):
        raise RuntimeError(
            f"{label} did not reach the end of the run; see {solver_log}"
        )


def read_profile(run_dir: Path, set_name: str, scale: float, stress: float) -> Profile:
    """Return a sampled set as (coordinate/scale, sigma_xx/stress).

    OpenFOAM writes the sets to postProcessing/sample and foam-extend to
    postProcessing/sets.
    """
    candidates = sorted(
        (run_dir / "postProcessing").glob(f"*/*/{set_name}_sigma.xy")
    )
    if not candidates:
        raise RuntimeError(f"no sampled {set_name} data found in {run_dir}")
    profile = []
    for line in candidates[-1].read_text().splitlines():
        fields = line.split()
        if len(fields) >= 7:
            profile.append((float(fields[0]) / scale, float(fields[1]) / stress))
    if len(profile) < 2:
        raise RuntimeError(f"fewer than two numeric rows in {candidates[-1]}")
    return profile


def read_reference_curve(path: Path) -> Profile:
    """Return a digitised curve as (r/W, sigma_xx/sigma_0), sorted by r/W.

    The columns can be separated by commas or by white space, and lines that
    start with # are comments. An empty or missing file gives an empty curve.
    """
    if not path.is_file():
        return []
    points = []
    for line in path.read_text().splitlines():
        if line.strip().startswith("#") or not line.strip():
            continue
        fields = line.replace(",", " ").split()
        points.append((float(fields[0]), float(fields[1])))
    return sorted(points)


def write_profile(path: Path, header: str, profile: Profile) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        f"# {header}\n"
        + "".join(f"{x:.12g} {value:.10g}\n" for x, value in profile)
    )


def interpolate_log(profile: Profile, r: float) -> float:
    """Interpolate a stress profile linearly in log(r)."""
    positive = [(x, value) for x, value in profile if x > 0.0]
    if r <= positive[0][0]:
        return positive[0][1]
    if r >= positive[-1][0]:
        return positive[-1][1]
    for (x0, v0), (x1, v1) in zip(positive, positive[1:]):
        if r <= x1:
            weight = (math.log(r) - math.log(x0)) / (math.log(x1) - math.log(x0))
            return v0 + weight * (v1 - v0)
    return positive[-1][1]


def mean_stress(profile: Profile) -> float:
    """Mean of sigma_xx/sigma_0 over the interface, 0 <= r/W <= 1.

    Equilibrium requires it to be one. The sample utility of foam-extend does
    not write the last point of the line, on the symmetry plane, so the stress
    is taken as uniform from the last sampled point to the symmetry plane,
    where its gradient is zero.
    """
    integral = sum(
        0.5 * (v0 + v1) * (x1 - x0)
        for (x0, v0), (x1, v1) in zip(profile, profile[1:])
    )
    return integral + profile[-1][1] * (1.0 - profile[-1][0])


def fitted_order(profile: Profile, window: list[float]) -> float:
    """Least-squares slope of -log(sigma) against log(r) in the window."""
    points = [
        (math.log(x), math.log(value))
        for x, value in profile
        if window[0] < x < window[1] and value > 0.0
    ]
    if len(points) < 2:
        return math.nan
    n = len(points)
    sx = sum(x for x, _ in points)
    sy = sum(y for _, y in points)
    sxx = sum(x * x for x, _ in points)
    sxy = sum(x * y for x, y in points)
    return -(n * sxy - sx * sy) / (n * sxx - sx * sx)


def analytical_order(material_a: dict, material_b: dict) -> float:
    """Order p of the singularity, sigma ~ r^(-p), in plane strain.

    It is the root lambda = 1 - p of the characteristic equation of Bogy
    (1970) for two bonded quarter planes, written with the Dundurs'
    parameters, found by bisection between lambda = 0.5 and the trivial root
    lambda = 1.
    """
    mu_a = material_a["E_pa"] / (2.0 * (1.0 + material_a["nu"]))
    mu_b = material_b["E_pa"] / (2.0 * (1.0 + material_b["nu"]))
    k = mu_b / mu_a
    m_a = 4.0 * (1.0 - material_a["nu"])
    m_b = 4.0 * (1.0 - material_b["nu"])
    alpha = (m_b - k * m_a) / (m_b + k * m_a)
    beta = ((m_b - 2.0) - k * (m_a - 2.0)) / (m_b + k * m_a)

    def f(lam: float) -> float:
        s = math.sin(0.5 * math.pi * lam) ** 2
        return (
            lam**2 * (lam**2 - 1.0) * alpha**2
            + 2.0 * lam**2 * (s - lam**2) * alpha * beta
            + (s - lam**2) ** 2 * beta**2
            + s * (1.0 - s)
        )

    low, high = 0.5, 0.9999
    for _ in range(60):
        middle = 0.5 * (low + high)
        if f(low) * f(middle) <= 0.0:
            high = middle
        else:
            low = middle
    return 1.0 - 0.5 * (low + high)


def profile_change(coarse: Profile, fine: Profile, window: list[float]) -> float:
    """RMS relative change between two meshes at the points of the coarser."""
    changes = [
        (interpolate_log(fine, x) - value) / value
        for x, value in coarse
        if window[0] <= x <= window[1]
    ]
    if not changes:
        return math.nan
    return math.sqrt(sum(value**2 for value in changes) / len(changes))


def interfaces(thickness: str) -> list[str]:
    return ["12"] if float(thickness) == 0.0 else ["13", "23"]


def run_case(
    variant_name: str,
    variant: dict,
    thickness: str,
    level: int,
    reference: dict,
    reuse: bool,
    cache: dict,
) -> dict:
    """Run one thickness on one mesh level, or return it if already run."""
    key = (variant_name, thickness, level)
    if key in cache:
        return cache[key]

    label = f"{variant_name} h/W = {thickness} mesh{level}"
    run_dir = WORK_DIR / variant_name / f"h{thickness}_mesh{level}"
    solver_log = run_dir / "log.solids4Foam"
    completed_case = (
        reuse
        and solver_log.exists()
        and re.search(r"^End\s*$", solver_log.read_text(), re.MULTILINE)
    )

    if completed_case:
        print(f"Reusing {label} in {run_dir}")
    else:
        if run_dir.exists():
            shutil.rmtree(run_dir)
        shutil.copytree(CASE_DIR, run_dir, ignore=ignored, symlinks=True)
        set_resolution(run_dir, level, int(reference["tutorial_level"]), thickness)
        print(f"Running {label} in {run_dir}", flush=True)
        with (run_dir / "log.Allverify").open("w") as log:
            completed = subprocess.run(
                ["./Allrun", thickness, variant["approach"]],
                cwd=run_dir,
                stdout=log,
                stderr=subprocess.STDOUT,
                check=False,
            )
        if completed.returncode:
            raise RuntimeError(f"{label} failed; see {run_dir / 'log.Allverify'}")
        if "PETSc is not installed" in (run_dir / "log.Allverify").read_text(
            errors="replace"
        ):
            raise RuntimeError(
                f"{label} needs PETSc, which is not available; select "
                "another variant with --variants segregated"
            )
        if not solver_log.exists():
            raise RuntimeError(f"solver log was not created in {run_dir}")

    check_solver_log(solver_log, label)

    half_width = float(reference["half_width_m"])
    stress = float(reference["applied_stress_pa"])
    cell_match = re.search(
        r"nCells:\s*(\d+)",
        (run_dir / "log.blockMesh").read_text(errors="replace"),
    )
    clock_matches = re.findall(
        r"ClockTime\s*=\s*([0-9.eE+-]+)", solver_log.read_text(errors="replace")
    )
    outcome = {
        "cells": int(cell_match.group(1)) if cell_match else 0,
        "clock_time_s": float(clock_matches[-1]) if clock_matches else math.nan,
        "interfaces": {
            name: read_profile(run_dir, f"interface{name}", half_width, stress)
            for name in interfaces(thickness)
        },
        "midLine": read_profile(run_dir, "midLine", half_width, stress),
        "midLineZoom": read_profile(run_dir, "midLineZoom", half_width, stress),
    }
    cache[key] = outcome
    return outcome


def interface_metrics(
    profile: Profile, interface: str, reference: dict
) -> dict[str, float]:
    materials = reference["materials"]
    order = fitted_order(profile, reference["order_fit_window_r_over_w"])
    analytical = analytical_order(materials[interface[0]], materials[interface[1]])
    return {
        "mean_stress_error": abs(mean_stress(profile) - 1.0),
        "fitted_order": order,
        "analytical_order": analytical,
        "order_relative_error": abs(order - analytical) / analytical,
    }


def reference_curve(reference: dict, interface: str, thickness: str) -> Profile:
    """Return the digitised curve of an interface, or an empty curve.

    Without a mid-layer, interfaces 13 and 23 are the same interface, so its
    curve can be stored under either name.
    """
    names = [interface]
    if float(thickness) == 0.0:
        names.extend(name for name in ("13", "23") if name != interface)
    for name in names:
        curve = read_reference_curve(
            REFERENCE_DIR
            / reference["thickness_study"]["reference_curve"].format(
                interface=name, thickness=thickness
            )
        )
        if curve:
            return curve
    return []


def write_reference_curves(reference: dict) -> list[str]:
    """Write the digitised curves found as white-space separated columns.

    The plot scripts read these copies. Returns the names of the curves.
    """
    found = []
    for thickness in reference["thickness_study"]["thicknesses"]:
        for interface in ("13", "23"):
            curve = reference_curve(reference, interface, thickness)
            if curve:
                name = f"reference_h{thickness}_interface{interface}.txt"
                write_profile(
                    PROFILE_DIR / name,
                    "r/W sigma_xx/sigma_0, digitised from Ioka et al. (2007)",
                    curve,
                )
                found.append(name)
    return found


def thickness_study(
    variant_name: str,
    variant: dict,
    thicknesses: list[str],
    reference: dict,
    args: argparse.Namespace,
    cache: dict,
) -> tuple[list[dict], int]:
    rows = []
    failures = 0
    level = int(reference["tutorial_level"])
    for thickness in thicknesses:
        try:
            outcome = run_case(
                variant_name, variant, thickness, level, reference, args.reuse, cache
            )
        except RuntimeError as error:
            print(f"ERROR: {error}", file=sys.stderr)
            failures += 1
            if not args.keep_going:
                raise
            continue
        for interface, profile in outcome["interfaces"].items():
            write_profile(
                PROFILE_DIR
                / f"thickness_{variant_name}_h{thickness}_interface{interface}.txt",
                "r/W sigma_xx/sigma_0",
                profile,
            )
            row = {
                "variant": variant_name,
                "thickness": thickness,
                "interface": interface,
                "cells": outcome["cells"],
                "clock_time_s": outcome["clock_time_s"],
            }
            row.update(interface_metrics(profile, interface, reference))
            rows.append(row)
    return rows, failures


def mesh_study(
    variant_name: str,
    variant: dict,
    levels: list[int],
    reference: dict,
    args: argparse.Namespace,
    cache: dict,
) -> tuple[list[dict], int]:
    rows = []
    failures = 0
    thickness = str(reference["mesh_study"]["thickness"])
    window = reference["mesh_study"]["change_window_r_over_w"]
    previous: dict[str, Profile] = {}
    for level in levels:
        try:
            outcome = run_case(
                variant_name, variant, thickness, level, reference, args.reuse, cache
            )
        except RuntimeError as error:
            print(f"ERROR: {error}", file=sys.stderr)
            failures += 1
            previous = {}
            if not args.keep_going:
                raise
            continue
        for name in ("midLine", "midLineZoom"):
            write_profile(
                PROFILE_DIR / f"mesh_{variant_name}_mesh{level}_{name}.txt",
                "x/W sigma_xx/sigma_0",
                outcome[name],
            )
        for interface, profile in outcome["interfaces"].items():
            write_profile(
                PROFILE_DIR
                / f"mesh_{variant_name}_mesh{level}_interface{interface}.txt",
                "r/W sigma_xx/sigma_0",
                profile,
            )
            row = {
                "variant": variant_name,
                "thickness": thickness,
                "level": level,
                "interface": interface,
                "cells": outcome["cells"],
                "clock_time_s": outcome["clock_time_s"],
            }
            row.update(interface_metrics(profile, interface, reference))
            row["profile_rms_relative_change"] = (
                profile_change(previous[interface], profile, window)
                if interface in previous
                else math.nan
            )
            rows.append(row)
        previous = outcome["interfaces"]
    return rows, failures


def write_csv(path: Path, rows: list[dict]) -> None:
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)


def summarise(
    thickness_rows: list[dict],
    mesh_rows: list[dict],
    variant_names: list[str],
    reference: dict,
    curves: list[str],
    quick: bool,
) -> bool:
    acceptance = reference["acceptance"]
    mean_tolerance = float(acceptance["mean_stress_error_tolerance"])
    order_tolerance = float(acceptance["bimaterial_order_relative_error_tolerance"])
    passed = True
    lines = [
        "# multiMaterialPlate verification summary",
        "",
        f"- Mode: {'quick smoke test' if quick else 'full verification'}",
        f"- Variants: {', '.join(variant_names)}",
    ]

    if thickness_rows:
        write_csv(POST_DIR / "thickness_study.csv", thickness_rows)
        lines.extend(
            [
                "",
                "## Thickness study",
                "",
                "| Variant | h/W | Interface | Cells | Mean stress error "
                "| Fitted order | Analytical order | Time (s) |",
                "|---|---:|---:|---:|---:|---:|---:|---:|",
            ]
        )
        for row in thickness_rows:
            lines.append(
                f"| {row['variant']} | {row['thickness']} | {row['interface']} "
                f"| {row['cells']} | {row['mean_stress_error']:.2e} "
                f"| {row['fitted_order']:.4f} | {row['analytical_order']:.4f} "
                f"| {row['clock_time_s']:.0f} |"
            )
            passed = passed and float(row["mean_stress_error"]) <= mean_tolerance
            if row["interface"] == "12":
                passed = (
                    passed and float(row["order_relative_error"]) <= order_tolerance
                )
        lines.extend(
            [
                "",
                "Digitised curves of Ioka et al. (2007) drawn in the plots, "
                "for a visual comparison only: "
                + (", ".join(curves) if curves else "none found"),
            ]
        )

    if mesh_rows:
        write_csv(POST_DIR / "mesh_convergence.csv", mesh_rows)
        lines.extend(
            [
                "",
                f"## Mesh study, h/W = {mesh_rows[0]['thickness']}",
                "",
                "| Variant | Level | Interface | Cells | Mean stress error "
                "| Fitted order | Analytical order | Profile change "
                "| Time (s) |",
                "|---|---:|---:|---:|---:|---:|---:|---:|---:|",
            ]
        )
        for row in mesh_rows:
            change = float(row["profile_rms_relative_change"])
            lines.append(
                f"| {row['variant']} | {row['level']} | {row['interface']} "
                f"| {row['cells']} | {row['mean_stress_error']:.2e} "
                f"| {row['fitted_order']:.4f} | {row['analytical_order']:.4f} "
                f"| {'-' if math.isnan(change) else f'{change:.2e}'} "
                f"| {row['clock_time_s']:.0f} |"
            )
        for name in variant_names:
            for interface in sorted({row["interface"] for row in mesh_rows}):
                series = [
                    row
                    for row in mesh_rows
                    if row["variant"] == name and row["interface"] == interface
                ]
                if not series:
                    continue
                passed = (
                    passed
                    and float(series[-1]["mean_stress_error"]) <= mean_tolerance
                )
                changes = [
                    float(row["profile_rms_relative_change"])
                    for row in series[1:]
                    if math.isfinite(float(row["profile_rms_relative_change"]))
                ]
                if not quick:
                    passed = (
                        passed
                        and len(changes) >= 2
                        and all(
                            right < left for left, right in zip(changes, changes[1:])
                        )
                    )

    lines.extend(["", f"- Result: {'PASS' if passed else 'FAIL'}"])
    summary = "\n".join(lines) + "\n"
    (POST_DIR / "verification_summary.md").write_text(summary)
    print(summary)
    return passed


def create_plots(
    variant_name: str,
    thicknesses: list[str],
    levels: list[int],
    reference: dict,
) -> None:
    if shutil.which("gnuplot") is None:
        return
    common = f"profiles='{PROFILE_DIR}'; variant='{variant_name}'; "
    jobs = []
    if thicknesses:
        for interface in ("13", "23"):
            jobs.append(
                (
                    "plotThicknessStudy.gnuplot",
                    f"thicknesses='{' '.join(thicknesses)}'; "
                    f"interface='{interface}'",
                    f"multiMaterialPlate_thickness_interface{interface}.pdf",
                )
            )
    if levels:
        thickness = str(reference["mesh_study"]["thickness"])
        for name in interfaces(thickness):
            jobs.append(
                (
                    "plotMeshStudy.gnuplot",
                    f"levels='{' '.join(str(level) for level in levels)}'; "
                    f"set='interface{name}'; thickness='{thickness}'",
                    f"multiMaterialPlate_mesh_interface{name}.pdf",
                )
            )
        jobs.append(
            (
                "plotMeshStudy.gnuplot",
                f"levels='{' '.join(str(level) for level in levels)}'; "
                f"set='midLineZoom'; thickness='{thickness}'",
                "multiMaterialPlate_mesh_midLine.pdf",
            )
        )
    for script, variables, output_name in jobs:
        output = POST_DIR / output_name
        completed = subprocess.run(
            [
                "gnuplot",
                "-e",
                f"{common}{variables}; output='{output}'",
                str(SCRIPT_DIR / script),
            ],
            cwd=VERIFY_DIR,
            text=True,
            check=False,
        )
        if completed.returncode:
            print(f"WARNING: gnuplot failed for {output_name}.")
        else:
            print(f"Plot: {output}")


def selection(value: str | None, default: list, available: list, what: str) -> list:
    chosen = [item.strip() for item in value.split(",")] if value else list(default)
    invalid = [item for item in chosen if item not in available]
    if not chosen or invalid or len(set(chosen)) != len(chosen):
        raise SystemExit(
            f"{what} must be a unique selection from: "
            f"{', '.join(str(item) for item in available)}"
        )
    return chosen


def main() -> int:
    args = parse_args()
    reference = json.loads(REFERENCE_FILE.read_text())

    studies = selection(args.studies, STUDIES, list(STUDIES), "studies")
    variant_names = selection(
        args.variants,
        reference["default_variants"],
        list(reference["variants"]),
        "variants",
    )
    all_thicknesses = reference["thickness_study"]["thicknesses"]
    thicknesses = selection(
        args.thicknesses, all_thicknesses, all_thicknesses, "thicknesses"
    )
    if args.levels:
        levels = [int(value) for value in args.levels.split(",")]
    else:
        levels = [int(value) for value in reference["mesh_study"]["default_levels"]]
    if args.quick:
        if not args.thicknesses:
            thicknesses = [all_thicknesses[0], all_thicknesses[-1]]
        levels = levels[:2]
    if "mesh" in studies:
        if len(levels) < 2 or any(level <= 0 for level in levels):
            raise SystemExit("at least two positive levels are required")
        if any(right <= left for left, right in zip(levels, levels[1:])):
            raise SystemExit("levels must be strictly increasing")

    missing = [
        command
        for command in ("blockMesh", "solids4Foam", "m4")
        if shutil.which(command) is None
    ]
    if missing:
        raise SystemExit(f"required command(s) not found: {', '.join(missing)}")

    WORK_DIR.mkdir(parents=True, exist_ok=True)
    POST_DIR.mkdir(parents=True, exist_ok=True)
    thickness_rows: list[dict] = []
    mesh_rows: list[dict] = []
    failures = 0
    cache: dict = {}
    try:
        for name in variant_names:
            variant = reference["variants"][name]
            if "thickness" in studies:
                rows, failed = thickness_study(
                    name, variant, thicknesses, reference, args, cache
                )
                thickness_rows.extend(rows)
                failures += failed
            if "mesh" in studies:
                rows, failed = mesh_study(
                    name, variant, levels, reference, args, cache
                )
                mesh_rows.extend(rows)
                failures += failed
    except RuntimeError:
        return 1

    if not thickness_rows and not mesh_rows:
        print("ERROR: no case completed", file=sys.stderr)
        return 1
    curves = write_reference_curves(reference)
    passed = summarise(
        thickness_rows, mesh_rows, variant_names, reference, curves, args.quick
    )
    create_plots(
        variant_names[0],
        thicknesses if "thickness" in studies else [],
        levels if "mesh" in studies else [],
        reference,
    )
    return 0 if passed and failures == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
