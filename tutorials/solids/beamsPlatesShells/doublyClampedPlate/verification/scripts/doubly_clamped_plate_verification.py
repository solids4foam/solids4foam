#!/usr/bin/env python3
"""Run the doublyClampedPlate mesh-convergence verification study.

The plate is clamped on two opposite edges and free on the other two, and
there is no closed-form Kirchhoff solution for it. The reference is therefore
a converged Kirchhoff plate finite element solution (12-DOF ACM rectangular
element) of the same problem, recorded in the reference JSON file. This study
refines the mesh, reads the deflection at the plate centre and at the midpoint
of a free edge from the written fields, and checks the observed order of
accuracy against that reference.
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
REFERENCE_FILE = REFERENCE_DIR / f"{CASE_DIR.name}_verification_references.json"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=f"Run the {CASE_DIR.name} verification study"
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
    parser.add_argument(
        "--coupled",
        action="store_true",
        help="use the block-coupled solver; results go under work/coupled",
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
        ignored_names.add("case.foam")
    if directory_path.name == "constant":
        ignored_names.add("polyMesh")
    if directory_path.name == "faMesh":
        ignored_names.update({"faceLabels", "faBoundary"})

    return ignored_names.intersection(names)


# ---------------------------------------------------------------------------
# Case preparation
# ---------------------------------------------------------------------------

def set_resolution(run_dir: Path, level: int, tutorial_level: int) -> int:
    """Scale the single square block of the tutorial mesh by a power of two.

    Only the two in-plane divisions are scaled; the single through-thickness
    cell of the plate is left alone.
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
    if divisions[0] != divisions[1]:
        raise RuntimeError(f"the block in {block_mesh_dict} is not square")
    block_mesh_dict.write_text(text)
    return divisions[0]


def read_divisions(run_dir: Path) -> int:
    """Read the block divisions already written into a case."""
    block_mesh_dict = run_dir / "system" / "blockMeshDict"
    match = re.search(
        r"^\s*hex\s*\([^)]*\)\s*\(\s*(\d+)\s+(\d+)\s+\d+\s*\)",
        block_mesh_dict.read_text(),
        re.MULTILINE,
    )
    if not match:
        raise RuntimeError(f"no block divisions found in {block_mesh_dict}")
    return int(match.group(1))


def set_clamped_relaxation(run_dir: Path, patches: list[str], value: float) -> None:
    """Set the clampedMoment relaxation factor of the clamped patches.

    The loop gain of the explicitly coupled clamped moment condition grows
    with refinement, so finer meshes need a smaller factor to converge.
    """
    field = run_dir / "0" / "finite-area" / "M"
    text = field.read_text()
    for patch in patches:
        pattern = re.compile(
            r"(^\s*" + re.escape(patch) + r"\s*\{[^}]*?relaxationFactor\s+)"
            r"[-+0-9.eE]+(\s*;)",
            re.MULTILINE | re.DOTALL,
        )
        text, count = pattern.subn(lambda m: f"{m.group(1)}{value:g}{m.group(2)}", text)
        if count != 1:
            raise RuntimeError(
                f"could not set the relaxation factor of patch {patch} in {field}"
            )
    field.write_text(text)


def set_coupled(run_dir: Path) -> None:
    """Select the block-coupled solver in constant/solidProperties."""
    path = run_dir / "constant" / "solidProperties"
    new_text, count = re.subn(
        r"^([ \t]*coupled\s+)\w+;",
        r"\1yes;",
        path.read_text(),
        count=1,
        flags=re.MULTILINE,
    )
    if count != 1:
        raise RuntimeError(f"could not find the coupled entry in {path}")
    path.write_text(new_text)


def check_solver_log(solver_log: Path, label: str) -> tuple[int, float]:
    """Reject a run that crashed or did not converge; return correctors and time."""
    text = solver_log.read_text(errors="replace")
    if re.search(
        r"FOAM FATAL|FOAM aborting|^ERROR$|\[stack trace\]", text, re.MULTILINE
    ):
        raise RuntimeError(f"{label} did not complete; see {solver_log}")
    if not re.search(r"^End\s*$", text, re.MULTILINE):
        raise RuntimeError(
            f"{label} did not reach the end of the run; see {solver_log}"
        )
    coupled = re.findall(r"^\s*Coupled iteration (\d+),", text, re.MULTILINE)
    if not coupled and "The residuals have converged" not in text:
        raise RuntimeError(
            f"{label} ran out of correctors before converging; see {solver_log}"
        )
    correctors = re.findall(r"^\s*(\d+), [-+0-9.eE]+, ", text, re.MULTILINE)
    if coupled:
        correctors = coupled
    clock = re.findall(r"ClockTime\s*=\s*([0-9.eE+-]+)", text)
    return (
        int(correctors[-1]) if correctors else -1,
        float(clock[-1]) if clock else math.nan,
    )


# ---------------------------------------------------------------------------
# Reading the OpenFOAM ascii files
# ---------------------------------------------------------------------------

def strip_header(path: Path) -> str:
    """Return the file contents after the FoamFile header, without comments."""
    text = path.read_text()
    text = text[text.index("}") + 1:]
    text = re.sub(r"/\*.*?\*/", "", text, flags=re.DOTALL)
    return re.sub(r"//.*", "", text)


def read_list_body(path: Path) -> str:
    """Return the text inside the outer parentheses of a list file."""
    text = strip_header(path)
    match = re.search(r"(\d+)\s*\((.*)\)", text, re.DOTALL)
    if not match:
        raise RuntimeError(f"no list found in {path}")
    return match.group(2)


def read_points(mesh_dir: Path) -> list[tuple[float, float, float]]:
    body = read_list_body(mesh_dir / "points")
    return [
        tuple(float(v) for v in m.split())  # type: ignore[misc]
        for m in re.findall(r"\(([^()]*)\)", body)
    ]


def read_faces(mesh_dir: Path) -> list[list[int]]:
    body = read_list_body(mesh_dir / "faces")
    return [
        [int(v) for v in m.split()] for m in re.findall(r"\d+\s*\(([^()]*)\)", body)
    ]


def read_labels(path: Path) -> list[int]:
    return [int(v) for v in read_list_body(path).split()]


def read_patch_range(mesh_dir: Path, patch: str) -> tuple[int, int]:
    text = strip_header(mesh_dir / "boundary")
    match = re.search(
        r"\b" + re.escape(patch) + r"\s*\{[^}]*?nFaces\s+(\d+)\s*;\s*startFace\s+(\d+)\s*;",
        text,
        re.DOTALL,
    )
    if not match:
        raise RuntimeError(f"patch {patch} not found in {mesh_dir / 'boundary'}")
    return int(match.group(2)), int(match.group(1))


def read_scalar_values(text: str, key: str, size: int | None = None) -> list[float]:
    match = re.search(
        key + r"\s+nonuniform\s+List<scalar>\s*(\d+)\s*\((.*?)\)\s*;", text, re.DOTALL
    )
    if match:
        return [float(v) for v in match.group(2).split()]
    match = re.search(key + r"\s+uniform\s+([-+0-9.eE]+)\s*;", text)
    if match and size is not None:
        return [float(match.group(1))] * size
    raise RuntimeError(f"no {key} entry found")


def read_area_field(path: Path, patch: str, patch_size: int) -> tuple[list[float], list[float]]:
    """Return the internal values and the boundary values on one patch."""
    text = strip_header(path)
    internal = read_scalar_values(text, "internalField")
    match = re.search(r"\n\s*" + re.escape(patch) + r"\s*\{(.*?)\n\s*\}", text, re.DOTALL)
    if not match:
        raise RuntimeError(f"patch {patch} not found in {path}")
    return internal, read_scalar_values(match.group(1), "value", patch_size)


def face_centre(points, face) -> tuple[float, float]:
    return (
        sum(points[i][0] for i in face) / len(face),
        sum(points[i][1] for i in face) / len(face),
    )


class Mesh:
    """The in-plane geometry of the one-layer finite area mesh."""

    def __init__(self, run_dir: Path):
        poly = run_dir / "constant" / "polyMesh"
        self.points = read_points(poly)
        self.faces = read_faces(poly)
        self.owner = read_labels(poly / "owner")
        self.face_labels = read_labels(
            run_dir / "constant" / "finite-area" / "faMesh" / "faceLabels"
        )
        self.poly = poly
        # Area face centres and the area face of each cell
        self.centres = [face_centre(self.points, self.faces[f]) for f in self.face_labels]
        self.area_face_of_cell = {
            self.owner[f]: i for i, f in enumerate(self.face_labels)
        }

    def patch_edges(self, patch: str) -> tuple[list[tuple[float, float]], list[int]]:
        """Edge centres of a patch and the area face adjacent to each edge.

        The finite area patch edges are built from the faces of the polyMesh
        patch in order, so the patch face centres give the edge centres.
        """
        start, size = read_patch_range(self.poly, patch)
        centres = [face_centre(self.points, self.faces[start + i]) for i in range(size)]
        adjacent = [self.area_face_of_cell[self.owner[start + i]] for i in range(size)]
        return centres, adjacent


def evaluate_metric(mesh: Mesh, run_dir: Path, time: str, metric: dict) -> float:
    field = run_dir / time / "finite-area" / "w"
    location = metric["location"]

    if location == "centre":
        internal, _ = read_area_field(field, "", 0) if False else (
            read_scalar_values(strip_header(field), "internalField"), None
        )
        xs = [c[0] for c in mesh.centres]
        ys = [c[1] for c in mesh.centres]
        cx = 0.5 * (min(xs) + max(xs))
        cy = 0.5 * (min(ys) + max(ys))
        d = [math.hypot(c[0] - cx, c[1] - cy) for c in mesh.centres]
        dmin = min(d)
        nearest = [i for i, di in enumerate(d) if di <= dmin * (1 + 1e-6)]
        return sum(internal[i] for i in nearest) / len(nearest)

    patch = metric["patch"]
    centres, adjacent = mesh.patch_edges(patch)
    internal, boundary = read_area_field(field, patch, len(centres))
    if len(boundary) != len(centres):
        raise RuntimeError(f"patch {patch} has {len(centres)} edges but {len(boundary)} values")

    # Guard the edge ordering assumption: each boundary value must be close
    # to the value of its adjacent cell, within the change of the deflection
    # over half a cell, which is a small part of its largest value
    scale = max(abs(v) for v in internal)
    for value, cell in zip(boundary, adjacent):
        if abs(value - internal[cell]) > 0.25 * scale:
            raise RuntimeError(
                f"the edge values of patch {patch} do not match their cells; "
                "the patch edge ordering assumption failed"
            )

    if location == "mean":
        return sum(boundary) / len(boundary)
    if location == "midpoint":
        xs = [c[0] for c in centres]
        ys = [c[1] for c in centres]
        along = 0 if (max(xs) - min(xs)) > (max(ys) - min(ys)) else 1
        s = [c[along] for c in centres]
        mid = 0.5 * (min(s) + max(s))
        d = [abs(si - mid) for si in s]
        dmin = min(d)
        nearest = [i for i, di in enumerate(d) if di <= dmin + 1e-9 * (max(s) - min(s))]
        return sum(boundary[i] for i in nearest) / len(nearest)
    raise RuntimeError(f"unknown metric location {location}")


# ---------------------------------------------------------------------------
# Running one level
# ---------------------------------------------------------------------------

def run_level(
    variant_name: str,
    variant: dict,
    level: int,
    reference: dict,
    reuse: bool,
    coupled: bool,
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
        divisions = read_divisions(run_dir)
    else:
        if run_dir.exists():
            shutil.rmtree(run_dir)
        shutil.copytree(CASE_DIR, run_dir, ignore=ignored, symlinks=True)
        divisions = set_resolution(run_dir, level, int(reference["tutorial_level"]))
        relaxation = reference["clamped_moment_relaxation"].get(str(level))
        if relaxation is None:
            raise RuntimeError(f"no clampedMoment relaxation factor given for level {level}")
        set_clamped_relaxation(run_dir, reference["clamped_patches"], float(relaxation))
        if coupled:
            set_coupled(run_dir)
        command = ["./Allrun"]
        if variant["approach"]:
            command.append(variant["approach"])
        print(
            f"Running {variant_name} mesh{level} ({divisions} x {divisions}) in {run_dir}",
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
                f"{variant_name} mesh{level} failed; see {run_dir / 'log.Allverify'}"
            )
        if not solver_log.exists():
            raise RuntimeError(f"solver log was not created in {run_dir}")

    correctors, clock_time = check_solver_log(solver_log, f"{variant_name} mesh{level}")

    time_dirs = sorted(
        (d for d in run_dir.iterdir() if re.fullmatch(r"[0-9.]+", d.name) and d.name != "0"),
        key=lambda d: float(d.name),
    )
    if not time_dirs:
        raise RuntimeError(f"no result time directory in {run_dir}")
    mesh = Mesh(run_dir)
    result: dict = {
        "variant": variant_name,
        "level": level,
        "divisions": divisions,
        "cells": len(mesh.face_labels),
        "spacing_m": float(reference["plate_length_m"]) / divisions,
        "correctors": correctors,
        "clock_time_s": clock_time,
    }
    for name, metric in reference["metrics"].items():
        value = evaluate_metric(mesh, run_dir, time_dirs[-1].name, metric)
        result[f"{name}_m"] = value
        result[f"{name}_error"] = (value - float(metric["reference_m"])) / float(
            metric["reference_m"]
        )
    if not all(
        math.isfinite(float(value))
        for key, value in result.items()
        if key != "variant"
    ):
        raise RuntimeError(f"non-finite result extracted from {run_dir}")
    return result


# ---------------------------------------------------------------------------
# Orders and reporting
# ---------------------------------------------------------------------------

def net_order(results: list[dict], metric: str) -> float:
    coarse = abs(float(results[0][metric]))
    fine = abs(float(results[-1][metric]))
    spacing_ratio = float(results[0]["spacing_m"]) / float(results[-1]["spacing_m"])
    if coarse <= 0 or fine <= 0 or spacing_ratio <= 1:
        return math.nan
    return math.log(coarse / fine) / math.log(spacing_ratio)


def successive_orders(results: list[dict], metric: str) -> list[float]:
    orders = []
    for left, right in zip(results, results[1:]):
        a = abs(float(left[metric]))
        b = abs(float(right[metric]))
        r = float(left["spacing_m"]) / float(right["spacing_m"])
        orders.append(math.log(a / b) / math.log(r) if a > 0 and b > 0 and r > 1 else math.nan)
    return orders


def richardson(results: list[dict], key: str) -> tuple[float, float]:
    """Observed order and extrapolated limit from the three finest levels.

    This uses the computed values alone, so it does not depend on the
    reference solution.
    """
    if len(results) < 3:
        return math.nan, math.nan
    w1, w2, w3 = (float(r[key]) for r in results[-3:])
    r = float(results[-3]["spacing_m"]) / float(results[-2]["spacing_m"])
    if (w1 - w2) * (w2 - w3) <= 0:
        return math.nan, math.nan
    order = math.log(abs(w1 - w2) / abs(w2 - w3)) / math.log(r)
    limit = w3 + (w3 - w2) / (r**order - 1)
    return order, limit


def write_results(
    results: list[dict],
    variant_names: list[str],
    reference: dict,
    quick: bool,
) -> bool:
    POST_DIR.mkdir(parents=True, exist_ok=True)
    acceptance = reference["acceptance"]
    metrics = list(reference["metrics"])
    primary = reference["primary_metric"]
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
            "# spacing_m " + " ".join(f"{m}_relative_error" for m in metrics) + "\n"
            + "".join(
                f"{float(row['spacing_m']):.10g} "
                + " ".join(f"{abs(float(row[f'{m}_error'])):.10g}" for m in metrics)
                + "\n"
                for row in grouped[name]
            )
        )

    passed = all(
        math.isfinite(float(row[f"{m}_error"])) for row in results for m in metrics
    )
    lines = [
        f"# {CASE_DIR.name} verification summary",
        "",
        f"- Mode: {'quick smoke test' if quick else 'full verification'}",
        f"- Variants: {', '.join(variant_names)}",
        f"- Reference: {reference['reference']['comment']}",
    ]
    for name in variant_names:
        rows = grouped[name]
        errors = [abs(float(row[f"{primary}_error"])) for row in rows]
        order = net_order(rows, f"{primary}_error")
        rich_order, rich_limit = richardson(rows, f"{primary}_m")
        if not quick:
            passed = (
                passed
                and all(right < left for left, right in zip(errors, errors[1:]))
                and errors[-1] <= float(acceptance["relative_tolerance"])
                and order >= float(acceptance["minimum_net_order"])
            )
        lines.extend(
            [
                "",
                f"## {name}",
                "",
                "| Level | Cells | Spacing (m) | Correctors | "
                + " | ".join(f"{m} (m) | {m} error" for m in metrics)
                + " |",
                "|---:|---:|---:|---:|" + "---:|---:|" * len(metrics),
            ]
        )
        for row in rows:
            lines.append(
                f"| {row['level']} | {row['cells']} | {float(row['spacing_m']):.4g} | "
                f"{row['correctors']} | "
                + " | ".join(
                    f"{float(row[f'{m}_m']):.6g} | {float(row[f'{m}_error']):+.3%}"
                    for m in metrics
                )
                + " |"
            )
        lines.extend(
            [
                "",
                f"- Successive orders ({primary}): "
                + ", ".join(f"{o:.2f}" for o in successive_orders(rows, f"{primary}_error")),
                f"- Net order ({primary}): {order:.3f}",
                f"- Richardson order from the three finest levels ({primary}): "
                f"{rich_order:.3f}, extrapolated value {rich_limit:.6g} m "
                f"({(rich_limit - float(reference['metrics'][primary]['reference_m'])) / float(reference['metrics'][primary]['reference_m']):+.3%} from the reference)"
                if math.isfinite(rich_order)
                else f"- Richardson order ({primary}): not available",
                "- Clock time per level (s): "
                + ", ".join(f"{float(row['clock_time_s']):.0f}" for row in rows),
            ]
        )
    lines.extend(["", f"- Result: {'PASS' if passed else 'FAIL'}"])
    summary = "\n".join(lines) + "\n"
    (POST_DIR / "verification_summary.md").write_text(summary)
    print(summary)
    print(f"CSV: {POST_DIR / 'mesh_convergence.csv'}")
    return passed


def create_plot(variant_names: list[str], reference: dict) -> None:
    plot_script = SCRIPT_DIR / "plotConvergence.gnuplot"
    if shutil.which("gnuplot") is None or not plot_script.is_file():
        return
    output = POST_DIR / f"{CASE_DIR.name}_convergence.pdf"
    metrics = list(reference["metrics"])
    completed = subprocess.run(
        [
            "gnuplot",
            "-e",
            f"profiles='{POST_DIR / 'profiles'}'; "
            f"variant='{variant_names[0]}'; "
            f"output='{output}'; "
            f"metric1='{metrics[0]}'; "
            f"metric2='{metrics[1] if len(metrics) > 1 else metrics[0]}'",
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
    global WORK_DIR, POST_DIR
    args = parse_args()
    if args.coupled:
        WORK_DIR = WORK_DIR / "coupled"
        POST_DIR = POST_DIR / "coupled"
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
        for command in ("blockMesh", "makeFaMesh", "solids4Foam")
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
                    run_level(
                        name, all_variants[name], level, reference, args.reuse, args.coupled
                    )
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
    create_plot(variant_names, reference)
    return 0 if passed and failures == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
