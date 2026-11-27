#!/usr/bin/env python3
"""Run the manufacturedSolutionPlate mesh-convergence verification study.

The plate is loaded so that the exact deflection is the manufactured solution
defined in writeManufacturedSolution.py of the case, which satisfies the
clamped, free and symmetry edge conditions of the tutorial. This study refines
the mesh for each mesh variant (uniform quadrilaterals, randomly perturbed
quadrilaterals and unstructured triangles), measures the errors of the
deflection, the moment sum and the rotation against the exact solution, and
checks the observed order of accuracy where the discretisation is expected
to converge.
"""

from __future__ import annotations

import argparse
import csv
import importlib.util
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
    REFERENCE_DIR / "manufacturedSolutionPlate_verification_references.json"
)

# The outward normals of the plate edges, by patch name
NORMALS = {
    "left": (-1.0, 0.0),
    "right": (1.0, 0.0),
    "top": (0.0, 1.0),
    "bottom": (0.0, -1.0),
}
FREE_PATCH = "right"
CLAMPED_PATCH = "left"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run the manufacturedSolutionPlate verification variants"
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
    ignored_names.update(name for name in names if name.endswith(".msh"))

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
    if directory_path.name == "finite-area" and directory_path.parent.name == "0":
        ignored_names.update({"p", "wExact", "MExact"})

    return ignored_names.intersection(names)


# ---------------------------------------------------------------------------
# The manufactured solution, imported from the case
# ---------------------------------------------------------------------------

def load_manufactured():
    """Import writeManufacturedSolution.py of the case as a module."""
    path = CASE_DIR / "writeManufacturedSolution.py"
    spec = importlib.util.spec_from_file_location("manufactured_plate", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"could not import {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


# ---------------------------------------------------------------------------
# Case preparation
# ---------------------------------------------------------------------------

def cells_per_side(level: int, reference: dict) -> int:
    return int(reference["cells_per_side_level_1"]) * 2 ** (level - 1)


def set_resolution(run_dir: Path, level: int, reference: dict, length: float) -> int:
    """Set the blockMesh divisions and the Gmsh spacing for the level."""
    cells = cells_per_side(level, reference)

    block_mesh_dict = run_dir / "system" / "blockMeshDict"
    pattern = re.compile(
        r"(^\s*hex\s*\([^)]*\)\s*)\(\s*(\d+)\s+(\d+)\s+(\d+)\s*\)", re.MULTILINE
    )
    text, count = pattern.subn(
        lambda m: f"{m.group(1)}({cells} {cells} {m.group(4)})",
        block_mesh_dict.read_text(),
    )
    if count != 1:
        raise RuntimeError(f"expected one block in {block_mesh_dict}, found {count}")
    block_mesh_dict.write_text(text)

    spacing_file = run_dir / "gmsh" / "meshSpacing.geo"
    text, count = re.subn(
        r"^dx\s*=\s*[-+0-9.eE]+\s*;",
        f"dx = {length / cells:.10g};",
        spacing_file.read_text(),
        flags=re.MULTILINE,
    )
    if count != 1:
        raise RuntimeError(f"could not set dx in {spacing_file}")
    spacing_file.write_text(text)
    return cells


def set_clamped_relaxation(run_dir: Path, patch: str, value: float) -> None:
    """Set the clampedMoment relaxation factor of the clamped patch.

    The loop gain of the explicitly coupled clamped moment condition grows
    with refinement, so finer meshes need a smaller factor to converge.
    """
    field = run_dir / "0" / "finite-area" / "M"
    text = field.read_text()
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


def read_vector_values(text: str, key: str) -> list[tuple[float, float, float]]:
    match = re.search(
        key + r"\s+nonuniform\s+List<vector>\s*(\d+)\s*\((.*?)\)\s*;", text, re.DOTALL
    )
    if not match:
        raise RuntimeError(f"no nonuniform {key} entry found")
    return [
        tuple(float(v) for v in m.split())  # type: ignore[misc]
        for m in re.findall(r"\(([^()]*)\)", match.group(2))
    ]


def read_patch_values(text: str, patch: str, size: int) -> list[float] | None:
    """Boundary values of a patch, or None for a patch without a value entry."""
    match = re.search(r"\n\s*" + re.escape(patch) + r"\s*\{(.*?)\n\s*\}", text, re.DOTALL)
    if not match:
        raise RuntimeError(f"patch {patch} not found")
    if "value" not in match.group(1):
        return None
    return read_scalar_values(match.group(1), "value", size)


def face_centre(points, face) -> tuple[float, float]:
    """Area-weighted centre of the polygon, as OpenFOAM computes it."""
    n = len(face)
    px = sum(points[i][0] for i in face) / n
    py = sum(points[i][1] for i in face) / n
    if n == 3:
        return px, py
    sum_a = sum_ax = sum_ay = 0.0
    for i in range(n):
        x1, y1 = points[face[i]][0], points[face[i]][1]
        x2, y2 = points[face[(i + 1) % n]][0], points[face[(i + 1) % n]][1]
        area = abs((x1 - px) * (y2 - py) - (x2 - px) * (y1 - py))
        sum_a += area
        sum_ax += area * (x1 + x2 + px) / 3.0
        sum_ay += area * (y1 + y2 + py) / 3.0
    return sum_ax / sum_a, sum_ay / sum_a


def face_area(points, face) -> float:
    area = 0.0
    n = len(face)
    for i in range(n):
        x1, y1 = points[face[i]][0], points[face[i]][1]
        x2, y2 = points[face[(i + 1) % n]][0], points[face[(i + 1) % n]][1]
        area += x1 * y2 - x2 * y1
    return 0.5 * abs(area)


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
        self.centres = [face_centre(self.points, self.faces[f]) for f in self.face_labels]
        self.areas = [face_area(self.points, self.faces[f]) for f in self.face_labels]
        self.area_face_of_cell = {
            self.owner[f]: i for i, f in enumerate(self.face_labels)
        }
        # Height of the plate surface, to pick the edge end points of a patch face
        self.z_top = max(self.points[i][2] for i in self.faces[self.face_labels[0]])

    def patch_edges(self, patch: str) -> tuple[list[tuple[float, float]], list[int]]:
        """Edge centres of a patch and the area face adjacent to each edge.

        The finite area patch edges are built from the faces of the polyMesh
        patch in order, so the midpoints of the plate-surface edges of the
        patch faces give the edge centres.
        """
        start, size = read_patch_range(self.poly, patch)
        centres = []
        adjacent = []
        for i in range(size):
            face = self.faces[start + i]
            top = [p for p in face if abs(self.points[p][2] - self.z_top) < 1e-9]
            if len(top) != 2:
                raise RuntimeError(f"patch face {start + i} has no plate-surface edge")
            centres.append(
                (
                    0.5 * (self.points[top[0]][0] + self.points[top[1]][0]),
                    0.5 * (self.points[top[0]][1] + self.points[top[1]][1]),
                )
            )
            adjacent.append(self.area_face_of_cell[self.owner[start + i]])
        return centres, adjacent


# ---------------------------------------------------------------------------
# Errors
# ---------------------------------------------------------------------------

def relative_norms(values, exact, areas) -> tuple[float, float]:
    """Area-weighted relative L2 error and relative L-infinity error."""
    num = den = 0.0
    max_err = max_exact = 0.0
    for v, e, a in zip(values, exact, areas):
        num += a * (v - e) ** 2
        den += a * e * e
        max_err = max(max_err, abs(v - e))
        max_exact = max(max_exact, abs(e))
    return math.sqrt(num / den), max_err / max_exact


def evaluate_errors(run_dir: Path, time: str, manufactured) -> dict:
    exact = manufactured.Manufactured(
        manufactured.L, manufactured.W, manufactured.NU, manufactured.W0
    )
    D = manufactured.D
    mesh = Mesh(run_dir)

    w_text = strip_header(run_dir / time / "finite-area" / "w")
    m_text = strip_header(run_dir / time / "finite-area" / "M")
    theta_text = strip_header(run_dir / time / "finite-area" / "theta")
    w = read_scalar_values(w_text, "internalField")
    M = read_scalar_values(m_text, "internalField")
    theta = read_vector_values(theta_text, "internalField")
    if len(w) != len(mesh.centres):
        raise RuntimeError(f"w has {len(w)} values for {len(mesh.centres)} faces")

    w_exact = [exact.w(x, y) for x, y in mesh.centres]
    m_exact = [exact.M(x, y, D) for x, y in mesh.centres]
    grad_exact = [exact.grad_w(x, y) for x, y in mesh.centres]

    w_l2, w_linf = relative_norms(w, w_exact, mesh.areas)
    m_l2, m_linf = relative_norms(M, m_exact, mesh.areas)

    # Rotation theta = -grad(w): relative L2 error of the vector
    num = den = 0.0
    for (tx, ty, _), (gx, gy), a in zip(theta, grad_exact, mesh.areas):
        num += a * ((tx + gx) ** 2 + (ty + gy) ** 2)
        den += a * (gx * gx + gy * gy)
    theta_l2 = math.sqrt(num / den)

    result = {
        "w_l2_error": w_l2,
        "w_linf_error": w_linf,
        "M_l2_error": m_l2,
        "M_linf_error": m_linf,
        "theta_l2_error": theta_l2,
    }

    # Edge values: the finite area boundary values sit at the foot of the
    # normal from the adjacent cell centre, which is where OpenFOAM's finite
    # area discretisation places them
    w_max = max(abs(v) for v in w_exact)
    m_max = max(abs(v) for v in m_exact)
    for patch, field_text, exact_fn, scale, key in (
        (FREE_PATCH, w_text, lambda x, y: exact.w(x, y), w_max, "free_edge_w_linf_error"),
        (FREE_PATCH, m_text, lambda x, y: exact.M(x, y, D), m_max, "free_edge_M_linf_error"),
        (CLAMPED_PATCH, m_text, lambda x, y: exact.M(x, y, D), m_max, "clamped_edge_M_linf_error"),
    ):
        centres, adjacent = mesh.patch_edges(patch)
        values = read_patch_values(field_text, patch, len(centres))
        if values is None:
            raise RuntimeError(f"patch {patch} has no boundary values")
        nx, ny = NORMALS[patch]
        internal = w if field_text is w_text else M
        max_err = 0.0
        for value, (cx, cy), cell in zip(values, centres, adjacent):
            px, py = mesh.centres[cell]
            # Guard the edge ordering assumption: the edge value must be
            # close to the value of its adjacent cell
            if abs(value - internal[cell]) > 0.5 * scale:
                raise RuntimeError(
                    f"the edge values of patch {patch} do not match their cells; "
                    "the patch edge ordering assumption failed"
                )
            dn = (cx - px) * nx + (cy - py) * ny
            foot = (px + dn * nx, py + dn * ny)
            max_err = max(max_err, abs(value - exact_fn(*foot)))
        result[key] = max_err / scale

    return result


# ---------------------------------------------------------------------------
# Running one level
# ---------------------------------------------------------------------------

def run_level(
    variant_name: str,
    variant: dict,
    level: int,
    reference: dict,
    manufactured,
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
    cells_side = cells_per_side(level, reference)

    if completed_case:
        print(f"Reusing {variant_name} mesh{level} in {run_dir}")
    else:
        if run_dir.exists():
            shutil.rmtree(run_dir)
        shutil.copytree(CASE_DIR, run_dir, ignore=ignored, symlinks=True)
        set_resolution(run_dir, level, reference, manufactured.L)
        relaxation = reference["clamped_moment_relaxation"].get(str(level))
        if relaxation is None:
            raise RuntimeError(f"no clampedMoment relaxation factor given for level {level}")
        set_clamped_relaxation(run_dir, CLAMPED_PATCH, float(relaxation))
        if coupled:
            set_coupled(run_dir)
        print(
            f"Running {variant_name} mesh{level} ({cells_side} cells per side) in {run_dir}",
            flush=True,
        )
        with (run_dir / "log.Allverify").open("w") as log:
            completed = subprocess.run(
                ["./Allrun", variant["mesh"]],
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
    errors = evaluate_errors(run_dir, time_dirs[-1].name, manufactured)
    n_cells = len(Mesh(run_dir).centres)
    result: dict = {
        "variant": variant_name,
        "level": level,
        "cells": n_cells,
        "spacing_m": math.sqrt(manufactured.L * manufactured.W / n_cells),
        "correctors": correctors,
        "clock_time_s": clock_time,
    }
    result.update(errors)
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


def write_results(
    results: list[dict],
    variant_names: list[str],
    reference: dict,
    quick: bool,
) -> bool:
    POST_DIR.mkdir(parents=True, exist_ok=True)
    acceptance = reference["acceptance"]
    variants = reference["variants"]
    grouped = {
        name: [row for row in results if row["variant"] == name]
        for name in variant_names
    }
    fieldnames: list[str] = []
    for row in results:
        for key in row:
            if key not in fieldnames:
                fieldnames.append(key)

    with (POST_DIR / "mesh_convergence.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(results)

    (POST_DIR / "profiles").mkdir(parents=True, exist_ok=True)
    for name in variant_names:
        (POST_DIR / "profiles" / f"{name}_errors.txt").write_text(
            "# spacing_m w_relative_l2_error M_relative_l2_error w_relative_linf_error\n"
            + "".join(
                f"{float(row['spacing_m']):.10g} "
                f"{float(row['w_l2_error']):.10g} "
                f"{float(row['M_l2_error']):.10g} "
                f"{float(row['w_linf_error']):.10g}\n"
                for row in grouped[name]
            )
        )

    passed = True
    lines = [
        "# manufacturedSolutionPlate verification summary",
        "",
        f"- Mode: {'quick smoke test' if quick else 'full verification'}",
        f"- Variants: {', '.join(variant_names)}",
        "- Reference: the manufactured solution of writeManufacturedSolution.py",
        "- Errors: relative L2 (area weighted) and L-infinity norms over the cells;"
        " edge errors are the largest error of the boundary values relative to"
        " the largest exact value",
    ]
    for name in variant_names:
        rows = grouped[name]
        check = bool(variants[name].get("check_orders", False))
        order_w = net_order(rows, "w_l2_error")
        order_m = net_order(rows, "M_l2_error")
        if not quick and check:
            w_errors = [float(row["w_l2_error"]) for row in rows]
            passed = (
                passed
                and all(right < left for left, right in zip(w_errors, w_errors[1:]))
                and w_errors[-1] <= float(acceptance["finest_w_l2_tolerance"])
                and order_w >= float(acceptance["minimum_net_order_w_l2"])
                and order_m >= float(acceptance["minimum_net_order_M_l2"])
            )
        lines.extend(
            [
                "",
                f"## {name}: {variants[name]['description']}",
                "",
                "| Level | Cells | Spacing (m) | Corr. | L2 w | Linf w | L2 M | Linf M |"
                " L2 theta | Free edge w | Free edge M | Clamped edge M |",
                "|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
            ]
        )
        for row in rows:
            lines.append(
                f"| {row['level']} | {row['cells']} | {float(row['spacing_m']):.4g} | "
                f"{row['correctors']} | "
                f"{float(row['w_l2_error']):.3%} | {float(row['w_linf_error']):.3%} | "
                f"{float(row['M_l2_error']):.3%} | {float(row['M_linf_error']):.3%} | "
                f"{float(row['theta_l2_error']):.3%} | "
                f"{float(row['free_edge_w_linf_error']):.3%} | "
                f"{float(row['free_edge_M_linf_error']):.3%} | "
                f"{float(row['clamped_edge_M_linf_error']):.3%} |"
            )
        lines.extend(
            [
                "",
                "- Successive orders (L2 w): "
                + ", ".join(f"{o:.2f}" for o in successive_orders(rows, "w_l2_error")),
                f"- Net order (L2 w): {order_w:.3f}",
                f"- Net order (Linf w): {net_order(rows, 'w_linf_error'):.3f}",
                f"- Net order (L2 M): {order_m:.3f}",
                f"- Net order (L2 theta): {net_order(rows, 'theta_l2_error'):.3f}",
                "- Orders checked against the acceptance criteria: "
                + ("yes" if check else "no"),
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


def create_plot(variant_names: list[str]) -> None:
    plot_script = SCRIPT_DIR / "plotConvergence.gnuplot"
    if shutil.which("gnuplot") is None or not plot_script.is_file():
        return
    output = POST_DIR / "manufacturedSolutionPlate_convergence.pdf"
    completed = subprocess.run(
        [
            "gnuplot",
            "-e",
            f"profiles='{POST_DIR / 'profiles'}'; "
            f"variants='{' '.join(variant_names)}'; "
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
    global WORK_DIR, POST_DIR
    args = parse_args()
    if args.coupled:
        WORK_DIR = WORK_DIR / "coupled"
        POST_DIR = POST_DIR / "coupled"
    reference = json.loads(REFERENCE_FILE.read_text())
    manufactured = load_manufactured()
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
    if any(all_variants[name]["mesh"] == "tri" for name in variant_names):
        if shutil.which("gmsh") is None:
            raise SystemExit("gmsh is required for the tri variant")

    WORK_DIR.mkdir(parents=True, exist_ok=True)
    POST_DIR.mkdir(parents=True, exist_ok=True)
    results: list[dict] = []
    failures = 0
    for name in variant_names:
        for level in levels:
            try:
                results.append(
                    run_level(
                        name,
                        all_variants[name],
                        level,
                        reference,
                        manufactured,
                        args.reuse,
                        args.coupled,
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
    create_plot(variant_names)
    return 0 if passed and failures == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
