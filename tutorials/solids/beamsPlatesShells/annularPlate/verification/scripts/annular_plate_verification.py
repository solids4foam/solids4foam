#!/usr/bin/env python3
"""Run the annularPlate mesh-convergence verification study.

The annular Kirchhoff plate has closed-form solutions for the three
configurations of the tutorial: clamped inner edge with a free outer edge
under a uniform load or under the load p = p0 cos(phi), and clamped outer
edge with a free inner edge under a uniform load. This study refines the mesh,
measures the error of the computed deflection and moment sum against those
solutions, and checks the observed order of accuracy.
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
REFERENCE_FILE = REFERENCE_DIR / "annularPlate_verification_references.json"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run the annularPlate verification variants"
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
        ignored_names.add("case.foam")
    if directory_path.name == "constant":
        ignored_names.add("polyMesh")
    if directory_path.name == "faMesh":
        ignored_names.update({"faceLabels", "faBoundary"})

    return ignored_names.intersection(names)


# ---------------------------------------------------------------------------
# Closed-form solutions
# ---------------------------------------------------------------------------

def solve_linear(matrix: list[list[float]], rhs: list[float]) -> list[float]:
    """Gaussian elimination with partial pivoting for the small systems here."""
    n = len(rhs)
    m = [row[:] + [rhs[i]] for i, row in enumerate(matrix)]
    for c in range(n):
        pivot = max(range(c, n), key=lambda r: abs(m[r][c]))
        m[c], m[pivot] = m[pivot], m[c]
        for r in range(n):
            if r != c:
                f = m[r][c] / m[c][c]
                for k in range(c, n + 1):
                    m[r][k] -= f * m[c][k]
    return [m[i][n] / m[i][i] for i in range(n)]


class Axisymmetric:
    """Uniformly loaded annulus clamped at r = rc and free at r = rf.

    w = C1 + C2 r^2 + C3 ln r + C4 r^2 ln r + p r^4 / (64 D), with the
    constants from w = dw/dr = 0 at the clamped edge and vanishing radial
    moment and shear force at the free edge.
    """

    def __init__(self, rc: float, rf: float, p: float, D: float, nu: float):
        self.p, self.D, self.nu = p, D, nu
        # The shear force condition at the free edge fixes C4 directly
        self.c4 = -p * rf**2 / (8.0 * D)

        def residuals(c):
            self.c = (c[0], c[1], c[2], self.c4)
            return [self.w(rc), self.dw(rc), self.mr(rf)]

        base = residuals((0.0, 0.0, 0.0))
        columns = []
        for j in range(3):
            unit = [0.0] * 3
            unit[j] = 1.0
            col = residuals(unit)
            columns.append([col[i] - base[i] for i in range(3)])
        matrix = [[columns[j][i] for j in range(3)] for i in range(3)]
        c1, c2, c3 = solve_linear(matrix, [-b for b in base])
        self.c = (c1, c2, c3, self.c4)

    def w(self, r: float, phi: float = 0.0) -> float:
        c1, c2, c3, c4 = self.c
        return c1 + c2 * r * r + c3 * math.log(r) + c4 * r * r * math.log(r) + self.p * r**4 / (64.0 * self.D)

    def dw(self, r: float) -> float:
        c1, c2, c3, c4 = self.c
        return 2 * c2 * r + c3 / r + c4 * (2 * r * math.log(r) + r) + self.p * r**3 / (16.0 * self.D)

    def d2w(self, r: float) -> float:
        c1, c2, c3, c4 = self.c
        return 2 * c2 - c3 / r**2 + c4 * (2 * math.log(r) + 3) + 3 * self.p * r**2 / (16.0 * self.D)

    def mr(self, r: float) -> float:
        """Radial bending moment divided by -D."""
        return self.d2w(r) + self.nu * self.dw(r) / r

    def M(self, r: float, phi: float = 0.0) -> float:
        """Moment sum M = -D laplacian(w)."""
        return -self.D * (self.d2w(r) + self.dw(r) / r)


class CosineLoad:
    """Annulus clamped at r = rc, free at r = rf, loaded by p0 cos(phi).

    w = f(r) cos(phi) with f = A r + B r^3 + C / r + E r ln r + p0 r^4 / (45 D),
    with the constants from f = f' = 0 at the clamped edge and vanishing
    radial moment and Kirchhoff shear force at the free edge.
    """

    def __init__(self, rc: float, rf: float, p: float, D: float, nu: float):
        self.p, self.D, self.nu = p, D, nu

        def residuals(c):
            self.c = c
            return [self.f(rc), self.df(rc), self.mr(rf), self.vr(rf)]

        base = residuals((0.0,) * 4)
        columns = []
        for j in range(4):
            unit = [0.0] * 4
            unit[j] = 1.0
            col = residuals(unit)
            columns.append([col[i] - base[i] for i in range(4)])
        matrix = [[columns[j][i] for j in range(4)] for i in range(4)]
        self.c = tuple(solve_linear(matrix, [-b for b in base]))

    def f(self, r: float) -> float:
        a, b, c, e = self.c
        return a * r + b * r**3 + c / r + e * r * math.log(r) + self.p * r**4 / (45.0 * self.D)

    def df(self, r: float) -> float:
        a, b, c, e = self.c
        return a + 3 * b * r**2 - c / r**2 + e * (math.log(r) + 1) + 4 * self.p * r**3 / (45.0 * self.D)

    def d2f(self, r: float) -> float:
        a, b, c, e = self.c
        return 6 * b * r + 2 * c / r**3 + e / r + 12 * self.p * r**2 / (45.0 * self.D)

    def laplacian_f(self, r: float) -> float:
        """f'' + f'/r - f/r^2."""
        a, b, c, e = self.c
        return 8 * b * r + 2 * e / r + self.p * r**2 / (3.0 * self.D)

    def d_laplacian_f(self, r: float) -> float:
        a, b, c, e = self.c
        return 8 * b - 2 * e / r**2 + 2 * self.p * r / (3.0 * self.D)

    def d_f_over_r(self, r: float) -> float:
        """d/dr (f'/r - f/r^2) = d/dr of (1/r) d/dphi of the twist term."""
        a, b, c, e = self.c
        return 2 * b * r - 2 * c / r**3 + e / r + self.p * r**2 / (15.0 * self.D)

    def mr(self, r: float) -> float:
        return self.d2f(r) + self.nu * (self.df(r) / r - self.f(r) / r**2)

    def vr(self, r: float) -> float:
        return self.d_laplacian_f(r) - (1.0 - self.nu) * self.d_f_over_r(r) / r

    def w(self, r: float, phi: float) -> float:
        return self.f(r) * math.cos(phi)

    def M(self, r: float, phi: float) -> float:
        return -self.D * self.laplacian_f(r) * math.cos(phi)

    def dMdr(self, r: float, phi: float) -> float:
        return -self.D * self.d_laplacian_f(r) * math.cos(phi)


def exact_solution(variant: dict, plate: dict):
    a = float(plate["inner_radius_m"])
    b = float(plate["outer_radius_m"])
    h = float(plate["thickness_m"])
    p = float(plate["pressure_Pa"])
    E = float(plate["youngs_modulus_Pa"])
    nu = float(plate["poissons_ratio"])
    D = E * h**3 / (12.0 * (1.0 - nu**2))
    rc, rf = (a, b) if variant["clamped_patch"] == "inner" else (b, a)
    if variant["load"] == "cosine":
        return CosineLoad(rc, rf, p, D, nu), rf
    return Axisymmetric(rc, rf, p, D, nu), rf


# ---------------------------------------------------------------------------
# Case preparation
# ---------------------------------------------------------------------------

def set_resolution(run_dir: Path, level: int, tutorial_level: int) -> tuple[int, int]:
    """Scale the radial and circumferential divisions of the four blocks."""
    block_mesh_dict = run_dir / "system" / "blockMeshDict"
    pattern = re.compile(
        r"(^\s*hex\s*\([^)]*\)\s*)\(\s*(\d+)\s+(\d+)\s+(\d+)\s*\)", re.MULTILINE
    )
    shift = level - tutorial_level
    numerator = 2**shift if shift > 0 else 1
    denominator = 2 ** (-shift) if shift < 0 else 1
    divisions: list[tuple[int, int]] = []

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
        divisions.append((scaled[0], scaled[1]))
        return f"{match.group(1)}({scaled[0]} {scaled[1]} {match.group(4)})"

    text, count = pattern.subn(scale, block_mesh_dict.read_text())
    if count != 4:
        raise RuntimeError(
            f"expected four blocks in {block_mesh_dict}, found {count}"
        )
    if len(set(divisions)) != 1:
        raise RuntimeError(f"the blocks in {block_mesh_dict} differ in resolution")
    block_mesh_dict.write_text(text)
    return divisions[0]


def read_divisions(run_dir: Path) -> tuple[int, int]:
    block_mesh_dict = run_dir / "system" / "blockMeshDict"
    match = re.search(
        r"^\s*hex\s*\([^)]*\)\s*\(\s*(\d+)\s+(\d+)\s+\d+\s*\)",
        block_mesh_dict.read_text(),
        re.MULTILINE,
    )
    if not match:
        raise RuntimeError(f"no block divisions found in {block_mesh_dict}")
    return int(match.group(1)), int(match.group(2))


def set_clamped_relaxation(run_dir: Path, field_file: str, patch: str, value: float) -> None:
    """Set the clampedMoment relaxation factor of the clamped patch.

    The loop gain of the explicitly coupled clamped moment condition grows
    with refinement, so finer meshes need a smaller factor to converge.
    """
    field = run_dir / "0" / "finite-area" / field_file
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
    if "The residuals have converged" not in text:
        raise RuntimeError(
            f"{label} ran out of correctors before converging; see {solver_log}"
        )
    correctors = re.findall(r"^\s*(\d+), [-+0-9.eE]+, ", text, re.MULTILINE)
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


# ---------------------------------------------------------------------------
# Errors
# ---------------------------------------------------------------------------

def evaluate_errors(run_dir: Path, time: str, variant: dict, plate: dict) -> dict:
    exact, rf = exact_solution(variant, plate)
    mesh = Mesh(run_dir)
    free = variant["free_patch"]
    cosine = variant["load"] == "cosine"

    w_file = run_dir / time / "finite-area" / "w"
    m_file = run_dir / time / "finite-area" / "M"
    centres, adjacent = mesh.patch_edges(free)
    w, wb = read_area_field(w_file, free, len(centres))
    M, Mb = read_area_field(m_file, free, len(centres))
    if len(w) != len(mesh.centres):
        raise RuntimeError(f"{w_file} has {len(w)} values for {len(mesh.centres)} faces")

    # Guard the edge ordering assumption: each boundary value must be close
    # to the value of its adjacent cell, within the change of the deflection
    # over half a cell, which is a small part of its largest value
    scale = max(abs(v) for v in w)
    for value, cell in zip(wb, adjacent):
        if abs(value - w[cell]) > 0.25 * scale:
            raise RuntimeError(
                f"the edge values of patch {free} do not match their cells; "
                "the patch edge ordering assumption failed"
            )

    # Relative L2 errors of w and M over the cells
    ew2 = em2 = sw2 = sm2 = 0.0
    for (x, y), wi, mi in zip(mesh.centres, w, M):
        r, phi = math.hypot(x, y), math.atan2(y, x)
        wa, ma = exact.w(r, phi), exact.M(r, phi)
        ew2 += (wi - wa) ** 2
        em2 += (mi - ma) ** 2
        sw2 += wa**2
        sm2 += ma**2

    # Free-edge deflection: the mean over the edge for the axisymmetric
    # cases, the amplitude of the cos(phi) mode for the cosine load
    phis = [math.atan2(y, x) for x, y in centres]

    def project(values: list[float]) -> float:
        if cosine:
            return 2.0 * sum(v * math.cos(p) for v, p in zip(values, phis)) / len(values)
        return sum(values) / len(values)

    w_edge = project(wb)
    w_edge_exact = project([exact.w(rf, p) for p in phis])

    result = {
        "free_edge_w_m": w_edge,
        "free_edge_w_exact_m": w_edge_exact,
        "free_edge_w_error": (w_edge - w_edge_exact) / w_edge_exact,
        "w_l2_error": math.sqrt(ew2 / sw2),
        "M_l2_error": math.sqrt(em2 / sm2),
    }

    if cosine:
        # Normal gradient of M on the free edge, which the freeEdgeMoment
        # condition sets through the curvature term on this curved edge
        d = [math.hypot(c[0] - mesh.centres[cell][0], c[1] - mesh.centres[cell][1])
             for c, cell in zip(centres, adjacent)]
        sign = 1.0 if free == "outer" else -1.0
        g = project([(mb - M[cell]) / di for mb, cell, di in zip(Mb, adjacent, d)])
        g_exact = project([sign * exact.dMdr(rf, p) for p in phis])
        result["free_edge_dMdn"] = g
        result["free_edge_dMdn_exact"] = g_exact
        result["free_edge_dMdn_error"] = (g - g_exact) / g_exact
    return result


# ---------------------------------------------------------------------------
# Running one level
# ---------------------------------------------------------------------------

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
        print(f"Reusing {variant_name} mesh{level} in {run_dir}")
        radial, circumferential = read_divisions(run_dir)
    else:
        if run_dir.exists():
            shutil.rmtree(run_dir)
        shutil.copytree(CASE_DIR, run_dir, ignore=ignored, symlinks=True)
        radial, circumferential = set_resolution(
            run_dir, level, int(reference["tutorial_level"])
        )
        relaxation = reference["clamped_moment_relaxation"].get(str(level))
        if relaxation is None:
            raise RuntimeError(f"no clampedMoment relaxation factor given for level {level}")
        set_clamped_relaxation(
            run_dir,
            "M." + ("clampedInner" if variant["clamped_patch"] == "inner" else "clampedOuter"),
            variant["clamped_patch"],
            float(relaxation),
        )
        print(
            f"Running {variant_name} mesh{level} ({radial} x {4 * circumferential}) in {run_dir}",
            flush=True,
        )
        with (run_dir / "log.Allverify").open("w") as log:
            completed = subprocess.run(
                ["./Allrun", variant["configuration"]],
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
    plate = reference["plate"]
    result: dict = {
        "variant": variant_name,
        "level": level,
        "radial_cells": radial,
        "circumferential_cells": 4 * circumferential,
        "cells": 4 * radial * circumferential,
        "spacing_m": (float(plate["outer_radius_m"]) - float(plate["inner_radius_m"])) / radial,
        "correctors": correctors,
        "clock_time_s": clock_time,
    }
    result.update(evaluate_errors(run_dir, time_dirs[-1].name, variant, plate))
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
            "# spacing_m free_edge_w_relative_error w_relative_l2_error M_relative_l2_error\n"
            + "".join(
                f"{float(row['spacing_m']):.10g} "
                f"{abs(float(row['free_edge_w_error'])):.10g} "
                f"{float(row['w_l2_error']):.10g} "
                f"{float(row['M_l2_error']):.10g}\n"
                for row in grouped[name]
            )
        )

    passed = all(
        math.isfinite(float(row[key]))
        for row in results
        for key in ("free_edge_w_error", "w_l2_error", "M_l2_error")
    )
    lines = [
        "# annularPlate verification summary",
        "",
        f"- Mode: {'quick smoke test' if quick else 'full verification'}",
        f"- Variants: {', '.join(variant_names)}",
        "- Reference: closed-form Kirchhoff solutions of the annular plate",
    ]
    for name in variant_names:
        rows = grouped[name]
        errors = [abs(float(row["free_edge_w_error"])) for row in rows]
        order = net_order(rows, "free_edge_w_error")
        if not quick:
            passed = (
                passed
                and all(right < left for left, right in zip(errors, errors[1:]))
                and errors[-1] <= float(acceptance["free_edge_relative_tolerance"])
                and order >= float(acceptance["minimum_net_order"])
            )
        cosine = "free_edge_dMdn_error" in rows[0]
        header = (
            "| Level | Cells | Spacing (m) | Correctors | Free-edge w (m) | Exact (m) | "
            "Error | L2 error w | L2 error M |"
        )
        rule = "|---:|---:|---:|---:|---:|---:|---:|---:|---:|"
        if cosine:
            header += " dM/dn error |"
            rule += "---:|"
        lines.extend(["", f"## {name}", "", header, rule])
        for row in rows:
            line = (
                f"| {row['level']} | {row['cells']} | {float(row['spacing_m']):.4g} | "
                f"{row['correctors']} | {float(row['free_edge_w_m']):.6g} | "
                f"{float(row['free_edge_w_exact_m']):.6g} | "
                f"{float(row['free_edge_w_error']):+.3%} | "
                f"{float(row['w_l2_error']):.3%} | {float(row['M_l2_error']):.3%} |"
            )
            if cosine:
                line += f" {float(row['free_edge_dMdn_error']):+.3%} |"
            lines.append(line)
        lines.extend(
            [
                "",
                "- Successive orders (free-edge w): "
                + ", ".join(f"{o:.2f}" for o in successive_orders(rows, "free_edge_w_error")),
                f"- Net order (free-edge w): {order:.3f}",
                f"- Net order (L2 of w): {net_order(rows, 'w_l2_error'):.3f}",
                f"- Net order (L2 of M): {net_order(rows, 'M_l2_error'):.3f}",
            ]
        )
        if cosine:
            lines.append(
                f"- Net order (free-edge dM/dn): {net_order(rows, 'free_edge_dMdn_error'):.3f}"
            )
        lines.append(
            "- Clock time per level (s): "
            + ", ".join(f"{float(row['clock_time_s']):.0f}" for row in rows)
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
    output = POST_DIR / "annularPlate_convergence.pdf"
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
