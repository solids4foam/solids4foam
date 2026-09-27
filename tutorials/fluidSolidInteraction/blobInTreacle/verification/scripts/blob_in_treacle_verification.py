#!/usr/bin/env python3
"""Run the blobInTreacle mesh and time-step verification studies.

The mesh study refines the tutorial fluid and solid meshes, runs each level
until the flow is steady, and compares the steady top-point displacement and
interface shape with the steady solution of Liu (arXiv:1401.0082), the
companion preprint of Liu et al. (2014). The time-step study repeats the
tutorial up to t = 1 s with successively halved time steps and checks the
observed temporal order of the top-point displacement, as Liu et al. (2014)
do for the same case.
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
REFERENCE_FILE = REFERENCE_DIR / "blobInTreacle_verification_references.json"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run the blobInTreacle verification studies"
    )
    parser.add_argument(
        "--study",
        choices=("all", "mesh", "time"),
        default="all",
        help="which study to run (default: all)",
    )
    parser.add_argument(
        "--levels",
        help="comma-separated mesh refinement factors relative to the "
        "tutorial mesh, e.g. 0.5,1,2",
    )
    parser.add_argument(
        "--quick",
        action="store_true",
        help="run a short smoke test of both studies",
    )
    parser.add_argument(
        "--reuse",
        action="store_true",
        help="reuse completed cases under verification/work",
    )
    parser.add_argument(
        "--cores",
        type=int,
        default=1,
        help="MPI ranks for the mesh levels finer than the tutorial mesh "
        "(default: 1, i.e. serial)",
    )
    parser.add_argument(
        "--keep-going",
        action="store_true",
        help="continue after an individual case fails",
    )
    args = parser.parse_args()
    if args.cores < 1:
        parser.error("--cores must be a positive integer")
    return args


# ---------------------------------------------------------------------------
# Case set-up
# ---------------------------------------------------------------------------


def ignored(directory: str, names: list[str]) -> set[str]:
    """Copy only the inputs a fresh run needs, never a previous result."""
    ignored_names = {"verification", "regressionTests", "postProcessing",
                     "images", "dynamicCode"}
    ignored_names.update(name for name in names if name.startswith("processor"))
    ignored_names.update(name for name in names if name.startswith("log."))
    ignored_names.update(name for name in names if name.endswith(".pdf"))
    directory_path = Path(directory)
    if directory_path == CASE_DIR:
        ignored_names.update(
            name
            for name in names
            if name != "0" and re.fullmatch(r"[0-9]+(?:\.[0-9]+)?(?:e-?[0-9]+)?", name)
        )
    if (
        directory_path.name in {"fluid", "solid"}
        and directory_path.parent.name == "constant"
    ):
        ignored_names.add("polyMesh")
    return ignored_names.intersection(names)


def refine_mesh(path: Path, factor: float) -> None:
    """Scale the in-plane block divisions, leaving the one-cell span alone."""
    pattern = re.compile(r"(hex\s*\([^)]*\)\s*)\(\s*(\d+)\s+(\d+)\s+(\d+)\s*\)")

    def scale(match: re.Match[str]) -> str:
        divisions = []
        for index in (2, 3):
            value = int(match.group(index)) * factor
            if abs(value - round(value)) > 1e-9 or round(value) < 1:
                raise RuntimeError(
                    f"refinement factor {factor:g} does not divide the block "
                    f"divisions in {path}"
                )
            divisions.append(int(round(value)))
        return f"{match.group(1)}({divisions[0]} {divisions[1]} {match.group(4)})"

    text, count = pattern.subn(scale, path.read_text())
    if not count:
        raise RuntimeError(f"no block divisions found in {path}")
    path.write_text(text)


def set_entry(path: Path, key: str, value: str) -> None:
    text, count = re.subn(
        rf"(^\s*{re.escape(key)}\s+)[^;]+;",
        rf"\g<1>{value};",
        path.read_text(),
        count=1,
        flags=re.MULTILINE,
    )
    if count != 1:
        raise RuntimeError(f"could not set {key} in {path}")
    path.write_text(text)


def set_cores(run_dir: Path, cores: int) -> None:
    for dictionary in (
        run_dir / "system" / "decomposeParDict",
        run_dir / "system" / "fluid" / "decomposeParDict",
        run_dir / "system" / "solid" / "decomposeParDict",
    ):
        set_entry(dictionary, "numberOfSubdomains", str(cores))


def check_solver_log(solver_log: Path, end_time: float, label: str) -> None:
    """Reject a run that crashed or stopped short of the end time.

    The tutorial Allrun returns zero even when the solver aborts, so a failed
    case would otherwise be read as a valid result.
    """
    text = solver_log.read_text(errors="replace")
    if re.search(r"FOAM FATAL|FOAM aborting|^ERROR$|\[stack trace\]", text,
                 re.MULTILINE):
        raise RuntimeError(f"{label} diverged or aborted; see {solver_log}")
    times = [float(match) for match in
             re.findall(r"^Time = ([0-9.eE+-]+)", text, re.MULTILINE)]
    if not times or times[-1] < end_time - 1.0e-8:
        reached = times[-1] if times else 0.0
        raise RuntimeError(
            f"{label} stopped at t = {reached:g} before the end time "
            f"{end_time:g}; see {solver_log}"
        )


def run_case(
    name: str,
    factor: float,
    delta_t: float,
    end_time: float,
    cores: int,
    reuse: bool,
) -> Path:
    run_dir = WORK_DIR / name
    solver_log = run_dir / "log.solids4Foam"
    if (
        reuse
        and solver_log.exists()
        and re.search(r"^End\s*$", solver_log.read_text(errors="replace"),
                      re.MULTILINE)
    ):
        print(f"Reusing {name} in {run_dir}")
    else:
        if run_dir.exists():
            shutil.rmtree(run_dir)
        shutil.copytree(CASE_DIR, run_dir, ignore=ignored, symlinks=True)
        if factor != 1:
            refine_mesh(run_dir / "system" / "fluid" / "blockMeshDict", factor)
            refine_mesh(run_dir / "system" / "solid" / "blockMeshDict", factor)
        control = run_dir / "system" / "controlDict"
        set_entry(control, "endTime", f"{end_time:.10g}")
        set_entry(control, "deltaT", f"{delta_t:.10g}")
        command = ["./Allrun"]
        if cores > 1:
            set_cores(run_dir, cores)
            command.append("parallel")
        print(
            f"Running {name} (mesh factor {factor:g}, deltaT = {delta_t:g} s, "
            f"{cores} core(s)) in {run_dir}",
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
            raise RuntimeError(f"{name} failed; see {run_dir / 'log.Allverify'}")
        if not solver_log.exists():
            raise RuntimeError(f"solver log was not created in {run_dir}")
    check_solver_log(solver_log, end_time, name)
    return run_dir


# ---------------------------------------------------------------------------
# Post-processing
# ---------------------------------------------------------------------------


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


def displacement_history(run_dir: Path) -> list[tuple[float, float, float]]:
    candidates = sorted(
        run_dir.glob("postProcessing/**/solidPointDisplacement_*.dat")
    )
    if not candidates:
        raise RuntimeError(f"no solidPointDisplacement data in {run_dir}")
    rows = numeric_rows(candidates[-1])
    if not rows:
        raise RuntimeError(f"no numeric rows in {candidates[-1]}")
    return [(row[0], row[1], row[2]) for row in rows if len(row) >= 3]


def value_at(history: list[tuple[float, float, float]], time: float):
    for row in history:
        if abs(row[0] - time) < 1e-8:
            return row
    raise RuntimeError(f"no displacement sample at t = {time:g}")


def steady_spread(history, window: float) -> float:
    """Relative spread of the x displacement over the closing window."""
    end_time = history[-1][0]
    values = [row[1] for row in history if row[0] >= end_time - window - 1e-9]
    return (max(values) - min(values)) / max(abs(history[-1][1]), 1e-30)


def foam_list(path: Path) -> str:
    text = path.read_text(errors="replace")
    text = re.sub(r"/\*.*?\*/", "", text, flags=re.DOTALL)
    text = re.sub(r"//[^\n]*", "", text)
    body = text[text.index("}", text.index("FoamFile")) + 1:]
    match = re.search(r"(\d+)\s*\(", body)
    if not match:
        raise RuntimeError(f"could not parse {path}")
    return body[match.end():]


def time_name(run_dir: Path, time: float) -> str:
    for candidate in run_dir.iterdir():
        try:
            if candidate.is_dir() and abs(float(candidate.name) - time) < 1e-8:
                return candidate.name
        except ValueError:
            continue
    raise RuntimeError(f"no time directory t = {time:g} in {run_dir}")


def interface_vertices(run_dir: Path, time: float) -> list[tuple]:
    """Interface vertices on the back plane as (x0, y0, x, y).

    The vertices are ordered by their undeformed angle about the centre of
    the half cylinder, from the upstream end.
    """
    mesh = run_dir / "constant" / "fluid" / "polyMesh"
    boundary = (mesh / "boundary").read_text()
    match = re.search(
        r"\binterface\s*\{[^}]*?nFaces\s+(\d+);\s*startFace\s+(\d+);", boundary
    )
    if not match:
        raise RuntimeError(f"no interface patch in {mesh / 'boundary'}")
    n_faces, start = int(match.group(1)), int(match.group(2))
    faces = re.findall(r"\d+\(([^()]*)\)", foam_list(mesh / "faces"))
    ids = set()
    for face in faces[start:start + n_faces]:
        ids.update(int(value) for value in face.split())

    def read_points(path: Path) -> list[tuple[float, ...]]:
        return [
            tuple(float(value) for value in point.split())
            for point in re.findall(r"\(([^()]*)\)", foam_list(path))
        ]

    original = read_points(mesh / "points")
    deformed = read_points(
        run_dir / time_name(run_dir, time) / "fluid" / "polyMesh" / "points"
    )
    z_min = min(original[index][2] for index in ids)
    selected = [
        original[index][:2] + deformed[index][:2]
        for index in ids
        if abs(original[index][2] - z_min) < 1e-9
    ]
    selected.sort(key=lambda row: -math.atan2(row[1] + 0.5, row[0] - 1.5))
    return selected


def read_curve(path: Path) -> list[tuple[float, float]]:
    """Read a published interface, closed by the two fixed base corners.

    The published curves stop short of the floor, so the fixed corners of
    the half cylinder are added to give the interface vertices next to the
    base a segment to be measured against.
    """
    points = [(1.0, -0.5)]
    for line in path.read_text().splitlines():
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        x, y = (float(value) for value in line.split(",")[:2])
        points.append((x, y))
    points.append((2.0, -0.5))
    return points


def distance_to_curve(point, curve) -> float:
    best = math.inf
    for (ax, ay), (bx, by) in zip(curve, curve[1:]):
        dx, dy = bx - ax, by - ay
        length = dx * dx + dy * dy
        s = 0.0 if length == 0 else max(
            0.0, min(1.0, ((point[0] - ax) * dx + (point[1] - ay) * dy) / length)
        )
        best = min(best, math.hypot(point[0] - ax - s * dx, point[1] - ay - s * dy))
    return best


def shape_difference(vertices, curve) -> dict:
    """Distance of the deformed interface vertices from a reference curve.

    The distance is normalised by the largest interface displacement, so the
    measure is the shape error relative to the deformation itself. The two
    fixed base corners, where both curves coincide, are left out.
    """
    inner = [row for row in vertices if row[1] > -0.5 + 1e-6]
    distances = [distance_to_curve(row[2:], curve) for row in inner]
    scale = max(math.hypot(row[2] - row[0], row[3] - row[1]) for row in inner)
    return {
        "max": max(distances) / scale,
        "rms": math.sqrt(sum(d * d for d in distances) / len(distances)) / scale,
        "scale_m": scale,
    }


def observed_order(coarse: float, medium: float, fine: float,
                   ratio: float = 2.0) -> float:
    first, second = abs(medium - coarse), abs(fine - medium)
    if first == 0 or second == 0:
        return math.nan
    return math.log(first / second) / math.log(ratio)


def solver_statistics(run_dir: Path) -> dict:
    log_text = (run_dir / "log.solids4Foam").read_text(errors="replace")
    clock = re.findall(r"ClockTime\s*=\s*([0-9.eE+-]+)", log_text)
    iterations: dict[str, int] = {}
    residual_file = run_dir / "postProcessing" / "fsiResiduals.dat"
    if residual_file.exists():
        for line in residual_file.read_text().splitlines()[1:]:
            fields = line.split()
            if len(fields) >= 2:
                iterations[fields[0]] = max(iterations.get(fields[0], 0),
                                            int(fields[1]))
    counts = list(iterations.values())
    return {
        "clock_time_s": float(clock[-1]) if clock else math.nan,
        "fsi_iterations_mean": sum(counts) / len(counts) if counts else math.nan,
        "fsi_iterations_max": max(counts) if counts else 0,
    }


def count_cells(run_dir: Path, region: str) -> int:
    owner = run_dir / "constant" / region / "polyMesh" / "owner"
    match = re.search(r"nCells:\s*(\d+)", owner.read_text(errors="replace"))
    if not match:
        raise RuntimeError(f"could not read the {region} cell count in {owner}")
    return int(match.group(1))


# ---------------------------------------------------------------------------
# Studies
# ---------------------------------------------------------------------------


def run_mesh_study(args, reference: dict, lines: list[str]) -> bool:
    spec = reference["mesh"]
    if args.levels:
        factors = [float(value) for value in args.levels.split(",")]
    elif args.quick:
        factors = [float(value) for value in spec["quick_factors"]]
    else:
        factors = [float(value) for value in spec["factors"]]
    # The observed order and the extrapolation assume a refinement ratio of 2
    if any(abs(right / left - 2.0) > 1e-9 for left, right in zip(factors, factors[1:])):
        raise SystemExit("successive mesh refinement factors must double, "
                         "e.g. 0.5,1,2,4")
    end_time = float(spec["end_time"])
    window = float(spec["steady_window"])
    liu = reference["liu"]
    curve_t1 = read_curve(REFERENCE_DIR / liu["interface_t1"])
    curve_steady = read_curve(REFERENCE_DIR / liu["interface_steady"])

    rows = []
    failed = False
    for factor in factors:
        delta_t = float(spec["delta_t"].get(f"{factor:g}", spec["delta_t"]["default"]))
        cores = args.cores if factor > 1 else 1
        name = f"mesh_x{factor:g}"
        try:
            run_dir = run_case(name, factor, delta_t, end_time, cores, args.reuse)
        except RuntimeError as error:
            print(f"ERROR: {error}", file=sys.stderr)
            if not args.keep_going:
                return False
            # A failed level fails the study, and leaves a gap in the
            # refinement sequence, so the remaining levels are reported only
            failed = True
            continue
        history = displacement_history(run_dir)
        _, dx_t1, dy_t1 = value_at(history, 1.0)
        _, dx, dy = history[-1]
        shape_t1 = shape_difference(interface_vertices(run_dir, 1.0), curve_t1)
        shape_steady = shape_difference(
            interface_vertices(run_dir, end_time), curve_steady
        )
        row = {
            "factor": factor,
            "cells_fluid": count_cells(run_dir, "fluid"),
            "cells_solid": count_cells(run_dir, "solid"),
            "delta_t_s": delta_t,
            "cores": cores,
            "dx_t1_m": dx_t1,
            "dy_t1_m": dy_t1,
            "dx_steady_m": dx,
            "dy_steady_m": dy,
            "steady_spread": steady_spread(history, window),
            "shape_t1_max": shape_t1["max"],
            "shape_t1_rms": shape_t1["rms"],
            "shape_steady_max": shape_steady["max"],
            "shape_steady_rms": shape_steady["rms"],
        }
        row.update(solver_statistics(run_dir))
        rows.append(row)

    if len(rows) < 2:
        print("ERROR: fewer than two mesh levels completed", file=sys.stderr)
        return False

    POST_DIR.mkdir(parents=True, exist_ok=True)
    with (POST_DIR / "mesh_convergence.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)
    write_interfaces(rows, end_time)

    acceptance = spec["acceptance"]
    passed = not failed and all(
        row["steady_spread"] <= float(spec["steady_tolerance"]) for row in rows
    )
    finest = rows[-1]
    dx_ref = float(liu["top_point_steady_m"][0])
    dy_ref = float(liu["top_point_steady_m"][1])
    dx_error = abs(finest["dx_steady_m"] - dx_ref) / abs(dx_ref)
    dy_error = abs(finest["dy_steady_m"] - dy_ref) / abs(dy_ref)
    order = math.nan
    extrapolated = math.nan
    changes = [abs(b["dx_steady_m"] - a["dx_steady_m"]) for a, b in zip(rows, rows[1:])]
    if len(rows) >= 3:
        order = observed_order(*(row["dx_steady_m"] for row in rows[-3:]))
        if math.isfinite(order) and order > 0:
            extrapolated = finest["dx_steady_m"] + (
                finest["dx_steady_m"] - rows[-2]["dx_steady_m"]
            ) / (2.0 ** order - 1.0)
    if not args.quick:
        passed = (
            passed
            and len(rows) >= 3
            and all(right < left for left, right in zip(changes, changes[1:]))
            and math.isfinite(order)
            and order >= float(acceptance["minimum_order"])
            and abs(finest["dx_steady_m"] - extrapolated) / abs(extrapolated)
            <= float(acceptance["finest_vs_extrapolated_tolerance"])
            and dx_error <= float(acceptance["top_dx_tolerance"])
            and dy_error <= float(acceptance["top_dy_tolerance"])
            and finest["shape_steady_max"]
            <= float(acceptance["shape_steady_max_tolerance"])
            and finest["shape_t1_max"] <= float(acceptance["shape_t1_max_tolerance"])
        )

    lines += [
        "## Mesh study",
        "",
        "- Levels (refinement factors): "
        + ", ".join(f"{row['factor']:g}" for row in rows),
        f"- End time: {end_time:g} s; steady spread of u_x over the last "
        f"{window:g} s: " + ", ".join(f"{r['steady_spread']:.2e}" for r in rows),
        "- Steady top-point u_x (m): "
        + ", ".join(f"{r['dx_steady_m']:.6g}" for r in rows),
        "- Change between successive levels (m): "
        + ", ".join(f"{change:.3g}" for change in changes),
        f"- Observed order: {order:.3f}; extrapolated u_x: {extrapolated:.6g} m",
        f"- Finest u_x = {finest['dx_steady_m']:.6g} m vs Liu {dx_ref:g} m "
        f"({100 * dx_error:.2f}%); u_y = {finest['dy_steady_m']:.6g} m vs "
        f"Liu {dy_ref:g} m ({100 * dy_error:.1f}%)",
        f"- Finest interface shape vs Liu, max distance / max displacement: "
        f"t = 1 s {finest['shape_t1_max']:.3f}, steady "
        f"{finest['shape_steady_max']:.3f}",
        "- FSI iterations per step (mean/max): "
        + ", ".join(
            f"{r['fsi_iterations_mean']:.1f}/{r['fsi_iterations_max']}" for r in rows
        ),
        "- Clock time (s): " + ", ".join(f"{r['clock_time_s']:.0f}" for r in rows),
        f"- Result: {'PASS' if passed else 'FAIL'}",
        "",
    ]
    return passed


def write_interfaces(rows, end_time) -> None:
    """Interface shapes for the comparison plot."""
    out = POST_DIR / "interfaces"
    out.mkdir(parents=True, exist_ok=True)
    for row in rows:
        run_dir = WORK_DIR / f"mesh_x{row['factor']:g}"
        for label, time in (("t1", 1.0), ("steady", end_time)):
            vertices = interface_vertices(run_dir, time)
            (out / f"mesh_x{row['factor']:g}_{label}.csv").write_text(
                "# x (m), y (m)\n"
                + "".join(f"{v[2]:.6f}, {v[3]:.6f}\n" for v in vertices)
            )


def run_time_study(args, reference: dict, lines: list[str]) -> bool:
    spec = reference["time"]
    delta_ts = [float(value) for value in
                (spec["quick_delta_t"] if args.quick else spec["delta_t"])]
    end_time = float(spec["end_time"])
    rows = []
    for delta_t in delta_ts:
        name = f"time_dt{delta_t:g}"
        try:
            run_dir = run_case(name, 1.0, delta_t, end_time, 1, args.reuse)
        except RuntimeError as error:
            print(f"ERROR: {error}", file=sys.stderr)
            return False
        _, dx, dy = value_at(displacement_history(run_dir), end_time)
        points = [row[2:] for row in interface_vertices(run_dir, end_time)]
        row = {"delta_t_s": delta_t, "dx_t1_m": dx, "dy_t1_m": dy}
        row.update(solver_statistics(run_dir))
        rows.append((row, points))

    # Successive differences need no reference solution; the interface
    # difference is the l2 norm over the interface vertices, as in Liu et al.
    table = []
    for index, (row, points) in enumerate(rows):
        if index:
            previous_row, previous_points = rows[index - 1]
            row["dx_change_m"] = abs(row["dx_t1_m"] - previous_row["dx_t1_m"])
            row["interface_change_m"] = math.sqrt(sum(
                (a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2
                for a, b in zip(points, previous_points)
            ))
        else:
            row["dx_change_m"] = math.nan
            row["interface_change_m"] = math.nan
        table.append(row)
    orders = {}
    for metric in ("dx_change_m", "interface_change_m"):
        values = [row[metric] for row in table[1:]]
        orders[metric] = [
            math.log(a / b) / math.log(2.0) if a > 0 and b > 0 else math.nan
            for a, b in zip(values, values[1:])
        ]

    POST_DIR.mkdir(parents=True, exist_ok=True)
    with (POST_DIR / "time_convergence.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=table[0].keys())
        writer.writeheader()
        writer.writerows(table)

    acceptance = spec["acceptance"]
    passed = True
    if not args.quick:
        order_window = orders["dx_change_m"][: int(acceptance["order_pairs"])]
        interface_window = orders["interface_change_m"][
            : int(acceptance["order_pairs"])
        ]
        passed = (
            len(order_window) == int(acceptance["order_pairs"])
            and all(
                float(acceptance["minimum_order"]) <= value
                <= float(acceptance["maximum_order"])
                for value in order_window + interface_window
            )
            and table[-1]["dx_change_m"] / abs(table[-1]["dx_t1_m"])
            <= float(acceptance["finest_change_tolerance"])
        )

    lines += [
        "## Time-step study",
        "",
        "- Time steps (s): " + ", ".join(f"{row['delta_t_s']:g}" for row in table),
        "- u_x at t = 1 s (m): " + ", ".join(f"{row['dx_t1_m']:.7g}" for row in table),
        "- Change in u_x (m): "
        + ", ".join(f"{row['dx_change_m']:.3g}" for row in table[1:]),
        "- Observed order of u_x: "
        + ", ".join(f"{value:.2f}" for value in orders["dx_change_m"]),
        "- Observed order of the interface position: "
        + ", ".join(f"{value:.2f}" for value in orders["interface_change_m"]),
        "- FSI iterations per step (mean/max): "
        + ", ".join(
            f"{r['fsi_iterations_mean']:.1f}/{r['fsi_iterations_max']}" for r in table
        ),
        f"- Result: {'PASS' if passed else 'FAIL'}",
        "",
    ]
    return passed


def create_plot() -> None:
    plot_script = SCRIPT_DIR / "plotInterfaces.gnuplot"
    if shutil.which("gnuplot") is None or not plot_script.is_file():
        return
    output = POST_DIR / "blobInTreacle_interfaces.png"
    completed = subprocess.run(
        [
            "gnuplot",
            "-e",
            f"interfaces='{POST_DIR / 'interfaces'}'; "
            f"reference='{REFERENCE_DIR}'; output='{output}'",
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
    missing = [
        command
        for command in ("blockMesh", "solids4Foam")
        if shutil.which(command) is None
    ]
    if missing:
        raise SystemExit(f"required command(s) not found: {', '.join(missing)}")

    WORK_DIR.mkdir(parents=True, exist_ok=True)
    POST_DIR.mkdir(parents=True, exist_ok=True)
    lines = [
        "# blobInTreacle verification summary",
        "",
        f"- Mode: {'quick smoke test' if args.quick else 'full verification'}",
        "",
    ]
    passed = True
    if args.study in ("all", "time"):
        passed = run_time_study(args, reference, lines) and passed
    if args.study in ("all", "mesh"):
        passed = run_mesh_study(args, reference, lines) and passed
        if (POST_DIR / "interfaces").is_dir():
            create_plot()
    lines.append(f"Overall: {'PASS' if passed else 'FAIL'}")
    summary = "\n".join(lines) + "\n"
    (POST_DIR / "verification_summary.md").write_text(summary)
    print(summary)
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
