#!/usr/bin/env python3
"""Run opt-in 3dTube verification studies in isolated case copies.

The mesh study refines the tutorial fluid and solid meshes uniformly, halves
the time step with each level, and follows the pressure pulse through the
tube: the radial and axial wall displacement at point A, the arrival of the
wave front at stations along the wall, and the resulting pulse-wave speed. The
coupling study compares the default Robin-Neumann coupling with the
Dirichlet-Neumann IQN-ILS (and optionally Aitken) coupling on the tutorial
mesh.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

SCRIPT = Path(__file__).resolve()
SCRIPT_DIR = SCRIPT.parent
VERIFICATION = SCRIPT.parents[1]
TUTORIAL = VERIFICATION.parent
REFERENCE_DIR = VERIFICATION / "reference"
REFERENCE_FILE = REFERENCE_DIR / "3dTube_verification_references.json"
WORK_ROOT = VERIFICATION / "work"
OUTPUT_ROOT = VERIFICATION / "postProcessing"

# Allrun argument and fluidSolidInterface type of each coupling; None keeps the
# type selected by the tutorial variant.
COUPLINGS = {
    "robin": {"allrun": [], "fsiProperties": "constant/fsiProperties.robin",
              "interface": None},
    "iqnils": {"allrun": ["dirichletNeumann"],
               "fsiProperties": "constant/fsiProperties.pimpleFluid",
               "interface": "IQNILS", "predictor": True},
    "aitken": {"allrun": ["dirichletNeumann"],
               "fsiProperties": "constant/fsiProperties.pimpleFluid",
               "interface": "Aitken", "predictor": False},
}


def fail(message: str) -> None:
    raise RuntimeError(message)


# ---------------------------------------------------------------------------
# Case preparation
# ---------------------------------------------------------------------------

def ignored(directory: str, names: list[str]) -> set[str]:
    """Copy only the inputs a fresh run needs, never a previous result."""
    skip = {"verification", "regressionTests", "postProcessing", "images",
            "case.foam"}
    skip.update(name for name in names if name.startswith("processor"))
    skip.update(name for name in names if name.startswith("log."))
    skip.update(name for name in names if name.endswith(".pdf"))
    directory_path = Path(directory)
    if directory_path == TUTORIAL:
        skip.update(
            name for name in names
            if name != "0" and re.fullmatch(r"[0-9.eE+-]+", name)
        )
    # Skip a generated mesh, but keep a foam-extend blockMeshDict
    if (directory_path.parent.name == "constant" and "polyMesh" in names
            and not (directory_path / "polyMesh" / "blockMeshDict").is_file()):
        skip.add("polyMesh")
    return skip.intersection(names)


def copy_case(name: str) -> Path:
    destination = WORK_ROOT / name
    if destination.exists():
        shutil.rmtree(destination)
    shutil.copytree(TUTORIAL, destination, symlinks=True, ignore=ignored)
    return destination


def replace_entry(path: Path, key: str, value: str,
                  top_level: bool = False) -> None:
    """Set one dictionary entry; top_level ignores indented sub-entries."""
    text = path.read_text()
    indent = "" if top_level else r"\s*"
    pattern = rf"^({indent}{re.escape(key)}\s+)[^;]+;"
    text, count = re.subn(pattern, rf"\g<1>{value};", text, flags=re.MULTILINE)
    if count != 1:
        fail(f"Expected one '{key}' entry in {path}, found {count}")
    path.write_text(text)


def block_mesh_dict(case: Path, region: str) -> Path:
    """Locate a region's blockMeshDict in the OpenFOAM or foam-extend layout."""
    for path in (case / "system" / region / "blockMeshDict",
                 case / "constant" / region / "polyMesh" / "blockMeshDict"):
        if path.is_file():
            return path
    fail(f"No blockMeshDict for region {region} in {case}")


def refine_mesh(path: Path, factor: int) -> None:
    """Multiply every block cell-count triple, keeping the block topology."""
    pattern = re.compile(
        r"(hex\s*\([^)]*\)\s*)\(\s*(\d+)\s+(\d+)\s+(\d+)\s*\)"
    )

    def scale(match: re.Match[str]) -> str:
        counts = [int(match.group(index)) * factor for index in (2, 3, 4)]
        return match.group(1) + "(" + " ".join(map(str, counts)) + ")"

    text, count = pattern.subn(scale, path.read_text())
    if not count:
        fail(f"No block cell counts found in {path}")
    path.write_text(text)


def cell_count(path: Path) -> int:
    pattern = re.compile(r"hex\s*\([^)]*\)\s*\(\s*(\d+)\s+(\d+)\s+(\d+)\s*\)")
    return sum(
        math.prod(int(value) for value in match.groups())
        for match in pattern.finditer(path.read_text())
    )


def configure_time_scheme(case: Path, scheme: str) -> None:
    """Apply one time scheme to the fluid and to both solid derivatives."""
    for path, sections in (
        (case / "system/fluid/fvSchemes", ("ddtSchemes",)),
        (case / "system/solid/fvSchemes", ("d2dt2Schemes", "ddtSchemes")),
    ):
        text = path.read_text()
        for section in sections:
            text, count = re.subn(
                rf"({section}\s*\{{\s*default\s+)[^;]+;",
                rf"\g<1>{scheme};",
                text,
                count=1,
            )
            if count != 1:
                fail(f"Expected one {section}/default entry in {path}")
        path.write_text(text)


def configure_solid_preconditioner(case: Path, preconditioner: str) -> None:
    """Replace the tutorial block-Jacobi LU preconditioner, if requested.

    The LU factorisation of the whole solid block grows quickly with the mesh,
    so the finer levels use an algebraic multigrid preconditioner. The
    preconditioner only changes how the solid linear system is solved, not
    the converged solution.
    """
    if preconditioner == "lu":
        return
    path = case / "system/solid/fvSolution"
    text, count = re.subn(
        r"pc_type\s+bjacobi\s*;\s*sub_pc_type\s+lu\s*;",
        "pc_type hypre;\n            pc_hypre_type boomeramg;\n"
        '            pc_hypre_boomeramg_max_iter "1";\n'
        '            pc_hypre_boomeramg_strong_threshold "0.7";',
        path.read_text(),
    )
    if count != 1:
        fail(f"Expected one block-Jacobi LU preconditioner in {path}")
    path.write_text(text)


def configure_tight_tolerances(case: Path, coupling: str) -> None:
    """Tighten every iterative tolerance, to bound the iterative error.

    A diagnostic, not the benchmark setup: the fluid linear solvers go from
    relTol 1e-3 to 1e-6 (and the pressure tolerance from 1e-6 to 1e-9), the
    PIMPLE residual control from 1e-5 to 1e-7, the solid Newton tolerance
    from 1e-6 to 1e-9 with a linear tolerance of 1e-8, and the FSI
    outerCorrTolerance from 1e-6 to 1e-8.
    """
    fluid = case / "system/fluid/fvSolution"
    text = fluid.read_text()
    text, count = re.subn(r"(relTol\s+)(1e-03|0\.001)\s*;", r"\g<1>1e-06;", text)
    if count != 3:
        fail(f"Expected three relTol 1e-3 entries in {fluid}, found {count}")
    text, count = re.subn(r"(tolerance\s+)1e-06\s*;", r"\g<1>1e-09;", text)
    if count != 1:
        fail(f"Expected one pressure tolerance 1e-6 in {fluid}")
    text, count = re.subn(r"((?:relTol|tolerance)\s+)1e-5\s*;", r"\g<1>1e-7;", text)
    if count != 2:
        fail(f"Expected the PIMPLE residualControl in {fluid}")
    fluid.write_text(text)
    solid = case / "system/solid/fvSolution"
    text, count = re.subn(
        r'snes_rtol\s+"1e-6"\s*;',
        'snes_rtol "1e-9";\n            ksp_rtol "1e-8";', solid.read_text()
    )
    if count != 1:
        fail(f"Expected one snes_rtol in {solid}")
    solid.write_text(text)
    properties = case / COUPLINGS[coupling]["fsiProperties"]
    text, count = re.subn(
        r"(outerCorrTolerance\s+)1e-6\s*;", r"\g<1>1e-8;",
        properties.read_text()
    )
    if count < 1:
        fail(f"Expected outerCorrTolerance 1e-6 in {properties}")
    properties.write_text(text)


def station_name(z: float) -> str:
    return f"verificationWallZ{round(z * 1.0e4):03d}"


def add_monitors(case: Path, reference: dict, probe_z_shift: float = 0.0
                 ) -> None:
    """Add wall-displacement and axis-pressure monitors along the tube."""
    geometry = reference["geometry"]
    radius = float(geometry["inner_radius_m"])
    entries = []
    for z in reference["stations_z_m"]:
        entries.append(
            f"    {station_name(z)}\n    {{\n"
            "        type            solidPointDisplacement;\n"
            "        region          solid;\n"
            f"        point           (0 {radius:.10g} {z:.10g});\n"
            "    }\n"
        )
    offset = float(reference["pressure_probe_offset_m"])
    locations = " ".join(
        f"({offset:.10g} {offset:.10g} {z + probe_z_shift:.10g})"
        for z in reference["stations_z_m"]
    )
    entries.append(
        "    verificationAxisPressure\n    {\n"
        "        type            probes;\n"
        '        libs            ("libsampling.so");\n'
        "        region          fluid;\n"
        "        writeControl    timeStep;\n"
        "        writeInterval   1;\n"
        "        fields          (p);\n"
        f"        probeLocations  ({locations});\n"
        "    }\n"
    )
    control = case / "system/controlDict"
    text, count = re.subn(
        r"(^functions\s*\{\s*\n)",
        lambda match: match.group(1) + "".join(entries),
        control.read_text(),
        count=1,
        flags=re.MULTILINE,
    )
    if count != 1:
        fail(f"Could not find the functions dictionary in {control}")
    control.write_text(text)


def configure_coupling(case: Path, coupling: str) -> None:
    """Select the interface type and make it write its residual history."""
    spec = COUPLINGS[coupling]
    if spec["interface"] is None:
        return
    path = case / spec["fsiProperties"]
    replace_entry(path, "fluidSolidInterface", spec["interface"], True)
    block = spec["interface"] + "Coeffs"
    text = path.read_text()
    match = re.search(rf"^{block}\s*\{{(.*?)^\}}", text,
                      flags=re.MULTILINE | re.DOTALL)
    if not match:
        fail(f"No {block} dictionary in {path}")
    # The fluid-solid interface predictor avoids the first-iterate added-mass
    # spike of IQN-ILS (issue #489). Aitken predicts the solid itself
    # (predictSolid), and the Robin-Neumann tutorial setup is used as it is.
    added = ""
    keys = ["writeResidualsToFile"]
    if spec["predictor"]:
        keys.append("predictor")
    for key in keys:
        if not re.search(rf"^\s*{key}\s", match.group(1), re.MULTILINE):
            added += f"    {key} yes;\n"
        elif not re.search(rf"^\s*{key}\s+(yes|on|true)\s*;",
                           match.group(1), re.MULTILINE):
            fail(f"{block}/{key} in {path} is not enabled")
    text = text[:match.end(1)] + added + text[match.end(1):]
    path.write_text(text)


def configure_parallel(case: Path, cores: int) -> None:
    """Run the solver of the copy on several ranks.

    The tutorial Allrun runs in serial only, so the copy's Allrun is changed
    to decompose both regions and run the solver in parallel. The monitored
    histories are written by the master rank, so no reconstruction is needed.
    """
    if cores <= 1:
        return
    for dictionary in (
        case / "system/decomposeParDict",
        case / "system/fluid/decomposeParDict",
        case / "system/solid/decomposeParDict",
    ):
        replace_entry(dictionary, "numberOfSubdomains", str(cores))
    allrun = case / "Allrun"
    serial = "solids4Foam::runApplication solids4Foam\n"
    text = allrun.read_text()
    if text.count(serial) != 1:
        fail(f"Expected one serial solver call in {allrun}")
    allrun.write_text(text.replace(
        serial,
        "solids4Foam::runApplication -s fluid decomposePar -region fluid\n"
        "solids4Foam::runApplication -s solid decomposePar -region solid\n"
        "solids4Foam::runParallel solids4Foam || exit 1\n",
    ))


def prepare_case(name: str, coupling: str, factor: int, delta_t: float,
                 end_time: float, args: argparse.Namespace,
                 reference: dict, cores: int) -> Path:
    case = copy_case(name)
    refine_mesh(block_mesh_dict(case, "fluid"), factor)
    refine_mesh(block_mesh_dict(case, "solid"), factor)
    control = case / "system/controlDict"
    replace_entry(control, "deltaT", f"{delta_t:.10g}", True)
    replace_entry(control, "endTime", f"{end_time:.10g}", True)
    # Fields are only needed at the end; the monitors write every time step.
    replace_entry(control, "writeInterval", f"{end_time:.10g}", True)
    configure_time_scheme(case, args.time_scheme)
    configure_solid_preconditioner(case, args.solid_preconditioner)
    configure_coupling(case, coupling)
    if args.tight_tolerances:
        configure_tight_tolerances(case, coupling)
    add_monitors(case, reference, args.probe_z_shift)
    configure_parallel(case, cores)
    return case


# ---------------------------------------------------------------------------
# Running
# ---------------------------------------------------------------------------

def solver_completed(case: Path) -> bool:
    log = case / "log.solids4Foam"
    return log.is_file() and bool(
        re.search(r"^End\s*$", log.read_text(errors="replace"), re.MULTILINE)
    )


def run_case(case: Path, label: str, coupling: str) -> None:
    allrun = case / "Allrun"
    if not allrun.is_file():
        fail(f"Missing {allrun}")
    log = case / "log.Allverify"
    print(f"Running {label} in {case}", flush=True)
    with log.open("w") as handle:
        result = subprocess.run(
            [str(allrun), *COUPLINGS[coupling]["allrun"]],
            cwd=case, stdout=handle, stderr=subprocess.STDOUT, text=True,
        )
    if result.returncode:
        fail(f"{label} failed; see {log}")
    check_solver_log(case, label)


def check_solver_log(case: Path, label: str) -> None:
    """Reject a run that aborted, since the Allrun returns zero regardless."""
    solver_log = case / "log.solids4Foam"
    if not solver_log.is_file():
        fail(f"{label} did not create {solver_log}")
    text = solver_log.read_text(errors="replace")
    if re.search(
        r"FOAM FATAL|FOAM aborting|PRTE ERROR|No sockets were able|^ERROR$",
        text, re.MULTILINE,
    ):
        fail(f"{label} failed inside Allrun; see {solver_log}")
    if not re.search(r"^End\s*$", text, re.MULTILINE):
        fail(f"{label} did not run to completion; see {solver_log}")


def solver_ranks(case: Path) -> int:
    """Number of MPI ranks the solver of a completed run actually used."""
    match = re.search(
        r"^nProcs\s*:\s*(\d+)",
        (case / "log.solids4Foam").read_text(errors="replace"), re.MULTILINE,
    )
    return int(match.group(1)) if match else 1


def clock_time(case: Path) -> float:
    matches = re.findall(
        r"ClockTime\s*=\s*([0-9.eE+-]+)",
        (case / "log.solids4Foam").read_text(errors="replace"),
    )
    return float(matches[-1]) if matches else math.nan


# ---------------------------------------------------------------------------
# Post-processing
# ---------------------------------------------------------------------------

def numeric_rows(path: Path) -> list[list[float]]:
    """Return every data row, rejecting malformed or non-finite ones.

    Comment lines and a single leading header line (as in fsiResiduals.dat)
    are skipped; any other row that is not entirely finite numbers, such as a
    row cut short by an aborted run, is an error rather than being dropped.
    """
    rows = []
    header_allowed = True
    for number, line in enumerate(
        path.read_text(errors="replace").splitlines(), 1
    ):
        fields = line.replace("(", " ").replace(")", " ").split()
        if not fields or fields[0].startswith("#"):
            continue
        try:
            values = [float(field) for field in fields]
        except ValueError:
            if header_allowed:
                header_allowed = False
                continue
            fail(f"{path}:{number} is not a numeric data row: {line.strip()}")
        header_allowed = False
        if not all(math.isfinite(value) for value in values):
            fail(f"{path}:{number} contains non-finite values")
        rows.append(values)
    return rows


def last_per_time(rows: list[list[float]]) -> list[list[float]]:
    """Keep the final row written for each time, in time order."""
    by_time: dict[float, list[float]] = {}
    for row in rows:
        by_time[row[0]] = row
    return [by_time[time] for time in sorted(by_time)]


def dictionary_scalar(path: Path, key: str) -> float | None:
    """Read a top-level scalar entry of an OpenFOAM dictionary."""
    match = re.search(
        rf"^{re.escape(key)}\s+([-+0-9.eE]+)\s*;",
        path.read_text(), flags=re.MULTILINE,
    )
    return float(match.group(1)) if match else None


def require_complete(path: Path, rows: list[list[float]], end_time: float,
                     columns: int) -> None:
    """Require one complete, finite row per time step up to the end time.

    The time step is read from the copy's controlDict, so a history with a
    gap or a missing tail is rejected rather than evaluated.
    """
    short = [row for row in rows if len(row) < columns]
    if short:
        fail(f"{path} has {len(short)} row(s) with fewer than {columns} "
             f"columns, starting at t={short[0][0]:g}")
    control = path
    while control.name != "postProcessing":
        control = control.parent
    delta_t = dictionary_scalar(control.parent / "system/controlDict", "deltaT")
    if delta_t is None or delta_t <= 0.0:
        fail(f"Could not read deltaT for {path}")
    expected = round(end_time / delta_t)
    times = [row[0] for row in rows]
    if not times or not math.isclose(
        times[-1], end_time, rel_tol=1e-6, abs_tol=1e-12
    ):
        reached = times[-1] if times else 0.0
        fail(f"{path} ends at t={reached:g}, not at the end time "
             f"t={end_time:g}")
    # The function objects write the time with the case's timePrecision
    # (6 significant digits), which does not resolve dt = 6.25e-6 s beyond
    # t = 10 ms; each time must therefore round to the next multiple of dt,
    # and is then replaced by that exact multiple, so that the histories
    # are equally spaced for the peak and arrival interpolation.
    steps = [round(time / delta_t) for time in times]
    if len(times) != expected or steps != list(
        range(steps[0], steps[0] + expected)
    ) or any(
        abs(time - step * delta_t) > 1e-5 * abs(time) + 1e-3 * delta_t
        for time, step in zip(times, steps)
    ):
        fail(f"{path} has {len(times)} time levels, expected {expected} "
             f"at dt = {delta_t:g}")
    for row, step in zip(rows, steps):
        row[0] = step * delta_t
    if any(not math.isfinite(value) for row in rows for value in row):
        fail(f"{path} contains non-finite values")


def displacement_history(case: Path, name: str, end_time: float
                         ) -> list[tuple[float, float, float]]:
    """Return (t, radial, axial) wall displacement of one monitor.

    Every monitor lies in the x = 0 symmetry plane, where the radial
    direction is y.
    """
    candidates = sorted(
        case.glob(f"postProcessing/**/solidPointDisplacement_{name}.dat")
    )
    if not candidates:
        fail(f"Displacement monitor '{name}' not found in {case}")
    rows = last_per_time(numeric_rows(candidates[0]))
    if not rows:
        fail(f"No numeric rows in {candidates[0]}")
    # Columns: time, then the x, y and z displacement and its magnitude
    require_complete(candidates[0], rows, end_time, 5)
    return [(row[0], row[2], row[3]) for row in rows]


def pressure_histories(case: Path, n_stations: int, end_time: float,
                       density: float) -> list[list[tuple[float, float]]]:
    """Return the axis pressure (Pa) history of every station."""
    candidates = sorted(case.glob("postProcessing/**/verificationAxisPressure/**/p"))
    if not candidates:
        fail(f"Axis pressure probes not found in {case}")
    rows = last_per_time(numeric_rows(candidates[0]))
    if not rows or any(len(row) != n_stations + 1 for row in rows):
        fail(f"Unexpected pressure probe columns in {candidates[0]}")
    require_complete(candidates[0], rows, end_time, n_stations + 1)
    # The fluid pressure is kinematic.
    return [
        [(row[0], density * row[index + 1]) for row in rows]
        for index in range(n_stations)
    ]


def peak(history: list[tuple[float, float]], sign: float = 1.0
         ) -> tuple[float, float]:
    """Return (time, value) of the extreme of sign*value.

    A parabola through the extreme sample and its neighbours refines both the
    time and the value, so the result is not quantised to the time step.
    """
    index = max(range(len(history)), key=lambda i: sign * history[i][1])
    time, value = history[index]
    if 0 < index < len(history) - 1:
        (t0, v0), (t1, v1), (t2, v2) = history[index - 1:index + 2]
        denominator = (v0 - 2.0 * v1 + v2)
        if denominator != 0.0 and math.isclose(t1 - t0, t2 - t1, rel_tol=1e-6):
            shift = 0.5 * (v0 - v2) / denominator
            time = t1 + shift * (t1 - t0)
            value = v1 - 0.25 * (v0 - v2) * shift
    return time, value


def first_peak(history: list[tuple[float, float]], fraction: float
               ) -> tuple[float, float]:
    """Return the first local maximum that exceeds fraction of the maximum.

    The wave front is followed by reflections from the clamped outlet, which
    can exceed the incident peak near the outlet; the incident peak is the
    first significant maximum.
    """
    maximum = max(value for _, value in history)
    threshold = fraction * maximum
    for index in range(1, len(history) - 1):
        value = history[index][1]
        if (value >= threshold and value >= history[index - 1][1]
                and value > history[index + 1][1]):
            return peak(history[max(0, index - 1):index + 2])
    return peak(history)


def arrival_time(history: list[tuple[float, float]], level: float) -> float:
    """Return the first time the signal rises through level (interpolated)."""
    for (t0, v0), (t1, v1) in zip(history, history[1:]):
        if v0 < level <= v1:
            return t0 + (level - v0) / (v1 - v0) * (t1 - t0)
    return math.nan


def front_arrival(history: list[tuple[float, float]], reference: dict
                  ) -> tuple[float, float, float]:
    """Return (arrival time, first-peak time, first-peak value)."""
    spec = reference["wave_front"]
    peak_time, peak_value = first_peak(history, float(spec["peak_fraction"]))
    leading = [(t, v) for t, v in history if t <= peak_time]
    arrival = arrival_time(leading, float(spec["arrival_fraction"]) * peak_value)
    return arrival, peak_time, peak_value


def fit_speed(positions: list[float], times: list[float]) -> float:
    """Least-squares slope dz/dt of the arrival times."""
    # Every station must have a front arrival; a partial fit is rejected.
    if len(times) < 2 or not all(math.isfinite(t) for t in times):
        return math.nan
    pairs = list(zip(positions, times))
    mean_z = sum(z for z, _ in pairs) / len(pairs)
    mean_t = sum(t for _, t in pairs) / len(pairs)
    s_tt = sum((t - mean_t) ** 2 for _, t in pairs)
    s_tz = sum((t - mean_t) * (z - mean_z) for z, t in pairs)
    return s_tz / s_tt if s_tt > 0.0 else math.nan


def extract(case: Path, reference: dict, end_time: float) -> dict:
    stations = [float(z) for z in reference["stations_z_m"]]
    point_a = float(reference["point_A_z_m"])
    fit_range = reference["wave_front"]["fit_z_range_m"]
    in_fit = [fit_range[0] - 1e-9 <= z <= fit_range[1] + 1e-9 for z in stations]

    wall = {z: displacement_history(case, station_name(z), end_time)
            for z in stations}
    radial = {z: [(t, r) for t, r, _ in wall[z]] for z in stations}
    axial_a = [(t, a) for t, _, a in wall[point_a]]
    pressures = pressure_histories(
        case, len(stations), end_time, float(reference["fluid"]["rho_kg_m3"])
    )

    row: dict[str, float] = {}
    row["ur_max_A_m"], row["t_ur_max_A_s"] = reversed(peak(radial[point_a]))
    # The incident trough of the axial displacement (about 4.7 ms)
    incident_end = float(reference["incident_axial_trough_before_s"])
    row["uz_min_A_m"], row["t_uz_min_A_s"] = reversed(peak(
        [(time, value) for time, value in axial_a if time <= incident_end], -1.0
    ))
    # The second, reflected trough of the radial displacement (about 17 ms)
    late_start = float(reference["late_minimum_after_s"])
    late = [(time, value) for time, value in radial[point_a] if time >= late_start]
    if len(late) > 2 and end_time >= float(reference["end_time_s"]) - 1e-12:
        row["ur_min_late_A_m"], row["t_ur_min_late_A_s"] = reversed(
            peak(late, -1.0)
        )
    else:
        row["ur_min_late_A_m"] = row["t_ur_min_late_A_s"] = math.nan
    wall_arrivals, pressure_arrivals = [], []
    for index, z in enumerate(stations):
        arrival, peak_time, peak_value = front_arrival(radial[z], reference)
        wall_arrivals.append(arrival)
        # The probes sample the value of the cell that contains them, which
        # steps as the fluid mesh moves; a fixed level, half the inlet
        # pressure, keeps the pressure arrival independent of those steps.
        _, _, p_peak = front_arrival(pressures[index], reference)
        p_arrival = arrival_time(
            pressures[index],
            float(reference["wave_front"]["arrival_fraction"])
            * float(reference["inlet_pressure_Pa"]),
        )
        pressure_arrivals.append(p_arrival)
        tag = f"z{round(z * 1.0e3):02d}mm"
        row[f"wall_arrival_{tag}_s"] = arrival
        row[f"wall_first_peak_{tag}_m"] = peak_value
        row[f"wall_first_peak_time_{tag}_s"] = peak_time
        row[f"pressure_arrival_{tag}_s"] = p_arrival
        row[f"pressure_first_peak_{tag}_Pa"] = p_peak
        if math.isclose(z, point_a):
            row["t_arrival_A_s"] = arrival
    fit_z = [z for z, keep in zip(stations, in_fit) if keep]
    row["wave_speed_wall_m_s"] = fit_speed(
        fit_z, [t for t, keep in zip(wall_arrivals, in_fit) if keep]
    )
    row["wave_speed_pressure_m_s"] = fit_speed(
        fit_z, [t for t, keep in zip(pressure_arrivals, in_fit) if keep]
    )
    row["_histories"] = {
        "radial_A": radial[point_a], "axial_A": axial_a,
        "pressure_A": pressures[stations.index(point_a)],
    }
    return row


def iteration_summary(case: Path, end_time: float) -> dict[str, float]:
    candidates = sorted(case.glob("postProcessing/**/fsiConvergenceData.dat"))
    if not candidates:
        fail(f"fsiConvergenceData.dat not found in {case}")
    rows = last_per_time(numeric_rows(candidates[0]))
    if not rows:
        fail(f"No numeric rows in {candidates[0]}")
    require_complete(candidates[0], rows, end_time, 2)
    counts = [row[1] for row in rows]
    return {
        "time_steps": len(counts),
        "total_fsi_iterations": sum(counts),
        "mean_fsi_iterations": sum(counts) / len(counts),
        "maximum_fsi_iterations": max(counts),
    }


def coefficient(path: Path, block: str, key: str) -> float:
    """Read one scalar from a coefficients sub-dictionary."""
    match = re.search(
        rf"^{re.escape(block)}\s*\{{(.*?)^\}}",
        path.read_text(), flags=re.MULTILINE | re.DOTALL,
    )
    if not match:
        fail(f"No {block} dictionary in {path}")
    value = re.search(
        rf"^\s*{re.escape(key)}\s+([-+0-9.eE]+)\s*;",
        match.group(1), flags=re.MULTILINE,
    )
    if not value:
        fail(f"Could not read {block}/{key} from {path}")
    return float(value.group(1))


def residual_rows(case: Path, end_time: float) -> tuple[list[str], list[list[float]]]:
    """Return the header and the final row of each time step of fsiResiduals."""
    path = case / "postProcessing/fsiResiduals.dat"
    if not path.is_file():
        fail(f"FSI residual data not found in {case}")
    header = path.read_text(errors="replace").splitlines()[0].split()
    rows = last_per_time(numeric_rows(path))
    if not rows or any(len(row) < 3 for row in rows):
        fail(f"FSI residual columns are missing from {path}")
    require_complete(path, rows, end_time, 3)
    return header, rows


def robin_residual_summary(case: Path, end_time: float) -> dict[str, float]:
    """Require every Robin time step to have met its convergence criteria.

    fsiResiduals.dat has the columns time, outerCorrector, residual,
    robinPressureResidual, robinFluxResidual and robinConvergenceState
    (0 not converged, 1 converged, 2 stalled within the stall tolerance).
    """
    header, rows = residual_rows(case, end_time)
    if "robinConvergenceState" not in header:
        fail(f"fsiResiduals.dat in {case} has no robinConvergenceState column")
    state = header.index("robinConvergenceState")
    if any(len(row) <= state for row in rows):
        fail(f"Robin residual columns are missing in {case}")
    unconverged = [row for row in rows if row[state] <= 0]
    if unconverged:
        fail(
            f"{len(unconverged)} Robin time step(s) in {case} reached "
            "nOuterCorr without satisfying all convergence criteria"
        )
    properties = case / COUPLINGS["robin"]["fsiProperties"]
    return {
        "unconverged_steps": 0,
        "robin_stalled_steps": sum(row[state] == 2 for row in rows),
        "n_outer_corr": coefficient(
            properties, "fixedRelaxationCoeffs", "nOuterCorr"
        ),
        "maximum_final_displacement_residual": max(row[2] for row in rows),
        "maximum_final_pressure_residual": max(row[3] for row in rows),
        "maximum_final_flux_residual": max(row[4] for row in rows),
    }


def dirichlet_neumann_summary(case: Path, coupling: str, end_time: float
                              ) -> dict[str, float]:
    """Require every Dirichlet-Neumann step to reach outerCorrTolerance."""
    _, rows = residual_rows(case, end_time)
    properties = case / COUPLINGS[coupling]["fsiProperties"]
    block = COUPLINGS[coupling]["interface"] + "Coeffs"
    tolerance = coefficient(properties, block, "outerCorrTolerance")
    unconverged = [row for row in rows if row[2] > tolerance]
    if unconverged:
        fail(
            f"{len(unconverged)} {coupling} time step(s) in {case} reached "
            f"nOuterCorr with the FSI residual above {tolerance:g}"
        )
    return {
        "unconverged_steps": 0,
        "n_outer_corr": coefficient(properties, block, "nOuterCorr"),
        "maximum_final_displacement_residual": max(row[2] for row in rows),
    }


def coupling_summary(case: Path, coupling: str, end_time: float
                     ) -> dict[str, float]:
    summary = iteration_summary(case, end_time)
    # The residual evaluated for each step must be that of its last
    # FSI iteration, as recorded independently in fsiConvergenceData.dat.
    _, residuals = residual_rows(case, end_time)
    delta_t = dictionary_scalar(case / "system/controlDict", "deltaT")
    iterations = {
        round(row[0] / delta_t): row[1] for row in last_per_time(numeric_rows(
            next(case.glob("postProcessing/**/fsiConvergenceData.dat"))
        ))
    }
    mismatched = [row[0] for row in residuals
                  if iterations.get(round(row[0] / delta_t)) != row[1]]
    if mismatched:
        fail(f"{len(mismatched)} step(s) in {case} have a final residual row "
             "that is not their last FSI iteration, starting at "
             f"t={mismatched[0]:g}")
    if coupling == "robin":
        summary.update(robin_residual_summary(case, end_time))
    else:
        summary.update(dirichlet_neumann_summary(case, coupling, end_time))
    return summary


def within(value: float, tolerance: float) -> bool:
    """A check passes only for a finite value within the tolerance."""
    return math.isfinite(value) and value <= tolerance


def relative_difference(value: float, reference: float) -> float:
    return abs(value - reference) / abs(reference) if reference else math.inf


def wave_speed_estimates(reference: dict) -> dict[str, float]:
    """Analytical pulse-wave speed estimates for the tube.

    mk_*: Moens-Korteweg, c = sqrt(E h / (2 rho_f R)), for a thin, inviscid,
    long-wave tube without wall inertia; with the inner or the mean radius,
    and with the 1/(1 - nu^2) factor of an axially tethered wall.
    thick_wall: Tukovic et al. (2018) Eq. (35), the thick-wall (Lame)
    distensibility of an incompressible fluid-filled tube.
    thick_wall_inertia: Tukovic et al. (2018) Eq. (36), Eq. (35) corrected for
    the axial stress waves in the wall.
    """
    geometry, solid = reference["geometry"], reference["solid"]
    rho = float(reference["fluid"]["rho_kg_m3"])
    rho_s = float(solid["rho_kg_m3"])
    young, poisson = float(solid["E_Pa"]), float(solid["nu"])
    h = float(geometry["wall_thickness_m"])
    r = float(geometry["inner_radius_m"])
    mean_radius = r + 0.5 * h
    thick_wall = math.sqrt(
        young * h / (2.0 * rho * r)
        / ((h / r) * (1.0 + poisson) + 2.0 * r / (2.0 * r + h))
    )
    return {
        "mk_inner_radius": math.sqrt(young * h / (2.0 * rho * r)),
        "mk_mean_radius": math.sqrt(young * h / (2.0 * rho * mean_radius)),
        "mk_inner_radius_tethered": math.sqrt(
            young * h / (2.0 * rho * r * (1.0 - poisson**2))
        ),
        "thick_wall": thick_wall,
        "thick_wall_inertia": thick_wall * math.sqrt(
            1.0 - poisson**2 / (1.0 - h * rho / (2.0 * r * rho_s))
        ),
    }


# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------

def write_histories(name: str, histories: dict) -> Path:
    """Write the point-A histories of one run for plotting."""
    path = OUTPUT_ROOT / "histories" / f"{name}.csv"
    path.parent.mkdir(parents=True, exist_ok=True)
    pressure = dict(histories["pressure_A"])
    axial = dict(histories["axial_A"])
    with path.open("w") as handle:
        handle.write("time_s,radial_displacement_m,axial_displacement_m,"
                     "axis_pressure_Pa\n")
        for time, radial in histories["radial_A"]:
            handle.write(
                f"{time:.10g},{radial:.10g},{axial.get(time, math.nan):.10g},"
                f"{pressure.get(time, math.nan):.10g}\n"
            )
    return path


def write_csv(path: Path, rows: list[dict]) -> None:
    columns: list[str] = []
    for row in rows:
        columns += [key for key in row if key not in columns and not key.startswith("_")]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def run_gnuplot(script: str, variables: dict[str, str]) -> None:
    if shutil.which("gnuplot") is None:
        print("gnuplot not found; skipping plots (the CSV results are valid)")
        return
    assignments = "; ".join(f"{key}='{value}'" for key, value in variables.items())
    for terminal, extension in (("png", "png"), ("pdf", "pdf")):
        # Paths are passed relative to the verification directory, so the
        # checkout path, which may contain spaces or quotes, never reaches
        # the gnuplot expressions.
        output = OUTPUT_ROOT / f"{variables['stem']}.{extension}"
        relative_output = output.relative_to(VERIFICATION)
        result = subprocess.run(
            ["gnuplot", "-e",
             f"{assignments}; term='{terminal}'; output='{relative_output}'",
             str(SCRIPT_DIR / script)],
            cwd=VERIFICATION, text=True, check=False,
        )
        if result.returncode:
            print(f"WARNING: gnuplot failed for {output}")
        else:
            print(f"Plot: {output}")


def format_value(value: float, digits: int = 5) -> str:
    return f"{value:.{digits}g}" if math.isfinite(value) else "n/a"


# ---------------------------------------------------------------------------
# Studies
# ---------------------------------------------------------------------------

def study_cores(requested: str, study: str, level: int, reference: dict) -> int:
    if requested != "auto":
        return int(requested)
    return int(reference[study]["auto_cores"][str(level)])


SETTINGS_FILE = "verification_settings.json"


def tutorial_fingerprint() -> str:
    """Hash of the tutorial inputs that a run copies."""
    digest = hashlib.sha256()
    for directory in ("0", "constant", "system"):
        for path in sorted((TUTORIAL / directory).rglob("*")):
            if path.is_symlink():
                digest.update(f"{path.relative_to(TUTORIAL)}->"
                              f"{os.readlink(path)}".encode())
            elif path.is_file() and (
                "polyMesh" not in path.parts or path.name == "blockMeshDict"
            ):
                digest.update(str(path.relative_to(TUTORIAL)).encode())
                digest.update(path.read_bytes())
    digest.update((TUTORIAL / "Allrun").read_bytes())
    return digest.hexdigest()


def run_settings(coupling: str, factor: int, delta_t: float, end_time: float,
                 args: argparse.Namespace, reference: dict, cores: int
                 ) -> dict:
    """Everything that defines a run, so a stale copy is never reused."""
    return {
        "coupling": coupling, "refinement": factor,
        "delta_t_s": f"{delta_t:.10g}", "end_time_s": f"{end_time:.10g}",
        "time_scheme": args.time_scheme,
        "solid_preconditioner": args.solid_preconditioner,
        "stations_z_m": reference["stations_z_m"],
        "pressure_probe_offset_m": reference["pressure_probe_offset_m"],
        "interface_predictor": COUPLINGS[coupling].get("predictor", False),
        "inner_radius_m": reference["geometry"]["inner_radius_m"],
        "cores": cores,
        "openfoam": os.environ.get("WM_PROJECT", "") + "-"
        + os.environ.get("WM_PROJECT_VERSION", ""),
        "tutorial_inputs": tutorial_fingerprint(),
        # Recorded only when used, so runs without a shift stay reusable
        **({"pressure_probe_z_shift_m": args.probe_z_shift}
           if args.probe_z_shift else {}),
        **({"tight_tolerances": True} if args.tight_tolerances else {}),
    }


def run_or_reuse(name: str, coupling: str, factor: int, delta_t: float,
                 end_time: float, args: argparse.Namespace, reference: dict,
                 cores: int, label: str) -> Path:
    case = WORK_ROOT / name
    settings = run_settings(coupling, factor, delta_t, end_time, args,
                            reference, cores)
    settings_file = case / SETTINGS_FILE
    if args.reuse and solver_completed(case):
        stored = (json.loads(settings_file.read_text())
                  if settings_file.is_file() else None)
        if stored == settings:
            print(f"Reusing {label} in {case}")
            check_solver_log(case, label)
            return case
        print(f"Not reusing {case}: it was run with other settings")
    case = prepare_case(name, coupling, factor, delta_t, end_time, args,
                        reference, cores)
    run_case(case, label, coupling)
    settings_file.write_text(json.dumps(settings, indent=2) + "\n")
    return case


def observed_order(values: list[float]) -> float:
    """Order from the last three members of a sweep with ratio two."""
    if len(values) < 3:
        return math.nan
    coarse, medium, fine = values[-3:]
    numerator, denominator = coarse - medium, medium - fine
    if numerator == 0.0 or denominator == 0.0 or numerator * denominator < 0.0:
        return math.nan
    return math.log(abs(numerator / denominator), 2.0)


def sweep_levels(args: argparse.Namespace, reference: dict, study: str
                 ) -> list[tuple[int, int, float]]:
    """Return (level, mesh refinement factor, time step) of each member."""
    spec = reference[study]
    levels = args.levels or list(spec["default_levels"])
    if args.quick:
        levels = levels[:1]
    members = []
    for level in levels:
        if study == "mesh":
            factor = 2 ** (level - 1)
            # --delta-t holds the time step fixed, so that the mesh alone is
            # refined; by default it is halved with the mesh spacing.
            delta_t = (args.delta_t if args.delta_t
                       else float(spec["base_delta_t_s"]) / factor)
            members.append((level, factor, delta_t))
        else:
            if level > len(spec["delta_t_s"]):
                fail(f"{study} level {level} is not defined")
            members.append((level, 1, float(spec["delta_t_s"][level - 1])))
    return members


def published_checks(rows: list[dict], reference: dict, study: str
                     ) -> list[tuple[str, bool]]:
    """Compare the finest member with the published values for this study."""
    checks = []
    finest = rows[-1]
    for source, quantities in reference["published_checks"].get(study, {}).items():
        citation = reference["published"][source]["citation"]
        for quantity, spec in quantities.items():
            value, tolerance = float(spec["value"]), float(spec["tolerance"])
            difference = relative_difference(finest[quantity], value)
            finest[f"{quantity}_vs_{source}"] = difference
            checks.append((
                f"{quantity} {format_value(finest[quantity], 4)} vs {citation} "
                f"{value:.4g} ({spec['quality']}): {100 * difference:.1f}% <= "
                f"{100 * tolerance:g}%",
                within(difference, tolerance),
            ))
    return checks


def history_differences(row: dict, reference: dict, study: str
                        ) -> list[str]:
    """Report the largest difference from each published radial history.

    The difference is normalised by the largest published |u_r| and taken
    over the published time range; it is reported, not checked.
    """
    lines = []
    history = row["_histories"]["radial_A"]
    times = [time for time, _ in history]
    for source in reference["published_histories"].get(study, []):
        spec = reference["published"][source]
        published = []
        for line in (REFERENCE_DIR / spec["history"]).read_text().splitlines():
            fields = line.split(",")
            if line.startswith("#") or fields[0] == "time_s":
                continue
            if fields[1] != "nan" and times[0] <= float(fields[0]) <= times[-1]:
                published.append((float(fields[0]), float(fields[1])))
        scale = max(abs(value) for _, value in published)
        worst = 0.0
        for time, value in published:
            index = next(i for i, t in enumerate(times) if t >= time - 1e-12)
            if index == 0:
                computed = history[0][1]
            else:
                (t0, v0), (t1, v1) = history[index - 1], history[index]
                computed = v0 + (v1 - v0) * (time - t0) / (t1 - t0)
            worst = max(worst, abs(computed - value))
        row[f"radial_history_difference_vs_{source}"] = worst / scale
        lines.append(
            f"- u_r(A) history vs {spec['citation']}: largest difference "
            f"{100 * worst / scale:.1f}% of the published peak (information)"
        )
    return lines


def run_sweep(args: argparse.Namespace, reference: dict) -> bool:
    """Run the mesh, time-step or published-discretisation study."""
    study = args.study
    if study == "literature":
        # The published finite element results use first-order implicit
        # Euler with dt = 1e-4 s, which damps the pulse appreciably.
        args.time_scheme = reference["literature"]["time_scheme"]
    end_time = float(reference["quick"]["end_time_s"] if args.quick
                     else reference["end_time_s"])
    suffix = "_quick" if args.quick else ""
    if args.delta_t:
        suffix = f"_dt{args.delta_t:g}" + suffix
    if args.probe_z_shift:
        suffix = f"_pz{args.probe_z_shift:g}" + suffix
    if args.tight_tolerances:
        suffix = "_tight" + suffix
    rows, history_files = [], []
    for level, factor, delta_t in sweep_levels(args, reference, study):
        cores = study_cores(args.cores, study, level, reference)
        name = f"{args.coupling}_{args.time_scheme}_{study}{level}{suffix}"
        case = run_or_reuse(
            name, args.coupling, factor, delta_t, end_time, args, reference,
            cores, f"{args.coupling} {study} level {level}",
        )
        row = {
            "study": study, "coupling": args.coupling,
            "time_scheme": args.time_scheme, "level": level,
            "refinement": factor,
            "fluid_cells": cell_count(block_mesh_dict(case, "fluid")),
            "solid_cells": cell_count(block_mesh_dict(case, "solid")),
            "delta_t_s": delta_t, "end_time_s": end_time,
            "cores": solver_ranks(case),
            "clock_time_s": clock_time(case),
        }
        row.update(extract(case, reference, end_time))
        row.update(coupling_summary(case, args.coupling, end_time))
        history_files.append(write_histories(name, row["_histories"]))
        rows.append(row)

    acceptance = reference["acceptance"]
    estimates = wave_speed_estimates(reference)
    wave_reference = estimates[reference["wave_speed_reference"]]
    for row in rows:
        for name, value in estimates.items():
            row[f"{name}_wave_speed_m_s"] = value
        row["wave_speed_pressure_vs_estimate"] = relative_difference(
            row["wave_speed_pressure_m_s"], wave_reference
        )
    tolerances = acceptance.get(f"{study}_self_convergence", {})
    for previous, row in zip(rows, rows[1:]):
        for quantity in tolerances:
            row[f"{quantity}_change"] = relative_difference(
                row[quantity], previous[quantity]
            )
    orders = {
        quantity: observed_order([row[quantity] for row in rows])
        for quantity in tolerances
    }
    for row in rows:
        for quantity, order in orders.items():
            row[f"{quantity}_observed_order"] = order

    checks: list[tuple[str, bool]] = [(
        f"all {len(rows)} run(s) completed, and every time step met the "
        f"{args.coupling} convergence criteria",
        True,
    )]
    if not args.quick:
        finest = rows[-1]
        if tolerances and len(rows) < 2:
            checks.append(("at least two levels are required", False))
        for quantity, tolerance in tolerances.items():
            if len(rows) < 2:
                break
            change = finest[f"{quantity}_change"]
            checks.append((
                f"{quantity}: change between the two finest levels "
                f"{100 * change:.2f}% <= {100 * tolerance:g}%",
                within(change, tolerance),
            ))
            if len(rows) >= 3:
                # A change below the resolution floor counts as converged:
                # below it, the order of successive changes is noise.
                floor = float(acceptance["change_resolution_floor"])
                changes = [row[f"{quantity}_change"] for row in rows[1:]]
                checks.append((
                    f"{quantity}: successive changes decrease or are below "
                    f"{100 * floor:g}% "
                    f"({', '.join(f'{100 * c:.3g}%' for c in changes)}; "
                    f"observed order {format_value(orders[quantity], 3)})",
                    all(math.isfinite(b) and (b < a or b <= floor)
                        for a, b in zip(changes, changes[1:])),
                ))
        tolerance = acceptance["wave_speed_vs_estimate"].get(study)
        if tolerance is not None:
            difference = finest["wave_speed_pressure_vs_estimate"]
            checks.append((
                "pressure-front speed "
                f"{format_value(finest['wave_speed_pressure_m_s'], 4)} m/s vs "
                f"the {reference['wave_speed_reference']} estimate "
                f"{wave_reference:.4g} m/s (approximate): "
                f"{100 * difference:.1f}% <= {100 * tolerance:g}%",
                within(difference, tolerance),
            ))
        checks += published_checks(rows, reference, study)
    passed = all(ok for _, ok in checks)
    information = history_differences(rows[-1], reference, study)

    stem = f"{study}_{args.coupling}_{args.time_scheme}{suffix}"
    csv_path = OUTPUT_ROOT / f"{stem}.csv"
    write_csv(csv_path, rows)
    title = {
        "mesh": "Mesh study (mesh and time step refined together",
        "timestep": "Time-step study (tutorial mesh",
        "literature": "Published-discretisation study (tutorial mesh",
    }[study]
    if args.delta_t:
        title = f"Mesh study (fixed Δt = {args.delta_t:g} s"
    lines = [
        f"## {title}; {args.coupling}, {args.time_scheme}"
        f"{', quick smoke test' if args.quick else ''})",
        "",
        f"- Result: {'PASS' if passed else 'FAIL'}",
        "- Analytical wave-speed estimates (m/s): " + ", ".join(
            f"{name} {value:.3f}" for name, value in estimates.items()
        ),
        "",
        "| Level | Fluid cells | Solid cells | Δt (s) | Cores "
        "| u_r,max(A) (mm) | t(u_r,max) (ms) | u_z,min(A) (mm) "
        "| u_r,min(A), t > 14 ms (mm) | t_arr(A) (ms) | c_wall (m/s) "
        "| c_p (m/s) | Mean FSI iter. | Clock (s) |",
        "|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in rows:
        lines.append(
            f"| {row['level']} | {row['fluid_cells']} | {row['solid_cells']} "
            f"| {row['delta_t_s']:.4g} | {row['cores']} "
            f"| {1e3 * row['ur_max_A_m']:.5f} | {1e3 * row['t_ur_max_A_s']:.3f} "
            f"| {1e3 * row['uz_min_A_m']:.5f} "
            f"| {format_value(1e3 * row['ur_min_late_A_m'], 4)} "
            f"| {1e3 * row['t_arrival_A_s']:.3f} "
            f"| {row['wave_speed_wall_m_s']:.3f} "
            f"| {row['wave_speed_pressure_m_s']:.3f} "
            f"| {row['mean_fsi_iterations']:.2f} | {row['clock_time_s']:.0f} |"
        )
    lines += ["", "Checks:", ""]
    lines += [f"- {'PASS' if ok else 'FAIL'}: {text}" for text, ok in checks]
    lines += information
    lines += ["", f"- Data: `{csv_path.name}`", ""]
    write_summary(stem, lines)
    create_history_plot(
        stem, history_files, [f"level_{row['level']}" for row in rows],
        reference, study,
    )
    return passed


def create_history_plot(stem: str, history_files: list[Path],
                        titles: list[str], reference: dict, study: str) -> None:
    """Plot the point-A histories with the published histories overlaid."""
    sources = reference["published_histories"].get(study, [])
    axial_sources = [
        source for source in sources
        if "axial_displacement_m" in (
            REFERENCE_DIR / reference["published"][source]["history"]
        ).read_text()
    ]
    run_gnuplot("plotPointAHistory.gnuplot", {
        "stem": stem,
        "files": " ".join(
            str(path.relative_to(VERIFICATION)) for path in history_files
        ),
        "titles": " ".join(titles),
        "references": " ".join(
            str((REFERENCE_DIR / reference["published"][source]["history"])
                .relative_to(VERIFICATION))
            for source in sources
        ),
        "referenceTitles": " ".join(
            reference["published"][source]["label"] for source in sources
        ),
        "axialReferences": " ".join(
            str((REFERENCE_DIR / reference["published"][source]["history"])
                .relative_to(VERIFICATION))
            for source in axial_sources
        ),
        "axialReferenceTitles": " ".join(
            reference["published"][source]["label"] for source in axial_sources
        ),
    })


def run_coupling_study(args: argparse.Namespace, reference: dict) -> bool:
    spec = reference["coupling"]
    end_time = float(reference["quick"]["end_time_s"] if args.quick
                     else reference["end_time_s"])
    delta_t = float(reference["mesh"]["base_delta_t_s"])
    cores = study_cores(args.cores, "coupling", 1, reference)
    couplings = args.couplings or list(spec["default_couplings"])
    if "robin" not in couplings or "iqnils" not in couplings:
        fail("the coupling study needs at least the robin and iqnils couplings")
    suffix = "_quick" if args.quick else ""
    rows, history_files, results = [], [], {}
    for coupling in couplings:
        name = f"{coupling}_{args.time_scheme}_coupling{suffix}"
        case = run_or_reuse(
            name, coupling, 1, delta_t, end_time, args, reference, cores,
            f"{coupling} coupling on the tutorial mesh",
        )
        row = {"study": "coupling", "coupling": coupling,
               "time_scheme": args.time_scheme, "delta_t_s": delta_t,
               "end_time_s": end_time, "cores": solver_ranks(case),
               "clock_time_s": clock_time(case)}
        row.update(extract(case, reference, end_time))
        row.update(coupling_summary(case, coupling, end_time))
        results[coupling] = row
        history_files.append(write_histories(name, row["_histories"]))
        rows.append(row)

    tolerance = float(spec["relative_tolerance"])
    reference_history = results["iqnils"]["_histories"]["radial_A"]
    scale = max(abs(value) for _, value in reference_history)
    checks = [(
        "every time step of every coupling met its convergence criteria "
        f"(Robin: {int(results['robin']['robin_stalled_steps'])} step(s) "
        "accepted as stalled within the stall tolerance)",
        True,
    )]
    for coupling in couplings:
        if coupling == "iqnils":
            continue
        row = results[coupling]
        history = dict(row["_histories"]["radial_A"])
        if any(time not in history for time, _ in reference_history):
            fail(f"{coupling} and iqnils histories have different time levels")
        row["radial_A_history_difference_vs_iqnils"] = max(
            abs(history[time] - value) for time, value in reference_history
        ) / scale
        difference = row["radial_A_history_difference_vs_iqnils"]
        checks.append((
            f"{coupling} vs iqnils: max |u_r(A) difference| / max |u_r(A)| "
            f"{100 * difference:.3f}% <= {100 * tolerance:g}%",
            within(difference, tolerance),
        ))
        for quantity in spec["compared_quantities"]:
            difference = relative_difference(
                row[quantity], results["iqnils"][quantity]
            )
            row[f"{quantity}_vs_iqnils"] = difference
            checks.append((
                f"{coupling} vs iqnils: {quantity} {100 * difference:.3f}% "
                f"<= {100 * tolerance:g}%",
                within(difference, tolerance),
            ))
    passed = all(ok for _, ok in checks)

    stem = f"coupling_{args.time_scheme}{suffix}"
    csv_path = OUTPUT_ROOT / f"{stem}.csv"
    write_csv(csv_path, rows)
    robin = results["robin"]
    lines = [
        f"## Coupling study (tutorial mesh, Δt = {delta_t:g} s, "
        f"{args.time_scheme}{', quick smoke test' if args.quick else ''})",
        "",
        f"- Result: {'PASS' if passed else 'FAIL'}",
        "",
        "| Coupling | u_r,max(A) (mm) | u_z,min(A) (mm) | t_arr(A) (ms) "
        "| c_p (m/s) | Total FSI iter. | Mean | Max | Cores | Clock (s) |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in rows:
        lines.append(
            f"| {row['coupling']} | {1e3 * row['ur_max_A_m']:.5f} "
            f"| {1e3 * row['uz_min_A_m']:.5f} "
            f"| {1e3 * row['t_arrival_A_s']:.3f} "
            f"| {row['wave_speed_pressure_m_s']:.3f} "
            f"| {int(row['total_fsi_iterations'])} "
            f"| {row['mean_fsi_iterations']:.2f} "
            f"| {int(row['maximum_fsi_iterations'])} "
            f"| {row['cores']} | {row['clock_time_s']:.0f} |"
        )
    lines += ["", "Checks:", ""]
    lines += [f"- {'PASS' if ok else 'FAIL'}: {text}" for text, ok in checks]
    lines += [
        "- Worst final Robin residuals: displacement "
        f"{robin['maximum_final_displacement_residual']:.3g}, pressure "
        f"{robin['maximum_final_pressure_residual']:.3g}, leakage flux "
        f"{robin['maximum_final_flux_residual']:.3g}",
        "", f"- Data: `{csv_path.name}`", "",
    ]
    write_summary(stem, lines)
    create_history_plot(stem, history_files, couplings, reference, "coupling")
    return passed


def write_summary(stem: str, lines: list[str]) -> None:
    """Store this study's section and rebuild the combined summary."""
    sections = OUTPUT_ROOT / "sections"
    sections.mkdir(parents=True, exist_ok=True)
    (sections / f"{stem}.md").write_text("\n".join(lines) + "\n")
    text = "# 3dTube verification summary\n\n" + "\n".join(
        path.read_text() for path in sorted(sections.glob("*.md"))
    )
    (OUTPUT_ROOT / "verification_summary.md").write_text(text)
    print("\n".join(lines))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--study", choices=("mesh", "timestep", "literature", "coupling"),
        default="mesh",
        help="mesh (default): mesh and time step refined together; "
             "timestep: time step only, tutorial mesh; literature: the "
             "published finite element time discretisation; coupling: "
             "couplings compared on the tutorial mesh",
    )
    parser.add_argument(
        "--coupling", choices=tuple(COUPLINGS), default="robin",
        help="coupling of the mesh and time-step studies (default: robin)",
    )
    parser.add_argument(
        "--couplings",
        help="comma-separated couplings of the coupling study "
             "(default: robin,iqnils)",
    )
    parser.add_argument(
        "--levels",
        help="comma-separated levels; mesh level 1 is the tutorial mesh",
    )
    parser.add_argument(
        "--delta-t", type=float,
        help="mesh study only: use this time step on every level instead of "
             "halving it with the mesh (separates the spatial error)",
    )
    parser.add_argument(
        "--tight-tolerances", action="store_true",
        help="diagnostic: tighten the fluid, solid and FSI tolerances to "
             "bound the iterative error",
    )
    parser.add_argument(
        "--probe-z-shift", type=float, default=0.0,
        help="shift the axis pressure probes axially by this distance (m); "
             "the stations lie on cell faces of every mesh level, where the "
             "containing cell is ambiguous, so a shift of half a level-3 "
             "cell (7.8125e-5) samples each probe inside one cell",
    )
    parser.add_argument(
        "--time-scheme", choices=("backward", "Euler"), default="backward",
        help="time scheme of the fluid and solid (default: backward; the "
             "tutorial uses Euler)",
    )
    parser.add_argument(
        "--solid-preconditioner", choices=("hypre", "lu"), default="hypre",
        help="solid PETSc preconditioner (default: hypre; the tutorial uses lu)",
    )
    parser.add_argument(
        "--cores", default="auto",
        help="MPI ranks per case: positive integer or auto (default)",
    )
    parser.add_argument(
        "--quick", action="store_true",
        help="smoke test: first level only, to a short end time, "
             "no accuracy checks",
    )
    parser.add_argument(
        "--reuse", action="store_true",
        help="re-evaluate completed runs under verification/work",
    )
    args = parser.parse_args()
    if args.cores != "auto" and (not args.cores.isdecimal() or int(args.cores) < 1):
        parser.error("--cores must be a positive integer or auto")
    if args.levels:
        args.levels = [int(value) for value in args.levels.split(",")]
        if any(level < 1 for level in args.levels) or any(
            b <= a for a, b in zip(args.levels, args.levels[1:])
        ):
            parser.error("--levels must be strictly increasing positive integers")
    if args.delta_t is not None and (
        args.study != "mesh" or not args.delta_t > 0.0
    ):
        parser.error("--delta-t takes a positive time step and applies only "
                     "to the mesh study")
    if args.couplings:
        args.couplings = [value.strip() for value in args.couplings.split(",")]
        invalid = [name for name in args.couplings if name not in COUPLINGS]
        if invalid or len(set(args.couplings)) != len(args.couplings):
            parser.error(
                f"--couplings must be unique names from {', '.join(COUPLINGS)}"
            )

    for executable in ("blockMesh", "solids4Foam"):
        if shutil.which(executable) is None:
            fail(f"Required executable '{executable}' is unavailable. "
                 "Source OpenFOAM and build solids4foam first.")
    reference = json.loads(REFERENCE_FILE.read_text())
    WORK_ROOT.mkdir(parents=True, exist_ok=True)
    OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)
    if args.study == "coupling":
        passed = run_coupling_study(args, reference)
    else:
        passed = run_sweep(args, reference)
    print(f"Summary: {OUTPUT_ROOT / 'verification_summary.md'}")
    return 0 if passed else 1


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (RuntimeError, subprocess.SubprocessError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        sys.exit(2)
