#!/usr/bin/env python3
"""Run the opt-in numerical verification of the hessenthalerFsi Phase I case.

The studies hold the physical model fixed (geometry, inflow, fluid
properties, neo-Hookean solid with one fixed shear modulus, nu = 0.45,
momentum stabilisation 0.01, buoyancy, ramp) and change only numerical
parameters: the solid mesh, the fluid mesh, the time step, the end time and
the coupling tolerance. The measured data are not used here; the comparison
with the experiment is the validation study in ../validation.

solid
    The flap alone under its net buoyancy (the zero-flow problem) on a
    nested solid mesh family S1-S4, with momentum stabilisation 0.01 and
    0.001. Gives the static tip deflection at the fixed shear modulus and the
    shear modulus each discretisation needs to reproduce the measured
    zero-flow deflection, with observed orders and Richardson estimates.

meshes
    Build the fluid meshes F1-F5 and the solid meshes S1-S4 and record their
    sizes and quality, without solving.

run
    One coupled Phase I run, given by --spec, for example
    "F2:S3", "F1:S2:dt=0.001", "F1:S2:tol=1e-6", "F1:S2:T=30",
    "F1:S2:coupling=robin" or "F1:S2:fluidtol=tight". The run is evaluated
    and its quantities of interest are written to
    postProcessing/runs/<name>.json.

analyse
    Collect the evaluated runs (postProcessing/runs/*.json and the solid
    study) into tables of successive differences, observed orders, Richardson
    estimates and numerical uncertainties, and evaluate every run against the
    measurements, kept separate from the verification.

Every run is a copy of the tutorial under verification/work; results are
written to verification/postProcessing.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import re
import shutil
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

SCRIPT = Path(__file__).resolve()
VERIFICATION = SCRIPT.parents[1]
TUTORIAL = VERIFICATION.parent
VALIDATION = TUTORIAL / "validation"
WORK_ROOT = VERIFICATION / "work"
OUTPUT_ROOT = VERIFICATION / "postProcessing"
RUNS_ROOT = OUTPUT_ROOT / "runs"
SETTINGS_FILE = "verification_settings.json"

# Reuse the file, run and sampling helpers of the validation driver
sys.path.insert(0, str(VALIDATION / "scripts"))
import hessenthaler_fsi_validation as vmod  # noqa: E402

fail = vmod.fail
numeric_rows = vmod.numeric_rows
last_per_time = vmod.last_per_time
read_csv = vmod.read_csv
write_csv = vmod.write_csv
replace_entry = vmod.replace_entry
replace_text = vmod.replace_text
# The verification runs keep their own settings file, which vmod.reusable reads
vmod.SETTINGS_FILE = SETTINGS_FILE

# ---------------------------------------------------------------------------
# Fixed model parameters of the verification problem
# ---------------------------------------------------------------------------

# Shear modulus held fixed for every mesh, time step and coupling setting.
# It is the mesh-converged zero-flow calibration of the standard solid (see
# the solid study and README.md); it is NOT recalibrated per mesh.
MU_VERIFICATION = 64200.0
NU = 0.45
STABILISATION = 0.01
DELTA_T = 0.002
END_TIME = 15.0
# Fraction of the run at its end over which quantities are time-averaged,
# and the preceding window of the same length used for the drift test
WINDOW_FRACTION = 0.2

# ---------------------------------------------------------------------------
# Mesh families
# ---------------------------------------------------------------------------

# Nested hexahedral flap meshes (cells in x, y, z) of the 11 x 2 x 65 mm
# flap; S2 is the tutorial mesh, and S2-S4 refine it by exactly 2 in every
# direction (S1 halves S2 with nz = 33 instead of 32.5)
SOLID_LEVELS = {
    "S1": (6, 4, 33),
    "S2": (12, 8, 65),
    "S3": (24, 16, 130),
    "S4": (48, 32, 260),
}
FLAP_MM = (11.0, 2.0, 65.0)

# cfMesh cell sizes (m) of fluid level F1, which is the tutorial's coarse mesh
# (system/meshDict.coarse). Level Fk scales EVERY size, including the
# octree root size maxCellSize, by 2^(-(k-1)/2), so that successive levels
# refine every region consistently by sqrt(2) and F1, F3, F5 refine it by 2.
# Each size is set 1 % above an octree level, as in the tutorial, so that
# rounding never selects the next finer level.
FLUID_F1_SIZES = {
    "maxCellSize": 0.004,
    "boundaryCellSize": 0.00202,
    "interface": 0.000505,
    "flapSweep": 0.00101,
    "merging": 0.00202,
}
FLUID_LEVELS = ("F1", "F2", "F3", "F4", "F5")
FLUID_RATIO = math.sqrt(2.0)

# Points of the velocity probes (m): the measured jet peaks on the two MRI
# planes, the shear layer between the jets, and the wake of the flap tip
PROBES = (
    (0.0, 0.025, 0.010), (0.0, -0.025, 0.010),
    (0.0, 0.025, 0.030), (0.0, -0.025, 0.030),
    (0.0, 0.010, 0.030), (0.0, 0.0, 0.080),
)

# Undeformed (material) z (mm) of the centreline monitors: every 1 mm, which
# are exact mesh nodes of S2-S4 (S1 uses its nearest nodes, 1.97 mm apart)
CENTRELINE_Z = [float(j) for j in range(1, 66)]
# Material stations at which displacements are compared between levels;
# exact nodes of S2-S4, interpolated linearly between the S1 nodes
COMPARE_Z = [5.0 * j for j in range(1, 14)]
# Time between the fixed-voxel velocity samples of the averaging window (s)
VOXEL_SAMPLE_INTERVAL = 0.1


def fluid_scale(level: str) -> float:
    k = int(level[1:])
    return FLUID_RATIO ** (-(k - 1))


def fluid_sizes(level: str) -> dict[str, float]:
    return {key: value * fluid_scale(level)
            for key, value in FLUID_F1_SIZES.items()}


def write_mesh_dict(case: Path, level: str) -> Path:
    """Write system/meshDict.<level> from the tutorial's coarse meshDict."""
    sizes = fluid_sizes(level)
    text = (TUTORIAL / "system" / "meshDict.coarse").read_text()
    replacements = (
        (r"^(maxCellSize\s+)[^;]+;", sizes["maxCellSize"]),
        (r"^(boundaryCellSize\s+)[^;]+;", sizes["boundaryCellSize"]),
        (r"(interface\s*\{\s*cellSize\s+)[^;]+;", sizes["interface"]),
        (r"(flapSweep\s*\{[^}]*?cellSize\s+)[^;]+;", sizes["flapSweep"]),
        (r"(merging\s*\{[^}]*?cellSize\s+)[^;]+;", sizes["merging"]),
    )
    for pattern, value in replacements:
        text, found = re.subn(pattern, rf"\g<1>{value:.8g};", text,
                              flags=re.MULTILINE | re.DOTALL)
        if found != 1:
            fail(f"Could not set '{pattern}' in the meshDict for {level}")
    text = re.sub(r"^// cfMesh cartesianMesh settings for the coarse.*$",
                  f"// Verification fluid level {level}: every cell size of "
                  f"meshDict.coarse scaled by {fluid_scale(level):.6g}",
                  text, count=1, flags=re.MULTILINE)
    path = case / "system" / f"meshDict.{level}"
    path.write_text(text)
    return path


def centreline_nodes(nz: int) -> list[float]:
    """Undeformed z (mm) of the mesh nodes nearest to every 1 mm.

    solidPointDisplacement reports the displacement of the closest mesh
    point, so the monitors are placed exactly on nodes, and their true
    undeformed positions are used.
    """
    dz = FLAP_MM[2] / nz
    nodes = sorted({round(round(z / dz) * dz, 9) for z in CENTRELINE_Z})
    return [z for z in nodes if z > 0]


def write_centreline_monitors(case: Path, nz: int) -> list[float]:
    nodes = centreline_nodes(nz)
    blocks = []
    for j, z in enumerate(nodes):
        blocks.append(f"""flapCentrelineN{j:02d}
{{
    type            solidPointDisplacement;
    region          solid;
    point           (0 0 {z * 1e-3:.10g});
}}
""")
    path = case / "system" / "flapCentreline"
    header = path.read_text()
    header = header[:header.index("// Solid points")]
    path.write_text(header + "// Solid points on the flap centreline at "
                    "mesh nodes, written by the verification driver\n"
                    + "\n".join(blocks))
    return nodes


def voxel_history(avg_start: float, end_time: float, delta_t: float) -> str:
    """Fixed-point velocity samples on the MRI voxels during the window.

    The fluid mesh moves, so the cell-based fieldAverage mean is not the time
    average at the fixed MRI voxels; these samples are, after averaging.
    """
    velocity = read_csv(vmod.REFERENCE_DIR / "phaseI_velocity.csv")
    points = "\n".join(f"                ({p[0]:.8g} {p[1]:.8g} {p[2]:.8g})"
                       for _, p in vmod.voxel_points(velocity))
    # At least ten samples in the window
    interval = min(VOXEL_SAMPLE_INTERVAL, (end_time - avg_start) / 10)
    every = max(1, round(interval / delta_t))
    return f"""FoamFile
{{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      voxelHistory;
}}

// Velocity at the sub-sampling points of every MRI voxel, every
// {every*delta_t:g} s of the averaging window (verification driver)
voxelHistory
{{
    type            sets;
    libs            (sampling);
    region          fluid;
    interpolationScheme cellPoint;
    setFormat       raw;
    fields          (U);
    timeStart       {avg_start:.10g};
    writeControl    timeStep;
    writeInterval   {every};
    sets
    {{
        voxels
        {{
            type    cloud;
            axis    xyz;
            points
            (
{points}
            );
        }}
    }}
}}
"""


def verification_monitors(avg_start: float) -> str:
    probes = "\n".join(f"            ({x:g} {y:g} {z:g})"
                       for x, y, z in PROBES)
    return f"""FoamFile
{{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      verificationMonitors;
}}

// Added by the verification driver (included in system/controlDict)

// Fluid force on the flap
flapForce
{{
    type            forces;
    libs            (forces);
    region          fluid;
    patches         (interface);
    rho             rhoInf;
    rhoInf          1163.3;
    CofR            (0 0 0);
    writeControl    timeStep;
    writeInterval   1;
}}

// Mean pressure at each inlet (the outlet pressure is zero)
upperInletPressure
{{
    type            surfaceFieldValue;
    libs            (fieldFunctionObjects);
    region          fluid;
    regionType      patch;
    name            upperInlet;
    operation       areaAverage;
    fields          (p);
    writeFields     false;
}}
lowerInletPressure
{{
    type            surfaceFieldValue;
    libs            (fieldFunctionObjects);
    region          fluid;
    regionType      patch;
    name            lowerInlet;
    operation       areaAverage;
    fields          (p);
    writeFields     false;
}}

// Velocity history at the jet peaks, the shear layer and the wake
velocityProbes
{{
    type            probes;
    libs            (sampling);
    region          fluid;
    fields          (U);
    probeLocations
    (
{probes}
    );
}}

// Time average of the fluid fields over the final window
fluidAverage
{{
    type            fieldAverage;
    libs            (fieldFunctionObjects);
    region          fluid;
    timeStart       {avg_start:.10g};
    // Write the means with the fields only
    writeControl    writeTime;
    restartOnRestart false;
    restartOnOutput false;
    fields
    (
        U
        {{
            mean        on;
            prime2Mean  off;
            base        time;
        }}
    );
}}
"""


# ---------------------------------------------------------------------------
# Small helpers
# ---------------------------------------------------------------------------

def environment(threads: int = 1) -> dict[str, str]:
    env = vmod.run_environment()
    env["OMP_NUM_THREADS"] = str(threads)
    return env


def run(command: list[str], case: Path, log_name: str,
        threads: int = 1) -> float:
    """Run a command, logging to case/log_name; return its wall time (s)."""
    start = time.time()
    with (case / log_name).open("w") as handle:
        result = subprocess.run(command, cwd=case, stdout=handle,
                                stderr=subprocess.STDOUT, text=True,
                                env=environment(threads))
    if result.returncode:
        fail(f"'{' '.join(command)}' failed in {case}; see {case / log_name}")
    return time.time() - start


def mpirun(cores: int, command: list[str]) -> list[str]:
    if cores == 1:
        return command
    return ["mpirun", "-np", str(cores)] + command + ["-parallel"]


def use_hypre(fv_solution: Path) -> None:
    """Replace the solid's LU preconditioner by hypre BoomerAMG."""
    text = fv_solution.read_text()
    hypre = ("pc_type hypre;\n"
             "            pc_hypre_type boomeramg;\n"
             "            pc_hypre_boomeramg_max_iter \"1\";\n"
             "            pc_hypre_boomeramg_strong_threshold \"0.7\";\n"
             "            pc_hypre_boomeramg_grid_sweeps_up \"1\";\n"
             "            pc_hypre_boomeramg_grid_sweeps_down \"1\";\n"
             "            pc_hypre_boomeramg_agg_nl \"1\";\n"
             "            pc_hypre_boomeramg_agg_num_paths \"1\";\n"
             "            pc_hypre_boomeramg_max_levels \"25\";\n"
             "            pc_hypre_boomeramg_coarsen_type HMIS;\n"
             "            pc_hypre_boomeramg_interp_type ext+i;\n"
             "            pc_hypre_boomeramg_P_max \"1\";\n"
             "            pc_hypre_boomeramg_truncfactor \"0.3\";")
    text, found = re.subn(
        r"pc_type\s+(bjacobi\s*;\s*sub_pc_type\s+lu|lu\s*;\s*"
        r"pc_factor_mat_solver_type\s+mumps)\s*;", hypre, text)
    if found != 1:
        fail(f"Expected one LU preconditioner in {fv_solution}")
    fv_solution.write_text(text)


def telescope_ranks(spec: dict, cores: int) -> int:
    """Ranks of the solid preconditioner's sub-communicator.

    The flap is small next to the fluid, so its LU preconditioner is applied
    on a few ranks (PCTELESCOPE): an exact LU distributed over hundreds of
    ranks is dominated by communication latency.
    """
    wanted = int(spec["pcranks"]) or (
        8 if math.prod(SOLID_LEVELS[spec["solid"]]) < 20000 else 16)
    wanted = max(1, min(wanted, cores))
    while cores % wanted:
        wanted -= 1
    return wanted


def use_telescope(fv_solution: Path, cores: int, ranks: int) -> None:
    """Apply the solid's LU (MUMPS) preconditioner on 'ranks' ranks."""
    text = fv_solution.read_text()
    telescope = ("pc_type telescope;\n"
                 f"            pc_telescope_reduction_factor \"{cores // ranks}\";\n"
                 "            telescope_ksp_type preonly;\n"
                 "            telescope_pc_type lu;\n"
                 "            telescope_pc_factor_mat_solver_type mumps;")
    text, found = re.subn(r"pc_type\s+bjacobi\s*;\s*sub_pc_type\s+lu\s*;",
                          telescope, text)
    if found != 1:
        fail(f"Expected one block-Jacobi LU preconditioner in {fv_solution}")
    fv_solution.write_text(text)


def set_material(mechanical: Path, mu: float, nu: float) -> None:
    bulk = 2.0 * mu * (1.0 + nu) / (3.0 * (1.0 - 2.0 * nu))
    text = mechanical.read_text()
    text, found = re.subn(r"(mu\s+mu\s+\[[^]]*\]\s+)[^;]+;",
                          rf"\g<1>{mu:.10g};", text)
    text, found2 = re.subn(r"(K\s+K\s+\[[^]]*\]\s+)[^;]+;",
                           rf"\g<1>{bulk:.10g};", text)
    if found != 1 or found2 != 1:
        fail(f"Could not set mu and K in {mechanical}")
    mechanical.write_text(text)


def window_stats(rows: list[list[float]], column: int, start: float,
                 end: float) -> dict[str, float]:
    """Statistics of one column over the half-open window (start, end].

    The histories have a constant time step, so the plain mean is the time
    average.
    """
    values = [r[column] for r in rows if start + 1e-9 < r[0] <= end + 1e-9]
    if not values:
        fail(f"No samples in the window ({start:g}, {end:g}]")
    mean = sum(values) / len(values)
    return {"mean": mean, "min": min(values), "max": max(values),
            "std": math.sqrt(sum((v - mean) ** 2 for v in values)
                             / len(values)), "samples": len(values)}


def history_rows(case: Path, name: str, filename: str, end_time: float,
                 delta_t: float, columns: int) -> list[list[float]]:
    """A complete postProcessing history of function object 'name'.

    Region function objects write to postProcessing/<region>/<name>/<t0>/ or
    postProcessing/<name>/<region>/<t0>/; restart segments are merged, and
    the history must cover every time step to the end time.
    """
    candidates = [path for path in case.glob(f"postProcessing/**/{filename}")
                  if name in path.parts]
    candidates.sort(key=lambda path: float(path.parent.name))
    if not candidates:
        fail(f"No '{name}' history ({filename}) in {case}")
    rows: list[list[float]] = []
    for path in candidates:
        start = float(path.parent.name)
        rows = [row for row in rows if row[0] <= start + 1e-12]
        rows.extend(numeric_rows(path))
    rows = last_per_time(rows)
    vmod.require_history(candidates[-1], rows, end_time, delta_t, columns)
    return rows


def settling_estimate(rows: list[list[float]], column: int, end_time: float,
                      scale: float = 1.0) -> dict:
    """Time-average of a QoI over the last third of the run, with its
    residual time (statistical and drift) uncertainty.

    After the ramp the flap does not settle to a strict steady state: it
    wanders at low frequency (periods of seconds) by a few hundredths of a
    millimetre, and on some meshes still drifts slowly. The QoI is the mean
    over (2T/3, T]; its time uncertainty is the larger of
      - the difference between the means of the two halves of that window
        (a drift that has not died out), and
      - twice the standard error of the means of its 1-s batches (the
        low-frequency wander; batches of 1 s are not independent, so this
        is indicative).
    """
    start = end_time * 2.0 / 3.0
    data = [(r[0], scale * r[column]) for r in rows
            if start + 1e-9 < r[0] <= end_time + 1e-9]
    mid = 0.5 * (start + end_time)
    first = [y for t, y in data if t <= mid + 1e-9]
    second = [y for t, y in data if t > mid + 1e-9]
    mean = sum(y for _, y in data) / len(data)
    batches = []
    edge = start
    while edge + 1.0 <= end_time + 1e-9:
        batch = [y for t, y in data if edge + 1e-9 < t <= edge + 1.0 + 1e-9]
        if batch:
            batches.append(sum(batch) / len(batch))
        edge += 1.0
    n = len(batches)
    se = (math.sqrt(sum((b - sum(batches) / n) ** 2 for b in batches)
                    / (n - 1) / n) if n > 1 else math.nan)
    drift = abs(sum(second) / len(second) - sum(first) / len(first))
    return {"mean": mean, "window_start_s": start, "half_difference": drift,
            "batch_standard_error": se, "batches": n,
            "uncertainty": max(drift, 2.0 * se if math.isfinite(se) else 0.0)}


def checkmesh_metrics(log: Path) -> dict:
    text = log.read_text(errors="replace")
    metrics = {}
    number = r"([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)"
    patterns = {
        "cells": r"^\s*cells:\s+(\d+)",
        "faces": r"^\s*faces:\s+(\d+)",
        "points": r"^\s*points:\s+(\d+)",
        "max_non_orthogonality_deg": r"Mesh non-orthogonality Max:\s*" + number,
        "mean_non_orthogonality_deg": r"Mesh non-orthogonality Max:\s*[0-9.eE+-]+\s+average:\s*" + number,
        "max_skewness": r"Max skewness\s*=\s*" + number,
        "max_aspect_ratio": r"Max aspect ratio\s*=\s*" + number,
        "total_volume_m3": r"Total volume\s*=\s*" + number,
        "min_volume_m3": r"Min volume\s*=\s*" + number,
    }
    for key, pattern in patterns.items():
        match = re.search(pattern, text, re.MULTILINE)
        if match:
            metrics[key] = float(match.group(1))
    metrics["checkmesh_ok"] = "Mesh OK." in text
    return metrics


def patch_faces(boundary: Path, patch: str) -> int:
    text = boundary.read_text()
    match = re.search(rf"\b{re.escape(patch)}\s*\{{[^}}]*?nFaces\s+(\d+);",
                      text, re.DOTALL)
    return int(match.group(1)) if match else -1


# ---------------------------------------------------------------------------
# Solid study: the zero-flow problem on the nested solid family
# ---------------------------------------------------------------------------

# Buoyancy load factors of the plateaus, relative to the Phase I value. The
# static deflection of a neo-Hookean solid with a fixed Poisson's ratio (and a
# stabilisation that scales with the elastic moduli) depends only on
# rho*g/mu, so load f at MU_REFERENCE is the static state at MU_REFERENCE/f:
# the plateaus span 48.8-66.3 kPa.
SOLID_LOADS = (0.92, 1.0, 1.1, 1.25)
# Ranks of the parallel solid-study runs (the recorded campaign)
SOLID_STUDY_CORES = {"S3": 32, "S4": 128}
SOLID_PLATEAU = 2.0
# Largest tip movement (mm) over the last 0.25 s of a settled plateau
# (0.02 mm is about 0.06 kPa in the calibrated modulus)
SETTLED_DRIFT_MM = 0.02
MU_REFERENCE = 61000.0


def quadratic_through(points: list[tuple[float, float]], x: float) -> float:
    """Lagrange interpolation through the given (x, y) points."""
    total = 0.0
    for i, (xi, yi) in enumerate(points):
        term = yi
        for j, (xj, _) in enumerate(points):
            if j != i:
                term *= (x - xj) / (xi - xj)
        total += term
    return total


def solve_monotone(func, lo: float, hi: float, target: float) -> float:
    """Bisection for func(x) = target on [lo, hi]."""
    flo, fhi = func(lo) - target, func(hi) - target
    if flo * fhi > 0:
        return math.nan
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        fmid = func(mid) - target
        if flo * fmid <= 0:
            hi, fhi = mid, fmid
        else:
            lo, flo = mid, fmid
    return 0.5 * (lo + hi)


def solid_interpolation(points: list[tuple[float, float]], target_tip: float,
                        mu_fixed: float) -> dict:
    """Calibrated mu and the tip at mu_fixed from (load, tip) plateaus.

    The tip is interpolated in the load factor f (proportional to 1/mu) by a
    quadratic through the three plateaus nearest the target; the linear and
    cubic interpolants give the interpolation uncertainty.
    """
    points = sorted(points)

    def nearest(x: float, n: int) -> list[tuple[float, float]]:
        return sorted(sorted(points, key=lambda p: abs(p[0] - x))[:n])

    result = {}
    f_fixed = MU_REFERENCE / mu_fixed
    for n, label in ((2, "linear"), (3, "quadratic"), (4, "cubic")):
        if len(points) < n:
            continue
        # The load factor at the target, then the tip at the fixed modulus
        guess = min(points, key=lambda p: abs(p[1] - target_tip))[0]
        local = nearest(guess, n)
        f_cal = solve_monotone(lambda f: quadratic_through(local, f),
                               points[0][0] - 0.3, points[-1][0] + 0.3,
                               target_tip)
        local = nearest(f_cal if math.isfinite(f_cal) else guess, n)
        f_cal = solve_monotone(lambda f: quadratic_through(local, f),
                               points[0][0] - 0.3, points[-1][0] + 0.3,
                               target_tip)
        result[f"mu_cal_{label}_Pa"] = (MU_REFERENCE / f_cal
                                        if math.isfinite(f_cal) else math.nan)
        result[f"tip_at_mu_fixed_{label}_mm"] = quadratic_through(
            nearest(f_fixed, n), f_fixed)
        result[f"cal_extrapolated_{label}"] = not (
            points[0][0] <= f_cal <= points[-1][0]) if math.isfinite(f_cal) \
            else True
    result["mu_cal_Pa"] = result["mu_cal_quadratic_Pa"]
    result["tip_at_mu_fixed_mm"] = result["tip_at_mu_fixed_quadratic_mm"]
    result["mu_cal_interp_uncertainty_Pa"] = max(
        abs(result.get(f"mu_cal_{k}_Pa", math.nan) - result["mu_cal_Pa"])
        for k in ("linear", "cubic"))
    result["tip_interp_uncertainty_mm"] = max(
        abs(result.get(f"tip_at_mu_fixed_{k}_mm", math.nan)
            - result["tip_at_mu_fixed_mm"]) for k in ("linear", "cubic"))
    return result


def solid_run_name(level: str, sf: float, cores: int) -> str:
    return f"solid_{level}_sf{sf:g}_np{cores}"


def run_solid(args: argparse.Namespace) -> bool:
    reference = json.loads(vmod.REFERENCE_FILE.read_text())
    target = float(reference["targets"]["zero_flow_tip_y_mm"])
    levels = args.levels or list(SOLID_LEVELS)
    stabilisations = [float(s) for s in args.sf.split(",")]
    # The validation driver builds the solid-only case; point it at this
    # study's work directory, loads and reference modulus
    vmod.WORK_ROOT = WORK_ROOT
    vmod.LOAD_LEVELS = SOLID_LOADS
    vmod.PLATEAU_TIME = SOLID_PLATEAU
    vmod.MU_REFERENCE = MU_REFERENCE
    if args.quick:
        vmod.LOAD_LEVELS = (1.0,)
        vmod.PLATEAU_TIME = 0.25
    loads = vmod.LOAD_LEVELS
    end_time = len(loads) * vmod.PLATEAU_TIME

    specs = []
    for level in levels:
        for sf in stabilisations:
            cells = math.prod(SOLID_LEVELS[level])
            cores = 1 if cells < 20000 else (
                args.cores_solid or SOLID_STUDY_CORES[level])
            specs.append({"level": level, "sf": sf, "cores": cores})

    def settings_of(spec: dict) -> dict:
        return {"study": "solid", **spec, "loads": list(loads),
                "plateau_time": vmod.PLATEAU_TIME,
                "mu_reference": MU_REFERENCE, "nu": NU,
                "delta_t": vmod.CAL_DELTA_T, "damping": vmod.CAL_DAMPING,
                "tutorial_inputs": vmod.tutorial_fingerprint(),
                "build": vmod.build_fingerprint()}

    def execute(spec: dict) -> None:
        name = solid_run_name(spec["level"], spec["sf"], spec["cores"])
        settings = settings_of(spec)
        if vmod.reusable(WORK_ROOT / name, settings, args, name):
            return
        vspec = {"solid": "standard", "model": "tutorial",
                 "mesh": SOLID_LEVELS[spec["level"]], "nu": NU,
                 "sf": spec["sf"], "cores": spec["cores"]}
        case = vmod.solid_case(name, vspec, end_time)
        if spec["cores"] > 1:
            # An exact LU of the finest flap is needlessly expensive
            use_hypre(case / "system" / "fvSolution")
        print(f"Running {name}", flush=True)
        vmod.run_solid_case(case, spec["cores"])
        vmod.check_solver_log(case, name)
        (case / SETTINGS_FILE).write_text(json.dumps(settings, indent=2)
                                          + "\n")

    if not args.evaluate_only:
        with ThreadPoolExecutor(max_workers=max(1, args.jobs)) as pool:
            for future in [pool.submit(execute, spec) for spec in specs]:
                future.result()

    rows = []
    passed = True
    for spec in specs:
        name = solid_run_name(spec["level"], spec["sf"], spec["cores"])
        case = WORK_ROOT / name
        history = vmod.tip_history(case, end_time, vmod.CAL_DELTA_T)
        if args.quick:
            print(f"  {name}: tip y = {1e3 * history[-1][2]:.4f} mm")
            continue
        points, drifts = [], []
        for k, load in enumerate(loads):
            end = (k + 1) * vmod.PLATEAU_TIME
            tail = [r for r in history if end - 0.25 <= r[0] <= end + 1e-9]
            drifts.append(1e3 * (max(r[2] for r in tail)
                                 - min(r[2] for r in tail)))
            points.append((load, 1e3 * tail[-1][2]))
        nx, ny, nz = SOLID_LEVELS[spec["level"]]
        row = {"level": spec["level"], "nx": nx, "ny": ny, "nz": nz,
               "cells": nx * ny * nz,
               "hx_mm": FLAP_MM[0] / nx, "hy_mm": FLAP_MM[1] / ny,
               "hz_mm": FLAP_MM[2] / nz,
               "h_eff_mm": (math.prod(FLAP_MM) / (nx * ny * nz)) ** (1 / 3),
               "stabilisation": spec["sf"], "ranks": spec["cores"],
               "max_plateau_drift_mm": max(drifts),
               "clock_time_s": vmod.clock_time(case)}
        for load, tip in points:
            row[f"tip_mm_load{load:g}"] = tip
            row[f"tip_mm_mu{MU_REFERENCE / load / 1e3:.2f}kPa"] = tip
        row.update(solid_interpolation(points, target, MU_VERIFICATION))
        # Accept a run only if every plateau has settled well below the
        # calibration resolution and the target lies within the plateaus
        row["accepted"] = (row["max_plateau_drift_mm"] < SETTLED_DRIFT_MM
                           and not row["cal_extrapolated_quadratic"])
        passed &= row["accepted"]
        rows.append(row)
        print(f"  {name}: mu_cal = {row['mu_cal_Pa']:.6g} Pa, tip at "
              f"{MU_VERIFICATION / 1e3:g} kPa = "
              f"{row['tip_at_mu_fixed_mm']:.4f} mm, plateau drift "
              f"{row['max_plateau_drift_mm']:.2e} mm", flush=True)
    if args.quick:
        return True
    OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)
    existing = []
    path = OUTPUT_ROOT / "solid_study.json"
    if path.is_file():
        keys = {(r["level"], r["stabilisation"]) for r in rows}
        existing = [r for r in json.loads(path.read_text())
                    if (r["level"], r["stabilisation"]) not in keys]
    rows = sorted(existing + rows,
                  key=lambda r: (r["stabilisation"], r["cells"]))
    path.write_text(json.dumps(rows, indent=2) + "\n")
    write_csv(OUTPUT_ROOT / "solid_study.csv", fmt_rows(rows))
    return passed


def fmt_rows(rows: list[dict]) -> list[dict]:
    keys: list[str] = []
    for row in rows:
        keys += [k for k in row if k not in keys]
    return [{k: (f"{r[k]:.7g}" if isinstance(r.get(k), float) else
                 r.get(k, "")) for k in keys} for r in rows]


# ---------------------------------------------------------------------------
# Coupled Phase I runs
# ---------------------------------------------------------------------------

DEFAULT_SPEC = {
    "fluid": "F1", "solid": "S2", "mu": MU_VERIFICATION, "sf": STABILISATION,
    "dt": DELTA_T, "T": END_TIME, "tol": 1e-4, "coupling": "iqnils",
    "fluidtol": "default", "pimple": 1.0, "pc": "auto",
    "lag": -2.0, "pcranks": 0.0, "solidranks": 0.0,
}


def parse_spec(text: str) -> dict:
    spec = dict(DEFAULT_SPEC)
    parts = text.split(":")
    for part in parts:
        if "=" not in part:
            if part in FLUID_LEVELS or part in ("coarse", "medium"):
                spec["fluid"] = part
            elif part in SOLID_LEVELS:
                spec["solid"] = part
            else:
                fail(f"Unknown spec item '{part}'")
            continue
        key, value = part.split("=", 1)
        if key not in DEFAULT_SPEC:
            fail(f"Unknown spec key '{key}'")
        spec[key] = (value if isinstance(DEFAULT_SPEC[key], str)
                     else float(value))
    if spec["coupling"] not in ("iqnils", "robin"):
        fail("coupling must be iqnils or robin")
    if spec["pc"] not in ("auto", "mumps", "hypre", "bjacobi", "telescope"):
        fail("pc must be auto, mumps, hypre, bjacobi or telescope")
    if spec["fluidtol"] not in ("default", "tight"):
        fail("fluidtol must be default or tight")
    return spec


def spec_name(spec: dict) -> str:
    name = (f"{spec['fluid']}_{spec['solid']}_mu{spec['mu'] / 1e3:g}"
            f"_dt{spec['dt']:g}_T{spec['T']:g}")
    if spec["sf"] != STABILISATION:
        name += f"_sf{spec['sf']:g}"
    if spec["tol"] != DEFAULT_SPEC["tol"]:
        name += f"_tol{spec['tol']:g}"
    if spec["coupling"] != "iqnils":
        name += f"_{spec['coupling']}"
    if spec["fluidtol"] != "default":
        name += f"_fluid{spec['fluidtol']}"
    if spec["pimple"] != 1:
        name += f"_pimple{spec['pimple']:g}"
    if spec["pc"] != "auto":
        name += f"_pc{spec['pc']}"
    if spec["lag"] != -2:
        name += f"_lag{spec['lag']:g}"
    if spec["pcranks"]:
        name += f"_pcranks{spec['pcranks']:g}"
    if spec["solidranks"]:
        name += f"_solidranks{spec['solidranks']:g}"
    return name


def ignored(directory: str, names: list[str]) -> set[str]:
    skip = vmod.ignored(directory, names)
    if Path(directory) == TUTORIAL:
        skip.update({"verification", "validation"}.intersection(names))
    return skip


ROBIN_FSI = """FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "constant";
    object      fsiProperties;
}

// Robin-Neumann coupling written by the verification driver: the fluid
// interface uses elasticWallPressure/elasticWallVelocity and unrelaxed
// fixed-point iterations
fluidSolidInterface    fixedRelaxation;

fixedRelaxationCoeffs
{
    solidPatch interface;
    fluidPatch interface;
    predictor            yes;
    predictSolid         yes;
    relaxationFactor     1.0;
    outerCorrTolerance   TOL;
    robinPressureTolerance TOL10;
    robinFluxTolerance     5e-3;
    nOuterCorr           100;
    coupled              yes;
    interfaceTransferMethod AMI;
    writeResidualsToFile yes;
}
"""


def prepare_run(spec: dict, cores: int, name: str) -> tuple[Path, dict]:
    case = WORK_ROOT / name
    if case.exists():
        shutil.rmtree(case)
    shutil.copytree(TUTORIAL, case, symlinks=True, ignore=ignored)
    info: dict = {}

    # Fluid mesh
    if spec["fluid"] in FLUID_LEVELS:
        write_mesh_dict(case, spec["fluid"])
        info["fluid_sizes_m"] = fluid_sizes(spec["fluid"])

    # Solid mesh and centreline monitors on its nodes
    nx, ny, nz = SOLID_LEVELS[spec["solid"]]
    block = case / "system" / "solid" / "blockMeshDict"
    for key, value in zip(("nx", "ny", "nz"), (nx, ny, nz)):
        replace_entry(block, key, str(value), True)
    info["centreline_nodes_mm"] = write_centreline_monitors(case, nz)

    # Material and stabilisation
    set_material(case / "constant" / "solid" / "mechanicalProperties",
                 spec["mu"], NU)
    replace_entry(case / "constant" / "solid" / "solidProperties",
                  "scaleFactor", f"{spec['sf']:g}")

    # Time
    end_time, delta_t = spec["T"], spec["dt"]
    control = case / "system" / "controlDict"
    replace_entry(control, "startFrom", "latestTime", True)
    replace_entry(control, "endTime", f"{end_time:.10g}", True)
    replace_entry(control, "deltaT", f"{delta_t:.10g}", True)
    # Fields every 5 s (restart points) and at the end
    write_every = end_time if end_time <= 10 else 5.0
    replace_entry(control, "writeInterval", f"{write_every:.10g}", True)
    replace_text(control, '    #include "flapCentreline"\n',
                 '    #include "flapCentreline"\n\n'
                 '    // Monitors of the verification study\n'
                 '    #include "verificationMonitors"\n'
                 '    #include "voxelHistory"\n')
    avg_start = end_time * (1.0 - WINDOW_FRACTION)
    (case / "system" / "verificationMonitors").write_text(
        verification_monitors(avg_start))
    (case / "system" / "voxelHistory").write_text(
        voxel_history(avg_start, end_time, delta_t))

    # Coupling
    fsi = case / "constant" / "fsiProperties"
    if spec["coupling"] == "robin":
        fsi.write_text(ROBIN_FSI.replace("TOL10", f"{10 * spec['tol']:g}")
                       .replace("TOL", f"{spec['tol']:g}"))
        p = case / "0" / "fluid" / "p"
        text = p.read_text()
        text, found = re.subn(r"(interface\s*\{\s*type\s+)zeroGradient;",
                              r"\1elasticWallPressure;\n        value"
                              r"           uniform 0;", text)
        if found != 1:
            fail(f"Could not set the Robin pressure condition in {p}")
        p.write_text(text)
        u = case / "0" / "fluid" / "U"
        text = u.read_text()
        text, found = re.subn(r"(interface\s*\{\s*type\s+)"
                              r"newMovingWallVelocity;",
                              r"\1elasticWallVelocity;", text)
        if found != 1:
            fail(f"Could not set the Robin velocity condition in {u}")
        u.write_text(text)
    else:
        replace_entry(fsi, "outerCorrTolerance", f"{spec['tol']:g}")
        if spec["tol"] < 1e-4:
            replace_entry(fsi, "nOuterCorr", "60")

    # Fluid linear-solver tolerances
    if spec["fluidtol"] == "tight":
        solution = case / "system" / "fluid" / "fvSolution"
        text = solution.read_text()
        text, n1 = re.subn(r'("p\|pFinal"\s*\{[^}]*?tolerance\s+)[^;]+;',
                           r"\g<1>1e-10;", text, flags=re.DOTALL)
        text, n2 = re.subn(r'("p\|pFinal"\s*\{[^}]*?relTol\s+)[^;]+;',
                           r"\g<1>0;", text, flags=re.DOTALL)
        text, n3 = re.subn(r'("U\|UFinal"\s*\{[^}]*?tolerance\s+)[^;]+;',
                           r"\g<1>1e-11;", text, flags=re.DOTALL)
        text, n4 = re.subn(r'("cellMotionU\|cellMotionUFinal"\s*\{[^}]*?'
                           r'tolerance\s+)[^;]+;', r"\g<1>1e-11;", text,
                           flags=re.DOTALL)
        if (n1, n2, n3, n4) != (1, 1, 1, 1):
            fail(f"Could not tighten the fluid tolerances in {solution}")
        solution.write_text(text)
    if spec["pimple"] != 1:
        solution = case / "system" / "fluid" / "fvSolution"
        replace_entry(solution, "nOuterCorrectors", f"{int(spec['pimple'])}")

    # Decomposition and the solid's parallel linear solver
    vmod.set_subdomains(case, cores, ("", "fluid", "solid"))
    info["solid_preconditioner"] = solid_preconditioner(spec, cores)
    if info["solid_preconditioner"] == "telescope":
        info["solid_preconditioner_ranks"] = telescope_ranks(spec, cores)
    solution = case / "system" / "solid" / "fvSolution"
    if info["solid_preconditioner"] == "mumps":
        vmod.use_parallel_lu(solution)
    elif info["solid_preconditioner"] == "hypre":
        use_hypre(solution)
    elif info["solid_preconditioner"].startswith("telescope"):
        use_telescope(solution, cores, telescope_ranks(spec, cores))
    if spec["lag"] != -2:
        # Rebuild the solid preconditioner every 'lag' Newton iterations
        # (the tutorial builds it once and keeps it for the whole run)
        replace_entry(solution, "snes_lag_preconditioner",
                      f'"{int(spec["lag"])}"')
    return case, info


def solid_preconditioner(spec: dict, cores: int) -> str:
    """The solid's linear preconditioner (it does not change the solution).

    Serial runs keep the tutorial's LU (block Jacobi with one block). In
    parallel, hypre BoomerAMG on all ranks: on 128-256 ranks it was 6-10
    times faster per time step than an exact MUMPS LU on all ranks or on a
    sub-communicator of 4-16 ranks (PCTELESCOPE), although it needs about
    1.6 times more Krylov iterations.
    """
    if spec["pc"] != "auto":
        return spec["pc"]
    if cores == 1:
        return "bjacobi"
    return PC_AUTO_RULE(spec, cores)


def PC_AUTO_RULE(spec: dict, cores: int) -> str:
    return "hypre"


def build_meshes(case: Path, spec: dict, cores: int) -> dict:
    """Build both meshes as the tutorial's Allrun does; record metrics."""
    timings = {}
    timings["surface_s"] = run(["python3", "makeFluidSurface.py"], case,
                               "log.python3")
    mesh_dict = case / "system" / "meshDict"
    if mesh_dict.exists() or mesh_dict.is_symlink():
        mesh_dict.unlink()
    mesh_dict.symlink_to(f"meshDict.{spec['fluid']}")
    # cartesianMesh is multi-threaded (OpenMP)
    timings["cartesianMesh_s"] = run(["cartesianMesh"], case,
                                     "log.cartesianMesh", threads=cores)
    (case / "constant" / "fluid").mkdir(exist_ok=True)
    shutil.move(str(case / "constant" / "polyMesh"),
                str(case / "constant" / "fluid" / "polyMesh"))
    run(["blockMesh", "-region", "solid"], case, "log.blockMesh")
    if spec.get("solidranks"):
        write_solid_decomposition(case, spec, cores)
    run(["checkMesh", "-region", "fluid"], case, "log.checkMesh.fluid")
    run(["checkMesh", "-region", "solid"], case, "log.checkMesh.solid")
    metrics = {
        "fluid": checkmesh_metrics(case / "log.checkMesh.fluid"),
        "solid": checkmesh_metrics(case / "log.checkMesh.solid"),
        "timings": timings,
    }
    metrics["fluid"]["interface_faces"] = patch_faces(
        case / "constant" / "fluid" / "polyMesh" / "boundary", "interface")
    for patch in ("upperInlet", "lowerInlet", "outlet"):
        metrics["fluid"][f"{patch}_faces"] = patch_faces(
            case / "constant" / "fluid" / "polyMesh" / "boundary", patch)
    metrics["solid"]["interface_faces"] = patch_faces(
        case / "constant" / "solid" / "polyMesh" / "boundary", "interface")
    fluid = metrics["fluid"]
    if fluid.get("total_volume_m3") and fluid.get("cells"):
        fluid["h_eff_mm"] = 1e3 * (fluid["total_volume_m3"]
                                   / fluid["cells"]) ** (1 / 3)
    return metrics


def write_solid_decomposition(case: Path, spec: dict, cores: int) -> None:
    """Put the whole flap on the first 'solidranks' ranks, in z slabs.

    The flap is small next to the fluid: spread over hundreds of ranks, the
    solid's matrix-free Krylov iterations are dominated by communication.
    The other ranks get no solid cells (the fluid still uses every rank).
    """
    ranks = int(spec["solidranks"])
    nx, ny, nz = SOLID_LEVELS[spec["solid"]]
    labels = [min(ranks - 1, (c // (nx * ny)) * ranks // nz)
              for c in range(nx * ny * nz)]
    (case / "constant" / "solid" / "cellDecomposition").write_text(
        "FoamFile\n{\n    version     2.0;\n    format      ascii;\n"
        "    class       labelList;\n    object      cellDecomposition;\n}\n\n"
        f"{len(labels)}\n(\n" + "\n".join(map(str, labels)) + "\n)\n")
    path = case / "system" / "solid" / "decomposeParDict"
    text = path.read_text()
    text, found = re.subn(r"^method\s+\w+;", "method          manual;\n\n"
                          "manualCoeffs\n{\n    dataFile    "
                          '"cellDecomposition";\n}', text, flags=re.MULTILINE)
    if found != 1:
        fail(f"Could not set the manual decomposition in {path}")
    path.write_text(text)


def run_meshes(args: argparse.Namespace) -> bool:
    """Build every requested fluid level with the S2 flap; no solve."""
    levels = args.levels or list(FLUID_LEVELS)
    rows = []
    for level in levels:
        spec = dict(DEFAULT_SPEC, fluid=level)
        name = f"mesh_{level}"
        case, info = prepare_run(spec, 1, name)
        print(f"Meshing {level}", flush=True)
        metrics = build_meshes(case, spec, args.threads)
        row = {"level": level, **{f"fluid_{k}": v for k, v in
                                  metrics["fluid"].items()},
               **{f"size_{k}_mm": 1e3 * v
                  for k, v in info["fluid_sizes_m"].items()},
               "cartesianMesh_s": metrics["timings"]["cartesianMesh_s"]}
        rows.append(row)
        print(json.dumps(row, indent=1), flush=True)
        # Keep only the logs and the boundary file
        for path in (case / "constant" / "fluid" / "polyMesh").iterdir():
            if path.name != "boundary":
                path.unlink()
    OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)
    path = OUTPUT_ROOT / "fluid_meshes.json"
    existing = json.loads(path.read_text()) if path.is_file() else []
    existing = [r for r in existing if r["level"] not in levels] + rows
    existing.sort(key=lambda r: r["level"])
    path.write_text(json.dumps(existing, indent=2) + "\n")
    write_csv(OUTPUT_ROOT / "fluid_meshes.csv", fmt_rows(existing))
    return True


def read_set_samples(path: Path, values: int) -> dict[tuple, list[float]]:
    """Samples of a raw cloud set, keyed by the point on a 10 um grid."""
    lookup = {}
    for row in numeric_rows(path):
        if len(row) != 3 + values:
            fail(f"{path}: expected {3 + values} columns, found {len(row)}")
        lookup[tuple(round(v * 1.0e5) for v in row[:3])] = row[3:]
    return lookup


def voxel_average(velocity: list[dict], points: list[tuple],
                  lookups: list[dict]) -> list[dict]:
    """Average samples over each voxel's points and over the sample times.

    A voxel is compared only if at least half of its points lie in the
    computed fluid (at every sample time).
    """
    per_voxel = math.prod(vmod.VOXEL_SAMPLES)
    sums = [[0.0, 0.0, 0.0, 0] for _ in velocity]
    counts = [0] * len(velocity)
    for lookup in lookups:
        found = [0] * len(velocity)
        for index, point in points:
            value = lookup.get(tuple(round(v * 1.0e5) for v in point))
            if value is not None:
                for d in range(3):
                    sums[index][d] += value[d]
                sums[index][3] += 1
                found[index] += 1
        for index, n in enumerate(found):
            if n >= 0.5 * per_voxel:
                counts[index] += 1
    voxels = []
    for row, total, count in zip(velocity, sums, counts):
        entry = dict(row)
        entry["fluid_fraction"] = total[3] / (per_voxel * len(lookups))
        valid = count == len(lookups) and total[3] > 0
        for d, comp in enumerate(("vx", "vy", "vz")):
            entry[f"{comp}_sim_mm_s"] = (1.0e3 * total[d] / total[3]
                                         if valid else None)
        voxels.append(entry)
    return voxels


def sample_field_at_voxels(case: Path, cores: int, velocity: list[dict],
                           time_name: str, field: str) -> list[dict]:
    """Voxel averages of one volVectorField at one written time."""
    points = vmod.voxel_points(velocity)
    function = f"sampleVoxels{field}"
    lines = "\n".join(f"                ({p[0]:.8g} {p[1]:.8g} {p[2]:.8g})"
                      for _, p in points)
    (case / "system" / function).write_text(f"""FoamFile
{{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      {function};
}}

functions
{{
    {function}
    {{
        type            sets;
        libs            (sampling);
        region          fluid;
        interpolationScheme cellPoint;
        setFormat       raw;
        fields          ({field});
        sets
        {{
            voxels
            {{
                type    cloud;
                axis    xyz;
                points
                (
{lines}
                );
            }}
        }}
    }}
}}
""")
    command = mpirun(cores, ["postProcess", "-region", "fluid", "-dict",
                             f"system/{function}", "-fields", f"({field})",
                             "-time", time_name])
    run(command, case, f"log.{function}")
    files = [path for path in case.glob(f"postProcessing/**/voxels_{field}.xy")
             if function in path.parts and path.parent.name == time_name]
    if len(files) != 1:
        fail(f"Expected one voxel sample file of {field} in {case}, "
             f"found {len(files)}")
    return voxel_average(velocity, points, [read_set_samples(files[0], 3)])


def voxel_history_average(case: Path, velocity: list[dict], start: float,
                          end: float) -> tuple[list[dict], int]:
    """Time average of the fixed-voxel velocity samples in (start, end]."""
    points = vmod.voxel_points(velocity)
    files = [path for path in case.glob("postProcessing/**/voxels_U.xy")
             if "voxelHistory" in path.parts
             and start + 1e-9 < float(path.parent.name) <= end + 1e-9]
    if len(files) < 5:
        fail(f"Only {len(files)} fixed-voxel velocity samples in the window "
             f"of {case}")
    lookups = [read_set_samples(path, 3) for path in sorted(files)]
    return voxel_average(velocity, points, lookups), len(files)


def jet_peaks(voxels: list[dict]) -> dict:
    """Largest voxel-averaged vz of each jet on each plane (mm/s)."""
    peaks = {}
    for plane in (10.0, 30.0):
        for jet, sign in (("upper", 1), ("lower", -1)):
            values = [v["vz_sim_mm_s"] for v in voxels
                      if abs(v["z_mm"] - plane) < 1e-6
                      and sign * v["y_mm"] > 0
                      and v["vz_sim_mm_s"] is not None]
            if values:
                peaks[f"peak_vz_{jet}_z{plane:g}_mm_s"] = max(values)
    return peaks


def solid_centreline_final(case: Path, cores: int, time_name: str) -> list:
    """Displacement at the material stations COMPARE_Z at the final time.

    Samples D (cell-point interpolation) on the undeformed centreline
    x = y = 0: a check of the node monitors that does not depend on where
    the nodes are.
    """
    points = "\n".join(f"                (0 0 {z * 1e-3:.10g})"
                       for z in COMPARE_Z)
    (case / "system" / "sampleCentreline").write_text(f"""FoamFile
{{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      sampleCentreline;
}}

functions
{{
    sampleCentreline
    {{
        type            sets;
        libs            (sampling);
        region          solid;
        interpolationScheme cellPoint;
        setFormat       raw;
        fields          (D);
        sets
        {{
            centreline
            {{
                type    cloud;
                axis    xyz;
                points
                (
{points}
                );
            }}
        }}
    }}
}}
""")
    command = mpirun(cores, ["postProcess", "-region", "solid", "-dict",
                             "system/sampleCentreline", "-fields", "(D)",
                             "-time", time_name])
    run(command, case, "log.sampleCentreline")
    files = [path for path in case.glob("postProcessing/**/centreline_D.xy")
             if "sampleCentreline" in path.parts
             and path.parent.name == time_name]
    if len(files) != 1:
        fail(f"Expected one centreline sample file in {case}")
    # A station on a processor boundary is sampled by both processors
    unique = {}
    for row in numeric_rows(files[0]):
        unique.setdefault(round(row[2] * 1e6), row)
    rows = [unique[key] for key in sorted(unique)]
    if len(rows) != len(COMPARE_Z):
        fail(f"{files[0]}: {len(rows)} of {len(COMPARE_Z)} stations sampled")
    # (z0, Dy, Dz) in mm
    return [[1e3 * r[2], 1e3 * r[4], 1e3 * r[5]] for r in rows]


def material_stations(nodes: list[float], dy: list[float],
                      dz: list[float]) -> list[list[float]]:
    """Displacement (Dy, Dz) at COMPARE_Z by linear interpolation in z0."""
    curve_y = list(zip([0.0] + nodes, [0.0] + dy))
    curve_z = list(zip([0.0] + nodes, [0.0] + dz))
    return [[z, vmod.interpolate(curve_y, z), vmod.interpolate(curve_z, z)]
            for z in COMPARE_Z]


def evaluate_run(case: Path, spec: dict, cores: int, info: dict,
                 meshes: dict) -> dict:
    end_time, delta_t = spec["T"], spec["dt"]
    window = WINDOW_FRACTION * end_time
    w2 = (end_time - window, end_time)
    w1 = (end_time - 2 * window, end_time - window)
    result = {"name": spec_name(spec), "spec": spec, "ranks": cores,
              "window_s": list(w2), "previous_window_s": list(w1)}
    result["meshes"] = meshes
    result["fluid_sizes_mm"] = {k: 1e3 * v for k, v in
                                info.get("fluid_sizes_m", {}).items()}
    nx, ny, nz = SOLID_LEVELS[spec["solid"]]
    result["solid_mesh"] = {"nx": nx, "ny": ny, "nz": nz,
                            "cells": nx * ny * nz,
                            "h_mm": [FLAP_MM[0] / nx, FLAP_MM[1] / ny,
                                     FLAP_MM[2] / nz]}

    # Tip
    tip = vmod.tip_history(case, end_time, delta_t)
    q = {"tip_y_final_mm": 1e3 * tip[-1][2],
         "tip_z_final_mm": 65.0 + 1e3 * tip[-1][3]}
    for label, (a, b) in (("w2", w2), ("w1", w1)):
        stats_y = window_stats(tip, 2, a, b)
        stats_z = window_stats(tip, 3, a, b)
        q[f"tip_y_{label}_mean_mm"] = 1e3 * stats_y["mean"]
        q[f"tip_y_{label}_std_mm"] = 1e3 * stats_y["std"]
        q[f"tip_y_{label}_range_mm"] = 1e3 * (stats_y["max"] - stats_y["min"])
        q[f"tip_z_{label}_mean_mm"] = 65.0 + 1e3 * stats_z["mean"]
    q["tip_y_window_change_mm"] = q["tip_y_w2_mean_mm"] - q["tip_y_w1_mean_mm"]
    tail = [r[2] for r in tip if r[0] > 0.9 * end_time + 1e-9]
    q["tip_drift_last_tenth_mm"] = 1e3 * (max(tail) - min(tail))
    # Primary steady-state values: means over the last third of the run,
    # with their time (wander and drift) uncertainty
    for label, column in (("y", 2), ("z", 3)):
        settle = settling_estimate(tip, column, end_time, 1e3)
        offset = 65.0 if label == "z" else 0.0
        q[f"tip_{label}_mean_mm"] = offset + settle["mean"]
        q[f"tip_{label}_time_uncertainty_mm"] = settle["uncertainty"]
        q[f"tip_{label}_half_difference_mm"] = settle["half_difference"]
        q[f"tip_{label}_batch_se_mm"] = settle["batch_standard_error"]
    result["tip_history"] = [[round(r[0], 6), 1e3 * r[2], 1e3 * r[3]]
                             for r in tip
                             if abs(r[0] * 50 - round(r[0] * 50)) < 1e-6]

    # Centreline from the node monitors: the window-mean and final
    # displacements at the nodes, in material coordinates, and the deformed
    # line y(z) for the validation
    names = sorted({p.name[len("solidPointDisplacement_"):-len(".dat")]
                    for p in case.glob("postProcessing/*/"
                                       "solidPointDisplacement_flapCentrelineN*.dat")})
    nodes = info["centreline_nodes_mm"]
    if len(names) != len(nodes):
        fail(f"Expected {len(nodes)} centreline monitors in {case}, "
             f"found {len(names)}")
    mean_dy, mean_dz, final_dy, final_dz = [], [], [], []
    for name in names:
        rows = vmod.monitor_history(case, name, end_time, delta_t)
        last_third = (end_time * 2.0 / 3.0, end_time)
        mean_dy.append(1e3 * window_stats(rows, 2, *last_third)["mean"])
        mean_dz.append(1e3 * window_stats(rows, 3, *last_third)["mean"])
        final_dy.append(1e3 * rows[-1][2])
        final_dz.append(1e3 * rows[-1][3])
    result["centreline_nodes_mm"] = nodes
    result["centreline_mean_D_mm"] = [mean_dy, mean_dz]
    result["stations_mean"] = material_stations(nodes, mean_dy, mean_dz)
    result["stations_final"] = material_stations(nodes, final_dy, final_dz)
    result["centreline_mean"] = sorted(
        [(0.0, 0.0)] + [(z + dz, dy) for z, dy, dz in
                        zip(nodes, mean_dy, mean_dz)])
    result["centreline_final"] = sorted(
        [(0.0, 0.0)] + [(z + dz, dy) for z, dy, dz in
                        zip(nodes, final_dy, final_dz)])
    time_name = vmod.latest_time(case, cores)
    if not math.isclose(float(time_name), end_time, rel_tol=1e-9):
        fail(f"The latest field time of {case} is {time_name}")
    sampled = solid_centreline_final(case, cores, time_name)
    result["stations_final_from_D"] = sampled
    q["centreline_monitor_vs_field_max_mm"] = max(
        math.hypot(a[1] - b[1], a[2] - b[2])
        for a, b in zip(result["stations_final"], sampled))

    # Fluid force on the flap and the inlet pressures
    force = history_rows(case, "flapForce", "force.dat", end_time, delta_t, 4)
    for comp, column in (("x", 1), ("y", 2), ("z", 3)):
        q[f"force_{comp}_w2_mean_N"] = window_stats(force, column, *w2)["mean"]
        q[f"force_{comp}_w1_mean_N"] = window_stats(force, column, *w1)["mean"]
    q["force_y_w2_std_N"] = window_stats(force, 2, *w2)["std"]
    for comp, column in (("y", 2), ("z", 3)):
        settle = settling_estimate(force, column, end_time)
        q[f"force_{comp}_mean_N"] = settle["mean"]
        q[f"force_{comp}_time_uncertainty_N"] = settle["uncertainty"]
    for inlet in ("upperInlet", "lowerInlet"):
        pressure = history_rows(case, f"{inlet}Pressure",
                                "surfaceFieldValue.dat", end_time, delta_t, 2)
        # Kinematic pressure: Pa = rho*p
        q[f"dp_{inlet}_w2_mean_Pa"] = 1163.3 * window_stats(
            pressure, 1, *w2)["mean"]
        q[f"dp_{inlet}_w1_mean_Pa"] = 1163.3 * window_stats(
            pressure, 1, *w1)["mean"]
        settle = settling_estimate(pressure, 1, end_time, 1163.3)
        q[f"dp_{inlet}_mean_Pa"] = settle["mean"]
        q[f"dp_{inlet}_time_uncertainty_Pa"] = settle["uncertainty"]
    probes = history_rows(case, "velocityProbes", "U", end_time, delta_t,
                          1 + 3 * len(PROBES))
    for i in range(len(PROBES)):
        column = 1 + 3 * i + 2
        q[f"probe{i}_vz_w2_mean_mm_s"] = 1e3 * window_stats(
            probes, column, *w2)["mean"]
        q[f"probe{i}_vz_w1_mean_mm_s"] = 1e3 * window_stats(
            probes, column, *w1)["mean"]
        q[f"probe{i}_vz_w2_std_mm_s"] = 1e3 * window_stats(
            probes, column, *w2)["std"]
        settle = settling_estimate(probes, column, end_time, 1e3)
        q[f"probe{i}_vz_mean_mm_s"] = settle["mean"]
        q[f"probe{i}_vz_time_uncertainty_mm_s"] = settle["uncertainty"]

    # Velocity on the MRI voxels: the time average of the fixed-voxel
    # samples over the window (primary), the final instantaneous field, and
    # the cell-based fieldAverage mean (check of the moving-mesh effect)
    velocity = read_csv(vmod.REFERENCE_DIR / "phaseI_velocity.csv")
    mean_voxels, samples = voxel_history_average(case, velocity, *w2)
    q["voxel_history_samples"] = samples
    final_voxels = sample_field_at_voxels(case, cores, velocity, time_name,
                                          "U")
    cell_mean_voxels = sample_field_at_voxels(case, cores, velocity,
                                              time_name, "UMean")
    result["voxels"] = {}
    for label, voxels in (("mean", mean_voxels), ("final", final_voxels),
                          ("cellmean", cell_mean_voxels)):
        result["voxels"][label] = [
            [v["x_mm"], v["y_mm"], v["z_mm"], v["vx_sim_mm_s"],
             v["vy_sim_mm_s"], v["vz_sim_mm_s"], v["fluid_fraction"]]
            for v in voxels]
        for key, value in jet_peaks(voxels).items():
            q[f"{key}_{label}"] = value
    check = velocity_difference(result["voxels"]["mean"],
                                result["voxels"]["cellmean"])
    q["voxel_mean_vs_cellmean_mean"] = check["mean"]
    q["voxel_mean_vs_cellmean_max"] = check["max"]

    # Deformed fluid mesh quality at the end
    run(mpirun(cores, ["checkMesh", "-region", "fluid", "-time", time_name]),
        case, "log.checkMesh.fluid.final")
    final_mesh = checkmesh_metrics(case / "log.checkMesh.fluid.final")
    q["fluid_final_max_non_orthogonality_deg"] = final_mesh.get(
        "max_non_orthogonality_deg", math.nan)
    q["fluid_final_max_skewness"] = final_mesh.get("max_skewness", math.nan)
    text = (case / "log.solids4Foam").read_text(errors="replace")
    ami = re.findall(r"interface-to-interface face error:\s*([0-9.eE+-]+)",
                     text)
    q["ami_face_error_last"] = float(ami[-1]) if ami else math.nan

    # Cost and coupling
    counts: dict[str, int] = {}
    clock = 0.0
    for log in sorted(case.glob("log.solids4Foam*")):
        text = log.read_text(errors="replace")
        for time_value, iteration in re.findall(
                r"^Time = ([0-9.eE+-]+), iteration: (\d+)", text,
                re.MULTILINE):
            counts[time_value] = max(counts.get(time_value, 0),
                                     int(iteration))
        times = re.findall(r"ClockTime\s*=\s*([0-9.eE+-]+)", text)
        clock += float(times[-1]) if times else 0.0
    q["fsi_iterations_mean"] = sum(counts.values()) / len(counts)
    q["fsi_iterations_max"] = max(counts.values())
    q["clock_time_h"] = clock / 3600.0
    q["restart_segments"] = len(list(case.glob("log.solids4Foam*")))
    q["fsi_max_iteration_hits"] = count_max_iterations(case, spec)
    result["qoi"] = q
    return result


def count_max_iterations(case: Path, spec: dict) -> int:
    """Time steps whose FSI loop stopped at nOuterCorr without converging."""
    text = (case / "log.solids4Foam").read_text(errors="replace")
    limit = 100 if spec["coupling"] == "robin" else (
        60 if spec["tol"] < 1e-4 else 30)
    counts: dict[str, int] = {}
    for time_value, iteration in re.findall(
            r"^Time = ([0-9.eE+-]+), iteration: (\d+)", text, re.MULTILINE):
        counts[time_value] = max(counts.get(time_value, 0), int(iteration))
    return sum(1 for value in counts.values() if value >= limit)


def run_coupled(args: argparse.Namespace) -> bool:
    if not args.spec:
        fail("--study run needs --spec")
    spec = parse_spec(args.spec)
    if args.quick:
        spec["T"] = 0.02
    cores = args.cores
    name = spec_name(spec)
    case = WORK_ROOT / name
    settings = {"study": "run", "spec": spec, "cores": cores,
                "driver": vmod.hashlib.sha256(SCRIPT.read_bytes()).hexdigest(),
                "tutorial_inputs": vmod.tutorial_fingerprint(),
                "build": vmod.build_fingerprint()}
    info_file = case / "verification_info.json"
    if args.restart:
        # Continue an interrupted run from its latest written time
        info = json.loads(info_file.read_text())
        meshes = info.pop("meshes")
        segment = 1
        while (case / f"log.solids4Foam.{segment}").exists():
            segment += 1
        shutil.move(str(case / "log.solids4Foam"),
                    str(case / f"log.solids4Foam.{segment}"))
        print(f"Restarting {name} on {cores} ranks", flush=True)
        run(mpirun(cores, ["solids4Foam"]), case, "log.solids4Foam")
        vmod.check_solver_log(case, name)
        (case / SETTINGS_FILE).write_text(json.dumps(settings, indent=2)
                                          + "\n")
    elif args.evaluate_only:
        info = json.loads(info_file.read_text())
        meshes = info.pop("meshes")
    elif vmod.reusable(case, settings, args, name):
        info = json.loads(info_file.read_text())
        meshes = info.pop("meshes")
    else:
        case, info = prepare_run(spec, cores, name)
        print(f"Building the meshes of {name}", flush=True)
        meshes = build_meshes(case, spec, cores)
        info_file.write_text(json.dumps({**info, "meshes": meshes},
                                        indent=2) + "\n")
        if cores > 1:
            run(["decomposePar", "-region", "fluid"], case,
                "log.decomposePar.fluid")
            run(["decomposePar", "-region", "solid"], case,
                "log.decomposePar.solid")
        print(f"Running {name} on {cores} ranks", flush=True)
        run(mpirun(cores, ["solids4Foam"]), case, "log.solids4Foam")
        vmod.check_solver_log(case, name)
        (case / SETTINGS_FILE).write_text(json.dumps(settings, indent=2)
                                          + "\n")
    if args.quick:
        print(f"{name}: completed the quick run")
        return True
    vmod.check_solver_log(case, name)
    result = evaluate_run(case, spec, cores, info, meshes)
    result["settings"] = settings
    RUNS_ROOT.mkdir(parents=True, exist_ok=True)
    (RUNS_ROOT / f"{name}.json").write_text(json.dumps(result, indent=1)
                                            + "\n")
    print(json.dumps(result["qoi"], indent=1))
    return True


# ---------------------------------------------------------------------------
# Analysis: successive differences, observed orders, uncertainties
# ---------------------------------------------------------------------------

def three_level(values: list[float], ratio: float, noise: float = 0.0,
                formal: float = 2.0) -> dict:
    """Observed order and Richardson estimate of a refinement triplet.

    values are ordered coarse, medium, fine with constant nominal ratio.
    'noise' is the resolution of the QoI (e.g. its residual transient
    error): differences below it are not interpreted. Monotone convergence
    with a convergence ratio in (0, 1) and an order in [0.5, 2*formal] is a
    screening criterion consistent with (not proof of) the asymptotic range;
    then a Richardson estimate and a fine-level GCI (Fs = 1.25) are given.
    Otherwise no order is claimed and a heuristic difference envelope
    (3 x the largest successive difference) is reported instead.
    """
    f1, f2, f3 = values
    e21, e32 = f2 - f1, f3 - f2
    out = {"coarse": f1, "medium": f2, "fine": f3,
           "diff_medium_coarse": e21, "diff_fine_medium": e32,
           "ratio": ratio, "resolution": noise, "order": math.nan,
           "richardson": math.nan, "gci_fine": math.nan,
           "envelope": 3.0 * max(abs(e21), abs(e32))}
    if max(abs(e21), abs(e32)) <= noise:
        out["status"] = "changes below the QoI resolution"
        return out
    if abs(e32) <= noise:
        out["status"] = ("fine-level change below the QoI resolution: "
                         "no order")
        return out
    if abs(e21) <= noise:
        out["status"] = "coarse-level change below resolution: no order"
        return out
    convergence_ratio = e32 / e21
    out["convergence_ratio"] = convergence_ratio
    if 0 < convergence_ratio < 1:
        order = math.log(1 / convergence_ratio) / math.log(ratio)
        out["order"] = order
        if 0.5 <= order <= 2 * formal:
            out["status"] = "monotone convergence"
            out["richardson"] = f3 + e32 / (ratio ** order - 1)
            out["gci_fine"] = 1.25 * abs(e32) / (ratio ** order - 1)
        else:
            out["status"] = "monotone, order outside [0.5, 2 x formal]"
    elif convergence_ratio < 0:
        out["status"] = "oscillatory (difference changes sign)"
    else:
        out["status"] = "divergent (differences grow)"
    return out


def load_runs() -> dict[str, dict]:
    runs = {}
    for path in sorted(RUNS_ROOT.glob("*.json")):
        data = json.loads(path.read_text())
        runs[data["name"]] = data
    return runs


def station_difference(a: list, b: list) -> tuple[float, float]:
    """RMS and max displacement-vector difference (mm) at COMPARE_Z."""
    diffs = [math.hypot(pa[1] - pb[1], pa[2] - pb[2])
             for pa, pb in zip(a, b)
             if math.isfinite(pa[1]) and math.isfinite(pb[1])]
    if len(diffs) != len(COMPARE_Z):
        fail("Centreline stations are missing")
    return (math.sqrt(sum(d * d for d in diffs) / len(diffs)), max(diffs))


def velocity_difference(a: list, b: list, scale: float = 630.0) -> dict:
    """Mean and max vector difference of two voxel sets, normalised.

    Only voxels valid in both sets are compared (a common mask).
    """
    distances = []
    for va, vb in zip(a, b):
        if va[3] is None or vb[3] is None:
            continue
        distances.append(math.sqrt(sum((va[k] - vb[k]) ** 2
                                       for k in (3, 4, 5))))
    if not distances:
        return {"mean": math.nan, "max": math.nan, "voxels": 0}
    return {"mean": sum(distances) / len(distances) / scale,
            "max": max(distances) / scale, "voxels": len(distances)}


def validation_metrics(run: dict, reference: dict) -> dict:
    """Comparison with the measurements (validation, not verification)."""
    measured_line = read_csv(vmod.REFERENCE_DIR / "phaseI_centreline.csv")
    velocity = read_csv(vmod.REFERENCE_DIR / "phaseI_velocity.csv")
    tip_target = float(reference["targets"]["phase_I_tip_y_mm"])
    q = run["qoi"]
    out = {"tip_error_final_mm": q["tip_y_final_mm"] - tip_target,
           "tip_error_mean_mm": q["tip_y_mean_mm"] - tip_target}
    for label in ("mean", "final"):
        line = [tuple(p) for p in run[f"centreline_{label}"]]
        errors = vmod.centreline_errors(line, measured_line)
        out[f"centreline_rms_{label}_mm"] = errors["rms_mm"]
        out[f"centreline_max_{label}_mm"] = errors["max_mm"]
        out[f"centreline_mean_{label}_mm"] = errors["mean_mm"]
    for label in ("mean", "final"):
        voxels = []
        for row, sim in zip(velocity, run["voxels"][label]):
            entry = dict(row)
            for k, comp in enumerate(("vx", "vy", "vz")):
                entry[f"{comp}_sim_mm_s"] = sim[3 + k]
            voxels.append(entry)
        errors = vmod.velocity_errors(voxels, 630.0)
        for key in ("voxels", "d_mean", "d_max", "d_mean_z10", "d_mean_z30",
                    "rms_vx_mm_s", "rms_vy_mm_s", "rms_vz_mm_s"):
            out[f"velocity_{key}_{label}"] = errors.get(key, math.nan)
    return out


# Principal QoIs: key, label, and the key of its resolution (residual
# transient error), if any
QOIS = (
    ("tip_y_mean_mm", "tip y, mean of last third (mm)",
     "tip_y_time_uncertainty_mm"),
    ("tip_z_mean_mm", "tip z, mean of last third (mm)",
     "tip_z_time_uncertainty_mm"),
    ("force_y_mean_N", "flap force y (N)", "force_y_time_uncertainty_N"),
    ("force_z_mean_N", "flap force z (N)", "force_z_time_uncertainty_N"),
    ("dp_upperInlet_mean_Pa", "upper inlet pressure (Pa)",
     "dp_upperInlet_time_uncertainty_Pa"),
    ("dp_lowerInlet_mean_Pa", "lower inlet pressure (Pa)",
     "dp_lowerInlet_time_uncertainty_Pa"),
    ("peak_vz_upper_z10_mm_s_mean", "upper jet peak vz, z = 10 mm (mm/s)",
     None),
    ("peak_vz_lower_z10_mm_s_mean", "lower jet peak vz, z = 10 mm (mm/s)",
     None),
    ("peak_vz_upper_z30_mm_s_mean", "upper jet peak vz, z = 30 mm (mm/s)",
     None),
    ("peak_vz_lower_z30_mm_s_mean", "lower jet peak vz, z = 30 mm (mm/s)",
     None),
)


def compare_pair(a: dict, b: dict) -> dict:
    rms, mx = station_difference(a["stations_mean"], b["stations_mean"])
    vel = velocity_difference(a["voxels"]["mean"], b["voxels"]["mean"])
    out = {"from": a["name"], "to": b["name"],
           "centreline_rms_change_mm": rms, "centreline_max_change_mm": mx,
           "velocity_mean_change": vel["mean"],
           "velocity_max_change": vel["max"],
           "velocity_compared_voxels": vel["voxels"]}
    out["tip_y_change_mm"] = (b["qoi"]["tip_y_mean_mm"]
                              - a["qoi"]["tip_y_mean_mm"])
    for key, _, _ in QOIS:
        out[f"{key}_change"] = b["qoi"].get(key, math.nan) - a["qoi"].get(
            key, math.nan)
    return out


def series_analysis(series: dict, runs: dict) -> dict:
    """Successive comparisons and explicit refinement triplets.

    series: {label: {"runs": [names, coarse to fine], "ratio": r,
                     "triplets": true}}. A triplet is analysed only when its
    three runs are all present and consecutive in the stated list.
    """
    out = {}
    for label, definition in series.items():
        names = definition["runs"]
        entry = {"runs": names, "ratio": definition.get("ratio"),
                 "missing": [n for n in names if n not in runs]}
        present = [n for n in names if n in runs]
        entry["successive"] = []
        for a, b in zip(names, names[1:]):
            if a in runs and b in runs:
                entry["successive"].append(compare_pair(runs[a], runs[b]))
        entry["values"] = {key: [runs[n]["qoi"].get(key, math.nan)
                                 if n in runs else None for n in names]
                           for key, _, _ in QOIS}
        entry["triplets"] = []
        ratio = definition.get("ratio")
        if ratio and definition.get("triplets", True):
            for i in range(len(names) - 2):
                trio = names[i:i + 3]
                if not all(n in runs for n in trio):
                    continue
                result = {"runs": trio}
                for key, _, noise_key in QOIS:
                    values = [runs[n]["qoi"].get(key, math.nan) for n in trio]
                    if not all(math.isfinite(v) for v in values):
                        continue
                    noise = max((runs[n]["qoi"].get(noise_key, 0.0)
                                 for n in trio), default=0.0) \
                        if noise_key else 0.0
                    result[key] = three_level(values, ratio, noise)
                # Norms of the field differences: convergence ratio only
                # (norms discard signs, so no order is claimed)
                s1 = compare_pair(runs[trio[0]], runs[trio[1]])
                s2 = compare_pair(runs[trio[1]], runs[trio[2]])
                for key in ("centreline_rms_change_mm",
                            "velocity_mean_change"):
                    result[f"{key}_ratio"] = (s2[key] / s1[key]
                                              if s1[key] > 0 else math.nan)
                entry["triplets"].append(result)
        out[label] = entry
    return out


def solid_limit_analysis() -> dict:
    """Richardson analysis of the solid study and the common-limit test."""
    path = OUTPUT_ROOT / "solid_study.json"
    if not path.is_file():
        return {}
    rows = json.loads(path.read_text())
    out = {}
    for sf in sorted({r["stabilisation"] for r in rows}):
        by_level = {r["level"]: r for r in rows if r["stabilisation"] == sf}
        entry = {}
        for trio in (("S1", "S2", "S3"), ("S2", "S3", "S4")):
            if all(level in by_level for level in trio):
                for key in ("mu_cal_Pa", "tip_at_mu_fixed_mm"):
                    noise = max(by_level[level].get(
                        "mu_cal_interp_uncertainty_Pa" if key == "mu_cal_Pa"
                        else "tip_interp_uncertainty_mm", 0.0)
                        for level in trio)
                    entry[f"{key}_{''.join(trio)}"] = three_level(
                        [by_level[level][key] for level in trio], 2.0, noise)
        out[f"sf{sf:g}"] = entry
    return out


def default_series() -> dict:
    """The refinement series of the verification campaign (README.md)."""
    def n(text: str) -> str:
        return spec_name(parse_spec(text))

    return {
        # F5 (11 M cells) was not affordable: the solid's matrix-free
        # Krylov iterations do not scale to the ~512 ranks it needs
        "fluid at S2": {"runs": [n(f"F{k}:S2") for k in range(1, 5)],
                        "ratio": FLUID_RATIO},
        "solid at F3": {"runs": [n(f"F3:S{k}") for k in range(1, 4)],
                        "ratio": 2.0},
        # Two levels of the r = 2 matched path; its third level, F5:S3,
        # was not affordable
        "matched": {"runs": [n("F1:S1"), n("F3:S2")], "ratio": 2.0},
        "solid at F1": {"runs": [n("F1:S2"), n("F1:S3")], "ratio": 2.0},
        "time": {"runs": [n("F2:S2"), n("F2:S2:dt=0.001"),
                          n("F2:S2:dt=0.0005")], "ratio": 2.0},
        # With the default fluid tolerances IQN-ILS stalls at a residual of
        # about 3e-5 (it fails at 1e-5 and 1e-6); with tight fluid
        # tolerances it reaches 1e-5 but stalls at 2.5-5e-6 (it fails at
        # 1e-6). The iterative study is therefore (default fluid, 1e-4),
        # (tight fluid, 1e-4) and (tight fluid, 1e-5).
        "coupling tolerance": {"runs": [n("F2:S2:fluidtol=tight"),
                                        n("F2:S2:tol=1e-5:fluidtol=tight")]},
        "all iterative tolerances": {
            "runs": [n("F2:S2"), n("F2:S2:tol=1e-5:fluidtol=tight")]},
        "fluid linear tolerance": {"runs": [n("F2:S2"),
                                            n("F2:S2:fluidtol=tight")]},
        "PIMPLE outer correctors": {"runs": [n("F2:S2"),
                                             n("F2:S2:pimple=3")]},
        "Robin-Neumann": {"runs": [n("F2:S2"), n("F2:S2:coupling=robin")]},
        "stabilisation": {"runs": [n("F2:S2"), n("F2:S2:sf=0.001")]},
        "end time": {"runs": [n("F2:S2"), n("F2:S2:T=40")]},
        "tutorial modulus": {"runs": [n("F1:S2:mu=58450"), n("F1:S2")]},
    }


def run_analyse(args: argparse.Namespace) -> bool:
    reference = json.loads(vmod.REFERENCE_FILE.read_text())
    runs = load_runs()
    table = []
    for name, data in runs.items():
        row = {"name": name, **{f"spec_{k}": v for k, v in data["spec"].items()},
               "ranks": data["ranks"],
               "fluid_cells": data["meshes"]["fluid"].get("cells"),
               "fluid_h_eff_mm": data["meshes"]["fluid"].get("h_eff_mm"),
               "fluid_interface_faces": data["meshes"]["fluid"].get(
                   "interface_faces"),
               "fluid_max_non_orthogonality_deg": data["meshes"]["fluid"].get(
                   "max_non_orthogonality_deg"),
               "solid_cells": data["solid_mesh"]["cells"],
               "solid_interface_faces": data["meshes"]["solid"].get(
                   "interface_faces"),
               **{f"fluid_size_{k}": v
                  for k, v in data.get("fluid_sizes_mm", {}).items()},
               **data["qoi"],
               **{f"validation_{k}": v
                  for k, v in validation_metrics(data, reference).items()}}
        table.append(row)
    OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)
    if table:
        write_csv(OUTPUT_ROOT / "verification_runs.csv", fmt_rows(table))
    series = (json.loads(Path(args.series).read_text()) if args.series
              else default_series())
    analysis = {"solid_study": solid_limit_analysis(),
                "series": series_analysis(series, runs)}
    (OUTPUT_ROOT / "verification_series.json").write_text(
        json.dumps(analysis, indent=2, default=str) + "\n")
    print(json.dumps(analysis["solid_study"], indent=1, default=str))
    return True


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--study", required=True,
                        choices=("solid", "meshes", "run", "analyse"))
    parser.add_argument("--levels",
                        help="comma-separated levels: S1-S4 for solid, "
                             "F1-F5 for meshes")
    parser.add_argument("--sf", default="0.01,0.001",
                        help="solid study stabilisation factors")
    parser.add_argument("--spec", help="run: the coupled run, e.g. F2:S3")
    parser.add_argument("--cores", type=int, default=1,
                        help="MPI ranks of a coupled run")
    parser.add_argument("--cores-solid", type=int, default=0,
                        help="ranks of the solid study's large meshes")
    parser.add_argument("--threads", type=int, default=1,
                        help="meshes: OpenMP threads of cartesianMesh")
    parser.add_argument("--jobs", type=int, default=1)
    parser.add_argument("--series",
                        help="analyse: JSON file naming the run series")
    parser.add_argument("--quick", action="store_true")
    parser.add_argument("--reuse", action="store_true")
    parser.add_argument("--evaluate-only", action="store_true",
                        help="run, solid: evaluate existing run directories")
    parser.add_argument("--restart", action="store_true",
                        help="run: continue an interrupted run")
    parser.add_argument("--work", help="work directory holding the runs "
                        "(default verification/work)")
    args = parser.parse_args()
    if args.levels:
        args.levels = [v.strip() for v in args.levels.split(",")]
    if args.work:
        global WORK_ROOT
        WORK_ROOT = Path(args.work).resolve()
    WORK_ROOT.mkdir(parents=True, exist_ok=True)
    OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)
    study = {"solid": run_solid, "meshes": run_meshes, "run": run_coupled,
             "analyse": run_analyse}[args.study]
    return 0 if study(args) else 1


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (RuntimeError, subprocess.SubprocessError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        sys.exit(2)
