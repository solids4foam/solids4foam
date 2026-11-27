#!/usr/bin/env python3
"""Verification of the oscillatingWallStokesLayer tutorial: the velocity and
the wall shear stress next to an immersed oscillating slab against the exact
solutions of the Stokes second problem, with the cutLink, penalty, legacy
incremental and (exploratory) ghostCell methods of immersedBoundaryForce.

Exact solutions (s = y - y_w, k = 1/delta, delta = sqrt(2 nu/omega), the slab
velocity U0 sin(omega t)):

  periodic half-space   u = U0 exp(-k s) sin(omega t - k s)
  wall shear stress     tau_w = rho nu du/ds(0) = -rho nu k U0 (sin + cos)
                        = rho nu k U0 sqrt(2) sin(omega t + 5 pi/4)
  finite height H       periodic solution with a slip lid at s = H:
                        u = Im{U0 e^{i omega t} cosh(l (H-s))/cosh(l H)},
                        l = (1 + i)/delta
  initial-value problem the exact solution of the simulated problem: fluid at
                        rest, slab velocity U0 sin(omega t) (zero at t = 0),
                        slip lid at s = H: the finite-height periodic
                        solution plus the decaying modes
                        sum a_n exp(-nu k_n^2 t) sin(k_n s),
                        k_n = (n + 1/2) pi/H, a_n from the initial condition
                        (N_MODES modes, Simpson integration with N_SIMPSON
                        intervals)

tau_w is the shear stress that the fluid exerts on the wall, i.e. the force on
the slab per unit wall area, per unit density when rho = 1. The solver
force (a force over the mesh depth) is divided by the wall area of the mesh.

The start-up transient has NOT decayed after four periods (the lowest lid
mode decays with rate nu k_1^2 = 0.14 1/s), so every metric is reported
against both the periodic solution (the official target) and the exact
initial-value solution of the simulated problem ('_ivp' columns).

Metrics over the fourth period 3 < t <= 4 s: the root mean square (RMS)
difference of the wall shear stress relative to the RMS of the exact shear;
the amplitude and phase of the fundamental of a least-squares fit of the
mean plus three harmonics (the amplitude difference relative to the exact
amplitude, the phase difference in degrees); the noise, the RMS of the
residual of that fit relative to the RMS of the exact shear; and, at the
eight written times 3.125 ... 4 s, the volume-weighted L1, L2 (RMS over the
snapshots) and Linf (maximum over the snapshots) norms of the velocity
difference relative to U0 in the fluid cells (cell centre above the wall):
'all', 'near' (s <= delta) and 'iface' (s < cell width).

The force estimators: for cutLink the momentum exchange (ME) and the surface
traction (ST) are two estimators of one run; the forcing estimator (FO) is a
separate run. For penalty, incremental and ghostCell the force is the
forcing plus the inertia of the fluid inside the slab (a first-order time
difference). A wall shear stress is never declared accurate from its
amplitude alone: the phase is reported next to it.

Each case is a copy of the tutorial under verification/work/cases/ with only
the method settings, the mesh level and the decomposition changed. The
tutorial itself is never modified. Standard library only. A manifest,
work/manifest.json, records every case with the SHA-256 of its inputs.

Usage:
    ./Allverify [run] [--runs ...]   generate and run locally (does not
                                     analyse or check: follow with the next
                                     two commands)
    ./Allverify unit                 exact tests that need no OpenFOAM
    ./Allverify dry-run --runs ...   list the cases and ranks
    ./Allverify generate --runs ...  generate the cases only
    ./Allverify status               state of the cases in the manifest
    ./Allverify analyse              analyse the completed cases
    ./Allverify check                PASS / FAIL / SKIPPED
    ./Allverify tables               the master table from the tracked data

The default run, 'smoke', is the cutLink case on Level 1 (20 x 2 cells, four
periods): a few seconds on one core.

Controls, all from this driver (the tracked rows of reference/results/controls
are what these named runs and options reproduce; see README.md):
    --runs dt-fine dt-L1 dt-6e-5 dt-3e-5   the time steps down to 3.125e-5 s
    --runs wall-face     cutLink on Levels 2-4 with the top vertices of the
                         slab moved onto the mesh face y = 0.125 m (--wall-y)
    --runs T8            the eight-period control (--end-time 8): cutLink
                         momentum exchange and penalty, Level 3
    --end-time T         endTime override of every selected case (tagged -T..)
    --wall-y Y           move the wall (top vertices of slab.stl) to y = Y
    --work-dir DIR       keep cases and manifest in DIR (independent jobs)
"""

import argparse
import cmath
import csv
import hashlib
import json
import math
import os
import re
import shutil
import subprocess
import sys
from datetime import datetime
from pathlib import Path

VERIFICATION_DIR = Path(__file__).resolve().parent.parent
TUTORIAL_DIR = VERIFICATION_DIR.parent
WORK_DIR = VERIFICATION_DIR / "work"
CASES_DIR = WORK_DIR / "cases"
MANIFEST = WORK_DIR / "manifest.json"
REFERENCE_DIR = VERIFICATION_DIR / "reference"
RESULTS_DIR = REFERENCE_DIR / "results"
CASES_CSV = RESULTS_DIR / "cases.csv"
CONTROLS_DIR = RESULTS_DIR / "controls"
CONTROLS_CSV = CONTROLS_DIR / "cases-timestep.csv"

LEVELS = (1, 2, 3, 4)
# (method, estimator): cutLink ME is the run that gives the ME and ST rows
VARIANTS = (("cutLink", "ME"), ("cutLink", "ST"), ("cutLink", "FO"),
            ("penalty", ""), ("incremental", ""), ("ghostCell", ""))
RUN_VARIANTS = (("cutLink", "ME"), ("cutLink", "FO"), ("penalty", ""),
                ("incremental", ""), ("ghostCell", ""))
ESTIMATOR_NAMES = {"ME": "momentumExchange", "ST": "surfaceTraction",
                   "FO": "forcing"}
N_OUTER = {"cutLink": 1, "penalty": 1, "incremental": 5, "ghostCell": 1}

# Settings written explicitly into every case: those of the stationary and
# oscillating cylinder verifications; ghostCell is documented for static
# bodies and is exploratory here (the slab does not change its position
# relative to the mesh, only its tangential velocity)
METHOD_SETTINGS = {
    "cutLink": ([], []),
    "penalty": ([("method", "penalty"), ("penaltyCoeff", "1e3"),
                 ("weighting", "volumeFraction"), ("surfaceRateCoeff", "3")],
                [("occupancy", "signedDistance")]),
    "incremental": ([("method", "incremental"), ("couplingCoeff", "0.2")],
                    [("occupancy", "vertexFraction")]),
    "ghostCell": ([("method", "ghostCell"), ("penaltyCoeff", "1e3"),
                   ("imageDistance", "1.5"),
                   ("imageInterpolation", "cellPoint")], []),
}

# Predefined analysis settings (fixed before the matrix was analysed)
PERIOD = (3.0, 4.0)                 # fourth period
NHARM = 3                           # harmonics of the noise definition
N_SNAPSHOTS = 8                     # written times in the period
N_MODES = 150                       # modes of the initial-value solution
N_SIMPSON = 1500                    # Simpson intervals of the mode integrals

# Predefined check tolerances, independent of the measured values
CHECK_SHEAR_AMP = 0.005             # |difference of the amplitude error|
CHECK_VELOCITY = 2e-4               # |difference of u_L2_all_ivp|
CHECK_CLOSED_FORM = 1e-9            # exact tau against rho nu du/ds
CHECK_SERIES = 2e-7                 # |shear history - tracked history| (Pa):
#                                     twice the rounding of the tracked files
CHECK_LONG_RUN = 1e-4               # |amplitude error, last period - fourth|
#                                     (the 0.01 % of the amplitude of
#                                     results.md, Exact solutions)

# Observed allocated rank-hours of the tracked runs are read from the tracked
# rows (observed_rank_hours)

# dt = 0 is the tutorial: the fixed time step 0.001 s (adjustTimeStep no). A
# positive dt sets deltaT in system/controlDict, the control that separates
# the time error from the mesh error. The steps divide the write interval
# 0.125 s, so that the fields are written on the phase grid of the analysis
# (0.002 s does not)
TUTORIAL_DT = 0.001
DT_VARIANTS = (("cutLink", "ME"), ("penalty", ""), ("incremental", ""),
               ("ghostCell", ""))
DT_STEPS = (0.0025, 0.0005, 0.00025)
DT_LEVELS = (2, 3, 4)

# Control runs of reference/results/controls. A run entry is (method, estimator,
# level, dt) or (method, estimator, level, dt, end_time, wall_y)
CONTROL_CUTLINK = (("cutLink", "ME"), ("cutLink", "FO"))
CONTROL_STEP_CASES = {
    # step: {level: variants}; the L1 steps of dt-L1 and the 1.25e-4 s step
    # of Levels 2-4 are all variants, the finer steps those tracked in
    # reference/results/controls/cases-timestep.csv (ST is a column of ME)
    0.000125: {lv: RUN_VARIANTS for lv in LEVELS},
    0.0000625: {3: (("penalty", ""), ("incremental", ""), ("ghostCell", "")),
                4: RUN_VARIANTS},
    0.00003125: {3: (("incremental", ""),), 4: (("incremental", ""),)},
}
WALL_FACE = 0.125                   # a mesh face on Levels 1-4 (0.5 m / 2^k)
T_LONG = 8.0                        # the eight-period control
RUNS = {
    "smoke": ("cutLink (momentum exchange and surface traction), Level 1",
              [("cutLink", "ME", 1, 0.0)]),
    "L1": ("all variants, Level 1", [(m, e, 1, 0.0) for m, e in RUN_VARIANTS]),
    "L2": ("all variants, Level 2", [(m, e, 2, 0.0) for m, e in RUN_VARIANTS]),
    "L3": ("all variants, Level 3", [(m, e, 3, 0.0) for m, e in RUN_VARIANTS]),
    "L4": ("all variants, Level 4", [(m, e, 4, 0.0) for m, e in RUN_VARIANTS]),
    "matrix": ("all variants, Levels 1-4 (the tracked matrix)",
               [(m, e, lv, 0.0) for lv in LEVELS for m, e in RUN_VARIANTS]),
    "dt": ("time-step control: deltaT 0.0025, 0.0005, 0.00025 s on Levels "
           "2-4, four variants",
           [(m, e, lv, dt) for lv in DT_LEVELS for dt in DT_STEPS
            for m, e in DT_VARIANTS]),
    "dt-L1": ("control: deltaT 0.0025, 0.0005, 0.00025 s on Level 1, all "
              "variants",
              [(m, e, 1, dt) for dt in DT_STEPS for m, e in RUN_VARIANTS]),
    "dt-fine": ("control: deltaT 1.25e-4 s, Levels 1-4, all variants",
                [(m, e, lv, 0.000125) for lv, vs in
                 CONTROL_STEP_CASES[0.000125].items() for m, e in vs]),
    "dt-6e-5": ("control: deltaT 6.25e-5 s (Level 3: penalty, incremental, "
                "ghostCell; Level 4: all variants)",
                [(m, e, lv, 0.0000625) for lv, vs in
                 CONTROL_STEP_CASES[0.0000625].items() for m, e in vs]),
    "dt-3e-5": ("control: deltaT 3.125e-5 s, incremental, Levels 3-4",
                [(m, e, lv, 0.00003125) for lv, vs in
                 CONTROL_STEP_CASES[0.00003125].items() for m, e in vs]),
    "wall-face": ("control: the wall on the mesh face y = 0.125 m, cutLink "
                  "(momentum exchange, forcing), Levels 2-4",
                  [(m, e, lv, 0.0, 0.0, WALL_FACE) for lv in (2, 3, 4)
                   for m, e in CONTROL_CUTLINK]),
    "T8": ("control: eight periods, cutLink momentum exchange and penalty, "
           "Level 3",
           [("cutLink", "ME", 3, 0.0, T_LONG, 0.0),
            ("penalty", "", 3, 0.0, T_LONG, 0.0)]),
}
RUNS["controls"] = ("all control runs of reference/results/controls",
                    [i for n in ("dt-L1", "dt-fine", "dt-6e-5", "dt-3e-5",
                                 "wall-face", "T8") for i in RUNS[n][1]])
DEFAULT_RUNS = ("smoke",)
LEVEL_RANKS = {1: 1, 2: 1, 3: 1, 4: 1}

FOAM_HEADER = (
    "/*--------------------------------*- C++ -*----------------------------"
    "------*\\\n| solids4foam: solid mechanics and fluid-solid interaction "
    "simulations        |\n\\*-----------------------------------------------"
    "----------------------------*/\nFoamFile\n{{\n    version     2.0;\n"
    "    format      ascii;\n    class       dictionary;\n"
    "    location    \"{location}\";\n    object      {object};\n}}\n"
    "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * "
    "* * * * //\n\n")


def fail(message: str) -> None:
    print(f"ERROR: {message}", file=sys.stderr)
    sys.exit(2)


def fmt(value: float) -> str:
    return f"{value:.10g}"


# ------------------------------------------------------------------ parameters
def dict_scalar(text: str, key: str) -> float:
    """First numerical entry 'key [key] [dimensions] number;' of an OpenFOAM
    dictionary (the name may be repeated, as in transportProperties)"""
    m = re.search(rf"^\s*{re.escape(key)}\s+(?:{re.escape(key)}\s+)?"
                  r"(?:\[[^\]]*\]\s*)?([-+0-9.eE]+)\s*;", text, re.MULTILINE)
    if not m:
        fail(f"entry '{key}' not found")
    return float(m.group(1))


def stl_top(path: Path) -> float:
    """Largest y of the vertices of the slab surface (its wall), rounded to
    1e-6 m: the vertices are single-precision numbers"""
    ys = [float(m) for m in re.findall(r"vertex\s+\S+\s+(\S+)\s+\S+",
                                       path.read_text())]
    if not ys:
        fail(f"no vertices in {path}")
    return round(max(ys), 6)


def parameters(tutorial: Path = TUTORIAL_DIR) -> dict:
    """The physical parameters of the tutorial, from its files: nu and rho
    from constant/transportProperties; the amplitude, period, phase and
    direction of the sinusoidalTranslation (position A sin(omega t + phase))
    from constant/fvOptions, hence U0 = A omega; the wall y_w from the slab
    surface; the mesh height, the cell-count rule, the cyclic length rule and
    the depth from system/blockMeshDict; the time step from the controlDict"""
    fv = (tutorial / "constant/fvOptions").read_text()
    tp = (tutorial / "constant/transportProperties").read_text()
    bm = (tutorial / "system/blockMeshDict").read_text()
    cd = (tutorial / "system/controlDict").read_text()
    motion = re.search(r"motion\s*\{(.*?)\}", fv, re.DOTALL).group(1)
    amplitude = dict_scalar(motion, "amplitude")
    period = dict_scalar(motion, "period")
    phase = dict_scalar(motion, "phase")
    direction = tuple(float(x) for x in re.search(
        r"direction\s+\(([^)]*)\)", motion).group(1).split())
    omega = 2.0 * math.pi / period
    n0 = int(re.search(r"round\((\d+)\*pow", bm).group(1))
    div = int(re.search(r"round\(\$ny/(\d+)\)", bm).group(1))
    height = float(re.search(r"lx\s+#eval\s+\"([0-9.]+)\*", bm).group(1))
    vertices = re.search(r"vertices\s*\((.*?)\n\);", bm, re.DOTALL).group(1)
    zs = sorted({float(m.split()[2]) for m in
                 re.findall(r"\(([-0-9.$a-z ]+)\)", vertices)
                 if len(m.split()) == 3 and not m.split()[2].startswith("$")})
    nu = dict_scalar(tp, "nu")
    return {"nu": nu, "rho": dict_scalar(tp, "rho"), "amplitude": amplitude,
            "period": period, "phase": phase, "direction": direction,
            "omega": omega, "U0": amplitude * omega,
            "yw": stl_top(tutorial / "constant/triSurface/slab.stl"),
            "height": height, "n0": n0, "div": div,
            "depth": zs[-1] - zs[0], "delta": math.sqrt(2.0 * nu / omega),
            "dt": dict_scalar(cd, "deltaT")}


def mesh_shape(p: dict, level: int) -> tuple:
    """ny, nx, cyclic length lx, wall area lx*depth of a mesh level"""
    ny = round(p["n0"] * 2 ** (level - 1))
    nx = round(ny / p["div"])
    lx = p["height"] * nx / ny
    return ny, nx, lx, lx * p["depth"]


# ------------------------------------------------------------ exact solutions
def slab_velocity(p: dict, t: float) -> float:
    """velocity of the slab: d/dt of A sin(omega t + phase)"""
    return p["amplitude"] * p["omega"] * math.cos(p["omega"] * t + p["phase"])


def exact_u(p: dict, s: float, t: float) -> float:
    """periodic half-space solution"""
    k = 1.0 / p["delta"]
    return p["U0"] * math.exp(-k * s) * math.sin(p["omega"] * t - k * s)


def exact_tau(p: dict, t: float) -> float:
    """shear stress of the fluid on the wall, rho nu du/ds(0) (Pa)"""
    k = 1.0 / p["delta"]
    return -p["rho"] * p["nu"] * k * p["U0"] * (
        math.sin(p["omega"] * t) + math.cos(p["omega"] * t))


def exact_amp_phase(p: dict) -> tuple:
    """tau = A sin(omega t + phi): A = rho nu k U0 sqrt(2), phi = 5 pi/4"""
    k = 1.0 / p["delta"]
    return (p["rho"] * p["nu"] * k * p["U0"] * math.sqrt(2.0),
            math.pi / 4.0 + math.pi)


def exact_u_finite(p: dict, s: float, t: float, H: float) -> float:
    """periodic solution with a slip lid at s = H"""
    lam = (1.0 + 1.0j) / p["delta"]
    f = cmath.cosh(lam * (H - s)) / cmath.cosh(lam * H)
    return (p["U0"] * cmath.exp(1j * p["omega"] * t) * f).imag


def exact_tau_finite(p: dict, t: float, H: float) -> float:
    lam = (1.0 + 1.0j) / p["delta"]
    return (p["rho"] * p["nu"] * p["U0"] * cmath.exp(1j * p["omega"] * t)
            * (-lam) * cmath.tanh(lam * H)).imag


def ivp_modes(p: dict, H: float, n: int = N_MODES,
              nint: int = N_SIMPSON) -> list:
    """(k_n, a_n) of the start-up transient from rest"""
    h = H / nint
    f0 = [exact_u_finite(p, i * h, 0.0, H) for i in range(nint + 1)]
    out = []
    for m in range(n):
        k = (m + 0.5) * math.pi / H
        acc = 0.0
        for i in range(nint + 1):
            wgt = 1 if i in (0, nint) else (4 if i % 2 else 2)
            acc += wgt * f0[i] * math.sin(k * i * h)
        out.append((k, -(2.0 / H) * acc * h / 3.0))
    return out


def exact_u_ivp(p: dict, s: float, t: float, H: float, modes: list) -> float:
    return exact_u_finite(p, s, t, H) + sum(
        a * math.exp(-p["nu"] * k * k * t) * math.sin(k * s)
        for k, a in modes)


def exact_tau_ivp(p: dict, t: float, H: float, modes: list) -> float:
    return exact_tau_finite(p, t, H) + p["rho"] * p["nu"] * sum(
        a * math.exp(-p["nu"] * k * k * t) * k for k, a in modes)


_MODES = {}


def cached_modes(p: dict) -> list:
    H = p["height"] - p["yw"]
    key = (p["nu"], p["omega"], p["U0"], H)
    if key not in _MODES:
        _MODES[key] = ivp_modes(p, H)
    return _MODES[key]


# --------------------------------------------------------- files and fields
def read_force_file(path: Path) -> dict:
    """Columns of postProcessing/immersedBoundary/0/<body>.dat by name"""
    header, rows = None, []
    for line in path.read_text().splitlines():
        if line.startswith("# Time"):
            header = line[2:].split("\t")
        elif line and not line.startswith("#"):
            rows.append([float(x) for x in line.split()])
    if header is None or not rows:
        fail(f"no data in {path}")
    out, col = {"t": [r[0] for r in rows]}, 1
    for name in header[1:]:
        key = ("force" if name.startswith("force (") else
               "torque" if name.startswith("torque") else
               "inertia" if name.startswith("inertia") else
               "momentum" if "momentum exchange" in name else "traction")
        out[key] = [r[col:col + 3] for r in rows]
        col += 3
    return out


def read_vector_field(path: Path) -> list:
    text = path.read_text()
    m = re.search(r"internalField\s+nonuniform\s+List<vector>\s*\n?(\d+)"
                  r"\s*\n?\(", text)
    if not m:
        fail(f"{path}: not a nonuniform vector field")
    body = text[m.end():]
    body = body[:body.index("\n)\n")]
    return [tuple(float(x) for x in line.strip("()").split())
            for line in body.strip().splitlines()]


def shear_series(columns: dict, method: str, estimator: str,
                 area: float) -> tuple:
    """Time and wall shear stress (force/area) of an estimator. cutLink: ME
    and FO are the selected estimator ('force' column), ST is the extra
    'traction' column (writeSurfaceTraction), and in an FO run the ME and ST
    columns are written after the inertia. Others: forcing plus the inertia
    of the fluid inside the slab"""
    if method == "cutLink":
        if estimator == "ST":
            key = "traction"
        else:
            key = "force"
        f = [r[0] for r in columns[key]]
    else:
        f = [a[0] + b[0] for a, b in zip(columns["force"], columns["inertia"])]
    return columns["t"], [x / area for x in f]


def aligned_shear_series(columns: dict, area: float) -> tuple:
    """Diagnostic for penalty, incremental and ghostCell: the inertia column
    is a backward difference, (M_n - M_(n-1))/dt, centred at t_n - dt/2,
    while the forcing is attributed to t_n. Here the inertia is averaged over
    two successive steps, which centres it at t_n (the last row keeps its
    own value). The solver output is unchanged; this is a post-processing
    variant reported next to the tutorial's (columns with the suffix _al)"""
    inertia = [b[0] for b in columns["inertia"]]
    centred = [0.5 * (inertia[i] + inertia[min(i + 1, len(inertia) - 1)])
               for i in range(len(inertia))]
    f = [a[0] + c for a, c in zip(columns["force"], centred)]
    return columns["t"], [x / area for x in f]


# ---------------------------------------------------------------- statistics
def read_scalar_field(path: Path) -> list:
    text = path.read_text()
    m = re.search(r"internalField\s+nonuniform\s+List<scalar>\s*\n?(\d+)"
                  r"\s*\n?\(", text)
    if not m:
        fail(f"{path}: not a nonuniform scalar field")
    body = text[m.end():]
    body = body[:body.index("\n)\n")]
    return [float(x) for x in body.strip().splitlines()]


def covered_volume_ratio(p: dict, level: int, occupancy: list) -> float:
    """Volume of the body in the mesh, sum(lambda V), over the exact volume
    y_w * lx * depth (equal cell volumes). It sets the inertia of the fluid
    inside the slab, the larger term of the force of penalty, incremental
    and ghostCell"""
    ny, nx, _, area = mesh_shape(p, level)
    return sum(occupancy) * (p["height"] / ny) * (area / nx) / (p["yw"] * area)


def solve_linear(a: list, b: list) -> list:
    """Gaussian elimination with partial pivoting"""
    n = len(b)
    m = [row[:] + [b[i]] for i, row in enumerate(a)]
    for c in range(n):
        piv = max(range(c, n), key=lambda r: abs(m[r][c]))
        m[c], m[piv] = m[piv], m[c]
        for r in range(c + 1, n):
            f = m[r][c] / m[c][c]
            for k in range(c, n + 1):
                m[r][k] -= f * m[c][k]
    x = [0.0] * n
    for r in range(n - 1, -1, -1):
        x[r] = (m[r][n] - sum(m[r][k] * x[k] for k in range(r + 1, n))) / m[r][r]
    return x


def harmonic_fit(t: list, y: list, omega: float, nh: int) -> tuple:
    """Least squares of mean + nh harmonics (columns 1, sin k w t, cos k w t)
    by the normal equations; returns the coefficients and the fit"""
    def row(x):
        out = [1.0]
        for k in range(1, nh + 1):
            out += [math.sin(k * omega * x), math.cos(k * omega * x)]
        return out
    rows = [row(x) for x in t]
    n = len(rows[0])
    ata = [[sum(r[i] * r[j] for r in rows) for j in range(n)] for i in range(n)]
    aty = [sum(r[i] * v for r, v in zip(rows, y)) for i in range(n)]
    c = solve_linear(ata, aty)
    return c, [sum(ci * ri for ci, ri in zip(c, r)) for r in rows]


def rms(values: list) -> float:
    return math.sqrt(sum(v * v for v in values) / len(values))


def shear_metrics(p: dict, t: list, tau: list, period=PERIOD) -> dict:
    a, b = period
    idx = [i for i, x in enumerate(t) if a < x <= b + 1e-12]
    tt = [t[i] for i in idx]
    yy = [tau[i] for i in idx]
    H = p["height"] - p["yw"]
    aex, _ = exact_amp_phase(p)
    out = {}
    for tag, fn in (("", lambda x: exact_tau(p, x)),
                    ("_ivp", lambda x: exact_tau_ivp(p, x, H, cached_modes(p)))):
        ref = [fn(x) for x in tt]
        rms_ref = rms(ref)
        cn, fit = harmonic_fit(tt, yy, p["omega"], NHARM)
        cr, _ = harmonic_fit(tt, ref, p["omega"], NHARM)
        an, ar = math.hypot(cn[1], cn[2]), math.hypot(cr[1], cr[2])
        pn, pr = math.atan2(cn[2], cn[1]), math.atan2(cr[2], cr[1])
        dphi = (pn - pr + math.pi) % (2 * math.pi) - math.pi
        out["shear_rms_err" + tag] = rms([x - r for x, r in zip(yy, ref)]) / rms_ref
        out["shear_amp_err" + tag] = (an - ar) / ar
        out["shear_phase_err_deg" + tag] = math.degrees(dphi)
        if not tag:
            out.update(shear_mean=cn[0], shear_amp=an, shear_amp_exact=aex,
                       shear_noise=rms([x - f for x, f in zip(yy, fit)]) / rms_ref)
    return out


def snapshot_times(period=PERIOD) -> list:
    return [period[0] + 0.125 * k for k in range(1, N_SNAPSHOTS + 1)]


def row_profile(p: dict, level: int, field: list) -> list:
    """u_x averaged over the cyclic direction for each row of cells"""
    ny, nx, _, _ = mesh_shape(p, level)
    if len(field) != nx * ny:
        fail(f"field has {len(field)} cells, expected {nx * ny}")
    return [sum(field[j * nx + i][0] for i in range(nx)) / nx
            for j in range(ny)]


def cell_heights(p: dict, level: int) -> tuple:
    ny, _, _, _ = mesh_shape(p, level)
    h = p["height"] / ny
    return h, [(j + 0.5) * h - p["yw"] for j in range(ny)]


def velocity_norms(p: dict, level: int, profiles: dict) -> dict:
    """Norms against the periodic and the initial-value solution from the
    row profiles {time: [u_x per row]}. Equal cell volumes: the
    volume-weighted norms are plain means over the fluid rows"""
    h, s = cell_heights(p, level)
    H = p["height"] - p["yw"]
    modes = cached_modes(p)
    sets = {"all": [j for j, x in enumerate(s) if x > 0],
            "near": [j for j, x in enumerate(s) if 0 < x <= p["delta"]],
            "iface": [j for j, x in enumerate(s) if 0 < x < h]}
    out = {}
    for tag, fn in (("", lambda si, tk: exact_u(p, si, tk)),
                    ("_ivp", lambda si, tk: exact_u_ivp(p, si, tk, H, modes))):
        acc = {k: [] for k in sets}
        for tk, ux in profiles.items():
            err = [ux[j] - (fn(s[j], tk) if s[j] > 0 else 0.0)
                   for j in range(len(s))]
            for k, rows in sets.items():
                e = [abs(err[j]) for j in rows]
                acc[k].append((sum(e) / len(e) / p["U0"],
                               rms(e) / p["U0"], max(e) / p["U0"]))
        for k in sets:
            out[f"u_L1_{k}{tag}"] = rms([a[0] for a in acc[k]])
            out[f"u_L2_{k}{tag}"] = rms([a[1] for a in acc[k]])
            out[f"u_Linf_{k}{tag}"] = max(a[2] for a in acc[k])
    return out


# ------------------------------------------------------------ case handling
def dt_tag(dt: float) -> str:
    return f"-dt{fmt(dt)}" if dt else ""


def case_name(method: str, est: str, level: int, dt: float = 0.0,
              end_time: float = 0.0, wall_y: float = 0.0,
              ranks_tag: str = "") -> str:
    """Tutorial cases have no tag; a changed step, end time, wall position
    or rank count is tagged -dt..., -T..., -yw..., -np..."""
    return (f"{method}{'-' + est if method == 'cutLink' else ''}-L{level}"
            f"{dt_tag(dt)}{'-T' + fmt(end_time) if end_time else ''}"
            f"{'-yw' + fmt(wall_y) if wall_y else ''}{ranks_tag}")


def load_manifest() -> dict:
    if MANIFEST.exists():
        return json.loads(MANIFEST.read_text())
    return {"cases": {}}


def save_manifest(manifest: dict) -> None:
    WORK_DIR.mkdir(parents=True, exist_ok=True)
    MANIFEST.write_text(json.dumps(manifest, indent=1, sort_keys=True) + "\n")


def plan(args) -> list:
    """The cases (method, estimator, level, name, ranks, dt, end_time,
    wall_y) of the named runs; --end-time and --wall-y override the options
    of every case, a case with a changed step, end time, wall position or
    rank count carries it in its name"""
    out, seen = [], set()
    end_arg = getattr(args, "end_time", 0.0) or 0.0
    wall_arg = getattr(args, "wall_y", 0.0) or 0.0
    for name in args.runs:
        if name not in RUNS:
            fail(f"unknown run '{name}'; choose from {', '.join(RUNS)}")
        for method, est, level, dt, *opt in RUNS[name][1]:
            if args.methods and method not in args.methods:
                continue
            if args.levels and level not in args.levels:
                continue
            end_time = end_arg or (opt[0] if opt else 0.0)
            wall_y = wall_arg or (opt[1] if len(opt) > 1 else 0.0)
            ranks = args.ranks or LEVEL_RANKS[level]
            key = (method, est, level, dt, end_time, wall_y, ranks)
            if key not in seen:
                seen.add(key)
                tag = f"-np{ranks}" if ranks != LEVEL_RANKS[level] else ""
                out.append((method, est, level,
                            case_name(method, est, level, dt, end_time,
                                      wall_y, tag), ranks, dt, end_time,
                            wall_y))
    return out


def move_wall(path: Path, wall_y: float) -> None:
    """Move the top vertices of the slab surface (its wall) to y = wall_y,
    leaving the other vertices (the control of the wall on a mesh face)"""
    text = path.read_text()
    top = stl_top(path)
    if not 0.0 < wall_y < parameters()["height"]:
        fail(f"--wall-y {wall_y} is outside the mesh")
    count = [0]

    def move(m):
        if abs(float(m.group(2)) - top) < 1e-6:
            count[0] += 1
            return f"{m.group(1)}{wall_y:.6f}{m.group(3)}"
        return m.group(0)
    new = re.sub(r"(vertex\s+\S+\s+)(\S+)(\s+\S+)", move, text)
    if not count[0]:
        fail(f"no top vertex found in {path}")
    path.write_text(new)


def replace_entry(path: Path, key: str, value: str) -> None:
    text = path.read_text()
    pattern = re.compile(rf"^(\s*{re.escape(key)}\s+)[^;\n]*;", re.MULTILINE)
    new, count = pattern.subn(lambda m: f"{m.group(1)}{value};", text)
    if count != 1:
        fail(f"found {count} '{key}' entries in {path}")
    path.write_text(new)


def write_fv_options(path: Path, method: str, est: str) -> None:
    """The tutorial's fvOptions with the method settings of METHOD_SETTINGS;
    the motion block is never rewritten. cutLink keeps
    'writeSurfaceTraction yes', so that the momentum exchange and the surface
    traction are two estimators of one run (the forcing estimator is a
    separate run); the other methods do not write it"""
    text = path.read_text()
    option, body = METHOD_SETTINGS[method]
    lines = [("method", "cutLink")] if method == "cutLink" else list(option)
    if method == "cutLink" and est != "ME":
        lines.append(("forceEstimator", ESTIMATOR_NAMES[est]))
    block = "".join(f"    {k.ljust(24)}{v};\n" for k, v in lines)
    new, count = re.subn(r"^    method\s+cutLink;\n", lambda m: block, text,
                         flags=re.MULTILINE)
    if count != 1:
        fail("the method entry of the tutorial fvOptions was not found "
             "exactly once")
    if method != "cutLink":
        new, count = re.subn(r"^    writeSurfaceTraction\s+yes;\n", "", new,
                             flags=re.MULTILINE)
        if count != 1:
            fail("writeSurfaceTraction entry not found")
    if body:
        extra = "".join(f"            {k.ljust(12)}{v};\n" for k, v in body)
        new, count = re.subn(r"(\n\s*surface\s+\"[a-z]+\.stl\";\n)",
                             lambda m: m.group(1) + extra, new)
        if count != 1:
            fail("expected the surface entry of the slab")
    path.write_text(new)


def write_decompose_par_dict(path: Path, ranks: int) -> None:
    """Slabs along y (the cyclic direction x has only a few cells)"""
    text = FOAM_HEADER.format(location="system", object="decomposeParDict")
    text += (f"numberOfSubdomains  {ranks};\n\nmethod          simple;\n\n"
             f"simpleCoeffs\n{{\n    n           (1 {ranks} 1);\n}}\n\n"
             "// **************************************************************"
             "*********** //\n")
    path.write_text(text)


RUN_SCRIPT = """\
#!/bin/bash
# Generated by verification/scripts/oscillating_wall_stokes_layer_verification.py
cd "${{0%/*}}" || exit 1
export OMP_NUM_THREADS=1
# Open MPI binds the ranks of every mpirun to the first cores by default; two
# simultaneous runs then share the same cores. No binding is requested here
export OMPI_MCA_hwloc_base_binding_policy=none
export MESH_LEVEL={level}
start=$(date +%s)
{body}
status=$?
echo "wall_s=$(( $(date +%s) - start )) ranks={ranks} allrun_status=$status" \\
    >> wallclock.txt
exit $status
"""

SERIAL_BODY = """\
(
    . "$WM_PROJECT_DIR/bin/tools/RunFunctions"
    source solids4FoamScripts.sh
    solids4Foam::caseDoesNotRunWithFoamExtend
    solids4Foam::caseDoesNotRunWithOpenFOAMOrg
    solids4Foam::runApplication blockMesh || exit 1
    solids4Foam::runApplication solids4Foam
)"""

PARALLEL_BODY = """\
(
    . "$WM_PROJECT_DIR/bin/tools/RunFunctions"
    source solids4FoamScripts.sh
    solids4Foam::caseDoesNotRunWithFoamExtend
    solids4Foam::caseDoesNotRunWithOpenFOAMOrg
    solids4Foam::runApplication blockMesh || exit 1
    solids4Foam::runApplication decomposePar || exit 1
    solids4Foam::runParallel solids4Foam || exit 1
    solids4Foam::runApplication reconstructPar -newTimes
)"""


def input_hashes(case: Path) -> dict:
    files = ["constant/fvOptions", "system/controlDict", "system/fvSolution",
             "system/fvSchemes", "system/blockMeshDict",
             "system/decomposeParDict", "0/U", "0/p"]
    return {f: hashlib.sha256((case / f).read_bytes()).hexdigest()
            for f in files if (case / f).exists()}


def generate_case(method: str, est: str, level: int, name: str, ranks: int,
                  force: bool = False, cases_dir: Path = None,
                  dt: float = 0.0, end_time: float = 0.0,
                  wall_y: float = 0.0) -> Path:
    cases_dir = CASES_DIR if cases_dir is None else cases_dir
    case = cases_dir / name
    if case.exists():
        if not force:
            fail(f"{case} exists; use --force to regenerate it")
        shutil.rmtree(case)
    skip = shutil.ignore_patterns("verification", "regressionTests", "log.*",
                                  "postProcessing", "processor*", "polyMesh",
                                  "[1-9]*", "0.[0-9]*", "*.pdf")
    shutil.copytree(TUTORIAL_DIR, case, ignore=skip)
    write_fv_options(case / "constant/fvOptions", method, est)
    replace_entry(case / "system/fvSolution", "nOuterCorrectors",
                  str(N_OUTER[method]))
    if dt:
        control = case / "system/controlDict"
        if not re.search(r"^adjustTimeStep\s+no;", control.read_text(),
                         re.MULTILINE):
            fail("the tutorial no longer has a fixed time step")
        replace_entry(control, "deltaT", fmt(dt))
    if end_time:
        replace_entry(case / "system/controlDict", "endTime", fmt(end_time))
    if wall_y:
        move_wall(case / "constant/triSurface/slab.stl", wall_y)
    write_decompose_par_dict(case / "system/decomposeParDict", ranks)
    script = case / "Allrun.verify"
    script.write_text(RUN_SCRIPT.format(
        level=level, ranks=ranks,
        body=PARALLEL_BODY if ranks > 1 else SERIAL_BODY))
    script.chmod(0o755)
    return case


def manifest_entry(method, est, level, dt, ranks, end_time, wall_y,
                   case: Path) -> dict:
    return {"method": method, "estimator": est, "level": level, "dt": dt,
            "end_time": end_time, "wall_y": wall_y, "ranks": ranks,
            "inputs_sha256": input_hashes(case),
            "generated": datetime.now().isoformat(timespec="seconds")}


def cmd_generate(args) -> int:
    manifest = load_manifest()
    for method, est, level, name, ranks, dt, end_time, wall_y in plan(args):
        case = generate_case(method, est, level, name, ranks, args.force,
                             dt=dt, end_time=end_time, wall_y=wall_y)
        manifest["cases"][name] = manifest_entry(
            method, est, level, dt, ranks, end_time, wall_y, case)
        print(f"{name}: generated in {case}")
    save_manifest(manifest)
    return 0


def cmd_dry_run(args) -> int:
    p = parameters()
    print(f"Tutorial: {TUTORIAL_DIR.name}; U0 = {p['U0']:.6g}, omega = "
          f"{p['omega']:.6g}, nu = {p['nu']:.8g}, delta = {p['delta']:.6g}, "
          f"y_w = {p['yw']}, dt = {p['dt']}")
    a, _ = exact_amp_phase(p)
    print(f"Exact wall shear amplitude: {a:.6e} Pa (unit density)")
    total = 0.0
    observed = observed_rank_hours()
    for method, est, level, name, ranks, dt, end_time, wall_y in plan(args):
        ny, nx, _, area = mesh_shape(p, level)
        obs = (None if end_time or wall_y else observed.get(
            (method, est, level, dt or TUTORIAL_DT)))
        total += obs or 0.0
        print(f"  {name:20s} {nx * ny:6d} cells, {ranks} rank(s), observed "
              f"{'%.4f' % obs if obs else 'n/a'} allocated rank-h")
    print(f"Observed total of the listed cases: {total:.3f} rank-h "
          "(single runs on a shared node, not a performance model)")
    return 0


def check_environment(ranks: int) -> str:
    """Empty if the cases can be run, else the reason they cannot. The rule
    of solids4FoamScripts.sh: OpenFOAM.com sets WM_PROJECT=OpenFOAM and a
    WM_PROJECT_VERSION that starts with 'v' (v2412); OpenFOAM.org sets the
    same WM_PROJECT with a number (11, 12)"""
    if "WM_PROJECT_DIR" not in os.environ:
        return "no OpenFOAM environment (WM_PROJECT_DIR is not set)"
    if os.environ.get("WM_PROJECT") != "OpenFOAM" or \
            not os.environ.get("WM_PROJECT_VERSION", "").startswith("v"):
        return "OpenFOAM.com (ESI) is required for libimmersedBoundary"
    for tool in ("blockMesh", "solids4Foam"):
        if not shutil.which(tool):
            return f"{tool} not found in PATH"
    if ranks > 1 and not shutil.which("mpirun"):
        return "mpirun not found"
    return ""


def solver_problem(case: Path) -> str:
    log = case / "log.solids4Foam"
    if not log.exists():
        return "no solver log"
    # the banner line "trapFpe: ... (FOAM_SIGFPE)" is not an error
    text = "\n".join(line for line in log.read_text(errors="replace").splitlines()
                     if "trapFpe" not in line and "FOAM_SIGFPE" not in line)
    if re.search(r"FOAM FATAL|sigFpe|\bnan\b", text, re.IGNORECASE):
        return "the solver log contains an error or a NaN"
    if not re.search(r"^End\s*$", text, re.MULTILINE):
        return "the solver did not finish"
    return ""


def cmd_run(args) -> int:
    manifest = load_manifest()
    todo = []
    for item in plan(args):
        entry = manifest["cases"].get(item[3], {})
        if entry.get("state") == "completed" and \
                (CASES_DIR / item[3]).exists() and not args.force:
            print(f"{item[3]}: completed, kept (use --force to rerun)")
        else:
            todo.append(item)
    for method, est, level, name, ranks, dt, end_time, wall_y in todo:
        case = generate_case(method, est, level, name, ranks, True, dt=dt,
                             end_time=end_time, wall_y=wall_y)
        manifest["cases"][name] = manifest_entry(
            method, est, level, dt, ranks, end_time, wall_y, case)
        print(f"{name}: generated in {case}")
    save_manifest(manifest)
    status = 0
    for method, est, level, name, ranks, dt, end_time, wall_y in todo:
        reason = check_environment(ranks)
        entry = manifest["cases"][name]
        if reason:
            print(f"SKIPPED {name}: {reason}")
            entry["state"] = "skipped"
            continue
        case = CASES_DIR / name
        print(f"Running {name} ({ranks} rank(s)) ...", flush=True)
        rc = subprocess.call([str(case / "Allrun.verify")])
        problem = solver_problem(case) if rc == 0 else f"Allrun status {rc}"
        if rc == 0 and problem == "no solver log":
            # Allrun.verify exits 0 without a solver log only through the
            # guards of solids4FoamScripts.sh (caseDoesNotRunWith...): the
            # case was not run, which is SKIPPED, never PASS and never FAIL
            print(f"SKIPPED {name}: Allrun.verify exited 0 without running "
                  "the solver (an OpenFOAM or solids4foam guard)")
            entry["state"] = "skipped"
            continue
        entry["state"] = "completed" if not problem else "failed"
        if problem:
            print(f"FAILED {name}: {problem}")
            status = 1
    save_manifest(manifest)
    return status


def cmd_status(args) -> int:
    for name, entry in sorted(load_manifest()["cases"].items()):
        print(f"{name:22s} {entry.get('state', 'generated')}")
    return 0


# ------------------------------------------------------------------ analysis
def wall_seconds(case: Path) -> float:
    f = case / "wallclock.txt"
    if not f.exists():
        return float("nan")
    m = re.findall(r"wall_s=(\d+)", f.read_text())
    return float(m[-1]) if m else float("nan")


def analysis_windows(end_time: float = 0.0) -> list:
    """The periods analysed: the fourth, and for an end-time control also
    the last (end_time - 1 < t <= end_time)"""
    if end_time and end_time > PERIOD[1] + 1e-9:
        return [PERIOD, (end_time - 1.0, end_time)]
    return [PERIOD]


def window_label(period: tuple) -> str:
    return f"{period[0]:g}-{period[1]:g}"


def analyse_case(case: Path, method: str, est: str, level: int,
                 ranks: int = 1, wall_y: float = 0.0,
                 end_time: float = 0.0) -> tuple:
    """The result rows of a completed case (two for the cutLink ME run: ME
    and ST) and the shear series {estimator: (t, tau)}. wall_y is the wall
    of a case with the slab moved (the exact solutions use it); a case that
    ran beyond four periods has the rows of the fourth and of the last
    period, tagged 'window'"""
    p = parameters(TUTORIAL_DIR)
    if wall_y:
        p["yw"] = wall_y
    ny, nx, _, area = mesh_shape(p, level)
    cols = read_force_file(case / "postProcessing/immersedBoundary/0/slab.dat")
    dt_text = (case / "system/controlDict").read_text()
    labels = ["ME", "ST"] if (method == "cutLink" and est == "ME") else \
        [est if method == "cutLink" else ""]
    windows = analysis_windows(end_time)
    rows, series = [], {}
    for period in windows:
        times = snapshot_times(period)
        profiles = {tk: row_profile(
            p, level, read_vector_field(case / f"{tk:g}" / "U"))
            for tk in times}
        vel = velocity_norms(p, level, profiles)
        lam_file = case / f"{times[-1]:g}" / "immersedBoundary:lambda"
        covered = (covered_volume_ratio(p, level, read_scalar_field(lam_file))
                   if lam_file.exists() else "")
        for lab in labels:
            t, tau = shear_series(cols, method, lab, area)
            row = {"method": method, "estimator": lab, "level": level,
                   "cells": nx * ny, "ranks": ranks,
                   "dt": dict_scalar(dt_text, "deltaT"), "steps": len(t),
                   "wall_s": wall_seconds(case), "covered_volume": covered}
            row.update(shear_metrics(p, t, tau, period))
            if method != "cutLink":
                ta, taua = aligned_shear_series(cols, area)
                al = shear_metrics(p, ta, taua, period)
                for key in ("shear_amp_err_ivp", "shear_phase_err_deg_ivp",
                            "shear_rms_err_ivp"):
                    row[key[:-4] + "_al"] = al[key]
            row.update(vel)
            if wall_y or len(windows) > 1:
                row["window"] = window_label(period)
                row["wall_y"] = wall_y or p["yw"]
            rows.append(row)
            series[lab] = (t, tau)
    return rows, series


def observed_rank_hours() -> dict:
    """Allocated rank-hours of the tracked runs, {(method, estimator, level,
    dt): rank_h}, from the tracked rows (reference/results/cases.csv and
    controls/cases-timestep.csv)"""
    out = {}
    for path in (CASES_CSV, CONTROLS_CSV):
        for r in read_results(path):
            try:
                out[row_key(r)] = float(r["rank_h"])
            except (KeyError, ValueError):
                pass
    return out


def row_dt(row: dict) -> float:
    """The time step of a row; blank is the tutorial step"""
    dt = row.get("dt", "")
    return float(dt) if dt not in ("", None) else TUTORIAL_DT


def row_key(row: dict) -> tuple:
    return (row["method"], row["estimator"], int(row["level"]), row_dt(row))


def sort_key(r: dict) -> tuple:
    dt = row_dt(r)
    return (dt != TUTORIAL_DT, -dt, variant_index(r), int(r["level"]))


FIELDS = (["method", "estimator", "level", "cells", "ranks", "dt", "steps",
           "shear_amp", "shear_amp_exact", "shear_amp_err",
           "shear_phase_err_deg", "shear_rms_err", "shear_mean",
           "shear_noise", "shear_amp_err_ivp", "shear_phase_err_deg_ivp",
           "shear_rms_err_ivp", "shear_amp_err_al", "shear_phase_err_deg_al",
           "shear_rms_err_al"]
          + [f"u_{n}_{s}" for s in ("all", "near", "iface")
             for n in ("L1", "L2", "Linf")]
          + [f"u_{n}_all_ivp" for n in ("L1", "L2", "Linf")]
          + ["covered_volume", "wall_s", "rank_h"])


def variant_index(r: dict) -> int:
    return VARIANTS.index((r["method"], r["estimator"]))


def write_results(path: Path, rows: list, fields=None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fields or FIELDS, lineterminator="\n",
                           extrasaction="ignore", restval="")
        w.writeheader()
        for r in sorted(rows, key=sort_key):
            w.writerow({k: (f"{v:.6g}" if isinstance(v, float) else v)
                        for k, v in r.items()})


def read_results(path: Path = CASES_CSV) -> list:
    if not path.exists():
        return []
    with path.open(newline="") as f:
        rows = list(csv.DictReader(f))
    for r in rows:
        r["level"] = int(r["level"])
    return rows


def completed_cases(manifest: dict) -> list:
    return [(n, e) for n, e in sorted(manifest["cases"].items())
            if e.get("state") == "completed"]


def cmd_analyse(args) -> int:
    manifest = load_manifest()
    rows = []
    for name, entry in completed_cases(manifest):
        new, _ = analyse_case(CASES_DIR / name, entry["method"],
                              entry["estimator"], entry["level"],
                              entry["ranks"], entry.get("wall_y", 0.0),
                              entry.get("end_time", 0.0))
        for row in new:
            print(f"{name:18s} {row['estimator'] or '-':2s} "
                  f"{row.get('window', ''):4s} amplitude "
                  f"{100 * row['shear_amp_err']:+.2f} %, phase "
                  f"{row['shear_phase_err_deg']:+.2f} deg, velocity L2 "
                  f"(initial-value) {row['u_L2_all_ivp']:.3e}")
        rows += new
    if not rows:
        print("No completed cases: run ./Allverify first")
        return 1
    out = WORK_DIR / "results.csv"
    extra = ["window", "wall_y"] if any("window" in r for r in rows) else []
    write_results(out, rows, FIELDS + extra)
    print(f"Results: {out}")
    return 0


# --------------------------------------------------------------------- check
def check_analytic(p: dict) -> list:
    out = []
    k = 1.0 / p["delta"]
    worst = 0.0
    for t in (0.0, 0.3, 0.77):
        # rho nu du/ds at s = 0 by a central difference of the solution
        h = 1e-6
        d = (exact_u(p, h, t) - exact_u(p, -h, t)) / (2 * h)
        worst = max(worst, abs(p["rho"] * p["nu"] * d - exact_tau(p, t)))
    out.append(("PASS" if worst < CHECK_CLOSED_FORM * p["U0"] else "FAIL",
                "tau_w = rho nu du/ds at the wall"))
    ok = all(abs(slab_velocity(p, t) - p["U0"] * math.sin(p["omega"] * t))
             < 1e-12 for t in (0.0, 0.2, 0.9)) and \
        abs(exact_u(p, 0.0, 0.37) - slab_velocity(p, 0.37)) < 1e-12
    out.append(("PASS" if ok else "FAIL",
                "the slab velocity is U0 sin(omega t) and the solution "
                "equals it at the wall"))
    a, _ = exact_amp_phase(p)
    out.append(("PASS" if abs(a - max(abs(exact_tau(p, i / 4000.0))
                                      for i in range(4000))) < 1e-8 else "FAIL",
                f"exact shear amplitude {a:.6e} Pa"))
    return out


def check_tracked(rows: list, p: dict) -> list:
    """Consistency of the tracked results with the exact solution; nothing is
    compared with a measured value"""
    out = []
    aex, _ = exact_amp_phase(p)
    bad = [r for r in rows
           if abs(float(r["shear_amp_exact"]) / aex - 1.0) > 1e-5 or
           abs(float(r["shear_amp"]) / aex - 1.0 - float(r["shear_amp_err"]))
           > 1e-4]
    out.append(("PASS" if not bad else "FAIL",
                f"{len(rows)} tracked rows: the exact amplitude and the "
                "amplitude error agree with the closed form"))
    bad = [r for r in rows
           if int(r["cells"]) != mesh_shape(p, int(r["level"]))[0]
           * mesh_shape(p, int(r["level"]))[1]]
    out.append(("PASS" if not bad else "FAIL",
                "cell counts follow the mesh rule of the tutorial"))
    return out


def check_case(row: dict, tracked: dict) -> list:
    out = []
    name = f"{row['method']} {row['estimator']} L{row['level']}"
    d = abs(float(row["shear_amp_err"]) - float(tracked["shear_amp_err"]))
    out.append(("PASS" if d <= CHECK_SHEAR_AMP else "FAIL",
                f"{name}: amplitude error {float(row['shear_amp_err']):+.4f} "
                f"against tracked {float(tracked['shear_amp_err']):+.4f}"))
    d = abs(float(row["u_L2_all_ivp"]) - float(tracked["u_L2_all_ivp"]))
    out.append(("PASS" if d <= CHECK_VELOCITY else "FAIL",
                f"{name}: velocity L2 {float(row['u_L2_all_ivp']):.3e} "
                f"against tracked {float(tracked['u_L2_all_ivp']):.3e}"))
    return out


def read_series_csv(path: Path) -> dict:
    with path.open(newline="") as f:
        r = list(csv.DictReader(f))
    return {k: [float(x[k]) for x in r] for k in r[0]}


def check_series(path: Path, column: str, t: list, tau: list,
                 name: str) -> tuple:
    """A re-run shear history against a tracked history (time, column) at the
    tracked times"""
    if not path.exists():
        return ("SKIPPED", f"{name}: no tracked history {path.name}")
    ref = read_series_csv(path)
    if column not in ref:
        return ("SKIPPED", f"{name}: no column {column} in {path.name}")
    at = {round(x, 9): y for x, y in zip(t, tau)}
    miss = [x for x in ref["time"] if round(x, 9) not in at]
    if miss:
        return ("FAIL", f"{name}: {len(miss)} tracked times are missing")
    worst = max(abs(at[round(x, 9)] - y)
                for x, y in zip(ref["time"], ref[column]))
    return ("PASS" if worst <= CHECK_SERIES else "FAIL",
            f"{name}: shear history against {path.name}, largest difference "
            f"{worst:.2e} Pa")


def check_long_run(rows: list, name: str) -> list:
    """The eighth period against the fourth of one end-time control: the
    error against the initial-value solution is the same"""
    by = {r["window"]: r for r in rows if "window" in r}
    wins = sorted(by, key=lambda w: float(w.split("-")[0]))
    if len(wins) < 2:
        return []
    out = []
    first, last = by[wins[0]], by[wins[-1]]
    d = abs(float(last["shear_amp_err_ivp"])
            - float(first["shear_amp_err_ivp"]))
    out.append(("PASS" if d <= CHECK_LONG_RUN else "FAIL",
                f"{name}: amplitude error against the initial-value solution "
                f"{float(first['shear_amp_err_ivp']):+.5f} (period {wins[0]})"
                f" and {float(last['shear_amp_err_ivp']):+.5f} (period "
                f"{wins[-1]})"))
    return out


def check_variant(name: str, entry: dict, rows: list, series: dict,
                  index: dict) -> list:
    """Checks of a completed case: the tutorial and time-step cases against
    the tracked rows of the same key; the wall control against the tracked
    shear histories; an end-time control against its fourth period (tracked
    row) and against itself in the last period"""
    out = []
    wall, end = entry.get("wall_y", 0.0), entry.get("end_time", 0.0)
    if wall:
        for lab, (t, tau) in series.items():
            if lab not in ("ME", "FO"):
                continue
            path = CONTROLS_DIR / f"shear-yw{round(wall * 1000):04d}-" \
                f"L{entry['level']}.csv"
            out.append(check_series(path, f"cutLink-{lab}", t, tau,
                                    f"{name} {lab}"))
        return out
    for row in rows:
        if "window" in row and row["window"] != window_label(PERIOD):
            continue
        ref = index.get(row_key(row))
        if ref is None:
            out.append(("SKIPPED", f"{name}: no tracked row"))
        else:
            out += check_case(row, ref)
    if end:
        out += check_long_run(rows, name)
    return out


def cmd_check(args) -> int:
    p = parameters()
    results = check_analytic(p)
    tracked = read_results()
    if tracked:
        results += check_tracked(tracked, p)
    else:
        results.append(("SKIPPED", "no tracked results"))
    manifest = load_manifest()
    for name, entry in sorted(manifest["cases"].items()):
        if entry.get("state") == "failed":
            results.append(("FAIL", f"{name}: the run failed"))
    done = completed_cases(manifest)
    if not done:
        results.append(("SKIPPED", "no completed case in work/ "
                        "(run ./Allverify)"))
    index = {row_key(r): r for r in tracked + read_results(CONTROLS_CSV)}
    for name, entry in done:
        rows, series = analyse_case(
            CASES_DIR / name, entry["method"], entry["estimator"],
            entry["level"], entry["ranks"], entry.get("wall_y", 0.0),
            entry.get("end_time", 0.0))
        results += check_variant(name, entry, rows, series, index)
    for state, message in results:
        print(f"{state}: {message}")
    n = {s: sum(1 for r in results if r[0] == s)
         for s in ("PASS", "FAIL", "SKIPPED")}
    print(f"{n['PASS']} PASS, {n['FAIL']} FAIL, {n['SKIPPED']} SKIPPED "
          "(SKIPPED is not PASS)")
    return 1 if n["FAIL"] else 0


# -------------------------------------------------------------------- tables
def method_label(r: dict) -> str:
    if r["method"] == "cutLink":
        return {"ME": "cutLink, momentum exchange",
                "ST": "cutLink, surface traction",
                "FO": "cutLink, forcing"}[r["estimator"]]
    return {"ghostCell": "ghostCell (exploratory)"}.get(r["method"],
                                                         r["method"])


def master_table(rows: list, target: str = "periodic") -> str:
    """The tutorial configuration (time step 0.001 s) only. target is the
    exact solution of the errors: 'periodic' (the official periodic solution
    of the Stokes second problem) or 'ivp' (the exact solution of the
    problem that is simulated: the fluid at rest at t = 0, with the lid)"""
    sfx = "" if target == "periodic" else "_ivp"
    rows = [r for r in rows if abs(row_dt(r) - TUTORIAL_DT) < 1e-12]
    head = ("| Method | Level | Cells | Exact shear amplitude (Pa) | "
            "Computed amplitude (Pa) | Amplitude error | Phase error (deg) | "
            "Shear RMS error | u L1 | u L2 | u Linf | Noise | Ranks | "
            "Rank-h |\n| --- | --- | --- | --- | --- | --- | --- | --- | "
            "--- | --- | --- | --- | --- | --- |\n")
    lines = []
    for r in sorted(rows, key=sort_key):
        rh = f"{float(r['rank_h']):.4f}" if r.get("rank_h") else "n/a"
        ref_amp = float(r["shear_amp_exact"])
        shown = float(r["shear_amp"])
        lines.append(
            f"| {method_label(r)} | {r['level']} | {r['cells']} | "
            f"{ref_amp:.4e} | {shown:.4e} | "
            f"{100 * float(r['shear_amp_err' + sfx]):+.2f} % | "
            f"{float(r['shear_phase_err_deg' + sfx]):+.2f} | "
            f"{100 * float(r['shear_rms_err' + sfx]):.2f} % | "
            f"{float(r['u_L1_all' + sfx]):.2e} | "
            f"{float(r['u_L2_all' + sfx]):.2e} | "
            f"{float(r['u_Linf_all' + sfx]):.2e} | "
            f"{float(r['shear_noise']):.1e} | {r['ranks']} | {rh} |")
    return head + "\n".join(lines) + "\n"


def time_step_table(rows: list) -> str:
    """Shear amplitude error and velocity L2 norm (both against the
    initial-value solution) at every time step, for the Levels and methods
    of the time-step control, with the tutorial step 0.001 s in its column"""
    steps = sorted({row_dt(r) for r in rows}, reverse=True)
    levels = sorted({int(r["level"]) for r in rows if row_dt(r) != TUTORIAL_DT})
    head = ("| Method | Level | " + " | ".join(f"dt = {fmt(d)} s" for d in steps)
            + " |\n| --- | --- | " + " | ".join("---" for _ in steps) + " |\n")
    lines = []
    for method, est in (DT_VARIANTS[:1] + (("cutLink", "ST"),)
                        + DT_VARIANTS[1:]):
        for level in levels:
            cells = []
            for d in steps:
                hit = [r for r in rows if (r["method"], r["estimator"]) ==
                       (method, est) and int(r["level"]) == level
                       and abs(row_dt(r) - d) < 1e-12]
                cells.append(
                    f"{100 * float(hit[0]['shear_amp_err_ivp']):+.2f} % / "
                    f"{float(hit[0]['u_L2_all_ivp']):.2e}" if hit else "n/a")
            lines.append(f"| {method_label({'method': method, 'estimator': est})}"
                         f" | {level} | " + " | ".join(cells) + " |")
    return head + "\n".join(lines) + "\n"


def aligned_table(rows: list) -> str:
    """Shear amplitude, phase and RMS error against the initial-value
    solution of penalty, incremental and ghostCell, as the solver writes the
    inertia and with the inertia averaged over two steps (diagnostic), at
    the tutorial step and at 0.00025 s"""
    head = ("| Method | Level | dt (s) | Amplitude error | Phase error (deg) | "
            "RMS error | Amplitude error, aligned | Phase error, aligned | "
            "RMS error, aligned |\n| --- | --- | --- | --- | --- | --- | --- "
            "| --- | --- |\n")
    lines = []
    for method in ("penalty", "incremental", "ghostCell"):
        for level in (2, 3, 4):
            for dt in (TUTORIAL_DT, 0.00025):
                hit = [r for r in rows if r["method"] == method
                       and int(r["level"]) == level
                       and abs(row_dt(r) - dt) < 1e-12 and r.get("shear_amp_err_al")]
                if not hit:
                    continue
                r = hit[0]
                lines.append(
                    f"| {method_label(r)} | {level} | {fmt(dt)} | "
                    f"{100 * float(r['shear_amp_err_ivp']):+.2f} % | "
                    f"{float(r['shear_phase_err_deg_ivp']):+.2f} | "
                    f"{100 * float(r['shear_rms_err_ivp']):.2f} % | "
                    f"{100 * float(r['shear_amp_err_al']):+.2f} % | "
                    f"{float(r['shear_phase_err_deg_al']):+.2f} | "
                    f"{100 * float(r['shear_rms_err_al']):.2f} % |")
    return head + "\n".join(lines) + "\n"


def cmd_tables(args) -> int:
    rows = read_results()
    if not rows:
        fail("no tracked results")
    print(master_table(rows, "periodic"))
    print(master_table(rows, "ivp"))
    print(time_step_table(rows))
    print(aligned_table(rows))
    return 0


# ---------------------------------------------------------------------- unit
def cmd_unit(args) -> int:
    import unittest
    here = Path(__file__).resolve().parent
    sys.path.insert(0, str(here))
    suite = unittest.defaultTestLoader.discover(
        str(here), pattern="test_oscillating_wall_stokes_layer_verification.py")
    result = unittest.TextTestRunner(verbosity=1, warnings="error").run(suite)
    return 0 if result.wasSuccessful() else 1


COMMANDS = ("run", "generate", "dry-run", "status", "analyse", "check",
            "tables", "unit")


def normalise_argv(argv: list) -> list:
    """No subcommand means 'run': ./Allverify --runs L1 L2 is
    ./Allverify run --runs L1 L2 (and ./Allverify alone is the smoke run)"""
    argv = list(argv)
    if not argv or argv[0] not in COMMANDS + ("-h", "--help"):
        argv = ["run"] + argv
    return argv


def main(argv=None) -> int:
    argv = normalise_argv(sys.argv[1:] if argv is None else argv)
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="command")
    for name in ("run", "generate", "dry-run"):
        s = sub.add_parser(name)
        s.add_argument("--runs", nargs="+", default=list(DEFAULT_RUNS),
                       help="named runs: " + ", ".join(RUNS))
        s.add_argument("--methods", nargs="+", default=[])
        s.add_argument("--levels", nargs="+", type=int, default=[])
        s.add_argument("--ranks", type=int, default=0,
                       help="MPI ranks (default: by level)")
        s.add_argument("--end-time", type=float, default=0.0,
                       help="endTime (s) of every selected case; the cases "
                       "are tagged -T<time>")
        s.add_argument("--wall-y", type=float, default=0.0,
                       help="y (m) of the wall: the top vertices of "
                       "slab.stl are moved there; the cases are tagged "
                       "-yw<y>")
        s.add_argument("--force", action="store_true",
                       help="regenerate existing cases")
        s.add_argument("--work-dir", default="",
                       help="cases and manifest of these runs (default "
                       "verification/work)")
    for name in ("status", "analyse", "check"):
        s = sub.add_parser(name)
        s.add_argument("--work-dir", default="",
                       help="work directory of the runs")
    for name in ("tables", "unit"):
        sub.add_parser(name)
    args = parser.parse_args(argv)
    work = getattr(args, "work_dir", "")
    if work:
        set_work_dir(Path(work).resolve())
    handlers = {"run": cmd_run, "generate": cmd_generate,
                "dry-run": cmd_dry_run, "status": cmd_status,
                "analyse": cmd_analyse, "check": cmd_check,
                "tables": cmd_tables, "unit": cmd_unit}
    return handlers[args.command](args)


def set_work_dir(path: Path) -> None:
    global WORK_DIR, CASES_DIR, MANIFEST
    WORK_DIR, CASES_DIR, MANIFEST = path, path / "cases", path / "manifest.json"


if __name__ == "__main__":
    sys.exit(main())
