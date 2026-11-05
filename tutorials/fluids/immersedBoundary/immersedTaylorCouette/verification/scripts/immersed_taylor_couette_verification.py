#!/usr/bin/env python3
"""Verification of the immersedTaylorCouette tutorial: the azimuthal velocity
and the torque on a rotating immersed cylinder against the exact solution of
the circular Couette flow, with the cutLink, penalty, legacy incremental and
(exploratory) ghostCell methods of immersedBoundaryForce.

The exact solution needs no reference solver. With R1 and R2 the radii of the
rotor and of the stator, omega the angular velocity of the rotor and L_z the
thickness of the mesh,

    u_theta(r) = A r + B / r,   A = -omega R1^2 / (R2^2 - R1^2),
                                B =  omega R1^2 R2^2 / (R2^2 - R1^2)
    tau_r,theta = mu r d(u_theta/r)/dr = -2 mu B / r^2
    T_rotor = 2 pi R1 L_z R1 tau(R1) = -4 pi mu B L_z

where T_rotor is the torque that the fluid exerts on the rotor about its
axis: it is negative because the fluid resists the rotation. The parameters
are read from the tutorial files (radii from the surfaces, omega and the
centre of rotation from constant/fvOptions, the viscosity from
constant/transportProperties, L_z from system/blockMeshDict); nothing is
copied from the README.

Four things are reported, and kept apart:

  A. the rotor torque against the exact torque, as a relative difference
     (T - T_exact)/T_exact, so that a positive value is a torque that is too
     large in magnitude. It is the mean over the window 8 <= t <= 10 s,
     with the variance in the same window. The flow is steady when the
     relative change of the torque over the window and the relative L2
     change of the velocity between t = 5 and 10 s are both below 1e-3
     (predefined); a case that is not steady is reported with its mean and
     variance and is never described as converged;
  B. the velocity at t = 10 s in the fluid cells (those with a cell centre
     in R1 < r < R2) against the exact velocity at the cell centres, as
     volume-weighted L1, L2 and Linf norms of the azimuthal difference and of
     the radial velocity (the exact value is zero), relative to omega R1.
     'far' are the cells at least two cell widths from both walls, and
     'iface' the others; both are reported next to the whole set;
  C. an independent check of the torque: the torque transmitted across
     circles inside the gap, from the viscous stress of the velocity field
     alone (it does not use the force output of the solver). In a steady flow
     it equals the torque on the rotor;
  D. the observed cost, with the conditions of the run.

Each case is a copy of the tutorial under verification/work/cases/ with only
the method settings, the mesh level and the decomposition changed. The
tutorial itself is never modified. Everything but the OpenFOAM commands is
plain Python 3 (standard library only). A manifest, work/manifest.json,
records every case with the SHA-256 of its inputs.

Usage:
    ./Allverify [run] [--runs ...]   generate and run locally (does not
                                     analyse or check: follow with the
                                     'analyse' and 'check' commands); 'run'
                                     may be omitted: ./Allverify --runs L1 L2
    ./Allverify unit                 exact tests that need no OpenFOAM
    ./Allverify dry-run --runs ...   list the cases and ranks
    ./Allverify generate --runs ...  generate the cases only
    (run, generate and dry-run take --methods, --levels, --ranks, --force
    and --end-time T, which writes endTime T into controlDict and tags the
    case name with -tT; the named runs 'dt-extended' and 'steady20' are the
    control runs of reference/results/controls/, see controls.md)
    ./Allverify status               state of the cases in the manifest
    ./Allverify analyse              analyse the completed cases
    ./Allverify check                PASS / FAIL / SKIPPED
    ./Allverify tables               the tables of results.md from the data

The default run, 'smoke', is the cutLink case on Level 1 (1024 cells, ten
seconds of flow): about 20 s on one core.
"""

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

LEVELS = (1, 2, 3, 4)
MESH_CELLS_ACROSS = {1: 32, 2: 64, 3: 128, 4: 256}
SIDE = 0.32                       # side of the square mesh
# (method, estimator) of the variants; the estimator selects the torque
# column of the cutLink force file
VARIANTS = (("cutLink", "ME"), ("cutLink", "ST"), ("cutLink", "FO"),
            ("penalty", ""), ("incremental", ""), ("ghostCell", ""))
ESTIMATOR_NAMES = {"ME": "momentumExchange", "ST": "surfaceTraction",
                   "FO": "forcing"}
N_OUTER = {"cutLink": 1, "penalty": 1, "incremental": 5, "ghostCell": 1}

# Settings written explicitly into every case: those of the stationary and
# oscillating cylinder verifications; ghostCell is documented for static
# bodies and is exploratory here (the rotor does not change its position)
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
WINDOW = (8.0, 10.0)              # steady-state window, all methods
STEADY_DRIFT = 1e-3               # relative change of T over the window
STEADY_DU = 1e-3                  # relative L2 change of u_theta, 5 -> 10 s
FIELD_FRACTIONS = (0.3, 0.4, 0.5, 0.6, 0.7)   # circles of the field torque

# Predefined check tolerances, independent of the measured values
CHECK_TORQUE = 0.5e-2             # |difference of the torque error| (smoke)
CHECK_VELOCITY = 1.0e-3           # |difference of the velocity L2| (smoke)
CHECK_EXACT_TORQUE = 1e-9         # relative, the closed form against the
                                  # direct integration of the viscous stress
CHECK_FIELD_TORQUE = 0.01         # |difference| of the torque from the field
                                  # and from the force output (Level >= 2)

def load_observed_rank_hours() -> dict:
    """Observed allocated rank-hours of the tracked runs (serial except Level
    4), {(method, estimator, level, dt_cap): rank_h}, from the rank_h column of
    reference/results/cases.csv and of controls/cases-timestep.csv"""
    out = {}
    for path in (CASES_CSV, RESULTS_DIR / "controls" / "cases-timestep.csv"):
        for r in read_results(path):
            if r.get("rank_h") not in ("", None):
                out[row_key(r)] = float(r["rank_h"])
    return out

# dt_cap = 0 is the tutorial: the time step is limited by the Courant
# number 0.5 and, because the cells halve with each level, is itself halved
# with each level. A positive dt_cap adds maxDeltaT = dt_cap to controlDict,
# the control that separates the time step from the mesh
DT_VARIANTS = (("cutLink", "ME"), ("cutLink", "ST"), ("penalty", ""),
               ("incremental", ""), ("ghostCell", ""))
DT_CONTROLS = [(lv, cap) for lv in (1, 2) for cap in (0.01, 0.005, 0.0025)] \
    + [(3, 0.01), (3, 0.005)]
# The control runs of reference/results/controls/cases-timestep.csv that are
# not in DT_CONTROLS: Level 4 (4 ranks) at the caps 0.01 and 0.005 (four
# variants, the momentum exchange among them) and 0.0025 (penalty and
# incremental), Level 3 at the cap 0.0025 (penalty, incremental, ghostCell) and
# the penalty method on Level 3 at 0.00125. The keys are exactly those of the
# CSV (a test compares them)
DT_EXTENDED = [
    ("cutLink", "ME", 4, 0.01), ("penalty", "", 4, 0.01),
    ("incremental", "", 4, 0.01), ("ghostCell", "", 4, 0.01),
    ("cutLink", "ME", 4, 0.005), ("penalty", "", 4, 0.005),
    ("incremental", "", 4, 0.005), ("ghostCell", "", 4, 0.005),
    ("penalty", "", 3, 0.0025), ("penalty", "", 4, 0.0025),
    ("incremental", "", 3, 0.0025), ("incremental", "", 4, 0.0025),
    ("ghostCell", "", 3, 0.0025), ("penalty", "", 3, 0.00125)]
# The 20 s runs of reference/results/controls/steady-extended.csv: (method,
# estimator, level, dt_cap, end_time); dt_cap 0 is the tutorial step
STEADY_END = 20.0
STEADY20 = [
    ("ghostCell", "", 1, 0.01, STEADY_END), ("ghostCell", "", 1, 0.005, STEADY_END),
    ("incremental", "", 1, 0.0, STEADY_END),
    ("incremental", "", 1, 0.01, STEADY_END),
    ("incremental", "", 1, 0.005, STEADY_END),
    ("incremental", "", 1, 0.0025, STEADY_END),
    ("incremental", "", 2, 0.0025, STEADY_END),
    ("penalty", "", 1, 0.0, STEADY_END), ("penalty", "", 1, 0.01, STEADY_END),
    ("penalty", "", 1, 0.005, STEADY_END), ("penalty", "", 1, 0.0025, STEADY_END),
    ("penalty", "", 2, 0.0, STEADY_END), ("penalty", "", 2, 0.01, STEADY_END),
    ("penalty", "", 2, 0.005, STEADY_END), ("penalty", "", 2, 0.0025, STEADY_END)]
RUNS = {
    "smoke": ("cutLink (momentum exchange), Level 1",
              [("cutLink", "ME", 1, 0.0)]),
    "L1": ("all variants, Level 1", [(m, e, 1, 0.0) for m, e in VARIANTS]),
    "L2": ("all variants, Level 2", [(m, e, 2, 0.0) for m, e in VARIANTS]),
    "L3": ("all variants, Level 3", [(m, e, 3, 0.0) for m, e in VARIANTS]),
    "matrix": ("all variants, Levels 1-3 (the tracked matrix)",
               [(m, e, lv, 0.0) for lv in (1, 2, 3) for m, e in VARIANTS]),
    "L4": ("the five variants of the tracked set, Level 4 (4 ranks)",
           [(m, e, 4, 0.0) for m, e in VARIANTS if e != "FO"]),
    "dt": ("time-step control: maxDeltaT 0.01, 0.005, 0.0025 s on Levels "
           "1-2 and 0.01, 0.005 s on Level 3, five variants",
           [(m, e, lv, cap) for lv, cap in DT_CONTROLS
            for m, e in DT_VARIANTS]),
    "dt-extended": ("time-step control, Levels 3-4 and the 0.00125 s cap: "
                    "the 14 rows of controls/cases-timestep.csv (Level 4 on "
                    "4 ranks; about 30 rank-h in all, the momentum exchange at Level 4 alone 8.5)", list(DT_EXTENDED)),
    "st-audit": ("cutLink surface traction, Levels 1-4 (Level 4 on 4 ranks): "
                 "the runs of the surface-traction audit, "
                 "scripts/taylor_couette_st_audit.py",
                 [("cutLink", "ST", lv, 0.0) for lv in LEVELS]),
    "steady20": ("the 15 runs to t = 20 s of controls/steady-extended.csv "
                 "(about 2 serial hours)", list(STEADY20)),
}
DEFAULT_RUNS = ("smoke",)
LEVEL_RANKS = {1: 1, 2: 1, 3: 1, 4: 4}

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


def stl_radii(path: Path) -> tuple:
    """Smallest and largest distance of the vertices of an ascii STL from
    the z axis, and the z range. The vertices are single-precision numbers
    (0.05 is stored as 0.0500000513), so the radii are rounded to 1e-6 m;
    the polygon of the surface is inscribed in the circle (256 sides: the
    mean radius is lower by 7.5e-5 relative)"""
    pts = [tuple(float(x) for x in m)
           for m in re.findall(r"vertex\s+(\S+)\s+(\S+)\s+(\S+)",
                               path.read_text())]
    if not pts:
        fail(f"no vertices in {path}")
    r = [math.hypot(x, y) for x, y, _ in pts]
    z = [p[2] for p in pts]
    return round(min(r), 6), round(max(r), 6), min(z), max(z)


def parameters(tutorial: Path = TUTORIAL_DIR) -> dict:
    """The physical parameters of the tutorial, from its files: the rotor
    and stator radii from their surfaces (the rotor is the cylinder and the
    stator the annulus), omega and the centre and axis of the rotation from
    constant/fvOptions, rho and nu from constant/transportProperties and the
    thickness of the mesh from system/blockMeshDict"""
    fv = (tutorial / "constant/fvOptions").read_text()
    tp = (tutorial / "constant/transportProperties").read_text()
    bm = (tutorial / "system/blockMeshDict").read_text()
    r1_min, r1_max, _, _ = stl_radii(tutorial / "constant/triSurface/rotor.stl")
    r2_min, r2_max, _, _ = stl_radii(tutorial / "constant/triSurface/stator.stl")
    if r1_min > 1e-9:
        fail("the rotor is not a solid cylinder")
    motion = re.search(r"motion\s*\{(.*?)\}", fv, re.DOTALL).group(1)
    origin = tuple(float(x) for x in re.search(
        r"origin\s+\(([^)]*)\)", motion).group(1).split())
    axis = tuple(float(x) for x in re.search(
        r"axis\s+\(([^)]*)\)", motion).group(1).split())
    cofrs = {tuple(float(x) for x in m.split())
             for m in re.findall(r"CofR\s+\(([^)]*)\)", fv)}
    if len(cofrs) != 1:
        fail("the bodies of the tutorial have different centres of rotation")
    cofr = cofrs.pop()
    vertices = re.search(r"vertices\s*\((.*?)\n\);", bm, re.DOTALL).group(1)
    pts = [tuple(float(x) for x in v.split())
           for v in re.findall(r"\(([-0-9. eE]+)\)", vertices)]
    zs = sorted({q[2] for q in pts})
    xs, ys = sorted({q[0] for q in pts}), sorted({q[1] for q in pts})
    cells = re.search(
        r'^\s*n\s+#eval\s+"round\(\s*([0-9.]+)\s*\*\s*pow\(\s*([0-9.]+)'
        r'\s*,\s*\$MESH_LEVEL\s*-\s*([0-9.]+)\s*\)\s*\)"', bm, re.MULTILINE)
    blocks = re.search(r"hex\s*\([^)]*\)\s*\(\s*\$n\s+\$n\s+1\s*\)", bm)
    return {"R1": r1_max, "R2": r2_min, "omega": dict_scalar(motion, "omega"),
            "origin": origin, "axis": axis, "CofR": cofr,
            "nu": dict_scalar(tp, "nu"), "rho": dict_scalar(tp, "rho"),
            "Lz": zs[-1] - zs[0], "side": SIDE,
            "mesh_x": (xs[0], xs[-1]), "mesh_y": (ys[0], ys[-1]),
            "cells_formula": tuple(float(g) for g in cells.groups())
            if cells else None,
            "square_block": blocks is not None}


# ------------------------------------------------------------ exact solution
def exact_AB(p: dict) -> tuple:
    d = p["R2"] ** 2 - p["R1"] ** 2
    return (-p["omega"] * p["R1"] ** 2 / d,
            p["omega"] * p["R1"] ** 2 * p["R2"] ** 2 / d)


def exact_utheta(p: dict, r: float) -> float:
    a, b = exact_AB(p)
    return a * r + b / r


def exact_shear(p: dict, r: float) -> float:
    """tau_r,theta = mu r d(u_theta/r)/dr = -2 mu B / r^2"""
    return -2.0 * p["rho"] * p["nu"] * exact_AB(p)[1] / (r * r)


def exact_rotor_torque(p: dict) -> float:
    """Torque of the fluid on the rotor about the axis, -4 pi mu B L_z"""
    return -4.0 * math.pi * p["rho"] * p["nu"] * exact_AB(p)[1] * p["Lz"]


def integrated_rotor_torque(p: dict) -> float:
    """The same torque from the stress on the surface of the rotor,
    tau(R1) * (2 pi R1 L_z) * R1: tests the closed form (sign included)"""
    r = p["R1"]
    return exact_shear(p, r) * 2.0 * math.pi * r * p["Lz"] * r


def exact_stator_torque(p: dict) -> float:
    """Steady angular momentum balance of the fluid: the torques on the two
    bodies are equal and opposite"""
    return -exact_rotor_torque(p)


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


def interpolate(t: list, y: list, x: float) -> float:
    if x <= t[0]:
        return y[0]
    if x >= t[-1]:
        return y[-1]
    lo, hi = 0, len(t) - 1
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if t[mid] <= x:
            lo = mid
        else:
            hi = mid
    return y[lo] + (y[hi] - y[lo]) * (x - t[lo]) / (t[hi] - t[lo])


def window_series(t: list, y: list, window=WINDOW) -> tuple:
    """The samples in the window, with the ends interpolated"""
    a, b = window
    tt = [a] + [x for x in t if a < x < b] + [b]
    return tt, [interpolate(t, y, x) for x in tt]


def trapezoid(t: list, f: list) -> float:
    return sum(0.5 * (f[i] + f[i + 1]) * (t[i + 1] - t[i])
               for i in range(len(t) - 1))


def window_stats(t: list, y: list, window=WINDOW) -> tuple:
    """Time-weighted mean and variance in the window (adaptive time step)"""
    tt, yy = window_series(t, y, window)
    span = window[1] - window[0]
    mean = trapezoid(tt, yy) / span
    var = trapezoid(tt, [(v - mean) ** 2 for v in yy]) / span
    return mean, var


# ------------------------------------------------------------ field analysis
def cell_centres(level: int) -> tuple:
    """Centres of the uniform mesh (x runs fastest, as in blockMesh)"""
    n = MESH_CELLS_ACROSS[level]
    h = SIDE / n
    xs = [-SIDE / 2 + (i + 0.5) * h for i in range(n)]
    return n, h, [(x, y) for y in xs for x in xs]


def velocity_errors(p: dict, level: int, u10: list, u5: list) -> dict:
    n, h, centres = cell_centres(level)
    scale = p["omega"] * p["R1"]
    vol = h * h * p["Lz"]
    sets = {"all": [], "far": [], "iface": []}
    ur_all, eth = [], []
    num5 = den5 = 0.0
    for (x, y), u, w in zip(centres, u10, u5):
        r = math.hypot(x, y)
        if not p["R1"] < r < p["R2"]:
            continue
        c, s = x / r, y / r
        ut = -s * u[0] + c * u[1]
        ur = c * u[0] + s * u[1]
        ut5 = -s * w[0] + c * w[1]
        e = ut - exact_utheta(p, r)
        far = (r > p["R1"] + 2 * h) and (r < p["R2"] - 2 * h)
        sets["all"].append((e, ur))
        sets["far" if far else "iface"].append((e, ur))
        num5 += (ut - ut5) ** 2
        den5 += ut ** 2

    def norms(items, comp):
        v = [abs(it[comp]) for it in items]
        # equal cell volumes: the volume-weighted norm is the plain mean
        return (sum(v) / len(v) / scale,
                math.sqrt(sum(x * x for x in v) / len(v)) / scale,
                max(v) / scale)
    out = {"du_5_10": math.sqrt(num5 / den5), "fluid_cells": len(sets["all"]),
           "cell_volume": vol}
    for name in ("all", "far", "iface"):
        for comp, label in ((0, "ut"), (1, "ur")):
            if name != "all" and label == "ur":
                continue
            l1, l2, li = norms(sets[name], comp)
            suffix = "" if name == "all" else f"_{name}"
            out[f"{label}_L1{suffix}"] = l1
            out[f"{label}_L2{suffix}"] = l2
            out[f"{label}_Linf{suffix}"] = li
    return out


def bilinear(level: int, comp: list, px: float, py: float) -> float:
    n = MESH_CELLS_ACROSS[level]
    h = SIDE / n
    fx = (px + SIDE / 2) / h - 0.5
    fy = (py + SIDE / 2) / h - 0.5
    i = min(max(int(math.floor(fx)), 0), n - 2)
    j = min(max(int(math.floor(fy)), 0), n - 2)
    a, b = fx - i, fy - j
    return ((1 - a) * (1 - b) * comp[j * n + i] + a * (1 - b) * comp[j * n + i + 1]
            + (1 - a) * b * comp[(j + 1) * n + i]
            + a * b * comp[(j + 1) * n + i + 1])


def field_torque(p: dict, level: int, u: list, rs: float, m: int = 720) -> float:
    """Torque transmitted across the circle r = rs from the viscous stress of
    the velocity field alone, rho nu (du_theta/dr - u_theta/r) 2 pi rs^2 L_z,
    with the derivative from a central difference over one cell width"""
    h = SIDE / MESH_CELLS_ACROSS[level]
    ux = [v[0] for v in u]
    uy = [v[1] for v in u]

    def ut(r, k):
        a = (k + 0.5) * 2.0 * math.pi / m
        px, py = r * math.cos(a), r * math.sin(a)
        return (-math.sin(a) * bilinear(level, ux, px, py)
                + math.cos(a) * bilinear(level, uy, px, py))
    d = sum(ut(rs + h, k) - ut(rs - h, k) for k in range(m)) / (2 * h * m)
    mean = sum(ut(rs, k) for k in range(m)) / m
    tau = p["rho"] * p["nu"] * (d - mean / rs)
    return tau * 2.0 * math.pi * rs * rs * p["Lz"]


def median(values: list) -> float:
    v = sorted(values)
    k = len(v) // 2
    return v[k] if len(v) % 2 else 0.5 * (v[k - 1] + v[k])


# ------------------------------------------------------------ case handling
def dt_tag(dt_cap: float) -> str:
    return f"-dt{fmt(dt_cap)}" if dt_cap else ""


def case_name(method: str, est: str, level: int, dt_cap: float = 0.0,
              end_time: float = 0.0) -> str:
    """-dt<cap> for a step cap and -t<end> for an end time other than the
    10 s of the tutorial"""
    return (f"{method}{'-' + est if method == 'cutLink' else ''}-L{level}"
            f"{dt_tag(dt_cap)}{f'-t{fmt(end_time)}' if end_time else ''}")


def load_manifest() -> dict:
    if MANIFEST.exists():
        return json.loads(MANIFEST.read_text())
    return {"cases": {}}


def save_manifest(manifest: dict) -> None:
    WORK_DIR.mkdir(parents=True, exist_ok=True)
    MANIFEST.write_text(json.dumps(manifest, indent=1, sort_keys=True) + "\n")


def plan(args) -> list:
    """(method, estimator, level, name, ranks, dt_cap, end_time) of the cases
    of the named runs; --end-time overrides the end time of every case"""
    out, seen = [], set()
    override = float(getattr(args, "end_time", 0.0) or 0.0)
    for name in args.runs:
        if name not in RUNS:
            fail(f"unknown run '{name}'; choose from {', '.join(RUNS)}")
        for method, est, level, cap, *rest in RUNS[name][1]:
            end = override or (rest[0] if rest else 0.0)
            if args.methods and method not in args.methods:
                continue
            if args.levels and level not in args.levels:
                continue
            key = (method, est, level, cap, end)
            if key not in seen:
                seen.add(key)
                ranks = args.ranks or LEVEL_RANKS[level]
                out.append((method, est, level,
                            case_name(method, est, level, cap, end), ranks,
                            cap, end))
    return out


def replace_entry(path: Path, key: str, value: str) -> None:
    text = path.read_text()
    pattern = re.compile(rf"^(\s*{re.escape(key)}\s+)[^;\n]*;", re.MULTILINE)
    new, count = pattern.subn(lambda m: f"{m.group(1)}{value};", text)
    if count != 1:
        fail(f"found {count} '{key}' entries in {path}")
    path.write_text(new)


def write_fv_options(path: Path, method: str, est: str) -> None:
    """The tutorial's fvOptions with the method settings of METHOD_SETTINGS;
    the motion blocks are never rewritten"""
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
    if body:
        extra = "".join(f"            {k.ljust(12)}{v};\n" for k, v in body)
        new, count = re.subn(r"(\n\s*surface\s+\"[a-z]+\.stl\";\n)",
                             lambda m: m.group(1) + extra, new)
        if count != 2:
            fail("expected the surface entries of the two bodies")
    path.write_text(new)


def write_decompose_par_dict(path: Path, ranks: int) -> None:
    """Slabs along y. The pressure reference point of the tutorial,
    (0 0.075 0.05), lies on the plane x = 0, which a decomposition along x
    cuts; setRefCell then finds it in two subdomains and stops"""
    text = FOAM_HEADER.format(location="system", object="decomposeParDict")
    text += (f"numberOfSubdomains  {ranks};\n\nmethod          simple;\n\n"
             f"simpleCoeffs\n{{\n    n           (1 {ranks} 1);\n}}\n\n"
             "// **************************************************************"
             "*********** //\n")
    path.write_text(text)


RUN_SCRIPT = """\
#!/bin/bash
# Generated by verification/scripts/immersed_taylor_couette_verification.py
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
                  dt_cap: float = 0.0, end_time: float = 0.0) -> Path:
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
    write_decompose_par_dict(case / "system/decomposeParDict", ranks)
    if dt_cap:
        control = case / "system/controlDict"
        text = control.read_text()
        new, count = re.subn(r"^(maxCo\s+[^;]+;)",
                             lambda m: f"{m.group(1)}\nmaxDeltaT      "
                             f"{fmt(dt_cap)};", text, flags=re.MULTILINE)
        if count != 1:
            fail("maxCo entry of the tutorial controlDict not found")
        control.write_text(new)
    if end_time:
        replace_entry(case / "system/controlDict", "endTime", fmt(end_time))
    script = case / "Allrun.verify"
    script.write_text(RUN_SCRIPT.format(
        level=level, ranks=ranks,
        body=PARALLEL_BODY if ranks > 1 else SERIAL_BODY))
    script.chmod(0o755)
    return case


def cmd_generate(args) -> int:
    manifest = load_manifest()
    for method, est, level, name, ranks, cap, end in plan(args):
        case = generate_case(method, est, level, name, ranks, args.force,
                             dt_cap=cap, end_time=end)
        manifest["cases"][name] = {
            "method": method, "estimator": est, "level": level,
            "dt_cap": cap, "end_time": end, "ranks": ranks,
            "inputs_sha256": input_hashes(case),
            "generated": datetime.now().isoformat(timespec="seconds")}
        print(f"{name}: generated in {case}")
    save_manifest(manifest)
    return 0


def cmd_dry_run(args) -> int:
    p = parameters()
    print(f"Tutorial: {TUTORIAL_DIR.name}; R1 = {p['R1']}, R2 = {p['R2']}, "
          f"omega = {p['omega']}, nu = {p['nu']}, rho = {p['rho']}, "
          f"L_z = {p['Lz']}")
    print(f"Exact rotor torque: {exact_rotor_torque(p):.6e} N m")
    total = 0.0
    observed = load_observed_rank_hours()
    for method, est, level, name, ranks, cap, end in plan(args):
        cells = MESH_CELLS_ACROSS[level] ** 2
        # the observed cost is that of the 10 s tutorial run of the same key;
        # a 20 s run costs about twice as much and has no tracked cost
        obs = None if end else observed.get((method, est, level, cap))
        total += obs or 0.0
        print(f"  {name:20s} {cells:6d} cells, {ranks} rank(s), observed "
              f"{'%.3f' % obs if obs else 'n/a'} allocated rank-h")
    print(f"Observed total of the listed cases: {total:.2f} rank-h "
          "(single runs on a shared node, not a performance model)")
    return 0


def is_openfoam_com(environ=None) -> bool:
    """OpenFOAM.com (ESI): WM_PROJECT is OpenFOAM and WM_PROJECT_VERSION
    starts with 'v' (v2412); OpenFOAM.org has a number or 'dev' (the rule of
    solids4FoamScripts.sh), foam-extend has WM_PROJECT foam"""
    environ = os.environ if environ is None else environ
    return environ.get("WM_PROJECT") == "OpenFOAM" and \
        environ.get("WM_PROJECT_VERSION", "").startswith("v")


def check_environment(ranks: int) -> str:
    """Empty if the cases can be run, else the reason they cannot"""
    if "WM_PROJECT_DIR" not in os.environ:
        return "no OpenFOAM environment (WM_PROJECT_DIR is not set)"
    if not is_openfoam_com():
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
        if entry.get("state") == "completed" and (CASES_DIR / item[3]).exists() \
                and not args.force:
            print(f"{item[3]}: completed, kept (use --force to rerun)")
        else:
            todo.append(item)
    for method, est, level, name, ranks, cap, end in todo:
        case = generate_case(method, est, level, name, ranks, True, dt_cap=cap,
                             end_time=end)
        manifest["cases"][name] = {
            "method": method, "estimator": est, "level": level,
            "dt_cap": cap, "end_time": end, "ranks": ranks,
            "inputs_sha256": input_hashes(case),
            "generated": datetime.now().isoformat(timespec="seconds")}
        print(f"{name}: generated in {case}")
    save_manifest(manifest)
    status = 0
    for method, est, level, name, ranks, cap, end in todo:
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


def analyse_case(case: Path, method: str, est: str, level: int,
                 ranks: int = 1, dt_cap: float = 0.0) -> tuple:
    """The result row of a completed case, and its torque series"""
    p = parameters()
    rotor = read_force_file(case / "postProcessing/immersedBoundary/0/rotor.dat")
    stator = read_force_file(case / "postProcessing/immersedBoundary/0/stator.dat")
    t = rotor["t"]
    tz = [v[2] for v in rotor["torque"]]
    exact = exact_rotor_torque(p)
    mean, var = window_stats(t, tz)
    a, b = WINDOW
    drift = abs(interpolate(t, tz, b) - interpolate(t, tz, a)) / abs(mean)
    stator_mean, _ = window_stats(stator["t"], [v[2] for v in stator["torque"]])
    u10 = read_vector_field(case / "10/U")
    u5 = read_vector_field(case / "5/U")
    err = velocity_errors(p, level, u10, u5)
    fields = [field_torque(p, level, u10,
                           p["R1"] + (p["R2"] - p["R1"]) * f)
              for f in FIELD_FRACTIONS]
    row = {
        "method": method, "estimator": est, "level": level,
        "dt_cap": fmt(dt_cap) if dt_cap else "",
        "cells": MESH_CELLS_ACROSS[level] ** 2, "ranks": ranks,
        "steps": len(t), "T_mean": mean, "T_exact": exact,
        "T_relerr": (mean - exact) / exact, "T_std_rel": math.sqrt(var) / abs(exact),
        "T_drift_rel": drift, "dU_5_10": err["du_5_10"],
        "steady": int(drift < STEADY_DRIFT and err["du_5_10"] < STEADY_DU),
        "T_field_rel": median(fields) / exact - 1.0,
        "T_stator_balance": (mean + stator_mean) / abs(exact),
        "wall_s": wall_seconds(case)}
    for k in ("ut_L1", "ut_L2", "ut_Linf", "ur_L1", "ur_L2", "ur_Linf",
              "ut_L2_far", "ut_Linf_far", "ut_L2_iface", "ut_Linf_iface"):
        row[k] = err[k]
    return row, (t, tz)


def row_key(row: dict) -> tuple:
    cap = row.get("dt_cap", "")
    return (row["method"], row["estimator"], int(row["level"]),
            float(cap) if cap not in ("", None) else 0.0)


FIELDS = ["method", "estimator", "level", "dt_cap", "cells", "ranks", "steps",
          "T_mean", "T_exact", "T_relerr", "T_std_rel", "T_drift_rel",
          "dU_5_10", "steady", "T_field_rel", "T_stator_balance",
          "ut_L1", "ut_L2", "ut_Linf", "ur_L1", "ur_L2", "ur_Linf",
          "ut_L2_far", "ut_Linf_far", "ut_L2_iface", "ut_Linf_iface",
          "wall_s", "rank_h"]


def sort_key(r: dict) -> tuple:
    cap = row_key(r)[3]
    return (cap > 0, -cap, VARIANTS.index((r["method"], r["estimator"])),
            int(r["level"]))


def write_results(path: Path, rows: list) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, FIELDS, lineterminator="\n",
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
        case = CASES_DIR / name
        row, series = analyse_case(case, entry["method"], entry["estimator"],
                                   entry["level"], entry["ranks"],
                                   entry.get("dt_cap", 0.0))
        rows.append(row)
        print(f"{name:20s} T/T_exact - 1 = {row['T_relerr']:+.4f}, "
              f"u_theta L2/(omega R1) = {row['ut_L2']:.3e}, "
              f"{'steady' if row['steady'] else 'NOT steady'}")
    if not rows:
        print("No completed cases: run ./Allverify first")
        return 1
    out = WORK_DIR / "results.csv"
    write_results(out, rows)
    print(f"Results: {out}")
    return 0


# --------------------------------------------------------------------- check
def check_analytic(p: dict) -> list:
    out = []
    closed, integrated = exact_rotor_torque(p), integrated_rotor_torque(p)
    ok = abs(closed - integrated) <= CHECK_EXACT_TORQUE * abs(closed)
    out.append(("PASS" if ok else "FAIL",
                f"exact torque {closed:.6e} N m, from the stress at the "
                f"rotor {integrated:.6e}"))
    a, b = exact_AB(p)
    wall = (exact_utheta(p, p["R1"]) - p["omega"] * p["R1"],
            exact_utheta(p, p["R2"]))
    ok = max(abs(x) for x in wall) < 1e-15
    out.append(("PASS" if ok else "FAIL",
                "u_theta(R1) = omega R1 and u_theta(R2) = 0"))
    ok = p["origin"] == p["CofR"] and p["axis"] == (0.0, 0.0, 1.0) and \
        p["omega"] > 0
    out.append(("PASS" if ok else "FAIL",
                "the torque axis is z, the centre of rotation of the body "
                "is the origin of the rotation, omega > 0"))
    # the field analysis measures r and the azimuthal direction from (0, 0)
    ok = tuple(p["origin"][:2]) == (0.0, 0.0)
    out.append(("PASS" if ok else "FAIL",
                "the rotation is centred at (0, 0), as the field analysis "
                f"assumes (origin {p['origin'][:2]})"))
    out.append(mesh_guard(p))
    return out


def mesh_guard(p: dict) -> tuple:
    """The field analysis uses SIDE and MESH_CELLS_ACROSS: they must be those
    of system/blockMeshDict (a square of side SIDE centred on the axis, one
    block of n x n x 1 cells with n = round(a pow(b, MESH_LEVEL - c)))"""
    problems = []
    for key in ("mesh_x", "mesh_y"):
        lo, hi = p[key]
        if abs(hi - lo - SIDE) > 1e-12 or abs(hi + lo) > 1e-12:
            problems.append(f"{key} = ({lo}, {hi}) is not a square of side "
                            f"{SIDE} centred on the axis")
    if not p.get("square_block"):
        problems.append("the block is not n x n x 1 cells")
    f = p.get("cells_formula")
    if f is None:
        problems.append("the cell count formula of blockMeshDict was not "
                        "recognised")
    else:
        a, b, c = f
        for level, n in MESH_CELLS_ACROSS.items():
            if round(a * b ** (level - c)) != n:
                problems.append(f"Level {level} has {round(a * b ** (level - c))}"
                                f" cells across, not {n}")
    return ("FAIL" if problems else "PASS",
            "the mesh of blockMeshDict is the one assumed by the field "
            f"analysis (side {SIDE}, cells across "
            f"{', '.join(str(n) for n in MESH_CELLS_ACROSS.values())})"
            + (": " + "; ".join(problems) if problems else ""))


def check_tracked(rows: list, p: dict) -> list:
    """Consistency of the tracked results with the exact solution and with
    the predefined criteria; nothing is compared with a measured value"""
    out = []
    exact = exact_rotor_torque(p)
    bad = [r for r in rows
           if abs(float(r["T_exact"]) / exact - 1.0) > 1e-5 or
           abs((float(r["T_mean"]) - exact) / exact - float(r["T_relerr"]))
           > 1e-5]
    out.append(("PASS" if not bad else "FAIL",
                f"{len(rows)} tracked rows: T_exact and T_relerr agree "
                "with the closed form"))
    flags = [r for r in rows if int(r["steady"]) != int(
        float(r["T_drift_rel"]) < STEADY_DRIFT and
        float(r["dU_5_10"]) < STEADY_DU)]
    out.append(("PASS" if not flags else "FAIL",
                "steady flags follow the predefined criteria"))
    unsteady = [r for r in rows if not int(r["steady"])]
    if unsteady:
        names = ", ".join(
            f"{r['method']}{'-' + r['estimator'] if r['estimator'] else ''} "
            f"L{r['level']}" + (f" cap {r['dt_cap']}" if r.get("dt_cap")
                                else "") for r in unsteady)
        out.append(("PASS", f"{len(unsteady)} tracked rows are not steady by "
                    "the predefined criteria and are reported with their mean "
                    f"and variance, not as converged: {names}"))
    resolved = [r for r in rows if r["T_field_rel"] != ""
                and r["method"] != "incremental" and r["estimator"] != "ST"]
    bad = [r for r in resolved
           if abs(float(r["T_field_rel"]) - float(r["T_relerr"])) >
           CHECK_FIELD_TORQUE and int(r["level"]) >= 2]
    out.append(("PASS" if not bad else "FAIL",
                "torque from the velocity field agrees with the force "
                f"output within {CHECK_FIELD_TORQUE:.2f} (Level >= 2; not "
                "applied to the incremental method and the surface "
                "traction, whose differences of 2-4 points are reported)"))
    return out


def check_case(row: dict, tracked: dict) -> list:
    out = []
    d = abs(float(row["T_relerr"]) - float(tracked["T_relerr"]))
    out.append(("PASS" if d <= CHECK_TORQUE else "FAIL",
                f"{row['method']} {row['estimator']} L{row['level']}: torque "
                f"error {float(row['T_relerr']):+.4f} against tracked "
                f"{float(tracked['T_relerr']):+.4f}"))
    d = abs(float(row["ut_L2"]) - float(tracked["ut_L2"]))
    out.append(("PASS" if d <= CHECK_VELOCITY else "FAIL",
                f"{row['method']} {row['estimator']} L{row['level']}: "
                f"velocity L2 {float(row['ut_L2']):.3e} against tracked "
                f"{float(tracked['ut_L2']):.3e}"))
    drift = float(tracked["T_drift_rel"])
    if 0.5 * STEADY_DRIFT < drift < 2.0 * STEADY_DRIFT:
        out.append(("SKIPPED", "steady-state flag: the tracked drift is within "
                    "a factor 2 of the threshold, the endpoint criterion is "
                    "too noisy to compare"))
    else:
        out.append(("PASS" if int(row["steady"]) == int(tracked["steady"])
                    else "FAIL", "steady-state flag"))
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
    for name, entry in done:
        row, _ = analyse_case(CASES_DIR / name, entry["method"],
                              entry["estimator"], entry["level"],
                              entry["ranks"], entry.get("dt_cap", 0.0))
        ref = {row_key(r): r for r in tracked + read_results(
            RESULTS_DIR / "controls" / "cases-timestep.csv")}.get(row_key(row))
        if ref is None:
            results.append(("SKIPPED", f"{name}: no tracked row"))
        else:
            results += check_case(row, ref)
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
    return {"ghostCell": "ghostCell (exploratory)"}.get(r["method"], r["method"])


def master_table(rows: list) -> str:
    """The tutorial configuration (Courant-limited time step) only"""
    rows = [r for r in rows if row_key(r)[3] == 0.0]
    head = ("| Method | Level | Cells | Rotor torque (N m) | Exact (N m) | "
            "Torque difference | u_theta L1 | u_theta L2 | u_theta Linf | "
            "Torque std. (relative) | Steady | Ranks | Rank-h |\n"
            "| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- "
            "| --- | --- | --- |\n")
    lines = []
    for r in sorted(rows, key=sort_key):
        rh = f"{float(r['rank_h']):.3f}" if r.get("rank_h") else "n/a"
        lines.append(
            f"| {method_label(r)} | {r['level']} | {r['cells']} | "
            f"{float(r['T_mean']):.4e} | {float(r['T_exact']):.4e} | "
            f"{100 * float(r['T_relerr']):+.2f} % | "
            f"{float(r['ut_L1']):.2e} | {float(r['ut_L2']):.2e} | "
            f"{float(r['ut_Linf']):.2e} | {float(r['T_std_rel']):.1e} | "
            f"{'yes' if int(r['steady']) else 'no'} | {r['ranks']} | {rh} |")
    return head + "\n".join(lines) + "\n"


def time_step_table(rows: list) -> str:
    """Torque difference from the exact torque and azimuthal velocity L2
    norm at the tutorial time step (Courant 0.5, the step halves with each
    level) and with the step capped at 0.01, 0.005 and 0.0025 s"""
    caps = [0.0, 0.01, 0.005, 0.0025]
    head = ("| Method | Level | tutorial step | cap 0.01 s | cap 0.005 s | "
            "cap 0.0025 s |\n| --- | --- | --- | --- | --- | --- |\n")
    lines = []
    for method, est in DT_VARIANTS:
        for level in (1, 2, 3):
            cells = []
            for cap in caps:
                hit = [r for r in rows if (r["method"], r["estimator"]) ==
                       (method, est) and int(r["level"]) == level
                       and row_key(r)[3] == cap]
                cells.append(
                    f"{100 * float(hit[0]['T_relerr']):+.2f} % / "
                    f"{float(hit[0]['ut_L2']):.2e}"
                    f"{'' if int(hit[0]['steady']) else ' *'}"
                    if hit else "n/a")
            lines.append(f"| {method_label({'method': method, 'estimator': est})}"
                         f" | {level} | " + " | ".join(cells) + " |")
    return head + "\n".join(lines) + "\n"


def cmd_tables(args) -> int:
    rows = read_results()
    if not rows:
        fail("no tracked results")
    print(master_table(rows))
    print(time_step_table(rows))
    return 0


# ---------------------------------------------------------------------- unit
def cmd_unit(args) -> int:
    import unittest
    here = Path(__file__).resolve().parent
    sys.path.insert(0, str(here))
    suite = unittest.defaultTestLoader.discover(
        str(here), pattern="test_immersed_taylor_couette_verification.py")
    result = unittest.TextTestRunner(verbosity=1, warnings="error").run(suite)
    return 0 if result.wasSuccessful() else 1


COMMANDS = ("run", "generate", "dry-run", "status", "analyse", "check",
            "tables", "unit")


def normalise_argv(argv: list) -> list:
    """'run' is the default command: it is inserted before parsing when the
    arguments are empty or do not start with a command or -h/--help, so that
    ./Allverify --runs L1 L2 is ./Allverify run --runs L1 L2"""
    argv = list(argv)
    if not argv or (argv[0] not in COMMANDS and argv[0] not in ("-h", "--help")):
        return ["run"] + argv
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
        s.add_argument("--force", action="store_true",
                       help="regenerate existing cases")
        s.add_argument("--end-time", type=float, default=0.0, metavar="T",
                       help="write endTime T into controlDict (default: the "
                       "10 s of the tutorial) and tag the case name with -tT")
    for name in COMMANDS[3:]:
        sub.add_parser(name)
    args = parser.parse_args(argv)
    command = args.command
    handlers = {"run": cmd_run, "generate": cmd_generate,
                "dry-run": cmd_dry_run, "status": cmd_status,
                "analyse": cmd_analyse, "check": cmd_check,
                "tables": cmd_tables, "unit": cmd_unit}
    return handlers[command](args)


if __name__ == "__main__":
    sys.exit(main())
