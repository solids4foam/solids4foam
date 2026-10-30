#!/usr/bin/env python3
"""Verification of the oscillatingCylinderInChannel tutorial: the cylinder
oscillating in a closed channel with the cutLink, ghostCell, penalty and
legacy incremental methods of immersedBoundaryForce, compared with two moving
body-fitted mesh solutions and, as a secondary comparison, the result of Wan
and Turek (2007) digitised from a figure of Erzincanli and Sahin (2013).

Each case is a copy of the tutorial under verification/work/cases/, with
only the method settings, the mesh level, a fixed time step, the end time,
the outer correctors and the decomposition changed. Serial cases are run
with the tutorial's Allrun; parallel cases with the solids4foam run
functions. The tutorial itself is never modified. A manifest,
work/manifest.json, records every case, its settings, its SLURM job (if
submitted) and its state, so that a campaign can be resumed.

The cutLink surface-traction estimator is run with the quadrature of the
tracked results, the greedy thinning of the candidate points, always written
explicitly (tractionQuadrature greedy), so that a change of the solver
default does not change the verification silently (see
--traction-quadrature).

Usage:
    ./Allverify [run] [--runs ...]   generate, run locally, analyse, check
    ./Allverify generate --runs ...  generate the cases only
    ./Allverify dry-run --runs ...   list the cases and their settings
    ./Allverify submit --runs ...    optional: one SLURM job per case
    ./Allverify status               state of the cases in the manifest
    ./Allverify analyse              analyse the completed cases
    ./Allverify check                compare them with the tracked results
    ./Allverify export               compact results to reference/results/
    ./Allverify plot                 figures from the results

The default run is the cutLink Level 1 case (about 1.5 minutes on one
core); './Allverify --runs smoke' runs the first 0.5 s only.
"""

import argparse
import csv
import hashlib
import json
import math
import os
import re
import shlex
import shutil
import subprocess
import sys
import time
from datetime import datetime
from pathlib import Path

VERIFICATION_DIR = Path(__file__).resolve().parent.parent
TUTORIAL_DIR = VERIFICATION_DIR.parent
WORK_DIR = VERIFICATION_DIR / "work"
CASES_DIR = WORK_DIR / "cases"
MANIFEST = WORK_DIR / "manifest.json"
OUTPUT_DIR = VERIFICATION_DIR / "postProcessing"
REFERENCE_DIR = VERIFICATION_DIR / "reference"

METHODS = ("cutLink", "ghostCell", "penalty", "incremental")
LEVELS = (1, 2, 3, 4)

# Common settings of the verification cases
END_TIME = 8.0
WRITE_INTERVAL = 0.5
N_CORRECTORS = 3
N_OUTER = {"cutLink": 1, "ghostCell": 1, "penalty": 1, "incremental": 5}
COUPLING_COEFF = 0.2

# Time step of each level, halved with the cell size: a body Courant number
# 2 pi A/T dt/h of about 0.1
LEVEL_DT = {1: 0.0025, 2: 0.00125, 3: 0.000625, 4: 0.0003125}

# Default MPI ranks of each level
LEVEL_RANKS = {1: 1, 2: 2, 3: 8, 4: 8}

# Optional SLURM requests of each level
SLURM_MEMORY = {1: "2G", 2: "4G", 3: "8G", 4: "16G"}
SLURM_TIME = {1: "06:00:00", 2: "12:00:00", 3: "24:00:00", 4: "47:00:00"}

# The runs of the campaign: name -> (method, level, time step, end time
# [, MPI ranks [, outer correctors]]); a run with its own ranks or outer
# correctors has them in its case name. R1-R3, R5: the level time step;
# R4: Level 2 at half the time step; N: outer-corrector sensitivity of the
# legacy method over the first 2.25 s; O4: serial run of R2A
RUNS = {
    "smoke": ("cutLink", 1, 0.0025, 0.5),
    "R1A": ("cutLink", 1, 0.0025, END_TIME),
    "R1B": ("penalty", 1, 0.0025, END_TIME),
    "R1C": ("incremental", 1, 0.0025, END_TIME),
    "R1D": ("ghostCell", 1, 0.0025, END_TIME),
    "R2A": ("cutLink", 2, 0.00125, END_TIME),
    "R2B": ("penalty", 2, 0.00125, END_TIME),
    "R2C": ("incremental", 2, 0.00125, END_TIME),
    "R2D": ("ghostCell", 2, 0.00125, END_TIME),
    "R3A": ("cutLink", 3, 0.000625, END_TIME),
    "R3B": ("penalty", 3, 0.000625, END_TIME),
    "R3C": ("incremental", 3, 0.000625, END_TIME),
    "R3D": ("ghostCell", 3, 0.000625, END_TIME),
    "R4A": ("cutLink", 2, 0.000625, END_TIME),
    "R4B": ("penalty", 2, 0.000625, END_TIME),
    "R4C": ("incremental", 2, 0.000625, END_TIME),
    "R4D": ("ghostCell", 2, 0.000625, END_TIME),
    "R5A": ("cutLink", 4, 0.0003125, END_TIME),
    "R5B": ("penalty", 4, 0.0003125, END_TIME),
    "R5D": ("ghostCell", 4, 0.0003125, END_TIME),
    "N5": ("incremental", 2, 0.00125, 2.25, 2, 5),
    "N10": ("incremental", 2, 0.00125, 2.25, 2, 10),
    "N20": ("incremental", 2, 0.00125, 2.25, 2, 20),
    "O4": ("cutLink", 2, 0.00125, END_TIME, 1),
}
DEFAULT_RUNS = ("R1A",)

# Analysis window (s), excluding the impulsive start
WINDOW = (0.25, 7.5)

# Settings of each method in constant/fvOptions, those of the static-cylinder
# verification (solids4foam pull request #595), written explicitly, and of the
# body. They are not all source defaults: couplingCoeff 0.2 is the value of
# #595, the default of immersedBoundaryForce is 0.8
METHOD_SETTINGS = {
    "cutLink": (
        [
            ("penaltyCoeff", "1e3"),
            ("pinnedRateCoeff", "100"),
            ("linkCorrection", "yes"),
            ("apertureCoupling", "yes"),
            ("forceEstimator", "surfaceTraction"),
        ],
        [],
    ),
    "ghostCell": (
        [
            ("penaltyCoeff", "1e3"),
            ("imageDistance", "1.5"),
            ("imageInterpolation", "cellPoint"),
        ],
        [],
    ),
    "penalty": (
        [
            ("penaltyCoeff", "1e3"),
            ("weighting", "volumeFraction"),
            ("surfaceRateCoeff", "3"),
        ],
        [("occupancy", "signedDistance")],
    ),
    "incremental": (
        [("couplingCoeff", "0.2")],
        [("occupancy", "vertexFraction")],
    ),
}

FOAM_HEADER = """\
/*--------------------------------*- C++ -*----------------------------------*\\
| solids4foam: oscillatingCylinderInChannel verification case                 |
\\*---------------------------------------------------------------------------*/
FoamFile
{{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "{location}";
    object      {object};
}}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

"""

MESH_FILES = ("points", "faces", "owner", "neighbour", "boundary")


def fail(message: str) -> None:
    raise RuntimeError(message)


def fmt(value: float) -> str:
    return f"{value:.12g}"


# -----------------------------------------------------------------------------
# Physical parameters, read from the tutorial (or a case)
# -----------------------------------------------------------------------------

def dict_scalar(text: str, key: str) -> float:
    """Last value of 'key ... value;' (dimensioned or not) in a dictionary"""
    matches = re.findall(rf"^\s*{key}\s+(?:{key}\s+\[[^\]]*\]\s+)?"
                         r"([-+0-9.eE]+)\s*;", text, re.MULTILINE)
    if not matches:
        fail(f"no '{key}' entry")
    return float(matches[-1])


def stl_vertices(path: Path) -> list:
    return [tuple(map(float, m)) for m in re.findall(
        r"vertex\s+(\S+)\s+(\S+)\s+(\S+)", path.read_text())]


def parameters(case: Path) -> dict:
    """Density, viscosity, motion, diameter and thickness of a case, and the
    force coefficient scale 2/(rho Uref^2 D Lz) with Uref = 2 pi A/T"""
    transport = (case / "constant/transportProperties").read_text()
    fv_options = (case / "constant/fvOptions").read_text()
    motion = fv_options[fv_options.index("motion"):]
    amplitude = dict_scalar(motion, "amplitude")
    period = dict_scalar(motion, "period")

    vertices = stl_vertices(case / "constant/triSurface/cylinder.stl")
    xc = sum(v[0] for v in vertices)/len(vertices)
    yc = sum(v[1] for v in vertices)/len(vertices)
    radius = max(math.hypot(v[0] - xc, v[1] - yc) for v in vertices)

    block = (case / "system/blockMeshDict").read_text()
    zs = [float(m) for m in re.findall(
        r"^\s*\(\s*[-+0-9.eE]+\s+[-+0-9.eE]+\s+([-+0-9.eE]+)\s*\)", block,
        re.MULTILINE)]
    thickness = max(zs) - min(zs)

    rho = dict_scalar(transport, "rho")
    u_ref = 2*math.pi*amplitude/period
    return {
        "rho": rho, "nu": dict_scalar(transport, "nu"),
        "amplitude": amplitude, "period": period, "U_ref": u_ref,
        "centre": (xc, yc), "D": 2*radius, "Lz": thickness,
        "coeff_scale": 2.0/(rho*u_ref**2*2*radius*thickness),
    }


def body_centre(params: dict, t: float) -> float:
    """x of the cylinder centre at time t (the STL is at the zero
    displacement position)"""
    return (params["centre"][0]
            + params["amplitude"]*math.sin(2*math.pi*t/params["period"]))


# -----------------------------------------------------------------------------
# Manifest
# -----------------------------------------------------------------------------

def load_manifest() -> dict:
    if MANIFEST.exists():
        return json.loads(MANIFEST.read_text())
    return {"created": datetime.now().isoformat(timespec="seconds"),
            "cases": {}}


def save_manifest(manifest: dict) -> None:
    WORK_DIR.mkdir(parents=True, exist_ok=True)
    tmp = MANIFEST.with_suffix(".json.tmp")
    tmp.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    tmp.replace(MANIFEST)


def case_name(method: str, level: int, dt: float, end: float,
              ranks: int = 0, n_outer: int = 0, quadrature: str = "") -> str:
    name = f"{method}-L{level}-dt{fmt(dt).replace('.', 'p')}"
    if end != END_TIME:
        name += f"-t{fmt(end).replace('.', 'p')}"
    if ranks:
        name += f"-np{ranks}"
    if n_outer and n_outer != N_OUTER[method]:
        name += f"-nOuter{n_outer}"
    if quadrature:
        name += f"-{quadrature}"
    return name


def selected(args) -> list:
    """(run, method, level, dt, end, name) of the selected runs"""
    out = []
    for run in args.runs:
        if run not in RUNS:
            fail(f"unknown run {run}; choose from {', '.join(RUNS)}")
        method, level, dt, end = RUNS[run][:4]
        tag = quadrature_tag(method, getattr(args, "traction_quadrature",
                                             HISTORICAL_QUADRATURE))
        out.append((run, method, level, dt, end,
                    case_name(method, level, dt, end, run_ranks(run),
                              run_n_outer(run), tag)))
    return out


def run_ranks(run: str) -> int:
    """The MPI ranks fixed by a run, 0 if none"""
    return RUNS[run][4] if len(RUNS[run]) > 4 else 0


def run_n_outer(run: str) -> int:
    """The outer correctors of a run: its own, or those of its method"""
    return RUNS[run][5] if len(RUNS[run]) > 5 else N_OUTER[RUNS[run][0]]


# -----------------------------------------------------------------------------
# Surface-traction quadrature of the cutLink estimator
# -----------------------------------------------------------------------------

# The quadrature of the tracked results: the greedy thinning of the
# candidate points. The cutLink cases always select it explicitly. A solver
# without the tractionQuadrature keyword ignores the entry, and its only
# quadrature is greedy; a solver with the keyword and a different default
# (clustered) runs greedy because of the entry; a solver that does not know
# the chosen name stops with an error
HISTORICAL_QUADRATURE = "greedy"
QUADRATURES = ("greedy", "none", "clustered")


def quadrature_tag(method: str, choice: str) -> str:
    """Case-name tag of a cutLink case run with a quadrature chosen on the
    command line; '' for the historical one"""
    if method != "cutLink" or choice == HISTORICAL_QUADRATURE:
        return ""
    return choice


def ranks_for(level: int, args, run: str = "") -> int:
    if run and run_ranks(run):
        return run_ranks(run)
    return int(args.ranks) if getattr(args, "ranks", None) else \
        LEVEL_RANKS[level]


# -----------------------------------------------------------------------------
# Case generation
# -----------------------------------------------------------------------------

def replace_entry(path: Path, key: str, value: str) -> None:
    text = path.read_text()
    pattern = re.compile(rf"^({re.escape(key)}\s+)[^;\n]*;", re.MULTILINE)
    new, count = pattern.subn(lambda m: f"{m.group(1)}{value};", text)
    if count != 1:
        fail(f"found {count} '{key}' entries in {path}")
    path.write_text(new)


def set_pimple_entry(path: Path, key: str, value: int) -> None:
    text = path.read_text()
    pattern = re.compile(rf"^(\s+{key}\s+)\d+;", re.MULTILINE)
    new, count = pattern.subn(lambda m: f"{m.group(1)}{value};", text)
    if count != 1:
        fail(f"found {count} '{key}' entries in {path}")
    path.write_text(new)


def write_fv_options(path: Path, method: str, motion: dict,
                     quadrature: str = "") -> None:
    option, body = METHOD_SETTINGS[method]
    if method == "cutLink" and quadrature:
        option = option + [("tractionQuadrature", quadrature)]
    width = max(len(k) for k, _ in option + [("method", "")]) + 4
    lines = [FOAM_HEADER.format(location="constant", object="fvOptions"),
             "// Generated by verification/scripts/"
             "oscillating_cylinder_verification.py\n",
             "immersedBoundary\n{\n",
             f"    {'type'.ljust(width)}immersedBoundaryForce;\n\n",
             f"    {'method'.ljust(width)}{method};\n"]
    lines += [f"    {k.ljust(width)}{v};\n" for k, v in option]
    lines.append("\n    bodies\n    {\n        cylinder\n        {\n")
    lines.append(f"            {'surface'.ljust(12)}\"cylinder.stl\";\n")
    lines += [f"            {k.ljust(12)}{v};\n" for k, v in body]
    lines.append("\n            motion\n            {\n"
                 "                type        sinusoidalTranslation;\n"
                 f"                amplitude   {fmt(motion['amplitude'])};\n"
                 f"                period      {fmt(motion['period'])};\n"
                 "                direction   (1 0 0);\n"
                 "            }\n        }\n    }\n}\n\n")
    lines.append("// ****************************************************"
                 "********************* //\n")
    path.write_text("".join(lines))


def write_decompose_par_dict(path: Path, ranks: int) -> None:
    text = FOAM_HEADER.format(location="system", object="decomposeParDict")
    text += (f"numberOfSubdomains  {ranks};\n\nmethod          simple;\n\n"
             f"simpleCoeffs\n{{\n    n           ({ranks} 1 1);\n}}\n\n"
             "// **************************************************************"
             "*********** //\n")
    path.write_text(text)


RUN_SCRIPT = """\
#!/bin/bash
# Generated by verification/scripts/oscillating_cylinder_verification.py
cd "${{0%/*}}" || exit 1
export OMP_NUM_THREADS=1
start=$(date +%s)
{body}
status=$?
# The validity of a run is judged from the solver log; the status of Allrun
# is that of its last command (the optional plot)
echo "wall_s=$(( $(date +%s) - start )) allrun_status=$status" >> wallclock.txt
exit $status
"""

SERIAL_BODY = "./Allrun"

PARALLEL_BODY = """\
(
    # As in the tutorial's Allrun, with decomposePar and runParallel
    . "$WM_PROJECT_DIR/bin/tools/RunFunctions"
    source solids4FoamScripts.sh
    solids4Foam::caseDoesNotRunWithFoamExtend
    solids4Foam::caseDoesNotRunWithOpenFOAMOrg
    solids4Foam::runApplication blockMesh || exit 1
    solids4Foam::runApplication decomposePar || exit 1
    solids4Foam::runParallel solids4Foam
)"""


def git_sha() -> str:
    """Revision of the tutorial tree: the sources of the cases, not
    necessarily the build of the solver that runs them"""
    return subprocess.run(["git", "-C", str(TUTORIAL_DIR), "rev-parse",
                           "HEAD"], capture_output=True, text=True,
                          check=False).stdout.strip()


# The solver of the tracked results. They record no per-case provenance:
# this is the build record of the runs that produced them, stated in
# README.md, and is not written into any case
BASELINE = {"openfoam": "v2412", "solids4foam": "f33b244a2"}
PROVENANCE_COLUMNS = ("openfoam_version", "openfoam_variant", "solver_build",
                      "solver_binary_sha256", "tutorial_revision")


def openfoam_variant(version: str, fork: str = "") -> str:
    """OpenFOAM variant from WM_PROJECT_VERSION (and WM_FORK): v2512 is
    OpenFOAM.com, a number such as 12 the OpenFOAM Foundation; '' if
    unknown"""
    if fork:
        return fork
    if re.fullmatch(r"v\d{4}", version):
        return "OpenFOAM.com"
    if re.fullmatch(r"\d+(\.\d+)?", version):
        return "openfoam.org"
    return ""


def solver_build(executable: str = "") -> dict:
    """What is known of the solids4foam build that runs the cases: the
    executable, its SHA-256 and, if a BUILD_INFO.txt with a 'SHA=' line is in
    the root of its build (<root>/platforms/<options>/bin/solids4Foam), that
    revision. Nothing is guessed: an unknown entry is ''"""
    executable = executable or shutil.which("solids4Foam") or ""
    out = {"solver_executable": "", "solver_binary_sha256": "",
           "solver_build": ""}
    if not executable or not Path(executable).is_file():
        return out
    path = Path(executable).resolve()
    out["solver_executable"] = str(path)
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    out["solver_binary_sha256"] = digest.hexdigest()
    info = path.parents[3] / "BUILD_INFO.txt" if len(path.parents) > 3 \
        else Path()
    if info.is_file():
        for line in info.read_text(errors="replace").splitlines():
            if line.startswith("SHA="):
                out["solver_build"] = line[4:].strip()
                break
    return out


def provenance(executable: str = "") -> dict:
    """OpenFOAM version and variant, solids4foam build and tutorial revision
    of the environment that generates a case"""
    version = os.environ.get("WM_PROJECT_VERSION", "")
    out = {"openfoam_version": version,
           "openfoam_variant": openfoam_variant(
               version, os.environ.get("WM_FORK", ""))}
    out.update(solver_build(executable))
    out["tutorial_revision"] = git_sha()
    return out


def generate_case(run: str, method: str, level: int, dt: float, end: float,
                  name: str, ranks: int, force: bool,
                  quadrature_choice: str = HISTORICAL_QUADRATURE) -> Path:
    case = CASES_DIR / name
    if case.exists():
        if not force:
            fail(f"{case} exists (use --force to regenerate)")
        shutil.rmtree(case)
    case.mkdir(parents=True)
    for item in ("0", "constant", "system", "Allrun"):
        source = TUTORIAL_DIR / item
        if source.is_dir():
            shutil.copytree(source, case / item, symlinks=True)
        else:
            shutil.copy2(source, case / item)

    replace_entry(case / "system/blockMeshDict", "MESH_LEVEL", str(level))
    control = case / "system/controlDict"
    replace_entry(control, "deltaT", fmt(dt))
    replace_entry(control, "adjustTimeStep", "no")
    replace_entry(control, "endTime", fmt(end))
    replace_entry(control, "writeControl", "runTime")
    replace_entry(control, "writeInterval", fmt(WRITE_INTERVAL))
    # Time names exact to the time step, as they are the force file times
    replace_entry(control, "timePrecision", "10")

    motion = parameters(TUTORIAL_DIR)
    quadrature = quadrature_choice if method == "cutLink" else ""
    write_fv_options(case / "constant/fvOptions", method, motion,
                     quadrature)
    fv_solution = case / "system/fvSolution"
    n_outer = run_n_outer(run) if run in RUNS else N_OUTER[method]
    set_pimple_entry(fv_solution, "nOuterCorrectors", n_outer)
    set_pimple_entry(fv_solution, "nCorrectors", N_CORRECTORS)
    write_decompose_par_dict(case / "system/decomposeParDict", ranks)

    script = case / "Allrun.verify"
    script.write_text(RUN_SCRIPT.format(
        body=SERIAL_BODY if ranks == 1 else PARALLEL_BODY))
    script.chmod(0o755)

    meta = {
        "run": run, "case": name, "method": method, "level": level,
        "cells_per_D": 10*2**(level - 1), "deltaT": dt, "end_time": end,
        "write_interval": WRITE_INTERVAL,
        "nOuterCorrectors": n_outer, "nCorrectors": N_CORRECTORS,
        "ranks": ranks, "decomposition": f"simple ({ranks} 1 1)",
        "fvOptions": dict(METHOD_SETTINGS[method][0]
                          + ([("tractionQuadrature", quadrature)]
                             if quadrature else [])),
        "traction_quadrature": quadrature,
        "body": dict(METHOD_SETTINGS[method][1]),
        "git_sha": git_sha(),
        "provenance": provenance(),
        "tracked_baseline": BASELINE,
        "generated": datetime.now().isoformat(timespec="seconds"),
    }
    (case / "case.json").write_text(json.dumps(meta, indent=2) + "\n")
    return case


# -----------------------------------------------------------------------------
# Case state
# -----------------------------------------------------------------------------

def solver_problem(case: Path) -> str:
    """Empty if the solver log shows a completed run, else the reason"""
    log = case / "log.solids4Foam"
    if not log.is_file():
        return "no log.solids4Foam"
    text = log.read_text(errors="replace")
    # (not the "Floating point exception trapping enabled" banner)
    if re.search(r"FOAM FATAL|FOAM aborting|sigFpe::sigHandler|"
                 r"Floating point exception \(core dumped\)|^ERROR$", text,
                 re.MULTILINE):
        return "the solver failed"
    if not re.search(r"^\s*(\[\d+,\d+\]<stdout>:)?End\s*$", text,
                     re.MULTILINE):
        return "the run did not complete"
    return ""


def queue_jobs() -> dict:
    """{job id: (state, cpus)} of the user's SLURM jobs ({} without SLURM)"""
    if shutil.which("squeue") is None:
        return {}
    out = subprocess.run(["squeue", "--me", "--noheader",
                          "--format=%i %T %C"], capture_output=True,
                         text=True, check=False).stdout
    jobs = {}
    for line in out.split("\n"):
        parts = line.split()
        if len(parts) == 3 and parts[2].isdigit():
            jobs[parts[0]] = (parts[1], int(parts[2]))
    return jobs


def case_state(name: str, entry: dict, jobs: dict) -> str:
    case = CASES_DIR / name
    job = str(entry.get("job_id", ""))
    if job and job in jobs:
        return jobs[job][0].lower()
    if not case.exists():
        return "missing"
    problem = solver_problem(case)
    if not problem:
        return "completed"
    if (case / "log.solids4Foam").exists():
        return "failed" if "failed" in problem else (
            "incomplete" if job or entry.get("ran_locally") else "running?")
    return "submitted?" if job else "generated"


# -----------------------------------------------------------------------------
# Commands: generate, dry-run, run, submit, status
# -----------------------------------------------------------------------------

def plan(args) -> list:
    rows = selected(args)
    print(f"{'run':5s} {'case':34s} {'ranks':>5s} {'dt':>10s} {'end':>5s}")
    for run, method, level, dt, end, name in rows:
        print(f"{run:5s} {name:34s} {ranks_for(level, args, run):5d} "
              f"{fmt(dt):>10s} {fmt(end):>5s}")
    return rows


def cmd_generate(args) -> int:
    manifest = load_manifest()
    for run, method, level, dt, end, name in plan(args):
        entry = manifest["cases"].get(name, {})
        if entry.get("job_id") and not args.force:
            fail(f"{name} has been submitted (job {entry['job_id']}); "
                 "use --force to regenerate it")
        generate_case(run, method, level, dt, end, name,
                      ranks_for(level, args, run), args.force,
                      getattr(args, "traction_quadrature",
                              HISTORICAL_QUADRATURE))
        manifest["cases"][name] = {
            "run": run, "method": method, "level": level, "deltaT": dt,
            "end_time": end, "ranks": ranks_for(level, args, run),
            "generated": datetime.now().isoformat(timespec="seconds")}
        print(f"{name}: generated in {CASES_DIR / name}")
    save_manifest(manifest)
    return 0


def cmd_dry_run(args) -> int:
    plan(args)
    params = parameters(TUTORIAL_DIR)
    print(f"U_ref = {params['U_ref']:.10g} m/s, D = {params['D']:.6g} m, "
          f"Lz = {params['Lz']:.6g} m, coefficient scale = "
          f"{params['coeff_scale']:.6g}")
    print("Dry run: nothing generated, run or submitted")
    return 0


def check_environment() -> None:
    for executable in ("blockMesh", "solids4Foam", "solids4FoamScripts.sh"):
        if shutil.which(executable) is None:
            fail(f"'{executable}' is unavailable: source OpenFOAM and "
                 "solids4foam first")


def cmd_run(args) -> int:
    if args.dry_run:
        return cmd_dry_run(args)
    check_environment()
    manifest = load_manifest()
    for run, method, level, dt, end, name in selected(args):
        case = CASES_DIR / name
        if case.exists() and not solver_problem(case) and not args.force:
            print(f"{name}: completed run reused")
            continue
        generate_case(run, method, level, dt, end, name,
                      ranks_for(level, args, run), True,
                      getattr(args, "traction_quadrature",
                              HISTORICAL_QUADRATURE))
        manifest["cases"][name] = {
            "run": run, "method": method, "level": level, "deltaT": dt,
            "end_time": end, "ranks": ranks_for(level, args, run),
            "ran_locally": True,
            "started": datetime.now().isoformat(timespec="seconds")}
        save_manifest(manifest)
        print(f"{name}: running on {ranks_for(level, args, run)} rank(s)",
              flush=True)
        with (case / "log.Allverify").open("w") as handle:
            subprocess.run(["./Allrun.verify"], cwd=case, stdout=handle,
                           stderr=subprocess.STDOUT, check=False)
        problem = solver_problem(case)
        manifest["cases"][name]["finished"] = \
            datetime.now().isoformat(timespec="seconds")
        save_manifest(manifest)
        if problem:
            fail(f"{name}: {problem}; see {case / 'log.solids4Foam'}")
    status = cmd_analyse(args)
    return cmd_check(args) or status


def cmd_submit(args) -> int:
    """Optional: one SLURM job per generated case, running its
    Allrun.verify; only a summary without --yes"""
    manifest = load_manifest()
    jobs = queue_jobs()
    queued_cpus = sum(c for _, c in jobs.values())
    todo = []
    for run, method, level, dt, end, name in selected(args):
        case = CASES_DIR / name
        entry = manifest["cases"].get(name)
        if entry is None or not (case / "Allrun.verify").exists():
            fail(f"{name} has not been generated (./Allverify generate)")
        state = case_state(name, entry, jobs)
        if state in ("pending", "running", "completed", "submitted?",
                     "running?") and not args.force:
            print(f"{name}: {state}, not submitted")
            continue
        todo.append((name, entry))
    requested = sum(e["ranks"] for _, e in todo)
    print(f"{len(todo)} job(s), {requested} CPU(s); {queued_cpus} CPU(s) "
          "already queued by the user")
    if args.max_cpus and queued_cpus + requested > args.max_cpus:
        fail(f"{queued_cpus} + {requested} CPUs would exceed --max-cpus "
             f"{args.max_cpus}")
    if not todo:
        return 0
    for name, entry in todo:
        level = entry["level"]
        print(f"  {name}: {entry['ranks']} rank(s), "
              f"{args.mem or SLURM_MEMORY[level]}, "
              f"{args.time or SLURM_TIME[level]}")
    if not args.yes:
        print("Not submitted: re-run with --yes to submit")
        return 0
    if shutil.which("sbatch") is None:
        fail("sbatch is unavailable: run the cases locally with ./Allverify")
    setup = (f"source {shlex.quote(args.env_script)} && "
             if args.env_script else "")
    for name, entry in todo:
        case = CASES_DIR / name
        level = entry["level"]
        command = (f"{setup}cd {shlex.quote(str(case))} && "
                   "./Allrun.verify > log.Allverify 2>&1")
        sbatch = ["sbatch", "--parsable", f"--job-name=ocv-{name}",
                  f"--ntasks={entry['ranks']}", "--cpus-per-task=1",
                  f"--mem={args.mem or SLURM_MEMORY[level]}",
                  f"--time={args.time or SLURM_TIME[level]}",
                  f"--output={case / 'log.slurm'}", f"--chdir={case}"]
        if args.partition:
            sbatch.append(f"--partition={args.partition}")
        sbatch += shlex.split(args.sbatch_args)
        sbatch += ["--wrap", f"bash -c {shlex.quote(command)}"]
        job = subprocess.run(sbatch, capture_output=True, text=True,
                             check=True).stdout.strip().split(";")[0]
        entry.update({"job_id": job, "partition": args.partition,
                      "submitted": datetime.now().isoformat(
                          timespec="seconds")})
        save_manifest(manifest)
        print(f"{name}: SLURM job {job}")
    return 0


def last_time(case: Path) -> float:
    rows = force_rows(case)
    return rows[-1][0] if rows else math.nan


def cmd_status(args) -> int:
    manifest = load_manifest()
    jobs = queue_jobs()
    print(f"{'case':34s} {'run':5s} {'ranks':>5s} {'job':>8s} "
          f"{'state':>11s} {'t':>7s}")
    for name, entry in sorted(manifest["cases"].items(),
                              key=lambda kv: kv[1].get("run", "")):
        state = case_state(name, entry, jobs)
        t = last_time(CASES_DIR / name) if (CASES_DIR / name).exists() \
            else math.nan
        print(f"{name:34s} {entry.get('run', ''):5s} "
              f"{entry.get('ranks', 0):5d} {str(entry.get('job_id', '')):>8s} "
              f"{state:>11s} {t:7.3f}")
    queued = sum(c for _, c in jobs.values())
    print(f"User SLURM jobs: {len(jobs)}, {queued} CPU(s)")
    return 0


# -----------------------------------------------------------------------------
# Force files
# -----------------------------------------------------------------------------

def parse_force_header(line: str) -> dict:
    """Column indices (0-based) of the vectors of a force file header line
    '# Time<TAB>force (x y z)<TAB>torque (x y z)<TAB>...'"""
    names = [n.strip() for n in line.lstrip("#").split("\t") if n.strip()]
    if not names or names[0] != "Time":
        fail(f"unexpected force file header: {line.strip()}")
    columns, index = {}, 1
    for title in names[1:]:
        key = {"force (x y z)": "force", "torque (x y z)": "torque",
               "inertia of the fluid inside (x y z)": "inertia",
               "force from the momentum exchange (x y z)":
                   "momentumExchange",
               "force from the surface traction (x y z)":
                   "surfaceTraction"}.get(title)
        if key is None:
            fail(f"unknown force file column '{title}'")
        columns[key] = index
        index += 3
    return columns


def read_force_file(path: Path) -> tuple:
    """(columns, rows) of a force file; rows are lists of floats"""
    columns, rows = None, []
    for line in path.read_text().splitlines():
        if line.startswith("# Time"):
            columns = parse_force_header(line)
        elif line.strip() and not line.startswith("#"):
            values = line.split()
            try:
                rows.append([float(v) for v in values])
            except ValueError:
                continue
    if columns is None:
        fail(f"{path}: no header")
    width = 1 + 3*len(columns)
    return columns, [r for r in rows if len(r) == width]


def force_files(case: Path) -> list:
    base = case / "postProcessing/immersedBoundary"
    files = [d / "cylinder.dat" for d in base.iterdir()
             if (d / "cylinder.dat").exists()] if base.exists() else []
    return sorted(files, key=lambda f: float(f.parent.name))


def merge_segments(segments: list) -> tuple:
    """Merge (start time, columns, rows) segments of a run and its restarts:
    a later segment replaces the rows from its start time; the last row of
    a duplicated time wins"""
    columns, data = None, {}
    for start, cols, rows in sorted(segments, key=lambda s: s[0]):
        if columns is not None and cols != columns:
            fail("the force file columns change between restarts")
        columns = cols
        for t in [t for t in data if t >= start - 1e-12]:
            del data[t]
        for row in rows:
            data[row[0]] = row
    return columns, [data[t] for t in sorted(data)]


def force_rows(case: Path) -> list:
    files = force_files(case)
    if not files:
        return []
    segments = []
    for f in files:
        cols, rows = read_force_file(f)
        segments.append((float(f.parent.name), cols, rows))
    return merge_segments(segments)[1]


def hydrodynamic_forces(method: str, columns: dict, rows: list,
                        estimator: str = "surfaceTraction") -> dict:
    """Time and the hydrodynamic force estimates (x, y) of a method:
    - cutLink: the force columns are the selected estimator (inertia
      included); the momentum exchange or surface traction columns, if
      written, are an alternative estimate of the same force: they are
      reported separately and never added;
    - penalty, incremental: the forcing plus the inertia of the fluid
      inside the body (separate columns)"""
    t = [r[0] for r in rows]

    def vec(key):
        i = columns[key]
        return [(r[i], r[i + 1]) for r in rows]

    out = {"t": t}
    if method == "cutLink":
        out[estimator] = vec("force")
        for key in ("momentumExchange", "surfaceTraction"):
            if key in columns and key != estimator:
                out[key] = vec(key)
        if "inertia" in columns and any(
                abs(r[columns["inertia"]]) + abs(r[columns["inertia"] + 1])
                > 0 for r in rows):
            fail("cutLink inertia columns are not zero")
        out["total"] = out[estimator]
    else:
        forcing, inertia = vec("force"), vec("inertia")
        out["forcing"] = forcing
        out["inertia"] = inertia
        out["total"] = [(a[0] + b[0], a[1] + b[1])
                        for a, b in zip(forcing, inertia)]
    return out


# -----------------------------------------------------------------------------
# Time series metrics
# -----------------------------------------------------------------------------

def clean_series(t: list, y: list) -> tuple:
    """Sorted by time, duplicated times averaged"""
    data = {}
    for ti, yi in zip(t, y):
        data.setdefault(ti, []).append(yi)
    ts = sorted(data)
    return ts, [sum(data[x])/len(data[x]) for x in ts]


def interpolate(t: list, y: list, x: float) -> float:
    """Linear interpolation of a sorted series; nan outside its range"""
    if not t or x < t[0] - 1e-12 or x > t[-1] + 1e-12:
        return math.nan
    lo, hi = 0, len(t) - 1
    while hi - lo > 1:
        mid = (lo + hi)//2
        if t[mid] <= x:
            lo = mid
        else:
            hi = mid
    if t[hi] == t[lo]:
        return y[lo]
    s = (x - t[lo])/(t[hi] - t[lo])
    return y[lo] + s*(y[hi] - y[lo])


def trapezoid(t: list, f: list) -> float:
    return sum(0.5*(f[i] + f[i + 1])*(t[i + 1] - t[i])
               for i in range(len(t) - 1))


def rms_error(t_sim, y_sim, t_ref, y_ref, window=WINDOW, shift=0.0) -> dict:
    """Time-weighted RMS difference between a simulation and a reference
    over the window, normalised by the time-weighted RMS of the reference.
    Both are interpolated linearly onto the union of their sample times
    within the window (the simulation shifted by 'shift', for diagnostics
    only) and integrated with the trapezoidal rule; nothing is
    extrapolated"""
    import numpy as np
    ts, ys = clean_series(t_ref, y_ref)
    tr, yr = np.asarray(ts, float), np.asarray(ys, float)
    tsim, ysim = np.asarray(t_sim, float), np.asarray(y_sim, float)
    order = np.argsort(tsim, kind="stable")
    tsim, ysim = tsim[order], ysim[order]
    lo = max(window[0], tr[0], tsim[0] + shift)
    hi = min(window[1], tr[-1], tsim[-1] + shift)
    eps = 1e-12
    grid = np.union1d(np.union1d(tr[(tr >= lo - eps) & (tr <= hi + eps)],
                                 (tsim + shift)[(tsim + shift >= lo - eps)
                                                & (tsim + shift <= hi + eps)]),
                      [lo, hi])
    grid = grid[(grid >= lo - eps) & (grid <= hi + eps)]
    if hi - lo <= 0 or len(grid) < 3:
        return {"rms": math.nan, "ref_rms": math.nan, "nrms": math.nan,
                "n": int(len(grid)), "t0": lo, "t1": hi}
    sim = np.interp(grid - shift, tsim, ysim)
    ref = np.interp(grid, tr, yr)
    span = grid[-1] - grid[0]

    def integral(f):
        return float(np.sum(0.5*(f[1:] + f[:-1])*np.diff(grid)))
    rms = math.sqrt(integral((sim - ref)**2)/span)
    ref_rms = math.sqrt(integral(ref**2)/span)
    return {"rms": rms, "ref_rms": ref_rms,
            "nrms": rms/ref_rms if ref_rms > 0 else math.nan,
            "n": int(len(grid)), "t0": float(grid[0]), "t1": float(grid[-1])}


def best_lag(t_sim, y_sim, t_ref, y_ref, max_lag=0.1, step=0.0005) -> float:
    """Shift of the simulation (s, positive = simulation lags) minimising the
    RMS difference: a diagnostic only, never applied to the primary error"""
    best = (math.inf, 0.0)
    n = int(round(max_lag/step))
    for k in range(-n, n + 1):
        lag = k*step
        e = rms_error(t_sim, y_sim, t_ref, y_ref, shift=-lag)["rms"]
        if e < best[0]:
            best = (e, lag)
    return best[1]


def extrema(t, y, t0, t1) -> tuple:
    """(t_max, y_max, t_min, y_min) over [t0, t1], refined by a parabola
    through the extreme sample and its neighbours"""
    idx = [i for i, x in enumerate(t) if t0 <= x <= t1]
    if len(idx) < 3:
        return (math.nan,)*4
    out = []
    for pick in (max, min):
        i = pick(idx, key=lambda j: y[j])
        if 0 < i < len(t) - 1:
            y0, y1, y2 = y[i - 1], y[i], y[i + 1]
            den = y0 - 2*y1 + y2
            d = 0.5*(y0 - y2)/den if den != 0 else 0.0
            d = max(-1.0, min(1.0, d))
            h = 0.5*(t[i + 1] - t[i - 1])
            out += [t[i] + d*h, y1 - 0.25*(y0 - y2)*d]
        else:
            out += [t[i], y[i]]
    return tuple(out)


def high_frequency_rms(t, y, cutoff, window=WINDOW, guard=0.25) -> float:
    """RMS of the part of y above 'cutoff' (Hz): a zero-phase 4th-order
    Butterworth high-pass filter (scipy) of the series from the start of the
    window, evaluated over the window shortened by 'guard' (s) at each end
    to exclude the filter's edge transients. The series must be uniformly
    sampled; nan if scipy is unavailable"""
    try:
        import numpy as np
        from scipy.signal import butter, sosfiltfilt
    except ImportError:
        return math.nan
    tt = np.array(t)
    yy = np.array(y)
    mask = (tt >= window[0] - 1e-12) & (tt <= window[1] + 1e-12)
    tt, yy = tt[mask], yy[mask]
    if len(tt) < 64:
        return math.nan
    steps = np.diff(tt)
    dt = float(np.median(steps))
    if np.max(np.abs(steps - dt)) > 1e-6*dt + 1e-12:
        fail("high_frequency_rms needs a uniformly sampled series")
    if cutoff >= 0.5/dt:
        return math.nan
    sos = butter(4, cutoff, btype="highpass", fs=1.0/dt, output="sos")
    filtered = sosfiltfilt(sos, yy - np.polyval(np.polyfit(tt, yy, 1), tt))
    inner = (tt >= window[0] + guard) & (tt <= window[1] - guard)
    return float(np.sqrt(np.mean(filtered[inner]**2)))


def low_pass(t, y, cutoff, window=WINDOW) -> tuple:
    """The series from the start of the window, filtered by a zero-phase
    4th-order Butterworth low-pass filter (scipy): for the separately
    labelled low-frequency diagnostic only, never for the primary error.
    ([], []) if scipy is unavailable or the series is too short"""
    try:
        import numpy as np
        from scipy.signal import butter, sosfiltfilt
    except ImportError:
        return [], []
    tt = np.array(t)
    yy = np.array(y)
    mask = tt >= window[0] - 1e-12
    tt, yy = tt[mask], yy[mask]
    if len(tt) < 64:
        return [], []
    dt = float(np.median(np.diff(tt)))
    if np.max(np.abs(np.diff(tt) - dt)) > 1e-6*dt + 1e-12:
        fail("low_pass needs a uniformly sampled series")
    sos = butter(4, cutoff, btype="lowpass", fs=1.0/dt, output="sos")
    return list(tt), list(sosfiltfilt(sos, yy))


# -----------------------------------------------------------------------------
# References
# -----------------------------------------------------------------------------

def read_series(path: Path, column: int = 1) -> tuple:
    t, y = [], []
    for line in path.read_text().splitlines():
        if line.strip() and not line.startswith("#"):
            parts = line.replace(",", " ").split()
            try:
                t.append(float(parts[0]))
                y.append(float(parts[column]))
            except (ValueError, IndexError):
                continue
    return clean_series(t, y)


def references() -> dict:
    """The body-fitted references A and B (primary) and the digitised
    literature histories (secondary), as {name: {'Cd': (t, y), 'Cl': ...}}"""
    refs = {}
    pairs = {"A": ("CdBodyFitted.dat", "ClBodyFitted.dat"),
             "digitised": ("Cd.dat", "Cl.dat")}
    for name, (cd, cl) in pairs.items():
        if (REFERENCE_DIR / cd).exists():
            refs[name] = {"Cd": read_series(REFERENCE_DIR / cd),
                          "Cl": read_series(REFERENCE_DIR / cl)}
    b = REFERENCE_DIR / "bodyFitted901k.csv"
    if b.exists():
        refs["B"] = {"Cd": read_series(b, 1), "Cl": read_series(b, 2)}
    return refs


# -----------------------------------------------------------------------------
# Log diagnostics
# -----------------------------------------------------------------------------

def log_stats(case: Path) -> dict:
    """Maximum Courant number, continuity errors, time steps and execution
    time from the solver log"""
    log = case / "log.solids4Foam"
    out = {"max_co": math.nan, "max_continuity_local": math.nan,
           "final_continuity_cumulative": math.nan, "n_steps": 0,
           "execution_s": math.nan}
    if not log.exists():
        return out
    co, local, cumulative, steps, exe = [], [], [], 0, math.nan
    for line in log.read_text(errors="replace").splitlines():
        line = re.sub(r"^\[\d+,\d+\]<stdout>:", "", line)
        if line.startswith("Courant Number mean:"):
            co.append(float(line.split()[5]))
        elif line.startswith("time step continuity errors"):
            m = re.search(r"sum local = ([-+0-9.eE]+).*cumulative = "
                          r"([-+0-9.eE]+)", line)
            if m:
                local.append(abs(float(m.group(1))))
                cumulative.append(float(m.group(2)))
        elif line.startswith("Time = "):
            steps += 1
        elif line.startswith("ExecutionTime"):
            m = re.search(r"ExecutionTime = ([-+0-9.eE]+)", line)
            if m:
                exe = float(m.group(1))
    out.update({"max_co": max(co) if co else math.nan,
                "max_continuity_local": max(local) if local else math.nan,
                "final_continuity_cumulative": cumulative[-1] if cumulative
                else math.nan, "n_steps": steps, "execution_s": exe})
    return out


def wall_seconds(case: Path) -> float:
    path = case / "wallclock.txt"
    if not path.exists():
        return math.nan
    walls = re.findall(r"wall_s=(\d+)", path.read_text())
    return float(walls[-1]) if walls else math.nan


def mesh_hash(case: Path) -> str:
    import hashlib
    h = hashlib.sha256()
    mesh = case / "constant/polyMesh"
    for name in MESH_FILES:
        path = mesh / name
        if not path.exists():
            return ""
        text = re.sub(r"^\s*(note|location)\s+[^;]*;\s*$", "",
                      path.read_text(), flags=re.MULTILINE)
        h.update(text.encode())
    return h.hexdigest()[:16]


# -----------------------------------------------------------------------------
# Analysis
# -----------------------------------------------------------------------------

def analyse_case(case: Path, refs: dict) -> tuple:
    """Integrity checks and metrics of a case, and its coefficient
    histories"""
    meta = json.loads((case / "case.json").read_text())
    params = parameters(case)
    scale = params["coeff_scale"]
    problem = solver_problem(case)
    files = force_files(case)
    row = {"case": meta["case"], "run": meta["run"],
           "method": meta["method"], "level": meta["level"],
           "cells_per_D": meta["cells_per_D"], "deltaT": meta["deltaT"],
           "end_time": meta["end_time"], "ranks": meta["ranks"],
           "nOuterCorrectors": meta["nOuterCorrectors"],
           "nCorrectors": meta["nCorrectors"], "solver": problem or "End",
           "mesh_hash": mesh_hash(case)}
    row.update({k: meta.get("provenance", {}).get(k, "")
                for k in PROVENANCE_COLUMNS})
    histories = {}
    checks = []
    if not files:
        checks.append("no force file")
    else:
        segments = []
        for f in files:
            cols, rows = read_force_file(f)
            segments.append((float(f.parent.name), cols, rows))
        columns, rows = merge_segments(segments)
        forces = hydrodynamic_forces(meta["method"], columns, rows)
        t = forces["t"]
        expected = int(round(meta["end_time"]/meta["deltaT"]))
        row.update({"t_final": t[-1] if t else math.nan,
                    "n_force_rows": len(t), "n_expected": expected})
        if not t or abs(t[-1] - meta["end_time"]) > 0.5*meta["deltaT"]:
            checks.append("force history does not reach the end time")
        if len(t) < expected:
            checks.append(f"{len(t)} force rows, {expected} expected")
        finite = all(math.isfinite(v) for r in rows for v in r)
        if not finite:
            checks.append("non-finite force")
        for key, vecs in forces.items():
            if key == "t":
                continue
            histories[f"Cd_{key}"] = [scale*v[0] for v in vecs]
            histories[f"Cl_{key}"] = [scale*v[1] for v in vecs]
        histories["t"] = t
        for coeff in ("Cd", "Cl"):
            y = histories[f"{coeff}_total"]
            for ref_name, ref in refs.items():
                e = rms_error(t, y, *ref[coeff])
                row[f"{coeff}_nrms_{ref_name}"] = e["nrms"]
                row[f"{coeff}_rms_{ref_name}"] = e["rms"]
            if "B" in refs:
                row[f"{coeff}_lag_B"] = best_lag(t, y, *refs["B"][coeff])
                # Diagnostic only: the error of the force below 5 Hz, over
                # the window shortened by 0.25 s at each end (filter edges)
                tl, yl = low_pass(t, y, 5.0)
                if tl:
                    row[f"{coeff}_nrms_B_lowpass5Hz"] = rms_error(
                        tl, yl, *refs["B"][coeff],
                        window=(WINDOW[0] + 0.25, WINDOW[1] - 0.25))["nrms"]
                row[f"{coeff}_hf_rms_above_5Hz"] = high_frequency_rms(
                    t, y, 5.0)
            for k, (t0, t1) in enumerate(((WINDOW[0], 4.0), (4.0, 8.0))):
                tmax, ymax, tmin, ymin = extrema(t, y, t0, t1)
                row.update({f"{coeff}_max_{k + 1}": ymax,
                            f"{coeff}_tmax_{k + 1}": tmax,
                            f"{coeff}_min_{k + 1}": ymin,
                            f"{coeff}_tmin_{k + 1}": tmin})
        u_h = params["U_ref"]*meta["cells_per_D"]/params["D"]
        row["crossing_frequency_max_Hz"] = u_h
        if meta["method"] == "cutLink" and "Cd_momentumExchange" in histories:
            for coeff in ("Cd", "Cl"):
                a = histories[f"{coeff}_surfaceTraction"]
                b = histories[f"{coeff}_momentumExchange"]
                if "B" in refs:
                    e = rms_error(t, b, *refs["B"][coeff])
                    row[f"{coeff}_ME_nrms_B"] = e["nrms"]
                    if "A" in refs:
                        row[f"{coeff}_ME_nrms_A"] = rms_error(
                            t, b, *refs["A"][coeff])["nrms"]
                    row[f"{coeff}_ME_hf_rms_above_5Hz"] = high_frequency_rms(
                        t, b, 5.0)
                diff = [x - y for x, y in zip(a, b)]
                sel = [d for x, d in zip(t, diff)
                       if WINDOW[0] <= x <= WINDOW[1]]
                row[f"{coeff}_ST_minus_ME_rms"] = math.sqrt(
                    sum(d*d for d in sel)/len(sel)) if sel else math.nan
        if meta["method"] != "cutLink":
            sel = [c for x, c in zip(t, histories["Cd_inertia"])
                   if WINDOW[0] <= x <= WINDOW[1]]
            row["Cd_inertia_max_abs"] = max(abs(c) for c in sel) \
                if sel else math.nan
            if "B" in refs:
                row["Cd_forcing_only_nrms_B"] = rms_error(
                    t, histories["Cd_forcing"], *refs["B"]["Cd"])["nrms"]
    stats = log_stats(case)
    row.update(stats)
    wall = wall_seconds(case)
    row["wall_s"] = wall
    row["cpu_hours"] = wall*meta["ranks"]/3600.0
    if math.isfinite(stats["max_co"]) and stats["max_co"] > 1.0:
        checks.append(f"maximum Courant number {stats['max_co']:.2f} > 1")
    if problem:
        checks.append(problem)
    row["integrity"] = "PASS" if not checks else "FAIL"
    row["integrity_notes"] = "; ".join(checks)
    return row, histories


def cmd_analyse(args) -> int:
    refs = references()
    manifest = load_manifest()
    rows, all_hist = [], {}
    for name in sorted(manifest["cases"]):
        case = CASES_DIR / name
        if not (case / "case.json").exists() or not force_files(case):
            continue
        row, hist = analyse_case(case, refs)
        rows.append(row)
        all_hist[name] = hist
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    if rows:
        keys = []
        for r in rows:
            keys += [k for k in r if k not in keys]
        with (OUTPUT_DIR / "cases.csv").open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=keys,
                                    lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
        for name, hist in all_hist.items():
            keys = [k for k in hist if k != "t"]
            with (OUTPUT_DIR / f"history-{name}.csv").open(
                    "w", newline="") as handle:
                w = csv.writer(handle, lineterminator="\n")
                w.writerow(["t"] + keys)
                for i, t in enumerate(hist["t"]):
                    w.writerow([f"{t:.10g}"] + [f"{hist[k][i]:.8g}"
                                                for k in keys])
    hashes = {}
    for r in rows:
        if r["mesh_hash"]:
            hashes.setdefault(r["level"], set()).add(r["mesh_hash"])
    for level, h in sorted(hashes.items()):
        state = "identical" if len(h) == 1 else "DIFFERENT"
        print(f"L{level} background meshes: {state} ({', '.join(sorted(h))})")
    for r in rows:
        print(f"{r['case']:34s} {r['integrity']:4s} "
              f"Cd nRMS A {r.get('Cd_nrms_A', math.nan):.4f} "
              f"B {r.get('Cd_nrms_B', math.nan):.4f} "
              f"Cl nRMS B {r.get('Cl_nrms_B', math.nan):.4f} "
              f"maxCo {r['max_co']:.3f} CPU-h {r['cpu_hours']:.3f} "
              f"{r['integrity_notes']}")
    failed = [r for r in rows if r["integrity"] != "PASS"]
    complete = [r for r in rows if r["integrity"] == "PASS"]
    print(f"Numerical integrity: {len(complete)} PASS, {len(failed)} FAIL"
          f"{'' if rows else ' (no cases analysed: SKIPPED)'}")
    print(f"Results: {OUTPUT_DIR / 'cases.csv'}")
    return 1 if failed else 0


# Reproduction check of an analysed case against the tracked result of the
# same case (reference/results/cases.csv): each metric must agree within
# CHECK_RTOL relative or CHECK_ATOL absolute (fractions). The tolerances
# allow for a different compiler, MPI library or decomposition (serial and
# parallel runs of the same case agree to about 1e-5 relative) but detect a
# change of the method, its settings or the force estimator; they are not
# accuracy targets
CHECK_METRICS = ("Cd_nrms_B", "Cl_nrms_B", "Cd_nrms_A", "Cl_nrms_A",
                 "Cd_hf_rms_above_5Hz")
# Metrics that need scipy (a Butterworth filter); without it they are nan and
# are skipped, never compared. The raw-force metrics do not need scipy
SCIPY_METRICS = ("Cd_hf_rms_above_5Hz",)
CHECK_RTOL = 0.05
CHECK_ATOL = 0.002


def have_scipy() -> bool:
    try:
        import scipy.signal  # noqa: F401
    except ImportError:
        return False
    return True


def cmd_check(args) -> int:
    """PASS/FAIL of the analysed cases (postProcessing/cases.csv) against
    the tracked results: numerical integrity and the reproduction of the
    normalised RMS differences and the fluctuation level; cases without a
    tracked result (e.g. another quadrature) are reported, not checked.
    The metrics that need scipy are skipped, and reported as skipped, only
    if scipy is not installed; a nonfinite value of any other metric (or of
    one of them with scipy installed) is a failure"""
    current = OUTPUT_DIR / "cases.csv"
    tracked = REFERENCE_DIR / "results" / "cases.csv"
    if not current.exists():
        fail("no postProcessing/cases.csv: run ./Allverify analyse first")
    with current.open() as handle:
        rows = list(csv.DictReader(handle))
    with tracked.open() as handle:
        reference = {r["case"]: r for r in csv.DictReader(handle)}
    failures, checked, skipped = 0, 0, 0
    recorded = {(r.get("openfoam_version", ""), r.get("solver_build", ""))
                for r in rows if r.get("openfoam_version")}
    for version, build in sorted(recorded):
        print(f"Solver of the analysed cases: OpenFOAM {version}, "
              f"solids4foam build {build or 'not recorded'}; tracked "
              f"results: OpenFOAM {BASELINE['openfoam']}, solids4foam "
              f"{BASELINE['solids4foam']}")
    scipy_installed = have_scipy()
    for r in rows:
        if r["integrity"] != "PASS":
            print(f"FAIL {r['case']}: numerical integrity "
                  f"({r['integrity_notes']})")
            failures += 1
            continue
        ref = reference.get(r["case"])
        if ref is None:
            print(f"---- {r['case']}: no tracked result, not checked")
            continue
        for key in CHECK_METRICS:
            if not ref.get(key) or not r.get(key):
                continue
            value, expected = float(r[key]), float(ref[key])
            if (key in SCIPY_METRICS and not scipy_installed
                    and not math.isfinite(value)):
                skipped += 1
                print(f"SKIP {r['case']}: {key} needs scipy, which is not "
                      "installed (tracked "
                      f"{expected:.5g}); not compared")
                continue
            ok = abs(value - expected) <= max(CHECK_RTOL*abs(expected),
                                              CHECK_ATOL)
            checked += 1
            if not ok:
                failures += 1
            print(f"{'PASS' if ok else 'FAIL'} {r['case']}: {key} "
                  f"{value:.5g} (tracked {expected:.5g})")
    print(f"Check: {checked} comparisons, {failures} failure(s)"
          f"{f', {skipped} skipped (scipy)' if skipped else ''}"
          f"{'' if checked else ' (nothing compared: SKIPPED)'}")
    return 1 if failures else 0


EXPORT_COLUMNS = (
    "case", "run", "method", "level", "cells_per_D", "deltaT", "end_time",
    "ranks", "nOuterCorrectors", "nCorrectors", "integrity", "t_final",
    "n_steps", "max_co", "max_continuity_local", "wall_s", "cpu_hours",
    "Cd_nrms_A", "Cd_nrms_B", "Cd_nrms_digitised", "Cl_nrms_A", "Cl_nrms_B",
    "Cl_nrms_digitised", "Cd_nrms_B_lowpass5Hz", "Cl_nrms_B_lowpass5Hz",
    "Cd_hf_rms_above_5Hz", "Cl_hf_rms_above_5Hz", "Cd_lag_B", "Cl_lag_B",
    "Cd_max_1", "Cd_tmax_1", "Cd_min_1", "Cd_tmin_1", "Cd_max_2",
    "Cd_tmax_2", "Cd_min_2", "Cd_tmin_2", "Cd_ME_nrms_A", "Cd_ME_nrms_B",
    "Cl_ME_nrms_A", "Cl_ME_nrms_B", "Cd_ME_hf_rms_above_5Hz",
    "Cd_ST_minus_ME_rms", "Cl_ST_minus_ME_rms", "Cd_inertia_max_abs",
    "Cd_forcing_only_nrms_B", "mesh_hash")
EXPORT_STEP = 0.02
# Run whose history is not exported: the serial repeat of R2A (identical
# metrics, kept in cases.csv)
HISTORY_NOT_TRACKED = ("O4",)
ZOOM = (1.8, 2.2)


def cmd_export(args) -> int:
    """Compact results of the analysed complete cases for reference/results/:
    cases.csv (selected columns, computed from the full histories, also of
    the shorter runs) and the 0.02 s histories of the 8 s runs that the
    figures draw (not the serial repeat O4, whose metrics are in cases.csv,
    nor the shorter runs). With --full-resolution, also a window of 1.8-2.2 s
    of the finest 8 s case of each method, for cell_crossing only; these
    files are not tracked. The figures other than cell_crossing and the check
    can be regenerated from the tracked files without work/"""
    src = OUTPUT_DIR / "cases.csv"
    if not src.exists():
        fail("no postProcessing/cases.csv: run ./Allverify analyse first")
    with src.open() as handle:
        rows = [r for r in csv.DictReader(handle)
                if r["integrity"] == "PASS"]
    out = REFERENCE_DIR / "results"
    out.mkdir(parents=True, exist_ok=True)
    with (out / "cases.csv").open("w", newline="") as handle:
        w = csv.writer(handle, lineterminator="\n")
        # Provenance columns only where the analysed cases recorded them
        columns = list(EXPORT_COLUMNS) + [
            k for k in PROVENANCE_COLUMNS if any(r.get(k) for r in rows)]
        w.writerow(columns)
        for r in rows:
            w.writerow([r.get(k, "") for k in columns])
    full = [r for r in rows if abs(float(r["end_time"]) - END_TIME) < 1e-9]
    # The full-resolution window of the finest case of each method only
    finest = {}
    for r in full:
        key = (int(r["level"]), -float(r["deltaT"]))
        if r["method"] not in finest or key > finest[r["method"]][0]:
            finest[r["method"]] = (key, r["case"])
    zoom_cases = {c for _, c in finest.values()} \
        if getattr(args, "full_resolution", False) else set()
    for old in list(out.glob("zoom-*.csv")) + list(out.glob("history-*.csv")):
        old.unlink()
    for r in full:
        if r["run"] in HISTORY_NOT_TRACKED:
            continue
        with (OUTPUT_DIR / f"history-{r['case']}.csv").open() as handle:
            hist = list(csv.reader(handle))
        keys, data = hist[0], hist[1:]
        dec = [d for d in data if abs(float(d[0])/EXPORT_STEP
                                      - round(float(d[0])/EXPORT_STEP))
               < 1e-6]
        zoom = [d for d in data if ZOOM[0] <= float(d[0]) <= ZOOM[1]]
        for name, sel in (("history", dec), ("zoom", zoom)):
            if name == "zoom" and r["case"] not in zoom_cases:
                continue
            with (out / f"{name}-{r['case']}.csv").open("w",
                                                        newline="") as h:
                w = csv.writer(h, lineterminator="\n")
                w.writerow(keys)
                w.writerows(sel)
    print(f"Exported {len(rows)} cases to {out}")
    return 0


def cmd_plot(args) -> int:
    script = Path(__file__).resolve().parent / "plot_verification.py"
    return subprocess.run([sys.executable, str(script)] + args.extra,
                          check=False).returncode


# -----------------------------------------------------------------------------
# Command line
# -----------------------------------------------------------------------------

def main(argv=None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    commands = ("run", "generate", "dry-run", "submit", "status", "analyse",
                "check", "export", "plot")
    command = argv.pop(0) if argv and argv[0] in commands else "run"
    parser = argparse.ArgumentParser(
        prog="Allverify" + ("" if command == "run" else f" {command}"),
        description=__doc__.split("\n\n")[0])
    if command in ("run", "generate", "dry-run", "submit"):
        parser.add_argument("--runs", nargs="+", default=list(DEFAULT_RUNS),
                            help="runs of the campaign: "
                            + ", ".join(RUNS))
        parser.add_argument("--ranks", type=int, default=0,
                            help="MPI ranks per case (default by level: "
                            "1, 2, 8, 8)")
        parser.add_argument("--force", action="store_true",
                            help="regenerate or resubmit existing cases")
        parser.add_argument("--traction-quadrature",
                            default=HISTORICAL_QUADRATURE,
                            choices=QUADRATURES,
                            help="cutLink surface-traction quadrature, "
                            "written in fvOptions: greedy (default) is that "
                            "of the tracked results; another choice tags "
                            "the case names and needs a solver with the "
                            "tractionQuadrature keyword")
    if command == "run":
        parser.add_argument("--dry-run", action="store_true",
                            help="list the cases; run nothing")
    if command == "submit":
        parser.add_argument("--env-script",
                            default=os.environ.get("S4F_ENV_SCRIPT", ""),
                            help="script sourced in each job (default "
                            "$S4F_ENV_SCRIPT)")
        parser.add_argument("--partition", default="")
        parser.add_argument("--max-cpus", type=int, default=0,
                            help="refuse to exceed this many queued CPUs")
        parser.add_argument("--time", default="")
        parser.add_argument("--mem", default="")
        parser.add_argument("--sbatch-args", default="",
                            help="further sbatch options")
        parser.add_argument("--yes", action="store_true",
                            help="submit (otherwise only a summary)")
    if command == "export":
        parser.add_argument("--full-resolution", action="store_true",
                            help="also write the 1.8-2.2 s window of the "
                            "finest case of each method (zoom-*.csv, for "
                            "cell_crossing)")
    if command == "plot":
        parser.add_argument("extra", nargs=argparse.REMAINDER)
    args = parser.parse_args(argv)
    if hasattr(args, "ranks") and args.ranks < 0:
        parser.error("--ranks must be positive")
    return {"run": cmd_run, "generate": cmd_generate,
            "dry-run": cmd_dry_run, "submit": cmd_submit,
            "status": cmd_status, "analyse": cmd_analyse,
            "check": cmd_check, "export": cmd_export,
            "plot": cmd_plot}[command](args)


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (RuntimeError, subprocess.SubprocessError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        sys.exit(2)
