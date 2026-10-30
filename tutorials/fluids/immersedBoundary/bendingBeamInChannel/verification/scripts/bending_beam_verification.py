#!/usr/bin/env python3
"""Verification of the bendingBeamInChannel tutorial: the beam bending
through a stationary mesh with the cutLink, penalty, legacy incremental and
(exploratory) ghostCell methods of immersedBoundaryForce.

Three things are checked, and kept apart:

  A. the exact volume balance. The prescribed bending changes the volume of
     the beam inside the channel at the analytic rate
     dV/dt = c 2 pi sin(4 pi t), so the net flux through the inlet and the
     outlet of an incompressible flow must equal it. It is a property of the
     formulation and needs no reference solution;
  B. the drag and lift coefficients against the tracked body-fitted series
     (reference/CdBodyFitted.dat), as a normalised RMS difference. The
     reference has unresolved provenance, so this is an agreement and not
     an error, and the lift is not used to rank methods;
  C. the observed cost (allocated rank-hours, steps, time per cell-step),
     reported with the conditions of the run and never as an intrinsic cost.

Each case is a copy of the tutorial under verification/work/cases/ with
only the method settings, the mesh level, the end time, the outer
correctors, the write interval and the decomposition changed. The tutorial
itself is never modified. Everything but a few Allrun-type commands is plain
Python 3 (standard library only); numpy, matplotlib and SLURM are not
required. A manifest, work/manifest.json, records every case with the
SHA-256 of its inputs.

Usage:
    ./Allverify [run] [--runs ...]   generate, run locally, analyse, check
    ./Allverify unit                 exact tests that need no OpenFOAM
    ./Allverify dry-run --runs ...   list the cases, ranks and estimated cost
    ./Allverify generate --runs ...  generate the cases only
    ./Allverify status               state of the cases in the manifest
    ./Allverify analyse              forces, flux and cost of completed cases
    ./Allverify check                PASS / FAIL / SKIPPED against the tracked
                                     results (reference/results/)

The default run, 'smoke', is a cutLink and a penalty case on Level 1 for
the first 0.5 s (the volume-rate period): a few minutes on one core.
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
REFERENCE_FILE = REFERENCE_DIR / "CdBodyFitted.dat"

METHODS = ("cutLink", "penalty", "incremental", "ghostCell")
LEVELS = (1, 2, 3, 4)
CELLS = {1: 20020, 2: 80080, 3: 320320, 4: 1281280}
CELLS_PER_LEVEL = 20020
N_CORRECTORS = 3
N_OUTER = {"cutLink": 1, "penalty": 1, "incremental": 5, "ghostCell": 1}
WRITE_INTERVAL = 0.125          # a quarter of the volume-rate period 0.5 s
WINDOW = (1.0, 8.0)             # comparison window of the force series
PERIOD = 1.0                    # bending period

# Settings written explicitly into every case. They are those of the
# stationary and oscillating cylinder verifications; ghostCell is
# documented for static bodies and is exploratory here.
METHOD_SETTINGS = {
    "cutLink": ([("forceEstimator", "surfaceTraction")], []),
    "penalty": ([("method", "penalty"), ("penaltyCoeff", "1e3"),
                 ("weighting", "volumeFraction"),
                 ("surfaceRateCoeff", "3")],
                [("occupancy", "signedDistance")]),
    "incremental": ([("method", "incremental"), ("couplingCoeff", "0.2")],
                    [("occupancy", "vertexFraction")]),
    "ghostCell": ([("method", "ghostCell"), ("penaltyCoeff", "1e3"),
                   ("imageDistance", "1.5"),
                   ("imageInterpolation", "cellPoint")], []),
}

# Observed allocated rank-hours of the tracked runs (eight seconds, adaptive
# time step, single runs on a shared node with SLURM or direct MPI): context
# for the estimate printed by dry-run, not a performance model
OBSERVED_RANK_HOURS = {("cutLink", 1): 0.057, ("cutLink", 2): 0.505,
                       ("cutLink", 3): 5.12, ("penalty", 1): 0.042,
                       ("penalty", 2): 0.282, ("penalty", 3): 2.60,
                       ("incremental", 1): 0.140, ("incremental", 2): 1.256,
                       ("ghostCell", 1): 0.046,
                       # measured for 0.5 s (1.709 allocated rank-hours on 4
                       # ranks), scaled to the eight-second equivalent
                       ("cutLink", 4): 27.3}
# Mean number of adaptive steps of an eight-second run (tracked runs)
STEPS_8S = {1: 1280, 2: 2559, 3: 5216, 4: 10400}
LEVEL_RANKS = {1: 1, 2: 2, 3: 4, 4: 4}

# name: (description, [(method, level, study, end_time)])
STUDY_A, STUDY_B = "A", "B"
RUNS = {
    "smoke": ("cutLink and penalty, Level 1, 0 to 0.5 s",
              [("cutLink", 1, STUDY_A, 0.5), ("penalty", 1, STUDY_A, 0.5)]),
    "flux": ("cutLink and penalty, Levels 1-2, 0 to 0.5 s",
             [(m, l, STUDY_A, 0.5) for m in ("cutLink", "penalty")
              for l in (1, 2)]),
    "L1": ("all four methods, Level 1, 0 to 8 s",
           [(m, 1, STUDY_A, 8.0) for m in METHODS]),
    "L2": ("cutLink, penalty and incremental, Level 2, 0 to 8 s",
           [(m, 2, STUDY_A, 8.0) for m in METHODS[:3]]),
    "L3": ("cutLink and penalty, Level 3, 0 to 8 s (expensive)",
           [(m, 3, STUDY_A, 8.0) for m in ("cutLink", "penalty")]),
    "studyB": ("cutLink, A_y = 0 (exact flux zero), Levels 1-3, 0 to 0.5 s",
               [("cutLink", l, STUDY_B, 0.5) for l in (1, 2, 3)]),
    "L4flux": ("cutLink, Level 4, 0 to 0.5 s (expensive: 1.28 M cells)",
               [("cutLink", 4, STUDY_A, 0.5)]),
}
DEFAULT_RUNS = ("smoke",)

# Tolerances of the check against the tracked results. Both are the thresholds
# registered in advance for the verification campaign (a change of 1.0 point
# of the Cd nRMS, or of 0.1 in the second harmonic, is material). They were
# not adjusted after looking at any result. The tracked values come from runs
# with the adaptive time step of the tutorial, which is what the generated
# cases use; the archived fixed-time-step penalty runs of the campaign differ
# from the tracked second harmonic by 0.06 to 0.18 and are not conforming
# reruns (see the README)
CHECK_NRMS = 1.0
CHECK_A2 = 0.1
CHECK_FLUX_RATIO = 0.01
ZERO_FLUX_FRACTION = 1e-6       # of the analytic peak
ANALYTIC_RTOL = 1e-4

FOAM_HEADER = """\
/*--------------------------------*- C++ -*----------------------------------*\\
| solids4foam: solid mechanics and fluid-solid interaction simulations        |
| Version:     v2.4                                                           |
| Web:         https://solids4foam.github.io                                  |
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


def fail(message: str) -> None:
    sys.exit(f"ERROR: {message}")


def fmt(value: float) -> str:
    return f"{value:.10g}"


# ---------------------------------------------------------------------------
# Tutorial parameters and the exact volume rate
# ---------------------------------------------------------------------------

def dict_scalar(text: str, key: str) -> float:
    match = re.search(rf"^\s*{re.escape(key)}\s+([-+0-9.eE]+)\s*;", text,
                      re.MULTILINE)
    if not match:
        fail(f"entry '{key}' not found")
    return float(match.group(1))


def parameters(tutorial: Path = TUTORIAL_DIR) -> dict:
    """Geometry, motion and normalisation read from the tutorial"""
    fv = (tutorial / "constant/fvOptions").read_text()
    block = (tutorial / "system/blockMeshDict").read_text()
    inlet = (tutorial / "0/U").read_text()
    z = [float(v) for v in re.findall(
        r"\(\s*[-+0-9.eE]+\s+[-+0-9.eE]+\s+([-+0-9.eE]+)\s*\)",
        block.split("vertices", 1)[1].split("blocks", 1)[0])]
    y_top = max(float(v) for v in re.findall(
        r"\(\s*[-+0-9.eE]+\s+([-+0-9.eE]+)\s+[-+0-9.eE]+\s*\)",
        block.split("vertices", 1)[1].split("blocks", 1)[0]))
    # parabolic inlet Ux = a1 y + a2 y^2 on 0 < y < y_top: mean velocity
    poly = re.findall(r"\(\(([-+0-9.eE]+)\s+0\s+0\)\s+\(([0-9]+)\s+0\s+0\)\)",
                      inlet)
    u_mean = sum(float(c) * y_top ** int(p) / (int(p) + 1)
                 for c, p in poly)
    params = {
        "Ax": dict_scalar(fv, "amplitudeX"),
        "Ay": dict_scalar(fv, "amplitudeY"),
        "period": dict_scalar(fv, "period"),
        "H": dict_scalar(fv, "yMax"),
        "xCenter": dict_scalar(fv, "xCenter"),
        "width": 0.2,
        "Lz": max(z) - min(z),
        "U_ref": round(u_mean, 8),
    }
    params["coeff_scale"] = 2.0 / (params["U_ref"] ** 2 * params["H"]
                                   * params["Lz"])
    return params


def profile_g(s: float) -> float:
    return 10 * s ** 3 - 15 * s ** 4 + 6 * s ** 5


def profile_gp(s: float) -> float:
    return 30 * s ** 2 - 60 * s ** 3 + 30 * s ** 4


def volume_rate_amplitude(params: dict, ay: float = None) -> float:
    """Peak of dV/dt (m^3/s): c omega, with c = (10/7) w Lz A_x A_y / H and
    omega = 2 pi/period. The volume of the beam in the channel changes as
    V(t) = V0 + (c/2)(1 - cos 2 omega t), i.e. at twice the bending frequency"""
    ay = params["Ay"] if ay is None else ay
    c = (10.0 / 7.0) * params["width"] * params["Lz"] * params["Ax"] * ay \
        / params["H"]
    return c * 2.0 * math.pi / params["period"]


def volume_rate(params: dict, t: float, ay: float = None) -> float:
    """Exact dV/dt of the beam in the channel, m^3/s, positive when the
    beam volume grows: the net outward flux of the fluid through the inlet
    and the outlet of an incompressible flow"""
    omega = 2.0 * math.pi / params["period"]
    return volume_rate_amplitude(params, ay) * math.sin(2.0 * omega * t)


def beam_area(params: dict, t: float, ay: float = None, n: int = 4001) -> float:
    """Area of the deformed beam outline (shoelace), m^2, from the motion
    of customProfileBend: dx = Ax g(xi) sin(omega t),
    dy = -(x0 - xc) Ay/H g'(xi) sin(omega t), xi = y0/H"""
    ay = params["Ay"] if ay is None else ay
    sn = math.sin(2.0 * math.pi * t / params["period"])
    left, right = [], []
    for i in range(n):
        s = i / (n - 1)
        for side, x0 in ((left, -params["width"] / 2),
                         (right, params["width"] / 2)):
            side.append((x0 + params["Ax"] * profile_g(s) * sn,
                         s * params["H"] - x0 * ay / params["H"]
                         * profile_gp(s) * sn))
    outline = left + right[::-1]
    return 0.5 * abs(sum(x1 * y2 - x2 * y1 for (x1, y1), (x2, y2)
                         in zip(outline, outline[1:] + outline[:1])))


def polygon_volume_rate(params: dict, t: float, ay: float = None) -> float:
    """dV/dt from the shoelace area by a central difference: an independent
    evaluation of volume_rate"""
    h = 1e-4
    return params["Lz"] * (beam_area(params, t + h, ay)
                           - beam_area(params, t - h, ay)) / (2 * h)


# ---------------------------------------------------------------------------
# Plan, ranks and cost estimate
# ---------------------------------------------------------------------------

def load_manifest() -> dict:
    if MANIFEST.exists():
        return json.loads(MANIFEST.read_text())
    return {"cases": {}}


def save_manifest(manifest: dict) -> None:
    WORK_DIR.mkdir(parents=True, exist_ok=True)
    MANIFEST.write_text(json.dumps(manifest, indent=1, sort_keys=True) + "\n")


def case_name(method: str, level: int, study: str, end: float) -> str:
    return f"{method}-L{level}-{study}-t{fmt(end).replace('.', 'p')}"


def plan(args) -> list:
    """[(method, level, study, end_time, name, ranks)] of the selected runs"""
    cases = []
    for run in args.runs:
        if run not in RUNS:
            fail(f"unknown run '{run}'; choose from {', '.join(RUNS)}")
        for method, level, study, end in RUNS[run][1]:
            ranks = args.ranks or LEVEL_RANKS[level]
            entry = (method, level, study, end,
                     case_name(method, level, study, end), ranks)
            if entry not in cases:
                cases.append(entry)
    return cases


def estimated_rank_hours(method: str, level: int, end: float):
    """Observed allocated rank-hours of the tracked eight-second run, scaled
    with the end time; None where no tracked run of the method exists"""
    observed = OBSERVED_RANK_HOURS.get((method, level))
    return None if observed is None else observed * end / 8.0


def replace_entry(path: Path, key: str, value: str) -> None:
    text = path.read_text()
    pattern = re.compile(rf"^(\s*{re.escape(key)}\s+)[^;\n]*;", re.MULTILINE)
    new, count = pattern.subn(lambda m: f"{m.group(1)}{value};", text)
    if count != 1:
        fail(f"found {count} '{key}' entries in {path}")
    path.write_text(new)


def write_fv_options(path: Path, method: str, study: str) -> None:
    """The tutorial's fvOptions with the method settings of METHOD_SETTINGS
    (and A_y = 0 for Study B); the motion block is never rewritten"""
    text = path.read_text()
    option, body = METHOD_SETTINGS[method]
    width = max(len(k) for k, _ in option) + 4
    lines = "".join(f"    {k.ljust(width)}{v};\n" for k, v in option)
    pattern = re.compile(r"    // cutLink \(the default\).*?"
                         r"forceEstimator\s+surfaceTraction;\n", re.DOTALL)
    new, count = pattern.subn(lambda m: f"    // {method}: generated by "
                              "verification/scripts/bending_beam_"
                              "verification.py\n" + lines, text)
    if count != 1:
        fail("the forceEstimator entry of the tutorial fvOptions was not "
             "found exactly once")
    if body:
        anchor = re.search(r"(\n\s*CofR\s+\([^)]*\);\n)", new)
        if not anchor:
            fail("CofR entry of the tutorial fvOptions not found")
        extra = "".join(f"            {k.ljust(12)}{v};\n" for k, v in body)
        new = new[:anchor.end()] + extra + new[anchor.end():]
    path.write_text(new)
    if study == STUDY_B:
        replace_entry(path, "amplitudeY", "0")


def write_decompose_par_dict(path: Path, ranks: int) -> None:
    text = FOAM_HEADER.format(location="system", object="decomposeParDict")
    text += (f"numberOfSubdomains  {ranks};\n\nmethod          simple;\n\n"
             f"simpleCoeffs\n{{\n    n           ({ranks} 1 1);\n}}\n\n"
             "// **************************************************************"
             "*********** //\n")
    path.write_text(text)


RUN_SCRIPT = """\
#!/bin/bash
# Generated by verification/scripts/bending_beam_verification.py
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
    # the fluxes through the inlet and the outlet are read from the
    # reconstructed fields
    solids4Foam::runApplication reconstructPar -newTimes -fields '(phi)'
)"""


def input_hashes(case: Path) -> dict:
    files = ["constant/fvOptions", "system/controlDict", "system/fvSolution",
             "system/fvSchemes", "system/blockMeshDict",
             "system/decomposeParDict", "0/U", "0/p"]
    return {f: hashlib.sha256((case / f).read_bytes()).hexdigest()
            for f in files if (case / f).exists()}


def generate_case(method: str, level: int, study: str, end: float,
                  name: str, ranks: int, force: bool = False,
                  cases_dir: Path = None) -> Path:
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
    write_fv_options(case / "constant/fvOptions", method, study)
    control = case / "system/controlDict"
    replace_entry(control, "endTime", fmt(end))
    replace_entry(control, "writeInterval", fmt(WRITE_INTERVAL))
    fv_solution = case / "system/fvSolution"
    replace_entry(fv_solution, "nOuterCorrectors", str(N_OUTER[method]))
    replace_entry(fv_solution, "nCorrectors", str(N_CORRECTORS))
    write_decompose_par_dict(case / "system/decomposeParDict", ranks)
    script = case / "Allrun.verify"
    script.write_text(RUN_SCRIPT.format(
        level=level, ranks=ranks,
        body=PARALLEL_BODY if ranks > 1 else SERIAL_BODY))
    script.chmod(0o755)
    return case


def cmd_generate(args) -> int:
    manifest = load_manifest()
    for method, level, study, end, name, ranks in plan(args):
        case = generate_case(method, level, study, end, name, ranks,
                             args.force)
        manifest["cases"][name] = {
            "method": method, "level": level, "study": study,
            "end_time": end, "ranks": ranks,
            "inputs_sha256": input_hashes(case),
            "generated": datetime.now().isoformat(timespec="seconds")}
        print(f"{name}: generated in {case}")
    save_manifest(manifest)
    return 0


def cmd_dry_run(args) -> int:
    params = parameters()
    print(f"H = {params['H']:.6g} m, A_x = {params['Ax']:.3g} m, "
          f"A_y = {params['Ay']:.3g} m, Lz = {params['Lz']:.3g} m, "
          f"U_mean = {params['U_ref']:.6g} m/s, coefficient scale = "
          f"{params['coeff_scale']:.6g}")
    print(f"exact volume-rate amplitude = "
          f"{volume_rate_amplitude(params):.6g} m^3/s (Study A), "
          f"0 (Study B)")
    total = 0.0
    for method, level, study, end, name, ranks in plan(args):
        est = estimated_rank_hours(method, level, end)
        total += est or 0.0
        print(f"  {name}: {CELLS[level]:,} cells, {ranks} rank(s), "
              f"nOuter {N_OUTER[method]}, "
              + (f"observed-cost estimate {est:.3f} allocated rank-hours"
                 if est is not None else "no observed cost to scale"))
    print(f"estimated total {total:.2f} allocated rank-hours (scaled from "
          "single runs on a shared node; not a performance model)")
    print("Dry run: nothing generated, run or submitted")
    return 0


def check_environment(ranks: int) -> str:
    """Return '' if the environment can run the case, else the reason"""
    for executable in ("blockMesh", "solids4Foam", "solids4FoamScripts.sh"):
        if shutil.which(executable) is None:
            return f"'{executable}' is unavailable: source OpenFOAM and " \
                   "solids4foam first"
    if ranks > 1:
        for executable in ("decomposePar", "reconstructPar", "mpirun"):
            if shutil.which(executable) is None:
                return f"'{executable}' is unavailable for {ranks} ranks"
    try:
        cores = len(os.sched_getaffinity(0))
    except AttributeError:
        cores = os.cpu_count() or 1
    if ranks > cores:
        return (f"{ranks} ranks requested but the process may use {cores} "
                "core(s): MPI ranks would share cores (check the CPU "
                "affinity of this shell)")
    return ""


def solver_problem(case: Path) -> str:
    """Empty if the solver log shows a normal end"""
    log = case / "log.solids4Foam"
    if not log.exists():
        return "no solver log"
    text = log.read_text(errors="ignore")
    if not re.search(r"^End\s*$", text, re.MULTILINE):
        return "the solver did not reach 'End'"
    # the start-up message "sigFpe : Enabling floating point exception
    # trapping" is not a failure
    text = re.sub(r"^.*trapping.*$", "", text, flags=re.MULTILINE)
    if re.search(r"\bnan\b|sigFpe caught|FOAM FATAL|Floating point exception",
                 text, re.IGNORECASE):
        return "nan or a floating point exception in the log"
    return ""


def cmd_run(args) -> int:
    if args.dry_run:
        return cmd_dry_run(args)
    manifest = load_manifest()
    for method, level, study, end, name, ranks in plan(args):
        case = CASES_DIR / name
        if case.exists() and not solver_problem(case) and not args.force:
            print(f"{name}: completed run reused")
            continue
        reason = check_environment(ranks)
        if reason:
            fail(f"{name}: {reason}")
        case = generate_case(method, level, study, end, name, ranks, True)
        load = os.getloadavg()[0] if hasattr(os, "getloadavg") else -1.0
        manifest["cases"][name] = {
            "method": method, "level": level, "study": study,
            "end_time": end, "ranks": ranks, "ran_locally": True,
            "inputs_sha256": input_hashes(case),
            "load_average_at_start": load,
            "cores": os.cpu_count(),
            "started": datetime.now().isoformat(timespec="seconds")}
        save_manifest(manifest)
        print(f"{name}: running on {ranks} rank(s)", flush=True)
        with (case / "log.Allverify").open("w") as handle:
            subprocess.run(["./Allrun.verify"], cwd=case, stdout=handle,
                           stderr=subprocess.STDOUT, check=False)
        manifest["cases"][name]["finished"] = \
            datetime.now().isoformat(timespec="seconds")
        save_manifest(manifest)
        problem = solver_problem(case)
        if problem:
            fail(f"{name}: {problem}; see {case / 'log.solids4Foam'}")
    status = cmd_analyse(args)
    return cmd_check(args) or status


def cmd_status(args) -> int:
    manifest = load_manifest()
    if not manifest["cases"]:
        print("no cases in the manifest")
    for name, entry in sorted(manifest["cases"].items()):
        case = CASES_DIR / name
        state = "generated"
        if (case / "log.solids4Foam").exists():
            state = solver_problem(case) or "completed"
        print(f"{name}: {state}")
    return 0


# ---------------------------------------------------------------------------
# Forces
# ---------------------------------------------------------------------------

def read_force_file(path: Path) -> tuple:
    """(column map, rows) of postProcessing/immersedBoundary/0/beam.dat.
    The header names the groups of three columns, separated by tabs"""
    header, rows = None, []
    for line in path.read_text().splitlines():
        if line.startswith("# Time"):
            header = line
        elif line.strip() and not line.startswith("#"):
            rows.append([float(v) for v in line.split()])
    if header is None or not rows:
        fail(f"no force data in {path}")
    groups = [g.strip() for g in header[1:].split("\t")][1:]
    columns = {}
    for i, group in enumerate(groups):
        name = group.split("(")[0].strip()
        columns[name] = 1 + 3 * i
    return columns, rows


def hydrodynamic_forces(columns: dict, rows: list, scale: float) -> dict:
    """Coefficient series per estimator, never mixed:

    cutLink writes the surface traction ('force') and the momentum exchange
    ('force from the momentum exchange'): both are the hydrodynamic force.
    The other methods write the forcing ('force') and the inertia of the fluid
    inside the body ('inertia of the fluid inside'), whose sum is the
    hydrodynamic force; the forcing alone is not
    """
    t = [r[0] for r in rows]
    force = columns["force"]
    me = columns.get("force from the momentum exchange")
    inertia = columns.get("inertia of the fluid inside")
    out = {"t": t}
    if me is not None:
        out["cutLink ST"] = ([scale * r[force] for r in rows],
                             [scale * r[force + 1] for r in rows])
        out["cutLink ME"] = ([scale * r[me] for r in rows],
                             [scale * r[me + 1] for r in rows])
    elif inertia is not None:
        out["total"] = ([scale * (r[force] + r[inertia]) for r in rows],
                        [scale * (r[force + 1] + r[inertia + 1])
                         for r in rows])
        out["forcing"] = ([scale * r[force] for r in rows],
                          [scale * r[force + 1] for r in rows])
        out["inertia"] = ([scale * r[inertia] for r in rows],
                          [scale * r[inertia + 1] for r in rows])
    else:
        fail("unrecognised force file columns: " + ", ".join(columns))
    return out


def read_series(path: Path) -> tuple:
    t, cd, cl = [], [], []
    for line in path.read_text().splitlines():
        if line.strip() and not line.startswith("#"):
            v = line.split()
            t.append(float(v[0]))
            cd.append(float(v[1]))
            cl.append(float(v[2]))
    return t, cd, cl


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
    w = (x - t[lo]) / (t[hi] - t[lo])
    return y[lo] + w * (y[hi] - y[lo])


def trapezoid(t: list, f: list) -> float:
    return sum(0.5 * (f[i] + f[i + 1]) * (t[i + 1] - t[i])
               for i in range(len(t) - 1))


def window_grid(step: float, window=WINDOW) -> list:
    n = int(round((window[1] - window[0]) / step))
    return [window[0] + i * step for i in range(1, n)]


def nrms(t_sim: list, y_sim: list, t_ref: list, y_ref: list,
         window=WINDOW) -> float:
    """Normalised RMS difference (per cent): the simulated series is
    interpolated on the reference times inside the window, the RMS over
    time by the trapezoid rule, normalised by the RMS of the reference. A
    normalised difference to a reference of uncertain accuracy, not an error"""
    idx = [i for i, x in enumerate(t_ref) if window[0] < x < window[1]]
    tt = [t_ref[i] for i in idx]
    diff2 = [(interpolate(t_sim, y_sim, t_ref[i]) - y_ref[i]) ** 2
             for i in idx]
    ref2 = [y_ref[i] ** 2 for i in idx]
    return 100.0 * math.sqrt(trapezoid(tt, diff2) / trapezoid(tt, ref2))


def solve_linear(a: list, b: list) -> list:
    n = len(b)
    m = [row[:] + [b[i]] for i, row in enumerate(a)]
    for c in range(n):
        p = max(range(c, n), key=lambda r: abs(m[r][c]))
        m[c], m[p] = m[p], m[c]
        for r in range(c + 1, n):
            f = m[r][c] / m[c][c]
            for k in range(c, n + 1):
                m[r][k] -= f * m[c][k]
    x = [0.0] * n
    for r in range(n - 1, -1, -1):
        x[r] = (m[r][n] - sum(m[r][k] * x[k] for k in range(r + 1, n))) \
            / m[r][r]
    return x


def harmonics(t: list, y: list, window=WINDOW, step=0.002) -> dict:
    """Mean, RMS and the amplitudes and phases (degrees, of
    a sin(2 pi f t) + b cos(2 pi f t) = A sin(2 pi f t + phi)) of the
    fundamental (1/period) and second harmonic, by a direct least-squares
    fit over the window, without phase averaging or alignment"""
    grid = window_grid(step, window)
    values = [interpolate(t, y, x) for x in grid]
    basis = [lambda x: 1.0] + [
        (lambda x, k=k, f=f: f(2.0 * math.pi * k * x / PERIOD))
        for k in (1, 2) for f in (math.sin, math.cos)]
    rows = [[b(x) for b in basis] for x in grid]
    ata = [[sum(r[i] * r[j] for r in rows) for j in range(5)]
           for i in range(5)]
    aty = [sum(r[i] * v for r, v in zip(rows, values)) for i in range(5)]
    c = solve_linear(ata, aty)
    # time average and RMS over the whole window, by the trapezoid rule
    full = [window[0]] + grid + [window[1]]
    fv = [interpolate(t, y, x) for x in full]
    span = window[1] - window[0]
    return {"mean": trapezoid(full, fv) / span,
            "rms": math.sqrt(trapezoid(full, [v * v for v in fv]) / span),
            "A1": math.hypot(c[1], c[2]),
            "ph1": math.degrees(math.atan2(c[2], c[1])),
            "A2": math.hypot(c[3], c[4]),
            "ph2": math.degrees(math.atan2(c[4], c[3]))}


def noise(t: list, y: list, width=0.0104, window=WINDOW) -> float:
    """RMS of the series minus its moving mean over about 0.01 s: a
    diagnostic of the high-frequency content, not a force error"""
    dt = (t[-1] - t[0]) / (len(t) - 1)
    k = max(3, int(round(width / dt)) | 1)
    half = k // 2
    acc, count = 0.0, 0
    for i in range(half, len(y) - half):
        if window[0] < t[i] < window[1]:
            acc += (y[i] - sum(y[i - half:i + half + 1]) / k) ** 2
            count += 1
    return math.sqrt(acc / count) if count else float("nan")


# ---------------------------------------------------------------------------
# Flux through the inlet and the outlet
# ---------------------------------------------------------------------------

def patch_flux(phi_file: Path, patches=("inlet", "outlet")) -> float:
    """Sum of the face fluxes phi (m^3/s) on the patches of an ASCII
    surfaceScalarField, outward positive"""
    text = phi_file.read_text(errors="ignore")
    text = text[text.rfind("boundaryField"):]
    total = 0.0
    for patch in patches:
        match = re.search(rf"\n\s*{patch}\s*\{{(.*?)\n\s*\}}", text,
                          re.DOTALL)
        if not match:
            fail(f"patch '{patch}' not found in {phi_file}")
        block = match.group(1)
        uniform = re.search(r"value\s+uniform\s+([-+0-9.eE]+);", block)
        if uniform:
            continue  # a uniform value on an empty patch: not expected here
        values = re.search(r"nonuniform\s+List<scalar>\s*\n?\s*(\d+)\s*\(\s*"
                           r"\n(.*?)\n\s*\)", block, re.DOTALL)
        if not values:
            fail(f"no face values for '{patch}' in {phi_file}")
        total += sum(float(v) for v in values.group(2).split())
    return total


def flux_report(case: Path, params: dict, study: str) -> dict:
    """Net exterior flux at the write times against the exact rate. The
    residual is relative to the analytic peak of Study A; for Study B (exact
    zero) only the absolute flux is meaningful"""
    ay = params["Ay"] if study == STUDY_A else 0.0
    peak = volume_rate_amplitude(params, ay)
    rows = []
    for d in case.iterdir():
        try:
            t = float(d.name)
        except ValueError:
            continue
        on_grid = abs(t / WRITE_INTERVAL - round(t / WRITE_INTERVAL)) < 1e-6
        if t > 0 and on_grid and (d / "phi").exists():
            rows.append((t, patch_flux(d / "phi"), volume_rate(params, t, ay)))
    rows.sort()
    out = {"snapshots": len(rows), "peak_analytic": peak}
    if not rows:
        return out
    out["max_abs_Q"] = max(abs(q) for _, q, _ in rows)
    out["max_abs_residual"] = max(abs(q - e) for _, q, e in rows)
    if peak > 0:
        extremes = [(q, e) for t, q, e in rows
                    if abs(abs(math.sin(4 * math.pi * t)) - 1) < 1e-6]
        out["peak_snapshots"] = len(extremes)
        if extremes:
            out["ratio_at_peaks"] = sum(q * e for q, e in extremes) \
                / sum(e * e for _, e in extremes)
            out["peak_residual_pct"] = 100.0 * sum(
                abs(q - e) for q, e in extremes) / len(extremes) / peak
    return out


# ---------------------------------------------------------------------------
# Cost and log
# ---------------------------------------------------------------------------

def log_stats(case: Path) -> dict:
    log = case / "log.solids4Foam"
    if not log.exists():
        return {}
    text = log.read_text(errors="ignore")
    co = [float(x) for x in re.findall(
        r"Courant Number mean: \S+ max: (\S+)", text)]
    times = re.findall(r"ExecutionTime = ([\d.]+) s\s+ClockTime = (\d+) s",
                       text)
    stats = {"steps": len(re.findall(r"^Time = ", text, re.MULTILINE)),
             "max_Co_in_step": max(co) if co else float("nan"),
             "normal_end": bool(re.search(r"^End\s*$", text, re.MULTILINE))}
    if times:
        execution, clock = float(times[-1][0]), float(times[-1][1])
        stats["cpu_over_wall"] = (execution / clock if clock > 0
                                  else float("nan"))
    return stats


def wall_seconds(case: Path) -> float:
    path = case / "wallclock.txt"
    if path.exists():
        found = re.findall(r"wall_s=(\d+)", path.read_text())
        if found:
            return float(found[-1])
    return float("nan")


def mesh_cells(case: Path, level: int) -> int:
    log = case / "log.blockMesh"
    if log.exists():
        found = re.search(r"nCells:\s*(\d+)", log.read_text())
        if found:
            return int(found.group(1))
    return CELLS[level]


def analyse_case(case: Path, method: str, level: int, study: str,
                 ranks: int, rank_hours: float = None, refs=None,
                 params: dict = None) -> list:
    """Rows (one per estimator and coefficient) of one completed case"""
    params = params or parameters()
    refs = refs if refs is not None else read_series(REFERENCE_FILE)
    t_ref, cd_ref, cl_ref = refs
    columns, rows = read_force_file(case / "postProcessing/immersedBoundary/"
                                    "0/beam.dat")
    series = hydrodynamic_forces(columns, rows, params["coeff_scale"])
    t = series["t"]
    stats = log_stats(case)
    wall = wall_seconds(case)
    if rank_hours is None and not math.isnan(wall):
        rank_hours = ranks * wall / 3600.0
    flux = flux_report(case, params, study)
    cells = mesh_cells(case, level)
    steps = stats.get("steps", len(t))
    cost_us = (rank_hours * 3600e6 / (cells * steps)
               if rank_hours and steps else float("nan"))
    out = []
    for estimator in [k for k in series if k not in ("t", "forcing",
                                                       "inertia")]:
        cd, cl = series[estimator]
        for coef, y, ref in (("Cd", cd, cd_ref), ("Cl", cl, cl_ref)):
            row = {"method": method, "level": level, "study": study,
                   "estimator": estimator, "coefficient": coef,
                   "cells": cells, "ranks": ranks, "steps": steps,
                   "max_Co_in_step": stats.get("max_Co_in_step", float("nan")),
                   "allocated_rank_hours": rank_hours or float("nan"),
                   "rank_us_per_cell_step": cost_us,
                   "cpu_over_wall": stats.get("cpu_over_wall", float("nan")),
                   "end_time": t[-1]}
            if t[-1] >= WINDOW[1] - 0.05 and study == STUDY_A:
                row["nrms_pct"] = nrms(t, y, t_ref, ref)
                row.update({k: v for k, v in harmonics(t, y).items()})
                row["noise"] = noise(t, y)
            row.update({f"flux_{k}": v for k, v in flux.items()})
            out.append(row)
    return out


RESULT_COLUMNS = ("method", "level", "study", "estimator", "coefficient",
                  "cells", "ranks", "steps", "end_time", "max_Co_in_step",
                  "nrms_pct", "mean", "rms", "A1", "ph1", "A2", "ph2",
                  "noise", "flux_snapshots", "flux_peak_analytic",
                  "flux_max_abs_Q", "flux_max_abs_residual",
                  "flux_ratio_at_peaks", "flux_peak_residual_pct",
                  "allocated_rank_hours", "rank_us_per_cell_step",
                  "cpu_over_wall")


def write_results(path: Path, rows: list) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=RESULT_COLUMNS,
                                extrasaction="ignore",
                                lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({k: (f"{v:.6g}" if isinstance(v, float) else v)
                             for k, v in row.items()})


def read_results(path: Path) -> list:
    if not path.exists():
        return []
    with path.open() as handle:
        return list(csv.DictReader(handle))


def completed_cases(manifest: dict) -> list:
    out = []
    for name, entry in sorted(manifest["cases"].items()):
        case = CASES_DIR / name
        if case.exists() and not solver_problem(case):
            out.append((name, entry, case))
    return out


def cmd_analyse(args) -> int:
    manifest = load_manifest()
    rows = []
    for name, entry, case in completed_cases(manifest):
        rows += analyse_case(case, entry["method"], entry["level"],
                             entry["study"], entry["ranks"])
    if not rows:
        print("no completed cases to analyse")
        return 0
    out = WORK_DIR / "results.csv"
    write_results(out, rows)
    shown = [r for r in rows if r["coefficient"] == "Cd"]
    for r in shown:
        nr = f"Cd nRMS {r['nrms_pct']:.2f} %, A2 {r['A2']:.3f}, " \
            if "nrms_pct" in r else ""
        flux = ""
        zero = r.get("flux_max_abs_Q", 1) <= ZERO_FLUX_FRACTION * r.get(
            "flux_peak_analytic", 0)
        if "flux_ratio_at_peaks" in r and not zero:
            flux = f"flux ratio {r['flux_ratio_at_peaks']:.4f}, "
        elif "flux_max_abs_Q" in r:
            flux = f"max |Q| {r['flux_max_abs_Q']:.3g} m^3/s, "
        cost = r["allocated_rank_hours"]
        cost = "cost not recorded" if cost != cost else \
            f"{cost:.3f} allocated rank-h"
        print(f"{r['method']} L{r['level']} {r['study']} "
              f"{r['estimator']}: {nr}{flux}{cost}, {r['steps']} steps")
    print(f"results written to {out}")
    return 0


# ---------------------------------------------------------------------------
# Check: PASS / FAIL / SKIPPED
# ---------------------------------------------------------------------------

PASS, FAIL, SKIPPED = "PASS", "FAIL", "SKIPPED"


def check_analytic(params: dict) -> list:
    """Exact, no OpenFOAM: the closed-form volume rate against the shoelace
    area of the deformed outline"""
    results = []
    peak = volume_rate_amplitude(params)
    worst = max(abs(polygon_volume_rate(params, t) - volume_rate(params, t))
                for t in (0.0625, 0.125, 0.25, 0.375, 0.5, 0.8125))
    results.append((PASS if worst <= ANALYTIC_RTOL * peak else FAIL,
                    "analytic volume rate",
                    f"closed form {peak:.6g} m^3/s; largest difference "
                    f"from the outline area {worst:.2e} m^3/s"))
    results.append((PASS if volume_rate_amplitude(params, 0.0) == 0.0
                    else FAIL, "analytic volume rate, Study B",
                    "A_y = 0 gives zero volume rate"))
    return results


def tracked(rows: list, method: str, level: int, study: str,
            estimator: str, coefficient: str):
    for r in rows:
        if (r["method"], int(r["level"]), r["study"], r["estimator"],
                r["coefficient"]) == (method, level, study, estimator,
                                      coefficient):
            return r
    return None


def check_case(row: dict, tracked_rows: list, strict_flux: bool) -> list:
    """Checks of one analysed row against the tracked results"""
    label = (f"{row['method']} L{row['level']} {row['study']} "
             f"{row['estimator']} {row['coefficient']}")
    out = []
    reference = tracked(tracked_rows, row["method"], row["level"],
                        row["study"], row["estimator"], row["coefficient"])
    if row["coefficient"] == "Cd" and "nrms_pct" in row:
        if reference is None or not reference.get("nrms_pct"):
            out.append((SKIPPED, label + " nRMS",
                        "no tracked result for this case"))
        else:
            d = row["nrms_pct"] - float(reference["nrms_pct"])
            d2 = row["A2"] - float(reference["A2"])
            ok = abs(d) <= CHECK_NRMS and abs(d2) <= CHECK_A2
            out.append((PASS if ok else FAIL, label + " nRMS / A2",
                        f"nRMS {row['nrms_pct']:.2f} (tracked "
                        f"{float(reference['nrms_pct']):.2f}), A2 "
                        f"{row['A2']:.3f} (tracked "
                        f"{float(reference['A2']):.3f}); tolerances "
                        f"{CHECK_NRMS} point and {CHECK_A2}"))
    elif row["coefficient"] == "Cd":
        out.append((SKIPPED, label + " nRMS",
                    f"run ends at {row['end_time']:.3g} s, before the "
                    "comparison window"))
    # the flux is checked once per case, on the Cd row
    if row["coefficient"] == "Cd" and row["estimator"] in (
            "cutLink ST", "total"):
        snaps = row.get("flux_snapshots", 0)
        if not snaps:
            out.append((SKIPPED, label + " flux",
                        "no phi fields at the write times (parallel runs "
                        "need reconstructPar)"))
        elif row["method"] != "cutLink":
            peak = row["flux_peak_analytic"]
            zero = row["flux_max_abs_Q"] <= ZERO_FLUX_FRACTION * peak
            out.append((PASS if zero else FAIL, label + " flux",
                        f"max |Q| {row['flux_max_abs_Q']:.2e} m^3/s "
                        f"(zero to the solver tolerance, "
                        f"{100 * row['flux_max_abs_Q'] / peak:.1e} % of "
                        "the analytic peak)"))
        elif row["study"] == STUDY_A and "flux_ratio_at_peaks" in row:
            ref = tracked(tracked_rows, "cutLink", row["level"], "A",
                          "cutLink ST", "Cd")
            if ref is None or not ref.get("flux_ratio_at_peaks"):
                out.append((SKIPPED, label + " flux", "no tracked ratio"))
            else:
                d = row["flux_ratio_at_peaks"] - float(
                    ref["flux_ratio_at_peaks"])
                out.append((PASS if abs(d) <= CHECK_FLUX_RATIO else FAIL,
                            label + " flux ratio",
                            f"{row['flux_ratio_at_peaks']:.4f} (tracked "
                            f"{float(ref['flux_ratio_at_peaks']):.4f}); "
                            f"tolerance {CHECK_FLUX_RATIO}"))
        elif row["study"] == STUDY_B:
            ref = tracked(tracked_rows, "cutLink", row["level"], "B",
                          "cutLink ST", "Cd")
            if ref is None or not ref.get("flux_max_abs_Q"):
                out.append((SKIPPED, label + " flux", "no tracked bias"))
            else:
                r = row["flux_max_abs_Q"] / float(ref["flux_max_abs_Q"])
                out.append((PASS if abs(r - 1) <= 0.05 else FAIL,
                            label + " flux bias",
                            f"max |Q| {row['flux_max_abs_Q']:.3e} m^3/s "
                            f"(tracked {float(ref['flux_max_abs_Q']):.3e})"))
    return out


def cmd_check(args) -> int:
    results = check_analytic(parameters())
    tracked_rows = read_results(RESULTS_DIR / "cases.csv")
    manifest = load_manifest()
    cases = completed_cases(manifest)
    if not cases:
        results.append((SKIPPED, "simulation checks",
                        "no completed cases; run ./Allverify first"))
    for name, entry, case in cases:
        for row in analyse_case(case, entry["method"], entry["level"],
                                entry["study"], entry["ranks"]):
            results += check_case(row, tracked_rows, True)
    counts = {PASS: 0, FAIL: 0, SKIPPED: 0}
    for status, label, detail in results:
        counts[status] += 1
        print(f"{status:7s} {label}: {detail}")
    print(f"{counts[PASS]} PASS, {counts[FAIL]} FAIL, "
          f"{counts[SKIPPED]} SKIPPED")
    if counts[FAIL] or (args.strict and counts[SKIPPED]):
        return 1
    return 0


def cmd_unit(args) -> int:
    results = check_analytic(parameters())
    counts = {PASS: 0, FAIL: 0, SKIPPED: 0}
    for status, label, detail in results:
        counts[status] += 1
        print(f"{status:7s} {label}: {detail}")
    return 1 if counts[FAIL] else 0


def main(argv=None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    commands = ("run", "generate", "dry-run", "status", "analyse", "check",
                "unit")
    command = argv.pop(0) if argv and argv[0] in commands else "run"
    parser = argparse.ArgumentParser(
        prog="Allverify" + ("" if command == "run" else f" {command}"),
        description=__doc__.split("\n\n")[0])
    if command in ("run", "generate", "dry-run"):
        parser.add_argument("--runs", nargs="+", default=list(DEFAULT_RUNS),
                            help="runs: " + ", ".join(
                                f"{k} ({v[0]})" for k, v in RUNS.items()))
        parser.add_argument("--ranks", type=int, default=0,
                            help="MPI ranks per case (default by level: "
                            + ", ".join(f"L{k} {v}"
                                        for k, v in LEVEL_RANKS.items())
                            + ")")
        parser.add_argument("--force", action="store_true",
                            help="regenerate and rerun existing cases")
        parser.add_argument("--dry-run", action="store_true",
                            help="list the cases only")
    if command in ("run", "check"):
        parser.add_argument("--strict", action="store_true",
                            help="treat SKIPPED as a failure")
    args = parser.parse_args(argv)
    args.runs = getattr(args, "runs", list(DEFAULT_RUNS))
    args.ranks = getattr(args, "ranks", 0)
    args.force = getattr(args, "force", False)
    args.dry_run = getattr(args, "dry_run", False)
    args.strict = getattr(args, "strict", False)
    handlers = {"run": cmd_run, "generate": cmd_generate,
                "dry-run": cmd_dry_run, "status": cmd_status,
                "analyse": cmd_analyse, "check": cmd_check, "unit": cmd_unit}
    return handlers[command](args)


if __name__ == "__main__":
    sys.exit(main())
