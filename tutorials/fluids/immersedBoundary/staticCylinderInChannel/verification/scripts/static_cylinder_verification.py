#!/usr/bin/env python3
"""Verification of the staticCylinderInChannel tutorial: the Schafer and
Turek (1996) 2D-1 benchmark with the four forcing methods of
immersedBoundaryForce on the four MESH_LEVELs of the tutorial.

Each case is a copy of the tutorial under verification/work/<study>/, with
only the settings of the study changed: the forcing method and its settings,
the mesh level, the time step, end and write times, the outer correctors and
the decomposition. Serial cases are run with the tutorial's Allrun; parallel
cases with the solids4foam run functions (solids4Foam::runApplication and
solids4Foam::runParallel). The tutorial itself is never modified.

By default the cases use the frozen settings of the final comparison
(reference/final_method_comparison.csv): deltaT = 0.001 s, endTime = 15 s,
nCorrectors 3, nOuterCorrectors 1 for cutLink, ghostCell and penalty, and
nOuterCorrectors 5 with couplingCoeff 0.2 and the vertexFraction occupancy for
the legacy incremental method. Overriding a setting adds a tag to the case
names, so that the sensitivity cases do not replace the final ones.

Usage:
    ./Allverify [options]            generate, run locally, analyse, check
    ./Allverify generate [options]   generate the cases only
    ./Allverify analyse [--study S]  analyse and check the completed cases
    ./Allverify submit [options]     optional: submit generated cases to SLURM

The default run is cutLink on levels 1-4; --methods all runs the four
methods (16 cases). --dry-run lists the cases and their estimated cost.
"""

import argparse
import csv
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

# -----------------------------------------------------------------------------
# Study settings
# -----------------------------------------------------------------------------

VERIFICATION_DIR = Path(__file__).resolve().parent.parent
TUTORIAL_DIR = VERIFICATION_DIR.parent
WORK_DIR = VERIFICATION_DIR / "work"
OUTPUT_DIR = VERIFICATION_DIR / "postProcessing"
REFERENCE_DIR = VERIFICATION_DIR / "reference"

METHODS = ("cutLink", "ghostCell", "penalty", "incremental")
LEVELS = (1, 2, 3, 4)
DEFAULT_STUDY = "final"
DEFAULT_METHODS = ("cutLink",)

# Frozen settings of the final comparison
DELTA_T = 0.001
END_TIME = 15.0
WRITE_INTERVAL = 5.0
N_CORRECTORS = 3
N_OUTER = {"cutLink": 1, "ghostCell": 1, "penalty": 1, "incremental": 5}
COUPLING_COEFF = 0.2

# MPI ranks of each level with --cores auto (the decomposition of the
# reference runs, capped at the number of available cores)
AUTO_RANKS = {1: 1, 2: 1, 3: 2, 4: 8}

# Optional SLURM requests of each level (the legacy incremental cases are
# the longest: about 6 h on 8 ranks for Level 4)
SLURM_MEMORY = {1: "2G", 2: "2G", 3: "4G", 4: "8G"}
SLURM_TIME = {1: "01:00:00", 2: "03:00:00", 3: "08:00:00", 4: "12:00:00"}

# Schafer and Turek (1996) 2D-1 benchmark ranges
CD_BOUNDS = (5.57, 5.59)
CL_BOUNDS = (0.0104, 0.0110)

# High-precision numerical reference values (FeatFlow; Nabh, 1998)
CD_REF = 5.57953523384
CL_REF = 0.010618948146

# Cd = 2 F/(rho Umean^2 D Lz), with rho = 1, Umean = 0.2, D = Lz = 0.1
COEFF_SCALE = 2.0/(1.0*0.2**2*0.1*0.1)

# The verification check: the drag of cutLink, with the frozen settings, on
# the finest level
CHECK_METHOD = "cutLink"
CHECK_LEVEL = 4

# Settings of each method in constant/fvOptions, written explicitly (the
# defaults of immersedBoundaryForce at f33b244a2), and of the body
METHOD_SETTINGS = {
    "cutLink": (
        [
            ("penaltyCoeff", "1e3"),
            ("pinnedRateCoeff", "100"),
            ("linkCorrection", "yes"),
            ("apertureCoupling", "yes"),
            ("forceEstimator", "momentumExchange"),
            ("writeSurfaceTraction", "yes"),
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
        [
            ("couplingCoeff", "0.2"),
        ],
        [("occupancy", "vertexFraction")],
    ),
}

FOAM_HEADER = """\
/*--------------------------------*- C++ -*----------------------------------*\\
| solids4foam: staticCylinderInChannel verification case                      |
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
    raise RuntimeError(message)


def fmt(value: float) -> str:
    return f"{value:.12g}"


# -----------------------------------------------------------------------------
# Case settings and names
# -----------------------------------------------------------------------------

def case_settings(method: str, args) -> dict:
    """Settings of a case: the frozen settings, with the overrides given"""
    return {
        "deltaT": args.deltaT if args.deltaT else DELTA_T,
        "end_time": args.end_time if args.end_time else END_TIME,
        "write_interval": args.write_interval or WRITE_INTERVAL,
        "nOuterCorrectors": args.n_outer if args.n_outer else N_OUTER[method],
        "nCorrectors": N_CORRECTORS,
        "couplingCoeff": ((args.coupling_coeff if args.coupling_coeff
                           else COUPLING_COEFF)
                          if method == "incremental" else None),
    }


def auto_tag(method: str, args) -> str:
    """Tag of the case names for the settings that differ from the frozen
    ones, e.g. n3-cc0p8-dt0p0005 (empty for the frozen settings)"""
    def short(value):
        return fmt(value).replace(".", "p")
    parts = []
    if args.n_outer and args.n_outer != N_OUTER[method]:
        parts.append(f"n{args.n_outer}")
    if (method == "incremental" and args.coupling_coeff
            and args.coupling_coeff != COUPLING_COEFF):
        parts.append(f"cc{short(args.coupling_coeff)}")
    if args.deltaT and args.deltaT != DELTA_T:
        parts.append(f"dt{short(args.deltaT)}")
    if args.end_time and args.end_time != END_TIME:
        parts.append(f"t{short(args.end_time)}")
    return "-".join(parts)


def case_name(method: str, level: int, tag: str) -> str:
    return f"{method}-L{level}{('-' + tag) if tag else ''}"


def selected_cases(args) -> list:
    """(method, level, name) of the selected cases"""
    cases = []
    for method in args.methods:
        tag = "-".join(t for t in (args.tag, auto_tag(method, args)) if t)
        for level in args.levels:
            cases.append((method, level, case_name(method, level, tag)))
    return cases


def ranks_of(level: int, cores: str) -> int:
    if cores == "auto":
        return max(1, min(AUTO_RANKS[level], os.cpu_count() or 1))
    return int(cores)


# -----------------------------------------------------------------------------
# Case generation
# -----------------------------------------------------------------------------

def replace_entry(path: Path, key: str, value: str) -> None:
    """Replace the value of a top-level 'key value;' entry, exactly once"""
    text = path.read_text()
    pattern = re.compile(rf"^({re.escape(key)}\s+)[^;\n]*;", re.MULTILINE)
    new, count = pattern.subn(lambda m: f"{m.group(1)}{value};", text)
    if count != 1:
        fail(f"found {count} '{key}' entries in {path}")
    path.write_text(new)


def set_pimple_entry(path: Path, key: str, value: int) -> None:
    """Set an entry of the PIMPLE dictionary of fvSolution, exactly once"""
    text = path.read_text()
    pattern = re.compile(rf"^(\s+{key}\s+)\d+;", re.MULTILINE)
    new, count = pattern.subn(lambda m: f"{m.group(1)}{value};", text)
    if count != 1:
        fail(f"found {count} '{key}' entries in {path}")
    path.write_text(new)


def pimple_entry(case: Path, key: str):
    """Value of an entry of the PIMPLE dictionary of fvSolution of a case"""
    m = re.search(rf"^\s+{key}\s+(\d+);",
                  (case / "system/fvSolution").read_text(), re.MULTILINE)
    return int(m.group(1)) if m else ""


def write_fv_options(path: Path, method: str, settings: dict) -> None:
    option, body = METHOD_SETTINGS[method]
    if method == "incremental":
        option = [(k, fmt(settings["couplingCoeff"])
                   if k == "couplingCoeff" else v) for k, v in option]
    width = max(len(k) for k, _ in option + [("method", "")]) + 4
    lines = [FOAM_HEADER.format(location="constant", object="fvOptions")]
    lines.append("// Generated by verification/scripts/"
                 "static_cylinder_verification.py\n")
    lines.append("immersedBoundary\n{\n")
    lines.append(f"    {'type'.ljust(width)}immersedBoundaryForce;\n\n")
    lines.append(f"    {'method'.ljust(width)}{method};\n")
    for key, value in option:
        lines.append(f"    {key.ljust(width)}{value};\n")
    lines.append("\n    bodies\n    {\n        cylinder\n        {\n")
    lines.append(f"            {'surface'.ljust(12)}\"cylinder.stl\";\n")
    for key, value in body:
        lines.append(f"            {key.ljust(12)}{value};\n")
    lines.append("        }\n    }\n}\n\n")
    lines.append("// ****************************************************"
                 "********************* //\n")
    # A regular file in place of the tutorial's link to fvOptions.cutLink,
    # which the case's Allrun, run without a method, keeps
    path.unlink(missing_ok=True)
    path.write_text("".join(lines))


def write_decompose_par_dict(path: Path, ranks: int) -> None:
    text = FOAM_HEADER.format(location="system", object="decomposeParDict")
    text += (
        f"numberOfSubdomains  {ranks};\n\n"
        "method          simple;\n\n"
        "simpleCoeffs\n{\n"
        f"    n           ({ranks} 1 1);\n"
        "}\n\n"
        "// **************************************************************"
        "*********** //\n"
    )
    path.write_text(text)


# Run script of each case: the tutorial's Allrun in serial, the solids4foam
# run functions in parallel; it records the wall time in wallclock.txt
RUN_SCRIPT = """\
#!/bin/bash
# Generated by verification/scripts/static_cylinder_verification.py
cd "${{0%/*}}" || exit 1
start=$(date +%s)
{body}
status=$?
echo "wall_s=$(( $(date +%s) - start )) status=$status" >> wallclock.txt
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


def generate_case(method: str, level: int, name: str, args) -> Path:
    settings = case_settings(method, args)
    ranks = ranks_of(level, args.cores)
    case = WORK_DIR / args.study / name
    if case.exists():
        if not args.force:
            fail(f"{case} exists: use --reuse to keep completed runs or "
                 "--force to regenerate it")
        shutil.rmtree(case)
    case.mkdir(parents=True)

    for item in ("0", "constant", "system", "Allrun", "forceCoeffs.gnuplot"):
        source = TUTORIAL_DIR / item
        if source.is_dir():
            shutil.copytree(source, case / item, symlinks=True)
        else:
            shutil.copy2(source, case / item)

    # Mesh level, written into the dictionary rather than taken from the
    # environment
    replace_entry(case / "system/blockMeshDict", "MESH_LEVEL", str(level))

    control = case / "system/controlDict"
    replace_entry(control, "deltaT", fmt(settings["deltaT"]))
    replace_entry(control, "adjustTimeStep", "no")
    replace_entry(control, "endTime", fmt(settings["end_time"]))
    replace_entry(control, "writeControl", "runTime")
    replace_entry(control, "writeInterval", fmt(settings["write_interval"]))

    write_fv_options(case / "constant/fvOptions", method, settings)
    fv_solution = case / "system/fvSolution"
    set_pimple_entry(fv_solution, "nOuterCorrectors",
                     settings["nOuterCorrectors"])
    set_pimple_entry(fv_solution, "nCorrectors", settings["nCorrectors"])
    write_decompose_par_dict(case / "system/decomposeParDict", ranks)

    script = case / "Allrun.verify"
    script.write_text(RUN_SCRIPT.format(
        body=SERIAL_BODY if ranks == 1 else PARALLEL_BODY))
    script.chmod(0o755)

    sha = subprocess.run(
        ["git", "-C", str(TUTORIAL_DIR), "rev-parse", "HEAD"],
        capture_output=True, text=True, check=False).stdout.strip()
    meta = {
        "study": args.study,
        "case": name,
        "method": method,
        "level": level,
        "cells_per_D": 10*2**(level - 1),
        "h_over_D": 0.1/2**(level - 1),
        "ranks": ranks,
        "decomposition": f"simple ({ranks} 1 1)",
        "deltaT": settings["deltaT"],
        "end_time": settings["end_time"],
        "write_interval": settings["write_interval"],
        "nOuterCorrectors": settings["nOuterCorrectors"],
        "nCorrectors": settings["nCorrectors"],
        "fvOptions": dict(METHOD_SETTINGS[method][0]),
        "body": dict(METHOD_SETTINGS[method][1]),
        "git_sha": sha,
        "generated": datetime.now().isoformat(timespec="seconds"),
    }
    if method == "incremental":
        meta["fvOptions"]["couplingCoeff"] = fmt(settings["couplingCoeff"])
    (case / "case.json").write_text(json.dumps(meta, indent=2) + "\n")
    return case


def completed(case: Path) -> bool:
    log = case / "log.solids4Foam"
    return log.is_file() and bool(re.search(
        r"^End\s*$", log.read_text(errors="replace"), re.MULTILINE))


def solver_problem(case: Path) -> str:
    """Empty if the solver log shows a completed run, else the reason"""
    log = case / "log.solids4Foam"
    if not log.is_file():
        return f"no {log.name}"
    text = log.read_text(errors="replace")
    if re.search(r"FOAM FATAL|FOAM aborting|^ERROR$", text, re.MULTILINE):
        return f"the solver failed; see {log}"
    if not re.search(r"^End\s*$", text, re.MULTILINE):
        return f"the run did not complete; see {log}"
    return ""


def estimate_cpu_hours(cases, args) -> dict:
    """CPU hours of the reference runs with the same method, level and
    settings, where available (machine-dependent)"""
    rows = read_reference_rows()
    estimates = {}
    for method, level, name in cases:
        settings = case_settings(method, args)
        match = find_reference(rows, method, level, settings)
        estimates[name] = float(match["cpu_hours"]) if match else math.nan
    return estimates


def prepare(args, generate_only: bool = False) -> list:
    """Generate the selected cases (or keep the completed ones with --reuse)
    and return (case, name) pairs; with --dry-run, only list them"""
    cases = selected_cases(args)
    estimates = estimate_cpu_hours(cases, args)
    total = sum(v for v in estimates.values() if not math.isnan(v))
    print(f"{'case':30s} {'ranks':>5s} {'est. CPU-h':>11s}")
    for method, level, name in cases:
        est = estimates[name]
        print(f"{name:30s} {ranks_of(level, args.cores):5d} "
              f"{'n/a' if math.isnan(est) else f'{est:.2f}':>11s}")
    print(f"{len(cases)} case(s) in work/{args.study}; about {total:.1f} "
          "CPU-h on the reference machine (2 x AMD EPYC 9684X) for the cases "
          "with a reference run")
    if args.dry_run:
        print("Dry run: nothing generated or run")
        return []
    prepared = []
    for method, level, name in cases:
        case = WORK_DIR / args.study / name
        if args.reuse and case.exists() and (completed(case)
                                             or generate_only):
            print(f"{name}: reusing {case}")
        else:
            generate_case(method, level, name, args)
            print(f"{name}: generated in {case}")
        prepared.append((case, name))
    return prepared


# -----------------------------------------------------------------------------
# Running
# -----------------------------------------------------------------------------

def check_environment() -> None:
    for executable in ("blockMesh", "solids4Foam", "solids4FoamScripts.sh"):
        if shutil.which(executable) is None:
            fail(f"'{executable}' is unavailable: source OpenFOAM and "
                 "solids4foam first")


def run_local(case: Path, name: str) -> None:
    if completed(case):
        print(f"{name}: completed run reused")
        return
    meta = json.loads((case / "case.json").read_text())
    print(f"{name}: running on {meta['ranks']} rank(s) "
          f"({datetime.now():%H:%M:%S})", flush=True)
    with (case / "log.Allverify").open("w") as handle:
        result = subprocess.run(["./Allrun.verify"], cwd=case, stdout=handle,
                                stderr=subprocess.STDOUT, text=True)
    problem = solver_problem(case)
    if result.returncode or problem:
        fail(f"{name}: {problem or 'Allrun.verify failed'}; "
             f"see {case / 'log.Allverify'}")


def cmd_run(args) -> int:
    if not args.dry_run:
        check_environment()
    prepared = prepare(args)
    if not prepared:
        return 0
    start = time.time()
    for case, name in prepared:
        run_local(case, name)
    print(f"All cases completed in {(time.time() - start)/3600:.2f} h")
    return analyse_and_check(args.study)


def cmd_generate(args) -> int:
    prepare(args, generate_only=True)
    return 0


def queue_cpus() -> int:
    """CPUs of the user's running and pending SLURM jobs"""
    out = subprocess.run(
        ["squeue", "--me", "--noheader", "--states=R,PD", "--format=%C"],
        capture_output=True, text=True, check=True).stdout
    return sum(int(x) for x in out.split() if x.isdigit())


def cmd_submit(args) -> int:
    """Optional: submit generated cases to SLURM, one job per case running
    its Allrun.verify; only the summary is printed without --yes"""
    jobs = []
    for method, level, name in selected_cases(args):
        case = WORK_DIR / args.study / name
        if not (case / "Allrun.verify").exists():
            fail(f"{case} has not been generated (./Allverify generate)")
        if completed(case):
            print(f"{name}: completed run, not submitted")
            continue
        meta = json.loads((case / "case.json").read_text())
        jobs.append((case, meta))
    requested = sum(m["ranks"] for _, m in jobs)
    print(f"{'case':30s} {'ranks':>5s} {'memory':>7s} {'time':>9s}")
    for case, m in jobs:
        print(f"{m['case']:30s} {m['ranks']:5d} "
              f"{args.mem or SLURM_MEMORY[m['level']]:>7s} "
              f"{args.time or SLURM_TIME[m['level']]:>9s}")
    print(f"{len(jobs)} job(s), {requested} CPU(s)")
    if args.max_cpus:
        if shutil.which("squeue") is None:
            fail("--max-cpus needs squeue")
        in_queue = queue_cpus()
        print(f"{in_queue} CPU(s) already queued; limit {args.max_cpus}")
        if in_queue + requested > args.max_cpus:
            fail("the CPU limit given by --max-cpus would be exceeded")
    if not args.yes:
        print("Not submitted: re-run with --yes to submit")
        return 0
    if shutil.which("sbatch") is None:
        fail("sbatch is unavailable: run the cases locally with ./Allverify")
    setup = (f"source {shlex.quote(args.env_script)} && "
             if args.env_script else "")
    for case, m in jobs:
        command = (f"{setup}cd {shlex.quote(str(case))} && "
                   "./Allrun.verify > log.Allverify 2>&1")
        sbatch = ["sbatch", "--parsable", f"--job-name=scv-{m['case']}",
                  f"--ntasks={m['ranks']}", "--cpus-per-task=1",
                  f"--mem={args.mem or SLURM_MEMORY[m['level']]}",
                  f"--time={args.time or SLURM_TIME[m['level']]}",
                  f"--output={case / 'log.slurm'}", f"--chdir={case}"]
        if args.partition:
            sbatch.append(f"--partition={args.partition}")
        sbatch += shlex.split(args.sbatch_args)
        sbatch += ["--wrap", f"bash -c {shlex.quote(command)}"]
        job = subprocess.run(sbatch, capture_output=True, text=True,
                             check=True).stdout.strip()
        print(f"{m['case']}: SLURM job {job}")
    return 0


# -----------------------------------------------------------------------------
# Analysis
# -----------------------------------------------------------------------------

def force_files(case: Path) -> list:
    """Force files of the run and its restarts, in order of start time"""
    base = case / "postProcessing/immersedBoundary"
    files = [d / "cylinder.dat" for d in base.iterdir()
             if (d / "cylinder.dat").exists()] if base.exists() else []
    return sorted(files, key=lambda f: float(f.parent.name))


def read_forces(files: list):
    """Time and force; on a restart, the rows of the later file replace those
    of the earlier files from its start time"""
    data = {}
    for path in files:
        start = float(path.parent.name)
        for key in [k for k in data if k >= start - 1e-12]:
            del data[key]
        for line in path.read_text().splitlines():
            if not line.strip() or line.startswith("#"):
                continue
            cols = line.split()
            data[float(cols[0])] = (float(cols[1]), float(cols[2]))
    t = sorted(data)
    return t, [data[k][0] for k in t], [data[k][1] for k in t]


def window_stats(t, v, t0):
    """Time-weighted (trapezoidal) mean, range and least-squares slope of v
    over t >= t0"""
    idx = [i for i, ti in enumerate(t) if ti >= t0 - 1e-12]
    if len(idx) < 3:
        return None
    area = duration = 0.0
    for a, b in zip(idx, idx[1:]):
        dt = t[b] - t[a]
        area += 0.5*(v[a] + v[b])*dt
        duration += dt
    mean = area/duration
    tm = sum(t[i] for i in idx)/len(idx)
    vm = sum(v[i] for i in idx)/len(idx)
    denominator = sum((t[i] - tm)**2 for i in idx)
    slope = (sum((t[i] - tm)*(v[i] - vm) for i in idx)/denominator
             if denominator > 0 else 0.0)
    values = [v[i] for i in idx]
    return mean, min(values), max(values), slope


def range_ratio(t, v):
    """Range of v over the last second divided by that over the second
    before (about 1 for an oscillation of stationary amplitude)"""
    end = t[-1]
    last = [x for ti, x in zip(t, v) if ti >= end - 1.0 - 1e-12]
    prev = [x for ti, x in zip(t, v) if end - 2.0 - 1e-12 <= ti < end - 1.0]
    if len(last) < 2 or len(prev) < 2:
        return math.nan
    r_prev = max(prev) - min(prev)
    return (max(last) - min(last))/r_prev if r_prev > 0 else math.nan


def classify(t, cd, cl, end_time, window=2.0, steady_tol=1e-5,
             cl_tol=1e-5) -> dict:
    """Statistics of Cd and Cl over the final window and the status of the
    run: steady, limit cycle (constant mean, stationary oscillation),
    drifting or incomplete"""
    t0 = t[-1] - window
    scd = window_stats(t, cd, t0)
    scl = window_stats(t, cl, t0)
    complete = t[-1] >= end_time - 1e-9
    cd_mean = scd[0] if scd else math.nan
    cd_rel_range = (scd[2] - scd[1])/abs(cd_mean) if scd else math.nan
    cd_rel_drift = abs(scd[3])*10.0/abs(cd_mean) if scd else math.nan
    cl_range = scl[2] - scl[1] if scl else math.nan
    steady = bool(complete and scd and scl and cd_rel_range <= steady_tol
                  and cd_rel_drift <= steady_tol and cl_range <= cl_tol)
    if steady:
        status = "steady"
    elif (complete and scd and scl and cd_rel_drift <= steady_tol
          and (cd_rel_range <= steady_tol
               or 0.9 <= range_ratio(t, cd) <= 1.1)
          and (cl_range <= cl_tol or 0.9 <= range_ratio(t, cl) <= 1.1)):
        status = "limit cycle"
    else:
        status = "drifting" if complete else "incomplete"
    return {
        "end_time": t[-1], "Cd_mean": cd_mean,
        "Cd_range": (scd[2] - scd[1]) if scd else math.nan,
        "Cd_rel_drift_10s": cd_rel_drift,
        "Cl_mean": scl[0] if scl else math.nan, "Cl_range": cl_range,
        "complete": complete, "status": status,
    }


def log_stats(case: Path) -> dict:
    """Maximum Courant number and pressure iterations per time step of the
    run and its restarts"""
    logs = sorted(case.glob("log.solids4Foam*"),
                  key=lambda f: (f.name != "log.solids4Foam", f.name))
    max_co, steps, p_iters, corr = 0.0, 0, [], []
    for log in logs:
        for line in log.read_text(errors="replace").splitlines():
            if line.startswith("Courant Number mean:"):
                max_co = max(max_co, float(line.split()[5]))
            elif line.startswith("Time = "):
                steps += 1
                if corr:
                    p_iters.append(sum(corr))
                corr = []
            elif "Solving for p," in line:
                m = re.search(r"No Iterations (\d+)", line)
                if m:
                    corr.append(int(m.group(1)))
    if corr:
        p_iters.append(sum(corr))
    return {"max_co": max_co, "n_steps": steps,
            "p_iters_per_step": (sum(p_iters)/len(p_iters)) if p_iters
            else math.nan}


def wall_seconds(case: Path) -> float:
    """Wall time of the run (and its restarts) from wallclock.txt"""
    path = case / "wallclock.txt"
    if not path.exists():
        return math.nan
    walls = re.findall(r"wall_s=(\d+)", path.read_text())
    return float(sum(int(w) for w in walls)) if walls else math.nan


def case_row(case: Path):
    meta_path = case / "case.json"
    files = force_files(case)
    if not meta_path.exists() or not files:
        return None
    meta = json.loads(meta_path.read_text())
    t, fx, fy = read_forces(files)
    stats = classify(t, [COEFF_SCALE*f for f in fx],
                     [COEFF_SCALE*f for f in fy], meta["end_time"])
    wall = wall_seconds(case)
    row = {
        "study": meta["study"], "case": meta["case"],
        "method": meta["method"], "level": meta["level"],
        "cells_per_D": meta["cells_per_D"], "ranks": meta["ranks"],
        "couplingCoeff": meta["fvOptions"].get("couplingCoeff", ""),
        "occupancy": meta["body"].get("occupancy", ""),
        "deltaT": meta["deltaT"],
        "nOuterCorrectors": pimple_entry(case, "nOuterCorrectors"),
        "nCorrectors": pimple_entry(case, "nCorrectors"),
        "end_time_target": meta["end_time"],
    }
    row["frozen_settings"] = is_frozen(row)
    row.update(stats)
    row.update(log_stats(case))
    row["Cd_rel_deviation"] = stats["Cd_mean"]/CD_REF - 1.0
    row["Cl_rel_deviation"] = stats["Cl_mean"]/CL_REF - 1.0
    row["wall_s"] = wall
    row["cpu_hours"] = wall*meta["ranks"]/3600.0
    return row


def is_frozen(row: dict) -> bool:
    """Whether a case used the frozen settings of the final comparison"""
    method = row["method"]
    try:
        return (int(row["nOuterCorrectors"]) == N_OUTER[method]
                and int(row["nCorrectors"]) == N_CORRECTORS
                and abs(float(row["deltaT"]) - DELTA_T) < 1e-12
                and abs(float(row["end_time_target"]) - END_TIME) < 1e-9
                and (method != "incremental"
                     or abs(float(row["couplingCoeff"]) - COUPLING_COEFF)
                     < 1e-12))
    except (KeyError, ValueError, TypeError):
        return False


def study_rows(study: str) -> list:
    rows = []
    for case in sorted((WORK_DIR / study).glob("*-L*")):
        row = case_row(case)
        if row:
            rows.append(row)
    return rows


def read_reference_rows() -> list:
    rows = []
    for path in sorted(REFERENCE_DIR.glob("*.csv")):
        with path.open() as handle:
            lines = [line for line in handle if not line.startswith("#")]
        for row in csv.DictReader(lines):
            if "level" in row and "method" in row:
                row["_file"] = path.name
                rows.append(row)
    return rows


def find_reference(rows: list, method: str, level: int, settings: dict):
    """The reference row with the same method, level and settings"""
    for row in rows:
        try:
            if (row["method"] == method and int(row["level"]) == level
                    and abs(float(row["deltaT"]) - settings["deltaT"]) < 1e-12
                    and int(row["nOuterCorrectors"])
                    == settings["nOuterCorrectors"]
                    and abs(float(row.get("end_time") or 0)
                            - settings["end_time"]) < 1e-9
                    and (method != "incremental"
                         or abs(float(row["couplingCoeff"])
                                - settings["couplingCoeff"]) < 1e-12)):
                return row
        except (KeyError, ValueError):
            continue
    return None


def verification_check(rows: list) -> tuple:
    """The verification check on the analysed cases: the time-averaged drag
    of cutLink, with the frozen settings, on Level 4 within the Schafer and
    Turek (1996) range. Returns (result, message), result PASS, FAIL or
    SKIPPED"""
    case = next((r for r in rows if r["method"] == CHECK_METHOD
                 and int(r["level"]) == CHECK_LEVEL
                 and r.get("frozen_settings")), None)
    label = f"{CHECK_METHOD} L{CHECK_LEVEL} (frozen settings)"
    if case is None:
        return "SKIPPED", f"{label} was not run"
    if not case["complete"] or math.isnan(case["Cd_mean"]):
        return "SKIPPED", f"{label} is incomplete"
    lo, hi = CD_BOUNDS
    cd = case["Cd_mean"]
    if lo <= cd <= hi:
        return "PASS", (f"{label}: Cd = {cd:.6f} is within the Schafer and "
                        f"Turek range {lo} <= Cd <= {hi}")
    return "FAIL", (f"{label}: Cd = {cd:.6f} is outside the Schafer and "
                    f"Turek range {lo} <= Cd <= {hi}")


def analyse_and_check(study: str) -> int:
    rows = study_rows(study)
    out_dir = OUTPUT_DIR / study
    out_dir.mkdir(parents=True, exist_ok=True)
    references = read_reference_rows()
    for row in rows:
        settings = {"deltaT": float(row["deltaT"]),
                    "nOuterCorrectors": int(row["nOuterCorrectors"]),
                    "end_time": float(row["end_time_target"]),
                    "couplingCoeff": float(row["couplingCoeff"] or 0)}
        ref = find_reference(references, row["method"], int(row["level"]),
                             settings)
        row["Cd_reference_run"] = float(ref["Cd_mean"]) if ref else math.nan
        row["reference_file"] = ref["_file"] if ref else ""
    if rows:
        with (out_dir / "cases.csv").open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
            writer.writeheader()
            writer.writerows(rows)
    result, message = verification_check(rows)

    lines = [f"# staticCylinderInChannel verification: {study}", "",
             f"- Result: {result}", f"- {message}", ""]
    if rows:
        lines += ["| Case | Cd | Cd dev. | Cl | Status | Ref. run Cd |",
                  "| --- | --- | --- | --- | --- | --- |"]
        for r in rows:
            ref = ("" if math.isnan(r["Cd_reference_run"])
                   else f"{r['Cd_reference_run']:.6f}")
            lines.append(f"| {r['case']} | {r['Cd_mean']:.6f} | "
                         f"{100*r['Cd_rel_deviation']:+.3f}% | "
                         f"{r['Cl_mean']:.6f} | {r['status']} | {ref} |")
        lines += ["", "Cd dev. is the relative deviation from the numerical "
                  f"reference Cd = {CD_REF}; Ref. run Cd is the tracked "
                  "reference run with the same settings, if any.", ""]
    (out_dir / "verification_summary.md").write_text("\n".join(lines))

    for r in rows:
        print(f"{r['case']:30s} Cd={r['Cd_mean']:.6f} "
              f"({100*r['Cd_rel_deviation']:+.3f}%) Cl={r['Cl_mean']:.6f} "
              f"{r['status']}")
    print(f"Verification: {result}: {message}")
    print(f"Results: {out_dir / 'cases.csv'}, "
          f"{out_dir / 'verification_summary.md'}")
    return 1 if result == "FAIL" else 0


def cmd_analyse(args) -> int:
    if not (WORK_DIR / args.study).is_dir():
        print(f"No cases in {WORK_DIR / args.study}")
    return analyse_and_check(args.study)


# -----------------------------------------------------------------------------
# Command line
# -----------------------------------------------------------------------------

def positive_cores(value: str) -> str:
    if value != "auto" and (not value.isdecimal() or int(value) < 1):
        raise argparse.ArgumentTypeError("a positive integer or auto")
    return value


def method_list(values: list) -> list:
    if values == ["all"]:
        return list(METHODS)
    for value in values:
        if value not in METHODS:
            raise SystemExit(f"--methods: unknown method {value}; choose from "
                             f"{', '.join(METHODS)} or all")
    return values


def selection_options(p, run_options: bool = True) -> None:
    p.add_argument("--study", default=DEFAULT_STUDY,
                   help=f"sub-directory of work/ (default {DEFAULT_STUDY})")
    p.add_argument("--methods", nargs="+", default=list(DEFAULT_METHODS),
                   metavar="METHOD",
                   help="cutLink, ghostCell, penalty, incremental or all "
                   "(default cutLink)")
    p.add_argument("--levels", nargs="+", type=int, default=list(LEVELS),
                   choices=LEVELS, help="mesh levels (default 1 2 3 4)")
    p.add_argument("--tag", default="", help="extra suffix of the case names")
    p.add_argument("--n-outer", type=int, default=0,
                   help="override nOuterCorrectors (frozen: 1, and 5 for "
                   "incremental)")
    p.add_argument("--coupling-coeff", type=float, default=0.0,
                   help="override couplingCoeff of incremental (frozen 0.2)")
    p.add_argument("--deltaT", type=float, default=0.0,
                   help=f"override the time step (frozen {DELTA_T} s)")
    p.add_argument("--end-time", type=float, default=0.0,
                   help=f"override the end time (frozen {END_TIME:g} s)")
    if run_options:
        p.add_argument("--write-interval", type=float, default=0.0,
                       help=f"write interval (default {WRITE_INTERVAL:g} s)")
        p.add_argument("--cores", type=positive_cores, default="auto",
                       help="MPI ranks per case: a positive integer or auto "
                       "(default: 1, 1, 2 and 8 on levels 1-4, at most the "
                       "available cores)")
        p.add_argument("--reuse", action="store_true",
                       help="keep the completed cases found in work/")
        p.add_argument("--force", action="store_true",
                       help="regenerate existing cases")
        p.add_argument("--dry-run", action="store_true",
                       help="list the cases and their cost; run nothing")


def main(argv=None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    commands = ("run", "generate", "analyse", "submit")
    command = argv.pop(0) if argv and argv[0] in commands else "run"
    parser = argparse.ArgumentParser(
        prog=f"Allverify {command}" if command != "run" else "Allverify",
        description=__doc__.split("\n\n")[0])
    if command in ("run", "generate"):
        selection_options(parser)
    elif command == "analyse":
        parser.add_argument("--study", default=DEFAULT_STUDY,
                            help="sub-directory of work/ (default "
                            f"{DEFAULT_STUDY})")
    else:
        selection_options(parser, run_options=False)
        parser.add_argument("--env-script",
                            default=os.environ.get("S4F_ENV_SCRIPT", ""),
                            help="script sourced in each job to set up "
                            "OpenFOAM and solids4foam (default "
                            "$S4F_ENV_SCRIPT; none: the job inherits the "
                            "submitting environment)")
        parser.add_argument("--partition", default="",
                            help="SLURM partition (default: the default)")
        parser.add_argument("--max-cpus", type=int, default=0,
                            help="refuse to submit beyond this number of "
                            "queued CPUs of the user (default: no limit)")
        parser.add_argument("--time", default="",
                            help="wall time of each job (default by level)")
        parser.add_argument("--mem", default="",
                            help="memory of each job (default by level)")
        parser.add_argument("--sbatch-args", default="",
                            help="further sbatch options, e.g. site settings")
        parser.add_argument("--yes", action="store_true",
                            help="submit (otherwise only the summary)")
    args = parser.parse_args(argv)
    if hasattr(args, "methods"):
        args.methods = method_list(args.methods)
    if command == "submit":
        args.write_interval = 0.0
        args.cores = "auto"
    return {"run": cmd_run, "generate": cmd_generate,
            "analyse": cmd_analyse, "submit": cmd_submit}[command](args)


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (RuntimeError, subprocess.SubprocessError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        sys.exit(2)
