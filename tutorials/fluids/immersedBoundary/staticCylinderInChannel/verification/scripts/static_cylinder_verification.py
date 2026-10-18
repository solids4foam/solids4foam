#!/usr/bin/env python3
"""
Verification driver for the staticCylinderInChannel tutorial: the
Schafer and Turek (1996) 2D-1 benchmark with the four forcing methods of
immersedBoundaryForce on the four MESH_LEVELs of the tutorial.

Each case is generated from the canonical tutorial under verification/work/,
with only the settings of the study changed: the forcing method (and the
settings it requires), the mesh level, a fixed time step common to all the
methods and levels, the end and write times, and the decomposition. The
tutorial itself is never modified.

Sub-commands:
    generate   create the cases and their SLURM scripts
    submit     submit the generated cases (after a resource summary and an
               explicit --yes), within the aggregate CPU limit
    analyse    extract Cd and Cl, the steady-state statistics, the Courant
               number and the cost of the cases into runs.csv; with
               --report, also convergence.csv, summary.md and the plots
    extend     restart cases from their end time with the same time step
               and ranks, and extend them (default by 5 s)
"""

import argparse
import csv
import json
import math
import os
import re
import shutil
import subprocess
import sys
from datetime import datetime
from pathlib import Path

# -----------------------------------------------------------------------------
# Study settings
# -----------------------------------------------------------------------------

VERIFICATION_DIR = Path(__file__).resolve().parent.parent
TUTORIAL_DIR = VERIFICATION_DIR.parent
WORK_DIR = VERIFICATION_DIR / "work"
OUTPUT_DIR = VERIFICATION_DIR / "postProcessing"

METHODS = ("cutLink", "ghostCell", "penalty", "incremental")
LEVELS = (1, 2, 3, 4)

# Common fixed time step (s), chosen by the Level 4 Courant number
# calibration (see README.md)
DELTA_T = 0.001

END_TIME = 10.0
WRITE_INTERVAL = 5.0

# MPI ranks at each mesh level, the same for all the methods, from the
# Level 4 scaling calibration (see README.md)
RANKS = {1: 1, 2: 1, 3: 2, 4: 8}

# SLURM requests at each mesh level: memory and wall time, with a margin of
# about four over the calibration estimate for 10 s at DELTA_T (Level 4 with
# 4 ranks: 0.28-0.42 s per step, about 0.8-1.2 h)
MEMORY = {1: "2G", 2: "2G", 3: "4G", 4: "8G"}
WALLTIME = {1: "01:00:00", 2: "01:00:00", 3: "03:00:00", 4: "06:00:00"}
# Site settings, from the environment:
#   S4F_VERIFY_PARTITION  SLURM partition (default: the default partition)
#   S4F_VERIFY_MAX_CPUS   aggregate CPU limit of the user's jobs (default 96)
#   S4F_VERIFY_ENV        script sourced by the jobs to set up OpenFOAM and
#                         solids4foam (default: none, the jobs inherit the
#                         environment of the submitting shell)
PARTITION = os.environ.get("S4F_VERIFY_PARTITION", "")
MAX_CPUS = int(os.environ.get("S4F_VERIFY_MAX_CPUS", "96"))
ENV_SH = os.environ.get("S4F_VERIFY_ENV", "")

# Benchmark reference values (Schafer and Turek, 1996, bounds)
CD_BOUNDS = (5.57, 5.59)
CL_BOUNDS = (0.0104, 0.0110)

# High-precision reference values of the benchmark (FEATFLOW)
CD_REF = 5.57953523384
CL_REF = 0.010618948146

# Cd = 2 F/(rho Umean^2 D Lz), with rho = 1, Umean = 0.2, D = Lz = 0.1
COEFF_SCALE = 2.0/(1.0*0.2**2*0.1*0.1)

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
            ("couplingCoeff", "0.8"),
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


def case_name(method: str, level: int, tag: str = "") -> str:
    return f"{method}-L{level}{('-' + tag) if tag else ''}"


def fmt(value: float) -> str:
    return f"{value:.12g}"


def replace_entry(path: Path, key: str, value: str) -> None:
    """Replace the value of a top-level 'key value;' entry, exactly once"""
    text = path.read_text()
    pattern = re.compile(rf"^({re.escape(key)}\s+)[^;\n]*;", re.MULTILINE)
    new, count = pattern.subn(lambda m: f"{m.group(1)}{value};", text)
    if count != 1:
        sys.exit(f"ERROR: found {count} '{key}' entries in {path}")
    path.write_text(new)


def method_settings(method: str, coupling_coeff: float = 0.0):
    """Settings of the method in fvOptions and of the body, with the
    coupling coefficient of the incremental method overridden if given"""
    option, body = METHOD_SETTINGS[method]
    if coupling_coeff > 0:
        if method != "incremental":
            sys.exit("ERROR: --coupling-coeff is only used by incremental")
        option = [(k, fmt(coupling_coeff) if k == "couplingCoeff" else v)
                  for k, v in option]
    return option, body


def set_pimple_entry(path: Path, key: str, value: int) -> None:
    """Set an entry of the PIMPLE dictionary of fvSolution, exactly once"""
    text = path.read_text()
    pattern = re.compile(rf"^(\s+{key}\s+)\d+;", re.MULTILINE)
    new, count = pattern.subn(lambda m: f"{m.group(1)}{value};", text)
    if count != 1:
        sys.exit(f"ERROR: found {count} '{key}' entries in {path}")
    path.write_text(new)


def pimple_entry(case: Path, key: str):
    """Value of an entry of the PIMPLE dictionary of fvSolution of a case"""
    m = re.search(rf"^\s+{key}\s+(\d+);",
                  (case / "system/fvSolution").read_text(), re.MULTILINE)
    return int(m.group(1)) if m else ""


def write_fv_options(path: Path, method: str,
                     coupling_coeff: float = 0.0) -> None:
    option, body = method_settings(method, coupling_coeff)
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


def slurm_header(case: Path, job_name: str, ranks: int, memory: str,
                 walltime: str, partition: str, env_sh: str) -> str:
    """SLURM directives and environment set-up of a job script"""
    lines = [
        "#!/bin/bash",
        f"#SBATCH --job-name={job_name}",
    ]
    if partition:
        lines.append(f"#SBATCH --partition={partition}")
    lines += [
        "#SBATCH --nodes=1",
        f"#SBATCH --ntasks={ranks}",
        "#SBATCH --cpus-per-task=1",
        f"#SBATCH --mem={memory}",
        f"#SBATCH --time={walltime}",
        # With an environment script the job starts from a clean
        # environment; otherwise it inherits that of the submitting shell
        f"#SBATCH --export={'NONE' if env_sh else 'ALL'}",
        f"#SBATCH --chdir={case}",
        "#SBATCH --output=slurm-%j.out",
        "#SBATCH --error=slurm-%j.err",
        "# Generated by verification/scripts/static_cylinder_verification.py",
        "set -euo pipefail",
        'export USER="${USER:-$(id -un)}"',
        'export HOME="${HOME:-$(getent passwd "$USER" | cut -d: -f6)}"',
        "",
    ]
    if env_sh:
        lines += [
            f"source {env_sh}",
            "if declare -F s4f_check_binaries > /dev/null",
            "then",
            "    s4f_check_binaries",
            "fi",
        ]
    lines += [
        'command -v solids4Foam > /dev/null || '
        '{ echo "solids4Foam not found" >&2; exit 1; }',
        'echo "solids4Foam=$(command -v solids4Foam)"',
        "",
    ]
    return "\n".join(lines) + "\n"


def mpi_run(ranks: int) -> str:
    # SLURM may not confine the jobs to their cores (TaskPlugin=task/none):
    # binding with mpirun would then place the ranks of concurrent jobs on
    # the same cores, so the placement is left to the kernel
    return "solids4Foam" if ranks == 1 else (
        f"mpirun -np {ranks} --bind-to none solids4Foam -parallel")


def write_run_script(case: Path, name: str, level: int, ranks: int,
                     memory: str, walltime: str, partition: str,
                     env_sh: str) -> None:
    decompose = "" if ranks == 1 else (
        "decomposePar -force > log.decomposePar 2>&1\n")
    script = slurm_header(case, f"scv-{name}", ranks, memory, walltime,
                          partition, env_sh) + f"""\
echo "job=${{SLURM_JOB_ID:-}} host=$(hostname) ranks={ranks} start=$(date -Is)"
START=$(date +%s)

blockMesh > log.blockMesh 2>&1
{decompose}{mpi_run(ranks)} > log.solids4Foam 2>&1

END=$(date +%s)
echo "end=$(date -Is) wall_s=$((END - START))" | tee wallclock.txt
"""
    path = case / "run.slurm"
    path.write_text(script)
    path.chmod(0o755)


def generate_case(method: str, level: int, args) -> Path:
    ranks = args.ranks if args.ranks else RANKS[level]
    name = case_name(method, level, args.tag)
    case = WORK_DIR / (args.study) / name
    if case.exists():
        if not args.force:
            sys.exit(f"ERROR: {case} exists (use --force to regenerate)")
        shutil.rmtree(case)
    case.mkdir(parents=True)

    for item in ("0", "constant", "system"):
        shutil.copytree(TUTORIAL_DIR / item, case / item, symlinks=True)

    # Mesh level, written into the dictionary rather than taken from the
    # environment
    replace_entry(case / "system/blockMeshDict", "MESH_LEVEL", str(level))

    # Fixed time step
    control = case / "system/controlDict"
    replace_entry(control, "deltaT", fmt(args.deltaT))
    replace_entry(control, "adjustTimeStep", "no")
    replace_entry(control, "endTime", fmt(args.end_time))
    replace_entry(control, "writeControl", "runTime")
    replace_entry(control, "writeInterval", fmt(args.write_interval))

    write_fv_options(case / "constant/fvOptions", method,
                     args.coupling_coeff)
    if args.n_outer:
        set_pimple_entry(case / "system/fvSolution", "nOuterCorrectors",
                         args.n_outer)
    write_decompose_par_dict(case / "system/decomposeParDict", ranks)

    write_run_script(
        case, name, level, ranks,
        args.memory or MEMORY[level], args.walltime or WALLTIME[level],
        args.partition, args.env_sh
    )

    sha = subprocess.run(
        ["git", "-C", str(TUTORIAL_DIR), "rev-parse", "HEAD"],
        capture_output=True, text=True, check=False
    ).stdout.strip()
    dirty = subprocess.run(
        ["git", "-C", str(TUTORIAL_DIR), "status", "--porcelain", "--", "."],
        capture_output=True, text=True, check=False
    ).stdout.strip()

    meta = {
        "study": args.study,
        "case": name,
        "method": method,
        "level": level,
        "cells_per_D": 10*2**(level - 1),
        "h_over_D": 0.1/2**(level - 1),
        "n_cells_expected": (60 + 50)*41*4**(level - 1),
        "ranks": ranks,
        "decomposition": f"simple ({ranks} 1 1)",
        "deltaT": args.deltaT,
        "nOuterCorrectors": pimple_entry(case, "nOuterCorrectors"),
        "nCorrectors": pimple_entry(case, "nCorrectors"),
        "end_time": args.end_time,
        "write_interval": args.write_interval,
        "fvOptions": dict(method_settings(method, args.coupling_coeff)[0]),
        "body": dict(method_settings(method, args.coupling_coeff)[1]),
        "partition": args.partition,
        "memory": args.memory or MEMORY[level],
        "walltime": args.walltime or WALLTIME[level],
        "env_sh": args.env_sh,
        "git_sha": sha,
        "tutorial_modified": bool(dirty),
        "generated": datetime.now().isoformat(timespec="seconds"),
    }
    (case / "case.json").write_text(json.dumps(meta, indent=2) + "\n")
    return case


def selected_cases(args) -> list:
    return [(m, lv) for m in args.methods for lv in args.levels]


def cmd_generate(args) -> None:
    if not args.env_sh and shutil.which("solids4Foam") is None:
        sys.exit("ERROR: solids4Foam is not on the PATH and no environment "
                 "script is given: source OpenFOAM and solids4foam, or set "
                 "S4F_VERIFY_ENV (or --env-sh)")
    if args.env_sh and not Path(args.env_sh).is_file():
        sys.exit(f"ERROR: environment script {args.env_sh} not found")
    for method, level in selected_cases(args):
        case = generate_case(method, level, args)
        meta = json.loads((case / "case.json").read_text())
        print(f"{case.relative_to(VERIFICATION_DIR)}: ranks={meta['ranks']}"
              f" deltaT={meta['deltaT']} endTime={meta['end_time']}")


def queue_cpus() -> int:
    out = subprocess.run(
        ["squeue", "-h", "-u", subprocess.run(
            ["id", "-un"], capture_output=True, text=True).stdout.strip(),
         "-o", "%C"],
        capture_output=True, text=True, check=True
    ).stdout.split()
    return sum(int(c) for c in out)


def cmd_submit(args) -> None:
    cases = []
    for method, level in selected_cases(args):
        case = WORK_DIR / args.study / case_name(method, level, args.tag)
        if not (case / "run.slurm").exists():
            sys.exit(f"ERROR: {case} has not been generated")
        meta = json.loads((case / "case.json").read_text())
        cases.append((case, meta))

    in_queue = queue_cpus()
    requested = sum(m["ranks"] for _, m in cases)
    print(f"{'case':28s} {'ranks':>5s} {'memory':>7s} {'walltime':>9s}")
    for case, m in cases:
        print(f"{m['case']:28s} {m['ranks']:5d} {m['memory']:>7s} "
              f"{m['walltime']:>9s}")
    print(f"jobs={len(cases)} cpus_requested={requested} "
          f"cpus_already_in_queue={in_queue} limit={MAX_CPUS}")

    if in_queue + requested > MAX_CPUS:
        sys.exit("ERROR: the aggregate CPU limit would be exceeded")
    if not args.yes:
        print("Not submitted: re-run with --yes to submit")
        return
    for case, _ in cases:
        out = subprocess.run(["sbatch", "run.slurm"], cwd=case,
                             capture_output=True, text=True, check=True)
        print(f"{case.name}: {out.stdout.strip()}")


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
    """Time, force and surface-traction drag; on a restart, the rows of the
    later file replace those of the earlier files from its start time"""
    data = {}
    for path in files:
        start = float(path.parent.name)
        for key in [k for k in data if k >= start - 1e-12]:
            del data[key]
        for line in path.read_text().splitlines():
            if not line.strip() or line.startswith("#"):
                continue
            cols = line.split()
            data[float(cols[0])] = (
                float(cols[1]), float(cols[2]),
                float(cols[10]) if len(cols) >= 13 else math.nan
            )
    t = sorted(data)
    return (t, [data[k][0] for k in t], [data[k][1] for k in t],
            [data[k][2] for k in t])


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
    slope = (sum((t[i] - tm)*(v[i] - vm) for i in idx)
             / sum((t[i] - tm)**2 for i in idx))
    values = [v[i] for i in idx]
    return mean, min(values), max(values), slope, len(idx)


def range_ratio(t, v):
    """Range of v over the last second divided by that over the second
    before"""
    end = t[-1]
    last = [x for ti, x in zip(t, v) if ti >= end - 1.0 - 1e-12]
    prev = [x for ti, x in zip(t, v) if end - 2.0 - 1e-12 <= ti < end - 1.0]
    if len(last) < 2 or len(prev) < 2:
        return math.nan
    r_prev = max(prev) - min(prev)
    return (max(last) - min(last))/r_prev if r_prev > 0 else math.nan


def log_stats(logs: list):
    """Courant number, time steps, pressure iterations and solver times of
    the run and its restarts"""
    max_co = 0.0
    n_steps = 0
    exec_time = clock_time = 0.0
    p_iters = []
    corr = []
    lines = []
    for log in logs:
        last = None
        for line in log.read_text(errors="replace").splitlines():
            lines.append(line)
            if line.startswith("ExecutionTime"):
                last = line
        if last:
            m = re.search(r"ExecutionTime = ([\d.e+-]+) s\s+ClockTime = "
                          r"([\d.e+-]+) s", last)
            if m:
                exec_time += float(m.group(1))
                clock_time += float(m.group(2))
    if not lines:
        return {}
    max_co_developed = 0.0
    for line in lines:
        if line.startswith("Courant Number mean:"):
            parts = line.split()
            max_co = max(max_co, float(parts[5]))
            # After the start-up (the first 100 time steps)
            if n_steps > 100:
                max_co_developed = max(max_co_developed, float(parts[5]))
        elif line.startswith("Time = "):
            n_steps += 1
            if corr:
                p_iters.append(sum(corr))
            corr = []
        elif "Solving for p," in line:
            m = re.search(r"No Iterations (\d+)", line)
            if m:
                corr.append(int(m.group(1)))
    if corr:
        p_iters.append(sum(corr))
    return {
        "n_steps": n_steps,
        "max_co": max_co,
        "max_co_developed": max_co_developed,
        "execution_time_s": exec_time,
        "clock_time_s": clock_time,
        "p_iters_per_step": (sum(p_iters)/len(p_iters)) if p_iters
        else math.nan,
    }


def case_rows(study: str, window: float, steady_tol: float,
              cl_tol: float) -> list:
    """One row per case with results: force statistics over the final
    window, steady-state flag, Courant number, iterations and cost"""
    rows = []
    for case in sorted((WORK_DIR / study).glob("*-L*")):
        meta_path = case / "case.json"
        files = force_files(case)
        if not meta_path.exists() or not files:
            continue
        meta = json.loads(meta_path.read_text())
        t, fx, fy, fxs = read_forces(files)
        cd = [COEFF_SCALE*f for f in fx]
        cl = [COEFF_SCALE*f for f in fy]
        t0 = t[-1] - window
        scd = window_stats(t, cd, t0)
        scl = window_stats(t, cl, t0)
        logs = sorted(case.glob("log.solids4Foam*"),
                      key=lambda f: (f.name != "log.solids4Foam", f.name))
        stats = log_stats(logs)
        wall = math.nan
        wc = case / "wallclock.txt"
        if wc.exists():
            walls = re.findall(r"wall_s=(\d+)", wc.read_text())
            wall = float(sum(int(w) for w in walls)) if walls else math.nan
        n_cells = math.nan
        bm = case / "log.blockMesh"
        if bm.exists():
            m = re.search(r"nCells:\s*(\d+)", bm.read_text())
            n_cells = int(m.group(1)) if m else math.nan
        cd_mean = scd[0] if scd else math.nan
        cd_rel_range = ((scd[2] - scd[1])/abs(cd_mean)) if scd else math.nan
        cd_rel_drift = (abs(scd[3])*10.0/abs(cd_mean)) if scd else math.nan
        cl_range = (scl[2] - scl[1]) if scl else math.nan
        complete = t[-1] >= meta["end_time"] - 1e-9
        cd_steady = bool(
            complete and scd is not None
            and cd_rel_range <= steady_tol
            and cd_rel_drift <= steady_tol
        )
        cl_steady = bool(complete and scl is not None and cl_range <= cl_tol)
        steady = cd_steady and cl_steady
        # Ratio of the ranges over the last and the previous second of the
        # window: about 1 for a stationary oscillation amplitude
        cd_ratio = range_ratio(t, cd)
        cl_ratio = range_ratio(t, cl)
        if steady:
            status = "steady"
        elif (complete and cd_rel_drift <= steady_tol
              and (cd_rel_range <= steady_tol
                   or 0.9 <= cd_ratio <= 1.1)
              and (cl_range <= cl_tol or 0.9 <= cl_ratio <= 1.1)):
            # Constant mean with a persistent oscillation of constant
            # amplitude
            status = "limit cycle"
        else:
            status = "drifting" if complete else "incomplete"
        # Fraction of successive differences of Cl over the window that
        # change sign: close to 1 for a period-2 (time step to time step)
        # oscillation
        wcl = [c for ti, c in zip(t, cl) if ti >= t0 - 1e-12]
        dcl = [b - a for a, b in zip(wcl, wcl[1:])]
        cl_alternation = (
            sum(1 for a, b in zip(dcl, dcl[1:]) if a*b < 0)/(len(dcl) - 1)
            if len(dcl) > 1 else math.nan
        )
        row = {
            "study": meta["study"], "case": meta["case"],
            "method": meta["method"], "level": meta["level"],
            "n_cells": n_cells, "cells_per_D": meta["cells_per_D"],
            "h_over_D": meta["h_over_D"], "ranks": meta["ranks"],
            "decomposition": meta["decomposition"],
            "partition": meta["partition"],
            "deltaT": meta["deltaT"], "dt_strategy": "fixed",
            "end_time": t[-1], "window_t0": t0,
            "Cd_last": cd[-1], "Cd_mean": cd_mean,
            "Cd_min": scd[1] if scd else math.nan,
            "Cd_max": scd[2] if scd else math.nan,
            "Cd_rel_range": cd_rel_range,
            "Cd_slope_per_s": scd[3] if scd else math.nan,
            "Cd_rel_drift_10s": cd_rel_drift,
            "Cl_last": cl[-1], "Cl_mean": scl[0] if scl else math.nan,
            "Cl_min": scl[1] if scl else math.nan,
            "Cl_max": scl[2] if scl else math.nan,
            "Cl_range": cl_range,
            "Cd_surface_traction_last": COEFF_SCALE*fxs[-1],
            "Cd_ref": CD_REF,
            "Cd_abs_err": cd_mean - CD_REF,
            "Cd_rel_err": (cd_mean - CD_REF)/CD_REF,
            "Cl_ref": CL_REF,
            "Cl_abs_err": (scl[0] - CL_REF) if scl else math.nan,
            "Cd_in_bounds": CD_BOUNDS[0] <= cd_mean <= CD_BOUNDS[1],
            "Cl_in_bounds": bool(scl and CL_BOUNDS[0] <= scl[0] <= CL_BOUNDS[1]),
            "complete": complete,
            "Cd_steady": cd_steady,
            "Cl_steady": cl_steady,
            "Cl_step_alternation": cl_alternation,
            "Cd_range_ratio": cd_ratio,
            "Cl_range_ratio": cl_ratio,
            "steady": steady,
            "status": status,
            "wall_s": wall,
            "cpu_hours": wall*meta["ranks"]/3600.0,
            "couplingCoeff": meta["fvOptions"].get("couplingCoeff", ""),
            "nOuterCorrectors": pimple_entry(case, "nOuterCorrectors"),
            "nCorrectors": pimple_entry(case, "nCorrectors"),
            "git_sha": meta["git_sha"],
        }
        row.update(stats)
        rows.append(row)
    return rows


def cmd_analyse(args) -> None:
    rows = case_rows(args.study, args.window, args.steady_tol, args.cl_tol)
    if not rows:
        sys.exit(f"No cases with results in {WORK_DIR / args.study}")
    out_dir = OUTPUT_DIR / args.study
    out_dir.mkdir(parents=True, exist_ok=True)
    out = out_dir / "runs.csv"
    with out.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    for r in rows:
        print(f"{r['case']:16s} t={r['end_time']:<6g} Cd={r['Cd_mean']:.6f} "
              f"relRange={r['Cd_rel_range']:.1e} "
              f"drift={r['Cd_rel_drift_10s']:.1e} Cl={r['Cl_mean']:.6f} "
              f"ClRange={r['Cl_range']:.1e} maxCo={r.get('max_co', math.nan):.3f} "
              f"steady={r['steady']}")
    print(f"Wrote {out}")

    if args.report:
        sys.path.insert(0, str(Path(__file__).resolve().parent))
        import convergence_report
        convergence_report.report(rows, out_dir, CD_REF, CD_BOUNDS,
                                  CL_REF, CL_BOUNDS, METHODS, LEVELS)


def cmd_extend(args) -> None:
    """Restart cases from their latest time with the same time step and
    ranks, adding args.add seconds to the end time"""
    cases = []
    for method, level in selected_cases(args):
        case = WORK_DIR / args.study / case_name(method, level, args.tag)
        meta = json.loads((case / "case.json").read_text())
        control = case / "system/controlDict"
        end = float(re.search(r"^endTime\s+([^;]+);", control.read_text(),
                              re.MULTILINE).group(1))
        times_dir = case / "processor0" if meta["ranks"] > 1 else case
        times = [float(d.name) for d in times_dir.iterdir()
                 if re.fullmatch(r"[0-9.e+-]+", d.name)]
        latest = max(times)
        if abs(latest - end) > 1e-9:
            sys.exit(f"ERROR: {case.name}: latest time {latest} is not the "
                     f"end time {end} (incomplete run?)")
        cases.append((case, meta, end))

    in_queue = queue_cpus()
    requested = sum(m["ranks"] for _, m, _ in cases)
    for case, meta, end in cases:
        print(f"{meta['case']:16s} ranks={meta['ranks']} restart t={end:g} "
              f"-> {end + args.add:g}")
    print(f"jobs={len(cases)} cpus_requested={requested} "
          f"cpus_already_in_queue={in_queue} limit={MAX_CPUS}")
    if in_queue + requested > MAX_CPUS:
        sys.exit("ERROR: the aggregate CPU limit would be exceeded")
    if not args.yes:
        print("Not extended: re-run with --yes to restart and submit")
        return

    for case, meta, end in cases:
        new_end = end + args.add
        control = case / "system/controlDict"
        replace_entry(control, "startFrom", "latestTime")
        replace_entry(control, "endTime", fmt(new_end))
        k = len(meta.get("extensions", [])) + 1
        ranks = meta["ranks"]
        script = slurm_header(
            case, f"scv-{meta['case']}-ext{k}", ranks, meta["memory"],
            meta["walltime"], meta["partition"], meta["env_sh"]
        ) + f"""\
echo "job=${{SLURM_JOB_ID:-}} host=$(hostname) ranks={ranks} restart from t={end:g} to t={new_end:g} start=$(date -Is)"
START=$(date +%s)

{mpi_run(ranks)} > log.solids4Foam.ext{k} 2>&1

END=$(date +%s)
echo "ext{k} end=$(date -Is) wall_s=$((END - START))" | tee -a wallclock.txt
"""
        path = case / f"run-ext{k}.slurm"
        path.write_text(script)
        path.chmod(0o755)
        meta.setdefault("extensions", []).append(
            {"from": end, "to": new_end, "script": path.name,
             "date": datetime.now().isoformat(timespec="seconds")})
        meta["end_time"] = new_end
        (case / "case.json").write_text(json.dumps(meta, indent=2) + "\n")
        out = subprocess.run(["sbatch", path.name], cwd=case,
                             capture_output=True, text=True, check=True)
        print(f"{case.name}: {out.stdout.strip()} ({path.name})")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    sub = parser.add_subparsers(dest="command", required=True)

    def common(p):
        p.add_argument("--study", default="mesh",
                       help="sub-directory of work/ (default mesh)")
        p.add_argument("--methods", nargs="+", default=list(METHODS),
                       choices=METHODS)
        p.add_argument("--levels", nargs="+", type=int, default=list(LEVELS),
                       choices=LEVELS)
        p.add_argument("--tag", default="", help="suffix of the case names")

    g = sub.add_parser("generate", help="generate the cases")
    common(g)
    g.add_argument("--deltaT", type=float, default=DELTA_T)
    g.add_argument("--end-time", type=float, default=END_TIME)
    g.add_argument("--write-interval", type=float, default=WRITE_INTERVAL)
    g.add_argument("--coupling-coeff", type=float, default=0.0,
                   help="override couplingCoeff of the incremental method")
    g.add_argument("--n-outer", type=int, default=0,
                   help="nOuterCorrectors of PIMPLE (default: the tutorial)")
    g.add_argument("--ranks", type=int, default=0,
                   help="override the rank policy for all the cases")
    g.add_argument("--memory", default="")
    g.add_argument("--walltime", default="")
    g.add_argument("--partition", default=PARTITION)
    g.add_argument("--env-sh", default=ENV_SH)
    g.add_argument("--force", action="store_true",
                   help="regenerate existing cases")
    g.set_defaults(func=cmd_generate)

    s = sub.add_parser("submit", help="submit the generated cases")
    common(s)
    s.add_argument("--yes", action="store_true",
                   help="submit (otherwise only the summary is printed)")
    s.set_defaults(func=cmd_submit)

    a = sub.add_parser("analyse", help="analyse the completed cases")
    a.add_argument("--study", default="mesh")
    a.add_argument("--window", type=float, default=2.0,
                   help="final time window for the statistics (s)")
    a.add_argument("--steady-tol", type=float, default=1e-5,
                   help="relative range and 10 s drift tolerance of Cd")
    a.add_argument("--cl-tol", type=float, default=1e-5,
                   help="range tolerance of Cl over the window")
    a.add_argument("--report", action="store_true",
                   help="also write convergence.csv, summary.md and plots")
    a.set_defaults(func=cmd_analyse)

    e = sub.add_parser("extend", help="restart cases and extend them")
    common(e)
    e.add_argument("--add", type=float, default=5.0,
                   help="time added to the end time (s)")
    e.add_argument("--yes", action="store_true",
                   help="restart and submit (otherwise only the summary)")
    e.set_defaults(func=cmd_extend)

    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
