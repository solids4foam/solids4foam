#!/usr/bin/env python3
"""
Robin coupling study driver for the elasticWallPressure boundary condition.

Copies an FSI tutorial, applies boundary-condition, fsiProperties, controlDict
and solid-property edits, runs it, and summarises the FSI iteration counts
from postProcessing/fsiResiduals.dat.

Examples
--------
Run a sweep of hsScale for the thickness-limited model on 3dTube, four runs
at a time:

    robin_hs_study.py run --case 3dTube --jobs 4 \
        --variant "tl_s050:bc.hsModel=thicknessLimited,bc.hsScale=0.5" \
        --variant "tl_s100:bc.hsModel=thicknessLimited"

Summarise all runs of a case:

    robin_hs_study.py summarize --case 3dTube

Variant settings are comma-separated key=value pairs:
    bc.<key>=<value>     entry in the Robin pressure patch (bc.constantHs=none
                         removes the entry)
    fsi.<key>=<value>    entry in the coupling coefficients of fsiProperties
                         (fsi.fluidSolidInterface=<type> sets the method)
    ctl.<key>=<value>    controlDict entry (endTime, deltaT, writeInterval)
    pimple.<key>=<value> entry in the fluid PIMPLE dictionary
    solid.E=<value>      Young's modulus in the solid mechanical properties
    args=<a b>           arguments passed to Allrun
"""

import argparse
import concurrent.futures
import csv
import json
import os
import re
import shutil
import signal
import subprocess
import sys
import time
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parents[3]
TUTORIALS = REPO_ROOT / "tutorials" / "fluidSolidInteraction"
DEFAULT_WORK = REPO_ROOT / "tutorialsTest-robinHs"
OPENFOAM_VERSION = os.environ.get("S4F_STUDY_OPENFOAM", "v2412")

# Case definitions: the files holding the Robin pressure BC, the FSI
# coupling settings and the control settings, relative to the case
CASES = {
    "3dTube": dict(
        source="3dTube",
        args=[],
        p_file="0/fluid.robin/p",
        patch="wall",
        fsi_file="constant/fsiProperties.robin",
        control_file="system/controlDict",
        mech_file="constant/solid/mechanicalProperties",
        end_time="0.0015",
        outputs=[
            "postProcessing/0/solidPointDisplacement_displacement.dat",
            "postProcessing/fluid/forces/0/force.dat",
        ],
    ),
    "fillingElasticContainer": dict(
        source="fillingElasticContainer",
        args=["iqnils"],
        p_file="0/fluid/p_rgh",
        patch="interface",
        fsi_file="constant/fsiProperties.iqnils",
        control_file="system/controlDict",
        mech_file="constant/solid/mechanicalProperties",
        end_time="0.5",
        outputs=["postProcessing/0/solidPointDisplacement_disp.dat"],
    ),
    "beamInCrossFlow": dict(
        source="beamInCrossFlow",
        args=["robin"],
        p_file="0/fluid/p.robin",
        patch="interface",
        fsi_file="constant/fsiProperties.robin",
        control_file="system/controlDict.iqnils",
        mech_file="constant/solid/mechanicalProperties",
        end_time="1",
        outputs=[
            "postProcessing/0/solidPointDisplacement_displacement.dat",
            "postProcessing/fluid/forces/0/force.dat",
        ],
    ),
    # Original (stiff) beam form of the verification study
    "beamInCrossFlowOriginal": dict(
        source="beamInCrossFlow",
        args=["robin"],
        p_file="0/fluid/p.robin",
        patch="interface",
        fsi_file="constant/fsiProperties.robin",
        control_file="system/controlDict.iqnils",
        mech_file="constant/solid/mechanicalProperties",
        end_time="1",
        edits=[
            ("0/fluid/U.robin", "maxVelocity", "0.2"),
            ("0/fluid/U.robin", "timeVaryingEndTime", "4.0"),
            ("0/fluid/U.robin", "timeAtMaxVelocity", "4.0"),
            ("constant/solid/mechanicalProperties", "E",
             "E [1 -1 -2 0 0 0 0] 1.4e6"),
        ],
        outputs=[
            "postProcessing/0/solidPointDisplacement_displacement.dat",
            "postProcessing/fluid/forces/0/force.dat",
        ],
    ),
    # Robin variant of the Turek-Hron FSI3 benchmark: flap wetted on both
    # sides, density ratio 1; coupling starts at t = 2 once the flow has
    # developed, as in the tutorial (an impulsive coupled start fails)
    "HronTurekFsi3Robin": dict(
        source="HronTurekFsi3",
        args=[],
        p_file="0/fluid/p",
        patch="plate",
        fsi_file="constant/fsiProperties",
        control_file="system/controlDict",
        mech_file="constant/solid/mechanicalProperties",
        end_time="2.3",
        patch_edits=[
            ("0/fluid/U", "plate", {"type": "elasticWallVelocity"}),
            ("0/fluid/p", "plate", {"type": "elasticWallPressure"}),
        ],
        fsi_edits={
            "fluidSolidInterface": "fixedRelaxation",
            "relaxationFactor": "1",
            "nOuterCorr": "100",
        },
        outputs=["postProcessing/0/solidPointDisplacement_pointDisp.dat"],
    ),
    # Robin variant of the flexible dam break: interFluid (p_rgh), obstacle
    # wetted on both sides, adjustable time step (maxCo 0.5)
    "flexibleDamBreakRobin": dict(
        source="flexibleDamBreak",
        args=[],
        p_file="0/fluid/p_rgh",
        patch="interface",
        fsi_file="constant/fsiProperties",
        control_file="system/controlDict",
        mech_file="constant/solid/mechanicalProperties",
        end_time="0.3",
        subs=[
            ("constant/fsiProperties", r"^AitkenCoeffs", "fixedRelaxationCoeffs"),
        ],
        patch_edits=[
            ("0/fluid/U", "interface", {"type": "elasticWallVelocity"}),
            ("0/fluid/p_rgh", "interface", {"type": "elasticWallPressure"}),
        ],
        fsi_edits={
            "fluidSolidInterface": "fixedRelaxation",
            "relaxationFactor": "1",
            "nOuterCorr": "100",
            "writeResidualsToFile": "yes",
        },
        outputs=["postProcessing/0/solidPointDisplacement_displacement.dat"],
    ),
    "cerebralAneurysm": dict(
        source="cerebralAneurysm",
        args=[],
        p_file="0/fluid/p",
        patch="wall",
        fsi_file="constant/fsiProperties",
        control_file="system/controlDict",
        mech_file="constant/solid/mechanicalProperties",
        end_time="0.0005",
        outputs=[
            "postProcessing/0/solidDisplacementsinnerWall.dat",
            "postProcessing/0/solidForcesinnerWall.dat",
        ],
    ),
}

EXCLUDE_NAMES = {
    "regressionTests", "verification", "images", "postProcessing",
    "__pycache__",
}


# ---------------------------------------------------------------------------
# Case preparation
# ---------------------------------------------------------------------------

def is_time_dir(path):
    try:
        float(path.name)
    except ValueError:
        return False
    return path.name != "0"


def copy_case(src, dst):
    if dst.exists():
        shutil.rmtree(dst)

    def ignore(dirname, names):
        ignored = []
        d = Path(dirname)
        for n in names:
            p = d / n
            if n in EXCLUDE_NAMES or n.startswith("processor"):
                ignored.append(n)
            elif n.startswith("log.") or n.endswith((".pdf", ".png")):
                ignored.append(n)
            elif d == src and p.is_dir() and is_time_dir(p) and n != "0":
                ignored.append(n)
            elif n == "polyMesh" and (p / "owner").exists():
                ignored.append(n)
        return ignored

    shutil.copytree(src, dst, symlinks=True, ignore=ignore)


def read_text(path):
    return Path(path).read_text()


def write_text(path, text):
    p = Path(path)
    if p.is_symlink():
        p.unlink()
    p.write_text(text)


def set_top_level_entry(text, key, value):
    pattern = re.compile(
        rf"^(\s*){re.escape(key)}(\s+)[^;]*;", re.MULTILINE
    )
    if pattern.search(text):
        return pattern.sub(rf"\g<1>{key}\g<2>{value};", text, count=1)
    return text.rstrip() + f"\n{key} {value};\n"


def set_coeffs_entry(text, key, value):
    """Set an entry inside the coupling coefficients sub-dictionary."""
    if key == "fluidSolidInterface":
        pattern = re.compile(r"^fluidSolidInterface\s+[^;]*;", re.MULTILINE)
        return pattern.sub(f"fluidSolidInterface    {value};", text, count=1)

    pattern = re.compile(
        rf"^(\s+){re.escape(key)}(\s+)[^;]*;", re.MULTILINE
    )
    if pattern.search(text):
        return pattern.sub(rf"\g<1>{key}\g<2>{value};", text, count=1)

    # Insert after outerCorrTolerance
    anchor = re.compile(r"^(\s+)outerCorrTolerance\s+[^;]*;", re.MULTILINE)
    m = anchor.search(text)
    if not m:
        raise RuntimeError(f"Cannot find outerCorrTolerance to add {key}")
    indent = m.group(1)
    return (
        text[:m.end()] + f"\n{indent}{key} {value};" + text[m.end():]
    )


def set_patch_entries(text, patch, entries):
    """Set entries in the boundaryField/<patch> block of a field file."""
    m = re.search(
        rf"^\s*\"?{re.escape(patch)}\"?\s*\n\s*\{{", text, re.MULTILINE
    )
    if not m:
        raise RuntimeError(f"Cannot find patch {patch}")

    start = m.end()
    depth = 1
    i = start
    while depth > 0:
        if text[i] == "{":
            depth += 1
        elif text[i] == "}":
            depth -= 1
        i += 1
    block = text[start:i - 1]

    for key, value in entries.items():
        pattern = re.compile(
            rf"^\s*{re.escape(key)}\s+[^;]*;\s*\n", re.MULTILINE
        )
        block = pattern.sub("", block)
        if value != "none":
            block = block.rstrip() + f"\n        {key:<15} {value};\n    "

    return text[:start] + block + text[i - 1:]


def prepare_case(case, run_dir, settings):
    cfg = CASES[case]
    copy_case(TUTORIALS / cfg["source"], run_dir)

    # Case-form edits: replace the first matching "key value;" entry
    for rel, key, value in cfg.get("edits", []):
        f = run_dir / rel
        text, n = re.subn(
            rf"^(\s*){re.escape(key)}\s+[^;]*;", rf"\g<1>{key} {value};",
            read_text(f), count=1, flags=re.MULTILINE,
        )
        if n != 1:
            raise RuntimeError(f"Cannot set {key} in {f}")
        write_text(f, text)

    # Case-level text substitutions, e.g. to rename a coefficients dictionary
    for rel, pattern, repl in cfg.get("subs", []):
        f = run_dir / rel
        text, n = re.subn(pattern, repl, read_text(f), flags=re.MULTILINE)
        if n == 0:
            raise RuntimeError(f"Pattern {pattern} not found in {f}")
        write_text(f, text)

    # Case-level patch and coupling edits, e.g. to switch to Robin BCs
    for rel, patch, entries in cfg.get("patch_edits", []):
        f = run_dir / rel
        write_text(f, set_patch_entries(read_text(f), patch, entries))

    bc = {k[3:]: v for k, v in settings.items() if k.startswith("bc.")}
    fsi = dict(cfg.get("fsi_edits", {}))
    fsi.update(
        {k[4:]: v for k, v in settings.items() if k.startswith("fsi.")}
    )
    ctl = {k[4:]: v for k, v in settings.items() if k.startswith("ctl.")}
    ctl.setdefault("endTime", cfg["end_time"])

    if bc:
        f = run_dir / cfg["p_file"]
        write_text(f, set_patch_entries(read_text(f), cfg["patch"], bc))

    if fsi:
        f = run_dir / cfg["fsi_file"]
        text = read_text(f)
        for k, v in fsi.items():
            text = set_coeffs_entry(text, k, v)
        write_text(f, text)

    f = run_dir / cfg["control_file"]
    text = read_text(f)
    for k, v in ctl.items():
        text = set_top_level_entry(text, k, v)
    write_text(f, text)

    # Fluid PIMPLE dictionary entries, in every fluid fvSolution variant
    pimple = {k[7:]: v for k, v in settings.items() if k.startswith("pimple.")}
    if pimple:
        for f in sorted((run_dir / "system" / "fluid").glob("fvSolution*")):
            if f.is_symlink():
                continue
            text = read_text(f)
            m = re.search(r"^PIMPLE\s*\n\s*\{", text, re.MULTILINE)
            if not m:
                continue
            for k, v in pimple.items():
                text = re.sub(
                    rf"^(\s+){re.escape(k)}\s+[^;]*;\s*\n", "", text,
                    flags=re.MULTILINE,
                )
                text = text[:m.end()] + f"\n    {k} {v};" + text[m.end():]
            write_text(f, text)

    if "solid.E" in settings:
        f = run_dir / cfg["mech_file"]
        text = read_text(f)
        text, n = re.subn(
            r"^(\s*E\s+(?:E\s+)?(?:\[[^\]]*\]\s*)?)[^;]*;",
            lambda m: f"{m.group(1)}{settings['solid.E']};",
            text, count=1, flags=re.MULTILINE,
        )
        if n != 1:
            raise RuntimeError("Cannot set E in " + str(f))
        write_text(f, text)


# ---------------------------------------------------------------------------
# Running
# ---------------------------------------------------------------------------

def process_group_alive(pgid):
    try:
        os.killpg(pgid, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        return True
    return True


def kill_process_group(proc, grace=30):
    """Terminate the whole process group of proc (launcher and descendants),
    escalating to SIGKILL if any member survives the grace period."""
    pgid = proc.pid
    for sig in (signal.SIGTERM, signal.SIGKILL):
        if not process_group_alive(pgid):
            break
        try:
            os.killpg(pgid, sig)
        except ProcessLookupError:
            break
        deadline = time.time() + grace
        while time.time() < deadline and process_group_alive(pgid):
            time.sleep(1)
    try:
        proc.wait(timeout=5)
    except subprocess.TimeoutExpired:
        pass


def run_variant(case, tag, settings, work_dir, timeout):
    run_dir = work_dir / case / tag
    prepare_case(case, run_dir, settings)

    args = CASES[case]["args"]
    if "args" in settings:
        args = settings["args"].split()

    meta = dict(case=case, tag=tag, settings=settings, args=args)
    (run_dir / "study.json").write_text(json.dumps(meta, indent=2))

    cmd = (
        f"source ~/bin/load-openfoam {OPENFOAM_VERSION} > /dev/null 2>&1 && "
        f"./Allrun {' '.join(args)}"
    )
    t0 = time.time()
    status = "ok"
    # Run in a new session so that the whole process tree (Allrun, meshers,
    # solver, MPI) can be terminated on timeout
    with open(run_dir / "log.study", "w") as log:
        proc = subprocess.Popen(
            ["bash", "-lc", cmd], cwd=run_dir, stdout=log,
            stderr=subprocess.STDOUT, start_new_session=True,
        )
        try:
            proc.wait(timeout=timeout)
            if proc.returncode != 0:
                status = f"exit{proc.returncode}"
        except subprocess.TimeoutExpired:
            status = "timeout"
            kill_process_group(proc)

    meta["wallTime"] = time.time() - t0
    meta["status"] = status
    (run_dir / "study.json").write_text(json.dumps(meta, indent=2))
    return summarize_run(run_dir)


# ---------------------------------------------------------------------------
# Post-processing
# ---------------------------------------------------------------------------

def read_fsi_settings(run_dir, case):
    text = read_text(run_dir / CASES[case]["fsi_file"])

    def get(key, default):
        m = re.search(rf"^\s+{key}\s+([^;]+);", text, re.MULTILINE)
        return float(m.group(1)) if m else default

    tol = get("outerCorrTolerance", 1e-6)
    return dict(
        tol=tol,
        nOuterCorr=int(get("nOuterCorr", 100)),
        pTol=get("robinPressureTolerance", tol),
        fluxTol=get("robinFluxTolerance", tol),
    )


def parse_residuals(path):
    steps = {}
    order = []
    with open(path) as f:
        header = f.readline().split()
        for line in f:
            parts = line.split()
            if len(parts) < 3:
                continue
            t = parts[0]
            vals = [float(x) for x in parts[1:]]
            if t not in steps:
                steps[t] = []
                order.append(t)
            steps[t].append(vals)
    return header, order, steps


def first_iter(rows, test):
    for row in rows:
        if test(row):
            return int(row[0])
    return None


def solver_time(run_dir):
    log = run_dir / "log.solids4Foam"
    if not log.exists():
        return None
    t = None
    with open(log, errors="replace") as f:
        for line in f:
            m = re.match(r"ExecutionTime = ([0-9.eE+-]+) s", line)
            if m:
                t = float(m.group(1))
    return t


def last_row(path):
    if not path.exists():
        return None
    last = None
    with open(path) as f:
        for line in f:
            if line.strip() and not line.startswith("#"):
                last = line.split()
    return last


def summarize_run(run_dir):
    run_dir = Path(run_dir)
    meta = json.loads(read_text(run_dir / "study.json"))
    case = meta["case"]
    out = dict(case=case, tag=meta["tag"], status=meta.get("status", "?"))
    out["settings"] = ",".join(f"{k}={v}" for k, v in meta["settings"].items())

    # Allrun does not always propagate solver failures
    log = run_dir / "log.solids4Foam"
    if log.exists():
        tail = log.read_bytes()[-4000:].decode(errors="replace")
        if "FOAM FATAL" in tail or "PETSC ERROR" in tail or "ERROR" in tail:
            out["status"] = "crashed"
        elif not re.search(r"^End\s*$", tail, re.MULTILINE):
            if out["status"] == "ok":
                out["status"] = "incomplete"

    res_file = run_dir / "postProcessing" / "fsiResiduals.dat"
    if not res_file.exists():
        out["status"] += ":noResiduals"
        return out

    fsi = read_fsi_settings(run_dir, case)
    header, order, steps = parse_residuals(res_file)
    has_robin = len(header) >= 5

    if not order:
        out["status"] += ":noResidualData"
        return out

    iters = [len(steps[t]) for t in order]
    disp_iters = []
    dp_iters = []
    final_flux = []
    rates = []
    iters_1e3 = []
    for t in order:
        rows = steps[t]

        # Tolerance-independent measures: the median contraction factor of
        # the early iterations (Robin pressure residual if present, which
        # has no inner-solver floor at this level), and the iterations to
        # reduce both residuals below 1e-3
        col = 2 if has_robin else 1
        r = [row[col] for row in rows[1:8]]
        ratios = sorted(
            r[i + 1]/r[i] for i in range(len(r) - 1)
            if r[i] > 1e-12 and r[i + 1] > 1e-12
        )
        if ratios:
            rates.append(ratios[len(ratios)//2])
        iters_1e3.append(
            first_iter(
                rows,
                lambda row: row[0] > 1 and row[1] <= 1e-3
                and (not has_robin or row[2] <= 1e-3),
            ) or len(rows)
        )

        # Iteration at which the displacement criterion is met
        disp_iters.append(
            first_iter(rows, lambda r: r[1] <= fsi["tol"] and r[0] > 1)
            or len(rows)
        )
        if has_robin:
            dp_iters.append(
                first_iter(
                    rows,
                    lambda r: r[1] <= fsi["tol"] and r[2] <= fsi["pTol"]
                    and r[0] > 1,
                ) or len(rows)
            )
            final_flux.append(rows[-1][3])

    n = len(order)
    out.update(
        steps=n,
        endTime=order[-1] if order else "",
        meanIter=sum(iters)/n,
        maxIter=max(iters),
        nMaxed=sum(1 for i in iters if i >= fsi["nOuterCorr"]),
        meanDispIter=sum(disp_iters)/n,
        maxDispIter=max(disp_iters),
        meanIter1e3=sum(iters_1e3)/n,
        medianRate=sorted(rates)[len(rates)//2] if rates else "",
    )
    if has_robin:
        out.update(
            meanDispPIter=sum(dp_iters)/n,
            maxDispPIter=max(dp_iters),
            meanFinalFlux=sum(final_flux)/n,
        )
    if len(header) >= 7:
        # Final Robin convergence state of each step (iterationError
        # criterion): leakage in column 6, state in column 7
        states = [int(steps[t][-1][5]) for t in order]
        out.update(
            meanLeakage=sum(steps[t][-1][4] for t in order)/n,
            nStalled=sum(1 for s in states if s == 2),
            nUnconverged=sum(1 for s in states if s == 0),
        )

    out["solverTime"] = solver_time(run_dir)

    for rel in CASES[case]["outputs"]:
        row = last_row(run_dir / rel)
        if row:
            out[Path(rel).stem] = " ".join(row[:4])

    diag = list((run_dir / "postProcessing").glob("robinCoefficient_*.dat"))
    if diag:
        rows = []
        with open(diag[0]) as f:
            for line in f:
                if not line.startswith("#"):
                    rows.append([float(x) for x in line.split()[1:]])
        if rows:
            out["alphaMeanFinal"] = rows[-1][1]
            ss = [r[7] for r in rows if r[7] > 0]
            sf = [r[8] for r in rows if r[8] > 0]
            if ss:
                ss.sort()
                out["solidImpedanceMedian"] = ss[len(ss)//2]
            if sf:
                sf.sort()
                out["fluidImpedanceMedian"] = sf[len(sf)//2]

    return out


def read_series(path):
    """Read a whitespace-separated time series, ignoring comments."""
    data = {}
    if not path.exists():
        return data
    with open(path) as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.replace("(", " ").replace(")", " ").split()
            try:
                vals = [float(x) for x in parts]
            except ValueError:
                continue
            data[vals[0]] = vals[1:]
    return data


def compare_runs(case_dir, case, ref_tag):
    """Maximum relative difference of each output series against ref_tag."""
    rows = []
    ref_dir = case_dir / ref_tag
    for d in sorted(case_dir.iterdir()):
        if not (d / "study.json").exists():
            continue
        row = dict(tag=d.name)
        for rel in CASES[case]["outputs"]:
            ref = read_series(ref_dir / rel)
            run = read_series(d / rel)
            times = [t for t in run if t in ref]
            if not times:
                row[Path(rel).stem] = "n/a"
                continue
            # Differences relative to the largest reference value of the
            # file, so near-zero components do not dominate
            ncol = min(len(ref[times[0]]), len(run[times[0]]))
            scale = max(
                abs(ref[t][c]) for t in ref for c in range(ncol)
            ) or 1.0
            worst = 0.0
            for c in range(ncol):
                diff = max(abs(run[t][c] - ref[t][c]) for t in times)
                worst = max(worst, diff/scale)
            row[Path(rel).stem] = f"{worst:.3g} ({len(times)} times)"
        rows.append(row)
    return rows


def write_summary(rows, path):
    keys = []
    for r in rows:
        for k in r:
            if k not in keys:
                keys.append(k)
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        for r in rows:
            w.writerow(r)


def print_table(rows):
    cols = [
        "tag", "status", "steps", "meanIter", "maxIter", "nMaxed",
        "meanDispIter", "meanDispPIter", "meanIter1e3", "medianRate",
        "meanFinalFlux", "meanLeakage", "nStalled", "nUnconverged",
        "solverTime",
        "alphaMeanFinal", "solidImpedanceMedian", "fluidImpedanceMedian",
    ]
    cols = [c for c in cols if any(c in r for r in rows)]
    print(" | ".join(cols))
    for r in rows:
        vals = []
        for c in cols:
            v = r.get(c, "")
            vals.append(f"{v:.4g}" if isinstance(v, float) else str(v))
        print(" | ".join(vals))


# ---------------------------------------------------------------------------
# Command line
# ---------------------------------------------------------------------------

def parse_variant(text):
    tag, _, body = text.partition(":")
    settings = {}
    if body:
        for item in body.split(","):
            k, _, v = item.partition("=")
            settings[k.strip()] = v.strip()
    return tag, settings


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawTextHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    r = sub.add_parser("run")
    r.add_argument("--case", required=True, choices=sorted(CASES))
    r.add_argument("--variant", action="append", required=True)
    r.add_argument("--common", default="",
                   help="settings applied to every variant")
    r.add_argument("--jobs", type=int, default=1)
    r.add_argument("--timeout", type=float, default=3600)
    r.add_argument("--work", default=str(DEFAULT_WORK))

    s = sub.add_parser("summarize")
    s.add_argument("--case", required=True, choices=sorted(CASES))
    s.add_argument("--work", default=str(DEFAULT_WORK))
    s.add_argument("--filter", default="")

    rp = sub.add_parser("report")
    rp.add_argument("--work", default=str(DEFAULT_WORK))
    rp.add_argument("--output", default="")

    c = sub.add_parser("compare")
    c.add_argument("--case", required=True, choices=sorted(CASES))
    c.add_argument("--ref", required=True, help="reference run tag")
    c.add_argument("--work", default=str(DEFAULT_WORK))

    args = ap.parse_args()
    work = Path(args.work)

    if args.cmd == "run":
        _, common = parse_variant("x:" + args.common) if args.common else ("", {})
        variants = []
        for v in args.variant:
            tag, settings = parse_variant(v)
            merged = dict(common)
            merged.update(settings)
            variants.append((tag, merged))

        rows = []
        with concurrent.futures.ThreadPoolExecutor(args.jobs) as ex:
            futs = {
                ex.submit(run_variant, args.case, tag, s, work, args.timeout):
                tag for tag, s in variants
            }
            for fut in concurrent.futures.as_completed(futs):
                try:
                    row = fut.result()
                except Exception as e:  # report and continue
                    row = dict(tag=futs[fut], status=f"error:{e}")
                rows.append(row)
                print(f"finished {row.get('tag')}: {row.get('status')} "
                      f"meanIter={row.get('meanIter')}", flush=True)
        rows.sort(key=lambda r: r.get("tag", ""))
        print_table(rows)

    elif args.cmd == "summarize":
        rows = []
        case_dir = work / args.case
        for d in sorted(case_dir.iterdir()) if case_dir.exists() else []:
            if (d / "study.json").exists() and args.filter in d.name:
                rows.append(summarize_run(d))
        write_summary(rows, case_dir / "summary.csv")
        print_table(rows)
        print(f"\nWritten {case_dir / 'summary.csv'}")

    elif args.cmd == "report":
        lines = []
        for case in CASES:
            case_dir = work / case
            if not case_dir.exists():
                continue
            rows = [
                summarize_run(d) for d in sorted(case_dir.iterdir())
                if (d / "study.json").exists()
            ]
            if not rows:
                continue
            lines.append(f"\n### {case}\n")
            lines.append(
                "| run | settings | status | steps | mean iter | max iter "
                "| iter to 1e-3 | median rate |"
            )
            lines.append("|---|---|---|---|---|---|---|---|")
            for r in rows:
                def fmt(k):
                    v = r.get(k, "")
                    return f"{v:.3g}" if isinstance(v, float) else str(v)
                lines.append(
                    f"| {r['tag']} | {r.get('settings', '')} | {r['status']} "
                    f"| {fmt('steps')} | {fmt('meanIter')} | {fmt('maxIter')} "
                    f"| {fmt('meanIter1e3')} | {fmt('medianRate')} |"
                )
        text = "\n".join(lines) + "\n"
        if args.output:
            Path(args.output).write_text(text)
            print(f"Written {args.output}")
        else:
            print(text)

    elif args.cmd == "compare":
        rows = compare_runs(work / args.case, args.case, args.ref)
        keys = list(rows[0].keys()) if rows else []
        print(" | ".join(keys))
        for r in rows:
            print(" | ".join(str(r.get(k, "")) for k in keys))


if __name__ == "__main__":
    sys.exit(main())
