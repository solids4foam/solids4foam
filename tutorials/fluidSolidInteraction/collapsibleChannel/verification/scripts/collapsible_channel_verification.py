#!/usr/bin/env python3
"""Run the collapsibleChannel verification studies.

Each case is a copy of the tutorial under verification/work/ with refined
fluid and solid meshes, a different time step or a different solid
discretisation. The history of the vertical wall displacement at 25%, 50% and
75% of the elastic segment is compared with a converged oomph-lib reference
(reference/collapsibleChannel_oomph_reference.csv).
"""

from __future__ import annotations

import argparse
import bisect
import csv
import json
import math
import re
import shutil
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent
VERIFY_DIR = SCRIPT_DIR.parent
CASE_DIR = VERIFY_DIR.parent
WORK_DIR = VERIFY_DIR / "work"
POST_DIR = VERIFY_DIR / "postProcessing"
REFERENCE_DIR = VERIFY_DIR / "reference"
REFERENCE_JSON = REFERENCE_DIR / "collapsibleChannel_verification_references.json"
REFERENCE_CSV = REFERENCE_DIR / "collapsibleChannel_oomph_reference.csv"


def reference_csv(study_spec: dict, delta_t: float) -> Path:
    """The converged reference, or the oomph-lib solution at the same time
    step, which isolates the spatial error of a case"""
    if study_spec.get("reference", "converged") == "converged":
        return REFERENCE_CSV
    return REFERENCE_DIR / f"collapsibleChannel_oomph_dt{1.0/delta_t:g}.csv"

MONITORS = ("wallQuarter", "wallMid", "wallThreeQuarter")

# Tutorial block divisions: fluid (upstream, elastic, downstream) x height and
# solid length x thickness. Fluid level n multiplies the fluid divisions by n.
FLUID_BASE = (40, 80, 80, 16)
SOLID_BASE = (80, 8)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run the collapsibleChannel verification studies"
    )
    parser.add_argument(
        "--study",
        default="all",
        help="comma-separated studies: mesh, time, solid or all (default)",
    )
    parser.add_argument(
        "--quick",
        action="store_true",
        help="smoke test: the two coarsest meshes to t = 0.5",
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
        help="number of cases run at the same time, each in serial (default 1)",
    )
    parser.add_argument(
        "--end-time",
        type=float,
        help="override the end time of every case",
    )
    parser.add_argument(
        "--cases",
        help="comma-separated case names to run, from the selected studies",
    )
    parser.add_argument(
        "--keep-going",
        action="store_true",
        help="continue after an individual case fails",
    )
    parser.add_argument(
        "--list",
        action="store_true",
        help="list the cases of the selected studies and exit",
    )
    return parser.parse_args()


def load_references() -> dict:
    return json.loads(REFERENCE_JSON.read_text())


def case_name(fluid: int, solid: tuple[int, int], dt: float, solid_type: str) -> str:
    return (
        f"F{fluid}_S{solid[0]}x{solid[1]}_dt{1.0/dt:g}_{solid_type}"
    )


def make_case(fluid: int, solid: tuple[int, int], dt: float, solid_type: str) -> dict:
    return {
        "name": case_name(fluid, solid, dt, solid_type),
        "fluid": fluid,
        "solid": solid,
        "dt": dt,
        "solidType": solid_type,
    }


def study_cases(study: str, refs: dict, quick: bool) -> list[dict]:
    spec = refs["studies"][study]
    cases = []
    if study == "mesh":
        # Fluid mesh study: the solid mesh and discretisation are fixed
        dt = spec["deltaT"]
        levels = spec["quickLevels"] if quick else spec["levels"]
        for level in levels:
            cases.append(
                make_case(level, tuple(spec["solidMesh"]), dt, spec["solidType"])
            )
    elif study == "time":
        level = spec["level"]
        steps = spec["quickDeltaT"] if quick else spec["deltaT"]
        for dt in steps:
            cases.append(
                make_case(level, tuple(spec["solidMesh"]), dt, spec["solidType"])
            )
    elif study == "static":
        meshes = spec["quickSolidMeshes"] if quick else spec["solidMeshes"]
        for solid_type in spec["solidTypes"]:
            for solid in meshes:
                cases.append({
                    "name": f"static_S{solid[0]}x{solid[1]}_{solid_type}",
                    "solid": tuple(solid),
                    "solidType": solid_type,
                    "static": True,
                })
    elif study == "solid":
        level = spec["fluidLevel"]
        meshes = spec["quickSolidMeshes"] if quick else spec["solidMeshes"]
        for solid_type in spec["solidTypes"]:
            for solid in meshes:
                cases.append(
                    make_case(level, tuple(solid), spec["deltaT"], solid_type)
                )
    else:
        raise RuntimeError(f"unknown study {study}")

    # Optional fsiProperties entries that differ from the tutorial; they are
    # part of the case name, so such cases are never shared with other studies
    overrides = spec.get("fsiProperties", {})
    if overrides:
        suffix = "_".join(f"{key}{value}" for key, value in overrides.items())
        for case in cases:
            case["fsi"] = overrides
            case["name"] += f"_{suffix}"
    return cases


def ignored(directory: str, names: list[str]) -> set[str]:
    """Copy only the inputs a fresh run needs, never a previous result."""
    ignored_names = {
        "verification", "regressionTests", "postProcessing", "images",
        "dynamicCode",
    }
    ignored_names.update(n for n in names if n.startswith("processor"))
    ignored_names.update(n for n in names if n.startswith("log."))
    directory_path = Path(directory)
    if directory_path == CASE_DIR:
        ignored_names.update(
            n for n in names
            if n != "0" and re.fullmatch(r"[0-9]+(?:\.[0-9]+)?(?:e-?[0-9]+)?", n)
        )
    if directory_path.name in {"fluid", "solid"} and directory_path.parent.name == "constant":
        ignored_names.add("polyMesh")
    return ignored_names.intersection(names)


def replace_once(path: Path, pattern: str, replacement: str) -> None:
    text, count = re.subn(pattern, replacement, path.read_text(), count=1,
                          flags=re.MULTILINE)
    if count != 1:
        raise RuntimeError(f"pattern {pattern!r} not found in {path}")
    path.write_text(text)


def set_block_divisions(path: Path, divisions: list[tuple[int, int]]) -> None:
    """Set the (nx ny 1) divisions of each hex block, in order."""
    pattern = re.compile(r"(hex\s*\([^)]*\)\s*)\(\s*\d+\s+\d+\s+(\d+)\s*\)")
    matches = list(pattern.finditer(path.read_text()))
    if len(matches) != len(divisions):
        raise RuntimeError(f"expected {len(divisions)} blocks in {path}")
    counter = iter(divisions)

    def scale(match: re.Match[str]) -> str:
        nx, ny = next(counter)
        return f"{match.group(1)}({nx} {ny} {match.group(2)})"

    path.write_text(pattern.sub(scale, path.read_text()))


def configure_case(run_dir: Path, case: dict, end_time: float) -> None:
    f = case["fluid"]
    set_block_divisions(
        run_dir / "system" / "fluid" / "blockMeshDict",
        [
            (FLUID_BASE[0]*f, FLUID_BASE[3]*f),
            (FLUID_BASE[1]*f, FLUID_BASE[3]*f),
            (FLUID_BASE[2]*f, FLUID_BASE[3]*f),
        ],
    )
    set_block_divisions(
        run_dir / "system" / "solid" / "blockMeshDict", [case["solid"]]
    )

    control = run_dir / "system" / "controlDict"
    replace_once(control, r"^(deltaT\s+)[^;]+;", rf"\g<1>{case['dt']:.10g};")
    replace_once(control, r"^(endTime\s+)[^;]+;", rf"\g<1>{end_time:.10g};")
    replace_once(control, r"^(writeInterval\s+)[^;]+;", rf"\g<1>{end_time:.10g};")

    # The initial IQN-ILS relaxation must shrink with the square of the time
    # step: the added-mass pressure of a given interface increment grows
    # with 1/deltaT^2 while the wall stiffness does not
    fsi = run_dir / "constant" / "fsiProperties"
    omega = float(
        re.search(r"^\s*relaxationFactor\s+([^;]+);", fsi.read_text(),
                  re.MULTILINE).group(1)
    )
    base_dt = float(
        re.search(r"^deltaT\s+([^;]+);",
                  (CASE_DIR / "system" / "controlDict").read_text(),
                  re.MULTILINE).group(1)
    )
    replace_once(
        fsi,
        r"^(\s*relaxationFactor\s+)[^;]+;",
        rf"\g<1>{omega*(case['dt']/base_dt)**2:.6g};",
    )

    for key, value in case.get("fsi", {}).items():
        replace_once(fsi, rf"^(\s*{key}\s+)[^;]+;", rf"\g<1>{value};")

    # Matching interface faces can be mapped directly; otherwise use AMI
    if case["solid"][0] != FLUID_BASE[1]*f:
        replace_once(
            run_dir / "constant" / "fsiProperties",
            r"^(\s*interfaceTransferMethod\s+)\w+;",
            r"\g<1>AMI;",
        )



def completed(run_dir: Path, end_time: float) -> bool:
    log = run_dir / "log.solids4Foam"
    if not log.is_file() or not re.search(
        r"^End\s*$", log.read_text(errors="replace"), re.MULTILINE
    ):
        return False
    try:
        history = read_history(run_dir, "wallMid")
    except RuntimeError:
        return False
    return bool(history) and history[-1][0] >= end_time - 1e-9


def run_case(case: dict, end_time: float, reuse: bool) -> Path:
    run_dir = WORK_DIR / case["name"]
    if reuse and completed(run_dir, end_time):
        print(f"  reusing {case['name']}")
        return run_dir
    if run_dir.exists():
        shutil.rmtree(run_dir)
    shutil.copytree(CASE_DIR, run_dir, ignore=ignored, symlinks=True)
    configure_case(run_dir, case, end_time)

    command = ["./Allrun"]
    if case["solidType"] == "linear":
        command.append("linear")
    print(f"  running {case['name']}", flush=True)
    with (run_dir / "log.Allverify").open("w") as handle:
        result = subprocess.run(
            command, cwd=run_dir, stdout=handle, stderr=subprocess.STDOUT
        )
    log = run_dir / "log.solids4Foam"
    text = log.read_text(errors="replace") if log.is_file() else ""
    if result.returncode or not re.search(r"^End\s*$", text, re.MULTILINE):
        raise RuntimeError(f"{case['name']} failed; see {run_dir}")
    return run_dir


def run_static_case(case: dict, refs: dict, reuse: bool) -> Path:
    """Solve the wall alone under the external pressure, raised in steps."""
    run_dir = WORK_DIR / case["name"]
    log = run_dir / "log.solids4Foam"
    if (
        reuse and log.is_file()
        and re.search(r"^End\s*$", log.read_text(errors="replace"), re.MULTILINE)
    ):
        print(f"  reusing {case['name']}")
        return run_dir
    if run_dir.exists():
        shutil.rmtree(run_dir)
    steps = refs["studies"]["static"]["loadSteps"]
    for sub in ("0", "constant", "system"):
        (run_dir / sub).mkdir(parents=True)
    for name in ("dynamicMeshDict", "g", "mechanicalProperties"):
        shutil.copy(CASE_DIR / "constant" / "solid" / name, run_dir / "constant")
    shutil.copy(
        CASE_DIR / "constant" / "solid" / f"solidProperties.{case['solidType']}",
        run_dir / "constant" / "solidProperties",
    )
    for name in ("blockMeshDict", "fvSchemes"):
        shutil.copy(CASE_DIR / "system" / "solid" / name, run_dir / "system")
    shutil.copy(
        CASE_DIR / "system" / "solid" / f"fvSolution.{case['solidType']}",
        run_dir / "system" / "fvSolution",
    )
    for path in (CASE_DIR / "0" / "solid").iterdir():
        shutil.copy(path, run_dir / "0")
    set_block_divisions(run_dir / "system" / "blockMeshDict", [case["solid"]])
    (run_dir / "constant" / "physicsProperties").write_text(
        "FoamFile { version 2.0; format ascii; class dictionary;"
        " object physicsProperties; }\ntype solid;\n"
    )
    (run_dir / "constant" / "pExt").write_text(
        f"(\n    (0 0)\n    ({steps} 200)\n)\n"
    )
    replace_once(
        run_dir / "0" / "D",
        r'"\$FOAM_CASE/constant/solid/pExt"',
        '"$FOAM_CASE/constant/pExt"',
    )
    monitors = "\n".join(
        f"    {name} {{ type solidPointDisplacement; point ({x} 1 0.05); }}"
        for name, x in zip(MONITORS, (7.5, 10, 12.5))
    )
    (run_dir / "system" / "controlDict").write_text(
        "FoamFile { version 2.0; format ascii; class dictionary;"
        " object controlDict; }\n"
        "application solids4Foam;\nstartFrom startTime;\nstartTime 0;\n"
        f"stopAt endTime;\nendTime {steps};\ndeltaT 1;\n"
        f"writeControl timeStep;\nwriteInterval {steps};\n"
        "writeFormat ascii;\nwritePrecision 10;\ntimeFormat general;\n"
        "timePrecision 6;\nrunTimeModifiable no;\n"
        f"functions\n{{\n{monitors}\n}}\n"
    )
    # Run through the tutorial run functions, which also restore the library
    # path that macOS strips from child processes
    allrun = run_dir / "Allrun"
    allrun.write_text(
        "#!/bin/bash\n"
        ". $WM_PROJECT_DIR/bin/tools/RunFunctions\n"
        "source solids4FoamScripts.sh\n"
        "solids4Foam::runApplication blockMesh\n"
        "solids4Foam::runApplication solids4Foam\n"
    )
    allrun.chmod(0o755)
    print(f"  running {case['name']}", flush=True)
    with (run_dir / "log.Allverify").open("w") as handle:
        subprocess.run(
            ["./Allrun"], cwd=run_dir, stdout=handle, stderr=subprocess.STDOUT
        )
    text = log.read_text(errors="replace") if log.is_file() else ""
    if not re.search(r"^End\s*$", text, re.MULTILINE):
        raise RuntimeError(f"{case['name']} failed; see {run_dir}")
    return run_dir


def evaluate_static(case: dict, run_dir: Path, refs: dict) -> dict:
    reference = refs["studies"]["static"]["reference"]
    row = {
        "case": case["name"],
        "solidCells": f"{case['solid'][0]}x{case['solid'][1]}",
        "solidType": case["solidType"],
    }
    for monitor in MONITORS:
        value = read_history(run_dir, monitor)[-1][1]
        row[monitor] = value
        # Negative: the wall deflects less than the beam, i.e. is too stiff
        row[f"{monitor}_relError"] = (value - reference[monitor])/reference[monitor]
    text = (run_dir / "log.solids4Foam").read_text(errors="replace")
    clock = re.findall(r"ClockTime = ([0-9.eE+-]+) s", text)
    row["clockTime"] = float(clock[-1]) if clock else math.nan
    return row


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


def read_history(run_dir: Path, monitor: str) -> list[tuple[float, float]]:
    candidates = sorted(
        run_dir.glob(f"postProcessing/**/solidPointDisplacement_{monitor}.dat")
    )
    if not candidates:
        raise RuntimeError(f"no {monitor} displacement data in {run_dir}")
    # A restart writes a second file; keep the latest value at each time
    values: dict[float, float] = {}
    for path in candidates:
        for row in numeric_rows(path):
            if len(row) >= 3:
                values[round(row[0], 10)] = row[2]
    return sorted(values.items())


def read_reference(path: Path) -> dict[str, list[tuple[float, float]]]:
    with path.open() as handle:
        rows = [
            row for row in csv.DictReader(
                line for line in handle if not line.startswith("#")
            )
        ]
    return {
        monitor: [(float(row["time"]), float(row[monitor])) for row in rows]
        for monitor in MONITORS
    }


def interpolate(curve: list[tuple[float, float]], time: float) -> float:
    times = [point[0] for point in curve]
    index = bisect.bisect_left(times, time)
    if index <= 0:
        return curve[0][1]
    if index >= len(curve):
        return curve[-1][1]
    (t0, v0), (t1, v1) = curve[index - 1], curve[index]
    return v0 + (v1 - v0)*(time - t0)/(t1 - t0)


def trough(curve: list[tuple[float, float]], window: tuple[float, float]) -> tuple[float, float]:
    """Minimum of the curve in a time window, refined by a parabola."""
    inside = [
        (i, point) for i, point in enumerate(curve)
        if window[0] <= point[0] <= window[1]
    ]
    index, (time, value) = min(inside, key=lambda item: item[1][1])
    if 0 < index < len(curve) - 1:
        (ta, va), (tb, vb), (tc, vc) = curve[index - 1:index + 2]
        denominator = (ta - tb)*(ta - tc)*(tb - tc)
        if abs(denominator) > 0:
            a = (tc*(vb - va) + tb*(va - vc) + ta*(vc - vb))/denominator
            b = (tc*tc*(va - vb) + tb*tb*(vc - va) + ta*ta*(vb - vc))/denominator
            c = (
                tb*tc*(tb - tc)*va + tc*ta*(tc - ta)*vb + ta*tb*(ta - tb)*vc
            )/denominator
            if a > 0:
                time = -b/(2*a)
                value = c - b*b/(4*a)
    return time, value


def solver_statistics(run_dir: Path) -> dict:
    text = (run_dir / "log.solids4Foam").read_text(errors="replace")
    iterations: dict[str, int] = {}
    for time, iteration in re.findall(
        r"^Time = (\S+), iteration: (\d+)", text, re.MULTILINE
    ):
        iterations[time] = max(iterations.get(time, 0), int(iteration))
    execution = re.findall(r"ExecutionTime = ([0-9.eE+-]+) s", text)
    clock = re.findall(r"ClockTime = ([0-9.eE+-]+) s", text)
    counts = list(iterations.values())
    return {
        "steps": len(counts),
        "meanFsiIterations": sum(counts)/len(counts) if counts else math.nan,
        "maxFsiIterations": max(counts) if counts else 0,
        "executionTime": float(execution[-1]) if execution else math.nan,
        "clockTime": float(clock[-1]) if clock else math.nan,
    }


def cell_counts(run_dir: Path) -> tuple[int, int]:
    counts = []
    for region in ("fluid", "solid"):
        log = run_dir / f"log.blockMesh.{region}"
        match = re.search(r"nCells:\s*(\d+)", log.read_text(errors="replace"))
        counts.append(int(match.group(1)) if match else 0)
    return counts[0], counts[1]


def evaluate(case: dict, run_dir: Path, reference: dict, refs: dict) -> dict:
    row = {
        "case": case["name"],
        "fluidLevel": case["fluid"],
        "solidCells": f"{case['solid'][0]}x{case['solid'][1]}",
        "deltaT": case["dt"],
        "solidType": case["solidType"],
    }
    row["fluidCells"], row["solidCellCount"] = cell_counts(run_dir)
    window = tuple(refs["troughWindow"])
    for monitor in MONITORS:
        history = [
            point for point in read_history(run_dir, monitor) if point[0] > 0
        ]
        ref_curve = reference[monitor]
        scale = max(abs(value) for _, value in ref_curve)
        errors = [abs(value - interpolate(ref_curve, t)) for t, value in history]
        row[f"{monitor}_maxError"] = max(errors)/scale
        row[f"{monitor}_rmsError"] = math.sqrt(
            sum(e*e for e in errors)/len(errors)
        )/scale
    history = read_history(run_dir, "wallMid")
    row["troughTime"], row["troughValue"] = trough(history, window)
    row.update(solver_statistics(run_dir))
    return row


def self_convergence(results: list[tuple[dict, Path, dict]]) -> list[str]:
    """Differences between successive time steps, and the observed order.

    The comparison with the converged reference also carries the spatial
    error of the fixed mesh, so the temporal convergence itself is measured
    between the solids4foam solutions: the difference of the wall-midpoint
    histories at the times of the coarser run, relative to the peak of the
    reference.
    """
    ordered = sorted(results, key=lambda item: -item[0]["dt"])
    scale = max(abs(v) for _, v in read_reference(REFERENCE_CSV)["wallMid"])
    lines, diffs = [], []
    for (case_a, dir_a, _), (case_b, dir_b, _) in zip(ordered, ordered[1:]):
        coarse = dict(read_history(dir_a, "wallMid"))
        fine = dict(read_history(dir_b, "wallMid"))
        common = [t for t in coarse if t in fine and t > 0]
        diff = max(abs(coarse[t] - fine[t]) for t in common)/scale
        diffs.append((case_a["dt"]/case_b["dt"], diff))
        lines.append(
            f"- dt {case_a['dt']:g} vs {case_b['dt']:g}: max difference"
            f" {diff:.3%}"
        )
    for (ratio, d1), (_, d2) in zip(diffs, diffs[1:]):
        if d2 > 0:
            lines.append(
                f"- observed temporal order {math.log(d1/d2)/math.log(ratio):.2f}"
            )
    return lines


def observed_order(errors: list[float], ratio: float) -> float | None:
    if len(errors) < 2 or errors[-1] <= 0 or errors[-2] <= 0:
        return None
    return math.log(errors[-2]/errors[-1])/math.log(ratio)


def write_rows(study: str, rows: list[dict]) -> Path:
    POST_DIR.mkdir(parents=True, exist_ok=True)
    path = POST_DIR / f"{study}_study.csv"
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
    return path


def describe(row: dict) -> str:
    if "wallMid_relError" in row:
        return (
            f"{row['case']}: wall midpoint {row['wallMid']:.6f},"
            f" error {row['wallMid_relError']:+.3%}, {row['clockTime']:.0f} s"
        )
    return (
        f"{row['case']}: max error {row['wallMid_maxError']:.3%},"
        f" trough {row['troughValue']:.5f} at {row['troughTime']:.4f},"
        f" {row['meanFsiIterations']:.1f} FSI iterations/step,"
        f" {row['clockTime']:.0f} s"
    )


def summary_table(study: str, rows: list[dict]) -> list[str]:
    if study == "static":
        lines = [
            "| Case | Wall midpoint (m) | Error | Clock (s) |",
            "|---|---:|---:|---:|",
        ]
        for row in rows:
            lines.append(
                f"| {row['case']} | {row['wallMid']:.6f} |"
                f" {row['wallMid_relError']:+.3%} | {row['clockTime']:.0f} |"
            )
        return lines
    lines = [
        "| Case | Max error (mid) | RMS error (mid) | Trough | FSI it./step | Clock (s) |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for row in rows:
        lines.append(
            f"| {row['case']} | {row['wallMid_maxError']:.3%} |"
            f" {row['wallMid_rmsError']:.3%} | {row['troughValue']:.5f} |"
            f" {row['meanFsiIterations']:.1f} | {row['clockTime']:.0f} |"
        )
    return lines


def check_study(study: str, rows: list[dict], refs: dict, quick: bool) -> list[str]:
    """Return the failed acceptance checks of a study."""
    failures = []
    criteria = refs["studies"][study]["acceptance"]
    if quick:
        return failures
    if study in ("mesh", "time"):
        errors = [row["wallMid_maxError"] for row in rows]
        if any(b >= a for a, b in zip(errors, errors[1:])):
            failures.append(f"{study}: wall-midpoint error does not decrease")
        if errors[-1] > criteria["finestMaxError"]:
            failures.append(
                f"{study}: finest wall-midpoint error {errors[-1]:.3g}"
                f" > {criteria['finestMaxError']}"
            )
    if study == "static":
        # The high-order wall must be accurate on every mesh; the linear one
        # only once its in-plane spacing is fine enough
        for row in rows:
            tolerance = criteria["maxRelError"].get(row["solidType"])
            if (
                tolerance is None
                and row["solidCells"] == criteria["linearFinestMesh"]
            ):
                tolerance = criteria["linearFinestRelError"]
            if tolerance is not None and abs(row["wallMid_relError"]) > tolerance:
                failures.append(
                    f"static: {row['case']} wall-midpoint error"
                    f" {row['wallMid_relError']:.3g} > {tolerance}"
                )
    if study == "solid":
        for row in rows:
            tolerance = criteria["maxError"].get(row["solidType"])
            if tolerance is not None and row["wallMid_maxError"] > tolerance:
                failures.append(
                    f"solid: {row['case']} wall-midpoint error"
                    f" {row['wallMid_maxError']:.3g} > {tolerance}"
                )
    return failures


def write_plot(study: str, cases: list[tuple[dict, Path]], reference: Path) -> None:
    if not shutil.which("gnuplot"):
        return
    lines = [
        'set terminal pngcairo enhanced size 1000,600 font ",11"',
        f'set output "{POST_DIR / (study + "_wallMid.png")}"',
        'set xlabel "t (s)"',
        'set ylabel "vertical displacement at x = 10 m (m)"',
        'set grid',
        'set key bottom right',
        'set datafile separator ","',
        f'plot "{reference}" every ::1 using 1:3 with lines lw 3 '
        'lc rgb "black" title "oomph-lib reference"',
    ]
    for case, run_dir in cases:
        data = sorted(
            run_dir.glob("postProcessing/**/solidPointDisplacement_wallMid.dat")
        )
        if data:
            lines[-1] += (
                f', "{data[-1]}" using 1:3 with lines lw 1.5 '
                f'title "{case["name"]}" noenhanced'
            )
    lines[-1] = lines[-1].replace(', "', ',\\\n     "')
    # The solver output is whitespace separated
    script = "\n".join(lines).replace(
        'set datafile separator ","\n', ""
    )
    script = script.replace(
        f'plot "{reference}"',
        f'plot "< grep -v ^# {reference} | tr , \' \'"',
    )
    path = POST_DIR / f"{study}_wallMid.gnuplot"
    path.write_text(script + "\n")
    subprocess.run(["gnuplot", str(path)], check=False)


def main() -> int:
    args = parse_args()
    refs = load_references()
    studies = (
        ["static", "solid", "mesh", "time"] if args.study == "all"
        else [s.strip() for s in args.study.split(",")]
    )
    end_time = args.end_time or (
        refs["quickEndTime"] if args.quick else refs["endTime"]
    )
    selected = set(args.cases.split(",")) if args.cases else None

    if args.list:
        for study in studies:
            for case in study_cases(study, refs, args.quick):
                print(f"{study}: {case['name']}")
        return 0

    for tool in ("blockMesh", "solids4Foam"):
        if not shutil.which(tool):
            print(f"ERROR: {tool} not found; source OpenFOAM and solids4foam",
                  file=sys.stderr)
            return 2

    all_failures = []
    # Cases shared between studies are run once per invocation
    finished: dict[str, Path] = {}
    summary = ["# collapsibleChannel verification summary", ""]
    for study in studies:
        print(f"Study: {study}", flush=True)
        cases = [
            case for case in study_cases(study, refs, args.quick)
            if not selected or case["name"] in selected
        ]

        spec = refs["studies"][study]

        def run_one(case: dict):
            try:
                if case.get("static"):
                    run_dir = run_static_case(case, refs, args.reuse)
                    return case, run_dir, evaluate_static(case, run_dir, refs)
                run_dir = finished.get(case["name"])
                if run_dir is None:
                    run_dir = run_case(case, end_time, args.reuse)
                reference = read_reference(reference_csv(spec, case["dt"]))
                return case, run_dir, evaluate(case, run_dir, reference, refs)
            except RuntimeError as error:
                return case, None, str(error)

        with ThreadPoolExecutor(max_workers=max(args.cores, 1)) as pool:
            outcomes = list(pool.map(run_one, cases))

        results = []
        for case, run_dir, row in outcomes:
            if run_dir is None:
                print(f"  FAILED: {row}", flush=True)
                all_failures.append(row)
                continue
            finished[case["name"]] = run_dir
            results.append((case, run_dir, row))
            print(f"  {describe(row)}", flush=True)
        if all_failures and not args.keep_going:
            return 1
        if not results:
            continue
        rows = [row for _, _, row in results]
        path = write_rows(study, rows)
        if study != "static":
            write_plot(
                study,
                [(case, run_dir) for case, run_dir, _ in results],
                reference_csv(spec, results[-1][0]["dt"]),
            )
        failures = check_study(study, rows, refs, args.quick)
        all_failures.extend(failures)
        summary += [f"## {study}", "", f"Results: `{path.name}`", ""]
        summary += summary_table(study, rows)
        if study == "time" and len(results) > 1:
            lines = self_convergence(results)
            summary += [""] + lines
            for line in lines:
                print(f"  {line[2:]}", flush=True)
        summary += [""] + [f"- FAIL: {f}" for f in failures] + [""]

    POST_DIR.mkdir(parents=True, exist_ok=True)
    (POST_DIR / "verification_summary.md").write_text("\n".join(summary) + "\n")
    if all_failures:
        print("Verification FAILED:")
        for failure in all_failures:
            print(f"  {failure}")
        return 1
    print("Verification PASSED")
    return 0


if __name__ == "__main__":
    sys.exit(main())
