#!/usr/bin/env python3
"""Run the manufactured-solution FSI convergence study and report the orders.

Each level copies the case to verification/work/<coupling>/n<nCells>, sets
system/caseParameters (nCells cells along the interface, square cells and
deltaT = 0.32/nCells, i.e. combined space-time refinement) and runs Allrun.
The errors written by the mmsFsiErrors function object at t = 0.2 and 0.4 are
compared with the reference errors, and the observed orders between
successive levels are reported.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
import shutil
import subprocess
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent
VERIFY_DIR = SCRIPT_DIR.parent
CASE_DIR = VERIFY_DIR.parent
WORK_DIR = VERIFY_DIR / "work"
POST_DIR = VERIFY_DIR / "postProcessing"
REFERENCE_FILE = (
    VERIFY_DIR
    / "reference"
    / "manufacturedSolutionFSI_verification_references.json"
)

# deltaT = DT_TIMES_NCELLS/nCells: 0.02 on the 16-cell mesh
DT_TIMES_NCELLS = 0.32


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run the manufacturedSolutionFSI verification study"
    )
    parser.add_argument(
        "--couplings",
        help="comma-separated couplings: aitken, robin",
    )
    parser.add_argument(
        "--levels",
        help="comma-separated cells along the interface, multiples of 4",
    )
    parser.add_argument(
        "--quick",
        action="store_true",
        help="run the first two levels as a smoke test",
    )
    parser.add_argument(
        "--reuse",
        action="store_true",
        help="reuse completed cases under verification/work",
    )
    return parser.parse_args()


def ignored(directory: str, names: list[str]) -> set[str]:
    ignored_names = {
        "verification",
        "regressionTests",
        "postProcessing",
        "derivation",
        "src",
    }
    ignored_names.update(name for name in names if name.startswith("processor"))
    ignored_names.update(name for name in names if name.startswith("log."))

    directory_path = Path(directory)
    if directory_path == CASE_DIR:
        ignored_names.update(
            name
            for name in names
            if name != "0" and re.fullmatch(r"[0-9]+(?:\.[0-9]+)?", name)
        )
    if directory_path.name in {"fluid", "solid"}:
        ignored_names.add("polyMesh")

    return ignored_names.intersection(names)


def set_parameters(run_dir: Path, n_cells: int) -> None:
    path = run_dir / "system" / "caseParameters"
    text = path.read_text()
    values = {
        "nCells": f"{n_cells}",
        "nCellsFluidY": f"{n_cells // 2}",
        "nCellsSolidY": f"{n_cells // 4}",
        "deltaT": f"{DT_TIMES_NCELLS / n_cells:.10g}",
    }
    for key, value in values.items():
        text, count = re.subn(
            rf"(?m)^({key}\s+)\S+;",
            rf"\g<1>{value};",
            text,
        )
        if count != 1:
            raise RuntimeError(f"could not set {key} in {path}")
    path.write_text(text)


def read_errors(run_dir: Path, times: list[str]) -> dict[str, dict[str, float]]:
    path = run_dir / "postProcessing" / "mmsFsiErrors" / "errors.dat"
    lines = path.read_text().splitlines()
    names = lines[0].lstrip("#").split()
    errors = {}
    for line in lines[1:]:
        values = [float(value) for value in line.split()]
        for time in times:
            if abs(values[0] - float(time)) < 1e-8:
                errors[time] = dict(zip(names[1:], values[1:]))
    missing = [time for time in times if time not in errors]
    if missing:
        raise RuntimeError(f"no errors at t = {', '.join(missing)} in {path}")
    return errors


def run_level(
    coupling: str, n_cells: int, times: list[str], reuse: bool
) -> dict[str, dict[str, float]]:
    run_dir = WORK_DIR / coupling / f"n{n_cells}"
    solver_log = run_dir / "log.solids4Foam"
    completed_case = (
        reuse
        and solver_log.exists()
        and re.search(r"^End\s*$", solver_log.read_text(), re.MULTILINE)
    )

    if completed_case:
        print(f"Reusing {coupling} n{n_cells} in {run_dir}")
    else:
        if run_dir.exists():
            shutil.rmtree(run_dir)
        shutil.copytree(CASE_DIR, run_dir, ignore=ignored, symlinks=True)
        # Share the tutorial library, so that it is compiled only once
        (run_dir / "src").symlink_to(CASE_DIR / "src")
        set_parameters(run_dir, n_cells)
        print(f"Running {coupling} n{n_cells} in {run_dir}", flush=True)
        with (run_dir / "log.Allverify").open("w") as log:
            completed = subprocess.run(
                ["./Allrun", coupling],
                cwd=run_dir,
                stdout=log,
                stderr=subprocess.STDOUT,
                check=False,
            )
        if completed.returncode:
            raise RuntimeError(
                f"{coupling} n{n_cells} failed; "
                f"see {run_dir / 'log.Allverify'}"
            )
        if not solver_log.exists():
            raise RuntimeError(
                f"solver log was not created in {run_dir}; "
                f"see {run_dir / 'log.Allverify'}"
            )

    log_text = solver_log.read_text()
    if not re.search(r"^End\s*$", log_text, re.MULTILINE):
        raise RuntimeError(f"incomplete solver log: {solver_log}")
    if "DIVERGED_" in log_text:
        raise RuntimeError(f"PETSc failed to converge: {solver_log}")

    return read_errors(run_dir, times)


def order(coarse: float, fine: float, ratio: float) -> float:
    if coarse <= 0 or fine <= 0:
        return math.nan
    return math.log(coarse / fine) / math.log(ratio)


def write_results(
    coupling: str,
    levels: list[int],
    errors: dict[int, dict[str, dict[str, float]]],
    reference: dict,
    quick: bool,
) -> bool:
    POST_DIR.mkdir(parents=True, exist_ok=True)
    times = [f"{time:g}" for time in reference["times"]]
    metrics = reference["metrics"]
    minimum_order = reference["acceptance"]["minimum_net_order"]
    tolerance = reference["acceptance"]["relative_tolerance"]
    reference_errors = reference["couplings"][coupling]["reference_errors"]

    with (POST_DIR / f"{coupling}_errors.csv").open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["nCells", "time", *metrics])
        for n_cells in levels:
            for time in times:
                writer.writerow(
                    [n_cells, time]
                    + [f"{errors[n_cells][time][m]:.10g}" for m in metrics]
                )

    passed = True
    lines = [f"## {coupling}", ""]
    for time in times:
        lines.extend(
            [
                f"### t = {time}",
                "",
                "| quantity | "
                + " | ".join(f"n{n}" for n in levels)
                + " | orders | net order | reference |",
                "|---|" + "---|" * len(levels) + "---|---|---|",
            ]
        )
        for metric in metrics:
            values = [errors[n][time][metric] for n in levels]
            orders = [
                order(values[i], values[i + 1], levels[i + 1] / levels[i])
                for i in range(len(levels) - 1)
            ]
            net = order(values[0], values[-1], levels[-1] / levels[0])

            status = []
            if not all(math.isfinite(v) and v > 0 for v in values):
                status.append("FAIL: non-positive error")
            if not quick and not (math.isfinite(net) and net >= minimum_order):
                status.append(f"FAIL: net order < {minimum_order:g}")
            for n_cells, value in zip(levels, values):
                ref = (
                    reference_errors.get(str(n_cells), {})
                    .get(time, {})
                    .get(metric)
                )
                if ref is not None and abs(value - ref) > tolerance*ref:
                    status.append(f"FAIL: n{n_cells} differs from {ref:.4e}")
            if status:
                passed = False

            lines.append(
                f"| {metric} | "
                + " | ".join(f"{v:.3e}" for v in values)
                + " | "
                + ", ".join(f"{o:.2f}" for o in orders)
                + f" | {net:.2f} | "
                + ("; ".join(status) if status else "PASS")
                + " |"
            )
        lines.append("")

    print("\n".join(lines))
    with (POST_DIR / "verification_summary.md").open("a") as stream:
        stream.write("\n".join(lines) + "\n")
    return passed


def main() -> int:
    args = parse_args()
    reference = json.loads(REFERENCE_FILE.read_text())
    all_couplings = reference["couplings"]
    if args.couplings:
        couplings = [value.strip() for value in args.couplings.split(",")]
    else:
        couplings = list(reference["default_couplings"])
    invalid = [name for name in couplings if name not in all_couplings]
    if not couplings or invalid or len(set(couplings)) != len(couplings):
        available = ", ".join(all_couplings)
        raise SystemExit(
            f"couplings must be a unique selection from: {available}"
        )

    if args.levels:
        levels = [int(value) for value in args.levels.split(",")]
    else:
        levels = [int(value) for value in reference["default_levels"]]
    if args.quick:
        levels = levels[:2]
    if len(levels) < 2 or any(level <= 0 or level % 4 for level in levels):
        raise SystemExit(
            "at least two levels, positive multiples of 4, are required"
        )
    if any(right <= left for left, right in zip(levels, levels[1:])):
        raise SystemExit("levels must be strictly increasing")

    times = [f"{time:g}" for time in reference["times"]]

    POST_DIR.mkdir(parents=True, exist_ok=True)
    (POST_DIR / "verification_summary.md").write_text(
        "# Manufactured-solution FSI verification summary\n\n"
        f"- Mode: {'quick smoke test' if args.quick else 'full verification'}\n"
        f"- Levels: {', '.join(str(n) for n in levels)} cells along the "
        f"interface, deltaT = {DT_TIMES_NCELLS:g}/nCells\n\n"
    )

    passed = True
    for coupling in couplings:
        errors = {}
        try:
            for n_cells in levels:
                errors[n_cells] = run_level(
                    coupling, n_cells, times, args.reuse
                )
        except RuntimeError as error:
            print(f"ERROR: {error}")
            passed = False
            continue
        passed = (
            write_results(coupling, levels, errors, reference, args.quick)
            and passed
        )

    print(f"Summary: {POST_DIR / 'verification_summary.md'}")
    print(f"Result: {'PASS' if passed else 'FAIL'}")
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
