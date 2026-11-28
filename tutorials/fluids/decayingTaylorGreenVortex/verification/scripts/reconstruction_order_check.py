#!/usr/bin/env python3
"""Mesh refinement study of the k-exact reconstruction of the decaying
Taylor-Green vortex.

For each mesh family (uniform, graded) and refinement level, the driver copies
the case to verification/work/<family>.<level>, creates the mesh with
system/blockMeshDict.<family>.<level> and runs
checkTaylorGreenVortexReconstruction, which sets the cell values of U and p to
the exact cell averages and measures the errors of the reconstructed values
and gradients for polynomialOrder k = 1, 2 and 3.

The observed order between two levels is computed from the L2 errors, with
the cell size h proportional to nCells^(-1/2). On the two finest meshes it
must be at least the k-exact order minus 0.2: k + 1 for the values and k for
the gradients. Higher orders pass, as symmetric stencils on smooth meshes can
cancel the leading error term. The boundary-face errors are reported but not
checked.
"""

from __future__ import annotations

import argparse
import csv
import math
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

SCRIPT = Path(__file__).resolve()
VERIFICATION = SCRIPT.parents[1]
TUTORIAL = VERIFICATION.parent
WORK_ROOT = VERIFICATION / "work"
OUTPUT_ROOT = VERIFICATION / "postProcessing"
UTILITY = "checkTaylorGreenVortexReconstruction"

POLYNOMIAL_ORDERS = (1, 2, 3)
FIELDS = ("U", "p")
MEASURES = (
    "cellValue",
    "cellGrad",
    "internalFaceValue",
    "internalFaceGrad",
    "boundaryFaceValue",
    "boundaryFaceGrad",
)
CHECKED_MEASURES = (
    "cellValue",
    "cellGrad",
    "internalFaceValue",
    "internalFaceGrad",
)
TOLERANCE = 0.2

NORM_LINE = re.compile(
    r"^\s+(U|p) (\w+) L2 = (\S+) LInf = (\S+)\s*$"
)


def openfoam_env() -> dict:
    """The environment of the OpenFOAM commands.

    macOS drops DYLD_LIBRARY_PATH when a protected interpreter such as
    /usr/bin/env starts this script; restore it as OpenFOAM does.
    """
    env = dict(os.environ)
    if sys.platform == "darwin" and not env.get("DYLD_LIBRARY_PATH") \
            and env.get("FOAM_LD_LIBRARY_PATH"):
        env["DYLD_LIBRARY_PATH"] = env["FOAM_LD_LIBRARY_PATH"]
    return env


def run(command: list[str], cwd: Path, log: Path) -> None:
    with log.open("w") as output:
        result = subprocess.run(
            command, cwd=cwd, stdout=output, stderr=subprocess.STDOUT,
            env=openfoam_env(),
        )
    if result.returncode != 0:
        raise RuntimeError(
            f"'{' '.join(command)}' failed in {cwd}; see {log}"
        )


def build() -> None:
    WORK_ROOT.mkdir(parents=True, exist_ok=True)
    run(["wmake", UTILITY], TUTORIAL, WORK_ROOT / "log.wmake")
    run(["./Allwmake"], TUTORIAL / "src", WORK_ROOT / "log.Allwmake")


def prepare_case(case: Path) -> None:
    if case.exists():
        shutil.rmtree(case)
    case.mkdir(parents=True)
    shutil.copytree(TUTORIAL / "constant", case / "constant", symlinks=True)
    shutil.copytree(TUTORIAL / "system", case / "system")
    shutil.copytree(TUTORIAL / "0.orig", case / "0")


def parse_log(log: Path) -> tuple[int, dict]:
    n_cells = None
    order = None
    errors = {}
    for line in log.read_text().splitlines():
        match = re.match(r"^Number of cells (\d+)", line)
        if match:
            n_cells = int(match.group(1))
            continue
        match = re.match(r"^polynomialOrder (\d+)", line)
        if match:
            order = int(match.group(1))
            continue
        match = NORM_LINE.match(line)
        if match and order is not None:
            field, measure, l2, linf = match.groups()
            errors[(order, field, measure)] = (float(l2), float(linf))
    if n_cells is None or len(errors) != (
        len(POLYNOMIAL_ORDERS) * len(FIELDS) * len(MEASURES)
    ):
        raise RuntimeError(f"Could not read the errors from {log}")
    return n_cells, errors


def run_level(family: str, level: int) -> tuple[int, dict]:
    case = WORK_ROOT / f"{family}.{level}"
    prepare_case(case)
    run(
        ["blockMesh", "-dict", f"system/blockMeshDict.{family}.{level}"],
        case,
        case / "log.blockMesh",
    )
    run([UTILITY], case, case / f"log.{UTILITY}")
    return parse_log(case / f"log.{UTILITY}")


def observed_orders(n_cells: list[int], errors: list[float]) -> list[float]:
    orders = []
    for i in range(1, len(errors)):
        if errors[i - 1] <= 0 or errors[i] <= 0:
            orders.append(math.nan)
            continue
        orders.append(
            math.log(errors[i - 1] / errors[i])
            / (0.5*math.log(n_cells[i] / n_cells[i - 1]))
        )
    return orders


def expected_order(polynomial_order: int, measure: str) -> int:
    return polynomial_order + 1 if measure.endswith("Value") else polynomial_order


def report(family: str, levels: list[int], results: dict) -> bool:
    passed = True
    n_cells = [results[level][0] for level in levels]
    print(f"\n{family} meshes, cells: {' '.join(str(n) for n in n_cells)}")

    for k in POLYNOMIAL_ORDERS:
        print(
            f"\n  polynomialOrder {k}: expected order at least {k + 1} for "
            f"values, {k} for gradients"
        )
        for field in FIELDS:
            for measure in MEASURES:
                l2 = [results[level][1][(k, field, measure)][0] for level in levels]
                orders = observed_orders(n_cells, l2)
                expected = expected_order(k, measure)
                if measure in CHECKED_MEASURES:
                    ok = orders[-1] >= expected - TOLERANCE
                    status = "PASS" if ok else "FAIL"
                    passed = passed and ok
                else:
                    status = "reported only"
                print(
                    f"    {field} {measure:18s} L2 "
                    + " ".join(f"{e:9.3e}" for e in l2)
                    + "  order "
                    + " ".join(f"{o:5.2f}" for o in orders)
                    + f"  [{status}]"
                )
    return passed


def write_csv(results: dict) -> Path:
    OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)
    path = OUTPUT_ROOT / "reconstruction_errors.csv"
    with path.open("w", newline="") as output:
        writer = csv.writer(output)
        writer.writerow(
            ["family", "level", "nCells", "polynomialOrder", "field",
             "measure", "L2", "LInf"]
        )
        for (family, level), (n_cells, errors) in sorted(results.items()):
            for (k, field, measure), (l2, linf) in sorted(errors.items()):
                writer.writerow(
                    [family, level, n_cells, k, field, measure, l2, linf]
                )
    return path


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "--families", nargs="+", default=["uniform", "graded"],
        choices=["uniform", "graded"], help="mesh families to run",
    )
    parser.add_argument(
        "--levels", nargs="+", type=int, default=[2, 3, 4, 5, 6],
        help="refinement levels to run (blockMeshDict.<family>.<level>)",
    )
    args = parser.parse_args()

    if len(args.levels) < 2:
        parser.error("at least two levels are needed for an order")

    levels = sorted(args.levels)
    build()

    results = {}
    for family in args.families:
        for level in levels:
            print(f"Running {family}.{level}", flush=True)
            results[(family, level)] = run_level(family, level)

    passed = True
    for family in args.families:
        family_results = {
            level: results[(family, level)] for level in levels
        }
        passed = report(family, levels, family_results) and passed

    print(f"\nErrors written to {write_csv(results)}")
    print("\nOverall result:", "PASSED" if passed else "FAILED")
    return 0 if passed else 1


if __name__ == "__main__":
    sys.exit(main())
