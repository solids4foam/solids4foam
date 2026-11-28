#!/usr/bin/env python3
"""Mesh refinement studies of the high-order newtonIcoFluid residual.

Both studies solve the Stokes equations (convection off) on the static
uniform meshes of the tutorial, with the alpha pressure stabilisation, and
compare the second-order residual with the high-order (kExactLeastSquares)
residual for polynomialOrder k = 1, 2 and 3. All fields start from, and are
compared with, exact cell averages written by
initialiseTaylorGreenVortex -cellAverage.

- steady: the forced Stokes problem of initialiseTaylorGreenVortex
  -steadyStokesSource, solved in one step of the steadyState time scheme.
  The exact velocity is the Taylor-Green velocity at that time and the exact
  pressure is the non-constant initial pressure. The velocity and pressure
  errors are checked (the pressure up to a constant).
- unsteady: the decaying Taylor-Green vortex of the Stokes equations, whose
  exact pressure is constant, with the backward time scheme and a time step
  proportional to h^2, so that the second-order time error falls like h^4.
  The velocity error is checked.

By default the steady study runs the second-order, k = 2 and k = 3 residuals,
and the unsteady study the second-order and k = 3 residuals.

The observed order of the L2 error between the two finest meshes must be at
least the expected order minus 0.2.
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
WORK_ROOT = VERIFICATION / "work" / "stokes"
OUTPUT_ROOT = VERIFICATION / "postProcessing"

NU = 0.1
STEADY_TIME = 0.1
UNSTEADY_END_TIME = 0.04
UNSTEADY_DT_COEFF = 0.2
TOLERANCE = 0.2

# Configurations: name -> polynomial order (0 for the second-order residual)
CONFIGS = {"secondOrder": 0, "k1": 1, "k2": 2, "k3": 3}

# Minimum orders: (velocity, pressure); None if the error is not checked.
# The k-exact residual gives k + 1 for odd k; for k = 2 it gives 2, as the
# solid manufactured solution does with a quadratic kExact reconstruction
EXPECTED = {
    "steady": {
        "secondOrder": (2, 2),
        "k1": (2, 2),
        "k2": (2, 2),
        "k3": (4, 3),
    },
    "unsteady": {
        "secondOrder": (2, None),
        "k1": (2, None),
        "k2": (2, None),
        "k3": (4, None),
    },
}

DEFAULT_LEVELS = {"steady": [2, 3, 4, 5], "unsteady": [2, 3, 4]}

# Configurations run by default. k = 1 is left out of both: its steady linear
# solves stall at the iteration limit (the solution still converges, in about
# 30 minutes on the finest mesh), and on the unsteady meshes it has not reached
# its asymptotic order (1.76 between 20 and 40 cells per side). k = 2 is
# second order and is left out of the unsteady study. Run them with --configs
DEFAULT_CONFIGS = {
    "steady": ["secondOrder", "k2", "k3"],
    "unsteady": ["secondOrder", "k3"],
}

HIGH_ORDER_COEFFS = """
    highOrderCoeffs
    {{
        highOrderResidual yes;

        velocity
        {{
            type                    kExactLeastSquares;
            polynomialOrder         {k};
            faceStencilExtraCells   20;
            weightFunctionCoeffs
            {{
                type    radiallySymmetricExponential;
                k       6;
            }}
        }}

        pressure
        {{
            type                    kExactLeastSquares;
            polynomialOrder         {k};
            faceStencilExtraCells   20;
            weightFunctionCoeffs
            {{
                type    radiallySymmetricExponential;
                k       6;
            }}
        }}
    }}
"""


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


def replace_once(path: Path, old: str, new: str) -> None:
    text = path.read_text()
    if text.count(old) != 1:
        raise RuntimeError(f"Expected one '{old.strip()}' in {path}")
    path.write_text(text.replace(old, new))


def set_entry(path: Path, key: str, value: str) -> None:
    text = path.read_text()
    text, count = re.subn(
        rf"^({re.escape(key)}\s+)[^;]+;", rf"\g<1>{value};", text,
        flags=re.MULTILINE,
    )
    if count != 1:
        raise RuntimeError(f"Expected one '{key}' entry in {path}")
    path.write_text(text)


def build() -> None:
    WORK_ROOT.mkdir(parents=True, exist_ok=True)
    run(["wmake", "initialiseTaylorGreenVortex"], TUTORIAL,
        WORK_ROOT / "log.wmake")
    run(["./Allwmake"], TUTORIAL / "src", WORK_ROOT / "log.Allwmake")


def prepare_case(case: Path, study: str, config: str, level: int) -> float:
    """Set up one case and return its end time."""
    if case.exists():
        shutil.rmtree(case)
    case.mkdir(parents=True)
    for name in ("constant", "system", "0.orig"):
        shutil.copytree(TUTORIAL / name, case / name, symlinks=True)
    shutil.copytree(case / "0.orig", case / "0")

    dynamic_mesh_dict = case / "constant" / "dynamicMeshDict"
    dynamic_mesh_dict.unlink()
    os.symlink("dynamicMeshDict.static", dynamic_mesh_dict)

    fluid_properties = case / "constant" / "fluidProperties"
    replace_once(
        fluid_properties, "            type    RhieChow;\n",
        "            type    alpha;\n",
    )
    extra = "\n    convection no;\n"
    if CONFIGS[config]:
        extra += HIGH_ORDER_COEFFS.format(k=CONFIGS[config])
    replace_once(
        fluid_properties, "    pRefValue 0.0;\n", "    pRefValue 0.0;\n" + extra
    )

    fv_solution = case / "system" / "fvSolution"
    replace_once(fv_solution, 'snes_rtol "1e-3";',
                 'snes_rtol "1e-10";\n            snes_atol "1e-13";')
    replace_once(fv_solution, 'snes_stol "1e-3";', 'snes_stol "1e-14";')
    replace_once(fv_solution, 'ksp_rtol "1e-2";', 'ksp_rtol "1e-10";')

    control_dict = case / "system" / "controlDict"
    h = 1.0/(5*2**(level - 1))
    if study == "steady":
        end_time = STEADY_TIME
        delta_t = STEADY_TIME
        fv_schemes = case / "system" / "fvSchemes"
        replace_once(fv_schemes, "    ddt(U)     backward;",
                     "    ddt(U)     steadyState;")
        replace_once(fv_schemes, "    ddt(p)     backward;",
                     "    ddt(p)     steadyState;")
    else:
        end_time = UNSTEADY_END_TIME
        delta_t = UNSTEADY_DT_COEFF*h*h
    set_entry(control_dict, "endTime", f"{end_time:.12g}")
    set_entry(control_dict, "deltaT", f"{delta_t:.12g}")
    set_entry(control_dict, "writeInterval", f"{end_time:.12g}")
    set_entry(control_dict, "writePrecision", "15")
    set_entry(control_dict, "timePrecision", "12")

    run(
        ["blockMesh", "-dict", f"system/blockMeshDict.uniform.{level}"],
        case, case / "log.blockMesh",
    )
    command = ["initialiseTaylorGreenVortex", "-cellAverage"]
    if study == "steady":
        command += ["-steadyStokesSource", f"{STEADY_TIME:.12g}"]
    run(command, case, case / "log.initialiseTaylorGreenVortex")
    return end_time


def read_internal(path: Path) -> list[list[float]]:
    text = path.read_text()
    start = text.index("internalField")
    open_paren = text.index("(", text.index("\n", start))
    count = int(text[start:open_paren].split()[-1])
    body = text[open_paren + 1:]
    if "List<vector>" in text[start:open_paren]:
        values = re.findall(r"\(([^()]*)\)", body)[:count]
        return [[float(v) for v in item.split()] for item in values]
    return [[float(v)] for v in body.split()[:count]]


def time_name(case: Path, time: float) -> Path:
    for entry in case.iterdir():
        try:
            if abs(float(entry.name) - time) < 1e-9*max(1.0, time):
                return entry
        except ValueError:
            continue
    raise RuntimeError(f"No time directory {time} in {case}")


def errors(case: Path, study: str, end_time: float) -> tuple[float, float]:
    """L2 errors of U and of p (up to a constant), uniform cell volumes."""
    u_exact0 = read_internal(case / "0" / "U")
    p_exact = read_internal(case / "0" / "p")
    result = time_name(case, end_time)
    u = read_internal(result / "U")
    p = read_internal(result / "p")

    scale = math.exp(-2*math.pi**2*NU*end_time)
    u_error = math.sqrt(
        sum(
            sum((a - scale*b)**2 for a, b in zip(ui, ue))
            for ui, ue in zip(u, u_exact0)
        )/len(u)
    )

    if study != "steady":
        return u_error, math.nan

    diff = [pi[0] - pe[0] for pi, pe in zip(p, p_exact)]
    mean = sum(diff)/len(diff)
    p_error = math.sqrt(sum((d - mean)**2 for d in diff)/len(diff))
    return u_error, p_error


def observed_orders(n_cells: list[int], values: list[float]) -> list[float]:
    orders = []
    for i in range(1, len(values)):
        orders.append(
            math.log(values[i - 1]/values[i])
            /(0.5*math.log(n_cells[i]/n_cells[i - 1]))
        )
    return orders


def run_study(study: str, configs: list[str], levels: list[int]) -> list[dict]:
    rows = []
    for config in configs:
        for level in levels:
            case = WORK_ROOT / study / f"{config}.{level}"
            print(f"Running {study} {config} level {level}", flush=True)
            end_time = prepare_case(case, study, config, level)
            run(["solids4Foam"], case, case / "log.solids4Foam")
            u_error, p_error = errors(case, study, end_time)
            rows.append(
                {
                    "study": study,
                    "config": config,
                    "level": level,
                    "nCells": len(read_internal(case / "0" / "U")),
                    "U_L2": u_error,
                    "p_L2": p_error,
                }
            )
    return rows


def report(study: str, configs: list[str], rows: list[dict]) -> bool:
    passed = True
    print(f"\n{study} Stokes study")
    for config in configs:
        selected = [r for r in rows if r["config"] == config]
        n_cells = [r["nCells"] for r in selected]
        expected_u, expected_p = EXPECTED[study][config]
        for field, expected in (("U", expected_u), ("p", expected_p)):
            if expected is None:
                continue
            values = [r[f"{field}_L2"] for r in selected]
            orders = observed_orders(n_cells, values)
            ok = orders[-1] >= expected - TOLERANCE
            passed = passed and ok
            print(
                f"  {config:11s} {field} L2 "
                + " ".join(f"{v:9.3e}" for v in values)
                + "  order " + " ".join(f"{o:5.2f}" for o in orders)
                + f"  [{'PASS' if ok else 'FAIL'}, expected at least "
                + f"{expected}]"
            )
    return passed


def write_csv(rows: list[dict]) -> Path:
    OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)
    path = OUTPUT_ROOT / "stokes_errors.csv"
    with path.open("w", newline="") as output:
        writer = csv.DictWriter(output, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    return path


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "--studies", nargs="+", default=["steady", "unsteady"],
        choices=["steady", "unsteady"], help="studies to run",
    )
    parser.add_argument(
        "--configs", nargs="+", choices=list(CONFIGS),
        help="residual configurations to run (default: per study)",
    )
    parser.add_argument(
        "--levels", nargs="+", type=int,
        help="uniform mesh levels (blockMeshDict.uniform.<level>)",
    )
    args = parser.parse_args()

    build()

    rows = []
    passed = True
    for study in args.studies:
        levels = sorted(args.levels or DEFAULT_LEVELS[study])
        if len(levels) < 2:
            parser.error("at least two levels are needed for an order")
        configs = args.configs or DEFAULT_CONFIGS[study]
        study_rows = run_study(study, configs, levels)
        rows += study_rows
        passed = report(study, configs, study_rows) and passed

    print(f"\nErrors written to {write_csv(rows)}")
    print("\nOverall result:", "PASSED" if passed else "FAILED")
    return 0 if passed else 1


if __name__ == "__main__":
    sys.exit(main())
