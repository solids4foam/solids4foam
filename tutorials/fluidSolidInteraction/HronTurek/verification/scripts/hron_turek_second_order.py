#!/usr/bin/env python3
"""Replay variants for the second-order study of the FSI3 flag load.

Builds a restart case of a trajectory replay with hron_turek_energy.py (same
energy/generalised-force function object as the energy-balance study), then
applies the variants of hron_turek_variants.py and the new ones below. Each
was checked against OpenFOAM v2512 / solids4foam:

  lu_unlim  div(phi,U): Gauss linearUpwind cellLimited leastSquares 1
            -> Gauss linearUpwind grad(U)   (no gradient limiter)
  linear    div(phi,U) -> Gauss linear      (central, unlimited)
  dispLap   mesh motion velocityLaplacian -> displacementLaplacian with the
            same quadratic inverseDistance diffusivity; the plate condition
            prescribes the replayed displacement itself (pointDisplacement),
            so the interior mesh is a function of the boundary position only
            (no path dependence / creep). Only valid from the undeformed
            rigid state (restart 25).
  tight     p, U tolerances 1e-6 -> 1e-9 (relTol 0), mesh 1e-6 -> 1e-9
  noc4      nOuterCorrectors 3 -> 4 and nCorrectors 3 -> 4
  conv      tight, plus PIMPLE nOuterCorrectors/nCorrectors/nNonOrthogonal
            3/3/1 -> 6/6/3 (converged segregated solve per step)

    python3 hron_turek_second_order.py --source <replay> --restart 25 --end 33.2 \
        --name s_1x_dispLap --variants dispLap,lu_unlim --work ~/ht_2nd/work
"""
from __future__ import annotations

import argparse
import re
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import hron_turek_variants as hv  # noqa: E402

LIMITED = "div(phi,U) Gauss linearUpwind cellLimited leastSquares 1;"


def tdirs(case: Path, restart: str):
    procs = sorted(case.glob("processor*"))
    return [p / restart for p in procs] or [case / restart]


def disp_lap(case: Path, restart: str) -> None:
    if float(restart) != 25.0:
        raise SystemExit("dispLap needs the undeformed rigid state (restart 25)")
    for td in tdirs(case, restart):
        f = td / "pointMotionU"
        s = f.read_text()
        s = hv.sub_once(s, "object      pointMotionU;", "object      pointDisplacement;", "dispLap")
        s = hv.sub_once(s, "dimensions      [0 1 -1 0 0 0 0];", "dimensions      [0 1 0 0 0 0 0];", "dispLap")
        s = hv.sub_once(s, "v[i] = (d1 - d0)/dt;", "v[i] = d1;", "dispLap")
        (td / "pointDisplacement").write_text(s)
        f.unlink()
    f = case / "constant" / "dynamicMeshDict"
    s = f.read_text()
    s = hv.sub_once(s, '"solver|motionSolver" velocityLaplacian;',
                    'motionSolverLibs (fvMotionSolvers);\n\nmotionSolver displacementLaplacian;', "dispLap")
    f.write_text(s)
    f = case / "system" / "fvSolution"
    s = f.read_text()
    s = hv.sub_once(s, "    cellMotionU\n", '    "cellMotionU|cellDisplacement"\n', "dispLap")
    f.write_text(s)
    f = case / "system" / "fvSchemes"
    s = f.read_text()
    s = hv.sub_once(s, "    laplacian(diffusivity,cellMotionU) Gauss linear corrected;",
                    "    laplacian(diffusivity,cellMotionU) Gauss linear corrected;\n"
                    "    laplacian(diffusivity,cellDisplacement) Gauss linear corrected;", "dispLap")
    f.write_text(s)


def apply(case: Path, restart: str, variant: str) -> None:
    fs = case / "system" / "fvSchemes"
    fv = case / "system" / "fvSolution"
    if variant == "lu_unlim":
        fs.write_text(hv.sub_once(fs.read_text(), LIMITED, "div(phi,U) Gauss linearUpwind grad(U);", variant))
    elif variant == "linear":
        fs.write_text(hv.sub_once(fs.read_text(), LIMITED, "div(phi,U) Gauss linear;", variant))
    elif variant == "dispLap":
        disp_lap(case, restart)
    elif variant == "tight":
        s = fv.read_text()
        n = s.count("tolerance       1e-06;")
        if n != 3:
            raise SystemExit("tight: expected three tolerances of 1e-06")
        fv.write_text(s.replace("tolerance       1e-06;", "tolerance       1e-09;"))
    elif variant == "noc4":
        s = fv.read_text()
        s = hv.sub_once(s, "nOuterCorrectors    3;", "nOuterCorrectors    4;", variant)
        s = hv.sub_once(s, "nCorrectors         3;", "nCorrectors         4;", variant)
        fv.write_text(s)
    elif variant == "conv":
        apply(case, restart, "tight")
        s = fv.read_text()
        s = hv.sub_once(s, "nOuterCorrectors    3;", "nOuterCorrectors    6;", variant)
        s = hv.sub_once(s, "nCorrectors         3;", "nCorrectors         6;", variant)
        s = hv.sub_once(s, "nNonOrthogonalCorrectors 1;", "nNonOrthogonalCorrectors 3;", variant)
        fv.write_text(s)
    else:
        hv.apply(case, restart, variant)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--source", required=True, type=Path)
    ap.add_argument("--restart", required=True)
    ap.add_argument("--end", required=True)
    ap.add_argument("--name", required=True)
    ap.add_argument("--variants", default="base", help="comma-separated")
    ap.add_argument("--work", type=Path, required=True)
    ap.add_argument("--amplitude-factor", type=float, default=1.0)
    ap.add_argument("--blend", type=float, default=0.3)
    a = ap.parse_args()
    subprocess.run([sys.executable, str(HERE / "hron_turek_energy.py"), "--source", str(a.source),
                    "--restart", a.restart, "--end", a.end, "--name", a.name, "--work", str(a.work),
                    "--amplitude-factor", repr(a.amplitude_factor), "--blend", repr(a.blend)],
                   check=True)
    case = (a.work / a.name).resolve()
    for v in a.variants.split(","):
        if v != "base":
            apply(case, a.restart, v)
    with open(case / "energy_setup.txt", "a") as fh:
        fh.write(f"variants {a.variants}\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
