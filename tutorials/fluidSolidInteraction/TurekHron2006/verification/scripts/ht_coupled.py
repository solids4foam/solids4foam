#!/usr/bin/env python3
"""Coupled FSI3 case of the production study (corrected fluid set-up).

Builds a partitioned FSI3 case (IQN-ILS, Dirichlet-Neumann) from the
TurekHron2006 tutorial with the corrections of the second-order replay study:

  fluid mesh   ring+box family of ht_mesh_family.py (--ring-n 22 --rows-n 16
               --tip-n 4 --aft-n 60 --wake-n 70), level --level;
  solid mesh   the tutorial's 105 x 6 plate times --solid-level;
  convection   div(phi,U) Gauss linear (central, unlimited);
  mesh motion  displacementLaplacian, same quadratic inverseDistance
               diffusivity, driven by the interface through pointDisplacement
               (solids4foam PR #554), so the interior mesh is a function of
               the interface position only (no creep);
  wall p       movingWallPressure on the flag (dp/dn = -n.a_wall);
  coupling     requireAllResidualMeasures yes (max of the two residual
               measures, issue #552), outerCorrTolerance --fsi-tol, no reuse;
  solvers      fluid absolute tolerances --fluid-tol (p, U, cellDisplacement),
               non-final p solves relTol 0.01, PIMPLE 3/3/1 as the tutorial;
               p by GAMG (DICGaussSeidel smoother) unless --p-solver gs.

The FSI3 protocol is that of the tutorial and of the earlier refinement
study: impulsive parabolic inflow (Um = 2 m/s), coupling switched on at 2 s,
run to --end. Binary output every --write-every seconds with purgeWrite 3
(checkpoints). Extra function object: forcesPlate (flag only).

    python3 ht_coupled.py --tutorial ~/ht_prod/s4f/tutorials/fluidSolidInteraction/TurekHron2006 \
        --level 2 --solid-level 4 --delta-t 5e-4 --cores 16 --out ~/ht_prod/work/c_2x
"""
from __future__ import annotations

import argparse
import re
import shutil
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import ht_mesh_family as mf  # noqa: E402

RINGBOX = dict(rows_n=16, rows_ratio=1.0, tip_n=4, tip_ratio=1.0, aft_n=60, aft_ratio=1.0,
               wake_n=70, wake_ratio=1.0, ring_n=22)


def sub1(text: str, old: str, new: str, what: str, regex: bool = False) -> str:
    if regex:
        out, n = re.subn(old, new, text, flags=re.M)
    else:
        n = text.count(old)
        out = text.replace(old, new)
    if n != 1:
        raise SystemExit(f"{what}: expected one match of {old!r}, found {n}")
    return out


def edit(path: Path, fn) -> None:
    path.write_text(fn(path.read_text()))


def p_solver(kind: str, tol: float) -> str:
    if kind == "gs":
        smoother = "GaussSeidel"
        coarse = 20
    else:
        smoother = "DICGaussSeidel"
        coarse = 200
    return f"""    pFinal
    {{
        solver          GAMG;
        tolerance       {tol:g};
        relTol          0;
        minIter         1;
        maxIter         1000;
        smoother        {smoother};
        nPreSweeps      0;
        nPostSweeps     2;
        nFinestSweeps   2;
        scaleCorrection true;
        directSolveCoarsest false;
        cacheAgglomeration true;
        nCellsInCoarsestLevel {coarse};
        agglomerator    faceAreaPair;
        mergeLevels     1;
    }}

    p
    {{
        $pFinal;
        relTol          0.01;
    }}
"""


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--tutorial", type=Path, required=True)
    ap.add_argument("--level", type=int, required=True, help="ring+box fluid level")
    ap.add_argument("--solid-level", type=int, required=True)
    ap.add_argument("--delta-t", type=float, required=True)
    ap.add_argument("--end", type=float, default=7.0)
    ap.add_argument("--write-every", type=float, default=0.25, help="checkpoint interval (s)")
    ap.add_argument("--fsi-tol", type=float, default=1e-6)
    ap.add_argument("--n-outer", type=int, default=50)
    ap.add_argument("--fluid-tol", type=float, default=1e-9)
    ap.add_argument("--p-solver", choices=("dic", "gs"), default="dic")
    ap.add_argument("--cores", type=int, default=1)
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args()

    tut, case = a.tutorial.resolve(), a.out.resolve()
    if case.exists():
        raise SystemExit(f"{case} exists")
    case.mkdir(parents=True)
    for d in ("0", "constant", "system"):
        shutil.copytree(tut / d, case / d, symlinks=False)
    for junk in ("0/fluid/U.robin", "0/fluid/p.robin", "0/fluid/U.dirichletNeumann",
                 "0/fluid/p.dirichletNeumann", "constant/fsiProperties.robin",
                 "constant/fsiProperties.iqnils"):
        (case / junk).unlink(missing_ok=True)

    # FSI3 parameters (the tutorial's stored values; checked, not assumed)
    s = (case / "0/fluid/U").read_text()
    if not re.search(r"maxValue\s+3;", s) or not re.search(r"transitionPeriod\s+0;", s):
        raise SystemExit("0/fluid/U is not the FSI3 inflow (maxValue 3, transitionPeriod 0)")
    s = (case / "constant/solid/mechanicalProperties").read_text()
    if not re.search(r"rho\s+rho \[[^]]*\] 1000;", s) or not re.search(r"E\s+E \[[^]]*\] 5.6e6;", s):
        raise SystemExit("mechanicalProperties is not FSI3 (rho 1000, E 5.6e6)")

    # Meshes
    base = (case / "system/fluid/blockMeshDict").read_text()
    (case / "system/fluid/blockMeshDict").write_text(mf.build(base, a.level, **RINGBOX))
    edit(case / "system/solid/blockMeshDict",
         lambda t: sub1(t, "(105 6 1)", f"({105 * a.solid_level} {6 * a.solid_level} 1)", "solid mesh"))

    # Fluid schemes, mesh motion, wall pressure
    def schemes(t: str) -> str:
        t = sub1(t, "div(phi,U) Gauss linearUpwind cellLimited leastSquares 1;", "div(phi,U) Gauss linear;",
                 "central")
        return sub1(t, "    laplacian(diffusivity,cellMotionU) Gauss linear corrected;",
                    "    laplacian(diffusivity,cellDisplacement) Gauss linear corrected;", "dispLap schemes")
    edit(case / "system/fluid/fvSchemes", schemes)
    edit(case / "constant/fluid/dynamicMeshDict",
         lambda t: sub1(t, '"solver|motionSolver" velocityLaplacian;',
                        "motionSolverLibs (fvMotionSolvers);\n\nmotionSolver displacementLaplacian;",
                        "dispLap"))
    pm = case / "0/fluid/pointMotionU"
    t = pm.read_text()
    t = sub1(t, "object      pointMotionU;", "object      pointDisplacement;", "pointDisplacement")
    t = sub1(t, "dimensions      [0 1 -1 0 0 0 0];", "dimensions      [0 1 0 0 0 0 0];", "pointDisplacement")
    (case / "0/fluid/pointDisplacement").write_text(t)
    pm.unlink()

    def pwall(t: str) -> str:
        return sub1(t, r"(plate\s*\{\s*)type\s+zeroGradient;", r"\1type            movingWallPressure;\n"
                    "        gradient        uniform 0;", "movingWallPressure", regex=True)
    edit(case / "0/fluid/p", pwall)

    # Fluid solvers
    def solvers(t: str) -> str:
        i = t.index('    "p|pFinal"')
        j = t.index("    cellMotionU")
        t = t[:i] + p_solver(a.p_solver, a.fluid_tol) + "\n" + t[j:]
        t = sub1(t, "    cellMotionU\n", '    "cellMotionU|cellDisplacement"\n', "motion solver")
        t, n = re.subn(r"(?m)^(\s*tolerance\s+)1e-0?6;", rf"\g<1>{a.fluid_tol:g};", t)
        if n != 2:
            raise SystemExit(f"fvSolution: expected two remaining tolerances, found {n}")
        return t
    edit(case / "system/fluid/fvSolution", solvers)

    # Coupling
    def fsi(t: str) -> str:
        t = sub1(t, r"^(\s*outerCorrTolerance\s+)[^;]*;", rf"\g<1>{a.fsi_tol:g};", "fsi tol", regex=True)
        t = sub1(t, r"^(\s*nOuterCorr\s+)[^;]*;", rf"\g<1>{a.n_outer};", "nOuterCorr", regex=True)
        t = sub1(t, r"^(\s*couplingStartTime\s+)[^;]*;", r"\g<1>2.0;", "couplingStartTime", regex=True)
        return sub1(t, r"^(\s*writeResidualsToFile\s+yes;)",
                    r"\1\n\n    // Strict convergence test (issue #552): both residual measures\n"
                    r"    requireAllResidualMeasures yes;", "strict", regex=True)
    edit(case / "constant/fsiProperties", fsi)
    if re.search(r"(?m)^\s*couplingReuse\s", (case / "constant/fsiProperties").read_text()):
        raise SystemExit("couplingReuse must stay off")

    # Time control: binary checkpoints with purgeWrite
    nwrite = max(1, round(a.write_every / a.delta_t))

    def control(t: str) -> str:
        t = sub1(t, r"^deltaT .*$", f"deltaT          {a.delta_t:g};", "deltaT", regex=True)
        t = sub1(t, r"^endTime .*$", f"endTime         {a.end:g};", "endTime", regex=True)
        t = sub1(t, r"^writeInterval .*$", f"writeInterval   {nwrite};", "writeInterval", regex=True)
        t = sub1(t, r"^purgeWrite .*$", "purgeWrite      3;", "purgeWrite", regex=True)
        t = sub1(t, r"^writeFormat .*$", "writeFormat     binary;", "writeFormat", regex=True)
        t = sub1(t, r"^writePrecision .*$", "writePrecision  12;", "writePrecision", regex=True)
        return sub1(t, r"^timePrecision .*$", "timePrecision   8;", "timePrecision", regex=True)
    edit(case / "system/controlDict", control)

    def functions(t: str) -> str:
        k = t.rindex("}")
        return t[:k] + """
   forcesPlate
   {
       type                forces;
       libs                ( "libforces.so" );
       writeControl        timeStep;
       writeInterval       1;
       patches             (plate);
       "pName|p"           p;
       "UName|U"           U;
       "rhoName|rho"       rhoInf;
       log                 false;
       rhoInf              1000;
       CofR                (0.5 0.1 0);

       region fluid;
   }
""" + t[k:]
    edit(case / "system/functions", functions)

    for dp in (case / "system/decomposeParDict", case / "system/fluid/decomposeParDict",
               case / "system/solid/decomposeParDict"):
        if dp.is_file():
            t = re.sub(r"(?m)^numberOfSubdomains\s+\d+;", f"numberOfSubdomains {a.cores};", dp.read_text())
            dp.write_text(re.sub(r"(?m)^method\s+\w+;", "method          scotch;", t))

    def run(cmd: list[str], log: str) -> None:
        with (case / log).open("w") as fh:
            if subprocess.run(cmd, cwd=case, stdout=fh, stderr=subprocess.STDOUT).returncode:
                raise SystemExit(f"{' '.join(cmd)} failed, see {case / log}")

    run(["blockMesh", "-region", "solid"], "log.blockMesh.solid")
    run(["blockMesh", "-region", "fluid"], "log.blockMesh.fluid")
    run(["checkMesh", "-region", "fluid", "-noFunctionObjects"], "log.checkMesh.fluid")
    if a.cores > 1:
        run(["decomposePar", "-region", "fluid"], "log.decomposePar.fluid")
        run(["decomposePar", "-region", "solid"], "log.decomposePar.solid")
    (case / "coupled_setup.txt").write_text(
        "".join(f"{k} {v}\n" for k, v in vars(a).items()) + "".join(f"{k} {v}\n" for k, v in RINGBOX.items()))
    print(case, "ready")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
