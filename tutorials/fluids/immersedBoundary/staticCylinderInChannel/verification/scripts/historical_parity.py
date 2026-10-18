#!/usr/bin/env python3
"""
Historical parity of the incremental method: the static cylinder cases of the
openHFDIB-DEM pimpleHFDIBFoam solver (the manuscript-era cases), run with the
incremental method of immersedBoundaryForce and the solids4foam pimpleFluid
fluid model, with the historical numerical settings.

The historical case (read only) gives the mesh (blockMeshDict, scaled by 10:
D = 1, channel 22 x 4.1, centre (2, 2)), the cylinder surface, the fluid
properties (rho = 1000, nu = 0.05, Re = 20), the initial and boundary
conditions (including its parabolic inlet, coded from the extreme inlet face
centres), fvSchemes and fvSolution (nOuterCorrectors 5, nCorrectors 3,
nNonOrthogonalCorrectors 2, relaxation U 0.7, p 0.3), and the time step and
end time (0.005 s, 20 s); constant/HFDIBDEMDict gives ibCouplingCoeff 0.2,
sdBasedLambda false (vertexFraction occupancy) and surfaceThreshold 1e-4.
The inlet may be replaced by the exact parabola over the channel height,
u = 6 Ubar y (H - y)/H^2 with Ubar = 1 and H = 4.1 (--inlet canonical).

Differences from the historical cases, which do not change the solution of a
static body: the mesh is static (staticFvMesh) rather than a motion solver
with zero displacement, the cache of grad(U) is removed from fvSolution (it
is not supported by the solids4foam fluid models with v2412, and only avoids
recomputing the gradient), the fields are written every 5 s, and the
decomposition is over fewer ranks.

Usage:
    historical_parity.py generate --old <historical case> --name <case>
        --inlet historical|canonical --ranks N
    historical_parity.py analyse
"""

import argparse
import csv
import json
import re
import shutil
import subprocess
import sys
from datetime import datetime
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import static_cylinder_verification as scv  # noqa: E402

STUDY = "final-B"
H = 4.1
UBAR = 1.0
# Cd = Fx/(0.5 rho Ubar^2 D Lz), with rho = 1000, D = 1, Lz = 0.1
FORCE_SCALE = 1.0/(0.5*1000.0*UBAR**2*1.0*0.1)

HEADER = scv.FOAM_HEADER


def write_dict(path: Path, location: str, obj: str, body: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(HEADER.format(location=location, object=obj) + body
                    + "\n// " + "*"*73 + " //\n")


def canonical_inlet(u_text: str) -> str:
    """Replace the coded inlet by the exact parabola over the height H"""
    a1 = 6.0*UBAR/H
    a2 = -6.0*UBAR/H**2
    new = (
        "    inlet\n    {\n"
        "        // Exact parabola u = 6 Ubar y (H - y)/H^2, Ubar = 1, H = 4.1\n"
        "        type            fixedProfile;\n"
        "        profile         polynomial\n        (\n"
        f"            (({a1:.15g} 0 0) (1 0 0))\n"
        f"            (({a2:.15g} 0 0) (2 0 0))\n"
        "        );\n"
        "        direction       (0 1 0);\n"
        "        origin          0;\n"
        "    }\n"
    )
    # The inlet entry runs from "inlet" to the closing brace at its indent
    m = re.search(r"^    inlet\s*\n    \{.*?^    \}\n", u_text,
                  re.MULTILINE | re.DOTALL)
    if not m:
        sys.exit("ERROR: inlet entry not found in 0/U")
    return u_text[:m.start()] + new + u_text[m.end():]


def generate(args) -> None:
    old = Path(args.old).resolve()
    case = scv.WORK_DIR / STUDY / args.name
    if case.exists():
        sys.exit(f"ERROR: {case} exists")
    case.mkdir(parents=True)

    # Historical files (copied, never modified in place)
    for rel in ("system/blockMeshDict", "system/fvSchemes",
                "system/fvSolution", "constant/transportProperties",
                "constant/triSurface/cylinder.stl", "0/U", "0/p"):
        (case / rel).parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(old / rel, case / rel)

    # The solids4foam fluid models do not support the caching of grad(U)
    # with OpenFOAM.com v2412; the cache only avoids recomputing the gradient
    # and does not change the solution, so the cache dictionary is removed
    fv_solution = case / "system/fvSolution"
    text, count = re.subn(r"\n[^\n]*\n?cache\s*\{[^}]*\}\s*\n", "\n",
                          fv_solution.read_text())
    if count != 1:
        sys.exit(f"ERROR: expected one cache dictionary in {fv_solution}")
    fv_solution.write_text(text)

    if args.inlet == "canonical":
        u = case / "0/U"
        u.write_text(canonical_inlet(u.read_text()))

    # solids4foam fluid-only case files from the tutorial
    for rel in ("constant/physicsProperties", "constant/fluidProperties",
                "constant/dynamicMeshDict", "constant/g",
                "constant/turbulenceProperties"):
        shutil.copy2(scv.TUTORIAL_DIR / rel, case / rel)

    # Time and output settings of the historical case
    old_control = (old / "system/controlDict").read_text()

    def old_entry(key):
        m = re.search(rf"^{key}\s+([^;]+);", old_control, re.MULTILINE)
        return m.group(1).strip()

    write_dict(case / "system/controlDict", "system", "controlDict", f"""\
libs            ("libfvMotionSolvers.so" "libmeshTools.so" immersedBoundary);

application     solids4Foam;

startFrom       startTime;

startTime       0;

stopAt          endTime;

endTime         {old_entry('endTime')};

deltaT          {old_entry('deltaT')};

writeControl    runTime;

writeInterval   5;

purgeWrite      0;

writeFormat     ascii;

writePrecision  8;

writeCompression off;

timeFormat      general;

timePrecision   8;

runTimeModifiable no;

adjustTimeStep  {old_entry('adjustTimeStep')};

maxCo           {old_entry('maxCo')};
""")

    hfdib = (old / "constant/HFDIBDEMDict").read_text()
    cc = re.search(r"^ibCouplingCoeff\s+([^;]+);", hfdib, re.MULTILINE)
    thr = re.search(r"^surfaceThreshold\s+([^;]+);", hfdib, re.MULTILINE)
    sd = re.search(r"^\s*sdBasedLambda\s+(\w+);", hfdib, re.MULTILINE)
    if not (cc and thr and sd) or sd.group(1) != "false":
        sys.exit("ERROR: cannot read ibCouplingCoeff, surfaceThreshold and "
                 "sdBasedLambda false from HFDIBDEMDict")
    write_dict(case / "constant/fvOptions", "constant", "fvOptions", f"""\
// Historical pimpleHFDIBFoam settings (constant/HFDIBDEMDict): ibCouplingCoeff,
// surfaceThreshold, sdBasedLambda false (vertex fraction occupancy)
immersedBoundary
{{
    type                immersedBoundaryForce;
    method              incremental;
    couplingCoeff       {cc.group(1).strip()};
    surfaceThreshold    {thr.group(1).strip()};

    bodies
    {{
        cylinder
        {{
            surface     "cylinder.stl";
            occupancy   vertexFraction;
        }}
    }}
}}
""")

    write_dict(case / "system/decomposeParDict", "system", "decomposeParDict",
               f"numberOfSubdomains  {args.ranks};\n\nmethod          "
               "scotch;\n")

    decompose = "" if args.ranks == 1 else (
        "decomposePar -force > log.decomposePar 2>&1\n")
    script = scv.slurm_header(
        case, f"scv-{args.name}", args.ranks, args.memory, args.walltime,
        scv.PARTITION, scv.ENV_SH
    ) + f"""\
echo "job=${{SLURM_JOB_ID:-}} host=$(hostname) ranks={args.ranks} start=$(date -Is)"
START=$(date +%s)

blockMesh > log.blockMesh 2>&1
{decompose}{scv.mpi_run(args.ranks)} > log.solids4Foam 2>&1

END=$(date +%s)
echo "end=$(date -Is) wall_s=$((END - START))" | tee wallclock.txt
"""
    (case / "run.slurm").write_text(script)
    (case / "run.slurm").chmod(0o755)

    sha = subprocess.run(["git", "-C", str(scv.TUTORIAL_DIR), "rev-parse",
                          "HEAD"], capture_output=True, text=True).stdout
    meta = {
        "study": STUDY, "case": args.name, "historical_case": str(old),
        "inlet": args.inlet, "ranks": args.ranks,
        "couplingCoeff": cc.group(1).strip(),
        "deltaT": float(old_entry("deltaT")),
        "end_time": float(old_entry("endTime")),
        "memory": args.memory, "walltime": args.walltime,
        "partition": scv.PARTITION, "env_sh": scv.ENV_SH,
        "git_sha": sha.strip(),
        "generated": datetime.now().isoformat(timespec="seconds"),
    }
    (case / "case.json").write_text(json.dumps(meta, indent=2) + "\n")
    print(f"{case}: inlet={args.inlet} ranks={args.ranks}")


def analyse(args) -> None:
    rows = []
    for case in sorted((scv.WORK_DIR / STUDY).glob("*")):
        files = scv.force_files(case)
        if not (case / "case.json").exists() or not files:
            continue
        meta = json.loads((case / "case.json").read_text())
        t, fx, fy, _ = scv.read_forces(files)
        cd = [FORCE_SCALE*f for f in fx]
        cl = [FORCE_SCALE*f for f in fy]
        t0 = t[-1] - args.window
        scd = scv.window_stats(t, cd, t0)
        scl = scv.window_stats(t, cl, t0)
        stats = scv.log_stats(sorted(case.glob("log.solids4Foam*")))
        wall = sum(int(w) for w in re.findall(
            r"wall_s=(\d+)", (case / "wallclock.txt").read_text())) \
            if (case / "wallclock.txt").exists() else float("nan")
        n_cells = re.search(r"nCells:\s*(\d+)",
                            (case / "log.blockMesh").read_text()).group(1)
        rows.append({
            "case": meta["case"], "inlet": meta["inlet"],
            "historical_case": Path(meta["historical_case"]).name,
            "n_cells": n_cells, "ranks": meta["ranks"],
            "couplingCoeff": meta["couplingCoeff"],
            "deltaT": meta["deltaT"], "end_time": t[-1],
            "Cd_mean": scd[0], "Cd_last": cd[-1],
            "Cd_rel_range": (scd[2] - scd[1])/abs(scd[0]),
            "Cl_mean": scl[0], "Cl_range": scl[2] - scl[1],
            "max_co": stats.get("max_co"), "wall_s": wall,
            "cpu_hours": wall*meta["ranks"]/3600.0,
        })
    out = scv.OUTPUT_DIR / STUDY
    out.mkdir(parents=True, exist_ok=True)
    with (out / "historical_parity_runs.csv").open("w", newline="") as h:
        w = csv.DictWriter(h, fieldnames=list(rows[0]), lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    for r in rows:
        print(f"{r['case']:28s} t={r['end_time']:g} Cd={r['Cd_mean']:.5f} "
              f"(last {r['Cd_last']:.5f}) relRange={r['Cd_rel_range']:.1e} "
              f"Cl={r['Cl_mean']:.5f}")


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    sub = parser.add_subparsers(dest="command", required=True)
    g = sub.add_parser("generate")
    g.add_argument("--old", required=True)
    g.add_argument("--name", required=True)
    g.add_argument("--inlet", choices=("historical", "canonical"),
                   required=True)
    g.add_argument("--ranks", type=int, default=1)
    g.add_argument("--memory", default="4G")
    g.add_argument("--walltime", default="12:00:00")
    g.set_defaults(func=generate)
    a = sub.add_parser("analyse")
    a.add_argument("--window", type=float, default=1.0,
                   help="final window (s); the historical value is the mean "
                        "over 19 < t < 20 s")
    a.set_defaults(func=analyse)
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
