#!/usr/bin/env python3
"""Turek-Hron CFD3 (rigid flag, Re = 200) on the HronTurek FSI3 fluid meshes.

CFD3 is the fluid-only part of the benchmark: the FSI3 geometry and fluid, with
the flag held rigid, so that the cylinder and flag form one no-slip obstacle.
The case is built from the HronTurek tutorial's fluid region (blockMeshDict,
schemes, solver settings, inlet and outlet conditions), refined by the same
factor as the FSI3 mesh study, with a static mesh and no solid or coupling.

    python3 scripts/hron_turek_cfd3.py --level 2 --delta-t 0.0005 --end-time 8 --cores 8

The run is written to verification/work/cfd3_<level>x_dt<dt>; analysis is done
by cfd3_analysis.py.
"""

from __future__ import annotations

import argparse
import re
import shutil
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import hron_turek_verification as driver  # noqa: E402

TUTORIAL = driver.TUTORIAL
WORK_ROOT = driver.WORK_ROOT


def case_name(level: int, delta_t: float, tag: str = "") -> str:
    return f"cfd3_{level}x_dt{delta_t:g}{tag}"


def build_case(name: str, level: int, delta_t: float, end_time: float,
               cores: int) -> Path:
    case = WORK_ROOT / name
    if case.exists():
        shutil.rmtree(case)
    for sub in ("0", "constant", "system"):
        (case / sub).mkdir(parents=True)

    def copy(source: Path, destination: Path) -> None:
        shutil.copyfile(source.resolve(), destination)

    # Single-region fluid case from the tutorial's fluid region
    for field in ("U", "p"):
        copy(TUTORIAL / "0/fluid" / f"{field}.dirichletNeumann", case / "0" / field)
    for entry in ("fluidProperties", "transportProperties", "turbulenceProperties"):
        copy(TUTORIAL / "constant/fluid" / entry, case / "constant" / entry)
    copy(TUTORIAL / "constant/g", case / "constant/g")
    for entry in ("blockMeshDict", "fvSchemes", "fvSolution"):
        copy(TUTORIAL / "system/fluid" / entry, case / "system" / entry)
    copy(TUTORIAL / "system/controlDict", case / "system/controlDict")
    copy(TUTORIAL / "system/functions", case / "system/functions")
    copy(TUTORIAL / "system/fluid/decomposeParDict", case / "system/decomposeParDict")
    (case / "constant/physicsProperties").write_text(
        "FoamFile\n{\n    version 2.0;\n    format ascii;\n    class dictionary;\n"
        "    object physicsProperties;\n}\n\ntype    fluid;\n")
    # Rigid flag: the plate is a fixed no-slip wall like the cylinder
    u = case / "0/U"
    text = u.read_text()
    text, count = re.subn(r"newMovingWallVelocity", "fixedValue", text)
    if count != 1:
        driver.fail("Could not make the plate a fixed wall")
    u.write_text(text)
    # Static mesh: no dynamicMeshDict. The tutorial's `region fluid` and the
    # solid point-displacement function do not apply to a single-region case.
    functions = case / "system/functions"
    text = functions.read_text()
    text = re.sub(r"\n\s*pointDisp\s*\{[^}]*\}", "", text)
    text = re.sub(r"\n\s*region\s+fluid;", "", text)
    functions.write_text(text)

    driver.refine_mesh(case / "system/blockMeshDict", level)
    control = case / "system/controlDict"
    driver.replace_entry(control, "deltaT", f"{delta_t:.8g}")
    driver.replace_entry(control, "endTime", f"{end_time:.8g}")
    # Fields are written every 5 s, for restart or mapping; histories every step
    driver.replace_entry(control, "writeInterval", str(round(5.0 / delta_t)))
    driver.replace_entry(functions, "patches", "(plate cylinder)")
    driver.replace_entry(functions, "rhoInf", "1000")
    driver.replace_entry(case / "system/decomposeParDict", "numberOfSubdomains", str(cores))
    return case


def run(case: Path, cores: int) -> None:
    def execute(command: list[str], log: str) -> None:
        with (case / log).open("w") as handle:
            result = subprocess.run(command, cwd=case, stdout=handle,
                                    stderr=subprocess.STDOUT, text=True)
        if result.returncode:
            driver.fail(f"{' '.join(command)} failed; see {case / log}")

    execute(["blockMesh"], "log.blockMesh")
    if cores > 1:
        decompose = case / "system/decomposeParDict"
        driver.replace_entry(decompose, "method", "scotch")
        execute(["decomposePar", "-force"], "log.decomposePar")
        execute(["mpirun", "-np", str(cores), "solids4Foam", "-parallel"],
                "log.solids4Foam")
        execute(["reconstructPar", "-latestTime"], "log.reconstructPar")
    else:
        execute(["solids4Foam"], "log.solids4Foam")
    text = (case / "log.solids4Foam").read_text(errors="replace")
    if re.search(r"FOAM FATAL|FOAM aborting", text) or not re.search(
            r"^End\s*$", text, re.MULTILINE):
        driver.fail(f"{case.name} did not complete; see log.solids4Foam")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--level", type=int, required=True, choices=(1, 2, 4, 8))
    parser.add_argument("--delta-t", type=float, required=True)
    parser.add_argument("--end-time", type=float, default=8.0)
    parser.add_argument("--cores", type=int, default=1)
    parser.add_argument("--tag", default="")
    parser.add_argument("--setup-only", action="store_true")
    args = parser.parse_args()
    name = case_name(args.level, args.delta_t, args.tag)
    case = build_case(name, args.level, args.delta_t, args.end_time, args.cores)
    print(f"{name}: case written to {case}")
    if not args.setup_only:
        run(case, args.cores)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
