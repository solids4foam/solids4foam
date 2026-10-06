#!/usr/bin/env python3
"""Run the Richter geometry with the beam held rigid."""

from __future__ import annotations

import argparse
import csv
import importlib.util
import math
import re
import shutil
import subprocess
import time
from pathlib import Path

SCRIPT = Path(__file__).resolve()
VERIFICATION = SCRIPT.parents[1]
MAIN_SCRIPT = SCRIPT.with_name("beam_in_cross_flow_verification.py")

spec = importlib.util.spec_from_file_location("beam_verify", MAIN_SCRIPT)
verify = importlib.util.module_from_spec(spec)
assert spec.loader is not None
spec.loader.exec_module(verify)


def copy_contents(source: Path, destination: Path) -> None:
    destination.mkdir(parents=True, exist_ok=True)
    for path in source.iterdir():
        target = destination / path.name
        if target.exists() or target.is_symlink():
            if target.is_dir() and not target.is_symlink():
                shutil.rmtree(target)
            else:
                target.unlink()
        if path.is_dir():
            shutil.copytree(path, target, symlinks=False)
        else:
            shutil.copy2(path.resolve(), target)


def prepare_case(factor: int, cores: int, delta_t: float, end_time: float) -> Path:
    case = verify.copy_case(f"richter_rigid_fluid_{factor}x")
    shutil.rmtree(case / "0/solid")
    copy_contents(case / "0/fluid", case / "0")
    shutil.rmtree(case / "0/fluid")
    copy_contents(case / "constant/fluid", case / "constant")
    copy_contents(case / "system/fluid", case / "system")

    verify.configure_graded_fluid_mesh(case / "system/blockMeshDict", factor)
    verify.replace_entry(case / "system/decomposeParDict", "numberOfSubdomains", str(cores))

    physics = case / "constant/physicsProperties"
    text = physics.read_text().replace("type    fluidSolidInteraction;", "type    fluid;")
    physics.write_text(text)

    velocity = case / "0/U"
    text = velocity.read_text()
    interface_match = re.search(r"interface\s*\{.*?\n\s*\}", text, re.DOTALL)
    if not interface_match:
        raise RuntimeError("could not find the fluid interface boundary condition")
    rigid = """interface
    {
        type            fixedValue;
        value           uniform (0 0 0);
    }"""
    text = text[:interface_match.start()] + rigid + text[interface_match.end():]
    velocity.write_text(text)

    dynamic_mesh = case / "constant/dynamicMeshDict"
    verify.replace_entry(dynamic_mesh, "dynamicFvMesh", "staticFvMesh")

    control = case / "system/controlDict"
    shutil.copy2((case / "system/controlDict.iqnils").resolve(), control)
    verify.replace_entry(control, "endTime", f"{end_time:.8g}")
    verify.replace_entry(control, "deltaT", f"{delta_t:.8g}")
    verify.replace_entry(control, "writeInterval", str(round(end_time / delta_t)))
    functions = case / "system/functions"
    functions.write_text(
        "functions\n{\n"
        "    forces\n    {\n"
        "        type forces;\n"
        "        libs (\"libforces.so\");\n"
        "        writeControl timeStep;\n"
        "        writeInterval 1;\n"
        "        patches (interface);\n"
        "        pName p;\n"
        "        UName U;\n"
        "        \"rhoName|rho\" rhoInf;\n"
        "        rhoInf 1000;\n"
        "        CofR (0.45 0.1 0);\n"
        "    }\n}\n"
    )
    return case


def run(case: Path, cores: int) -> float:
    start = time.monotonic()
    commands = [["blockMesh"]]
    if cores > 1:
        commands.extend(
            [
                ["decomposePar", "-force"],
                ["mpirun", "-np", str(cores), "solids4Foam", "-parallel"],
                ["reconstructPar"],
            ]
        )
    else:
        commands.append(["solids4Foam"])
    for index, command in enumerate(commands):
        log = case / f"log.rigid.{index}"
        with log.open("w") as handle:
            result = subprocess.run(
                command,
                cwd=case,
                stdout=handle,
                stderr=subprocess.STDOUT,
                text=True,
            )
        if result.returncode:
            raise RuntimeError(f"{' '.join(command)} failed; see {log}")
    return time.monotonic() - start


def fluid_cells(case: Path) -> int:
    return sum(
        math.prod(int(value) for value in match.group(1).split())
        for match in re.finditer(
            r"hex\s+\([^)]*\)\s+\(([^()]+)\)\s+simpleGrading",
            (case / "system/blockMeshDict").read_text(),
        )
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--levels", default="1,2,4,8")
    parser.add_argument("--cores", type=int, default=1)
    parser.add_argument("--delta-t", type=float, default=0.00625)
    parser.add_argument("--end-time", type=float, default=8.0)
    args = parser.parse_args()
    factors = [int(value) for value in args.levels.split(",")]
    rows = []
    for level, factor in enumerate(factors):
        case = prepare_case(factor, args.cores, args.delta_t, args.end_time)
        elapsed = run(case, args.cores)
        _, force = verify.force_at_or_before(verify.find_force(case), args.end_time)
        _, force_start = verify.force_at_or_before(
            verify.find_force(case), args.end_time - 1.0
        )
        rows.append(
            {
                "level": level,
                "refinement_factor": factor,
                "fluid_cells": fluid_cells(case),
                "near_body_cell_size": verify.mesh_characteristics("graded", factor)["near_body_cell_size"],
                "delta_t": args.delta_t,
                "ranks": args.cores,
                "wall_time_seconds": elapsed,
                "fx": force[0],
                "fy": force[1],
                "fz": force[2],
                "fx_last_second_relative_change": verify.relative_error(force[0], force_start[0]),
            }
        )
    for quantity in ("fx", "fy", "fz"):
        for index, row in enumerate(rows):
            row[f"{quantity}_relative_change"] = (
                "" if index == 0 else verify.relative_error(row[quantity], rows[index - 1][quantity])
            )
            row[f"{quantity}_observed_order"] = (
                "" if index < 2 else verify.observed_order(
                    rows[index - 2][quantity], rows[index - 1][quantity], row[quantity]
                ) or ""
            )
    output = VERIFICATION / "postProcessing/richter_rigid_fluid.csv"
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)
    print(output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
