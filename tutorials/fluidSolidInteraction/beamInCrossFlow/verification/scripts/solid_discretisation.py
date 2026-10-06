#!/usr/bin/env python3
"""Run the Richter beam's static solid-discretisation diagnostic."""

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


def prepare_case(factor: int, cores: int, traction: float) -> Path:
    case = verify.copy_case(f"richter_solid_static_{factor}x")
    shutil.rmtree(case / "0/fluid")
    copy_contents(case / "0/solid", case / "0")
    shutil.rmtree(case / "0/solid")
    copy_contents(case / "constant/solid", case / "constant")
    copy_contents(case / "system/solid", case / "system")

    verify.refine_mesh(case / "system/blockMeshDict", factor)
    verify.replace_entry(case / "system/decomposeParDict", "numberOfSubdomains", str(cores))

    physics = case / "constant/physicsProperties"
    text = physics.read_text().replace("type    fluidSolidInteraction;", "type    solid;")
    physics.write_text(text)

    displacement = case / "0/D"
    text = displacement.read_text()
    text = re.sub(
        r"(interface\s*\{.*?traction\s+uniform\s*)\([^)]*\)",
        rf"\g<1>({traction:.12g} 0 0)",
        text,
        count=1,
        flags=re.DOTALL,
    )
    displacement.write_text(text)

    schemes = case / "system/fvSchemes"
    text = schemes.read_text()
    text = re.sub(
        r"(d2dt2Schemes\s*\{\s*default\s+)[^;]+;",
        r"\g<1>steadyState;",
        text,
        count=1,
        flags=re.DOTALL,
    )
    text = re.sub(
        r"(ddtSchemes\s*\{\s*default\s+)[^;]+;",
        r"\g<1>steadyState;",
        text,
        count=1,
        flags=re.DOTALL,
    )
    schemes.write_text(text)

    control = case / "system/controlDict"
    shutil.copy2((case / "system/controlDict.iqnils").resolve(), control)
    verify.replace_entry(control, "endTime", "1")
    verify.replace_entry(control, "deltaT", "1")
    verify.replace_entry(control, "writeInterval", "1")
    functions = case / "system/functions"
    functions.write_text(
        "functions\n{\n"
        "    displacement\n    {\n"
        "        type solidPointDisplacement;\n"
        "        point (0.45 0.15 -0.15);\n"
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
        log = case / f"log.static.{index}"
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


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--levels", default="1,2,4,8")
    parser.add_argument("--cores", type=int, default=1)
    parser.add_argument("--traction", type=float, default=100.0)
    args = parser.parse_args()
    factors = [int(value) for value in args.levels.split(",")]
    rows = []
    for level, factor in enumerate(factors):
        case = prepare_case(factor, args.cores, args.traction)
        elapsed = run(case, args.cores)
        _, displacement = verify.vector_at_or_before(
            verify.find_displacement(case), 1.0
        )
        cells = sum(
            math.prod(int(value) for value in match.group(1).split())
            for match in re.finditer(
                r"hex\s+\([^)]*\)\s+\(([^()]+)\)\s+simpleGrading",
                (case / "system/blockMeshDict").read_text(),
            )
        )
        rows.append(
            {
                "level": level,
                "refinement_factor": factor,
                "solid_cells": cells,
                "cells_across_thickness": 4 * factor,
                "interface_cells_y": 8 * factor,
                "interface_cells_z": 8 * factor,
                "traction_x_pa": args.traction,
                "ranks": args.cores,
                "wall_time_seconds": elapsed,
                "ux_A": displacement[0],
                "uy_A": displacement[1],
                "uz_A": displacement[2],
            }
        )
    for index, row in enumerate(rows):
        row["ux_relative_change"] = (
            "" if index == 0 else verify.relative_error(row["ux_A"], rows[index - 1]["ux_A"])
        )
        row["ux_observed_order"] = (
            "" if index < 2 else verify.observed_order(
                rows[index - 2]["ux_A"], rows[index - 1]["ux_A"], row["ux_A"]
            ) or ""
        )
    output = VERIFICATION / "postProcessing/richter_solid_discretisation.csv"
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)
    print(output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
