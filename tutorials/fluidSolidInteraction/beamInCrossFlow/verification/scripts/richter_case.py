#!/usr/bin/env python3
"""Build one Richter beamInCrossFlow case for the verification campaign.

The case is a copy of the tutorial with every setting that the campaign varies
written explicitly, so that a run directory documents itself.  The script only
prepares the case; meshing, decomposition and the solver are run by the batch
scripts (see ``../campaign``), which keeps one code path for MeluXina and
workstation runs.

Mesh family
    ``--factor f`` multiplies every block cell count of the graded family
    (fixed total expansion ratios, see ``beam_in_cross_flow_verification``).
    ``--thickness-multiplier m`` additionally multiplies the through-thickness
    (x) count of the solid and of the three fluid blocks above and beside the
    beam, so that the fluid-solid interface stays conformal (direct mapping).
    The solid then has ``4 m f`` cells through the thickness and ``8 f`` cells
    in y and z.

Examples
    richter_case.py work/L1_m4 --factor 2 --thickness-multiplier 4 \\
        --dt 0.00625 --end-time 8 --ranks 32
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import re
import shutil
import subprocess
from pathlib import Path

SCRIPT = Path(__file__).resolve()
MAIN_SCRIPT = SCRIPT.with_name("beam_in_cross_flow_verification.py")
spec = importlib.util.spec_from_file_location("beam_verify", MAIN_SCRIPT)
verify = importlib.util.module_from_spec(spec)
assert spec.loader is not None
spec.loader.exec_module(verify)

TUTORIAL = verify.TUTORIAL

# Graded fluid family: (counts, total expansion ratios) per block at f = 1.
# Blocks 4, 5 and 6 span the beam thickness (0.4 < x < 0.5).
FLUID_BLOCKS = [
    ((12, 8, 8), (0.125, 1, 1)),
    ((12, 8, 6), (0.125, 1, 1 / 6)),
    ((12, 6, 8), (0.125, 6, 1)),
    ((12, 6, 6), (0.125, 6, 1 / 6)),
    ((4, 8, 6), (1, 1, 1 / 6)),
    ((4, 6, 6), (1, 6, 1 / 6)),
    ((4, 6, 8), (1, 6, 1)),
    ((24, 8, 8), (8, 1, 1)),
    ((24, 8, 6), (8, 1, 1 / 6)),
    ((24, 6, 8), (8, 6, 1)),
    ((24, 6, 6), (8, 6, 1 / 6)),
]
THICKNESS_BLOCKS = (4, 5, 6)


def fluid_counts(factor: int, multiplier: int) -> list[tuple[int, int, int]]:
    counts = []
    for index, (base, _) in enumerate(FLUID_BLOCKS):
        nx, ny, nz = (value * factor for value in base)
        if index in THICKNESS_BLOCKS:
            nx *= multiplier
        counts.append((nx, ny, nz))
    return counts


def solid_counts(factor: int, multiplier: int) -> tuple[int, int, int]:
    return (4 * factor * multiplier, 8 * factor, 8 * factor)


def write_fluid_mesh(path: Path, factor: int, multiplier: int) -> None:
    counts = fluid_counts(factor, multiplier)
    index = 0

    def replace(match: re.Match[str]) -> str:
        nonlocal index
        cells = " ".join(str(value) for value in counts[index])
        grading = " ".join(f"{value:.12g}" for value in FLUID_BLOCKS[index][1])
        index += 1
        return f"{match.group(1)}({cells}) simpleGrading ({grading})"

    text = re.sub(
        r"(hex\s+\([^)]*\)\s+)\([^()]+\)\s+simpleGrading\s+\([^)]*\)",
        replace,
        path.read_text(),
    )
    if index != len(FLUID_BLOCKS):
        raise RuntimeError(f"expected {len(FLUID_BLOCKS)} fluid blocks in {path}")
    path.write_text(text)


def write_solid_mesh(path: Path, factor: int, multiplier: int) -> None:
    cells = " ".join(str(value) for value in solid_counts(factor, multiplier))
    text, count = re.subn(
        r"(hex\s+\([^)]*\)\s+)\([^()]+\)\s+simpleGrading",
        rf"\g<1>({cells}) simpleGrading",
        path.read_text(),
    )
    if count != 1:
        raise RuntimeError(f"expected one solid block in {path}")
    path.write_text(text)


def set_entry(path: Path, key: str, value: str) -> None:
    verify.replace_entry(path, key, value)


def set_in_block(path: Path, block: str, key: str, value: str) -> None:
    """Set ``key`` inside the first dictionary named ``block``."""
    text = path.read_text()
    match = re.search(rf"^(\s*){re.escape(block)}\s*\n\s*\{{", text, re.MULTILINE)
    if not match:
        raise RuntimeError(f"no '{block}' dictionary in {path}")
    start = match.end()
    depth, end = 1, start
    while depth and end < len(text):
        depth += {"{": 1, "}": -1}.get(text[end], 0)
        end += 1
    body = text[start:end - 1]
    pattern = rf"^(\s*{re.escape(key)}\s+)[^;]+;"
    body, count = re.subn(pattern, rf"\g<1>{value};", body, count=1, flags=re.MULTILINE)
    if count == 0:
        body = f"\n        {key} {value};" + body
    path.write_text(text[:start] + body + text[end - 1:])


def set_petsc_option(path: Path, key: str, value: str) -> None:
    text = path.read_text()
    text, count = re.subn(
        rf'^(\s*{re.escape(key)}\s+)"[^"]*";',
        rf'\g<1>"{value}";',
        text,
        flags=re.MULTILINE,
    )
    if count != 1:
        raise RuntimeError(f"expected one PETSc option {key} in {path}")
    path.write_text(text)


def set_scheme(path: Path, section: str, scheme: str) -> None:
    text, count = re.subn(
        rf"({section}\s*\{{\s*default\s+)[^;]+;",
        rf"\g<1>{scheme};",
        path.read_text(),
        count=1,
        flags=re.DOTALL,
    )
    if count != 1:
        raise RuntimeError(f"expected {section}/default in {path}")
    path.write_text(text)


def materialise(case: Path, link: str, target: str) -> Path:
    """Replace the tutorial's selector symlink by a real copy of target."""
    path = case / link
    source = path.parent / target
    if path.is_symlink() or path.exists():
        path.unlink()
    shutil.copy2(source, path)
    return path


def build(args: argparse.Namespace) -> Path:
    case = Path(args.case).resolve()
    if case.exists():
        if not args.force:
            raise RuntimeError(f"{case} exists; use --force to replace it")
        shutil.rmtree(case)
    ignore = shutil.ignore_patterns(
        "verification", "postProcessing", "processor*", "log.*", "*.pdf",
        "images",
    )
    shutil.copytree(TUTORIAL, case, symlinks=True, ignore=ignore)

    # Fixed selector choices: IQN-ILS, Dirichlet-Neumann interface
    fsi = materialise(case, "constant/fsiProperties", "fsiProperties.iqnils")
    solid_props = materialise(
        case,
        "constant/solid/solidProperties",
        "solidProperties.iqnils" + (".highOrder" if args.high_order_solid else ""),
    )
    control = materialise(case, "system/controlDict", "controlDict.iqnils")
    fluid_solution = materialise(case, "system/fluid/fvSolution", "fvSolution.iqnils")
    materialise(case, "0/fluid/U", "U.dirichletNeumann")
    materialise(case, "0/fluid/p", "p.dirichletNeumann")
    for name in ("fsiProperties.aitken", "fsiProperties.robin", "fsiProperties.iqnils"):
        (case / "constant" / name).unlink()
    for path in (case / "constant/solid").glob("solidProperties.*"):
        path.unlink()
    for path in (case / "system").glob("controlDict.*"):
        path.unlink()
    for path in (case / "system/fluid").glob("fvSolution.*"):
        path.unlink()
    for path in (case / "0/fluid").glob("[Up].*"):
        path.unlink()

    # Richter definition (inlet peak, ramp, material)
    u_file = case / "0/fluid/U"
    set_entry(u_file, "maxVelocity", f"{args.inlet_peak:.12g}")
    set_entry(u_file, "timeVaryingEndTime", f"{args.ramp:.12g}")
    set_entry(u_file, "timeAtMaxVelocity", f"{args.ramp:.12g}")
    set_entry(
        case / "constant/solid/mechanicalProperties", "E", "E [1 -1 -2 0 0 0 0] 1.4e6"
    )

    # Meshes
    write_fluid_mesh(case / "system/fluid/blockMeshDict", args.factor,
                     args.thickness_multiplier)
    write_solid_mesh(case / "system/solid/blockMeshDict", args.factor,
                     args.thickness_multiplier)

    # Time control and output
    set_entry(control, "deltaT", f"{args.dt:.12g}")
    set_entry(control, "endTime", f"{args.end_time:.12g}")
    set_entry(control, "writeInterval", str(args.write_interval))
    set_entry(control, "purgeWrite", str(args.purge_write))
    set_entry(control, "writeFormat", "binary" if args.binary else "ascii")
    set_entry(control, "writePrecision", str(args.write_precision))
    set_entry(control, "timePrecision", "10")
    set_entry(case / "system/functions", "writeInterval", "1")
    for scheme_file in (case / "system/fluid/fvSchemes",):
        set_scheme(scheme_file, "ddtSchemes", args.fluid_ddt)
    set_scheme(case / "system/solid/fvSchemes", "d2dt2Schemes", args.solid_d2dt2)
    set_scheme(case / "system/solid/fvSchemes", "ddtSchemes", args.solid_ddt)

    # Coupling
    set_entry(fsi, "outerCorrTolerance", f"{args.coupling_tol:.12g}")
    set_entry(fsi, "nOuterCorr", str(args.n_outer_corr))
    set_entry(fsi, "requireAllResidualMeasures", "yes" if args.strict else "no")
    set_entry(fsi, "relaxationFactor", f"{args.relaxation:.12g}")
    text = fsi.read_text()
    if "allowUnconvergedCoupling" not in text:
        text = text.replace(
            "    requireAllResidualMeasures",
            "    allowUnconvergedCoupling "
            + ("yes" if args.allow_unconverged else "no")
            + ";\n\n    requireAllResidualMeasures",
            1,
        )
        fsi.write_text(text)

    # Fluid algorithm and linear solvers
    set_in_block(fluid_solution, "PIMPLE", "nOuterCorrectors", str(args.pimple_outer))
    text = fluid_solution.read_text()
    for field, tolerance in (
        ('"p|pFinal"', args.p_tol), ('"U|UFinal"', args.u_tol),
        ("cellMotionU", args.motion_tol),
    ):
        block = re.search(rf'{re.escape(field)}\s*\{{[^}}]*\}}', text)
        if not block:
            raise RuntimeError(f"no {field} solver in {fluid_solution}")
        new = re.sub(r"(tolerance\s+)[^;]+;", rf"\g<1>{tolerance:.12g};", block.group(0))
        text = text[:block.start()] + new + text[block.end():]
    fluid_solution.write_text(text)
    if args.pimple_residual_control is not None:
        value = f"{args.pimple_residual_control:.12g}"
        text = fluid_solution.read_text()
        text = re.sub(r"(relTol\s+)1e-5;", rf"\g<1>{value};", text)
        text = re.sub(r"(tolerance\s+)1e-5;", rf"\g<1>{value};", text)
        fluid_solution.write_text(text)
    solid_solution = case / "system/solid/fvSolution"
    set_petsc_option(solid_solution, "snes_rtol", f"{args.snes_rtol:.12g}")
    set_petsc_option(solid_solution, "snes_stol", f"{args.snes_stol:.12g}")
    set_entry(solid_props, "predictor", "yes")
    text = solid_props.read_text()
    if "restart" not in text:
        # Total-displacement St Venant-Kirchhoff needs no extra restart
        # fields; continuations happen only near the steady state (the
        # missing old-old level of issue #555 then implies the correct, zero,
        # old velocity).
        text = text.replace("    predictor yes;", "    predictor yes;\n    restart no;", 1)
        solid_props.write_text(text)

    for dictionary in (
        case / "system/decomposeParDict",
        case / "system/fluid/decomposeParDict",
        case / "system/solid/decomposeParDict",
    ):
        if dictionary.is_file():
            set_entry(dictionary, "numberOfSubdomains", str(args.ranks))

    fluid = fluid_counts(args.factor, args.thickness_multiplier)
    solid = solid_counts(args.factor, args.thickness_multiplier)
    metadata = dict(vars(args))
    metadata.update({
        "fluid_cells": sum(a * b * c for a, b, c in fluid),
        "solid_cells": solid[0] * solid[1] * solid[2],
        "solid_counts_xyz": solid,
        "case_commit": subprocess.run(
            ["git", "-C", str(SCRIPT.parent), "rev-parse", "HEAD"],
            capture_output=True, text=True,
        ).stdout.strip(),
    })
    if args.init_fluid:
        initialise_fluid(case, Path(args.init_fluid))
    if args.rigid:
        make_rigid(case)
    (case / "campaign_case.json").write_text(json.dumps(metadata, indent=2) + "\n")
    return case


def internal_field(path: Path) -> str:
    """Return the internalField entry (keyword to semicolon) of an ascii file."""
    text = path.read_text()
    start = text.index("internalField")
    depth, index = 0, start
    while True:
        char = text[index]
        if char == "(":
            depth += 1
        elif char == ")":
            depth -= 1
        elif char == ";" and depth == 0:
            return text[start:index + 1]
        index += 1


def initialise_fluid(case: Path, source: Path) -> None:
    """Start the fluid from a converged rigid-beam flow on the same mesh.

    Only the internal fields of U and p are taken; the boundary conditions of
    the coupled case are kept.  The inlet is then applied at full strength
    from t = 0 (no ramp), because the flow is already developed.
    """
    for name in ("U", "p"):
        target = case / "0/fluid" / name
        text = target.read_text()
        old = internal_field(target)
        text = text.replace(old, internal_field(source / name), 1)
        target.write_text(text)
    u_file = case / "0/fluid/U"
    set_entry(u_file, "timeVaryingEndTime", "0")
    set_entry(u_file, "timeAtMaxVelocity", "1")


def move_contents(source: Path, destination: Path) -> None:
    for path in source.iterdir():
        target = destination / path.name
        if target.is_symlink() or target.is_file():
            target.unlink()
        elif target.is_dir():
            shutil.rmtree(target)
        shutil.move(str(path), str(target))
    source.rmdir()


def make_rigid(case: Path) -> None:
    """Turn the FSI case into the fluid-only case with a rigid beam.

    The fluid mesh, schemes, solvers, inlet and time control are those of the
    coupled case; the interface is a fixed no-slip wall and the mesh is static.
    """
    shutil.rmtree(case / "0/solid")
    shutil.rmtree(case / "constant/solid")
    shutil.rmtree(case / "system/solid")
    for name in ("decomposeParDict",):
        (case / "system" / name).unlink()
    move_contents(case / "0/fluid", case / "0")
    move_contents(case / "constant/fluid", case / "constant")
    move_contents(case / "system/fluid", case / "system")
    (case / "0/pointMotionU").unlink()
    (case / "constant/fsiProperties").unlink()
    physics = case / "constant/physicsProperties"
    physics.write_text(
        physics.read_text().replace("type    fluidSolidInteraction;", "type    fluid;")
    )
    velocity = case / "0/U"
    text = velocity.read_text()
    match = re.search(r"interface\s*\{.*?\n\s*\}", text, re.DOTALL)
    if not match:
        raise RuntimeError("no interface condition in 0/U")
    text = (
        text[:match.start()]
        + "interface\n    {\n        type            noSlip;\n    }"
        + text[match.end():]
    )
    velocity.write_text(text)
    (case / "constant/dynamicMeshDict").write_text(
        "FoamFile\n{\n    version 2.0;\n    format ascii;\n    class dictionary;\n"
        "    object dynamicMeshDict;\n}\n\ndynamicFvMesh staticFvMesh;\n"
    )
    functions = case / "system/functions"
    text = functions.read_text()
    text = re.sub(r"\n\s*displacement\s*\{[^}]*\}", "", text)
    text = re.sub(r"\n\s*region\s+fluid;", "", text)
    functions.write_text(text)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("case")
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--factor", type=int, required=True)
    parser.add_argument("--thickness-multiplier", type=int, default=1)
    parser.add_argument("--high-order-solid", action="store_true")
    parser.add_argument("--inlet-peak", type=float, default=0.3)
    parser.add_argument("--ramp", type=float, default=1.0)
    parser.add_argument("--dt", type=float, default=0.00625)
    parser.add_argument("--end-time", type=float, default=8.0)
    parser.add_argument("--write-interval", type=int, default=80,
                        help="time steps between writes")
    parser.add_argument("--purge-write", type=int, default=2)
    parser.add_argument("--binary", action="store_true")
    parser.add_argument("--write-precision", type=int, default=12)
    parser.add_argument("--fluid-ddt", default="backward")
    parser.add_argument("--solid-d2dt2", default="backward")
    parser.add_argument("--solid-ddt", default="backward")
    parser.add_argument("--coupling-tol", type=float, default=1e-6)
    parser.add_argument("--strict", action=argparse.BooleanOptionalAction, default=True,
                        help="requireAllResidualMeasures (max of r1, r2)")
    parser.add_argument("--n-outer-corr", type=int, default=200)
    parser.add_argument("--allow-unconverged", action=argparse.BooleanOptionalAction,
                        default=False)
    parser.add_argument("--relaxation", type=float, default=0.1)
    parser.add_argument("--pimple-outer", type=int, default=1)
    parser.add_argument("--pimple-residual-control", type=float)
    parser.add_argument("--p-tol", type=float, default=1e-6)
    parser.add_argument("--u-tol", type=float, default=1e-6)
    parser.add_argument("--motion-tol", type=float, default=1e-6,
                        help="cellMotionU (mesh-motion velocity) solver tolerance")
    parser.add_argument("--snes-rtol", type=float, default=1e-6)
    parser.add_argument("--snes-stol", type=float, default=1e-6)
    parser.add_argument("--ranks", type=int, default=1)
    parser.add_argument("--init-fluid",
                        help="time directory of a reconstructed rigid-beam run "
                        "on the same fluid mesh; its U and p start the fluid")
    parser.add_argument("--rigid", action="store_true",
                        help="fluid-only case with the beam held rigid")
    args = parser.parse_args()
    case = build(args)
    print(case)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
