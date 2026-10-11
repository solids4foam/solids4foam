#!/usr/bin/env python3
"""Fluid-side variants of the FSI3 trajectory replay (diagnostic study).

Builds a restart case with hron_turek_energy.py (same function objects and
analysis as the energy-balance study) and then applies one or more variants,
each a change to a single fluid treatment. Every entry was checked against the
source (solids4foam src and OpenFOAM v2512):

  pwall    p on `plate`: zeroGradient -> movingWallPressure (solids4foam):
           gradient = -n.a_wall, with a_wall the BDF2 wall acceleration stored
           by newMovingWallVelocity; pimpleFluid adds rAU*snGrad(p)*|Sf| to
           phiHbyA on such patches, so the wall flux stays the mesh flux.
  uwall    U on `plate`: newMovingWallVelocity -> movingWallVelocity (OpenFOAM
           v2512: Euler face-centre velocity, normal part from meshPhi)
  gradS4f  grad(U), grad(p): leastSquares -> leastSquaresS4f (solids4foam; needs the
           useBoundaryFaceValues lists that fluidModel creates for U and p only)
  gradS4fL as gradS4f, and the limiter gradient of div(phi,U) too
  mesh_inv diffusivity: quadratic inverseDistance -> inverseDistance
  mesh_exp diffusivity: quadratic inverseDistance -> exponential 1.0 inverseDistance
  gcl      mesh flux consistent with backward ddt: htVariants library,
           dynamicMotionSolverBDF2FvMesh (see platform/htVariants)
  euler    ddt backward -> Euler (not GCL-consistent comparison; for reference)

    python3 hron_turek_variants.py --source <replay> --restart 25 --end 33.2 \
        --name v2x_pwall --variants pwall --work ~/ht_variants/work
"""
from __future__ import annotations
import argparse, re, subprocess, sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
PLATE = re.compile(r"(?ms)^    plate\n    \{\n.*?^    \}\n")
PBC = ("    plate\n    {\n        type            movingWallPressure;\n"
       "        gradient        uniform 0;\n        value           uniform 0;\n    }\n")


def sub_once(text: str, old: str, new: str, what: str) -> str:
    if text.count(old) != 1:
        raise SystemExit(f"variant {what}: expected exactly one occurrence of {old!r}")
    return text.replace(old, new)


def tdirs(case: Path):
    d = sorted(case.glob("processor*"))
    return [x for x in d] or [case]


def apply(case: Path, restart: str, variant: str) -> None:
    for pd in tdirs(case):
        td = pd / restart
        if variant == "pwall":
            for name in ("p", "p_0"):
                f = td / name
                if f.is_file():
                    s = f.read_text()
                    assert len(PLATE.findall(s)) == 1, f
                    f.write_text(PLATE.sub(PBC, s))
        elif variant == "uwall":
            for name in ("U", "U_0"):
                f = td / name
                if f.is_file():
                    s = f.read_text()
                    blk = PLATE.search(s)
                    assert blk, f
                    if name == "U_0" and "newMovingWallVelocity" not in blk.group(0):
                        continue
                    nb = sub_once(blk.group(0), "newMovingWallVelocity", "movingWallVelocity", variant)
                    nb = re.sub(r"(?ms)^        oldFaceCentres.*?^        \)\n;\n|^        oldFaceCentres.*?\n\)\n;\n", "", nb)
                    s = s[:blk.start()] + nb + s[blk.end():]
                    f.write_text(s)
    sysd = case / "system"
    if variant in ("gradS4f", "gradS4fL"):
        f = sysd / "fvSchemes"
        # cellMotionU has no useBoundaryFaceValues list, so name the fields
        s = sub_once(f.read_text(), "    default leastSquares;",
                     "    default leastSquares;\n    grad(U) leastSquaresS4f;\n    grad(p) leastSquaresS4f;", variant)
        if variant == "gradS4fL":
            s = sub_once(s, "cellLimited leastSquares 1", "cellLimited leastSquaresS4f 1", variant)
        f.write_text(s)
    elif variant == "euler":
        f = sysd / "fvSchemes"
        f.write_text(sub_once(f.read_text(), "default backward;", "default Euler;", variant))
    elif variant in ("mesh_inv", "mesh_exp"):
        f = case / "constant" / "dynamicMeshDict"
        new = ("diffusivity inverseDistance 2(plate cylinder);" if variant == "mesh_inv"
               else "diffusivity exponential 1.0 inverseDistance 2(plate cylinder);")
        f.write_text(sub_once(f.read_text(), "diffusivity quadratic inverseDistance 2(plate cylinder);", new, variant))
    elif variant == "gcl":
        f = case / "constant" / "dynamicMeshDict"
        f.write_text(sub_once(f.read_text(), "dynamicFvMesh dynamicMotionSolverFvMesh;",
                              "dynamicFvMesh dynamicMotionSolverBDF2FvMesh;", variant))
        f = sysd / "controlDict"
        s = f.read_text()
        if "htVariants" not in s:
            s = sub_once(s, "application     solids4Foam;",
                         'application     solids4Foam;\n\nlibs            ( "libhtVariants.so" );', variant)
        f.write_text(s)
    elif variant not in ("base", "pwall", "uwall"):
        raise SystemExit(f"unknown variant {variant}")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--source", required=True, type=Path)
    ap.add_argument("--restart", required=True)
    ap.add_argument("--end", required=True)
    ap.add_argument("--name", required=True)
    ap.add_argument("--variants", default="base", help="comma-separated")
    ap.add_argument("--work", type=Path, required=True)
    a = ap.parse_args()
    subprocess.run([sys.executable, str(HERE / "hron_turek_energy.py"), "--source", str(a.source),
                    "--restart", a.restart, "--end", a.end, "--name", a.name, "--work", str(a.work)],
                   check=True)
    case = (a.work / a.name).resolve()
    for v in a.variants.split(","):
        apply(case, a.restart, v)
    with open(case / "energy_setup.txt", "a") as fh:
        fh.write(f"variants {a.variants}\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
