#!/usr/bin/env python3
"""Replay case on the smooth graded mesh family (ht_mesh_smooth.py).

Builds a ready-to-run, decomposed trajectory-replay case of the second-order
protocol (central convection, displacementLaplacian, movingWallPressure,
fastp; restart at 25 s from the developed rigid state, analysis window
ending at 30.2 s) from:

  --template   an existing replay case of the corrected set-up whose system/
               and constant/ are copied unchanged except for the mesh, the
               replay table and the run controls (here the 8x' ring+box case
               rb_8x_fastp_mlx: GAMG/DICGaussSeidel p solver, PIMPLE 3/3/1,
               tolerances 1e-6, energy and force function objects);
  --map-from   a case with the developed rigid state reconstructed at 25 s
               (U, p) and its constant/polyMesh: mapped with mapFields
               -consistent -mapMethod cellPointInterpolate;
  --point-template  the pointDisplacement file of a replay case at 25 s (the
               plate condition returns the replayed displacement; processor
               entries are replaced by a "procBoundary.*" entry);
  --table      the replay table of the coupled 2x motion, interpolated to the
               plate points of the new mesh by replay_table.py.

    python3 ht_smooth_case.py --template rb_8x_fastp_mlx --map-from src4x \\
        --point-template src4x/processor0/25/pointDisplacement --level 2 \\
        --delta-t 5e-4 --cores 16 --out work/sm_2x
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
import ht_mesh_smooth as ms  # noqa: E402

NONUNIFORM = re.compile(r"nonuniform\s+List<(\w+)>\s*\d*\s*\(.*?\)\s*;", re.S)
ZERO = {"scalar": "0", "vector": "(0 0 0)"}


def uniformise(text: str) -> str:
    return NONUNIFORM.sub(lambda m: f"uniform {ZERO[m.group(1)]};", text)


def strip_internal(text: str) -> str:
    """internalField -> uniform zero (the mapped values replace it)."""
    return re.sub(r"internalField\s+nonuniform\s+List<(\w+)>\s*\d*\s*\(.*?\)\s*;",
                  lambda m: f"internalField   uniform {ZERO[m.group(1)]};", text, count=1, flags=re.S)


def point_template(text: str) -> str:
    text = re.sub(r"\n    procBoundary\w+\s*\{[^{}]*\}", "", text)
    text = uniformise(text)
    i = text.rindex("}")  # end of boundaryField
    return text[:i] + "    \"procBoundary.*\"\n    {\n        type            processor;\n    }\n" + text[i:]


def run(cmd: list[str], case: Path, log: str) -> None:
    with (case / log).open("w") as fh:
        r = subprocess.run(cmd, cwd=case, stdout=fh, stderr=subprocess.STDOUT)
    if r.returncode:
        raise SystemExit(f"{' '.join(cmd)} failed, see {case / log}")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--template", type=Path, required=True)
    ap.add_argument("--map-from", type=Path, required=True)
    ap.add_argument("--point-template", type=Path, required=True)
    ap.add_argument("--base-mesh", type=Path, default=HERE.parent.parent / "system/fluid/blockMeshDict")
    ap.add_argument("--level", type=int, required=True)
    ap.add_argument("--delta-t", type=float, required=True)
    ap.add_argument("--cores", type=int, default=1)
    ap.add_argument("--end", type=float, default=30.2)
    ap.add_argument("--write-every", type=float, default=2.6, help="write interval (s of simulated time)")
    ap.add_argument("--table", type=Path, default=HERE.parent / "reference/energy_balance/plateMotion_2x.tab")
    ap.add_argument("--out", type=Path, required=True)
    for k, v in ms.DEFAULTS.items():
        ap.add_argument("--" + k.replace("_", "-"), type=float, default=v)
    a = ap.parse_args()
    params = {k: getattr(a, k) for k in ms.DEFAULTS}
    tpl, src, case = a.template.resolve(), a.map_from.resolve(), a.out.resolve()
    if case.exists():
        raise SystemExit(f"{case} exists")
    case.mkdir(parents=True)
    shutil.copytree(tpl / "system", case / "system")
    shutil.copytree(tpl / "constant", case / "constant",
                    ignore=shutil.ignore_patterns("polyMesh", "plateMotion.tab"))
    (case / "system/blockMeshDict").write_text(ms.build(a.base_mesh.read_text(), a.level, params))
    cd = case / "system/controlDict"
    s = cd.read_text()
    steps = round(a.write_every / a.delta_t)
    for k, v in (("deltaT", repr(a.delta_t)), ("startFrom", "latestTime"), ("startTime", "25"),
                 ("endTime", f"{a.end:g}"), ("writeControl", "timeStep"), ("writeInterval", str(steps)),
                 ("purgeWrite", "2")):
        s, n = re.subn(rf"(?m)^{k}\s+.*$", f"{k:15s} {v};", s)
        if n != 1:
            raise SystemExit(f"controlDict: {k} not found once")
    cd.write_text(s)
    run(["blockMesh"], case, "log.blockMesh")
    run(["checkMesh", "-noFunctionObjects"], case, "log.checkMesh")
    t25 = case / "25"
    t25.mkdir()
    for f in ("U", "p"):
        (t25 / f).write_text(strip_internal(uniformise((src / "25" / f).read_text())))
    run(["mapFields", str(src), "-consistent", "-sourceTime", "25", "-mapMethod", "cellPointInterpolate",
         "-noFunctionObjects"], case, "log.mapFields")
    for f in ("U_0", "p_0", "phi", "phi_0"):
        if (t25 / f).exists():
            (t25 / f).unlink()
    subprocess.run([sys.executable, str(HERE / "replay_table.py"), "--source", str(a.table),
                    "--mesh", str(case / "constant/polyMesh"), "--out", str(case / "constant/plateMotion.tab")],
                   check=True)
    dp = case / "system/decomposeParDict"
    s = re.sub(r"(?m)^numberOfSubdomains\s+\d+;", f"numberOfSubdomains {a.cores};", dp.read_text())
    s = re.sub(r"(?m)^method\s+\w+;", "method          scotch;", s)
    dp.write_text(s)
    if a.cores > 1:  # U and p only: the coded plate condition is compiled by the solver on the run host
        run(["decomposePar", "-time", "25", "-noFunctionObjects"], case, "log.decomposePar")
    (t25 / "pointDisplacement").write_text(point_template(a.point_template.read_text()))
    if a.cores > 1:
        for p in sorted(case.glob("processor*")):
            shutil.copy(t25 / "pointDisplacement", p / "25" / "pointDisplacement")
    setup = (tpl / "energy_setup.txt").read_text()
    setup = re.sub(r"(?m)^source .*$", f"source {src} (rigid 25 s state, mapped); mesh ht_mesh_smooth.py level {a.level}", setup)
    setup = re.sub(r"(?m)^nprocs .*$", f"nprocs {a.cores}", setup)
    setup = re.sub(r"(?m)^end .*$", f"end {a.end:g}", setup)
    (case / "energy_setup.txt").write_text(setup)
    (case / "mesh_family.txt").write_text(
        "generator ht_mesh_smooth.py\n" + "".join(f"{k} {v:g}\n" for k, v in params.items())
        + f"level {a.level}\ndelta_t {a.delta_t:g}\ncores {a.cores}\nmap_from {src}\ntemplate {tpl}\n")
    print(case, "ready")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
