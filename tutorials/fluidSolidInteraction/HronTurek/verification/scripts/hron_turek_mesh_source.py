#!/usr/bin/env python3
"""Replay source case (rigid developed state at 25 s) on a new fluid mesh.

Builds a case with the same constant/system as an existing replay source of a
standard level, the fluid mesh of ht_mesh_family.py (any level and
clustering), the developed rigid CFD3 flow mapped from that source at 25 s
(mapFields -consistent, cellPointInterpolate), and the replay table of the coupled
2x motion interpolated to the new plate points (replay_table.py). The result
can be used as --source of hron_turek_second_order.py with --restart 25. The
mapped state is not a converged rigid state of the new mesh; the replay starts
with the 1 s motion ramp from it, and the analysis windows are taken after
the flow has become periodic (check with qin_slide.py).

    python3 hron_turek_mesh_source.py --map-from ~/ht_variants/src/ale_2x_dt0.0005_replay \
        --level 2 --rows-n 16 --rows-ratio 4 --tip-n 6 --tip-ratio 3 --cores 8 \
        --delta-t 0.0005 --out ~/ht_2nd/src/cl_2x
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

NONUNIFORM = re.compile(r"nonuniform\s+List<(\w+)>\s*\d*\s*\(.*?\)\s*;", re.S)
ZERO = {"scalar": "0", "vector": "(0 0 0)", "symmTensor": "(0 0 0 0 0 0)", "tensor": "(0 0 0 0 0 0 0 0 0)"}


def uniformise(text: str) -> str:
    return NONUNIFORM.sub(lambda m: f"uniform {ZERO[m.group(1)]};", text)


def run(cmd: list[str], case: Path, log: str) -> None:
    with (case / log).open("w") as fh:
        r = subprocess.run(cmd, cwd=case, stdout=fh, stderr=subprocess.STDOUT)
    if r.returncode:
        raise SystemExit(f"{' '.join(cmd)} failed, see {case / log}")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--map-from", type=Path, required=True, help="replay source case with a 25 s state")
    ap.add_argument("--base-mesh", type=Path, help="base blockMeshDict (default: <map-from>/system, level-1 tutorial)")
    ap.add_argument("--level", type=int, required=True)
    ap.add_argument("--rows-n", type=int, default=8)
    ap.add_argument("--rows-ratio", type=float, default=1.0)
    ap.add_argument("--tip-n", type=int, default=2)
    ap.add_argument("--tip-ratio", type=float, default=1.0)
    ap.add_argument("--aft-n", type=int, default=30)
    ap.add_argument("--aft-ratio", type=float, default=1.0)
    ap.add_argument("--wake-n", type=int, default=35)
    ap.add_argument("--wake-ratio", type=float, default=1.0)
    ap.add_argument("--ring-n", type=int, default=11)
    ap.add_argument("--delta-t", type=float, required=True)
    ap.add_argument("--cores", type=int, default=1)
    ap.add_argument("--table", type=Path, default=Path.home() / "ht_2nd/plateMotion_2x.tab")
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args()
    src, case = a.map_from.resolve(), a.out.resolve()
    if case.exists():
        shutil.rmtree(case)
    case.mkdir(parents=True)
    shutil.copytree(src / "system", case / "system")
    shutil.copytree(src / "constant", case / "constant", ignore=shutil.ignore_patterns("polyMesh"))
    base = a.base_mesh or (HERE.parent.parent / "system/fluid/blockMeshDict")
    if not base.is_file():
        base = HERE.parent / "blockMeshDict.base"
    (case / "system/blockMeshDict").write_text(
        mf.build(base.read_text(), a.level, a.rows_n, a.rows_ratio, a.tip_n, a.tip_ratio,
                 a.aft_n, a.aft_ratio, a.wake_n, a.wake_ratio, a.ring_n))
    cd = case / "system/controlDict"
    s = cd.read_text()
    s = re.sub(r"(?m)^deltaT .*$", f"deltaT          {a.delta_t:g};", s)
    s = re.sub(r"(?m)^startFrom .*$", "startFrom       startTime;", s)
    s = re.sub(r"(?m)^startTime .*$", "startTime       25;", s)
    cd.write_text(s)
    run(["blockMesh"], case, "log.blockMesh")
    run(["checkMesh", "-noFunctionObjects"], case, "log.checkMesh")
    # Field templates at 25 s: the source's conditions with uniform values
    srct = src / "25" if (src / "25").is_dir() else src / "processor0" / "25"
    (case / "25").mkdir()
    for f in ("U", "p", "pointMotionU"):
        (case / "25" / f).write_text(uniformise((srct / f).read_text()))
    # Map the developed rigid state (reconstructed source needed)
    if not (src / "25" / "U").is_file():
        raise SystemExit(f"{src}/25 is not reconstructed")
    run(["mapFields", str(src), "-consistent", "-sourceTime", "25", "-mapMethod", "cellPointInterpolate",
         "-noFunctionObjects"], case, "log.mapFields")
    # mapFields leaves the time-0 style names; drop any old-time fields it may have written
    for f in ("U_0", "p_0", "phi", "phi_0"):
        p = case / "25" / f
        if p.exists():
            p.unlink()
    # mapFields moves the (uniform) point field it cannot map aside; restore it
    for p in (case / "25").glob("*.unmapped"):
        p.rename(p.with_suffix(""))
    subprocess.run([sys.executable, str(HERE / "replay_table.py"), "--source", str(a.table),
                    "--mesh", str(case / "constant/polyMesh"), "--out", str(case / "constant/plateMotion.tab")],
                   check=True)
    dp = case / "system/decomposeParDict"
    s = re.sub(r"(?m)^numberOfSubdomains\s+\d+;", f"numberOfSubdomains {a.cores};", dp.read_text())
    s = re.sub(r"(?m)^method\s+\w+;", "method          scotch;", s)
    dp.write_text(s)
    if a.cores > 1:
        run(["decomposePar", "-time", "25", "-noFunctionObjects"], case, "log.decomposePar")
        shutil.rmtree(case / "25")
    (case / "mesh_family.txt").write_text(
        f"map_from {src}\nlevel {a.level}\nrows_n {a.rows_n}\nrows_ratio {a.rows_ratio}\n"
        f"tip_n {a.tip_n}\ntip_ratio {a.tip_ratio}\naft_n {a.aft_n}\naft_ratio {a.aft_ratio}\n"
        f"wake_n {a.wake_n}\nwake_ratio {a.wake_ratio}\nring_n {a.ring_n}\ndelta_t {a.delta_t}\ncores {a.cores}\n")
    print(case, "ready")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
