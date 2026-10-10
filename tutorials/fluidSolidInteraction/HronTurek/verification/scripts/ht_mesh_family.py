#!/usr/bin/env python3
"""Fluid blockMeshDict of a HronTurek mesh family with near-flag clustering.

The tutorial's fluid blockMeshDict is the base (level 1x of the standard
family). A clustered family keeps the block topology and changes, at every
level identically, only the divisions and gradings of the block rows that
touch the flag:

  row below the flag (y 0.12-0.19) and above it (y 0.21-0.28): --rows-n cells
      (base 8) with a geometric grading of total ratio --rows-ratio
      (last/first > 1), the small cells at the flag faces;
  the columns along the aft flag (x 0.28-0.6, --aft-n, base 30) and behind
      the tip (x 0.6-1.0, --wake-n, base 35), graded with total ratio
      --aft-ratio / --wake-ratio so that the small cells are at x = 0.6;
  the flag row behind the tip (y 0.19-0.21, x > 0.6): --tip-n cells (base 2)
      graded symmetrically towards both flag faces with total ratio
      --tip-ratio per half.

Then the in-plane divisions of every block are multiplied by --level, with
the gradings kept, so that the levels form a smooth family of ratio 2 (the
same fixed mapping refined uniformly), as in the standard study.

    python3 ht_mesh_family.py --base system/fluid/blockMeshDict --level 2 \
        --rows-n 16 --rows-ratio 4 --tip-n 6 --tip-ratio 3 --out blockMeshDict
"""
from __future__ import annotations

import argparse
import re
from pathlib import Path

# Blocks of the rows that touch the flag, as written in the tutorial, and the
# local direction (0 = x, 1 = y) that crosses the row, with the sense of the
# grading: +1 if the local axis runs towards the flag (small cells last,
# ratio < 1), -1 if it runs away from it (small cells first, ratio > 1)
BELOW = {  # y 0.12 -> 0.19, towards the lower flag face
    "7 8 14 13": (1, +1),
    "8 65 66 14": (1, +1),
    "65 9 15 66": (1, +1),
    "7 13 12 11": (0, +1),
}
ABOVE = {  # y 0.21 -> 0.28, away from the upper flag face
    "19 20 25 24": (1, -1),
    "20 67 68 25": (1, -1),
    "67 21 26 68": (1, -1),
    "17 18 19 24": (0, +1),  # local x runs 17->18 and 24->19, i.e. towards the flag
}
TIP = ("14 66 67 20", "66 15 21 67")  # local y 0.19 -> 0.21
# Columns along the aft flag (x 0.28-0.6, local x towards the tip) and behind
# the tip (x 0.6-1.0, local x away from the tip)
# Ring of blocks around the cylinder: the radial direction (base 11 cells,
# graded towards the cylinder; in the two flag-junction blocks the same
# direction also runs along the first 0.05 m of the flag faces)
RING = {"6 10 16 23": 0, "6 7 11 10": 1, "7 13 12 11": 1, "16 17 24 23": 1, "17 18 19 24": 1}
COL_AFT = ("2 3 8 7", "7 8 14 13", "19 20 25 24", "24 25 30 29")
COL_WAKE = ("3 64 65 8", "8 65 66 14", "14 66 67 20", "20 67 68 25", "25 68 69 30")


def fmt(g):
    return f"{g:.10g}" if not isinstance(g, str) else g


def build(base: str, level: int, rows_n: int, rows_ratio: float, tip_n: int, tip_ratio: float,
          aft_n: int = 30, aft_ratio: float = 1.0, wake_n: int = 35, wake_ratio: float = 1.0,
          ring_n: int = 11) -> str:
    pat = re.compile(r"hex\s+\(([\d ]+)\)\s+\((\d+) (\d+) (\d+)\)\s+simpleGrading\s+\(([^()]+)\)")

    def repl(m: re.Match) -> str:
        verts = " ".join(m.group(1).split()[:4])
        n = [int(m.group(2)), int(m.group(3)), int(m.group(4))]
        g: list = [float(v) for v in m.group(5).split()]
        for table in (BELOW, ABOVE):
            if verts in table:
                d, sense = table[verts]
                n[d] = rows_n
                g[d] = 1.0 / rows_ratio if sense > 0 else rows_ratio
        if verts in RING:
            n[RING[verts]] = ring_n
        if verts in COL_AFT:
            n[0] = aft_n
            g[0] = 1.0 / aft_ratio
        if verts in COL_WAKE:
            n[0] = wake_n
            g[0] = wake_ratio
        if verts in TIP:
            n[1] = tip_n
            g[1] = f"((0.5 0.5 {tip_ratio:.10g}) (0.5 0.5 {1.0 / tip_ratio:.10g}))"
        n = [n[0] * level, n[1] * level, n[2]]
        return (f"hex ({m.group(1)}) ({n[0]} {n[1]} {n[2]}) "
                f"simpleGrading ({' '.join(fmt(v) for v in g)})")

    text, count = pat.subn(repl, base)
    if count != 24:
        raise SystemExit(f"expected 24 blocks, found {count}")
    for v in list(BELOW) + list(ABOVE) + list(TIP):
        if f"hex ({v}" not in text and not re.search(rf"hex \({v} ", text):
            raise SystemExit(f"block {v} not found")
    return text


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--base", type=Path, required=True)
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
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args()
    a.out.write_text(build(a.base.read_text(), a.level, a.rows_n, a.rows_ratio, a.tip_n, a.tip_ratio,
                           a.aft_n, a.aft_ratio, a.wake_n, a.wake_ratio, a.ring_n))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
