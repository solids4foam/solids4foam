#!/usr/bin/env python3
"""Replay QoIs of the smooth family against (1/L)^2, with the earlier families.

    python3 smooth_family_plot.py smooth_family.json --out smooth_family.png
"""
from __future__ import annotations

import argparse
import cmath
import csv
import json
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

REF = Path(__file__).resolve().parent.parent / "reference"
EARLIER = {  # family -> (runs file, case per level)
    "standard": (REF / "second_order/runs.csv",
                 {1: "s3_1x_linear_dispLap_pwall", 2: "s3_2x_linear_dispLap_pwall",
                  3: "s4_3x_linear_dispLap_pwall", 4: "s4_4x_linear_dispLap_pwall"}),
    "ring+box": (REF / "production/replay_8x_runs.csv",
                 {1: "rb_1x_fastp", 2: "rb_2x_fastp", 3: "rb_3x_fastp", 4: "rb_4x_fastp", 8: "rb_8x_fastp_mlx"}),
}


def earlier():
    out = {}
    for fam, (path, cases) in EARLIER.items():
        rows = {}
        for f in (path, REF / "second_order/runs.csv"):
            if f.is_file():
                for r in csv.DictReader(f.open()):
                    rows.setdefault(r["case"], r)
        vals = {}
        for lv, c in cases.items():
            r = rows.get(c)
            if r and r.get("lift1_amp"):
                vals[lv] = {k: float(r[k]) for k in ("Q_in", "lift1_amp", "lift1_re", "lift1_im")}
        out[fam] = vals
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("json", type=Path)
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args()
    d = json.loads(a.json.read_text())["runs"]
    fams = earlier()
    fams["smooth"] = {int(t[0]): d[t] for t in ("1x", "2x", "3x", "4x") if t in d}
    fig, ax = plt.subplots(1, 3, figsize=(13, 4))
    marks = {"standard": "s--", "ring+box": "^--", "smooth": "o-"}
    for fam, vals in fams.items():
        lv = sorted(vals)
        x = [1.0 / v ** 2 for v in lv]
        ax[0].plot(x, [vals[v]["Q_in"] for v in lv], marks[fam], label=fam)
        ax[1].plot(x, [vals[v]["lift1_amp"] for v in lv], marks[fam], label=fam)
        ax[2].plot(x, [math.degrees(cmath.phase(complex(vals[v]["lift1_re"], vals[v]["lift1_im"]))) for v in lv],
                   marks[fam], label=fam)
    for axi, t in zip(ax, ("Q_in (N/m)", "lift fundamental amplitude (N/m)", "lift fundamental phase (deg)")):
        axi.set_xlabel("(1/L)^2  (L = level; 8x' at 1/64)")
        axi.set_title(t)
        axi.grid(alpha=0.3)
    ax[0].legend()
    fig.tight_layout()
    fig.savefig(a.out, dpi=130)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
