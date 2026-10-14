#!/usr/bin/env python3
"""Compare coupled FSI3 runs step by step after a common time t0.

Used for (a) the solid-preconditioner equivalence check (bjacobi/LU against hypre BoomerAMG,
both restarted from the same write) and (b) the restart check (continuous against restarted).
For every case against the reference, over the common steps after t0:

  - tip displacement u_x, u_y (solidPointDisplacement): max |difference|, relative to the u_y
    amplitude of the reference over the window, and the least-squares drift of the difference;
  - total drag and lift (forces, cylinder + flag, per unit depth): max |difference|, relative to
    the lift amplitude, drift;
  - FSI iterations per step (mean, max, steps at nOuterCorr) and the final combined residual;
  - the solid SNES norm at acceptance (last norm of the last solid solve of each step) and
    the first SNES norm of the predictor after t0.

    python3 ht_branch_compare.py --t0 4.5 --ref rst_r4.5 --case eq_hypre_r4.5 [--case ...] [--json out.json]
"""
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import numpy as np

DEPTH = 0.015


def table(files: list[Path], ncol: int) -> np.ndarray:
    rows = []
    for f in files:
        for line in f.read_text().splitlines():
            line = line.replace("(", " ").replace(")", " ")
            if not line.strip() or line.lstrip().startswith("#"):
                continue
            try:
                rows.append([float(v) for v in line.split()[:ncol]])
            except ValueError:
                pass
    a = np.array(rows)
    _, idx = np.unique(a[::-1, 0], return_index=True)
    return a[::-1][idx]


def series(case: Path) -> dict:
    pp = case / "postProcessing"
    d = table(sorted(pp.rglob("solidPointDisplacement*.dat")), 3)
    f = table(sorted(pp.glob("fluid/forces/*/force.dat")), 3)
    return {"disp": d, "force": f}


def coupling(case: Path) -> dict:
    pp = case / "postProcessing"
    files = sorted(pp.glob("fsiResiduals.dat.upto_*"), key=lambda p: float(p.name.split("upto_")[1]))
    last = {}
    for res in files + [pp / "fsiResiduals.dat"]:
        if res.is_file():
            for line in res.read_text().splitlines()[1:]:
                p = line.split()
                try:
                    last[round(float(p[0]), 7)] = (int(p[1]), float(p[2]))
                except (ValueError, IndexError):
                    pass
    # solid SNES norms: first of each step (predictor) and last of each step (acceptance)
    first, acc, t, cur = {}, {}, None, None
    pat = re.compile(r"^Time = ([0-9.eE+-]+)\s*$")
    for log in sorted(case.glob("log.solids4Foam*")):
        with log.open(errors="replace") as fh:
            for line in fh:
                m = pat.match(line)
                if m:
                    t = round(float(m.group(1)), 7)
                elif "SNES Function norm" in line and t is not None:
                    v = float(line.split()[-1])
                    first.setdefault(t, v)
                    acc[t] = v
    return {"fsi": last, "snes_first": first, "snes_acc": acc}


def compare(ref: dict, cas: dict, t0: float, nout: int) -> dict:
    out = {}
    for key, cols, names, scale_col in (("disp", (1, 2), ("ux", "uy"), 2), ("force", (1, 2), ("drag", "lift"), 2)):
        r, c = ref[key], cas[key]
        r, c = r[r[:, 0] > t0 + 1e-9], c[c[:, 0] > t0 + 1e-9]
        rm = {round(t, 7): i for i, t in enumerate(r[:, 0])}
        pairs = [(rm[round(t, 7)], j) for j, t in enumerate(c[:, 0]) if round(t, 7) in rm]
        if not pairs:
            continue
        ir, ic = map(np.array, zip(*pairs))
        tt = r[ir, 0]
        amp = 0.5 * (r[ir, scale_col].max() - r[ir, scale_col].min())
        f = 1.0 / DEPTH if key == "force" else 1.0
        out[f"{key}_steps"] = int(len(tt))
        out[f"{key}_t"] = [float(tt[0]), float(tt[-1])]
        for col, name in zip(cols, names):
            dd = (c[ic, col] - r[ir, col]) * f
            out[f"max_abs_d_{name}"] = float(np.abs(dd).max())
            out[f"rel_d_{name}"] = float(np.abs(dd).max() / (amp * f))
            out[f"drift_d_{name}_per_s"] = float(np.polyfit(tt, dd, 1)[0]) if len(tt) > 2 else 0.0
            out[f"rms_d_{name}"] = float(np.sqrt(np.mean(dd ** 2)))
            out[f"first_nonzero_d_{name}_step"] = int(np.argmax(np.abs(dd) > 0)) if np.any(dd) else -1
        out[f"amp_{names[1]}_ref"] = float(amp * f)
    for label, cp in (("ref", ref["cpl"]), ("case", cas["cpl"])):
        ts = sorted(t for t in cp["fsi"] if t > t0 + 1e-9)
        common = [t for t in ts if t in cas["cpl"]["fsi"] and t in ref["cpl"]["fsi"]]
        it = np.array([cp["fsi"][t][0] for t in common])
        res = np.array([cp["fsi"][t][1] for t in common])
        acc = np.array([cp["snes_acc"][t] for t in common if t in cp["snes_acc"]])
        out[f"{label}_iter_mean"] = float(it.mean()) if len(it) else float("nan")
        out[f"{label}_iter_max"] = int(it.max()) if len(it) else 0
        out[f"{label}_steps_at_max"] = int((it >= nout).sum())
        out[f"{label}_fsi_res_median"] = float(np.median(res)) if len(res) else float("nan")
        out[f"{label}_snes_acc_median"] = float(np.median(acc)) if len(acc) else float("nan")
        out[f"{label}_snes_acc_max"] = float(acc.max()) if len(acc) else float("nan")
        t1 = common[0] if common else None
        out[f"{label}_first_step_snes_predictor"] = cp["snes_first"].get(t1, float("nan")) if t1 else float("nan")
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--t0", type=float, required=True)
    ap.add_argument("--ref", type=Path, required=True)
    ap.add_argument("--case", type=Path, action="append", required=True)
    ap.add_argument("--n-outer", type=int, default=30)
    ap.add_argument("--json", type=Path)
    a = ap.parse_args()
    ref = series(a.ref)
    ref["cpl"] = coupling(a.ref)
    res = {}
    for c in a.case:
        cas = series(c)
        cas["cpl"] = coupling(c)
        res[c.name] = r = compare(ref, cas, a.t0, a.n_outer)
        print(f"== {c.name} against {a.ref.name}, after t0 = {a.t0}")
        for k, v in r.items():
            print(f"  {k:32s} {v}")
    if a.json:
        a.json.write_text(json.dumps({"ref": str(a.ref), "t0": a.t0, "cases": res}, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
