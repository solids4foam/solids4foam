#!/usr/bin/env python3
"""Benchmark quantities, stationarity and coupling statistics of coupled FSI3 runs.

For each case (built by ht_coupled.py):

- periods are delimited by upward crossings of the window mean of the tip
  u_y (hysteresis band 0.25 of the range, as hron_turek_verification.py);
- the analysis window is the last --periods full u_y periods before --t-end
  (default five, about 0.9 s), the same definition on every level;
- per period: (max+min)/2 and (max-min)/2 of u_x, u_y, drag and lift (drag
  and lift = cylinder + flag, per unit depth 0.015 m, as the benchmark), the
  time mean, and the flag-only (forcesPlate) drag and lift; the quantities are
  the averages over the window's periods (the benchmark definition applied
  per period), the frequency is periods/duration;
- stationarity: the same window slid back period by period, with the
  least-squares drift (per second) of every quantity over --drift-from to
  --t-end;
- coupling: iterations per coupled step, the final combined residual
  (fsiResiduals.dat) and, from the logs, the residual1 and residual2 norms at
  acceptance; steps that hit nOuterCorr are counted.

    python3 ht_coupled_analysis.py --case c_1x=~/ht_prod/work/c_1x ... --json out.json --csv out.csv
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import re
from pathlib import Path

import numpy as np

DEPTH = 0.015


def rows(path: Path) -> np.ndarray:
    out = []
    for line in path.read_text().splitlines():
        line = line.replace("(", " ").replace(")", " ")
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        try:
            out.append([float(v) for v in line.split()])
        except ValueError:
            continue
    n = min(len(r) for r in out)
    return np.array([r[:n] for r in out])


def unique_time(a: np.ndarray) -> np.ndarray:
    """Keep the last row of every time (restarts append or rewrite)."""
    _, idx = np.unique(a[::-1, 0], return_index=True)
    return a[::-1][idx]


def history(case: Path) -> dict:
    pp = case / "postProcessing"
    disp = sorted(pp.rglob("solidPointDisplacement*.dat")) or sorted(pp.rglob("*pointDisp*.dat"))
    d = unique_time(np.vstack([rows(p) for p in disp]))
    out = {"t": d[:, 0], "ux": d[:, 1], "uy": d[:, 2]}
    for name in ("forces", "forcesPlate"):
        files = sorted(pp.rglob(f"{name}/*/force.dat"))
        if not files:
            continue
        f = unique_time(np.vstack([rows(p) for p in files]))
        key = "" if name == "forces" else "_plate"
        out["tf" + key] = f[:, 0]
        out["drag" + key] = f[:, 1] / DEPTH
        out["lift" + key] = f[:, 2] / DEPTH
    return out


def crossings(t: np.ndarray, y: np.ndarray) -> np.ndarray:
    mean = y.mean()
    band = 0.25 * (y.max() - y.min())
    out, armed = [], y[0] < mean - band
    for k in range(1, len(y)):
        if y[k] < mean - band:
            armed = True
        elif armed and y[k] >= mean:
            a, b = y[k - 1] - mean, y[k] - mean
            out.append(t[k - 1] + a / (a - b) * (t[k] - t[k - 1]))
            armed = False
    return np.array(out)


def period_stats(t, y, bounds):
    res = []
    for a, b in zip(bounds[:-1], bounds[1:]):
        m = (t >= a) & (t <= b)
        seg, ts = y[m], t[m]
        res.append((0.5 * (seg.max() + seg.min()), 0.5 * (seg.max() - seg.min()),
                    np.trapezoid(seg, ts) / (ts[-1] - ts[0])))
    return np.array(res)


QTY = ("ux", "uy", "drag", "lift", "drag_plate", "lift_plate")


def window_values(h: dict, bounds: np.ndarray) -> dict:
    v = {}
    for q in QTY:
        if q not in h:
            continue
        t = h["t"] if q in ("ux", "uy") else h["tf" + ("_plate" if q.endswith("_plate") else "")]
        ps = period_stats(t, h[q], bounds)
        v[f"{q}_midrange"] = float(ps[:, 0].mean())
        v[f"{q}_midrange_spread"] = float(ps[:, 0].max() - ps[:, 0].min())
        v[f"{q}_amplitude"] = float(ps[:, 1].mean())
        v[f"{q}_timemean"] = float(ps[:, 2].mean())
        v[f"{q}_amplitude_spread"] = float(ps[:, 1].max() - ps[:, 1].min())
        # 4 ms centred moving average (attenuates the 5.5 Hz fundamental by
        # < 0.1%): removes step-to-step coupling noise from the extrema
        dt = float(np.median(np.diff(t)))
        w = max(3, int(round(0.004 / dt)) | 1)
        sm = np.convolve(h[q], np.ones(w) / w, mode="same")
        v[f"{q}_amplitude_smoothed"] = float(period_stats(t, sm, bounds)[:, 1].mean())
    v["frequency"] = float((len(bounds) - 1) / (bounds[-1] - bounds[0]))
    return v


def coupling(case: Path, t0: float, t1: float, n_outer: int) -> dict:
    res = case / "postProcessing/fsiResiduals.dat"
    last = {}
    for line in res.read_text().splitlines()[1:]:
        p = line.split()
        if len(p) >= 3:
            try:
                last[float(p[0])] = (int(p[1]), float(p[2]))
            except ValueError:
                pass
    steps = {t: v for t, v in last.items() if t0 <= t <= t1 + 1e-9}
    it = np.array([v[0] for v in steps.values()])
    r = np.array([v[1] for v in steps.values()])
    out = {"steps": int(len(it)), "iter_mean": float(it.mean()) if len(it) else math.nan,
           "iter_max": int(it.max()) if len(it) else 0,
           "steps_at_max": int((it >= n_outer).sum()),
           "final_residual_max": float(r.max()) if len(r) else math.nan,
           "final_residual_median": float(np.median(r)) if len(r) else math.nan}
    # residual1/residual2 at acceptance, from the logs (last iteration of each step)
    r1, r2, cur, t = {}, {}, None, None
    pat_t = re.compile(r"^Time = ([0-9.eE+-]+)")
    for log in sorted(case.glob("log.solids4Foam*")):
        with log.open(errors="replace") as fh:
            for line in fh:
                if line.startswith("Time = "):
                    m = pat_t.match(line)
                    t = float(m.group(1).rstrip(","))
                elif line.startswith("FSI relative residual1 norm"):
                    r1[t] = float(line.split()[-1])
                elif line.startswith("FSI residual2 norm"):
                    r2[t] = float(line.split()[-1])
    a1 = np.array([v for k, v in r1.items() if k is not None and t0 <= k <= t1 + 1e-9])
    a2 = np.array([v for k, v in r2.items() if k is not None and t0 <= k <= t1 + 1e-9])
    if len(a1):
        out.update({"r1_accept_median": float(np.median(a1)), "r1_accept_max": float(a1.max()),
                    "r2_accept_median": float(np.median(a2)), "r2_accept_max": float(a2.max())})
    return out


def analyse(name: str, case: Path, t_end: float, periods: int, drift_from: float, n_outer: int) -> dict:
    h = history(case)
    m = h["t"] <= t_end + 1e-9
    for k in list(h):
        if k in ("t", "ux", "uy"):
            h[k] = h[k][m] if len(h[k]) == len(m) else h[k]
    for key in ("", "_plate"):
        if "tf" + key in h:
            mf = h["tf" + key] <= t_end + 1e-9
            for q in ("tf", "drag", "lift"):
                h[q + key] = h[q + key][mf]
    late = h["t"] >= drift_from - 0.5
    cr = crossings(h["t"][late], h["uy"][late])
    if len(cr) < periods + 1:
        raise SystemExit(f"{name}: too few periods before {t_end}")
    bounds = cr[-(periods + 1):]
    res = {"name": name, "case": str(case), "t_end": float(h["t"][-1]),
           "window": [float(bounds[0]), float(bounds[-1])], **window_values(h, bounds)}
    # sliding windows (one period apart) for the drift
    slide = []
    for j in range(len(cr) - periods - 1, -1, -1):
        b = cr[j:j + periods + 1]
        if b[0] < drift_from:
            break
        slide.append({"t_end": float(b[-1]), **window_values(h, b)})
    slide = slide[::-1]
    res["sliding"] = slide
    if len(slide) >= 3:
        te = np.array([s["t_end"] for s in slide])
        for k in slide[0]:
            if k == "t_end" or k.endswith("_spread"):
                continue
            y = np.array([s[k] for s in slide])
            res[f"drift_{k}_per_s"] = float(np.polyfit(te, y, 1)[0])
            res[f"range_{k}"] = float(y.max() - y.min())
    setup = {}
    st = case / "coupled_setup.txt"
    if st.is_file():
        for line in st.read_text().splitlines():
            k, _, v = line.partition(" ")
            setup[k] = v
    res["setup"] = setup
    res["coupling_window"] = coupling(case, bounds[0], bounds[-1], n_outer)
    res["coupling_all"] = coupling(case, 2.0 + 1e-9, t_end, n_outer)
    return res


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--case", action="append", required=True, help="name=path")
    ap.add_argument("--t-end", type=float, default=7.0)
    ap.add_argument("--periods", type=int, default=5)
    ap.add_argument("--drift-from", type=float, default=5.0)
    ap.add_argument("--n-outer", type=int, default=50)
    ap.add_argument("--json", type=Path)
    ap.add_argument("--csv", type=Path)
    a = ap.parse_args()
    out = []
    for spec in a.case:
        name, _, path = spec.partition("=")
        out.append(analyse(name, Path(path).expanduser(), a.t_end, a.periods, a.drift_from, a.n_outer))
    cols = ["ux_timemean", "ux_midrange", "ux_amplitude", "uy_amplitude", "frequency", "drag_timemean",
            "drag_amplitude", "lift_timemean", "lift_amplitude"]
    print("name " + " ".join(cols) + " iter_mean r1_med")
    for r in out:
        sc = {c: r[c] * (1e3 if c.startswith("u") else 1) for c in cols}
        print(r["name"], " ".join(f"{sc[c]:.4f}" for c in cols),
              f"{r['coupling_all']['iter_mean']:.2f}", f"{r['coupling_all'].get('r1_accept_median', math.nan):.2e}")
    if a.json:
        a.json.write_text(json.dumps(out, indent=1))
    if a.csv:
        keys = [k for k in out[0] if not isinstance(out[0][k], (dict, list))]
        with a.csv.open("w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(keys + ["iter_mean_all", "steps_at_max_all", "r1_accept_median", "r1_accept_max"])
            for r in out:
                c = r["coupling_all"]
                w.writerow([r.get(k) for k in keys]
                           + [c["iter_mean"], c["steps_at_max"], c.get("r1_accept_median"), c.get("r1_accept_max")])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
