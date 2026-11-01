#!/usr/bin/env python3
"""Controls of the oscillatingWallStokesLayer verification: time-step
convergence, the momentum-exchange wall offset and the cell-average comparison.
Everything is computed from the tracked data in reference/results/ (no CFD, no
numpy), with the exact solutions of the verification driver.

    stokes_controls.py temporal    time-step convergence table
    stokes_controls.py offset      momentum exchange minus forcing against the
                                   acceleration of the slab and the offset
                                   between the wall and the nearest mesh face
    stokes_controls.py cellavg     velocity L2 error against the exact value at
                                   the cell centre and against the exact cell
                                   average
    stokes_controls.py all         the three, written to reference/results/controls

Thresholds, fixed before the data were interpreted (report 12 of the
campaign): a series of time steps is converged to a bound only if three
successive halvings give monotone differences whose ratio lies in [1.6, 2.5]
(first order) or [3, 5] (second order); two agreeing values are not
convergence.
"""

import csv
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import oscillating_wall_stokes_layer_verification as v  # noqa: E402

RESULTS = Path(__file__).resolve().parent.parent / "reference" / "results"
CONTROLS = RESULTS / "controls"
RATIO_FIRST = (1.6, 2.5)
RATIO_SECOND = (3.0, 5.0)
STEP_FACTOR = 2.0
EXACT_RMS = None            # set in main: RMS of the exact periodic wall shear
QUANTITIES = ("u_L2_all_ivp", "shear_amp_err_ivp", "shear_rms_err_ivp")


# ------------------------------------------------------------ time-step series
def read_rows(path: Path) -> list:
    with path.open() as f:
        return list(csv.DictReader(f))


def merged_cases(results: Path = RESULTS) -> list:
    """The tutorial-baseline rows and the control rows, each key once"""
    rows = read_rows(results / "cases.csv")
    extra = results / "controls" / "cases-timestep.csv"
    if extra.exists():
        seen = {row_key(r) for r in rows}
        rows += [r for r in read_rows(extra) if row_key(r) not in seen]
    return rows


def row_key(row: dict) -> tuple:
    return (row["method"], row["estimator"], int(row["level"]),
            round(float(row["dt"]), 12))


def halving_series(rows: list, method: str, estimator: str, level: int,
                   quantity: str, dt0: float = 0.001) -> list:
    """[(dt, value)] of the unbroken sequence dt0, dt0/2, dt0/4, ... that
    exists in the rows (the sequence stops at the first missing step)"""
    by_dt = {round(float(r["dt"]), 12): float(r[quantity]) for r in rows
             if r["method"] == method and r["estimator"] == estimator
             and int(r["level"]) == level and r[quantity] != ""}
    out, dt = [], dt0
    while round(dt, 12) in by_dt:
        out.append((dt, by_dt[round(dt, 12)]))
        dt /= STEP_FACTOR
    return out


def assess(values: list) -> dict:
    """Convergence assessment of the values at successive halvings of the
    step. Returns the differences, the last two ratios, the verdict, the
    Richardson limit and a bound on the remaining change."""
    n = len(values)
    d = [values[i + 1] - values[i] for i in range(n - 1)]
    res = {"n": n, "last_diff": d[-1] if d else float("nan"),
           "ratio": float("nan"), "verdict": "insufficient", "limit": "",
           "bound": ""}
    if n < 4:
        return res
    d3 = d[-3:]
    scale = max(abs(values[-1]), 1e-30)
    if n >= 5 and abs(d3[1]) <= 0.01 * scale and abs(d3[2]) <= 0.01 * scale:
        # the last two changes are below 1 % of the value (needs five steps)
        res["verdict"] = "saturated"
        res["limit"] = values[-1]
        res["bound"] = abs(d3[2])
        return res
    r1 = d3[0] / d3[1] if d3[1] != 0 else float("inf")
    r2 = d3[1] / d3[2] if d3[2] != 0 else float("inf")
    res["ratio"] = math.sqrt(abs(r1 * r2)) if math.isfinite(r1 * r2) else float("inf")
    same_sign = d3[0] * d3[1] > 0 and d3[1] * d3[2] > 0
    if not same_sign:
        res["verdict"] = "non-monotone"
        return res
    in_first = all(RATIO_FIRST[0] <= r <= RATIO_FIRST[1] for r in (r1, r2))
    in_second = all(RATIO_SECOND[0] <= r <= RATIO_SECOND[1] for r in (r1, r2))
    if in_first or in_second:
        r = res["ratio"]
        res["verdict"] = "first-order" if in_first else "second-order"
        res["limit"] = values[-1] + d3[2] / (r - 1.0)
        res["bound"] = abs(d3[2] / (r - 1.0))
    else:
        res["verdict"] = "ratio outside the thresholds"
    return res


def temporal_table(rows: list) -> list:
    out = []
    keys = sorted({(r["method"], r["estimator"], int(r["level"])) for r in rows})
    for method, est, level in keys:
        for q in QUANTITIES:
            s = halving_series(rows, method, est, level, q)
            if len(s) < 3:
                continue
            a = assess([x[1] for x in s])
            out.append({"method": method, "estimator": est, "level": level,
                        "quantity": q, "points": a["n"],
                        "dt_first": s[0][0], "dt_last": s[-1][0],
                        "value_first": s[0][1], "value_last": s[-1][1],
                        "last_difference": a["last_diff"],
                        "ratio": a["ratio"], "verdict": a["verdict"],
                        "limit": a["limit"], "bound": a["bound"]})
    return out


# ------------------------------------------------- momentum-exchange wall offset
def wall_offset(p: dict, level: int) -> float:
    """Signed distance, in cell heights, from the wall to the nearest mesh
    face: round(y_w/h) - y_w/h. The cutLink momentum exchange acts at the
    nearest face, the forcing estimator at the surface"""
    ny, _, _, _ = v.mesh_shape(p, level)
    q = p["yw"] / (p["height"] / ny)
    return round(q) - q


def slab_acceleration(p: dict, t: float) -> float:
    a, w = p["amplitude"], p["omega"]
    return -a * w * w * math.sin(w * t + p["phase"])


def acceleration_coefficient(p: dict, t: list, diff: list) -> tuple:
    """Least-squares c of diff(t) = c a(t) and the RMS residual relative to
    the RMS of the difference (near 1 when the difference itself is noise)"""
    acc = [slab_acceleration(p, x) for x in t]
    c = sum(d * a for d, a in zip(diff, acc)) / sum(a * a for a in acc)
    res = v.rms([d - c * a for d, a in zip(diff, acc)])
    return c, res / max(v.rms(diff), 1e-30)


def read_shear(path: Path) -> dict:
    with path.open() as f:
        r = list(csv.DictReader(f))
    return {k: [float(x[k]) for x in r] for k in r[0]}


def exact_shear_rms(p: dict) -> float:
    amp, _ = v.exact_amp_phase(p)
    return amp / math.sqrt(2.0)


def offset_table(p: dict, results: Path = RESULTS) -> list:
    """Fitted rho*dH/h of (ME - FO) and of (ST - FO) against the predicted
    offset, at the tutorial step (tracked shear histories of the fourth
    period), and for the control with the wall on a mesh face"""
    rows = []
    for level in v.LEVELS:
        path = results / f"shear-L{level}.csv"
        if not path.exists():
            continue
        s = read_shear(path)
        h = p["height"] / v.mesh_shape(p, level)[0]
        entry = {"case": "tutorial", "yw": p["yw"], "level": level,
                 "predicted_offset_over_h": wall_offset(p, level)}
        for lab, col in (("ME", "cutLink-ME"), ("ST", "cutLink-ST")):
            c, r = acceleration_coefficient(
                p, s["time"], [a - b for a, b in zip(s[col], s["cutLink-FO"])])
            entry[f"fitted_{lab}_over_h"] = c / h
            entry[f"residual_{lab}"] = r
            entry[f"diff_rms_{lab}"] = v.rms([a - b for a, b in zip(s[col], s["cutLink-FO"])]) / EXACT_RMS
        rows.append(entry)
    for level in v.LEVELS:
        path = results / "controls" / f"shear-yw0125-L{level}.csv"
        if not path.exists():
            continue
        s = read_shear(path)
        q = dict(p)
        q["yw"] = 0.125
        h = p["height"] / v.mesh_shape(q, level)[0]
        c, r = acceleration_coefficient(
            q, s["time"], [a - b for a, b in zip(s["cutLink-ME"], s["cutLink-FO"])])
        m_me = v.shear_metrics(q, s["time"], s["cutLink-ME"])
        m_fo = v.shear_metrics(q, s["time"], s["cutLink-FO"])
        rows.append({"case": "wall on a mesh face", "yw": 0.125, "level": level,
                     "predicted_offset_over_h": wall_offset(q, level),
                     "fitted_ME_over_h": c / h, "residual_ME": r,
                     "diff_rms_ME": v.rms([a - b for a, b in zip(s["cutLink-ME"], s["cutLink-FO"])]) / EXACT_RMS,
                     "amp_err_ME": m_me["shear_amp_err_ivp"],
                     "amp_err_FO": m_fo["shear_amp_err_ivp"]})
    return rows


# --------------------------------------------------------- cell-average metric
GAUSS = ((0.0, 0.5688888888888889),
         (0.5384693101056831, 0.4786286704993665),
         (-0.5384693101056831, 0.4786286704993665),
         (0.9061798459386640, 0.2369268850561891),
         (-0.9061798459386640, 0.2369268850561891))


def physical_velocity(p: dict, s: float, t: float, h_fluid: float, modes: list) -> float:
    """exact initial-value solution for s > 0, rigid slab velocity inside"""
    if s > 0:
        return v.exact_u_ivp(p, s, t, h_fluid, modes)
    return v.slab_velocity(p, t)


def cell_average(p: dict, s: float, h: float, t: float, modes: list) -> float:
    """average of the physical field over [s - h/2, s + h/2] with the kink at
    the wall integrated separately"""
    a, b = s - h / 2, s + h / 2
    parts = [(a, 0.0), (0.0, b)] if a < 0 < b else [(a, b)]
    hf = p["height"] - p["yw"]
    total = 0.0
    for x, y in parts:
        total += 0.5 * (y - x) * sum(
            w * physical_velocity(p, 0.5 * (x + y) + 0.5 * (y - x) * g, t, hf, modes)
            for g, w in GAUSS)
    return total / (b - a)


def cellavg_table(p: dict, results: Path = RESULTS) -> list:
    """L2 error against the exact value at the cell centre and against the
    exact cell average, over the fluid rows (centre above the wall) and over
    the rows wholly in the fluid, from the tracked profiles of the fourth
    period at the tutorial step"""
    modes = v.cached_modes(p)
    hf = p["height"] - p["yw"]
    out = []
    for level in v.LEVELS:
        path = results / f"profiles-L{level}.csv"
        if not path.exists():
            continue
        h = p["height"] / v.mesh_shape(p, level)[0]
        data = read_rows(path)
        names = [c for c in data[0] if c not in ("phase", "s", "exact_periodic", "exact_ivp")]
        phases = sorted({float(r["phase"]) for r in data})
        for name in names:
            acc = {"centre": [], "average": [], "centre_whole": [], "average_whole": []}
            for ph in phases:
                rr = [r for r in data if float(r["phase"]) == ph]
                t = 3.0 + ph
                e = {k: [] for k in acc}
                for r in rr:
                    s = float(r["s"])
                    u = float(r[name])
                    ec = abs(u - v.exact_u_ivp(p, s, t, hf, modes))
                    ea = abs(u - cell_average(p, s, h, t, modes))
                    e["centre"].append(ec)
                    e["average"].append(ea)
                    if s >= h / 2:
                        e["centre_whole"].append(ec)
                        e["average_whole"].append(ea)
                for k in acc:
                    acc[k].append(v.rms(e[k]) / p["U0"])
            out.append({"level": level, "method": name,
                        **{k: v.rms(acc[k]) for k in acc}})
    return out


# --------------------------------------------------------------------- output
def write_csv(path: Path, rows: list) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    keys = []
    for r in rows:
        for k in r:
            if k not in keys:
                keys.append(k)
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys, restval="", lineterminator="\n")
        w.writeheader()
        for r in rows:
            w.writerow({k: (f"{x:.6g}" if isinstance(x, float) else x)
                        for k, x in r.items()})


def main(argv=None) -> int:
    cmd = (argv or sys.argv[1:] or ["all"])[0]
    p = v.parameters()
    global EXACT_RMS
    EXACT_RMS = exact_shear_rms(p)
    if cmd in ("temporal", "all"):
        write_csv(CONTROLS / "temporal-convergence.csv", temporal_table(merged_cases()))
    if cmd in ("offset", "all"):
        write_csv(CONTROLS / "offset-audit.csv", offset_table(p))
    if cmd in ("cellavg", "all"):
        write_csv(CONTROLS / "cell-average.csv", cellavg_table(p))
    if cmd not in ("temporal", "offset", "cellavg", "all"):
        print(__doc__)
        return 2
    print(f"written to {CONTROLS}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
