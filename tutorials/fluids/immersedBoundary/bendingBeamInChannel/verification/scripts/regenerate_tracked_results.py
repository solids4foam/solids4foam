#!/usr/bin/env python3
"""Regenerate the tracked files of reference/results/ from the complete saved
runs (standard library only; no solver and no SLURM job is run).

The runs, their fields and their logs are not part of this repository (about
45 GB). The script documents how every tracked file is computed and lets a
holder of the runs reproduce it byte for byte:

    python3 -I scripts/regenerate_tracked_results.py \\
        --runs DIR [DIR ...]          directories that contain the run
                                      directories (searched by name)
        --ledger-advanced FILE        ledger_advanced.csv of the sensitivity
                                      campaign (cost tables)
        --ledger-completion FILE      ledger.csv of the Level 3 completion
                                      campaign (rank-hours of those runs)
        --bodyfitted-archive DIR      the archived body-fitted case
        --out DIR                     where to write (default: a new
                                      directory, never the tracked one)
        [--compare DIR]               report whether each file is identical
                                      to the one in DIR (default: the tracked
                                      reference/results/)
        [--only NAME ...]             cases histories fv cost bodyfitted
                                      diagnostics decomposition courant
                                      inputs

Outputs, with the run names they need:

    cases.csv, history-*.csv      16 runs, analyse_case of the driver
    fv-sensitivity.csv            the 20 control/test pairs of FV_PAIRS
    cost-observed.csv,
    cost-matched-pairs.csv        ledger_advanced.csv and the solver logs
    bodyfitted-R1-*.csv           the archive and the two R1 runs
    diagnostics.csv               filtered (moving mean) nRMS of Section 2.2
    decomposition-scatter.csv     serial against parallel force differences
    bodyfitted-courant.csv        maximum Courant number of the body-fitted
                                  logs (cell-sum definition)
    provenance-inputs.csv         SHA-256 of every raw input file used

The same functions of bending_beam_verification.py compute every metric.
"""

import argparse
import csv
import datetime as dt
import hashlib
import importlib.util
import json
import math
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
TRACKED = HERE.parent / "reference" / "results"

_spec = importlib.util.spec_from_file_location(
    "bbv", HERE / "bending_beam_verification.py")
bbv = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(bbv)

# Runs of cases.csv and history-*.csv: method, level, study, run, ranks and the
# allocated rank-hours (a number as recorded at the time, or the name of the
# run in the completion ledger, or None when the wall time is not recorded)
CASES = [
    ("cutLink", 1, "A", "N1-cutLink-L1-A-wi0125-f33b244a2-v2412", 1, 0.057),
    ("cutLink", 2, "A", "N2L2-cutLink-L2-A-wi0125-f33b244a2-v2412", 1, 0.505),
    ("cutLink", 3, "A", "CA3-cutLink-A-L3-np2", 2, 5.12),
    ("penalty", 1, "A", "PA-penalty-L1-A-wi0125-f33b244a2-v2412", 1, 0.042),
    ("penalty", 2, "A", "PA2-penalty-L2-A-wi0125-f33b244a2-v2412", 1, 0.282),
    ("penalty", 3, "A", "PA3-penalty-A-L3-np2", 2, 2.60),
    ("incremental", 1, "A", "IA-incremental-L1-A-wi0125-f33b244a2-v2412", 1,
     0.140),
    ("incremental", 2, "A", "IA2-incremental-A-L2-np2", 2, 1.256),
    ("incremental", 3, "A", "I3-incremental-L3-A-np8", 8, "ledger"),
    ("ghostCell", 1, "A", "G1c-ghostCell-L1-A-wi0125", 1, "ledger"),
    ("ghostCell", 2, "A", "G2-ghostCell-L2-A-np4", 4, "ledger"),
    ("ghostCell", 3, "A", "G3-ghostCell-L3-A-np8", 8, "ledger"),
    ("cutLink", 1, "B", "N2-cutLink-L1-B0-wi0125-f33b244a2-v2412", 1, None),
    ("cutLink", 2, "B", "CB2-cutLink-B-L2-np2", 2, None),
    ("cutLink", 3, "B", "F-cutLink-B-L3-0.5s", 2, 0.1733),
    ("cutLink", 4, "A", "F-cutLink-A-L4-0.5s", 4, 1.7089),
]

# Finite-volume sensitivity pairs (question, comparison, method, estimator,
# control run, test run); the thresholds are 1.0 point in the nRMS and 0.1 in
# the second harmonic (0.5 and 0.05 for question Q5)
FV_PAIRS = [
    ("Q1", "dt: N48 vs N96 (Euler)", "cutLink", "ST",
     "T-cutLink-Euler-N96", "T-cutLink-Euler-N48"),
    ("Q1", "dt: N48 vs N96 (Euler)", "cutLink", "ME",
     "T-cutLink-Euler-N96", "T-cutLink-Euler-N48"),
    ("Q1", "dt: N48 vs N96 (Euler)", "penalty", "total",
     "T-penalty-Euler-N96", "T-penalty-Euler-N48"),
    ("Q2", "backward vs Euler, N96", "cutLink", "ST",
     "T-cutLink-Euler-N96", "T-cutLink-backward-N96"),
    ("Q2", "backward vs Euler, N96", "cutLink", "ME",
     "T-cutLink-Euler-N96", "T-cutLink-backward-N96"),
    ("Q2", "backward vs Euler, N96", "penalty", "total",
     "T-penalty-Euler-N96", "T-penalty-backward-N96"),
    ("Q2", "backward vs Euler, N48", "cutLink", "ST",
     "T-cutLink-Euler-N48", "T-cutLink-backward-N48"),
    ("Q2", "backward vs Euler, N48", "cutLink", "ME",
     "T-cutLink-Euler-N48", "T-cutLink-backward-N48"),
    ("Q2", "backward vs Euler, N48", "penalty", "total",
     "T-penalty-Euler-N48", "T-penalty-backward-N48"),
    ("Q4", "limitedLinear vs linearUpwind (N96)", "cutLink", "ST",
     "T-cutLink-Euler-N96", "S-cutLink-limitedLinear-N96"),
    ("Q4", "limitedLinear vs linearUpwind (N96)", "cutLink", "ME",
     "T-cutLink-Euler-N96", "S-cutLink-limitedLinear-N96"),
    ("Q4", "limitedLinear vs linearUpwind (N96)", "penalty", "total",
     "T-penalty-Euler-N96", "S-penalty-limitedLinear-N96"),
    ("Q5", "nOuter 2 vs 1 (N96)", "cutLink", "ST",
     "T-cutLink-Euler-N96", "P-cutLink-nOuter2-r2"),
    ("Q5", "nOuter 2 vs 1 (N96)", "cutLink", "ME",
     "T-cutLink-Euler-N96", "P-cutLink-nOuter2-r2"),
    ("Q5", "nOuter 2 vs 1 (N96)", "penalty", "total",
     "T-penalty-Euler-N96", "P-penalty-nOuter2"),
    ("Q5", "nCorr 6 vs 3 (N96)", "penalty", "total",
     "T-penalty-Euler-N96", "P-penalty-nCorr6"),
    ("Q5", "pressure relTol 1e-3 vs 1e-2 (N96)", "penalty", "total",
     "T-penalty-Euler-N96", "P-penalty-prelTol0.001"),
    ("Q4", "upwind vs linearUpwind (diagnostic, N96)", "penalty", "total",
     "T-penalty-Euler-N96", "S-penalty-upwind-N96"),
]
ESTIMATOR = {"ST": "cutLink ST", "ME": "cutLink ME", "total": "total"}
PENALTY_COEFFICIENT = [("S3-penalty-coeff1e2-L2-N96", "1e2"),
                       ("S3-penalty-coeff1e4-L2-N96", "1e4")]

# Matched cutLink/penalty cost pairs
COST_PAIRS = [
    ("Euler N96", "T-cutLink-Euler-N96", "T-penalty-Euler-N96"),
    ("Euler N48", "T-cutLink-Euler-N48", "T-penalty-Euler-N48"),
    ("backward N96", "T-cutLink-backward-N96", "T-penalty-backward-N96"),
    ("backward N48", "T-cutLink-backward-N48", "T-penalty-backward-N48"),
    ("limitedLinear N96", "S-cutLink-limitedLinear-N96",
     "S-penalty-limitedLinear-N96"),
    ("nOuter 2 N96", "P-cutLink-nOuter2-r2", "P-penalty-nOuter2")]

# Local time at which the CPU-affinity defect of the MPI launch (all ranks of a
# job bound to the first cores of the node) was corrected in the campaign
AFFINITY_FIX = dt.datetime(2026, 10, 30, 16, 21)

# Series of the three body-fitted histories: label, directory key
BODYFITTED = {"A": "archive", "B": "R1-bodyfitted-A24720-dt0.0025",
              "C": "R1C-bodyfitted-A24720-dt0.01"}

# Serial against parallel comparisons: (name, reference run, parallel run);
# the window is 0.5 to 1 s, the end of the shortest of these runs
DECOMPOSITION = [
    ("ghostCell L1, 1 against 2 ranks", "PC-G1-np1", "PC-G1-np2"),
    ("ghostCell L1, 1 against 4 ranks", "PC-G1-np1", "PC-G1-np4"),
    ("penalty L2, 2 against 4 ranks", "SM2-penalty-A-L2-np2-1s",
     "SM4-penalty-A-L2-np4-1s"),
    ("incremental L1, 1 against 2 ranks",
     "IA-incremental-L1-A-wi0125-f33b244a2-v2412",
     "SMI2-incremental-A-L1-np2-1s"),
]
DECOMPOSITION_WINDOW = (0.5, 1.0)


class Runs:
    """Finds run directories by name under the given parent directories"""

    def __init__(self, parents):
        self.parents = [Path(p) for p in parents]

    def find(self, name: str) -> Path:
        for parent in self.parents:
            if (parent / name).is_dir():
                return parent / name
        sys.exit(f"run '{name}' not found under "
                 + ", ".join(str(p) for p in self.parents))


def write_csv(path: Path, header: list, rows: list) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)


def force_series(run: Path, params: dict) -> dict:
    columns, rows = bbv.read_force_file(
        run / "postProcessing/immersedBoundary/0/beam.dat")
    return bbv.hydrodynamic_forces(columns, rows, params["coeff_scale"])


def ledger_rank_hours(path: Path, name: str) -> float:
    with path.open() as handle:
        rows = [r for r in csv.DictReader(handle) if r["run"] == name]
    if not rows:
        sys.exit(f"{name} is not in {path}")
    return round(sum(float(r["rank_h"]) for r in rows), 3)


# ---------------------------------------------------------------------------
# cases.csv and history-*.csv
# ---------------------------------------------------------------------------

def make_cases(args, runs, out, params, refs):
    rows, hist = [], {}
    for method, level, study, name, ranks, rank_h in CASES:
        if rank_h == "ledger":
            rank_h = ledger_rank_hours(Path(args.ledger_completion), name)
        run = runs.find(name)
        rows += bbv.analyse_case(run, method, level, study, ranks,
                                 rank_hours=rank_h, refs=refs, params=params)
        if study == "A" and level <= 3:
            hist[(method, level)] = force_series(run, params)
    # the tracked file has no solver columns: the baseline is stated in README
    columns = [c for c in bbv.RESULT_COLUMNS
               if c not in bbv.PROVENANCE_COLUMNS]
    with (out / "cases.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns,
                                extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({k: (f"{v:.6g}" if isinstance(v, float) else v)
                             for k, v in row.items()})
    return rows, hist


def make_histories(hist, out):
    grid = [round(0.02 * i, 2) for i in range(1, 401)]
    for (method, level), series in hist.items():
        keys = [k for k in series if k not in ("t", "forcing", "inertia")]
        header = ["t"] + [f"{c}_{k.replace(' ', '_')}"
                          for k in keys for c in ("Cd", "Cl")]
        body = [[f"{t:.2f}"] + [
            f"{bbv.interpolate(series['t'], series[k][i], t):.5g}"
            for k in keys for i in (0, 1)] for t in grid]
        write_csv(out / f"history-{method}-L{level}.csv", header, body)


# ---------------------------------------------------------------------------
# fv-sensitivity.csv
# ---------------------------------------------------------------------------

def fv_metrics(run, estimator, params, refs):
    t_ref, cd_ref, cl_ref = refs
    series = force_series(run, params)
    t = series["t"]
    cd, cl = series[estimator]
    harm = bbv.harmonics(t, cd)
    return {"nrms": bbv.nrms(t, cd, t_ref, cd_ref), "A2": harm["A2"],
            "noise": bbv.noise(t, cd),
            "nrms_cl": bbv.nrms(t, cl, t_ref, cl_ref)}


def fv_row(question, comparison, method, estimator, control, test, thr):
    dn = test["m"]["nrms"] - control["m"]["nrms"]
    da = test["m"]["A2"] - control["m"]["A2"]
    c, t = control["m"], test["m"]
    return [question, comparison, method, estimator, control["name"],
            test["name"], round(c["nrms"], 3), round(t["nrms"], 3),
            round(dn, 3), round(c["A2"], 4), round(t["A2"], 4), round(da, 4),
            round(t["noise"] / c["noise"], 3),
            round(t["nrms_cl"] - c["nrms_cl"], 3), thr[0], thr[1],
            abs(dn) >= thr[0] or abs(da) >= thr[1]]


def make_fv(runs, out, params, refs):
    header = ["question", "comparison", "method", "estimator", "control_run",
              "test_run", "Cd_nrms_control", "Cd_nrms_test", "d_Cd_nrms",
              "Cd_A2_control", "Cd_A2_test", "d_Cd_A2", "noise_ratio",
              "d_Cl_nrms", "thr_nrms", "thr_A2", "material"]
    body = []
    for q, comparison, method, short, ctrl, test in FV_PAIRS:
        est = ESTIMATOR[short]
        thr = (0.5, 0.05) if q == "Q5" else (1.0, 0.1)
        c = {"name": ctrl, "m": fv_metrics(runs.find(ctrl), est, params, refs)}
        t = {"name": test, "m": fv_metrics(runs.find(test), est, params, refs)}
        body.append(fv_row(q, comparison, method, est, c, t, thr))
    ref_name = "T-penalty-Euler-N96"
    c = {"name": ref_name,
         "m": fv_metrics(runs.find(ref_name), "total", params, refs)}
    for name, label in PENALTY_COEFFICIENT:
        t = {"name": name,
             "m": fv_metrics(runs.find(name), "total", params, refs)}
        body.append(fv_row("Q6", f"penaltyCoeff {label} vs 1e3 (N96)",
                           "penalty", "total", c, t, (1.0, 0.1)))
    write_csv(out / "fv-sensitivity.csv", header, body)


# ---------------------------------------------------------------------------
# cost-observed.csv and cost-matched-pairs.csv
# ---------------------------------------------------------------------------

def solver_log(run: Path) -> str:
    for name in ("log.solids4Foam", "log.pimpleFoam"):
        if (run / name).exists():
            return (run / name).read_text(errors="ignore")
    sys.exit(f"no solver log in {run}")


def cost_audit(runs, ledger: Path):
    """Per accepted run: wall time per step in three windows of the run, the
    CPU/wall ratio, the share of the window before the affinity correction
    and the time-averaged number of the other ranks of the campaign"""
    with ledger.open() as handle:
        all_rows = list(csv.DictReader(handle))
    interval = {}
    for r in all_rows:
        end = dt.datetime.fromisoformat(r["end"])
        interval[r["case"]] = (end - dt.timedelta(seconds=float(r["wall_s"])),
                               end, int(r["ranks"]))

    def concurrency(t0, t1, me):
        total = 0.0
        for key, (s, e, n) in interval.items():
            if key == me:
                continue
            overlap = (min(e, t1) - max(s, t0)).total_seconds()
            total += max(0.0, overlap) * n
        return total / max(1.0, (t1 - t0).total_seconds())

    accepted = [r for r in all_rows if r["status"] == "ACCEPTED"
                and not r["case"].startswith("W0")]
    out = []
    pattern = re.compile(r"ExecutionTime = ([\d.]+) s\s+ClockTime = (\d+) s")
    for r in accepted:
        text = solver_log(runs.find(Path(r["result"]).name))
        found = pattern.findall(text)
        execs = [float(m[0]) for m in found]
        clock = [float(m[1]) for m in found]
        n = len(execs)
        start, _end, ranks = interval[r["case"]]
        offset = float(r["wall_s"]) - clock[-1]

        def window(a, b):
            i, j = int(a * n), int(b * n) - 1
            steps = j - i
            wall = (clock[j] - clock[i]) / steps
            cpu = (execs[j] - execs[i]) / steps
            t0 = start + dt.timedelta(seconds=offset + clock[i])
            t1 = start + dt.timedelta(seconds=offset + clock[j])
            pinned = max(0.0, (min(t1, AFFINITY_FIX) - t0).total_seconds()) \
                / max(1.0, (t1 - t0).total_seconds())
            return {"wall_per_step": wall, "cpu_over_wall": cpu / wall,
                    "t0": t0.strftime("%H:%M"), "t1": t1.strftime("%H:%M"),
                    "after_fix": t0 >= AFFINITY_FIX,
                    "pinned_overlap": pinned,
                    "other_ranks": concurrency(t0, t1, r["case"])}
        out.append({"case": r["case"], "method": r["method"], "ranks": ranks,
                    "steps": n, "allocated_rank_h": float(r["rank_h"]),
                    "first": window(0, 0.2), "mid": window(0.2, 0.6),
                    "last": window(0.6, 1.0)})
    return out


def make_cost(runs, out, ledger):
    audit = cost_audit(runs, Path(ledger))
    by_case = {c["case"]: c for c in audit}
    body = []
    for c in audit:
        if c["method"] not in ("cutLink", "penalty") \
                or c["case"].startswith(("F-", "PC1")):
            continue
        pinned = (0.2 * c["first"]["pinned_overlap"]
                  + 0.4 * c["mid"]["pinned_overlap"]
                  + 0.4 * c["last"]["pinned_overlap"])
        steady = c["ranks"] * c["steps"] * c["last"]["wall_per_step"] / 3600
        last = c["last"]
        body.append([c["case"], c["method"], c["ranks"], c["steps"],
                     pinned > 0.01, round(pinned, 2), c["allocated_rank_h"],
                     round(steady, 3), round(last["wall_per_step"], 4),
                     round(last["cpu_over_wall"], 3),
                     round(last["other_ranks"], 1),
                     last["t0"] + "-" + last["t1"]])
    write_csv(out / "cost-observed.csv",
              ["run", "method", "ranks", "steps", "affinity_defect_affected",
               "share_of_steps_run_pinned", "allocated_rank_h",
               "steady_rank_h_at_last40_rate", "s_per_step_last40",
               "cpu_over_wall_last40", "other_campaign_ranks_active_last40",
               "window_last40"], body)
    pairs = []
    for label, cut, pen in COST_PAIRS:
        a, q = by_case[cut]["last"], by_case[pen]["last"]
        clean = (a["after_fix"] and q["after_fix"]
                 and a["pinned_overlap"] == 0 and q["pinned_overlap"] == 0)
        pairs.append([label, cut, pen, round(a["wall_per_step"], 4),
                      round(q["wall_per_step"], 4),
                      round(a["wall_per_step"] / q["wall_per_step"], 2),
                      round(a["other_ranks"], 1), round(q["other_ranks"], 1),
                      bool(clean)])
    write_csv(out / "cost-matched-pairs.csv",
              ["settings", "cutLink_run", "penalty_run",
               "cutLink_s_per_step", "penalty_s_per_step", "ratio",
               "cutLink_other_ranks", "penalty_other_ranks",
               "both_windows_after_affinity_fix"], pairs)


# ---------------------------------------------------------------------------
# Body-fitted control
# ---------------------------------------------------------------------------

def read_coefficients(case: Path) -> tuple:
    """(t, Cd, Cl) of postProcessing/forceCoeffs/0/coefficient.dat"""
    t, cd, cl = [], [], []
    path = case / "postProcessing/forceCoeffs/0/coefficient.dat"
    for line in path.read_text().splitlines():
        if line.strip() and line[0] != "#":
            v = [float(x) for x in line.split()[:8]]
            t.append(v[0])
            cd.append(v[1])
            cl.append(v[4])
    return t, {"Cd": cd, "Cl": cl}


def trms(t, y):
    return math.sqrt(bbv.trapezoid(t, [v * v for v in y]) / (t[-1] - t[0]))


def grid(step):
    """Sample times 1 + step, 1 + 2 step, ... below 8 s"""
    out, i = [], 1
    while 1 + i * step < 8 - 1e-9:
        out.append(round(1 + i * step, 6))
        i += 1
    return out


def lsq_harmonics(t, y):
    """Direct 5-term least-squares fit (constant, two harmonics of the period
    of the beam) of y(t): fundamental amplitude, second harmonic amplitude and
    phase (degrees)"""
    rows = [[1.0, math.sin(2 * math.pi * x), math.cos(2 * math.pi * x),
             math.sin(4 * math.pi * x), math.cos(4 * math.pi * x)]
            for x in t]
    ata = [[sum(r[i] * r[j] for r in rows) for j in range(5)]
           for i in range(5)]
    aty = [sum(r[i] * v for r, v in zip(rows, y)) for i in range(5)]
    c = bbv.solve_linear(ata, aty)
    return (math.hypot(c[1], c[2]), math.hypot(c[3], c[4]),
            math.degrees(math.atan2(c[4], c[3])))


def make_bodyfitted(runs, out, archive, tracked_ref):
    data = {}
    for key, name in BODYFITTED.items():
        case = Path(archive) if name == "archive" else runs.find(name)
        data[key] = read_coefficients(case)
    t_ref, cd_ref, cl_ref = tracked_ref
    ref = {"Cd": cd_ref, "Cl": cl_ref}

    def on(key, q, g):
        t, s = data[key]
        return [bbv.interpolate(t, s[q], x) for x in g]

    g01 = grid(0.1)
    hist_header = ["t", "Cd_A_dt0p01_archive", "Cl_A_dt0p01_archive",
                   "Cd_B_dt0p0025", "Cl_B_dt0p0025", "Cd_C_dt0p01_control",
                   "Cl_C_dt0p01_control"]
    g = [round(0.1 * i, 1) for i in range(1, 81)]
    cols = [on(k, q, g) for k in "ABC" for q in ("Cd", "Cl")]
    write_csv(out / "bodyfitted-R1-histories.csv", hist_header,
              [[f"{t:.1f}"] + [f"{c[i]:.5g}" for c in cols]
               for i, t in enumerate(g)])

    def native(key, q):
        t, s = data[key]
        idx = [i for i, x in enumerate(t) if 1 + 1e-9 < x < 8 - 1e-9]
        tt, yy = [t[i] for i in idx], [s[q][i] for i in idx]
        return (bbv.trapezoid(tt, yy) / (tt[-1] - tt[0]), trms(tt, yy))

    def difference(k1, k2, q, step):
        gg = grid(step)
        a, b = on(k1, q, gg), on(k2, q, gg)
        return 100 * trms(gg, [x - y for x, y in zip(a, b)]) / trms(gg, b)

    def to_reference(key, q, step):
        gg = grid(step)
        r = [bbv.interpolate(t_ref, ref[q], x) for x in gg]
        d = [x - y for x, y in zip(on(key, q, gg), r)]
        return 100 * trms(gg, d) / trms(gg, r)

    body = []
    for q in ("Cd", "Cl"):
        stats = {k: native(k, q) for k in "ABC"}
        harm = {k: lsq_harmonics(g01, on(k, q, g01)) for k in "ABC"}
        body.append(["mean (1<t<8, native samples)", q]
                    + [round(stats[k][0], 4) for k in "ABC"])
        body.append(["RMS (1<t<8, native samples)", q]
                    + [round(stats[k][1], 4) for k in "ABC"])
        body.append(["fundamental amplitude (direct LSQ, 0.1 s grid)", q]
                    + [round(harm[k][0], 4) for k in "ABC"])
        body.append(["second harmonic amplitude (direct LSQ, 0.1 s grid)", q]
                    + [round(harm[k][1], 4) for k in "ABC"])
        body.append(["second harmonic phase (deg, direct LSQ)", q]
                    + [round(harm[k][2], 2) for k in "ABC"])
        body.append(["nRMS to tracked reference (%), same sample times, no "
                     "interpolation of the force", q]
                    + [round(to_reference(k, q, 0.1), 3) for k in "ABC"])
        body.append(["nRMS to tracked reference (%), 0.01 s grid", q, ""]
                    + [round(to_reference(k, q, 0.01), 3) for k in "BC"])
        for a, b in (("B", "C"), ("A", "C"), ("A", "B")):
            body.append([f"nRMS difference {a}-{b} (%), 0.1 s grid, "
                         "normalised by the second", q,
                         round(difference(a, b, q, 0.1), 3), "", ""])
        body.append(["nRMS difference B-C (%), 0.01 s grid", q,
                     round(difference("B", "C", q, 0.01), 3), "", ""])
    write_csv(out / "bodyfitted-R1-control.csv",
              ["quantity", "coefficient", "A_archive_dt0.01",
               "B_R1_dt0.0025", "C_control_dt0.01"], body)


def make_courant(runs, out, archive):
    """Maximum of the Courant number printed by pimpleFoam (OpenFOAM cell-sum
    definition) in the three body-fitted logs"""
    body = []
    for key, name in BODYFITTED.items():
        case = Path(archive) if name == "archive" else runs.find(name)
        text = (case / "log.pimpleFoam").read_text(errors="ignore")
        co = [float(x) for x in re.findall(
            r"Courant Number mean: \S+ max: (\S+)", text)]
        steps = len(re.findall(r"^Time = ", text, re.M))
        body.append([key, "archive" if name == "archive" else name,
                     steps, max(co)])
    write_csv(out / "bodyfitted-courant.csv",
              ["series", "run", "steps", "max_Courant_cell_sum"], body)


# ---------------------------------------------------------------------------
# Section 2.2: filtered diagnostic
# ---------------------------------------------------------------------------

def moving_mean(t, y, width=0.0104):
    """Centred moving mean over round(width / mean dt) samples (odd, at least
    3), shortened at the ends of the series"""
    dt_mean = (t[-1] - t[0]) / (len(t) - 1)
    k = max(3, int(round(width / dt_mean)) | 1)
    h = k // 2
    out = []
    for i in range(len(y)):
        part = y[max(0, i - h):i + h + 1]
        out.append(sum(part) / len(part))
    return out


def make_diagnostics(rows, hist, out, params, refs):
    t_ref, cd_ref, cl_ref = refs
    order = [("cutLink", "cutLink ST"), ("cutLink", "cutLink ME"),
             ("penalty", "total"), ("incremental", "total"),
             ("ghostCell", "total")]
    steps = {(r["method"], r["level"]): r["steps"] for r in rows
             if r["study"] == "A"}
    body = []
    for method, estimator in order:
        for level in (1, 2, 3):
            series = hist[(method, level)]
            t = series["t"]
            row = [method, estimator, level, steps[(method, level)]]
            for y, ref in zip(series[estimator], (cd_ref, cl_ref)):
                row += [bbv.nrms(t, y, t_ref, ref),
                        bbv.nrms(t, moving_mean(t, y), t_ref, ref)]
            sel = [k for k, x in enumerate(t) if 1 < x < 8 and k + 1 < len(t)]
            dts = [t[k + 1] - t[k] for k in sel]
            row += [sum(dts) / len(dts), min(dts), max(dts)]
            body.append([f"{v:.6g}" if isinstance(v, float) else v
                         for v in row])
    write_csv(out / "diagnostics.csv",
              ["method", "estimator", "level", "steps", "Cd_nrms_raw",
               "Cd_nrms_movmean_0.0104s", "Cl_nrms_raw",
               "Cl_nrms_movmean_0.0104s", "mean_dt_1_8", "min_dt",
               "max_dt"], body)


# ---------------------------------------------------------------------------
# Decomposition scatter
# ---------------------------------------------------------------------------

def make_decomposition(runs, out, params):
    body = []
    lo, hi = DECOMPOSITION_WINDOW
    for label, ref_name, run_name in DECOMPOSITION:
        a = force_series(runs.find(ref_name), params)
        b = force_series(runs.find(run_name), params)
        idx = [i for i, x in enumerate(a["t"]) if lo <= x <= hi + 1e-12]
        same = (len(a["t"]) == len(b["t"])
                and all(abs(x - y) < 1e-12 for x, y in zip(a["t"], b["t"])))
        for name, i in (("Cd", 0), ("Cl", 1)):
            ya = [a["total"][i][k] for k in idx]
            yb = [bbv.interpolate(b["t"], b["total"][i], a["t"][k])
                  for k in idx]
            d = [x - y for x, y in zip(ya, yb)]
            rms = math.sqrt(sum(x * x for x in ya) / len(ya))
            mean = sum(ya) / len(ya)
            body.append([label, ref_name, run_name, name, len(ya), same,
                         round(100 * math.sqrt(sum(x * x for x in d)
                                               / len(d)) / rms, 3),
                         round(100 * max(abs(x) for x in d) / rms, 3),
                         round(100 * abs(sum(d) / len(d)) / rms, 3),
                         round(100 * abs(sum(d) / len(d)) / abs(mean), 3)])
    write_csv(out / "decomposition-scatter.csv",
              ["comparison", "reference_run", "parallel_run", "coefficient",
               "samples_0.5_1_s", "same_times", "rms_diff_pct_of_rms",
               "max_diff_pct_of_rms", "mean_diff_pct_of_rms",
               "mean_diff_pct_of_mean"], body)


# ---------------------------------------------------------------------------
# Input hashes
# ---------------------------------------------------------------------------

def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def make_inputs(args, runs, out):
    names = [c[3] for c in CASES]
    names += sorted({n for p in FV_PAIRS for n in p[4:]})
    names += [n for n, _ in PENALTY_COEFFICIENT]
    names += sorted({n for _, c, p in COST_PAIRS for n in (c, p)})
    names += [n for d in DECOMPOSITION for n in d[1:]]
    names += [BODYFITTED["B"], BODYFITTED["C"]]
    seen, body = set(), []
    for name in names:
        if name in seen:
            continue
        seen.add(name)
        run = runs.find(name)
        for rel in ("postProcessing/immersedBoundary/0/beam.dat",
                    "postProcessing/forceCoeffs/0/coefficient.dat",
                    "log.solids4Foam", "log.pimpleFoam"):
            path = run / rel
            if path.exists():
                body.append([name, rel, path.stat().st_size, sha256(path)])
    for label, path in (("ledger_advanced.csv", args.ledger_advanced),
                        ("ledger_completion.csv", args.ledger_completion)):
        if path:
            body.append(["campaign", label, Path(path).stat().st_size,
                         sha256(Path(path))])
    if args.bodyfitted_archive:
        arch = Path(args.bodyfitted_archive)
        for rel in ("postProcessing/forceCoeffs/0/coefficient.dat",
                    "log.pimpleFoam"):
            if (arch / rel).exists():
                body.append(["body-fitted archive", rel,
                             (arch / rel).stat().st_size,
                             sha256(arch / rel)])
    write_csv(out / "provenance-inputs.csv",
              ["run", "file", "bytes", "sha256"], body)


# ---------------------------------------------------------------------------

STEPS = ("cases", "histories", "fv", "cost", "bodyfitted", "diagnostics",
         "decomposition", "courant", "inputs")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--runs", nargs="+", required=True)
    ap.add_argument("--ledger-advanced")
    ap.add_argument("--ledger-completion")
    ap.add_argument("--bodyfitted-archive")
    ap.add_argument("--out", required=True)
    ap.add_argument("--compare", default=str(TRACKED))
    ap.add_argument("--only", nargs="+", choices=STEPS)
    args = ap.parse_args(argv)
    todo = args.only or STEPS
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    runs = Runs(args.runs)
    params = bbv.parameters()
    refs = bbv.read_series(bbv.REFERENCE_FILE)
    rows = hist = None
    if {"cases", "histories", "diagnostics"} & set(todo):
        rows, hist = make_cases(args, runs, out, params, refs)
    if "histories" in todo:
        make_histories(hist, out)
    if "diagnostics" in todo:
        make_diagnostics(rows, hist, out, params, refs)
    if "fv" in todo:
        make_fv(runs, out, params, refs)
    if "cost" in todo:
        make_cost(runs, out, args.ledger_advanced)
    if "bodyfitted" in todo:
        make_bodyfitted(runs, out, args.bodyfitted_archive, refs)
    if "courant" in todo:
        make_courant(runs, out, args.bodyfitted_archive)
    if "decomposition" in todo:
        make_decomposition(runs, out, params)
    if "inputs" in todo:
        make_inputs(args, runs, out)
    compare = Path(args.compare)
    status = 0
    for path in sorted(out.glob("*.csv")):
        other = compare / path.name
        if not other.exists():
            print(f"NEW        {path.name}")
        elif other.read_bytes() == path.read_bytes():
            print(f"IDENTICAL  {path.name}")
        else:
            print(f"DIFFERENT  {path.name}")
            status = 1
    return status


if __name__ == "__main__":
    sys.exit(main())
