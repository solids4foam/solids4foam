#!/usr/bin/env python3
"""Layer A: NUMERICAL ACCEPTANCE of one radialTorsionalRing run, from the solver log and diagnostics.

This layer answers "did the run behave" and nothing else. A passing layer A is NOT scientific
verification (see verify.py layer B). Statuses: PASS, CONDITIONAL, FAIL, INCOMPLETE, INFO.
A truncated or killed log without a solver error is INCOMPLETE, never PASS.

The solver (fluidSolidInterface.C) accepts an FSI iteration when, with requireAllResidualMeasures
unset, min(residual1, residual2) <= outerCorrTolerance, or when the absolute point mismatch is below
outerCorrAbsoluteTolerance ("accepted by outerCorrAbsoluteTolerance"). The criterion that ended
each time step is recorded from the LAST iteration of the step (the residual evaluation printed before
"Time = t, iteration: 1" belongs to the solid predictor, not to an iteration of the loop).
"""
import sys as _sys
_sys.dont_write_bytecode = True
import json
import math
import os
import re

import numpy as np

import status as ST  # noqa: E402

SEV = ST.ORDER
F = r"([-+0-9.eE]+|nan|-?inf)"
HDR = re.compile(r"^Time = " + F + r"\s*$")
ITR = re.compile(r"^Time = " + F + r", iteration: (\d+)")
SOL = re.compile(r"Solving for (\w+), Initial residual = " + F + r", Final residual = " + F + r", No Iterations (\d+)")
CON = re.compile(r"time step continuity errors : sum local = " + F)
R1 = re.compile(r"FSI relative residual1 norm for interface 0: " + F)
R2 = re.compile(r"FSI residual2 norm for interface 0: " + F)
AB = re.compile(r"FSI absolute residual \(max point\) for interface 0: " + F)
ABACC = re.compile(r"accepted by outerCorrAbsoluteTolerance: " + F + r" < " + F)
SNES = re.compile(r"Nonlinear solve (converged|did not converge) due to (\w+)")
NAN = re.compile(r"(?<![A-Za-z0-9_])-?nan(?![A-Za-z0-9_])")


def worst(statuses):
    return ST.worst(statuses)


def read_dict_value(path, key, cast=float):
    try:
        with open(path, errors="replace") as fh:
            m = re.search(r"^\s*" + re.escape(key) + r"\s+([^;\s]+)\s*;", fh.read(), re.M)
        return cast(m.group(1)) if m else None
    except OSError:
        return None


def case_settings(case_dir, thr):
    """Settings read from the case dictionaries when present, else the documented defaults."""
    A = thr["layerA"]
    fsi = os.path.join(case_dir, "constant", "fsiProperties")
    ctl = os.path.join(case_dir, "system", "controlDict")
    fv = os.path.join(case_dir, "system", "fluid", "fvSolution")
    s = dict(
        tol=read_dict_value(fsi, "outerCorrTolerance") or A["outerCorrTolerance_default"]["value"],
        absTol=read_dict_value(fsi, "outerCorrAbsoluteTolerance"),
        nOuterCorr=read_dict_value(fsi, "nOuterCorr", int) or A["nOuterCorr_default"]["value"],
        deltaT=read_dict_value(ctl, "deltaT"), endTime=read_dict_value(ctl, "endTime"),
        nNonOrth=read_dict_value(fv, "nNonOrthogonalCorrectors", int),
        settings_from_case=os.path.isfile(fsi) and os.path.isfile(ctl))
    s["p_per_cont"] = (s["nNonOrth"] + 1) if s["nNonOrth"] is not None else 3
    try:
        with open(os.path.join(case_dir, "reference.json")) as fh:
            ref = json.load(fh)["args"]
        s["radial_flow"] = abs(float(ref.get("P0", 0.0))) > 0
    except (OSError, ValueError, KeyError):
        s["radial_flow"] = None
    return s


def _residual(rx, line, cur):
    """Value of a printed FSI residual, NaN when the line is malformed, truncated or the number unreadable
    (the step is then recorded in cur["bad"] and the evidence judged INCOMPLETE)."""
    m = rx.search(line)
    try:
        v = float(m.group(1)) if m else math.nan
    except ValueError:
        v = math.nan
    if not math.isfinite(v):
        cur["bad"].append(line.strip()[:120])
    return v


def parse_log(path, thr, cfg):
    """One pass over the log. Returns the raw record used by the checks."""
    steps = {}                    # time key -> dict
    order = []
    cur = None
    k = 0
    nps = 0
    rec = dict(has_end=False, fatal=[], nan=0, snes_bad=[], nsolves_bad=[], p_final_max=0.0, U_final_max=0.0,
               n_cont=0, n_p=0, lines=0, last_line="", p_viol_t=[], U_viol_t=[], p_maxiter=0)
    pm_, um_ = thr["layerA"]["p_final_residual_max"]["value"], thr["layerA"]["U_final_residual_max"]["value"]
    sl_ = 1 + thr["layerA"]["final_residual_slack"]["value"]
    fh = open(path, errors="replace")
    for line in fh:
        rec["lines"] += 1
        if line.startswith("End"):
            rec["has_end"] = True
        # cheap substring tests before regexes
        if "FATAL" in line or "MPI_ABORT" in line or "sigFpe::sigHandler" in line or "core dumped" in line \
                or "did not converge after" in line or "Segmentation" in line:
            if len(rec["fatal"]) < 5:
                rec["fatal"].append(line.strip()[:200])
        if "nan" in line and NAN.search(line):
            rec["nan"] += 1
        if line.startswith("Time = "):
            m = ITR.match(line)
            if m:
                k = int(m.group(2))
                cur["iters_seen"] = max(cur["iters_seen"], k)
                cur["cont_by_iter"].setdefault(k, 0.0)
                cur["ux_by_iter"].setdefault(k, [])
                continue
            m = HDR.match(line)
            if m:
                key = m.group(1)
                cur = dict(t=float(key), key=key, closed=False, trip=[], bad=[], iters_seen=0, cont_by_iter={0: 0.0}, ux_by_iter={0: []})
                order.append(cur)
                steps.setdefault(key, []).append(cur)
                k = 0
                nps = 0
            continue
        if cur is None:
            continue
        if line.startswith("ExecutionTime"):
            cur["closed"] = True
            continue
        if "Solving for" in line:
            m = SOL.search(line)
            if m:
                f, ini, fin, nit = m.group(1), float(m.group(2)), float(m.group(3)), int(m.group(4))
                if f == "p":
                    nps += 1
                    rec["n_p"] += 1
                    rec["p_final_max"] = max(rec["p_final_max"], fin)
                    if fin > pm_ * sl_:
                        rec["p_viol_t"].append(cur["t"])
                        if nit >= 1000:
                            rec["p_maxiter"] += 1
                elif f in ("Ux", "Uy"):
                    rec["U_final_max"] = max(rec["U_final_max"], fin)
                    if fin > um_ * sl_:
                        rec["U_viol_t"].append(cur["t"])
                    if f == "Ux" and k > 0:
                        cur["ux_by_iter"][k].append(ini)
            continue
        if "continuity errors" in line:
            m = CON.search(line)
            if m:
                rec["n_cont"] += 1
                v = abs(float(m.group(1)))
                cur["cont_by_iter"][k] = max(cur["cont_by_iter"].get(k, 0.0), v)
                if nps != cfg["p_per_cont"] and len(rec["nsolves_bad"]) < 20:
                    rec["nsolves_bad"].append((cur["t"], k, nps))
                nps = 0
            continue
        if "Nonlinear solve" in line:
            m = SNES.search(line)
            if m and m.group(1) != "converged" and len(rec["snes_bad"]) < 10:
                rec["snes_bad"].append((cur["t"], line.strip()[:120]))
            continue
        if "FSI relative residual1" in line:
            v = _residual(R1, line, cur)
            cur["trip"].append(dict(r1=v, r2=math.nan, ab=math.nan, acc=False, accval=math.nan))
            continue
        if "FSI residual2" in line:
            v = _residual(R2, line, cur)
            if not cur["trip"] or not math.isnan(cur["trip"][-1]["r2"]):
                cur["bad"].append("residual2 without its residual1 line (residual1 line missing)")
            if cur["trip"]:
                cur["trip"][-1]["r2"] = v
            continue
        if "FSI absolute residual" in line:
            v = _residual(AB, line, cur)
            if not cur["trip"] or not math.isnan(cur["trip"][-1]["ab"]):
                cur["bad"].append("absolute residual without its residual1 line (residual1 line missing)")
            if cur["trip"]:
                cur["trip"][-1]["ab"] = v
            continue
        if "accepted by outerCorrAbsoluteTolerance" in line:
            m = ABACC.search(line)
            if m and cur["trip"]:
                cur["trip"][-1]["acc"] = True
                cur["trip"][-1]["accval"] = float(m.group(1))
    fh.close()
    rec["steps"] = order
    rec["by_key"] = steps
    return rec


def classify_step(st, tol, absTol):
    """Criterion that ended the step (last iteration). Returns (criterion, flags)."""
    if not st["trip"]:
        return "none", ["no FSI residual printed"]
    if not st["closed"]:
        return "none", ["step never completed (no ExecutionTime line)"]
    L = st["trip"][-1]
    a1, a2 = L["r1"] <= tol, L["r2"] <= tol
    flags = []
    if a1 and a2:
        c = "both"
    elif a1:
        c = "r1-only"
    elif a2:
        c = "r2-only"
    elif L["acc"]:
        c = "absolute"
        if absTol is not None and not (L["accval"] < absTol):
            flags.append("absolute residual not below outerCorrAbsoluteTolerance")
    else:
        c = "none"
    if c in ("r1-only", "r2-only"):
        flags.append("single measure; the other residual is above tolerance")
    if c == "absolute":
        flags.append("both relative residuals above tolerance (accepted by absolute tolerance)")
    return c, flags


def contraction(rec, thr):
    """Per-step PIMPLE contraction (audit definition): geometric mean of the Ux initial-residual ratio
    over the passes of the last three FSI iterations of the step."""
    T, R = [], []
    for st in rec["steps"]:
        its = st["ux_by_iter"]
        L = max(its) if its else 0
        rc = math.nan
        if L >= 4:
            seq = its[L - 2] + its[L - 1] + its[L]
            if len(seq) > 1 and seq[0] > 0 and seq[-1] > 0:
                rc = (seq[-1] / seq[0]) ** (1 / (len(seq) - 1))
        T.append(st["t"])
        R.append(rc)
    return np.array(T), np.array(R)


WINDOWS = [(0.5, 1.0), (1.0, 1.5), (1.5, 1.6), (1.6, 1.7), (1.7, 1.75), (1.75, 1.78), (1.78, 1.79),
           (1.79, 1.8), (1.8, 1.81), (1.81, 1.9), (1.9, 2.0)]


def gm(x):
    x = x[np.isfinite(x) & (x > 0)]
    return (float(np.exp(np.log(x).mean())), int(len(x))) if len(x) else (math.nan, 0)


def _run_layer_a(path, thr, cfg_override=None, use_diag=True):
    """path: case directory (log.solids4Foam inside) or a log file. Returns the layer-A report dict."""
    A = thr["layerA"]
    if os.path.isdir(path):
        case_dir, log = path, os.path.join(path, "log.solids4Foam")
    else:
        case_dir, log = os.path.dirname(os.path.abspath(path)), path
        use_diag = use_diag and os.path.basename(path) == "log.solids4Foam"
    cfg = case_settings(case_dir, thr)
    if cfg_override:
        cfg.update(cfg_override)
    out = dict(case=os.path.basename(os.path.abspath(case_dir)), log=os.path.basename(log), settings=cfg, checks={})
    chk = out["checks"]

    def add(name, status, detail, **extra):
        chk[name] = dict(status=status, detail=detail, **extra)

    if not os.path.isfile(log):
        add("log_present", "INCOMPLETE", "no solver log")
        out["status"] = "INCOMPLETE"
        return out
    rec = parse_log(log, thr, cfg)
    S = rec["steps"]
    # a log that stops without End and without an error is truncated: its last step is in progress and is
    # not judged (the run is INCOMPLETE through the termination check)
    trunc = (not rec["has_end"]) and not (rec["fatal"] or rec["nan"] or rec["snes_bad"])
    last_started = S[-1]["t"] if S else 0.0
    if trunc and S:
        S = S[:-1]
        rec["steps"] = S
    tol, absTol, nOuter = cfg["tol"], cfg["absTol"], cfg["nOuterCorr"]
    ts = np.array([s["t"] for s in S])
    # ---- A1 termination
    end = cfg["endTime"] if cfg["endTime"] is not None else float(os.environ.get("RT_END_TIME", 2.0))
    last_t = float(ts[-1]) if len(ts) else 0.0
    if trunc:
        last_t = float(last_started)
    if rec["fatal"] or rec["nan"] or rec["snes_bad"]:
        why = []
        if rec["fatal"]:
            why.append("solver error: " + rec["fatal"][0])
        if rec["nan"]:
            why.append(f"{rec['nan']} NaN lines")
        if rec["snes_bad"]:
            why.append(f"SNES failure at t={rec['snes_bad'][0][0]}: {rec['snes_bad'][0][1]}")
        add("termination", "FAIL", "; ".join(why) + f"; last step started t={last_t:g}",
            fatal=rec["fatal"], nan_lines=rec["nan"], snes_failures=rec["snes_bad"], last_step_t=last_t)
    elif not rec["has_end"]:
        add("termination", "INCOMPLETE", f"no 'End' and no solver error: truncated or killed log; last step started t={last_t:g}",
            last_step_t=last_t)
    else:
        target = end if end is not None else last_t
        if abs(last_t - target) > 1e-6:
            add("termination", "INCOMPLETE", f"End reached but last step t={last_t:g} is short of endTime {target:g}", last_step_t=last_t)
        else:
            add("termination", "PASS", f"End reached at t={last_t:g} (endTime {target:g}), no solver error, no NaN, no SNES failure", last_step_t=last_t)
    # ---- A2 steps monotone and unique
    if len(ts) == 0:
        add("steps_monotone_unique", "INCOMPLETE", "no time step in the log")
    else:
        dup = len(S) - len({s_["key"] for s_ in S})
        dts = np.diff(ts)
        dt0 = cfg["deltaT"] or (ts[0] if len(ts) else 0)
        bad_dt = int(np.sum(np.abs(dts - dt0) > 1e-9 * max(1.0, abs(dt0)) + 1e-12)) if len(dts) else 0
        mono = bool(np.all(dts > 0))
        st = "PASS" if (dup == 0 and mono and bad_dt == 0) else "FAIL"
        add("steps_monotone_unique", st, f"{len(ts)} steps, duplicates {dup}, non-increasing {int(np.sum(dts <= 0))}, steps with dt != {dt0:g}: {bad_dt}")
    # ---- A2b the printed FSI residuals of every judged, completed step must be readable and finite
    bad_ev = [(s_["t"], (s_["bad"] or ["missing residual"])[0]) for s_ in S if s_["closed"] and (s_["bad"] or (
        s_["trip"] and not all(math.isfinite(s_["trip"][-1][k_]) for k_ in ("r1", "r2", "ab"))))]
    if bad_ev:
        add("fsi_residual_evidence", "INCOMPLETE", f"{len(bad_ev)} completed step(s) with a missing, malformed, truncated, NaN or infinite FSI "
            f"residual (first t={bad_ev[0][0]:g}: {bad_ev[0][1]}); the acceptance criterion cannot be evaluated", steps=[t_ for t_, _ in bad_ev[:10]])
    # ---- A3 FSI acceptance criterion per step (final iteration)
    crit = {"both": 0, "r1-only": 0, "r2-only": 0, "absolute": 0, "none": 0}
    per_step = []
    abs_t, single_r2_above, none_t, abs_incons = [], [], [], []
    for s in S:
        c, fl = classify_step(s, tol, absTol)
        crit[c] += 1
        L = s["trip"][-1] if s["trip"] else dict(r1=math.nan, r2=math.nan, ab=math.nan)
        per_step.append(dict(t=s["t"], iterations=len(s["trip"]) - 1 if s["trip"] else 0, criterion=c, r1=L["r1"], r2=L["r2"],
                             abs_residual=L["ab"], flags="; ".join(fl)))
        if c == "absolute":
            abs_t.append(s["t"])
            if any("not below" in f for f in fl):
                abs_incons.append(s["t"])
        if c == "r1-only":
            single_r2_above.append(s["t"])
        if c == "none":
            none_t.append(s["t"])
    out["per_step"] = per_step
    # two different rejection fractions (kept apart): (A) steps failing the both-relative-measures test (absolute steps count as failing);
    # (B) steps the solver's strict mode (requireAllResidualMeasures, which still accepts absolute-tolerance steps) would reject at the accepted iterate
    nA = nB = nSingle = nRescued = 0
    for st_ in S:
        if not st_["trip"]:
            continue
        T_ = st_["trip"][-1]
        both_ = T_["r1"] <= tol and T_["r2"] <= tol
        nA += 0 if both_ else 1
        single_ = (T_["r1"] <= tol) != (T_["r2"] <= tol)
        nSingle += 1 if single_ else 0
        resc_ = absTol is not None and T_["ab"] < absTol
        nRescued += 1 if (single_ and resc_) else 0
        nB += 1 if (max(T_["r1"], T_["r2"]) > tol and not resc_) else 0
    nS = max(len(S), 1)
    out["rejection"] = dict(steps=len(S), A_fail_both_measures=nA, A_pct=100.0 * nA / nS, B_strict_mode_reject=nB, B_pct=100.0 * nB / nS,
                            single_measure=nSingle, single_pct=100.0 * nSingle / nS, single_measure_rescued_by_absolute=nRescued)
    out["criterion_counts"] = crit
    out["abs_accepted_times"] = abs_t
    out["r1_only_steps_with_r2_above_tol"] = dict(count=len(single_r2_above), first=single_r2_above[:3], last=single_r2_above[-3:])
    abs_dis = [t for t in abs_t if t >= A["abs_accept_startup_end_t"]["value"]]
    if S and (none_t or abs_incons):
        add("fsi_step_acceptance", "FAIL", f"{len(none_t)} step(s) ended without an accepting criterion (first t={none_t[:3]}); "
            f"{len(abs_incons)} absolute acceptance(s) inconsistent with outerCorrAbsoluteTolerance", counts=crit, steps_without_criterion=none_t)
    elif not S:
        add("fsi_step_acceptance", "INCOMPLETE", "no steps")
    else:
        flagged = crit["r1-only"] + crit["r2-only"] + crit["absolute"]
        st = "CONDITIONAL" if (flagged or abs_dis) else "PASS"
        add("fsi_step_acceptance", st,
            f"every step accepted by an explicit criterion: both {crit['both']}, r1-only {crit['r1-only']} (r2 above {tol:g}), "
            f"r2-only {crit['r2-only']} (r1 above {tol:g}), absolute {crit['absolute']}"
            + (f"; {len(abs_dis)} absolute acceptance(s) after t={A['abs_accept_startup_end_t']['value']:g}" if abs_dis else "")
            + ("; CONDITIONAL because single-measure and absolute acceptances are flagged explicitly (legacy minimum test)" if flagged else ""),
            counts=crit, abs_accepted_times_min_max=[min(abs_t), max(abs_t)] if abs_t else None,
            abs_accepted_after_startup=abs_dis)
    # ---- A4 iterations
    its = np.array([len(s["trip"]) - 1 if s["trip"] else 0 for s in S])
    if len(its):
        lim = A["max_fsi_iterations"]["value"]
        t0 = A["max_fsi_iterations_from_t"]["value"]
        late = its[ts >= t0]
        over_nouter = int(np.sum(its > nOuter))
        over = [float(t) for t, i in zip(ts, its) if t >= t0 and i > lim]
        st = "FAIL" if (over_nouter or over) else "PASS"
        add("fsi_iterations", st, f"max {int(its.max())} overall (at t={ts[int(its.argmax())]:g}), max {int(late.max()) if len(late) else 0} for t>={t0:g} "
            f"(limit {lim}); steps beyond {lim}: {len(over)}; steps beyond nOuterCorr {nOuter}: {over_nouter}; mean {its.mean():.2f}",
            max_overall=int(its.max()), max_after_startup=int(late.max()) if len(late) else 0, steps_over_limit=over)
    # ---- A5 p solves per continuity
    if rec["n_cont"]:
        st = "PASS" if not rec["nsolves_bad"] else "FAIL"
        add("p_solves_per_continuity", st, f"expected {cfg['p_per_cont']} (nNonOrthogonalCorrectors+1), "
            f"mean {rec['n_p'] / rec['n_cont']:.4f}, deviating groups {len(rec['nsolves_bad'])}" +
            (f" (first {rec['nsolves_bad'][:2]})" if rec["nsolves_bad"] else ""))
    else:
        add("p_solves_per_continuity", "INCOMPLETE", "no continuity line in the log")
    # ---- A6 linear finals
    pm, um = A["p_final_residual_max"]["value"], A["U_final_residual_max"]["value"]
    if rec["n_p"]:
        t_su = A["linear_finals_startup_end_t"]["value"]
        viol = sorted(set(rec["p_viol_t"] + rec["U_viol_t"]))
        late_v = [t for t in viol if t >= t_su]
        st = "PASS" if not viol else ("FAIL" if late_v else "CONDITIONAL")
        add("linear_finals", st, f"max final p {rec['p_final_max']:.3e} (<= {pm:g}), max final U {rec['U_final_max']:.3e} (<= {um:g}); "
            f"{len(rec['p_viol_t'])} p solve(s) and {len(rec['U_viol_t'])} U solve(s) above tolerance"
            + (f" in {len(viol)} step(s) between t={viol[0]:g} and t={viol[-1]:g} ({rec['p_maxiter']} p solves stopped at maxIter 1000); "
               f"start-up (t<{t_su:g}) exceedances are CONDITIONAL, later ones FAIL" if viol else ""),
            p_over_tol=len(rec["p_viol_t"]), U_over_tol=len(rec["U_viol_t"]), steps_with_exceedance=len(viol))
    # ---- A7 continuity (informational)
    if rec["n_cont"]:
        allmax = max((max(s["cont_by_iter"].values()) if s["cont_by_iter"] else 0.0) for s in S) if S else 0.0
        accmax = max((s["cont_by_iter"].get(max(s["cont_by_iter"]), 0.0)) for s in S) if S else 0.0
        radial = cfg["radial_flow"]
        lim = A["continuity_sum_local_max_radial_flow" if radial else "continuity_sum_local_max_torsion"]["value"]
        add("continuity_sum_local", "INFO", f"max over all iterations {allmax:.2e}, max in the accepted (last) iteration {accmax:.2e}; "
            f"documented limit {lim:g} for {'radial-flow (rebased)' if radial else 'torsion-only (or unknown)'} runs, applied to the accepted iterate: "
            f"{'within' if accmax <= lim else 'above'}", all_iterations=allmax, accepted_iterate=accmax)
    # ---- A8 contraction
    if len(S) > 100:
        T, R = contraction(rec, thr)
        c0, c1 = A["contraction_start_t"]["value"], A["contraction_audit_start_t"]["value"]
        gmax = A["contraction_window_gm_max"]["value"]
        rows, bad_hard, bad_soft = [], [], []
        for a, b in WINDOWS:
            if b <= c0 + 1e-9 or a >= (T[-1] if len(T) else 0):
                continue
            m = (T > a + 1e-9) & (T <= b + 1e-9)
            g, n = gm(R[m])
            rows.append(dict(window=f"{a:g}-{b:g}", gm=g, n=n, steps_above_1=int(np.sum(R[m] > 1))))
            if np.isfinite(g) and g > gmax:
                (bad_hard if a >= c1 - 1e-9 else bad_soft).append(f"{a:g}-{b:g}: {g:.3f}")
        runs, cur_run = [], 0
        for t, r in zip(T, R):
            if t >= c1 - 1e-9 and np.isfinite(r) and r > 1:
                cur_run += 1
            else:
                if cur_run:
                    runs.append(cur_run)
                cur_run = 0
        if cur_run:
            runs.append(cur_run)
        longest = max(runs, default=0)
        st = "FAIL" if (bad_hard or longest > A["contraction_run_max"]["value"]) else ("CONDITIONAL" if bad_soft else "PASS")
        add("contraction", st, f"windows from {c0:g} s (start-up excluded); GM limit {gmax:g} enforced from {c1:g} s; "
            f"max window GM {max([r['gm'] for r in rows if np.isfinite(r['gm'])], default=math.nan):.3f}; windows above limit "
            f"from {c1:g} s: {bad_hard or 'none'}, between {c0:g} and {c1:g} s: {bad_soft or 'none'}; longest run of steps >1 from {c1:g} s: {longest}",
            windows=rows)
    # ---- A9 diagnostics-based checks
    if use_diag and os.path.isfile(os.path.join(case_dir, "rtDiagSolid.dat")) and os.path.isfile(os.path.join(case_dir, "rtDiagFluid.dat")):
        try:
            Sd = np.loadtxt(os.path.join(case_dir, "rtDiagSolid.dat"))
            Fd = np.loadtxt(os.path.join(case_dir, "rtDiagFluid.dat"))
            td = Sd[:, 0]
            okrows = (len(Sd) == len(Fd)) and bool(np.all(np.diff(td) > 0)) and (trunc or len(S) == 0 or abs(td[-1] - last_t) < 1e-6)
            add("diagnostics_rows", "PASS" if okrows else ("INCOMPLETE" if trunc else "FAIL"), f"{len(Sd)} solid and {len(Fd)} fluid rows, monotone unique, last t={td[-1]:g}")
            Jmin = float(np.min(Sd[:, 17]))
            add("solid_J", "PASS" if Jmin > A["solid_Jmin_min"]["value"] else "FAIL",
                f"min J {Jmin:.5f} (> {A['solid_Jmin_min']['value']}), max J {float(np.max(Sd[:, 18])):.4f}, min principal stretch {float(np.min(Sd[:, 16])):.4f}")
        except (OSError, ValueError, IndexError) as e:
            add("diagnostics_rows", "FAIL", f"cannot read diagnostics: {e}")
    else:
        add("diagnostics_rows", "NOT-VERIFIED", "rtDiag*.dat not available for this log")
    add("mesh_Vmin_GCL_ALE_modes", "NOT-VERIFIED", "field-level checks (minimum cell volume, swept-volume GCL, ALE, A_mode growth) need field scripts that are not in the tutorial")
    out["status"] = ST.decide([v["status"] for v in chk.values()])
    return out


def run_layer_a(path, thr, cfg_override=None, use_diag=True):
    """Layer A of one run. Never raises: an unreadable or garbled log or diagnostics file gives INCOMPLETE."""
    try:
        return _run_layer_a(path, thr, cfg_override, use_diag)
    except Exception as e:
        name = os.path.basename(os.path.abspath(path))
        return dict(case=name, log="", settings={}, status="INCOMPLETE",
                    checks={"input": dict(status="INCOMPLETE", detail=f"cannot be analysed: {type(e).__name__}: {e}")})
