#!/usr/bin/env python3
"""radialTorsionalRing verification driver: two separate layers, explicit statuses.

LAYER A  NUMERICAL ACCEPTANCE (per run, from the solver log and diagnostics; acceptance.py).
LAYER B  SCIENTIFIC reference-error expectations (per family of three meshes; errors against the
         semi-analytical reference, observed and self-convergence orders, mixed terms).

A passing layer A is NOT scientific verification. Layer B is computed only from runs whose layer A is
PASS or CONDITIONAL; otherwise it is INVALID-INPUT. Statuses: PASS, CONDITIONAL, FAIL, INCOMPLETE,
NOT-VERIFIED (layer A), INVALID-INPUT (layer B). Exit code 1 if any layer-A run is FAIL or INCOMPLETE,
2 on usage errors, else 0.

Usage:
  verify.py WORK_DIR [--out OUT_DIR] [--thresholds FILE] [--log PATH ...]
    WORK_DIR   contains <config>_<mesh> case directories (config torsion|radial|combined, mesh M1|M2|M3),
               each with log.solids4Foam, rtDiagFluid.dat, rtDiagSolid.dat, reference.json
               (directories or symbolic links; read only, nothing is written there)
    --log PATH extra log file or case directory checked by layer A only (repeatable)
    --out      default WORK_DIR/verification_report  (verification_report.json, layerA_runs.csv,
               layerA_steps_<run>.csv, layerB_errors.csv, summary.txt)
"""
import sys
sys.dont_write_bytecode = True
import csv
import json
import re
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import acceptance as ACC  # noqa: E402
import status as ST  # noqa: E402
import family as FAM  # noqa: E402
import analyse_case as AC  # noqa: E402
import convergence as CV  # noqa: E402

CONFIGS = ("torsion", "radial", "combined")
MESHES = ("M1", "M2", "M3")

NOT_VERIFIED = [
    "radial-only M1-M3 family: INCOMPLETE (at nNonOrthogonalCorrectors 2 only M1 completed, M2 stalled at t = 0.11 s, M3 not run); no radial convergence claim",
    "three-way radial term d_u3 = u_C - u_R - u_T (needs the radial family)",
    "asymptotic convergence (three meshes cannot prove an asymptotic order)",
    "temporal order (a single time step 0.0025 s; temporal contamination 5-25% of the M3 error)",
    "origin of the combined reference-error floor (only partially explained)",
    "effect of the absolute coupling tolerance on M3 (quantified on M2 only; 1e-12 on M1/M2 and 1e-10 on M3)",
    "field-level checks: swept-volume GCL, ALE consecutive-save, A_mode growth, minimum cell volume (need field scripts that are not in the tutorial)",
    "independent validation of the reference (it is only numerically resolved)",
]


def band_status(v, lo, hi, margin, noise_lo=0.0, noise_hi=0.0):
    """PASS only if v lies inside [lo + noise_lo, hi - noise_hi] (the deterministic band shrunk by the
    documented numerical-noise allowance). Inside the band but closer to an edge than the allowance:
    CONDITIONAL (marginal). Outside the band by at most `margin` (fraction): CONDITIONAL; beyond: FAIL."""
    if lo + noise_lo <= v <= hi - noise_hi:
        return "PASS"
    if lo <= v <= hi:
        return "CONDITIONAL"
    d = (lo - v) / abs(lo) if v < lo else (v - hi) / abs(hi)
    return "CONDITIONAL" if d <= margin else "FAIL"


def worst(sts):
    return ST.worst(sts)


def cap(status, family_status):
    """A family with documented differences (CONDITIONAL) cannot give a layer-B PASS."""
    return ST.cap_pass(status, "CONDITIONAL") if family_status == "CONDITIONAL" else status


def mixed_allowance(noise, max_ref, signal):
    """Allowance of a mixed term in error/signal units: noise (fraction of max|ref|) x max|reference quantity of the combined
    run| / max|mixed-term reference signal| (the reference quantity and the signal are in the same units)."""
    return noise * max_ref / signal


def noise_limited(iv, bound):
    """True if the interval half-width exceeds the documented bound (or the interval is unresolved)."""
    return iv is None or (iv[3] - iv[2]) / 2.0 > bound


def interval_text(iv, bound):
    if iv is None:
        return "unresolved (noise >= error)"
    hw = (iv[3] - iv[2]) / 2.0
    return f"{iv[2]:.2f}-{iv[3]:.2f}" + (" noise-limited" if hw > bound else (" sharp" if hw < 0.02 else ""))


def interval_status(iv, lo, hi):
    """Downgrade only: an unresolved interval (noise >= error) or an interval leaving the band is CONDITIONAL; otherwise PASS (neutral)."""
    if iv is None:
        return "CONDITIONAL"
    return "CONDITIONAL" if (iv[2] < lo or iv[3] > hi) else "PASS"


def mixed_claim(mixed, thr):
    """Status and detail of the mixed-term claim on M3 from CV.mixed_terms output; the allowance is computed with mixed_allowance."""
    B, N = thr["layerB_deterministic"], thr["numerical_noise"]
    margin = B["band_conditional_margin"]["value"]
    sts, det = [], []
    for wn, key in (("0-2", "mixed_term_max_error_over_signal_M3_0_2"), ("0-1.8", "mixed_term_max_error_over_signal_M3_0_1p8")):
        mm = (mixed or {}).get("M3", {}).get(wn, {})
        for k, lim in B[key]["value"].items():
            if k in mm and mm[k]["error_over_signal"] is not None:
                v = mm[k]["error_over_signal"]
                nzq = N["noise_E_by_qoi_combined"]["value"][N["mixed_term_allowance"]["qoi_of_term"][k]]
                a_ = mixed_allowance(nzq, mm[k]["max_abs_ref_combined"], mm[k]["signal_max_abs"])
                st_ = band_status(v, 0.0, lim, margin, 0.0, a_)
                sts.append(st_)
                det.append(f"{wn} {k} {100 * v:.2f}% (limit {100 * lim:.1f}%, allowance {100 * a_:.3g}% = {nzq:g} x {mm[k]['max_abs_ref_combined']:.4g} / {mm[k]['signal_max_abs']:.4g}, {st_})")
    return (worst(sts) if sts else "NOT-VERIFIED"), "error/signal of the mixed terms on M3: " + "; ".join(det)


def layer_b(work, runs_a, thr):
    """runs_a: {(config, mesh): layerA report}. Returns the layer-B report.
    Order of the gates: layer A of every member (PASS or CONDITIONAL) -> family consistency
    (INCOMPLETE / INVALID-INPUT stop the assessment; CONDITIONAL caps every claim at CONDITIONAL)
    -> errors, orders, mixed terms."""
    B = thr["layerB_deterministic"]
    N = thr["numerical_noise"]
    margin = B["band_conditional_margin"]["value"]
    prim = B["primary_qoi"]["value"]
    rep = dict(families={}, claims={}, errors=[])
    CV.OVERRIDE.clear()
    cases = {}
    for key, a in runs_a.items():
        if key[0] in ("torsion", "combined") and a["status"] in ("PASS", "CONDITIONAL"):
            try:
                t, num, ref, _ = AC.load(os.path.join(work, f"{key[0]}_{key[1]}"))
                cases[key] = dict(t=t, num=num, ref=ref, dt=float(t[1] - t[0]))
            except Exception:          # truncated, garbled or missing series: the run is not usable
                pass
    for k, v in cases.items():
        CV.OVERRIDE[k] = v
    for (cfg, mesh), S in sorted(cases.items()):
        for wn, (a, b) in CV.WINDOWS.items():
            m = CV.select(S, a, b)
            for q in CV.QOI[cfg]:
                if q in S["num"] and q in S["ref"] and np.max(np.abs(S["ref"][q][m])) > 0:
                    sc = float(np.max(np.abs(S["ref"][q][m])))
                    e = S["num"][q][m] - S["ref"][q][m]
                    rep["errors"].append(dict(config=cfg, mesh=mesh, window=wn, qoi=q, E_max=float(np.max(np.abs(e))) / sc,
                                              E_L2=float(np.sqrt(np.mean(e ** 2))) / sc))
    for cfg in ("torsion", "combined"):
        need = [(cfg, m) for m in MESHES]
        missing = [k for k in need if k not in runs_a]
        invalid = [k for k in need if k in runs_a and k not in cases]
        if missing or invalid:
            why = []
            if missing:
                why.append("missing runs " + ",".join(f"{k[1]}" for k in missing))
            if invalid:
                why.append("layer A not PASS/CONDITIONAL for " + ",".join(f"{k[1]}({runs_a[k]['status']})" for k in invalid))
            rep["families"][cfg] = dict(status="INVALID-INPUT" if invalid else "NOT-VERIFIED", detail="; ".join(why))
            continue
        fam = FAM.check_family(work, cfg, thr)
        if fam["status"] in ("INCOMPLETE", "INVALID-INPUT"):
            bad = [f"{f['check']}: {f['detail']}" for f in fam["findings"] if f["status"] in ("INCOMPLETE", "INVALID-INPUT")]
            rep["families"][cfg] = dict(status=fam["status"], detail="family consistency: " + "; ".join(bad), consistency=fam)
            continue
        CV.NOISE["torsion"] = dict(N["noise_E_by_qoi"]["value"]["torsion"])
        CV.NOISE["combined"] = dict(N["noise_E_by_qoi"]["value"]["combined"])
        try:
            conv = CV.convergence(work, cfg)
            mix = {m: CV.mixed_terms(work, m) for m in MESHES} if cfg == "combined" else None
        except Exception as e:     # empty or inconsistent series
            rep["families"][cfg] = dict(status="INCOMPLETE", detail=f"errors could not be computed: {type(e).__name__}: {e}", consistency=fam)
            continue
        rep["families"][cfg] = dict(status="COMPUTED", family_status=fam["status"], consistency=fam, convergence=conv, mixed=mix)

    def claim(name, status, detail, family=None, **kw):
        if family is not None:
            status = cap(status, family)
        rep["claims"][name] = dict(status=status, detail=detail, **kw)

    tf = rep["families"].get("torsion")
    if tf and tf["status"] == "COMPUTED":
        lo, hi = B["torsion_order_band"]["value"]
        hl, hh = B["torsion_order_headline_band"]["value"]
        sts, outside, allp, l2p, unc, unresolved = [], [], [], [], [], []
        for wn, res in tf["convergence"].items():
            for q in prim:
                r = res.get(q)
                if not r:
                    continue
                o = r["observed_order"]
                if not isinstance(o, list):
                    sts.append("FAIL")
                    outside.append(f"{wn} {q}: non-monotone errors")
                    continue
                iv = r["order_interval_max_norm"]
                for k, p_ in zip(("p12", "p23"), o):
                    sts.append(band_status(p_, lo, hi, margin))
                    allp.append(p_)
                    if not (hl <= p_ <= hh):
                        outside.append(f"{wn} {q} {p_:.2f}")
                    sts.append(interval_status(iv.get(k), lo, hi))       # unresolved or leaving the band: downgrade only
                    if iv.get(k):
                        unc.append((wn, q, k, iv[k][2], iv[k][3]))
                    else:
                        unresolved.append(f"{wn} {q} {k}")
                o2 = r.get("observed_order_L2")
                if isinstance(o2, list):
                    for p_ in o2:
                        sts.append(band_status(p_, lo, hi, margin))
                        l2p.append(p_)
        nl_b = N["order_uncertainty"]["noise_limited_halfwidth"]["value"]
        t_unc = [u for u in unc if u[0] == "0-2" and u[2] == "p23"]
        claim("torsion_primary_orders", worst(sts),
              f"observed orders of the 9 primary rotational QoIs, windows 0-2 and 0-1.8 s: {min(allp):.2f}-{max(allp):.2f} in max-norm and "
              f"{min(l2p):.2f}-{max(l2p):.2f} in L2-in-time (band {lo}-{hi}); max-norm values outside 1.8-2.0 (listed individually, not rounded): "
              f"{'; '.join(outside) or 'none'}. At the documented per-quantity noise (2.4e-4 of max|ref| for the torque-type quantities, much less for alpha_B and H_f) the second-pair orders carry first-order uncertainty intervals, "
              f"0-2 s: " + "; ".join(f"{q} p23 {interval_text([0, 0, a, b], nl_b)}" for (w, q, k, a, b) in t_unc) +
              (f"; unresolved (noise >= error, a downgrade): {', '.join(unresolved)}" if unresolved else "") + ". Intervals that leave the band, and unresolved ones, downgrade, never upgrade. "
              "Near-second-order spatial convergence observed, an asymptotic order is not established; no temporal order claim.",
              family=tf["family_status"], order_intervals=unc)
        for wn, res in tf["convergence"].items():
            for q in B["secondary_qoi_excluded_from_orders"]["value"]:
                if q in res:
                    o = res[q]["observed_order"]
                    rep.setdefault("secondary_info", []).append(f"torsion {wn} {q}: orders {o if not isinstance(o, list) else [round(x, 2) for x in o]} (excluded from order claims)")
    else:
        claim("torsion_primary_orders", tf["status"] if tf else "NOT-VERIFIED", tf["detail"] if tf else "no torsion runs")
    cf = rep["families"].get("combined")
    if cf and cf["status"] == "COMPUTED":
        fam = cf["family_status"]
        fb = B["combined_floor_band_M3"]
        lo, hi = fb["value"]
        nz = N["combined_max_norm_noise_of_E"]["value"]
        conv0 = cf["convergence"]["0-2"]
        # per-quantity floor status: the floor statement holds only for the five listed quantities
        per_q, sts, det = {}, [], []
        for q in prim + ["u_B", "sig_rr", "p_b"]:
            r = conv0.get(q)
            if not r:
                continue
            e1, e2, e3 = r["E"]
            so = r["self_order"]
            if q in fb["qoi"]:
                st_ = band_status(e3, lo, hi, margin, nz, nz)
                st_ = "DISCLOSED" if st_ == "PASS" else st_          # the floor is never a verification PASS
                sts.append(st_)
                det.append(f"{q} E_M3 {e3:.2e} [{st_}]")
            else:
                st_ = "ABOVE-FLOOR" if e3 > hi else "CONDITIONAL"
                det.append(f"{q} E_M3 {e3:.2e} [{st_}: outside the floor statement]")
            per_q[q] = dict(status=st_, E_M1=e1, E_M2=e2, E_M3=e3, self_order=so, in_floor_statement=q in fb["qoi"])
        rep["floor_by_quantity"] = per_q
        stf = worst(sts)
        measured = [conv0[q]["E"][2] for q in fb["qoi"] if q in conv0]
        claim("combined_reference_error_floor", stf,
              f"floor statement (about 2.5-5.2e-3 of max|ref|) applies ONLY to {', '.join(fb['qoi'])}; measured M3 errors {min(measured):.2e}-{max(measured):.2e} "
              f"(2.5e-3 is the lower edge of the band, not a measured value; the lowest measured is {min(measured):.2e}). By quantity (noise allowance {nz:g}): " + "; ".join(det) +
              ". A disclosed floor against the reference, not clean second-order convergence to it; origin only partly explained. Never a verification PASS.", family=fam)
        so = B["combined_self_order_band"]["value"]
        nl_b = N["order_uncertainty"]["noise_limited_halfwidth"]["value"]
        sts, robust, limited, others = [], [], [], []
        for q in prim:
            r = conv0.get(q)
            if not r or r["self_order"] is None:
                continue
            iv = r.get("self_order_interval")
            if q in B["combined_second_order_qoi"]["value"]:
                st_ = band_status(r["self_order"], so[0], so[1], margin)
                st_ = worst([st_, interval_status(iv, so[0] - 0.0, so[1] + 0.0)]) if st_ == "PASS" else st_
                sts.append(st_)
                txt = f"{q} {r['self_order']:.2f} (interval {interval_text(iv, nl_b)})"
                (limited if noise_limited(iv, nl_b) else robust).append(txt)
            else:
                others.append(f"{q} {r['self_order']:.2f} (E_M3 {r['E'][2]:.2e})")
        claim("combined_self_convergence", worst(sts),
              f"self-difference orders 0-2 s (band {so[0]}-{so[1]}). Robustly resolved near second order: " + (", ".join(robust) or "none") +
              ". NOISE-LIMITED (nominal order near 2 but the noise-propagated interval is wide; these do not independently demonstrate a precise second-order rate): " +
              (", ".join(limited) or "none") + ". No second-order self-convergence (informational): " + ", ".join(others) +
              " (H_s 1.63, w_wall 1.66 and P_wall 0.70 are well below 2; KE_s nominal self-order 2.17 is noise-limited, interval 1.42-2.92, and its reference error is above the floor). The floor origin is only partially explained "
              "(about one third of the rotation excess follows from the radial wall-traction discrepancy; not a proven unique mechanism)", family=fam)
        st_m, det_m = mixed_claim(cf["mixed"], thr)
        claim("mixed_terms_M3", st_m, det_m, family=fam)
        extra = []
        for wn in ("0-2", "0-1.8"):
            x = (cf["mixed"] or {}).get("M3", {}).get(wn, {}).get("d_Hf")
            if x and x["error_over_signal"] is not None:
                a_ = N["noise_E_by_qoi_combined"]["value"]["H_f"] * x["max_abs_ref_combined"] / x["signal_max_abs"]
                extra.append(f"{wn} d_Hf {100 * x['error_over_signal']:.2f}% (informational, no band; allowance {100 * a_:.3f}%)")
        rep.setdefault("secondary_info", []).extend(extra)
    else:
        for n in ("combined_reference_error_floor", "combined_self_convergence", "mixed_terms_M3"):
            claim(n, cf["status"] if cf else "NOT-VERIFIED", cf["detail"] if cf else "no combined runs")
    claim("radial_only_family", "NOT-VERIFIED", NOT_VERIFIED[0])
    return rep


def write_outputs(out, runs, layerb, extra):
    os.makedirs(out, exist_ok=True)
    rj = dict(statement="Layer A is numerical acceptance only. It is not scientific verification. Layer B is the scientific layer.",
              layerA={n: {k: v for k, v in r.items() if k != "per_step"} for n, r in runs.items()},
              layerA_extra_logs={n: {k: v for k, v in r.items() if k != "per_step"} for n, r in extra.items()},
              layerB=layerb, not_verified=NOT_VERIFIED)
    with open(os.path.join(out, "verification_report.json"), "w") as fh:
        json.dump(rj, fh, indent=1, default=lambda x: x.item() if hasattr(x, "item") else str(x))
    with open(os.path.join(out, "layerA_runs.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["run", "status", "check", "check_status", "detail"])
        for n, r in list(runs.items()) + list(extra.items()):
            for c, v in r["checks"].items():
                w.writerow([n, r["status"], c, v["status"], v["detail"]])
    for n, r in list(runs.items()) + list(extra.items()):
        if r.get("per_step"):
            with open(os.path.join(out, f"layerA_steps_{n}.csv"), "w", newline="") as fh:
                w = csv.writer(fh)
                w.writerow(["t", "iterations", "criterion", "r1", "r2", "abs_residual", "flags"])
                for s in r["per_step"]:
                    w.writerow([f"{s['t']:.6g}", s["iterations"], s["criterion"], s["r1"], s["r2"], s["abs_residual"], s["flags"]])
    with open(os.path.join(out, "layerB_errors.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["config", "mesh", "window", "qoi", "E_max", "E_L2"])
        for e in layerb["errors"]:
            w.writerow([e["config"], e["mesh"], e["window"], e["qoi"], f"{e['E_max']:.4e}", f"{e['E_L2']:.4e}"])
    return rj


def summary(runs, extra, layerb):
    L = []
    L.append("=" * 78)
    L.append("radialTorsionalRing verification")
    L.append("=" * 78)
    L.append("LAYER A - NUMERICAL ACCEPTANCE (per run; this is NOT scientific verification)")
    for n, r in list(runs.items()) + list(extra.items()):
        L.append(f"  [{r['status']:<11}] {n}")
        for c, v in r["checks"].items():
            if v["status"] != "PASS" or c in ("termination", "fsi_step_acceptance"):
                L.append(f"      {v['status']:<12} {c}: {v['detail']}")
        rj = r.get("rejection")
        if rj and rj["steps"]:
            L.append(f"      rejection at 1e-6: (A) fail the both-relative-measures test {rj['A_pct']:.2f}% ({rj['A_fail_both_measures']}); "
                     f"(B) rejected by the solver's strict mode (absolute-tolerance steps accepted) {rj['B_pct']:.2f}% ({rj['B_strict_mode_reject']}); "
                     f"single-measure share {rj['single_pct']:.1f}%")
        cc = r.get("criterion_counts")
        if cc:
            ab = r["abs_accepted_times"]
            L.append(f"      criterion per step: {cc}; absolute-accepted times {ab[0]:g}..{ab[-1]:g} s ({len(ab)} steps)" if ab else f"      criterion per step: {cc}")
    L.append("")
    L.append("LAYER B - SCIENTIFIC reference-error expectations (only from runs with layer A PASS/CONDITIONAL)")
    for cfg, f in layerb["families"].items():
        cons = f.get("consistency")
        if cons:
            L.append(f"  family consistency {cfg}: {cons['status']}")
            for x in cons["findings"]:
                if x["status"] != "PASS":
                    L.append(f"      {x['status']:<13} {x['check']}: {x['detail']}")
        if f["status"] != "COMPUTED":
            L.append(f"  [{f['status']}] {cfg} family: {f['detail']}")
    for n, c in layerb["claims"].items():
        L.append(f"  [{c['status']:<13}] {n}: {c['detail']}")
    fq = layerb.get("floor_by_quantity")
    if fq:
        L.append("  combined floor, per quantity (0-2 s; errors of M1, M2, M3 and the self-order):")
        for q, v in fq.items():
            so = "n/a" if v["self_order"] is None else f"{v['self_order']:.2f}"
            L.append(f"      {v['status']:<12} {q:<8} {v['E_M1']:.2e} {v['E_M2']:.2e} {v['E_M3']:.2e}  self {so}" +
                     ("" if v["in_floor_statement"] else "   (outside the floor statement)"))
    for x in layerb.get("secondary_info", [])[:12]:
        L.append("      info: " + x)
    L.append("")
    L.append("NOT-VERIFIED (printed in every run):")
    for x in NOT_VERIFIED:
        L.append("  - " + x)
    L.append("")
    return "\n".join(L)


def run_label(path, taken=()):
    """Sanitised label of an extra run: parent directory and file stem for a log file, directory name for a case directory.
    Never contains a path separator, so output file names stay under the output directory; made unique if it is already taken."""
    ap = os.path.abspath(path)
    fname = os.path.basename(ap)
    stem = os.path.splitext(fname)[0] if os.path.splitext(fname)[1] in (".log", ".txt", ".out") else fname
    base = os.path.basename(ap) if os.path.isdir(ap) else (os.path.basename(os.path.dirname(ap)) + "_" + stem)
    lab = re.sub(r"[^A-Za-z0-9_.-]+", "_", base).strip("._") or "run"
    out, k = lab, 2
    while out in taken:
        out, k = f"{lab}_{k}", k + 1
    return out


def main(argv):
    if not argv or argv[0].startswith("-"):
        print(__doc__)
        return 2
    work = argv[0]
    out = os.path.join(work, "verification_report")
    thr_file = os.path.join(HERE, "..", "thresholds.json")
    logs = []
    i = 1
    while i < len(argv):
        if argv[i] in ("--out", "--thresholds", "--log"):
            if i + 1 >= len(argv) or argv[i + 1].startswith("--"):
                print(f"error: option {argv[i]} requires a value", file=sys.stderr); return 2
            if argv[i] == "--out":
                out = argv[i + 1]
            elif argv[i] == "--thresholds":
                thr_file = argv[i + 1]
            else:
                logs.append(argv[i + 1])
            i += 2
        else:
            print("unknown option", argv[i]); return 2
    if not out.strip():
        print("error: option --out requires a non-empty directory name", file=sys.stderr); return 2
    try:
        with open(thr_file) as fh_:
            thr = json.load(fh_)
    except (OSError, ValueError) as e_:
        print(f"error: cannot read the thresholds file {thr_file!r}: {e_}", file=sys.stderr); return 2
    runs, keyed = {}, {}
    for cfg in CONFIGS:
        for mesh in MESHES:
            d = os.path.join(work, f"{cfg}_{mesh}")
            if os.path.isdir(d):
                r = ACC.run_layer_a(d, thr)
                runs[f"{cfg}_{mesh}"] = r
                keyed[(cfg, mesh)] = r
    extra = {}
    for p in logs:
        r = ACC.run_layer_a(p, thr)
        extra[run_label(p, extra)] = r
    layerb = layer_b(work, keyed, thr)
    write_outputs(out, runs, layerb, extra)
    txt = summary(runs, extra, layerb)
    with open(os.path.join(out, "summary.txt"), "w") as fh:
        fh.write(txt)
    print(txt)
    bad = [n for n, r in list(runs.items()) + list(extra.items()) if r["status"] in ST.BLOCKING]
    if not runs and not extra:
        print("no runs found")
        return 2
    if bad:
        print("EXIT 1: numerical acceptance FAIL/INCOMPLETE for: " + ", ".join(bad))
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
