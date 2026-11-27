#!/usr/bin/env python3
"""Generators of the tracked data of reference/results/ from completed
cases of the driver (standard library only; the numerics are those of the
driver, analyse_case and the exact solutions).

    export_results.py cases   --work-dir D [D ...] [--out DIR]
        cases.csv, controls/cases-timestep.csv, shear-L<k>.csv and
        profiles-L<k>.csv from the completed cases of the work directories
        (the runs 'matrix', 'dt' and the controls 'dt-L1', 'dt-fine',
        'dt-6e-5', 'dt-3e-5')
    export_results.py wall    --work-dir D [D ...] [--out DIR]
        controls/shear-yw0125-L<k>.csv from the run 'wall-face'
    export_results.py scaling --work-dir D [D ...] [--out FILE]
        scaling.csv from the runs of Level 4 penalty with --end-time 1 on
        1, 2, 4, 8 ranks

The default output is verification/work/export, never the tracked files:
copy or pass --out reference/results to regenerate them. The column rank_h
is ranks x the wall time of Allrun.verify (blockMesh included), where the
tracked column is the allocated rank-hours of the scheduler job; the rows of
the other columns do not depend on it.
"""

import argparse
import json
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import oscillating_wall_stokes_layer_verification as v  # noqa: E402

ORDER = ("cutLink-ME", "cutLink-ST", "cutLink-FO", "penalty", "incremental",
         "ghostCell")
PROFILE_ORDER = ("cutLink", "penalty", "incremental", "ghostCell")
PHASES = ((0.0, 3.0), (0.25, 3.25), (0.5, 3.5), (0.75, 3.75))


def completed(work_dirs: list) -> list:
    """[(name, entry, case directory)] of the completed cases of the work
    directories (a name must not occur twice)"""
    out, seen = [], set()
    for d in work_dirs:
        manifest = Path(d) / "manifest.json"
        if not manifest.exists():
            v.fail(f"{manifest} not found")
        for name, entry in sorted(json.loads(manifest.read_text())
                                  ["cases"].items()):
            if entry.get("state") != "completed":
                continue
            if name in seen:
                v.fail(f"case {name} occurs in two work directories")
            seen.add(name)
            out.append((name, entry, Path(d) / "cases" / name))
    return out


def rank_hours(case: Path, ranks: int) -> float:
    wall = v.wall_seconds(case)
    return ranks * wall / 3600.0


def is_baseline(entry: dict) -> bool:
    return not entry.get("end_time") and not entry.get("wall_y")


def in_cases_csv(row: dict) -> bool:
    """reference/results/cases.csv: the tutorial step on every level and the
    steps of the run 'dt' on Levels 2-4; the rest is controls/"""
    dt = v.row_dt(row)
    return dt == v.TUTORIAL_DT or (
        int(row["level"]) >= 2 and any(abs(dt - s) < 1e-12
                                       for s in v.DT_STEPS))


def export_cases(cases: list, out: Path) -> dict:
    rows, runs = [], []
    for name, e, case in cases:
        if not is_baseline(e):
            continue
        new, series = v.analyse_case(case, e["method"], e["estimator"],
                                     e["level"], e["ranks"])
        for r in new:
            r["rank_h"] = f"{rank_hours(case, e['ranks']):.4f}"
            rows.append(r)
        dt = e["dt"] or v.TUTORIAL_DT
        if abs(dt - v.TUTORIAL_DT) < 1e-12:
            runs.append((e, case, series))
    main_rows = [r for r in rows if in_cases_csv(r)]
    ctrl_rows = [r for r in rows if not in_cases_csv(r)]
    if main_rows:
        v.write_results(out / "cases.csv", main_rows)
    if ctrl_rows:
        v.write_results(out / "controls" / "cases-timestep.csv", ctrl_rows)
    p = v.parameters()
    H = p["height"] - p["yw"]
    modes = v.cached_modes(p)
    for level in v.LEVELS:
        sel = [x for x in runs if x[0]["level"] == level]
        if not sel:
            continue
        cols = {}
        for e, case, series in sel:
            for lab, (t, tau) in series.items():
                vid = e["method"] + ("-" + lab if e["method"] == "cutLink"
                                     else "")
                cols[vid] = (t, tau)
        names = [o for o in ORDER if o in cols]
        t0 = cols[names[0]][0]
        with (out / f"shear-L{level}.csv").open("w") as f:
            f.write("time,exact_periodic,exact_ivp," + ",".join(names) + "\n")
            for i, tt in enumerate(t0):
                if not 3.0 - 1e-9 <= tt <= 4.0 + 1e-9:
                    continue
                f.write(f"{tt:.6g},{v.exact_tau(p, tt):.6e},"
                        f"{v.exact_tau_ivp(p, tt, H, modes):.6e},"
                        + ",".join(f"{cols[n][1][i]:.5e}" for n in names)
                        + "\n")
        h, s = v.cell_heights(p, level)
        flows = {}
        for e, case, _ in sel:
            if e["method"] == "cutLink" and e["estimator"] != "ME":
                continue
            flows[e["method"]] = case
        pnames = [o for o in PROFILE_ORDER if o in flows]
        with (out / f"profiles-L{level}.csv").open("w") as f:
            f.write("phase,s,exact_periodic,exact_ivp," + ",".join(pnames)
                    + "\n")
            for phase, tk in PHASES:
                prof = {n: v.row_profile(p, level, v.read_vector_field(
                    flows[n] / f"{tk:g}" / "U")) for n in pnames}
                for j in range(len(s)):
                    if s[j] <= 0:
                        continue
                    f.write(f"{phase},{s[j]:.6g},{v.exact_u(p, s[j], tk):.6e},"
                            f"{v.exact_u_ivp(p, s[j], tk, H, modes):.6e},"
                            + ",".join(f"{prof[n][j]:.5e}" for n in pnames)
                            + "\n")
    return {"rows": len(rows), "cases.csv": len(main_rows),
            "controls": len(ctrl_rows)}


def export_wall(cases: list, out: Path) -> list:
    """controls/shear-yw<y>-L<k>.csv: time and the shear histories of the
    cutLink momentum exchange and forcing runs, 3 <= t <= 4"""
    groups = {}
    for name, e, case in cases:
        if not e.get("wall_y") or e.get("end_time") or e["method"] != "cutLink":
            continue
        _, series = v.analyse_case(case, e["method"], e["estimator"],
                                   e["level"], e["ranks"], e["wall_y"])
        for lab in ("ME", "FO"):
            if lab in series:
                groups.setdefault((e["wall_y"], e["level"]), {})[
                    "cutLink-" + lab] = series[lab]
    written = []
    for (wall, level), cols in sorted(groups.items()):
        names = [n for n in ("cutLink-ME", "cutLink-FO") if n in cols]
        path = out / "controls" / f"shear-yw{round(wall * 1000):04d}-L{level}.csv"
        path.parent.mkdir(parents=True, exist_ok=True)
        t0 = cols[names[0]][0]
        with path.open("w") as f:
            f.write("time," + ",".join(names) + "\n")
            for i, tt in enumerate(t0):
                if 3.0 - 1e-9 <= tt <= 4.0 + 1e-9:
                    f.write(f"{tt:.6g}," + ",".join(
                        f"{cols[n][1][i]:.5e}" for n in names) + "\n")
        written.append(path)
    return written


def export_scaling(cases: list, out: Path) -> str:
    rows = []
    for name, e, case in cases:
        if not (e["method"] == "penalty" and e["level"] == 4
                and e.get("end_time") == 1.0 and not e.get("wall_y")):
            continue
        log = (case / "log.solids4Foam").read_text(errors="replace")
        clock = re.findall(r"ClockTime = (\d+) s", log)
        steps = len(re.findall(r"^Time = ", log, re.MULTILINE))
        rows.append((e["ranks"], steps, float(clock[-1])))
    rows.sort()
    if not rows:
        v.fail("no completed Level 4 penalty case with --end-time 1")
    t1 = next((w for n, _, w in rows if n == 1), rows[0][2] * rows[0][0])
    text = "ranks,steps,solve_wall_s,speedup,efficiency\n" + "".join(
        f"{n},{s},{w:.0f},{t1 / w:.2f},{t1 / w / n:.2f}\n"
        for n, s, w in rows)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(text)
    return text


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("what", choices=("cases", "wall", "scaling"))
    ap.add_argument("--work-dir", nargs="+", required=True)
    ap.add_argument("--out", default="")
    a = ap.parse_args(argv)
    cases = completed(a.work_dir)
    default = v.WORK_DIR / "export"
    if a.what == "scaling":
        out = Path(a.out) if a.out else default / "scaling.csv"
        print(export_scaling(cases, out), end="")
        print(f"written {out}")
    else:
        out = Path(a.out) if a.out else default
        out.mkdir(parents=True, exist_ok=True)
        if a.what == "cases":
            print(export_cases(cases, out))
        else:
            for p in export_wall(cases, out):
                print(f"written {p}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
