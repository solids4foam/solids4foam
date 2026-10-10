#!/usr/bin/env python3
"""Summarise replay-variant runs: Q_in, Q_quad, W, force amplitude and cost.

Runs energy_analysis.analyse on each case (same last-5-periods window as the
energy-balance study) and adds the wall-clock core-hours from run.out.

    python3 variants_summary.py <case> [<case> ...] --csv out.csv --json out.json
"""
from __future__ import annotations
import argparse, csv, json, sys
from datetime import datetime
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import energy_analysis as ea

FIELDS = ["case", "variants", "restart", "t_end", "window0", "window1", "Q_in", "Q_quad", "Q_in_p", "W", "W_p",
          "W_v", "plate_lift_amp1", "plate_drag_mean", "gross_abs_work", "core_hours"]


def core_hours(case: Path) -> float:
    f = case / "run.out"
    if not f.is_file():
        return float("nan")
    lines = f.read_text().split()
    try:
        s = next(l for l in f.read_text().splitlines() if l.startswith("start"))
        e = next(l for l in f.read_text().splitlines() if l.startswith("exit"))
        n = int(s.split("n=")[1])
        t0 = datetime.fromisoformat(s.split()[1])
        t1 = datetime.fromisoformat(e.split()[-1])
        return n * (t1 - t0).total_seconds() / 3600.0
    except Exception:
        return float("nan")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cases", nargs="+", type=Path)
    ap.add_argument("--csv", type=Path)
    ap.add_argument("--json", type=Path)
    ap.add_argument("--periods", type=int, default=5)
    a = ap.parse_args()
    rows, full = [], []
    for c in a.cases:
        r = ea.analyse(c, a.periods, None)
        s = ea.setup(c)
        r["variants"] = s.get("variants", "base")
        r["restart"] = s["restart"]
        r["core_hours"] = core_hours(c)
        full.append(r)
        rows.append({"case": r["case"], "variants": r["variants"], "restart": r["restart"],
                     "t_end": r["t_last"], "window0": r["window"][0], "window1": r["window"][1],
                     **{k: r[k] for k in ("Q_in", "Q_quad", "Q_in_p", "W", "W_p", "W_v", "plate_lift_amp1",
                                          "plate_drag_mean", "gross_abs_work", "core_hours")}})
    for r in rows:
        print("%-22s %-14s Q_in %8.3f Q_quad %7.3f W %7.3f lift %8.3f drag %8.3f ch %5.1f" % (
            r["case"], r["variants"], r["Q_in"], r["Q_quad"], r["W"], r["plate_lift_amp1"], r["plate_drag_mean"], r["core_hours"]))
    if a.csv:
        with open(a.csv, "w", newline="") as fh:
            w = csv.DictWriter(fh, FIELDS)
            w.writeheader()
            for r in rows:
                w.writerow({k: (round(v, 6) if isinstance(v, float) else v) for k, v in r.items()})
    if a.json:
        a.json.write_text(json.dumps(full, indent=1))


if __name__ == "__main__":
    main()
