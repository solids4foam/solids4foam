#!/usr/bin/env python3
"""Figures of the oscillatingCylinderInChannel verification, from the
analysis results (cases.csv and history-<case>.csv, written by
'./Allverify analyse' to postProcessing/, or the tracked copies in
reference/results/) and the reference histories in reference/. No CFD is
run.

    force_histories        Cd(t) and Cl(t) of the methods on the finest
                           level run, against both body-fitted references
                           and the digitised literature histories
    error_vs_resolution    normalised RMS difference from reference B
                           against D/h, for each time step
    error_vs_cost          the same against the allocated CPU hours of
                           each run (MPI ranks x wall time, on the machine
                           used; not a scaling study)
    cutlink_estimators     cutLink surface traction and momentum exchange
                           from the same run
    cell_crossing          a zoom of Cd(t) showing the force fluctuations as
                           the cylinder crosses the cells

Usage:
    python3 scripts/plot_verification.py [--data DIR] [--out DIR]
                                         [--dpi N] [--png-only]

The previews in images/ are written with '--out images --dpi 110
--png-only'.
"""

import argparse
import csv
import importlib.util
import math
import sys
from pathlib import Path

VERIFICATION_DIR = Path(__file__).resolve().parent.parent
spec = importlib.util.spec_from_file_location(
    "ocv", Path(__file__).resolve().parent
    / "oscillating_cylinder_verification.py")
ocv = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ocv)

STYLE = {
    "cutLink": dict(color="#0072B2", marker="o", linestyle="-"),
    "ghostCell": dict(color="#E69F00", marker="D", linestyle=":"),
    "penalty": dict(color="#D55E00", marker="s", linestyle="--"),
    "incremental": dict(color="#009E73", marker="^", linestyle="-."),
}
LABEL = {"cutLink": "cutLink", "ghostCell": "ghostCell",
         "penalty": "penalty", "incremental": "legacy incremental"}


def read_cases(data: Path) -> list:
    with (data / "cases.csv").open() as handle:
        rows = list(csv.DictReader(handle))
    for r in rows:
        for k, v in list(r.items()):
            try:
                r[k] = float(v)
            except (TypeError, ValueError):
                pass
    return [r for r in rows if r.get("integrity") == "PASS"
            and abs(r["end_time"] - ocv.END_TIME) < 1e-9]


def read_history(data: Path, case: str) -> dict:
    with (data / f"history-{case}.csv").open() as handle:
        rows = list(csv.reader(handle))
    keys = rows[0]
    cols = list(zip(*[[float(x) for x in r] for r in rows[1:]]))
    return dict(zip(keys, cols))


def setup(dpi=200):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({
        "font.family": "serif", "font.serif": ["STIXGeneral", "DejaVu Serif"],
        "mathtext.fontset": "stix", "font.size": 9, "axes.linewidth": 0.7,
        "xtick.direction": "in", "ytick.direction": "in",
        "xtick.top": True, "ytick.right": True, "legend.frameon": False,
        "legend.fontsize": 8, "lines.linewidth": 1.1,
        "savefig.dpi": dpi, "savefig.bbox": "tight",
        "svg.hashsalt": "ocv", "pdf.fonttype": 42,
    })
    return plt


PNG_ONLY = False

# Footnotes printed under the figures
SOLVER_NOTE = ("solids4foam f33b244a2 (unmodified); immersed boundary: Euler "
               "(first order); references: backward (second order)")
RAW_NOTE = ("raw force histories; < 5 Hz values are a filtered diagnostic, "
            "not a force estimator")
COST_NOTE = ("allocated CPU hours = MPI ranks x wall time (1, 2, 8, 8 ranks "
             "on Levels 1-4) on one shared node; no scaling study, not a "
             "speed-up measurement")


def footnote(fig, *lines):
    fig.text(0.01, -0.02, "\n".join(lines), fontsize=6.5, va="top",
             ha="left", color="0.3")


def reference_spread(refs: dict) -> dict:
    """Normalised RMS difference between references A and B, per
    coefficient (the resolution of the comparison)"""
    if "A" not in refs or "B" not in refs:
        return {}
    return {c: ocv.rms_error(*refs["A"][c], *refs["B"][c])["nrms"]
            for c in ("Cd", "Cl")}


def save(fig, out: Path, name: str, written: list) -> None:
    if not PNG_ONLY:
        fig.savefig(out / f"{name}.pdf", metadata={"CreationDate": None})
    fig.savefig(out / f"{name}.png", metadata={"Software": None})
    written.append(out / name)


def finest(cases: list, method: str):
    pick = [c for c in cases if c["method"] == method]
    if not pick:
        return None
    return max(pick, key=lambda c: (c["level"], -c["deltaT"]))


def force_histories(plt, cases, data, refs, out, written):
    fig, axes = plt.subplots(2, 1, figsize=(7.0, 5.2), sharex=True)
    for ax, coeff in zip(axes, ("Cd", "Cl")):
        if "B" in refs:
            ax.plot(*refs["B"][coeff], color="black", linewidth=1.6,
                    label="body-fitted B (901k cells)")
        if "A" in refs:
            ax.plot(*refs["A"][coeff], color="0.55", linestyle=":",
                    linewidth=1.4, label="body-fitted A (282k cells)")
        if "digitised" in refs:
            ax.plot(*refs["digitised"][coeff], linestyle="none", marker=".",
                    markersize=2.5, color="0.35",
                    label="digitised literature (secondary)")
        for method in ocv.METHODS:
            c = finest(cases, method)
            if c is None:
                continue
            h = read_history(data, c["case"])
            s = STYLE[method]
            ax.plot(h["t"], h[f"{coeff}_total"], color=s["color"],
                    linestyle=s["linestyle"], linewidth=0.9,
                    label=f"{LABEL[method]}, L{int(c['level'])}, "
                    f"$\\Delta t$ = {c['deltaT']:g} s")
        ax.set_ylabel(f"$C_{coeff[1]}$")
        ax.set_xlim(0, 8)
    axes[0].set_ylim(-4.5, 4.5)
    axes[1].set_ylim(-0.06, 0.06)
    axes[1].set_xlabel("$t$ (s)")
    axes[0].legend(ncol=2, loc="upper center", bbox_to_anchor=(0.5, 1.42))
    fig.tight_layout()
    footnote(fig, SOLVER_NOTE, "cutLink: surface traction with the greedy "
             "quadrature; legacy incremental not run on Level 4; histories "
             "sampled every 0.02 s")
    save(fig, out, "force_histories", written)
    plt.close(fig)


def error_vs(plt, cases, out, written, xkey, xlabel, name, logx,
             spread=None):
    fig, axes = plt.subplots(1, 2, figsize=(7.0, 2.8))
    for ax, coeff in zip(axes, ("Cd", "Cl")):
        if spread and coeff in spread:
            ax.axhline(100*spread[coeff], color="0.5", linestyle="--",
                       linewidth=0.8, label="A-B difference")
        for method in ocv.METHODS:
            s = STYLE[method]
            for dt_level in (True, False):
                pts = sorted((c[xkey], c[f"{coeff}_nrms_B"]) for c in cases
                             if c["method"] == method
                             and (abs(c["deltaT"] - ocv.LEVEL_DT[int(
                                 c["level"])]) < 1e-12) == dt_level)
                if not pts:
                    continue
                ax.plot([p[0] for p in pts], [100*p[1] for p in pts],
                        color=s["color"], marker=s["marker"],
                        linestyle=s["linestyle"] if dt_level else "none",
                        markerfacecolor=s["color"] if dt_level else "white",
                        label=f"{LABEL[method]}"
                        + ("" if dt_level else ", other $\\Delta t$"))
        ax.set_xlabel(xlabel)
        ax.set_ylabel(f"normalised RMS difference in $C_{coeff[1]}$ (%)")
        ax.set_yscale("log")
        if logx:
            ax.set_xscale("log")
        if xkey == "cells_per_D":
            ax.set_xticks([10, 20, 40, 80])
            ax.set_xticklabels(["10", "20", "40", "80"])
            ax.minorticks_off()
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, ncol=4, loc="upper center",
               bbox_to_anchor=(0.5, 0.0))
    fig.tight_layout()
    notes = [SOLVER_NOTE, RAW_NOTE + "; legacy incremental not run on "
             "Level 4"]
    if xkey == "cpu_hours":
        notes.append(COST_NOTE)
    fig.text(0.01, -0.27, "\n".join(notes), fontsize=6.5, va="top",
             ha="left", color="0.3")
    save(fig, out, name, written)
    plt.close(fig)


def cutlink_estimators(plt, cases, data, refs, out, written):
    """Cl from the two cutLink estimates of the same runs at D/h = 40, where
    the default (greedy) thinning of the surface-traction points acts, and
    at the finest level run (at D/h = 80 every candidate point is kept)"""
    runs = sorted((c for c in cases if c["method"] == "cutLink"
                   and int(c["level"]) >= 3
                   and abs(c["deltaT"] - ocv.LEVEL_DT[int(c["level"])])
                   < 1e-12), key=lambda c: c["level"])
    if not runs:
        return
    fig, axes = plt.subplots(len(runs), 1, figsize=(7.0, 2.2*len(runs) + 0.4),
                             sharex=True, squeeze=False)
    for ax, c in zip(axes[:, 0], runs):
        h = read_history(data, c["case"])
        if "Cl_momentumExchange" not in h:
            continue
        if "B" in refs:
            ax.plot(*refs["B"]["Cl"], color="black", linewidth=1.6,
                    label="body-fitted B")
        ax.plot(h["t"], h["Cl_surfaceTraction"], color="#0072B2",
                label="cutLink surface traction")
        ax.plot(h["t"], h["Cl_momentumExchange"], color="#CC79A7",
                linestyle="--", label="cutLink momentum exchange")
        ax.set_ylabel("$C_l$")
        ax.set_xlim(0, 8)
        ax.set_ylim(-0.08, 0.08)
        ax.set_title(f"$D/h$ = {int(c['cells_per_D'])}, $\\Delta t$ = "
                     f"{c['deltaT']:g} s", fontsize=9, loc="left")
    axes[-1, 0].set_xlabel("$t$ (s)")
    axes[0, 0].legend(ncol=3, loc="upper center",
                      bbox_to_anchor=(0.5, 1.45))
    fig.tight_layout()
    save(fig, out, "cutlink_estimators", written)
    plt.close(fig)


def read_zoom(data: Path, case: str) -> dict:
    """The full-resolution window of a history if exported, else the
    history"""
    path = data / f"zoom-{case}.csv"
    if not path.exists():
        return read_history(data, case)
    with path.open() as handle:
        rows = list(csv.reader(handle))
    cols = list(zip(*[[float(x) for x in r] for r in rows[1:]]))
    return dict(zip(rows[0], cols))


def cell_crossing(plt, cases, data, refs, out, written, window=(1.8, 2.2)):
    fig, ax = plt.subplots(figsize=(7.0, 3.0))
    if "B" in refs:
        t, y = refs["B"]["Cd"]
        sel = [(a, b) for a, b in zip(t, y) if window[0] <= a <= window[1]]
        ax.plot([a for a, _ in sel], [b for _, b in sel], color="black",
                linewidth=1.6, label="body-fitted B")
    for method in ocv.METHODS:
        c = finest(cases, method)
        if c is None:
            continue
        h = read_zoom(data, c["case"])
        sel = [(a, b) for a, b in zip(h["t"], h["Cd_total"])
               if window[0] <= a <= window[1]]
        s = STYLE[method]
        ax.plot([a for a, _ in sel], [b for _, b in sel], color=s["color"],
                linestyle=s["linestyle"], linewidth=0.8,
                label=f"{LABEL[method]}, L{int(c['level'])}")
    ax.set_xlabel("$t$ (s)")
    ax.set_ylabel("$C_d$")
    ax.set_xlim(*window)
    ax.legend(ncol=5, loc="upper center", bbox_to_anchor=(0.5, 1.18))
    fig.tight_layout()
    save(fig, out, "cell_crossing", written)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    default = VERIFICATION_DIR / "reference" / "results"
    if not (default / "cases.csv").exists():
        default = VERIFICATION_DIR / "postProcessing"
    parser.add_argument("--data", default=str(default),
                        help="directory of cases.csv and history-*.csv")
    parser.add_argument("--out", default=str(VERIFICATION_DIR
                                             / "postProcessing" / "figures"))
    parser.add_argument("--dpi", type=int, default=200,
                        help="resolution of the PNG figures")
    parser.add_argument("--png-only", action="store_true",
                        help="write the PNG figures only (e.g. previews)")
    args = parser.parse_args()
    global PNG_ONLY
    PNG_ONLY = args.png_only
    data, out = Path(args.data), Path(args.out)
    if not (data / "cases.csv").exists():
        sys.exit(f"No cases.csv in {data}: run ./Allverify analyse first")
    out.mkdir(parents=True, exist_ok=True)
    cases = read_cases(data)
    if not cases:
        sys.exit("No complete 8 s cases with numerical integrity PASS")
    refs = ocv.references()
    plt = setup(args.dpi)
    written = []
    force_histories(plt, cases, data, refs, out, written)
    spread = reference_spread(refs)
    error_vs(plt, cases, out, written, "cells_per_D", "$D/h$",
             "error_vs_resolution", True, spread)
    if any(math.isfinite(c.get("cpu_hours", math.nan)) for c in cases):
        error_vs(plt, cases, out, written, "cpu_hours",
                 "allocated CPU hours (MPI ranks $\\times$ wall time)",
                 "error_vs_cost", True, spread)
    cutlink_estimators(plt, cases, data, refs, out, written)
    cell_crossing(plt, cases, data, refs, out, written)
    for path in written:
        print(f"Wrote {path}.pdf and .png")


if __name__ == "__main__":
    main()
