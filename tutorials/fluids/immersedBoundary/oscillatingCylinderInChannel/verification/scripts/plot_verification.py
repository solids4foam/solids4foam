#!/usr/bin/env python3
"""Figures of the oscillatingCylinderInChannel verification, from the
analysis results (cases.csv and history-<case>.csv, written by
'./Allverify analyse' to postProcessing/, or the tracked copies in
reference/results/) and the reference histories in reference/. No CFD is
run, and nothing is filtered or shifted.

Main figures:

    force_histories        (a, b) Cd(t) and Cl(t) of the references: the
                           body-fitted solutions A and B and the digitised
                           Wan and Turek (2007) histories; (c, d) the four
                           methods on the finest mesh run by all of them
                           (D/h = 40), against B and the digitised histories
    error_vs_resolution    normalised RMS difference of the raw Cd and Cl
                           from reference B against D/h
    cutlink_estimators     cutLink Cl from its surface traction (greedy
                           quadrature) and from its momentum exchange, same
                           runs, at D/h = 40 and 80

Further diagnostics:

    error_vs_cost          the errors against the allocated CPU hours of
                           each run (MPI ranks x wall time, on the machine
                           used; not a scaling study)
    cell_crossing          Cd(t) at full time resolution over 0.4 s,
                           showing the force fluctuations as the cylinder
                           crosses the cells

The compact histories of reference/results/ are sampled every 0.02 s
(Nyquist frequency 25 Hz). They are used for display only: the cell-
crossing fluctuations (39-314 Hz at the maximum velocity on D/h = 10-80)
are aliased in them, so no frequency content or fluctuation level is read
from them. Every error in the figures is read from cases.csv, computed by
'./Allverify analyse' from the complete histories, and cell_crossing uses
the full-resolution extracts (zoom-<case>.csv).

The style follows the staticCylinderInChannel verification: one colour,
marker and line pattern per method (Okabe-Ito colours), the references in
black and grey. Only matplotlib is needed (no LaTeX).

Usage:
    python3 scripts/plot_verification.py [--data DIR] [--out DIR]
                                         [--dpi N] [--png-only]

The previews in images/ are written with '--out images --dpi 110
--png-only'. Without --png-only, a PDF and an SVG are also written.
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

# One identity per method, the same as in staticCylinderInChannel
STYLE = {
    "cutLink": dict(color="#0072B2", marker="o", linestyle="-"),
    "ghostCell": dict(color="#D55E00", marker="s", linestyle="--"),
    "penalty": dict(color="#009E73", marker="^", linestyle="-."),
    "incremental": dict(color="#CC79A7", marker="D", linestyle=":"),
}
LABEL = {"cutLink": "cutLink", "ghostCell": "ghostCell",
         "penalty": "penalty", "incremental": "legacy incremental"}
# Draw the noisiest histories first, so that they do not hide the others
DRAW_ORDER = ("incremental", "ghostCell", "penalty", "cutLink")

# References: B a wide grey band, A black dashed, digitised open circles
REF_STYLE = {
    "B": dict(color="#b0b0b0", linewidth=3.2, solid_capstyle="butt",
              zorder=1),
    "A": dict(color="black", linewidth=0.8, linestyle=(0, (3.5, 1.8)),
              zorder=6),
    "digitised": dict(color="black", linestyle="none", marker="o",
                      markersize=2.6, markerfacecolor="none",
                      markeredgewidth=0.5, zorder=7),
}
REF_LABEL = {"B": "body-fitted B (901 120 cells)",
             "A": "body-fitted A (282 624 cells)",
             "digitised": "Wan and Turek (2007), digitised"}
COMMON_LEVEL = 3                  # finest level run by all four methods
WINDOW = (0.25, 7.5)


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
        "font.family": "serif",
        "font.serif": ["cmr10", "Computer Modern Roman", "DejaVu Serif"],
        "mathtext.fontset": "cm", "axes.formatter.use_mathtext": True,
        "axes.unicode_minus": True,
        "font.size": 9, "axes.labelsize": 9.5, "legend.fontsize": 8,
        "xtick.labelsize": 8.5, "ytick.labelsize": 8.5,
        "axes.linewidth": 0.8, "lines.linewidth": 1.0,
        "xtick.direction": "in", "ytick.direction": "in",
        "xtick.top": True, "ytick.right": True,
        "xtick.major.size": 3.5, "ytick.major.size": 3.5,
        "legend.frameon": True, "legend.fancybox": False,
        "legend.edgecolor": "black", "legend.framealpha": 1.0,
        "axes.grid": True, "grid.linestyle": (0, (3, 2)),
        "grid.linewidth": 0.4, "grid.color": "#c8c8c8",
        "savefig.dpi": dpi, "savefig.bbox": "tight",
        "svg.hashsalt": "ocv", "svg.fonttype": "path", "pdf.fonttype": 42,
    })
    return plt


PNG_ONLY = False

# Footnotes printed under the figures
SOLVER_NOTE = ("solids4foam f33b244a2 (unmodified); immersed boundary: Euler "
               "(first order); references: backward (second order)")
SAMPLING_NOTE = ("histories drawn from the compact files (every 0.02 s): "
                 "fluctuations above 25 Hz are aliased in the drawing; the "
                 "errors are computed from the complete histories")
COST_NOTE = ("allocated CPU hours = MPI ranks x wall time (1, 2, 8, 8 ranks "
             "on Levels 1-4) on one shared node; no scaling study, not a "
             "speed-up measurement")


def footnote(fig, *lines, y=-0.01):
    fig.text(0.01, y, "\n".join(lines), fontsize=6.5, va="top",
             ha="left", color="0.3")


def panel(ax, text):
    ax.text(0.015, 0.97, text, transform=ax.transAxes, va="top", ha="left",
            fontsize=8, zorder=20,
            bbox=dict(boxstyle="square,pad=0.15", fc="white", ec="none"))


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
        fig.savefig(out / f"{name}.svg", metadata={"Date": None})
    fig.savefig(out / f"{name}.png", metadata={"Software": None})
    written.append(out / name)


def at_level(cases: list, method: str, level: int):
    """The run of a method on a level with the time step halved with the
    cell size (the R1-R3, R5 runs)"""
    for c in cases:
        if (c["method"] == method and int(c["level"]) == level
                and abs(c["deltaT"] - ocv.LEVEL_DT[level]) < 1e-12):
            return c
    return None


def window_limits(series, pad=0.06, head=0.0):
    """y range of every series for 0.25 <= t <= 8 s (the impulsive start
    is excluded), so that nothing drawn after the start is clipped"""
    lo = min(min(y for t, y in zip(*s) if WINDOW[0] <= t <= 8)
             for s in series)
    hi = max(max(y for t, y in zip(*s) if WINDOW[0] <= t <= 8)
             for s in series)
    d = hi - lo
    return lo - pad*d, hi + (pad + head)*d


def plot_refs(ax, refs, coeff, names):
    hs = []
    for name in names:
        if name in refs:
            hs += ax.plot(*refs[name][coeff], label=REF_LABEL[name],
                          **REF_STYLE[name])
    return hs


def force_histories(plt, cases, data, refs, out, written):
    """2x2: the references, then the four methods at D/h = 40"""
    fig, axes = plt.subplots(2, 2, figsize=(7.4, 5.4), sharex=True)
    runs = {m: at_level(cases, m, COMMON_LEVEL) for m in ocv.METHODS}
    for j, coeff in enumerate(("Cd", "Cl")):
        a, b = axes[0, j], axes[1, j]
        plot_refs(a, refs, coeff, ("B", "A", "digitised"))
        a.set_ylim(*window_limits([refs[n][coeff] for n in refs]))
        plot_refs(b, refs, coeff, ("B", "digitised"))
        series = [refs[n][coeff] for n in ("B", "digitised") if n in refs]
        for method in DRAW_ORDER:
            c = runs[method]
            if c is None:
                continue
            h = read_history(data, c["case"])
            s = STYLE[method]
            b.plot(h["t"], h[f"{coeff}_total"], color=s["color"],
                   linewidth=0.7, zorder=2 + DRAW_ORDER.index(method))
            series.append((h["t"], h[f"{coeff}_total"]))
        b.set_ylim(*window_limits(series))
        for ax in (a, b):
            ax.axvspan(0, WINDOW[0], color="0.9", linewidth=0, zorder=0)
            ax.set_xlim(0, 8)
            ax.set_ylabel(f"$C_{coeff[1]}$")
        b.set_xlabel("$t$ (s)")
        panel(a, f"({'ab'[j]}) references")
        panel(b, f"({'cd'[j]}) four methods, $D/h$ = 40")
    handles = [plt.Line2D([], [], **REF_STYLE[n]) for n in
               ("B", "A", "digitised") if n in refs]
    labels = [REF_LABEL[n] for n in ("B", "A", "digitised") if n in refs]
    for method in ocv.METHODS:
        s = STYLE[method]
        handles.append(plt.Line2D([], [], color=s["color"],
                                  linestyle=s["linestyle"], linewidth=1.2))
        labels.append(LABEL[method] + (" (surface traction, greedy "
                                       "quadrature)" if method == "cutLink"
                                       else ""))
    handles.append(plt.Line2D([], [], color="0.9", linewidth=6))
    labels.append("impulsive start ($t < 0.25$ s, off scale, not in "
                  "the errors)")
    fig.tight_layout(rect=(0, 0.13, 1, 1))
    fig.legend(handles, labels, ncol=3, loc="lower center",
               bbox_to_anchor=(0.5, 0.0))
    footnote(fig, SOLVER_NOTE, "D/h = 40, time step 0.000625 s, for all "
             "four methods (the finest mesh run by all of them)",
             SAMPLING_NOTE)
    save(fig, out, "force_histories", written)
    plt.close(fig)


def error_vs(plt, cases, out, written, xkey, xlabel, name, logx,
             spread=None):
    fig, axes = plt.subplots(1, 2, figsize=(7.4, 3.0))
    for k, (ax, coeff) in enumerate(zip(axes, ("Cd", "Cl"))):
        if spread and coeff in spread:
            ax.axhspan(1e-3, 100*spread[coeff], color="0.88", linewidth=0,
                       zorder=0, label="below the A-B difference")
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
                        markersize=4.5 if dt_level else 6,
                        label=LABEL[method] if dt_level else None)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(f"$C_{coeff[1]}$: normalised RMS difference "
                      "from B (%)")
        ax.set_yscale("log")
        ax.set_ylim(0.3, 300)
        if logx:
            ax.set_xscale("log")
        if xkey == "cells_per_D":
            ax.set_xlim(8, 100)
            ax.set_xticks([10, 20, 40, 80])
            ax.set_xticklabels(["10", "20", "40", "80"])
            ax.minorticks_off()
        panel(ax, f"({'ab'[k]}) {'drag' if coeff == 'Cd' else 'lift'}")
    handles, labels = axes[0].get_legend_handles_labels()
    handles.append(plt.Line2D([], [], color="0.3", marker="o",
                              markerfacecolor="white", linestyle="none",
                              markersize=6))
    labels.append("open: D/h = 20 at the D/h = 40 time step (R4)")
    fig.tight_layout(rect=(0, 0.12, 1, 1))
    fig.legend(handles, labels, ncol=3, loc="lower center",
               bbox_to_anchor=(0.5, 0.0))
    notes = [SOLVER_NOTE, "raw forces; filled markers: time step halved with "
             "the cell size (space and time refined together); legacy "
             "incremental not run on D/h = 80"]
    if xkey == "cpu_hours":
        notes.append(COST_NOTE)
    footnote(fig, *notes)
    save(fig, out, name, written)
    plt.close(fig)


def cutlink_estimators(plt, cases, data, refs, out, written):
    """Cl from the two cutLink estimates of the same runs at D/h = 40, where
    the default (greedy) thinning of the surface-traction points acts, and
    at D/h = 80, where every candidate point is kept"""
    runs = [c for c in (at_level(cases, "cutLink", 3),
                        at_level(cases, "cutLink", 4)) if c is not None]
    if not runs:
        return
    fig, axes = plt.subplots(len(runs), 1, figsize=(7.4, 2.0*len(runs) + 0.9),
                             sharex=True, squeeze=False)
    for ax, c in zip(axes[:, 0], runs):
        h = read_history(data, c["case"])
        if "Cl_momentumExchange" not in h:
            continue
        ax.axvspan(0, WINDOW[0], color="0.9", linewidth=0, zorder=0)
        plot_refs(ax, refs, "Cl", ("B",))
        ax.plot(h["t"], h["Cl_surfaceTraction"], color=STYLE["cutLink"][
            "color"], linewidth=1.0, zorder=4)
        ax.plot(h["t"], h["Cl_momentumExchange"], color="black",
                linestyle=(0, (4, 1.6)), linewidth=0.9, zorder=5)
        ax.set_ylabel("$C_l$")
        ax.set_xlim(0, 8)
        series = [refs["B"]["Cl"], (h["t"], h["Cl_surfaceTraction"]),
                  (h["t"], h["Cl_momentumExchange"])]
        ax.set_ylim(*window_limits(series, head=0.45))
        panel(ax, f"({'ab'[runs.index(c)]}) $D/h$ = {int(c['cells_per_D'])}"
              f", $\\Delta t$ = {c['deltaT']:g} s\n"
              f"Cl vs B: {100*c['Cl_nrms_B']:.1f}% (surface traction), "
              f"{100*c['Cl_ME_nrms_B']:.1f}% (momentum exchange)\n"
              f"Cd vs B: {100*c['Cd_nrms_B']:.2f}% (surface traction)")
    axes[-1, 0].set_xlabel("$t$ (s)")
    handles = [plt.Line2D([], [], **REF_STYLE["B"]),
               plt.Line2D([], [], color=STYLE["cutLink"]["color"],
                          linewidth=1.0),
               plt.Line2D([], [], color="black", linestyle=(0, (4, 1.6)),
                          linewidth=0.9)]
    labels = [REF_LABEL["B"], "cutLink surface traction (greedy "
              "quadrature; the estimator of the tracked results)",
              "cutLink momentum exchange (same run; reported separately, "
              "never added)"]
    fig.tight_layout(rect=(0, 0.14, 1, 1))
    fig.legend(handles, labels, ncol=1, loc="lower center",
               bbox_to_anchor=(0.5, 0.0))
    footnote(fig, SOLVER_NOTE, "the greedy quadrature keeps 128 of the 256 "
             "candidate surface points on D/h = 40 and all 256 on D/h = 80; "
             + SAMPLING_NOTE.split(":")[0])
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


def finest(cases: list, method: str):
    pick = [c for c in cases if c["method"] == method]
    if not pick:
        return None
    return max(pick, key=lambda c: (c["level"], -c["deltaT"]))


def cell_crossing(plt, cases, data, refs, out, written, window=(1.8, 2.2)):
    fig, ax = plt.subplots(figsize=(7.4, 3.0))
    if "B" in refs:
        t, y = refs["B"]["Cd"]
        sel = [(a, b) for a, b in zip(t, y) if window[0] <= a <= window[1]]
        ax.plot([a for a, _ in sel], [b for _, b in sel],
                label=REF_LABEL["B"], **REF_STYLE["B"])
    for method in DRAW_ORDER:
        c = finest(cases, method)
        if c is None:
            continue
        h = read_zoom(data, c["case"])
        sel = [(a, b) for a, b in zip(h["t"], h["Cd_total"])
               if window[0] <= a <= window[1]]
        s = STYLE[method]
        ax.plot([a for a, _ in sel], [b for _, b in sel], color=s["color"],
                linestyle=s["linestyle"], linewidth=0.8,
                label=f"{LABEL[method]}, D/h = {int(c['cells_per_D'])}")
    ax.set_xlabel("$t$ (s)")
    ax.set_ylabel("$C_d$")
    ax.set_xlim(*window)
    handles, labels = ax.get_legend_handles_labels()
    fig.tight_layout(rect=(0, 0.12, 1, 1))
    fig.legend(handles, labels, ncol=3, loc="lower center",
               bbox_to_anchor=(0.5, 0.0))
    footnote(fig, "every time step (full-resolution extracts zoom-*.csv); "
             "each method on its finest run (legacy incremental: D/h = 40)")
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
        print(f"Wrote {path}" + (".png" if PNG_ONLY else
                                 ".pdf, .svg and .png"))


if __name__ == "__main__":
    main()
