#!/usr/bin/env python3
"""Figures of the oscillatingWallStokesLayer verification from the tracked
results in reference/results/ (or those of './Allverify analyse' with
--data). No CFD is run and nothing is filtered or shifted.

    shear_histories      wall shear stress against time over the fourth
                         period, Levels 2 and 4, with the periodic solution
                         and the exact initial-value solution
    shear_vs_mesh        amplitude and phase error of the shear against the
                         mesh level (against the initial-value solution)
    velocity_profiles    velocity profiles at t/T = 0, 0.25, 0.5 and 0.75 of
                         the fourth period, Level 3, with the exact solutions
    velocity_convergence L1, L2 and Linf norms of the velocity difference
                         against the cell width
    shear_noise          fluctuation of the shear about its mean and three
                         harmonics, against the mesh
    transient            the difference between the exact initial-value
                         solution and the periodic solution in the shear,
                         and the numerical shear (forcing estimator) minus
                         the periodic solution
    time_step            amplitude error and velocity norm against the time
                         step on Levels 2-4 (separates the step from the mesh)
    cost                 observed cost against the number of cells

Style of the immersed-boundary verifications: one colour, marker and line
pattern per method (Okabe-Ito colours), references in black and grey. Only
matplotlib is needed (no LaTeX).

Usage:
    python3 scripts/plot_verification.py [--data DIR] [--out DIR]
                                         [--dpi N] [--png-only]

The previews in images/ are written with '--out images --dpi 110
--png-only'. Without --png-only, a PDF is also written.
"""

import argparse
import csv
import math
import sys
from pathlib import Path

VERIFICATION_DIR = Path(__file__).resolve().parent.parent
REFERENCE_DIR = VERIFICATION_DIR / "reference"

KEYS = (("cutLink", "ME"), ("cutLink", "ST"), ("penalty", ""),
        ("incremental", ""), ("ghostCell", ""))
COLOR = {("cutLink", "ME"): "#0072B2", ("cutLink", "ST"): "#56B4E9",
         ("cutLink", "FO"): "#7F7F7F", ("penalty", ""): "#009E73",
         ("incremental", ""): "#CC79A7", ("ghostCell", ""): "#D55E00"}
MARKER = {("cutLink", "ME"): "o", ("cutLink", "ST"): "v",
          ("cutLink", "FO"): "x", ("penalty", ""): "^",
          ("incremental", ""): "D", ("ghostCell", ""): "s"}
DASH = {("cutLink", "ME"): "-", ("cutLink", "ST"): (0, (3, 1.2)),
        ("cutLink", "FO"): (0, (1, 1)), ("penalty", ""): (0, (5, 1.5, 1.2, 1.5)),
        ("incremental", ""): (0, (1.2, 1.2)), ("ghostCell", ""): (0, (4, 1.6))}
LABEL = {("cutLink", "ME"): "cutLink, momentum exchange",
         ("cutLink", "ST"): "cutLink, surface traction",
         ("cutLink", "FO"): "cutLink, forcing",
         ("penalty", ""): "penalty", ("incremental", ""): "legacy incremental",
         ("ghostCell", ""): "ghostCell (exploratory)"}
EXACT = dict(color="black", lw=0.95, zorder=1.5)
EXACT_IVP = dict(color="0.5", lw=0.95, ls=(0, (3, 2)), zorder=1.4)
TUTORIAL_DT = "0.001"


def read_csv(path: Path) -> list:
    with path.open() as handle:
        return list(csv.DictReader(handle))


def num(row: dict, key: str) -> float:
    return float(row[key]) if row.get(key, "") not in ("", "nan") else float("nan")


def style(plt):
    plt.rcParams.update({
        "font.family": "serif", "font.size": 8.5, "axes.labelsize": 9,
        "axes.titlesize": 9, "legend.fontsize": 7.2, "axes.linewidth": 0.8,
        "xtick.direction": "in", "ytick.direction": "in",
        "xtick.top": True, "ytick.right": True,
        "legend.frameon": True, "legend.edgecolor": "black",
        "legend.fancybox": False, "legend.framealpha": 1.0,
        "savefig.bbox": "tight"})


def panel(ax, label):
    ax.text(0.02, 0.97, label, transform=ax.transAxes, va="top", ha="left",
            fontweight="bold")
    ax.grid(True, ls="--", lw=0.4, color="0.75")


def key(r):
    return (r["method"], r["estimator"])


def tutorial(rows):
    return [r for r in rows if r["dt"] in ("", TUTORIAL_DT)]


def series(rows, k, ykey, levels=(1, 2, 3, 4)):
    sel = sorted((int(r["level"]), num(r, ykey)) for r in tutorial(rows)
                 if key(r) == k and int(r["level"]) in levels)
    return [a for a, _ in sel], [b for _, b in sel]


def h_of(level):
    return 0.5 / (20 * 2 ** (level - 1))


def save(plt, fig, out, name, dpi, png_only):
    out.mkdir(parents=True, exist_ok=True)
    fig.savefig(out / f"{name}.png", dpi=dpi)
    if not png_only:
        fig.savefig(out / f"{name}.pdf")
    plt.close(fig)


def vid(k):
    return k[0] + ("-" + k[1] if k[0] == "cutLink" else "")


def shear_histories(plt, data, out, dpi, png_only):
    fig, axes = plt.subplots(1, 2, figsize=(7.4, 3.1), sharey=True)
    for ax, level, lab in zip(axes, (2, 4), "ab"):
        rows = read_csv(data / f"shear-L{level}.csv")
        t = [num(r, "time") for r in rows]
        ax.plot(t, [1e3 * num(r, "exact_periodic") for r in rows], **EXACT,
                label="exact, periodic")
        ax.plot(t, [1e3 * num(r, "exact_ivp") for r in rows], **EXACT_IVP,
                label="exact, initial-value")
        for k in KEYS:
            if vid(k) not in rows[0]:
                continue
            ax.plot(t, [1e3 * num(r, vid(k)) for r in rows], color=COLOR[k],
                    ls=DASH[k], lw=0.8, label=LABEL[k])
        ax.set_xlabel("Time (s)")
        ax.set_title(f"Level {level}: {2 * 2 ** (level - 1) * 10 // 4} cells "
                     f"across $\\delta$" if False else f"Level {level}",
                     fontsize=8.5)
        panel(ax, lab)
    axes[0].set_ylabel(r"Wall shear stress, $\tau_w$ (mPa)")
    handles, labels = axes[1].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=4,
               bbox_to_anchor=(0.5, -0.17))
    save(plt, fig, out, "shear_histories", dpi, png_only)


def shear_vs_mesh(plt, data, out, dpi, png_only):
    rows = read_csv(data / "cases.csv")
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.1))
    for ax, ykey, ylabel, lab in (
            (axes[0], "shear_amp_err_ivp", "Amplitude error (%)", "a"),
            (axes[1], "shear_phase_err_deg_ivp", "Phase error (deg)", "b")):
        for k in KEYS + (("cutLink", "FO"),):
            lv, y = series(rows, k, ykey)
            scale = 100.0 if "amp" in ykey else 1.0
            ax.plot([2 ** (a - 1) * 2 for a in lv], [scale * x for x in y],
                    color=COLOR[k], ls=DASH[k], marker=MARKER[k], ms=4.5,
                    lw=0.9, label=LABEL[k], mfc=COLOR[k])
        ax.axhline(0, **EXACT)
        ax.set_xscale("log", base=2)
        ax.set_xticks([2, 4, 8, 16])
        ax.set_xticklabels(["2", "4", "8", "16"])
        ax.set_xlabel(r"Cells across $\delta$")
        ax.set_ylabel(ylabel)
        ax.set_ylim(-23, 8) if "amp" in ykey else ax.set_ylim(-23, 8)
        panel(ax, lab)
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=3,
               bbox_to_anchor=(0.5, -0.16))
    save(plt, fig, out, "shear_vs_mesh", dpi, png_only)


def velocity_profiles(plt, data, out, dpi, png_only):
    level = 3
    rows = read_csv(data / f"profiles-L{level}.csv")
    fig, axes = plt.subplots(2, 4, figsize=(7.6, 5.6), sharey="row")
    names = [n for n in ("cutLink", "penalty", "incremental", "ghostCell")
             if n in rows[0]]
    kmap = {"cutLink": ("cutLink", "ME"), "penalty": ("penalty", ""),
            "incremental": ("incremental", ""), "ghostCell": ("ghostCell", "")}
    for col, (phase, lab) in enumerate(zip(("0.0", "0.25", "0.5", "0.75"),
                                           "abcd")):
        sel = [r for r in rows if r["phase"] == phase]
        s = [num(r, "s") for r in sel]
        ax = axes[0][col]
        ax.plot([1e3 * num(r, "exact_periodic") for r in sel], s, **EXACT,
                label="exact, periodic")
        ax.plot([1e3 * num(r, "exact_ivp") for r in sel], s, **EXACT_IVP,
                label="exact, initial-value")
        for n in names:
            k = kmap[n]
            ax.plot([1e3 * num(r, n) for r in sel], s, color=COLOR[k],
                    ls="none", marker=MARKER[k], ms=2.6, mfc="none", mew=0.6,
                    label=LABEL[k] if n != "cutLink" else "cutLink")
        ax.set_title(f"$t/T$ = {phase}", fontsize=8.5)
        ax.set_xlabel(r"$u_x$ (mm/s)")
        ax.set_ylim(0, 0.3)
        panel(ax, lab)
        ax = axes[1][col]
        for n in names:
            k = kmap[n]
            ax.plot([1e3 * (num(r, n) - num(r, "exact_ivp")) for r in sel], s,
                    color=COLOR[k], ls=DASH[k], marker=MARKER[k], ms=2.4,
                    lw=0.7)
        ax.plot([1e3 * (num(r, "exact_periodic") - num(r, "exact_ivp"))
                 for r in sel], s, **EXACT)
        ax.axvline(0, color="0.6", lw=0.5)
        ax.set_xlabel("$u_x$ - exact IVP (mm/s)")
        ax.set_ylim(0, 0.3)
        panel(ax, "")
    axes[0][0].set_ylabel("Distance from the wall (m)")
    axes[1][0].set_ylabel("Distance from the wall (m)")
    axes[0][0].legend(loc="upper left", fontsize=5.8)
    save(plt, fig, out, "velocity_profiles", dpi, png_only)


def velocity_convergence(plt, data, out, dpi, png_only):
    rows = read_csv(data / "cases.csv")
    fig, axes = plt.subplots(1, 3, figsize=(7.4, 2.9), sharey=True)
    for ax, norm, lab in zip(axes, ("u_L1_all_ivp", "u_L2_all_ivp",
                                    "u_Linf_all_ivp"), "abc"):
        for k in KEYS[:1] + KEYS[2:]:
            lv, e = series(rows, k, norm)
            ax.plot([h_of(a) * 1e3 for a in lv], e, color=COLOR[k],
                    ls=DASH[k], marker=MARKER[k], ms=4.2, lw=0.9,
                    label=LABEL[k] if k != ("cutLink", "ME") else "cutLink")
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("Cell width (mm)")
        ax.set_title({"u_L1_all_ivp": "L1", "u_L2_all_ivp": "L2",
                      "u_Linf_all_ivp": "Linf"}[norm], fontsize=8.5)
        panel(ax, lab)
    axes[0].set_ylabel(r"norm of $u_x - u_{x,\mathrm{exact}}$ over $U_0$")
    axes[0].legend(loc="lower right", fontsize=6.5)
    save(plt, fig, out, "velocity_convergence", dpi, png_only)


def shear_noise(plt, data, out, dpi, png_only):
    rows = read_csv(data / "cases.csv")
    fig, ax = plt.subplots(figsize=(3.9, 3.0))
    for k in KEYS:
        lv, y = series(rows, k, "shear_noise")
        ax.plot([2 ** (a - 1) * 2 for a in lv], y, color=COLOR[k],
                ls=DASH[k], marker=MARKER[k], ms=4.2, lw=0.9, label=LABEL[k])
    ax.set_xscale("log", base=2)
    ax.set_yscale("log")
    ax.set_xticks([2, 4, 8, 16])
    ax.set_xticklabels(["2", "4", "8", "16"])
    ax.set_xlabel(r"Cells across $\delta$")
    ax.set_ylabel("RMS fluctuation of the shear over the exact RMS")
    ax.legend(loc="upper right", fontsize=6.3)
    panel(ax, "")
    save(plt, fig, out, "shear_noise", dpi, png_only)


def transient(plt, data, out, dpi, png_only):
    rows = read_csv(data / "shear-L4.csv")
    t = [num(r, "time") for r in rows]
    fig, ax = plt.subplots(figsize=(4.6, 3.1))
    d = [1e3 * (num(r, "exact_ivp") - num(r, "exact_periodic")) for r in rows]
    ax.plot(t, d, **EXACT, label="exact initial-value minus periodic")
    k = ("cutLink", "FO")
    if vid(k) in rows[0]:
        ax.plot(t, [1e3 * (num(r, vid(k)) - num(r, "exact_periodic"))
                    for r in rows], color=COLOR[k], ls="none", marker="x",
                ms=3, label="cutLink forcing, Level 4, minus periodic")
        ax.plot(t, [1e3 * (num(r, vid(k)) - num(r, "exact_ivp"))
                    for r in rows], color="#D55E00", ls="none", marker=".",
                ms=3, label="cutLink forcing, Level 4, minus initial-value")
    ax.axhline(0, color="0.6", lw=0.5)
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Difference in the wall shear (mPa)")
    ax.set_ylim(-0.03, 0.13)
    ax.legend(loc="upper right", fontsize=6.3)
    panel(ax, "")
    save(plt, fig, out, "transient", dpi, png_only)


def time_step(plt, data, out, dpi, png_only):
    rows = read_csv(data / "cases.csv")
    fig, axes = plt.subplots(2, 3, figsize=(7.4, 4.8), sharex=True)
    methods = (("cutLink", "ME"), ("cutLink", "ST"), ("penalty", ""),
               ("incremental", ""), ("ghostCell", ""))
    for col, level in enumerate((2, 3, 4)):
        for row, (ykey, scale, ylabel) in enumerate((
                ("shear_amp_err_ivp", 100.0, "Amplitude error (%)"),
                ("u_L2_all_ivp", 1.0, "Velocity L2 over $U_0$"))):
            ax = axes[row][col]
            for k in methods:
                pts = sorted((num(r, "dt"), scale * num(r, ykey))
                             for r in rows if key(r) == k
                             and int(r["level"]) == level)
                if not pts:
                    continue
                ax.plot([p[0] for p in pts], [p[1] for p in pts],
                        color=COLOR[k], ls=DASH[k], marker=MARKER[k], ms=4,
                        lw=0.8, label=LABEL[k])
            ax.set_xscale("log")
            if row == 1:
                ax.set_yscale("log")
                ax.set_xlabel("Time step (s)")
            else:
                ax.axhline(0, **EXACT)
                ax.set_title(f"Level {level}", fontsize=8.5)
            if col == 0:
                ax.set_ylabel(ylabel)
            panel(ax, "")
    axes[0][0].legend(loc="lower right", fontsize=5.8)
    save(plt, fig, out, "time_step", dpi, png_only)


def cost(plt, data, out, dpi, png_only):
    rows = read_csv(data / "cases.csv")
    fig, ax = plt.subplots(figsize=(3.9, 3.0))
    for k in KEYS[:1] + KEYS[2:]:
        pts = sorted((num(r, "cells"), num(r, "rank_h")) for r in tutorial(rows)
                     if key(r) == k and r["rank_h"])
        ax.plot([p[0] for p in pts], [p[1] for p in pts], color=COLOR[k],
                ls=DASH[k], marker=MARKER[k], ms=4.2, lw=0.9,
                label=LABEL[k] if k != ("cutLink", "ME") else "cutLink")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Cells")
    ax.set_ylabel("Allocated rank-hours")
    ax.set_title("observed, single serial runs on a shared node", fontsize=7.5)
    ax.legend(loc="upper left", fontsize=6.5)
    panel(ax, "")
    save(plt, fig, out, "cost", dpi, png_only)


FIGURES = {"shear_histories": shear_histories, "shear_vs_mesh": shear_vs_mesh,
           "velocity_profiles": velocity_profiles,
           "velocity_convergence": velocity_convergence,
           "shear_noise": shear_noise, "transient": transient,
           "time_step": time_step, "cost": cost}


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--data", type=Path, default=REFERENCE_DIR / "results")
    parser.add_argument("--out", type=Path, default=Path("figures"))
    parser.add_argument("--dpi", type=int, default=300)
    parser.add_argument("--png-only", action="store_true")
    parser.add_argument("--only", nargs="+", default=list(FIGURES))
    args = parser.parse_args(argv)
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib is not installed: no figures written")
        return 0
    style(plt)
    for name in args.only:
        FIGURES[name](plt, args.data, args.out, args.dpi, args.png_only)
        print(f"{name}: written to {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
