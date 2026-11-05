#!/usr/bin/env python3
"""Figures of the immersedTaylorCouette verification from the tracked results
in reference/results/ (or those of './Allverify analyse' with --data). No CFD
is run and nothing is filtered or shifted.

    torque_histories     rotor torque against time, Levels 2 and 3, with the
                         exact torque and the window used for the mean
    torque_vs_mesh       torque difference from the exact torque against
                         the mesh level (tutorial time step), and with
                         the time step capped at 0.01 s
    torque_vs_dt         torque difference against the time-step cap on
                         Levels 1 and 2 (separates the time step from the mesh)
    velocity_profiles    azimuthal velocity against the exact velocity in the
                         gap (azimuthal averages of the cell centres) and the
                         difference near the rotor, Level 3
    velocity_convergence L1, L2 and Linf norms of the azimuthal velocity
                         difference against the cell width
    torque_variance      standard deviation of the torque in the window
    cost                 observed cost against the number of cells

Style of the immersed-boundary verifications: one colour, marker and line
pattern per method (Okabe-Ito colours), references in black. Only matplotlib
is needed (no LaTeX).

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

# (method, estimator) of a variant -> style
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
R1, R2 = 0.05, 0.1


def read_csv(path: Path) -> list:
    with path.open() as handle:
        return list(csv.DictReader(handle))


def num(row: dict, key: str) -> float:
    return float(row[key]) if row[key] not in ("", "nan") else float("nan")


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


def series(rows, k, cap="", ykey="T_relerr", levels=(1, 2, 3, 4)):
    sel = sorted((int(r["level"]), num(r, ykey)) for r in rows
                 if key(r) == k and r["dt_cap"] == cap
                 and int(r["level"]) in levels)
    return [a for a, _ in sel], [b for _, b in sel]


def h_of(level):
    return 0.32 / (32 * 2 ** (level - 1))


def save(plt, fig, out, name, dpi, png_only):
    out.mkdir(parents=True, exist_ok=True)
    fig.savefig(out / f"{name}.png", dpi=dpi)
    if not png_only:
        fig.savefig(out / f"{name}.pdf")
    plt.close(fig)


def torque_histories(plt, data, out, dpi, png_only):
    rows = read_csv(data / "cases.csv")
    exact = float(rows[0]["T_exact"])
    fig, axes = plt.subplots(2, 2, figsize=(7.2, 5.4))
    for col, (level, lab) in enumerate(zip((2, 3), "ab")):
        for k in KEYS:
            vid = k[0] + ("-" + k[1] if k[0] == "cutLink" else "")
            f = data / f"history-{vid}-L{level}.csv"
            if not f.exists():
                continue
            h = read_csv(f)
            t = [num(r, "time") for r in h]
            tz = [num(r, "rotor_torque_z") for r in h]
            axes[0][col].plot(t, [1e6 * x for x in tz], color=COLOR[k],
                              ls=DASH[k], lw=0.9, label=LABEL[k])
            axes[1][col].plot(t, [100 * (x / exact - 1) for x in tz],
                              color=COLOR[k], ls=DASH[k], lw=0.9)
        ax = axes[0][col]
        ax.axhline(1e6 * exact, **EXACT, label="exact")
        ax.axvspan(8, 10, color="0.9", zorder=0)
        ax.set_title(f"Level {level}, tutorial time step", fontsize=8.5)
        ax.set_xlim(0, 10)
        ax.set_ylim(-4.6, 0.2)
        ax.set_xlabel("Time (s)")
        panel(ax, lab)
        ax = axes[1][col]
        ax.axhline(0, **EXACT)
        ax.axvspan(8, 10, color="0.9", zorder=0)
        ax.set_xlim(2, 10)
        ax.set_ylim(-18, 4) if level == 2 else ax.set_ylim(-12, 4)
        ax.set_xlabel("Time (s)")
        panel(ax, "")
    axes[0][0].set_ylabel(r"Rotor torque, $T_z$ ($\mu$N m)")
    axes[1][0].set_ylabel(r"$(T - T_{\mathrm{exact}})/T_{\mathrm{exact}}$ (%)")
    handles, labels = axes[0][1].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=3,
               bbox_to_anchor=(0.5, -0.06))
    save(plt, fig, out, "torque_histories", dpi, png_only)


def torque_vs_mesh(plt, data, out, dpi, png_only):
    rows = read_csv(data / "cases.csv")
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.3), sharey=True)
    for ax, cap, title, lab in ((axes[0], "", "tutorial time step "
                                 "(Courant 0.5)", "a"),
                                (axes[1], "0.005", "time step $\\leq$ 0.005 s",
                                 "b")):
        keys = KEYS + (("cutLink", "FO"),) if cap == "" else KEYS
        for k in keys:
            lv, e = series(rows, k, cap)
            if not lv:
                continue
            ax.plot([10 * 2 ** (a - 1) for a in lv], [100 * abs(x) for x in e],
                    color=COLOR[k], ls=DASH[k], marker=MARKER[k], ms=4.5,
                    lw=0.9, label=LABEL[k], mfc=COLOR[k])
        ax.set_xscale("log", base=2)
        ax.set_yscale("log")
        ax.set_xlabel("Cells across the rotor")
        ax.set_xticks([10, 20, 40, 80])
        ax.set_xticklabels(["10", "20", "40", "80"])
        ax.set_ylim(0.03, 40)
        ax.set_title(title, fontsize=8.5)
        panel(ax, lab)
    axes[0].set_ylabel(r"$|T - T_{\mathrm{exact}}|/|T_{\mathrm{exact}}|$ (%)")
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=3,
               bbox_to_anchor=(0.5, -0.17))
    save(plt, fig, out, "torque_vs_mesh", dpi, png_only)


def torque_vs_dt(plt, data, out, dpi, png_only):
    rows = read_csv(data / "cases.csv")
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.3), sharey=True)
    methods = (("cutLink", "ME"), ("penalty", ""), ("incremental", ""),
               ("ghostCell", ""))
    for ax, level, lab in zip(axes, (1, 2), "ab"):
        for k in methods:
            pts = []
            for r in rows:
                if key(r) == k and int(r["level"]) == level:
                    cap = num(r, "dt_cap")
                    if cap != cap:                    # tutorial step
                        pts.append((10.0 / num(r, "steps"),
                                    num(r, "T_relerr"), True))
                    else:
                        pts.append((cap, num(r, "T_relerr"), False))
            pts.sort()
            ax.plot([p[0] for p in pts], [100 * p[1] for p in pts],
                    color=COLOR[k], ls=DASH[k], marker=MARKER[k], ms=4.5,
                    lw=0.9, label=LABEL[k])
            ax.plot([p[0] for p in pts if p[2]], [100 * p[1] for p in pts
                                                  if p[2]],
                    ls="none", marker=MARKER[k], ms=8, mfc="none",
                    color=COLOR[k], mew=0.8)
        ax.axhline(0, **EXACT)
        ax.set_xscale("log")
        ax.set_xlabel("Time step cap, or mean step (circled) (s)")
        ax.set_title(f"Level {level}", fontsize=8.5)
        panel(ax, lab)
    axes[0].set_ylabel(r"$(T - T_{\mathrm{exact}})/T_{\mathrm{exact}}$ (%)")
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=4,
               bbox_to_anchor=(0.5, -0.1))
    save(plt, fig, out, "torque_vs_dt", dpi, png_only)


def velocity_profiles(plt, data, out, dpi, png_only):
    level = 3
    prof = read_csv(data / f"profiles-L{level}.csv")
    h = h_of(level)
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.1))
    r_ex = [R1 + (R2 - R1) * i / 400 for i in range(401)]
    a_ = -1.0 * R1 ** 2 / (R2 ** 2 - R1 ** 2)
    b_ = R1 ** 2 * R2 ** 2 / (R2 ** 2 - R1 ** 2)
    ax = axes[0]
    ax.plot(r_ex, [a_ * r + b_ / r for r in r_ex], **EXACT, label="exact")
    for k, vid in ((("cutLink", "ME"), "cutLink"), (("penalty", ""), "penalty"),
                   (("incremental", ""), "incremental"),
                   (("ghostCell", ""), "ghostCell")):
        rr = [r for r in prof if r["variant"] == vid]
        ax.plot([num(r, "r") for r in rr], [num(r, "u_theta") for r in rr],
                color=COLOR[k], ls="none", marker=MARKER[k], ms=3.2,
                mfc="none", label=LABEL[k] if k != ("cutLink", "ME")
                else "cutLink")
    ax.set_xlabel("Radius (m)")
    ax.set_ylabel(r"$u_\theta$ (m/s)")
    ax.set_title("Level 3, azimuthal averages", fontsize=8.5)
    ax.legend(loc="upper right")
    panel(ax, "a")
    ax = axes[1]
    for k, vid in ((("cutLink", "ME"), "cutLink"), (("penalty", ""), "penalty"),
                   (("incremental", ""), "incremental"),
                   (("ghostCell", ""), "ghostCell")):
        rr = [r for r in prof if r["variant"] == vid and
              num(r, "r") < R1 + 6 * h]
        ax.plot([(num(r, "r") - R1) / h for r in rr],
                [1e3 * (num(r, "u_theta") - num(r, "u_theta_exact"))
                 for r in rr],
                color=COLOR[k], ls=DASH[k], marker=MARKER[k], ms=3.8, lw=0.9)
    ax.axhline(0, **EXACT)
    ax.set_xlabel("Distance from the rotor surface (cell widths)")
    ax.set_ylabel(r"$u_\theta - u_{\theta,\mathrm{exact}}$ (mm/s)")
    ax.set_title("Level 3, near the rotor", fontsize=8.5)
    panel(ax, "b")
    save(plt, fig, out, "velocity_profiles", dpi, png_only)


def velocity_convergence(plt, data, out, dpi, png_only):
    rows = read_csv(data / "cases.csv")
    fig, axes = plt.subplots(1, 3, figsize=(7.4, 2.9), sharey=True)
    for ax, norm, lab in zip(axes, ("ut_L1", "ut_L2", "ut_Linf"), "abc"):
        for k in KEYS[:1] + KEYS[2:]:
            lv, e = series(rows, k, "", norm)
            ax.plot([h_of(a) * 1e3 for a in lv], e, color=COLOR[k],
                    ls=DASH[k], marker=MARKER[k], ms=4.2, lw=0.9,
                    label=LABEL[k] if k != ("cutLink", "ME") else "cutLink")
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("Cell width (mm)")
        ax.set_title({"ut_L1": "L1", "ut_L2": "L2", "ut_Linf": "Linf"}[norm],
                     fontsize=8.5)
        panel(ax, lab)
    axes[0].set_ylabel(r"norm of $u_\theta - u_{\theta,\mathrm{exact}}$ "
                       r"over $\omega R_1$")
    axes[0].legend(loc="lower right")
    save(plt, fig, out, "velocity_convergence", dpi, png_only)


def torque_variance(plt, data, out, dpi, png_only):
    rows = read_csv(data / "cases.csv")
    fig, ax = plt.subplots(figsize=(3.8, 3.0))
    for k in KEYS:
        lv, e = series(rows, k, "", "T_std_rel")
        ax.plot([10 * 2 ** (a - 1) for a in lv], e, color=COLOR[k],
                ls=DASH[k], marker=MARKER[k], ms=4.2, lw=0.9, label=LABEL[k])
    ax.set_xscale("log", base=2)
    ax.set_yscale("log")
    ax.set_xticks([10, 20, 40, 80])
    ax.set_xticklabels(["10", "20", "40", "80"])
    ax.set_xlabel("Cells across the rotor")
    ax.set_ylabel(r"std. of $T_z$ in $8 \leq t \leq 10$ s over $|T_{exact}|$")
    ax.legend(loc="lower left", fontsize=6.5)
    panel(ax, "")
    save(plt, fig, out, "torque_variance", dpi, png_only)


def cost(plt, data, out, dpi, png_only):
    rows = read_csv(data / "cases.csv")
    fig, ax = plt.subplots(figsize=(3.9, 3.0))
    for k in KEYS:
        pts = sorted((num(r, "cells"), num(r, "rank_h")) for r in rows
                     if key(r) == k and r["dt_cap"] == "" and r["rank_h"])
        ax.plot([p[0] for p in pts], [p[1] for p in pts], color=COLOR[k],
                ls=DASH[k], marker=MARKER[k], ms=4.2, lw=0.9, label=LABEL[k])
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Cells")
    ax.set_ylabel("Allocated rank-hours")
    ax.set_title("observed, single runs on a shared node", fontsize=7.5)
    ax.legend(loc="upper left", fontsize=6.5)
    panel(ax, "")
    save(plt, fig, out, "cost", dpi, png_only)


FIGURES = {"torque_histories": torque_histories,
           "torque_vs_mesh": torque_vs_mesh, "torque_vs_dt": torque_vs_dt,
           "velocity_profiles": velocity_profiles,
           "velocity_convergence": velocity_convergence,
           "torque_variance": torque_variance, "cost": cost}


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
