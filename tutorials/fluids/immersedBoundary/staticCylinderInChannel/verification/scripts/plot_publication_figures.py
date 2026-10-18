#!/usr/bin/env python3
"""
Publication figures and tables of the staticCylinderInChannel verification,
from the compact data in verification/reference/ only: no CFD is run and
work/ is not needed. The data are validated before anything is drawn.

Layout variants (--variant, default manuscript; "all" writes all four):

    manuscript     Figure_A_legacy_incremental_sensitivity (175 mm, three
                       panels in a row) and Figure_A_..._single_column
                       (88 mm, stacked): (a) Cd vs n_outer, (b) Cd vs C_IB,
                       (c) Cd vs Delta t, legacy incremental, L2 and L4
                   Figure_B_method_comparison (175 mm) and
                       Figure_B_..._single_column (88 mm): (a) Cd and
                       (b) CPU-hours vs cells/D, four formulations;
                       Figure_B_method_comparison_inset: the same with a
                       zoom of the benchmark range in (a)
                   captions.md, historical_parity_table.{md,tex}
    expanded       fig_outer_corrector_sensitivity,
                   fig_coupling_coefficient_sensitivity,
                   fig_timestep_sensitivity, fig_four_method_Cd,
                   fig_four_method_CPU_hours (88 mm each), the tables
    presentation   talk_* figures with larger text
    review         contact_sheet.{png,pdf} of the manuscript and expanded
                   figures (requires them, written by --variant all)

Each figure is written as vector PDF and 600 dpi PNG to
verification/postProcessing/figures/<variant>/ (--out to change the root).
Text uses Matplotlib MathText with STIX fonts; --usetex renders it with an
external LaTeX installation instead.

Usage:
    python3 scripts/plot_publication_figures.py
    python3 scripts/plot_publication_figures.py --variant all
"""

import argparse
import csv
import sys
from pathlib import Path

VERIFICATION_DIR = Path(__file__).resolve().parent.parent
REFERENCE = VERIFICATION_DIR / "reference"

# Benchmark constants: the Schafer and Turek (1996) 2D-1 range and the
# high-precision reference value (FeatFlow; John and Matthies, 2001; Nabh,
# 1998)
CD_RANGE = (5.57, 5.59)
CD_REF = 5.57953523384

LEGACY = "legacy incremental"
BASELINE = "common-settings incremental baseline"
FORMULATIONS = ("cutLink", "ghostCell", "penalty", LEGACY)
# Frozen settings of the final comparison
FROZEN_OUTER = {"cutLink": 1, "ghostCell": 1, "penalty": 1, LEGACY: 5}
LEGACY_CC = 0.2
LEGACY_OUTER = 5
COMMON_DT = 0.001

# Okabe-Ito colours, with a marker and line style for each formulation so
# that the figures remain readable in greyscale
STYLE = {
    "cutLink": dict(color="#0072B2", marker="o", linestyle="-"),
    "ghostCell": dict(color="#D55E00", marker="s", linestyle="--"),
    "penalty": dict(color="#009E73", marker="^", linestyle="-."),
    LEGACY: dict(color="#CC79A7", marker="D", linestyle=":"),
}
LEVEL_STYLE = {
    2: dict(color="#000000", marker="o", linestyle="-",
            label=r"L2 (20 cells/$D$)"),
    4: dict(color="#7f7f7f", marker="s", linestyle="--",
            label=r"L4 (80 cells/$D$)"),
}
BAND = "#d9d9d9"
REF = "#555555"

MM = 1.0/25.4
FULL = 175*MM                  # full page width
COLUMN = 88*MM                 # single column width

VARIANT = {
    "manuscript": dict(font=8, tick=7.5, legend=7, small=6.8, line=1.0,
                       marker=4.2),
    "expanded": dict(font=8, tick=7.5, legend=7, small=6.8, line=1.0,
                     marker=4.2),
    "presentation": dict(font=14, tick=13, legend=12, small=12, line=1.8,
                         marker=7.5),
}

# Conditions of each sensitivity panel: what is held fixed
CONDITIONS = {
    "outer": r"$C_{\mathrm{IB}} = 0.2$, $\Delta t = 0.001$ s",
    "coupling": r"$n_{\mathrm{outer}} = 5$, $\Delta t = 0.001$ s",
    "timestep": r"$C_{\mathrm{IB}} = 0.2$, $n_{\mathrm{outer}} = 5$",
}


# -----------------------------------------------------------------------------
# Data
# -----------------------------------------------------------------------------

def read(name):
    path = REFERENCE / name
    with path.open() as handle:
        lines = [line for line in handle if not line.startswith("#")]
    return list(csv.DictReader(lines))


def num(row, key):
    return float(row[key])


def load_and_validate():
    """Read the final datasets and check that each one holds the intended
    cases; exit with a message if not"""
    outer_all = read("outer_corrector_study.csv")
    data = {
        "outer": [r for r in outer_all if r["method"] == "incremental"],
        "coupling": read("legacy_incremental_coupling.csv"),
        "timestep": read("legacy_incremental_timestep.csv"),
        "comparison": read("final_method_comparison.csv"),
        "parity": read("historical_parity.csv"),
    }
    errors = []

    def check(condition, message):
        if not condition:
            errors.append(message)

    def close(row, key, value):
        return abs(num(row, key) - value) < 1e-12

    outer = data["outer"]
    for level in (2, 4):
        have = sorted(int(r["nOuterCorrectors"]) for r in outer
                      if int(r["level"]) == level)
        check(have == [1, 3, 5, 7],
              f"outer correctors, L{level}: n_outer {have}, not [1, 3, 5, 7]")
    check(all(close(r, "couplingCoeff", LEGACY_CC)
              and close(r, "deltaT", COMMON_DT) for r in outer),
          "outer correctors: not all C_IB = 0.2, deltaT = 0.001")
    check(all(r["formulation"] == (BASELINE if int(r["nOuterCorrectors"]) == 1
                                   else LEGACY) for r in outer),
          "outer correctors: n_outer 1 must be the common-settings baseline")

    coupling = data["coupling"]
    check(all(r["formulation"] == LEGACY
              and int(r["nOuterCorrectors"]) == LEGACY_OUTER
              and close(r, "deltaT", COMMON_DT) for r in coupling),
          "coupling: not all legacy incremental, n_outer 5, deltaT 0.001")
    for level in (2, 4):
        have = sorted(num(r, "couplingCoeff") for r in coupling
                      if int(r["level"]) == level)
        check(have == [0.05, 0.2, 0.8],
              f"coupling, L{level}: C_IB {have}, not [0.05, 0.2, 0.8]")

    timestep = data["timestep"]
    check(all(r["formulation"] == LEGACY
              and int(r["nOuterCorrectors"]) == LEGACY_OUTER
              and close(r, "couplingCoeff", LEGACY_CC) for r in timestep),
          "timestep: not all legacy incremental, n_outer 5, C_IB 0.2")
    for level, want in ((2, [0.0005, 0.001, 0.004]),
                        (4, [0.00025, 0.0005, 0.001])):
        have = sorted(num(r, "deltaT") for r in timestep
                      if int(r["level"]) == level)
        check(have == want, f"timestep, L{level}: deltaT {have}, not {want}")

    comparison = data["comparison"]
    check(all(r["formulation"] in FORMULATIONS for r in comparison),
          "comparison: unexpected formulation (common-settings baseline?)")
    for f in FORMULATIONS:
        rows = [r for r in comparison if r["formulation"] == f]
        check(sorted(int(r["level"]) for r in rows) == [1, 2, 3, 4],
              f"comparison, {f}: levels not 1-4")
        check(all(int(r["nOuterCorrectors"]) == FROZEN_OUTER[f]
                  and close(r, "deltaT", COMMON_DT) for r in rows),
              f"comparison, {f}: not the frozen n_outer/deltaT")
        if f == LEGACY:
            check(all(close(r, "couplingCoeff", LEGACY_CC) for r in rows),
                  "comparison, legacy incremental: C_IB not 0.2")
    check(sorted(int(r["cells_per_D"]) for r in comparison)
          == sorted([10, 20, 40, 80]*4),
          "comparison: meshes other than the canonical levels")

    parity = data["parity"]
    check(len(parity) == 4 and {r["inlet"] for r in parity}
          == {"historical", "canonical"},
          "historical parity: expected 2 meshes x 2 inlets")
    check("cells_per_D" not in parity[0],
          "historical parity: has canonical mesh columns")

    if errors:
        sys.exit("Data validation failed:\n  " + "\n  ".join(errors))
    print("Data validated: outer correctors (n_outer 1 from the "
          "common-settings incremental baseline, 3-7 legacy incremental, "
          "C_IB 0.2, deltaT 0.001), coupling (legacy incremental, n_outer 5), "
          "time step (legacy incremental, C_IB 0.2, n_outer 5), comparison "
          "(frozen settings, L1-L4), historical parity (separate)")
    return data


# -----------------------------------------------------------------------------
# Drawing
# -----------------------------------------------------------------------------

def setup(variant, usetex):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    v = VARIANT.get(variant, VARIANT["manuscript"])
    plt.rcParams.update({
        "text.usetex": usetex,
        "font.family": "serif",
        "font.serif": ["STIXGeneral", "DejaVu Serif"],
        "mathtext.fontset": "stix",
        "font.size": v["font"],
        "axes.labelsize": v["font"],
        "axes.linewidth": 0.6 if variant != "presentation" else 1.0,
        "axes.grid": False,
        "xtick.direction": "in",
        "ytick.direction": "in",
        "xtick.top": True,
        "ytick.right": True,
        "xtick.major.width": 0.6,
        "ytick.major.width": 0.6,
        "xtick.minor.width": 0.4,
        "ytick.minor.width": 0.4,
        "xtick.major.size": 3.0,
        "ytick.major.size": 3.0,
        "xtick.minor.size": 1.6,
        "ytick.minor.size": 1.6,
        "xtick.labelsize": v["tick"],
        "ytick.labelsize": v["tick"],
        "lines.linewidth": v["line"],
        "lines.markersize": v["marker"],
        "legend.fontsize": v["legend"],
        "legend.frameon": False,
        "legend.handlelength": 2.2,
        "pdf.fonttype": 42,
        "savefig.dpi": 600,
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.03,
    })
    return plt, v


def benchmark(ax, labels=True):
    ax.axhspan(*CD_RANGE, color=BAND, alpha=0.6, linewidth=0, zorder=0,
               label="Schäfer–Turek benchmark range" if labels else None)
    ax.axhline(CD_REF, color=REF, linewidth=0.6, linestyle=(0, (4, 2)),
               zorder=1, label="reference value" if labels else None)


def plot_series(ax, xs, ys, statuses, style, label):
    ax.plot(xs, ys, color=style["color"], linestyle=style["linestyle"],
            zorder=3, label=label, marker=style["marker"],
            markerfacecolor=style["color"], markeredgecolor=style["color"],
            markeredgewidth=0.8)
    # Limit cycles: hollow markers (the value is the mean over the final
    # averaging interval)
    for x, y, s in zip(xs, ys, statuses):
        if s != "steady":
            ax.plot(x, y, linestyle="none", marker=style["marker"],
                    markerfacecolor="white", markeredgecolor=style["color"],
                    markeredgewidth=0.8, zorder=4)


def panel_label(ax, text, condition=None, small=None):
    ax.text(0.03, 0.95, text, transform=ax.transAxes, va="top", ha="left")
    if condition:
        ax.text(0.03, 0.85, condition, transform=ax.transAxes, va="top",
                ha="left", fontsize=small, color="#333333")


def level_series(ax, rows, xkey):
    for level in (2, 4):
        pts = sorted((num(r, xkey), num(r, "Cd_mean"), r["status"])
                     for r in rows if int(r["level"]) == level)
        plot_series(ax, [p[0] for p in pts], [p[1] for p in pts],
                    [p[2] for p in pts], LEVEL_STYLE[level],
                    LEVEL_STYLE[level]["label"])


def limit_cycle_proxy(ax, label="limit cycle (hollow marker)"):
    ax.plot([], [], linestyle="none", marker="o", markerfacecolor="white",
            markeredgecolor="#000000", markeredgewidth=0.8, label=label)


def log_axis(ax, values):
    import matplotlib.ticker as ticker
    ax.set_xscale("log")
    ax.xaxis.set_major_locator(ticker.FixedLocator(values))
    ax.xaxis.set_major_formatter(
        ticker.FixedFormatter([f"{v:g}" for v in values]))
    ax.xaxis.set_minor_locator(ticker.NullLocator())
    ax.set_xlim(min(values)/1.5, max(values)*1.5)


def draw_outer(ax, data):
    benchmark(ax)
    level_series(ax, data["outer"], "nOuterCorrectors")
    ax.set_xticks([1, 3, 5, 7])
    ax.set_xlim(0.4, 7.6)
    ax.set_xlabel(r"$n_{\mathrm{outer}}$")
    ax.set_ylabel(r"$C_d$")


def draw_coupling(ax, data):
    benchmark(ax)
    level_series(ax, data["coupling"], "couplingCoeff")
    log_axis(ax, sorted({num(r, "couplingCoeff") for r in data["coupling"]}))
    ax.set_xlabel(r"$C_{\mathrm{IB}}$")
    ax.set_ylabel(r"$C_d$")


def draw_timestep(ax, data):
    benchmark(ax)
    level_series(ax, data["timestep"], "deltaT")
    log_axis(ax, sorted({num(r, "deltaT") for r in data["timestep"]}))
    ax.set_xlabel(r"$\Delta t$ (s)")
    ax.set_ylabel(r"$C_d$")


def method_series(ax, rows, ykey):
    for f in FORMULATIONS:
        pts = sorted((num(r, "cells_per_D"), num(r, ykey), r["status"])
                     for r in rows if r["formulation"] == f)
        plot_series(ax, [p[0] for p in pts], [p[1] for p in pts],
                    [p[2] for p in pts], STYLE[f], f)


def cells_axis(ax, xlim=(8, 100)):
    import matplotlib.ticker as ticker
    ax.set_xscale("log", base=2)
    ax.xaxis.set_major_locator(ticker.FixedLocator([10, 20, 40, 80]))
    ax.xaxis.set_major_formatter(
        ticker.FixedFormatter(["10", "20", "40", "80"]))
    ax.xaxis.set_minor_locator(ticker.NullLocator())
    ax.set_xlim(*xlim)
    ax.set_xlabel(r"cells/$D$")


def draw_cd(ax, data, headroom=0.12):
    import matplotlib.ticker as ticker
    benchmark(ax)
    method_series(ax, data["comparison"], "Cd_mean")
    cells_axis(ax)
    lo, hi = ax.get_ylim()
    ax.set_ylim(lo, hi + headroom*(hi - lo))
    ax.yaxis.set_minor_locator(ticker.AutoMinorLocator())
    ax.set_ylabel(r"$C_d$")


def draw_cpu(ax, data):
    method_series(ax, data["comparison"], "cpu_hours")
    cells_axis(ax)
    ax.set_yscale("log")
    ax.set_ylabel("CPU-hours")


def zoom_inset(ax, data, bounds):
    """Inset of the Cd panel around the benchmark range, at the full x range
    of the panel; the main axes are not changed"""
    import matplotlib.ticker as ticker
    ins = ax.inset_axes(bounds)
    benchmark(ins, labels=False)
    method_series(ins, data["comparison"], "Cd_mean")
    cells_axis(ins)
    ins.set_xlabel("")
    ins.set_ylim(5.555, 5.607)
    ins.yaxis.set_major_locator(ticker.MultipleLocator(0.02))
    ins.yaxis.set_minor_locator(ticker.MultipleLocator(0.005))
    ins.tick_params(labelsize=6, length=2.0, pad=1.5)
    ins.set_facecolor("white")
    for line in ins.get_lines():
        line.set_markersize(3.0)
    return ins


def legacy_ylim(axes, data, keys=("outer", "coupling", "timestep"),
                top=0.42):
    """The same C_d range on the given sensitivity panels (by default that of
    all three studies, so that the panels of Figure A share a scale), with
    headroom for the panel labels and conditions"""
    import matplotlib.ticker as ticker
    values = [num(r, "Cd_mean") for key in keys for r in data[key]]
    lo = min(values + [CD_RANGE[0]])
    hi = max(values + [CD_RANGE[1]])
    span = hi - lo
    for ax in axes:
        ax.set_ylim(lo - 0.06*span, hi + top*span)
        ax.yaxis.set_minor_locator(ticker.AutoMinorLocator())


def top_legend(fig, ax, ncol):
    handles, labels = ax.get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=ncol,
               bbox_to_anchor=(0.5, 1.0), columnspacing=1.4)


def save(fig, out, name, written):
    for ext in ("pdf", "png"):
        fig.savefig(out / f"{name}.{ext}")
    written.append(out / name)


# -----------------------------------------------------------------------------
# Variants
# -----------------------------------------------------------------------------

def figure_a(plt, v, data, layout):
    if layout == "row":
        fig, axes = plt.subplots(1, 3, figsize=(FULL, 2.35), sharey=True)
    else:
        fig, axes = plt.subplots(3, 1, figsize=(COLUMN, 6.2), sharey=True)
    draw_outer(axes[0], data)
    draw_coupling(axes[1], data)
    draw_timestep(axes[2], data)
    legacy_ylim(axes, data)
    for ax, label, key in zip(axes, ("(a)", "(b)", "(c)"),
                              ("outer", "coupling", "timestep")):
        panel_label(ax, label, CONDITIONS[key], v["small"])
    if layout == "row":
        for ax in axes[1:]:
            ax.set_ylabel("")
    limit_cycle_proxy(axes[0])
    fig.tight_layout(w_pad=0.8, h_pad=0.6)
    top_legend(fig, axes[0], ncol=5 if layout == "row" else 2)
    return fig


def figure_b(plt, v, data, layout, inset=False):
    if layout == "row":
        fig, axes = plt.subplots(1, 2, figsize=(FULL, 2.7))
    else:
        fig, axes = plt.subplots(2, 1, figsize=(COLUMN, 5.0))
    draw_cd(axes[0], data)
    draw_cpu(axes[1], data)
    for ax, label in zip(axes, ("(a)", "(b)")):
        panel_label(ax, label)
    if inset:
        # The y-axis of (a) is extended downwards (the scale and the data
        # are unchanged), and the inset is placed below the lowest point, so
        # that it hides no data
        lo, hi = axes[0].get_ylim()
        new_lo = lo - 0.55*(hi - lo)
        axes[0].set_ylim(new_lo, hi)
        lowest = min(num(r, "Cd_mean") for r in data["comparison"])
        top = (lowest - new_lo)/(hi - new_lo) - 0.05
        zoom_inset(axes[0], data, (0.30, 0.05, 0.67, top - 0.05))
    limit_cycle_proxy(axes[0])
    fig.tight_layout(w_pad=1.5, h_pad=0.6)
    top_legend(fig, axes[0], ncol=4 if layout == "row" else 2)
    return fig


def manuscript(plt, v, data, out, written):
    for layout, suffix in (("row", ""), ("column", "_single_column")):
        fig = figure_a(plt, v, data, layout)
        save(fig, out, f"Figure_A_legacy_incremental_sensitivity{suffix}",
             written)
        plt.close(fig)
        fig = figure_b(plt, v, data, layout)
        save(fig, out, f"Figure_B_method_comparison{suffix}", written)
        plt.close(fig)
    fig = figure_b(plt, v, data, "row", inset=True)
    save(fig, out, "Figure_B_method_comparison_inset", written)
    plt.close(fig)


def expanded(plt, v, data, out, written):
    size = (COLUMN, 2.55)
    for name, draw, key in (
        ("fig_outer_corrector_sensitivity", draw_outer, "outer"),
        ("fig_coupling_coefficient_sensitivity", draw_coupling, "coupling"),
        ("fig_timestep_sensitivity", draw_timestep, "timestep"),
    ):
        fig, ax = plt.subplots(figsize=size)
        draw(ax, data)
        legacy_ylim([ax], data, (key,), top=0.22)
        ax.text(0.03, 0.95, CONDITIONS[key], transform=ax.transAxes,
                va="top", ha="left", fontsize=v["small"], color="#333333")
        limit_cycle_proxy(ax)
        fig.tight_layout()
        top_legend(fig, ax, ncol=2)
        save(fig, out, name, written)
        plt.close(fig)
    for name, draw in (("fig_four_method_Cd", draw_cd),
                       ("fig_four_method_CPU_hours", draw_cpu)):
        fig, ax = plt.subplots(figsize=size)
        draw(ax, data)
        limit_cycle_proxy(ax)
        fig.tight_layout()
        top_legend(fig, ax, ncol=2 if draw is draw_cpu else 3)
        save(fig, out, name, written)
        plt.close(fig)


def presentation(plt, v, data, out, written):
    def simple_legend(fig, ax, ncol):
        # Only the series, the benchmark range and the limit cycles
        handles, labels = ax.get_legend_handles_labels()
        keep = [(h, lab) for h, lab in zip(handles, labels)
                if lab != "reference value"]
        fig.legend([h for h, _ in keep], [lab for _, lab in keep],
                   loc="lower center", ncol=ncol, bbox_to_anchor=(0.5, 1.0))

    fig, ax = plt.subplots(figsize=(7.2, 4.6))
    draw_outer(ax, data)
    legacy_ylim([ax], data, ("outer",), top=0.2)
    ax.text(0.03, 0.95, "legacy incremental, " + CONDITIONS["outer"],
            transform=ax.transAxes, va="top", ha="left", fontsize=v["small"])
    limit_cycle_proxy(ax, "limit cycle")
    fig.tight_layout()
    simple_legend(fig, ax, ncol=2)
    save(fig, out, "talk_outer_corrector_sensitivity", written)
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(12.0, 4.6), sharey=True)
    draw_coupling(axes[0], data)
    draw_timestep(axes[1], data)
    legacy_ylim(axes, data, ("coupling", "timestep"), top=0.3)
    axes[1].set_ylabel("")
    for ax, label, key in zip(axes, ("(a)", "(b)"), ("coupling", "timestep")):
        panel_label(ax, label, CONDITIONS[key], v["small"])
    limit_cycle_proxy(axes[0], "limit cycle")
    fig.tight_layout(w_pad=1.0)
    simple_legend(fig, axes[0], ncol=4)
    save(fig, out, "talk_coupling_and_timestep", written)
    plt.close(fig)

    for name, draw in (("talk_four_method_Cd", draw_cd),
                       ("talk_four_method_CPU_hours", draw_cpu)):
        fig, ax = plt.subplots(figsize=(7.2, 4.6))
        draw(ax, data)
        limit_cycle_proxy(ax, "limit cycle")
        fig.tight_layout()
        simple_legend(fig, ax, ncol=3)
        save(fig, out, name, written)
        plt.close(fig)


def review(plt, root, written):
    """Contact sheet of the manuscript and expanded figures"""
    import matplotlib.image as mpimg
    panels = [
        ("manuscript", "Figure_A_legacy_incremental_sensitivity"),
        ("manuscript", "Figure_B_method_comparison"),
        ("manuscript", "Figure_B_method_comparison_inset"),
        ("manuscript", "Figure_A_legacy_incremental_sensitivity_single_column"),
        ("manuscript", "Figure_B_method_comparison_single_column"),
        ("expanded", "fig_outer_corrector_sensitivity"),
        ("expanded", "fig_coupling_coefficient_sensitivity"),
        ("expanded", "fig_timestep_sensitivity"),
        ("expanded", "fig_four_method_Cd"),
        ("expanded", "fig_four_method_CPU_hours"),
    ]
    missing = [f"{d}/{n}.png" for d, n in panels
               if not (root / d / f"{n}.png").exists()]
    if missing:
        sys.exit("review: run --variant all first; missing " + ", ".join(
            missing))
    fig = plt.figure(figsize=(16, 22))
    grid = fig.add_gridspec(6, 6, height_ratios=(1.0, 1.1, 1.1, 2.2, 1.0,
                                                 1.0))
    slots = [grid[0, :], grid[1, :], grid[2, :], grid[3, 0:3], grid[3, 3:6],
             grid[4, 0:2], grid[4, 2:4], grid[4, 4:6], grid[5, 0:3],
             grid[5, 3:6]]
    for (d, n), slot in zip(panels, slots):
        ax = fig.add_subplot(slot)
        ax.imshow(mpimg.imread(root / d / f"{n}.png"))
        ax.set_axis_off()
        ax.set_title(f"{d}/{n}", fontsize=9, family="monospace")
    fig.tight_layout()
    out = root / "review"
    out.mkdir(parents=True, exist_ok=True)
    fig.savefig(out / "contact_sheet.png", dpi=150)
    fig.savefig(out / "contact_sheet.pdf")
    plt.close(fig)
    written.append(out / "contact_sheet")


# -----------------------------------------------------------------------------
# Tables and captions
# -----------------------------------------------------------------------------

def captions(out, written):
    hollow = ("Hollow markers denote bounded limit-cycle results, for which "
              "the plotted value is the mean over the final averaging "
              "interval (2 s); filled markers are steady results.")
    band = (f"The grey band is the Schäfer–Turek benchmark range "
            f"{CD_RANGE[0]} ≤ C_d ≤ {CD_RANGE[1]}, and the dashed line the "
            f"high-precision reference value C_d = {CD_REF} (FeatFlow; "
            f"John and Matthies, 2001).")
    text = [
        "# Figure captions (draft)", "",
        "**Figure A.** Sensitivity of the drag coefficient of the legacy "
        "incremental formulation (vertexFraction occupancy, n_corr = 3) on "
        "L2 (20 cells/D) and L4 (80 cells/D) to (a) the number of outer "
        "correctors n_outer, with C_IB = 0.2 and Δt = 0.001 s; (b) the "
        "coupling coefficient C_IB, with n_outer = 5 and Δt = 0.001 s; and "
        "(c) the time step Δt, with C_IB = 0.2 and n_outer = 5. In (a), the "
        "values at n_outer = 1 are from the common-settings incremental "
        "baseline with the same C_IB and Δt. " + band + " " + hollow, "",
        "**Figure B.** (a) Drag coefficient and (b) CPU-hours (allocated "
        "ranks × wall time) against the mesh resolution for cutLink, "
        "ghostCell and penalty (n_outer = 1) and legacy incremental "
        "(C_IB = 0.2, vertexFraction, n_outer = 5), all with n_corr = 3 and "
        "Δt = 0.001 s; 1, 1, 2 and 8 MPI ranks on 10, 20, 40 and 80 cells/D. "
        "The CPU-hours of legacy incremental include its five outer "
        "correctors. " + band + " " + hollow + " In the inset variant, the "
        "inset of (a) magnifies the benchmark range at the same mesh "
        "resolutions; the C_d axis of (a) is extended downwards to make room "
        "for it below the data, with its scale unchanged.", "",
        "**Table (historical parity).** Drag coefficient of the "
        "manuscript-era pimpleHFDIBFoam static-cylinder cases (D = 1, "
        "ρ = 1000, Ū = 1) as saved, and as computed by the current "
        "solids4foam legacy incremental port with the historical numerical "
        "settings (C_IB = 0.2, vertexFraction, Δt = 0.005 s, n_outer = 5, "
        "n_corr = 3, two non-orthogonal correctors, relaxation 0.7/0.3), "
        "with the historical inlet and with the exact parabolic inlet.", "",
    ]
    (out / "captions.md").write_text("\n".join(text))
    written.append(out / "captions")


def parity_table(data, out, written):
    parity = data["parity"]
    meshes = sorted({r["historical_case"] for r in parity},
                    key=lambda m: int(m.split("Mesh")[-1]))
    md = [
        "| Mesh | Saved manuscript $C_d$ | Current port, historical inlet "
        "| Difference | Current port, canonical inlet | Inlet-induced change |",
        "| --- | --- | --- | --- | --- | --- |",
    ]
    tex = [
        r"\begin{tabular}{lccccc}", r"\hline",
        r"Mesh & Saved $C_d$ & Port, historical inlet & Difference & "
        r"Port, canonical inlet & Inlet-induced change \\", r"\hline",
    ]
    for mesh in meshes:
        hist = next(r for r in parity if r["historical_case"] == mesh
                    and r["inlet"] == "historical")
        canon = next(r for r in parity if r["historical_case"] == mesh
                     and r["inlet"] == "canonical")
        saved = num(hist, "Cd_historical")
        h = num(hist, "Cd_current")
        c = num(canon, "Cd_current")
        cells_tex = f"{int(hist['n_cells']):,}".replace(",", r"\,")
        cells_md = f"{int(hist['n_cells']):,}".replace(",", " ")
        md.append(f"| {cells_md} cells | {saved:.5f} | {h:.6f} | "
                  f"{h - saved:+.1e} | {c:.6f} | {100*(c - h)/h:+.2f}% |")
        tex.append(f"{cells_tex} cells & {saved:.5f} & {h:.6f} & "
                   f"${h - saved:+.1e}$ & {c:.6f} & "
                   f"${100*(c - h)/h:+.2f}$\\% " r"\\")
    tex += [r"\hline", r"\end{tabular}"]
    (out / "historical_parity_table.md").write_text("\n".join(md) + "\n")
    (out / "historical_parity_table.tex").write_text("\n".join(tex) + "\n")
    written.append(out / "historical_parity_table")


def comparison_table(data, out, written):
    md = ["| Formulation | $n_{outer}$ | Level | cells/$D$ | $C_d$ | $C_l$ | "
          "Status | Ranks | Wall time (s) | CPU-hours |",
          "| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |"]
    tex = [r"\begin{tabular}{llrrcclrrr}", r"\hline",
           r"Formulation & $n_{\mathrm{outer}}$ & Level & cells/$D$ & $C_d$ & "
           r"$C_l$ & Status & Ranks & Wall time (s) & CPU-hours \\",
           r"\hline"]
    for f in FORMULATIONS:
        for r in sorted((r for r in data["comparison"]
                         if r["formulation"] == f),
                        key=lambda r: int(r["level"])):
            cells = (f"{f} & {r['nOuterCorrectors']} & L{r['level']} & "
                     f"{r['cells_per_D']} & {num(r, 'Cd_mean'):.5f} & "
                     f"{num(r, 'Cl_mean'):.5f} & {r['status']} & "
                     f"{r['ranks']} & {num(r, 'wall_s'):.0f} & "
                     f"{num(r, 'cpu_hours'):.3f}")
            md.append("| " + cells.replace(" & ", " | ") + " |")
            tex.append(cells + r" \\")
    tex += [r"\hline", r"\end{tabular}"]
    (out / "final_comparison_table.md").write_text("\n".join(md) + "\n")
    (out / "final_comparison_table.tex").write_text("\n".join(tex) + "\n")
    written.append(out / "final_comparison_table")


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--variant", default="manuscript",
                        choices=("manuscript", "expanded", "presentation",
                                 "review", "all"))
    parser.add_argument("--out", default=str(VERIFICATION_DIR
                                             / "postProcessing" / "figures"))
    parser.add_argument("--usetex", action="store_true",
                        help="render text with an external LaTeX "
                             "installation (default: MathText)")
    args = parser.parse_args()
    root = Path(args.out)

    data = load_and_validate()
    variants = (("manuscript", "expanded", "presentation", "review")
                if args.variant == "all" else (args.variant,))
    written = []
    for variant in variants:
        plt, v = setup(variant, args.usetex)
        if variant == "review":
            review(plt, root, written)
            continue
        out = root / variant
        out.mkdir(parents=True, exist_ok=True)
        {"manuscript": manuscript, "expanded": expanded,
         "presentation": presentation}[variant](plt, v, data, out, written)
        if variant in ("manuscript", "expanded"):
            parity_table(data, out, written)
            comparison_table(data, out, written)
        if variant == "manuscript":
            captions(out, written)
    for path in written:
        print(f"Wrote {path}")


if __name__ == "__main__":
    main()
