#!/usr/bin/env python3
"""
Summary figures of the staticCylinderInChannel verification, from the
tracked data in verification/reference/ only: no CFD is run and work/ is not
needed. The data are validated before anything is drawn.

    legacy_incremental_sensitivity   Cd of the legacy incremental method on
                                     levels 2 and 4 against (a) the outer
                                     correctors, (b) the coupling coefficient
                                     and (c) the time step
    four_method_comparison           (a) Cd and (b) CPU hours of the four
                                     methods, with their frozen settings,
                                     against the cells per diameter

Each figure is written as PDF and PNG to
verification/postProcessing/figures/ (--out to change the directory). Text
uses Matplotlib MathText with STIX fonts.

Usage:
    python3 scripts/plot_verification.py
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

# Text and line sizes
SIZES = dict(font=8, tick=7.5, legend=7, small=6.8, line=1.0, marker=4.2)

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

def setup():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    v = SIZES
    plt.rcParams.update({
        "font.family": "serif",
        "font.serif": ["STIXGeneral", "DejaVu Serif"],
        "mathtext.fontset": "stix",
        "font.size": v["font"],
        "axes.labelsize": v["font"],
        "axes.linewidth": 0.6,
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
        "savefig.dpi": 300,
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
# Figures
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


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--out", default=str(VERIFICATION_DIR
                                             / "postProcessing" / "figures"),
                        help="output directory")
    args = parser.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    data = load_and_validate()
    plt, v = setup()
    written = []
    fig = figure_a(plt, v, data, "row")
    save(fig, out, "legacy_incremental_sensitivity", written)
    plt.close(fig)
    fig = figure_b(plt, v, data, "row", inset=True)
    save(fig, out, "four_method_comparison", written)
    plt.close(fig)
    for path in written:
        print(f"Wrote {path}.pdf and .png")


if __name__ == "__main__":
    main()
