#!/usr/bin/env python3
"""Figures of the controls of the oscillatingWallStokesLayer verification
(time-step convergence, matched-step mesh sequence, shear reconstruction,
ghostCell interface error, momentum-exchange wall offset, cell-average
comparison) from the tracked results in reference/results/ and
reference/results/controls/. Nothing is smoothed, filtered or shifted.

    controls_dt_velocity   velocity L2 error against the time step, every level
    controls_dt_shear      shear amplitude, phase and history error against the
                           time step on Level 4
    controls_matched_dt    velocity L2 error against the cell width at two
                           fixed steps (tutorial 0.001 s and 0.000125 s)
    controls_alignment     raw and inertia-aligned (post-processing) shear
                           amplitude error against the mesh at dt 0.001 and 0.000125 s
    controls_ghost_regions ghostCell velocity error in all fluid rows, within
                           one boundary-layer thickness and in the interface row
    controls_offset        (ME - FO) shear against the wall offset, and the
                           amplitude error with the wall on a mesh face
    controls_cellavg       velocity error against the cell-centre value and the
                           exact cell average

Usage: python3 scripts/plot_controls.py [--out DIR] [--dpi N] [--png-only]
"""

import argparse
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import plot_verification as pv  # noqa: E402
import stokes_controls as sc  # noqa: E402

METHODS = (("cutLink", "ME"), ("penalty", ""), ("incremental", ""), ("ghostCell", ""))
DT_REF = 0.000125


def merged():
    return sc.merged_cases()


def at(rows, k, level, dt):
    for r in rows:
        if pv.key(r) == k and int(r["level"]) == level and abs(float(r["dt"]) - dt) < 1e-12:
            return r
    return None


def tidy_dt_axis(ax):
    from matplotlib.ticker import NullFormatter
    ax.set_xticks([1e-3, 1e-4, 1e-5])
    ax.set_xticklabels(["$10^{-3}$", "$10^{-4}$", "$10^{-5}$"])
    ax.xaxis.set_minor_formatter(NullFormatter())
    ax.set_xlim(2e-3, 2e-5)


def dt_values(rows, k, level, ykey):
    pts = sorted((float(r["dt"]), pv.num(r, ykey)) for r in rows
                 if pv.key(r) == k and int(r["level"]) == level
                 and float(r["dt"]) <= 0.001 + 1e-12 and r[ykey] != "")
    return [a for a, _ in pts], [b for _, b in pts]


def dt_velocity(plt, rows, out, dpi, png_only):
    fig, axes = plt.subplots(1, 4, figsize=(7.6, 2.7), sharey=True)
    for ax, level, lab in zip(axes, (1, 2, 3, 4), "abcd"):
        for k in METHODS:
            x, y = dt_values(rows, k, level, "u_L2_all_ivp")
            if x:
                ax.loglog(x, y, color=pv.COLOR[k], marker=pv.MARKER[k], ms=3.2,
                          ls=pv.DASH[k], lw=0.8, label=pv.LABEL[k])
        ax.set_xlabel("$\\Delta t$ (s)")
        ax.set_title(f"Level {level}", fontsize=8.5)
        pv.panel(ax, lab)
        ax.texts[-1].set_position((0.02, 0.97))
        tidy_dt_axis(ax)
    axes[0].set_ylabel("$L_2$ velocity error / $U_0$")
    h, l = axes[0].get_legend_handles_labels()
    fig.legend(h, l, loc="lower center", ncol=4, bbox_to_anchor=(0.5, -0.2), fontsize=6.8)
    pv.save(plt, fig, out, "controls_dt_velocity", dpi, png_only)


def dt_shear(plt, rows, out, dpi, png_only):
    fig, axes = plt.subplots(1, 3, figsize=(7.6, 2.8))
    spec = (("shear_amp_err_ivp", "amplitude error (fraction)", True),
            ("shear_phase_err_deg_ivp", "phase error (degrees)", False),
            ("shear_rms_err_ivp", "history RMS error (fraction)", True))
    for ax, (ykey, label, absval), lab in zip(axes, spec, "abc"):
        for k in METHODS + (("cutLink", "FO"),):
            x, y = dt_values(rows, k, 4, ykey)
            if x:
                y = [abs(v) for v in y] if absval else y
                ax.plot(x, y, color=pv.COLOR[k], marker=pv.MARKER[k], ms=3.2,
                        ls=pv.DASH[k], lw=0.8, label=pv.LABEL[k])
        ax.set_xscale("log")
        if absval:
            ax.set_yscale("log")
        tidy_dt_axis(ax)
        ax.set_xlabel("$\\Delta t$ (s)")
        ax.set_ylabel(("|" + label.split(" (")[0] + "| " if absval else "") +
                      ("(fraction)" if absval else label) if False else label)
        pv.panel(ax, lab)
    axes[0].set_title("Level 4, against the initial-value solution", fontsize=8.5, loc="left")
    h, l = axes[0].get_legend_handles_labels()
    fig.legend(h, l, loc="lower center", ncol=5, bbox_to_anchor=(0.5, -0.17), fontsize=6.6)
    pv.save(plt, fig, out, "controls_dt_shear", dpi, png_only)


def matched_dt(plt, rows, out, dpi, png_only):
    fig, axes = plt.subplots(1, 2, figsize=(7.0, 3.0), sharey=True)
    for ax, dt, lab in zip(axes, (0.001, DT_REF), "ab"):
        hs = [pv.h_of(lv) for lv in (1, 2, 3, 4)]
        for k in METHODS:
            ys = [pv.num(at(rows, k, lv, dt), "u_L2_all_ivp") if at(rows, k, lv, dt) else float("nan")
                  for lv in (1, 2, 3, 4)]
            ax.loglog(hs, ys, color=pv.COLOR[k], marker=pv.MARKER[k], ms=3.4,
                      ls=pv.DASH[k], lw=0.8, label=pv.LABEL[k])
        for slope, ls in ((1, ":"), (2, "-.")):
            ax.loglog(hs, [2e-2 * (h / hs[0]) ** slope for h in hs], color="0.55", ls=ls, lw=0.7,
                      label=f"slope {slope} (guide)")
        from matplotlib.ticker import NullFormatter
        ax.set_xticks(hs)
        ax.set_xticklabels([f"{h * 1e3:.2f}" for h in hs], fontsize=7)
        ax.xaxis.set_minor_formatter(NullFormatter())
        ax.set_xlabel("Cell height $h$ (mm)")
        ax.set_title(f"$\\Delta t$ = {dt:g} s", fontsize=8.5)
        pv.panel(ax, lab)
    axes[0].set_ylabel("$L_2$ velocity error / $U_0$")
    axes[0].legend(loc="lower right", fontsize=6.4)
    pv.save(plt, fig, out, "controls_matched_dt", dpi, png_only)


def alignment(plt, rows, out, dpi, png_only):
    fig, axes = plt.subplots(1, 2, figsize=(7.0, 3.0), sharey=True)
    for ax, dt, lab in zip(axes, (0.001, DT_REF), "ab"):
        for k in METHODS[1:]:
            lv = [1, 2, 3, 4]
            raw = [pv.num(at(rows, k, x, dt), "shear_amp_err_ivp") for x in lv]
            al = [pv.num(at(rows, k, x, dt), "shear_amp_err_al") for x in lv]
            ax.plot(lv, [100 * a for a in raw], color=pv.COLOR[k], marker=pv.MARKER[k], ms=3.4, ls="-",
                    lw=0.8, label=pv.LABEL[k] + ", raw solver output")
            ax.plot(lv, [100 * a for a in al], color=pv.COLOR[k], marker=pv.MARKER[k], ms=3.4,
                    mfc="white", ls="--", lw=0.8, label=pv.LABEL[k] + ", inertia aligned (post-processing)")
        ax.axhline(0, color="black", lw=0.6)
        ax.set_xticks([1, 2, 3, 4])
        ax.set_xlabel("Mesh level")
        ax.set_title(f"$\\Delta t$ = {dt:g} s", fontsize=8.5)
        pv.panel(ax, lab)
    axes[0].set_ylabel("Shear amplitude error (%)")
    h, l = axes[0].get_legend_handles_labels()
    fig.legend(h, l, loc="lower center", ncol=2, bbox_to_anchor=(0.5, -0.17), fontsize=6.4)
    pv.save(plt, fig, out, "controls_alignment", dpi, png_only)


def ghost_regions(plt, rows, out, dpi, png_only):
    reg = sc.read_rows(sc.CONTROLS / "velocity-regions.csv")
    fig, axes = plt.subplots(1, 2, figsize=(6.8, 3.0), sharey=True)
    for ax, dt, lab in zip(axes, (0.001, DT_REF), "ab"):
        for col, name, mk in (("L2_all", "all fluid rows", "o"),
                              ("L2_interface_row", "interface row (centre within $h$ of the wall)", "^"),
                              ("L2_near", "within one boundary-layer thickness", "s"),
                              ("L2_bulk", "bulk, beyond one thickness", "v")):
            sel = sorted((int(r["level"]), float(r[col])) for r in reg
                         if r["method"] == "ghostCell" and abs(float(r["dt"]) - dt) < 1e-12)
            ax.semilogy([a for a, _ in sel], [b for _, b in sel], marker=mk, ms=3.4, lw=0.8, label=name)
        ax.set_xticks([1, 2, 3, 4])
        ax.set_xlabel("Mesh level")
        ax.set_title(f"ghostCell, $\\Delta t$ = {dt:g} s", fontsize=8.5)
        pv.panel(ax, lab)
    axes[0].set_ylabel("$L_2$ velocity error / $U_0$")
    h, l = axes[0].get_legend_handles_labels()
    fig.legend(h, l, loc="lower center", ncol=2, bbox_to_anchor=(0.5, -0.17), fontsize=6.4)
    pv.save(plt, fig, out, "controls_ghost_regions", dpi, png_only)


def offset(plt, out, dpi, png_only):
    rows = sc.read_rows(sc.CONTROLS / "offset-audit.csv")
    fig, axes = plt.subplots(1, 2, figsize=(7.0, 3.0))
    ax = axes[0]
    for case, mk, name in (("tutorial", "o", "tutorial wall $y_w$ = 0.1234 m"),
                           ("wall on a mesh face", "s", "wall on a mesh face, $y_w$ = 0.125 m")):
        sel = [r for r in rows if r["case"] == case]
        ax.plot([float(r["predicted_offset_over_h"]) for r in sel],
                [float(r["fitted_ME_over_h"]) for r in sel], mk, ms=4.5, label=name)
        for r in sel:
            ax.annotate(f"L{r['level']}", (float(r["predicted_offset_over_h"]), float(r["fitted_ME_over_h"])),
                        textcoords="offset points", xytext=(4, 3), fontsize=6.5)
    ax.plot([-0.6, 0.4], [-0.6, 0.4], color="black", lw=0.6)
    ax.set_xlabel("Nearest mesh face minus wall, $/h$ (predicted)")
    ax.set_ylabel("Fitted $\\rho\\,\\Delta H/h$ of (ME $-$ forcing)")
    ax.legend(fontsize=6.2, loc="upper left")
    pv.panel(ax, "a")
    ax = axes[1]
    cases = {r["level"]: r for r in rows if r["case"] == "wall on a mesh face"}
    base = {int(r["level"]): r for r in sc.read_rows(sc.RESULTS / "cases.csv")
            if r["method"] == "cutLink" and r["estimator"] == "ME" and r["dt"] in ("", pv.TUTORIAL_DT)}
    base_fo = {int(r["level"]): r for r in sc.read_rows(sc.RESULTS / "cases.csv")
               if r["method"] == "cutLink" and r["estimator"] == "FO" and r["dt"] in ("", pv.TUTORIAL_DT)}
    lv = [2, 3, 4]
    w = 0.2
    series = (("ME, $y_w$ = 0.1234", [100 * pv.num(base[x], "shear_amp_err_ivp") for x in lv], "#0072B2"),
              ("forcing, $y_w$ = 0.1234", [100 * pv.num(base_fo[x], "shear_amp_err_ivp") for x in lv], "#7F7F7F"),
              ("ME, wall on a face", [100 * float(cases[str(x)]["amp_err_ME"]) for x in lv], "#56B4E9"),
              ("forcing, wall on a face", [100 * float(cases[str(x)]["amp_err_FO"]) for x in lv], "#CCCCCC"))
    for i, (name, ys, colr) in enumerate(series):
        ax.bar([x + (i - 1.5) * w for x in lv], ys, w, color=colr, edgecolor="black", lw=0.5, label=name)
    ax.axhline(0, color="black", lw=0.6)
    ax.set_xticks(lv)
    ax.set_xlabel("Mesh level")
    ax.set_ylabel("Shear amplitude error (%)")
    ax.set_ylim(-6.5, 4.0)
    ax.legend(fontsize=6, loc="lower center", ncol=2)
    pv.panel(ax, "b")
    pv.save(plt, fig, out, "controls_offset", dpi, png_only)


def cellavg(plt, out, dpi, png_only):
    rows = sc.read_rows(sc.CONTROLS / "cell-average.csv")
    fig, axes = plt.subplots(1, 2, figsize=(7.0, 2.9), sharey=True)
    names = ["cutLink", "penalty", "incremental", "ghostCell"]
    for ax, level, lab in zip(axes, (3, 4), "ab"):
        sel = {r["method"]: r for r in rows if int(r["level"]) == level}
        w = 0.2
        for i, (col, name, colr) in enumerate((("centre", "exact value at the cell centre, fluid rows", "#0072B2"),
                                               ("average", "exact cell average, fluid rows", "#D55E00"),
                                               ("centre_whole", "centre value, rows wholly in the fluid", "#56B4E9"),
                                               ("average_whole", "cell average, rows wholly in the fluid", "#E69F00"))):
            ax.bar([j + (i - 1.5) * w for j in range(4)], [float(sel[n][col]) for n in names], w,
                   color=colr, edgecolor="black", lw=0.5, label=name)
        ax.set_yscale("log")
        ax.set_ylim(3e-5, 3e-2)
        ax.set_xticks(range(4))
        ax.set_xticklabels(["cutLink", "penalty", "incremental", "ghostCell"], fontsize=7)
        ax.set_title(f"Level {level}, $\\Delta t$ = 0.001 s", fontsize=8.5)
        pv.panel(ax, lab)
    axes[0].set_ylabel("$L_2$ velocity error / $U_0$")
    h, l = axes[0].get_legend_handles_labels()
    fig.legend(h, l, loc="lower center", ncol=2, bbox_to_anchor=(0.5, -0.14), fontsize=6.4)
    pv.save(plt, fig, out, "controls_cellavg", dpi, png_only)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--out", type=Path, default=pv.VERIFICATION_DIR / "images")
    parser.add_argument("--dpi", type=int, default=200)
    parser.add_argument("--png-only", action="store_true")
    args = parser.parse_args(argv)
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib is required for the figures", file=sys.stderr)
        return 2
    pv.style(plt)
    rows = merged()
    for fn in (dt_velocity, dt_shear, matched_dt, alignment, ghost_regions):
        fn(plt, rows, args.out, args.dpi, args.png_only)
    offset(plt, args.out, args.dpi, args.png_only)
    cellavg(plt, args.out, args.dpi, args.png_only)
    print(f"figures written to {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
