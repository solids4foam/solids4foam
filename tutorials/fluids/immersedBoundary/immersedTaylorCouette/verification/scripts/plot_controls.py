#!/usr/bin/env python3
"""Figures of the controls of the immersedTaylorCouette verification (mesh
and time step separated, estimators, steady state) from the tracked results in
reference/results/ and reference/results/controls/. Nothing is smoothed,
filtered or shifted.

    controls_torque_vs_mesh   torque error against the mesh at each step cap
                              (the tutorial step halves with the mesh level)
    controls_torque_vs_dt     torque error against the step cap on each mesh
    controls_velocity_mesh    azimuthal-velocity L2 error against the mesh at
                              two step caps and at the tutorial step
    controls_estimators       the cutLink estimators, the torque from the
                              velocity field and the audit recipe, per level
    controls_steady           torque windows and drift of the t = 20 s runs
    controls_histories        torque histories to 20 s of the cases that are
                              not steady at 10 s

Usage: python3 scripts/plot_controls.py [--out DIR] [--dpi N] [--png-only]
"""

import argparse
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import plot_verification as pv  # noqa: E402
import taylor_couette_controls as tc  # noqa: E402

METHODS = (("cutLink", "ME"), ("cutLink", "ST"), ("penalty", ""), ("incremental", ""), ("ghostCell", ""))
CAPS = (0.01, 0.005, 0.0025, 0.00125)
CAP_STYLE = {0.0: ("0.0", "x", "-"), 0.01: ("#0072B2", "o", "-"), 0.005: ("#009E73", "s", "--"),
             0.0025: ("#D55E00", "^", "-."), 0.00125: ("#CC79A7", "D", ":")}


def cap_label(cap):
    return "tutorial step" if cap == 0 else f"cap {cap:g} s"


def rows_of(rows, k, cap):
    out = {}
    for r in rows:
        if pv.key(r) == k and abs(float(r["dt_cap"] or 0) - cap) < 1e-12:
            out[int(r["level"])] = r
    return out


def torque_vs_mesh(plt, rows, out, dpi, png_only):
    fig, axes = plt.subplots(2, 3, figsize=(7.6, 4.6), sharex=True)
    for ax, k, lab in zip(axes.flat, METHODS, "abcde"):
        for cap in (0.0,) + CAPS:
            d = rows_of(rows, k, cap)
            if len(d) < 2:
                continue
            lv = sorted(d)
            colr, mk, ls = CAP_STYLE[cap]
            ax.plot(lv, [100 * pv.num(d[x], "T_relerr") for x in lv], color=colr, marker=mk, ms=3.4,
                    ls=ls, lw=0.8, label=cap_label(cap))
        ax.axhline(0, color="black", lw=0.6)
        ax.set_xticks([1, 2, 3, 4])
        ax.set_title(pv.LABEL[k], fontsize=8)
        pv.panel(ax, lab)
    for ax in axes[1]:
        ax.set_xlabel("Mesh level")
    for ax in axes[:, 0]:
        ax.set_ylabel("Rotor torque error (%)")
    axes[0, 0].set_ylim(-1, 3.5)
    axes[0, 1].set_ylim(-4.0, 0.5)
    axes[0, 2].set_ylim(-26, 2)
    axes[1, 0].set_ylim(-8, 9)
    axes[1, 1].set_ylim(-2, 3)
    h, l = axes[0, 2].get_legend_handles_labels()
    axes[1, 2].axis("off")
    axes[1, 2].legend(h, l, loc="center", fontsize=7)
    pv.save(plt, fig, out, "controls_torque_vs_mesh", dpi, png_only)


def torque_vs_dt(plt, rows, out, dpi, png_only):
    fig, axes = plt.subplots(2, 3, figsize=(7.6, 4.6), sharex=True)
    lcol = {1: "#0072B2", 2: "#009E73", 3: "#D55E00", 4: "#CC79A7"}
    for ax, k, lab in zip(axes.flat, METHODS, "abcde"):
        for lv in (1, 2, 3, 4):
            xs, ys = [], []
            for cap in CAPS:
                r = rows_of(rows, k, cap).get(lv)
                if r:
                    xs.append(cap)
                    ys.append(100 * pv.num(r, "T_relerr"))
            if xs:
                ax.plot(xs, ys, color=lcol[lv], marker="o", ms=3.4, lw=0.8, label=f"Level {lv}")
        ax.axhline(0, color="black", lw=0.6)
        ax.set_xscale("log")
        ax.set_title(pv.LABEL[k], fontsize=8)
        from matplotlib.ticker import NullFormatter
        ax.set_xticks([0.01, 0.005, 0.0025, 0.00125])
        ax.set_xticklabels(["0.01", "0.005", "0.0025", "0.00125"], fontsize=6.5, rotation=30)
        ax.xaxis.set_minor_formatter(NullFormatter())
        ax.invert_xaxis()
        pv.panel(ax, lab)
    for ax in axes[1]:
        ax.set_xlabel("maxDeltaT cap (s)")
    for ax in axes[:, 0]:
        ax.set_ylabel("Rotor torque error (%)")
    h, l = axes[0, 0].get_legend_handles_labels()
    axes[1, 2].axis("off")
    axes[1, 2].legend(h, l, loc="center", fontsize=7)
    pv.save(plt, fig, out, "controls_torque_vs_dt", dpi, png_only)


def velocity_mesh(plt, rows, out, dpi, png_only):
    fig, axes = plt.subplots(1, 3, figsize=(7.6, 2.9), sharey=True)
    hs = [pv.h_of(lv) for lv in (1, 2, 3, 4)]
    for ax, cap, lab in zip(axes, (0.0, 0.01, 0.005), "abc"):
        for k in METHODS:
            if k == ("cutLink", "ST"):
                continue
            d = rows_of(rows, k, cap)
            lv = sorted(d)
            if len(lv) < 2:
                continue
            ax.loglog([hs[x - 1] for x in lv], [pv.num(d[x], "ut_L2") for x in lv], color=pv.COLOR[k],
                      marker=pv.MARKER[k], ms=3.4, ls=pv.DASH[k], lw=0.8, label=pv.LABEL[k])
        from matplotlib.ticker import NullFormatter
        ax.set_xticks(hs)
        ax.set_xticklabels([f"{h * 1e3:.2f}" for h in hs], fontsize=6.5)
        ax.xaxis.set_minor_formatter(NullFormatter())
        ax.set_xlabel("Cell width $h$ (mm)")
        ax.set_title(cap_label(cap), fontsize=8.5)
        pv.panel(ax, lab)
    axes[0].set_ylabel("$L_2$ azimuthal-velocity error / $\\omega R_1$")
    h, l = axes[0].get_legend_handles_labels()
    fig.legend(h, l, loc="lower center", ncol=4, bbox_to_anchor=(0.5, -0.1), fontsize=6.6)
    pv.save(plt, fig, out, "controls_velocity_mesh", dpi, png_only)


def estimators(plt, rows, out, dpi, png_only):
    audit = tc.read_rows(tc.CONTROLS / "st-audit.csv")
    fig, ax = plt.subplots(figsize=(5.6, 3.0))
    lv = [1, 2, 3, 4]
    me, st = rows_of(rows, ("cutLink", "ME"), 0.0), rows_of(rows, ("cutLink", "ST"), 0.0)
    series = (("momentum exchange (solver)", [100 * pv.num(me[x], "T_relerr") for x in lv], "#0072B2", "o", "-"),
              ("torque from the velocity field", [100 * pv.num(me[x], "T_field_rel") for x in lv], "#009E73", "s", "--"),
              ("surface traction (solver)", [100 * pv.num(st[x], "T_relerr") for x in lv], "#56B4E9", "v", "-"),
              ("surface-traction recipe, numerical field", [float(r["recipe_numerical_field_pct"]) for r in audit], "#56B4E9", "x", ":"),
              ("surface-traction recipe, exact field", [float(r["recipe_exact_field_circle_wall_pct"]) for r in audit], "#D55E00", "^", "-."),
              ("exact field in the first fluid layer only", [float(r["first_fluid_layer_replaced_by_exact_pct"]) for r in audit], "#CC79A7", "D", ":"))
    for name, ys, colr, mk, ls in series:
        ax.plot(lv, ys, color=colr, marker=mk, ms=3.6, ls=ls, lw=0.8, label=name)
    ax.axhline(0, color="black", lw=0.6)
    ax.set_xticks(lv)
    ax.set_xlabel("Mesh level (tutorial step)")
    ax.set_ylabel("Torque difference from exact (%)")
    pv.panel(ax, "")
    h, l = ax.get_legend_handles_labels()
    fig.legend(h, l, loc="upper center", ncol=2, bbox_to_anchor=(0.5, -0.02), fontsize=6.4)
    pv.save(plt, fig, out, "controls_estimators", dpi, png_only)


def steady(plt, out, dpi, png_only):
    rows = tc.read_rows(tc.CONTROLS / "steady-extended.csv")
    rows.sort(key=lambda r: r["case"])
    fig, axes = plt.subplots(1, 2, figsize=(7.6, 4.0), sharey=True)
    names = [r["case"].replace("TC-", "").replace("-np1", "").replace("-t20", "") for r in rows]
    y = list(range(len(rows)))
    for ax, key, title, lab in ((axes[0], "drift_13_15_to_18_20", "mean torque, 13-15 s against 18-20 s", "a"),
                                (axes[1], "dU_15_20", "velocity change, 15 s to 20 s ($/\\omega R_1$)", "b")):
        vals = [float(r[key]) for r in rows]
        colr = ["#D55E00" if r["steady_by_criterion"] == "0" else "#0072B2" for r in rows]
        ax.barh(y, vals, color=colr, edgecolor="black", lw=0.4)
        ax.axvline(1e-3, color="black", ls="--", lw=0.8)
        ax.set_xscale("log")
        ax.set_xlabel(title)
        pv.panel(ax, lab)
    axes[0].set_yticks(y)
    axes[0].set_yticklabels(names, fontsize=5.8)
    axes[0].invert_yaxis()
    axes[1].text(1.2e-3, 0.5, "threshold 1e-3", fontsize=6.4, rotation=90, va="top")
    pv.save(plt, fig, out, "controls_steady", dpi, png_only)


def histories(plt, out, dpi, png_only):
    rows = tc.read_rows(tc.CONTROLS / "history-extended.csv")
    cases = sorted({r["case"] for r in rows})
    fig, axes = plt.subplots(len(cases), 1, figsize=(7.0, 1.7 * len(cases)), sharex=True)
    exact = float(tc.read_rows(tc.RESULTS / "cases.csv")[0]["T_exact"])
    for ax, c in zip(axes, cases):
        sel = [r for r in rows if r["case"] == c]
        ax.plot([float(r["time"]) for r in sel], [float(r["rotor_torque_z"]) / exact for r in sel],
                color="#0072B2", lw=0.5, label=c.replace("TC-", "").replace("-np1", "").replace("-t20", ""))
        ax.axhline(1.0, color="black", lw=0.7)
        ax.legend(loc="upper right", fontsize=6.4)
        ax.grid(True, ls="--", lw=0.4, color="0.75")
        ax.set_ylabel("$T/T_{exact}$")
    axes[-1].set_xlabel("Time (s)")
    pv.save(plt, fig, out, "controls_histories", dpi, png_only)


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
    rows = tc.merged_cases()
    for fn in (torque_vs_mesh, torque_vs_dt, velocity_mesh, estimators):
        fn(plt, rows, args.out, args.dpi, args.png_only)
    steady(plt, args.out, args.dpi, args.png_only)
    if (tc.CONTROLS / "history-extended.csv").exists():
        histories(plt, args.out, args.dpi, args.png_only)
    print(f"figures written to {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
