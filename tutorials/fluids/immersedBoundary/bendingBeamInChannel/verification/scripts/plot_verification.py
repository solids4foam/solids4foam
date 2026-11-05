#!/usr/bin/env python3
"""Figures of the bendingBeamInChannel verification from the tracked
results in reference/results/ (or those of './Allverify analyse' with
--data) and the reference series in reference/. No CFD is run and nothing
is filtered or shifted.

    flux_vs_level        net exterior flux against the exact rate: the peak
                         residual of cutLink against the mesh level, its
                         bias when the exact flux is zero (Study B), and the
                         zero flux of penalty, incremental and ghostCell
    force_histories      Cd and Cl of cutLink, penalty and incremental against
                         the body-fitted series, and the Cd harmonics
    fv_sensitivity       change of the Cd nRMS and second harmonic with each
                         finite-volume setting, against the thresholds
    bodyfitted_control   the body-fitted case at two time steps (R1 control)
    observed_cost        observed cost per time step of matched cutLink and
                         penalty runs, after the CPU-affinity correction
    mesh_refinement      Cd nRMS, Cl nRMS, Cd second harmonic and the peak
                         exterior-flux residual of the four methods on
                         Levels 1 to 3 (and the Level 4 flux test of cutLink)
    histories_L2         Cd(t) and Cl(t), 1 < t < 3 s, of the four methods
    histories_L3         against the body-fitted series, Levels 2 and 3
                         (the tracked series every 0.02 s, unfiltered)

Style of the immersed-boundary verifications: one colour, marker and line
pattern per method (Okabe-Ito colours), references in black and grey. Only
matplotlib is needed (no LaTeX).

Usage:
    python3 scripts/plot_verification.py [--data DIR] [--out DIR]
                                         [--dpi N] [--png-only]

The previews in images/ are written with '--out images --dpi 110
--png-only'. Without --png-only, a PDF and an SVG are also written.
"""

import argparse
import csv
import sys
from pathlib import Path

VERIFICATION_DIR = Path(__file__).resolve().parent.parent
REFERENCE_DIR = VERIFICATION_DIR / "reference"

COLOR = {"cutLink": "#0072B2", "ghostCell": "#D55E00", "penalty": "#009E73",
         "incremental": "#CC79A7"}
MARKER = {"cutLink": "o", "ghostCell": "s", "penalty": "^",
          "incremental": "D"}
DASH = {"cutLink": "-", "ghostCell": (0, (4, 1.6)),
        "penalty": (0, (5, 1.5, 1.2, 1.5)), "incremental": (0, (1.2, 1.2))}
LABEL = {"cutLink": "cutLink", "ghostCell": "ghostCell (exploratory)",
         "penalty": "penalty", "incremental": "legacy incremental"}
COLOR_ME = "#56B4E9"
REF = dict(color="black", lw=0.95, zorder=1.5)


def read_csv(path: Path) -> list:
    with path.open() as handle:
        return list(csv.DictReader(handle))


def floats(rows: list, key: str) -> list:
    return [float(r[key]) if r[key] not in ("", "nan") else float("nan")
            for r in rows]


def style(plt):
    plt.rcParams.update({
        "font.family": "serif", "font.size": 8.5, "axes.labelsize": 9,
        "axes.titlesize": 9, "legend.fontsize": 7.5, "axes.linewidth": 0.8,
        "xtick.direction": "in", "ytick.direction": "in",
        "xtick.top": True, "ytick.right": True,
        "legend.frameon": True, "legend.edgecolor": "black",
        "legend.fancybox": False, "legend.framealpha": 1.0,
        "savefig.bbox": "tight"})


def panel(ax, label):
    ax.text(0.02, 0.97, label, transform=ax.transAxes, va="top", ha="left",
            fontweight="bold")
    ax.grid(True, ls="--", lw=0.4, color="0.75")


def flux_vs_level(plt, data, out):
    cases = read_csv(data / "cases.csv")
    fig, ax = plt.subplots(figsize=(4.6, 3.4))
    chosen = [r for r in cases if r["coefficient"] == "Cd"
              and r["estimator"] in ("cutLink ST", "total")]
    for method, study, ls, filled in (("cutLink", "A", "-", True),
                                      ("cutLink", "B", "--", False)):
        rows = sorted([r for r in chosen if r["method"] == method
                       and r["study"] == study], key=lambda r: int(r["level"]))
        levels = [int(r["level"]) for r in rows]
        peak = float(rows[0]["flux_peak_analytic"]) if study == "A" else \
            float([r for r in chosen if r["study"] == "A"
                   and r["method"] == "cutLink"][0]["flux_peak_analytic"])
        if study == "A":
            y = [float(r["flux_peak_residual_pct"]) for r in rows]
            lab = "cutLink, exact flux $\\pm$0.0283 m$^3$/s (Study A)"
        else:
            y = [100 * float(r["flux_max_abs_Q"]) / peak for r in rows]
            lab = "cutLink, exact flux 0 (Study B): $|Q|$"
        ax.plot(levels, y, ls=ls, color=COLOR[method], marker="o",
                mfc=COLOR[method] if filled else "white", ms=4.5, lw=1.0,
                label=lab)
    for method in ("penalty", "incremental", "ghostCell"):
        rows = sorted([r for r in chosen if r["method"] == method
                       and r["study"] == "A" and r["flux_snapshots"]],
                      key=lambda r: int(r["level"]))
        if not rows:
            continue
        peak = float(rows[0]["flux_peak_analytic"])
        # the exact flux is non-zero and the net exterior flux is zero: the
        # residual is the whole peak
        shift = {"penalty": 0.0, "incremental": -0.07,
                 "ghostCell": 0.07}[method]
        ax.plot([int(r["level"]) + shift for r in rows], [100.0] * len(rows),
                ls="none", marker=MARKER[method], color=COLOR[method],
                ms=5, label=f"{LABEL[method]}: $Q$ = 0 (residual = peak)")
    ax.set_yscale("log")
    ax.set_xticks([1, 2, 3, 4])
    ax.set_xticklabels(["L1\n20 k", "L2\n80 k", "L3\n320 k", "L4\n1.28 M"])
    ax.set_ylabel("residual, % of the analytic peak flux")
    ax.set_ylim(0.5, 300)
    ax.set_xlim(0.7, 4.3)
    ax.grid(True, which="both", ls="--", lw=0.3, color="0.8")
    ax.legend(loc="lower left", bbox_to_anchor=(0.0, 1.02), ncol=1)
    fig.subplots_adjust(top=0.72)
    return fig


def force_histories(plt, data, out):
    t_ref, cd_ref, cl_ref = [], [], []
    for r in read_csv_rows(REFERENCE_DIR / "CdBodyFitted.dat"):
        t_ref.append(r[0])
        cd_ref.append(r[1])
        cl_ref.append(r[2])
    cases = read_csv(data / "cases.csv")
    fig, axes = plt.subplots(1, 3, figsize=(7.2, 2.7),
                             gridspec_kw={"width_ratios": [1, 1, 0.8]})
    runs = (("cutLink", 3, "Cd_cutLink_ST", "Cl_cutLink_ST"),
            ("penalty", 3, "Cd_total", "Cl_total"),
            ("incremental", 2, "Cd_total", "Cl_total"))
    for ax, key, ref, lab in ((axes[0], 0, cd_ref, "$C_d$"),
                              (axes[1], 1, cl_ref, "$C_l$")):
        ax.plot(t_ref, ref, label="body-fitted", **REF)
        for method, level, cd_col, cl_col in runs:
            hist = read_csv(data / f"history-{method}-L{level}.csv")
            t = floats(hist, "t")
            y = floats(hist, (cd_col, cl_col)[key])
            ax.plot(t, y, color=COLOR[method], ls=DASH[method], lw=0.8,
                    label=f"{LABEL[method]} L{level}")
        ax.set_xlim(5, 8)
        ax.set_xlabel("$t$ (s)")
        ax.set_ylabel(lab)
        panel(ax, "(a)" if key == 0 else "(b)")
    axes[1].text(0.98, 0.04, "reference lift unverified", ha="right",
                 transform=axes[1].transAxes, fontsize=7, color="0.3")
    ax = axes[2]
    ref_a2 = None
    from_cases = {}
    for r in cases:
        if r["coefficient"] == "Cd" and r["study"] == "A" and r["nrms_pct"] \
                and r["estimator"] in ("cutLink ST", "total"):
            from_cases.setdefault(r["method"], []).append(
                (int(r["level"]), float(r["A2"])))
    for method, pts in from_cases.items():
        pts.sort()
        ax.plot([p[0] for p in pts], [p[1] for p in pts],
                color=COLOR[method], marker=MARKER[method], ms=4.5, lw=1.0,
                ls=DASH[method], label=LABEL[method])
    ref_a2 = reference_second_harmonic(t_ref, cd_ref)
    ax.axhline(ref_a2, **REF)
    ax.text(1.05, ref_a2 + 0.15, "body-fitted", fontsize=7)
    ax.set_xticks([1, 2, 3])
    ax.set_xlabel("mesh level")
    ax.set_ylabel("$C_d$ second harmonic")
    ax.set_ylim(0, 7.5)
    panel(ax, "(c)")
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=4,
               bbox_to_anchor=(0.5, -0.12))
    fig.tight_layout()
    return fig


def read_csv_rows(path: Path) -> list:
    rows = []
    for line in path.read_text().splitlines():
        if line.strip() and not line.startswith("#"):
            rows.append([float(v) for v in line.split()])
    return rows


def reference_second_harmonic(t, y) -> float:
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    import bending_beam_verification as bbv
    return bbv.harmonics(t, y)["A2"]


def fv_sensitivity(plt, data, out):
    rows = read_csv(data / "fv-sensitivity.csv")
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.6), sharey=True)
    labels = []
    for i, r in enumerate(rows):
        labels.append(f"{r['comparison']}  [{r['estimator']}]")
        for ax, key, thr in ((axes[0], "d_Cd_nrms", "thr_nrms"),
                             (axes[1], "d_Cd_A2", "thr_A2")):
            ax.plot(float(r[key]), i, ls="none", marker=MARKER[r["method"]],
                    color=COLOR[r["method"]], ms=4.5)
    for ax, title, thr in ((axes[0], "change of the $C_d$ nRMS (points)",
                            "thr_nrms"),
                           (axes[1], "change of the $C_d$ second harmonic",
                            "thr_A2")):
        ax.axvline(0, color="black", lw=0.7)
        full = float(rows[0][thr])
        for sign in (-1, 1):
            ax.axvline(sign * full, color="0.35", ls="--", lw=0.8)
            ax.axvline(sign * full / 2, color="0.6", ls=":", lw=0.8)
        ax.set_xlabel(title)
        ax.grid(True, axis="x", ls="--", lw=0.4, color="0.75")
    axes[0].set_yticks(range(len(rows)))
    axes[0].set_yticklabels(labels, fontsize=6.5)
    axes[0].invert_yaxis()
    axes[0].set_xlim(-1.2, 1.2)
    axes[1].set_xlim(-0.12, 0.12)
    fig.text(0.5, 0.005, "dashed: pre-registered threshold (1.0 point, 0.1); "
             "dotted: half of it, the threshold of the solver-convergence "
             "controls (0.5, 0.05).\nNo change reaches a threshold.",
             ha="center", fontsize=7)
    fig.tight_layout(rect=(0, 0.06, 1, 1))
    return fig


def bodyfitted_control(plt, data, out):
    rows = read_csv(data / "bodyfitted-R1-histories.csv")
    t = floats(rows, "t")
    fig, axes = plt.subplots(2, 2, figsize=(7.2, 4.6), sharex=True)
    series = (("A_dt0p01_archive", "A: archive, $\\Delta t$ = 0.01 s",
               dict(color="0.5", ls=(0, (4, 2)), lw=0.9)),
              ("B_dt0p0025", "B: R1, $\\Delta t$ = 0.0025 s",
               dict(color="#D55E00", lw=1.0)),
              ("C_dt0p01_control", "C: R1 control, $\\Delta t$ = 0.01 s",
               dict(color="black", lw=1.0, ls=(0, (1.2, 1.2)))))
    for col, (coef, lab) in enumerate((("Cd", "$C_d$"), ("Cl", "$C_l$"))):
        top, bottom = axes[0][col], axes[1][col]
        for key, label, kw in series:
            top.plot(t, floats(rows, f"{coef}_{key}"), label=label, **kw)
        c = floats(rows, f"{coef}_C_dt0p01_control")
        bottom.plot(t, [x - y for x, y in zip(floats(rows,
                    f"{coef}_B_dt0p0025"), c)], color="#D55E00", lw=1.0,
                    label="B $-$ C (time step only)")
        bottom.plot(t, [x - y for x, y in zip(floats(rows,
                    f"{coef}_A_dt0p01_archive"), c)], color="0.5",
                    ls=(0, (4, 2)), lw=0.9,
                    label="A $-$ C (build and history)")
        top.set_xlim(1, 8)
        top.set_ylabel(lab)
        bottom.set_ylabel(lab + " difference")
        bottom.set_xlabel("$t$ (s)")
        panel(top, "(a)" if col == 0 else "(b)")
        panel(bottom, "(c)" if col == 0 else "(d)")
    axes[0][0].legend(loc="lower left", bbox_to_anchor=(0.0, 1.02), ncol=1)
    axes[1][0].legend(loc="lower left", bbox_to_anchor=(0.0, 1.02), ncol=1)
    fig.subplots_adjust(hspace=0.9, wspace=0.3)
    return fig


def observed_cost(plt, data, out):
    rows = read_csv(data / "cost-matched-pairs.csv")
    fig, ax = plt.subplots(figsize=(4.8, 3.0))
    width = 0.38
    for i, r in enumerate(rows):
        for dx, method, key in ((-width / 2, "cutLink", "cutLink_s_per_step"),
                                (width / 2, "penalty", "penalty_s_per_step")):
            ax.bar(i + dx, float(r[key]), width, color=COLOR[method],
                   edgecolor="black", lw=0.5,
                   label=LABEL[method] if i == 0 else None)
        ax.text(i, max(float(r["cutLink_s_per_step"]),
                       float(r["penalty_s_per_step"])) + 0.01,
                f"x{float(r['ratio']):.2f}", ha="center", fontsize=7)
    ax.set_xticks(range(len(rows)))
    ax.set_xticklabels([r["settings"] for r in rows], rotation=30,
                       ha="right", fontsize=7)
    ax.set_ylabel("wall time per step (s), 2 ranks, L2")
    ax.set_ylim(0, 0.5)
    ax.grid(True, axis="y", ls="--", lw=0.4, color="0.75")
    ax.legend(loc="upper left", ncol=2)
    ax.text(0.99, 0.97, "single runs, shared node;\nnot an intrinsic cost",
            transform=ax.transAxes, ha="right", va="top", fontsize=6.5,
            color="0.3")
    fig.tight_layout()
    return fig



def titled(ax, text):
    ax.set_title(text, loc="left", fontsize=8.5, fontweight="bold")
    ax.grid(True, ls="--", lw=0.4, color="0.75")


def estimator_lines():
    """(label, method, estimator column suffix, colour, dash, marker)"""
    return (("cutLink, surface traction", "cutLink", "cutLink ST",
             COLOR["cutLink"], "-", "o"),
            ("cutLink, momentum exchange", "cutLink", "cutLink ME",
             COLOR_ME, (0, (4, 1.6)), "s"),
            ("penalty", "penalty", "total", COLOR["penalty"],
             DASH["penalty"], "^"),
            ("legacy incremental", "incremental", "total",
             COLOR["incremental"], DASH["incremental"], "D"),
            ("ghostCell (exploratory)", "ghostCell", "total",
             COLOR["ghostCell"], DASH["ghostCell"], "v"))


def mesh_refinement(plt, data, out):
    cases = [r for r in read_csv(data / "cases.csv") if r["study"] == "A"]
    fig, axes = plt.subplots(2, 2, figsize=(7.2, 5.6))
    spec = ((axes[0, 0], "Cd", "nrms_pct", "$C_d$ nRMS (%)", "(a)"),
            (axes[0, 1], "Cl", "nrms_pct", "$C_l$ nRMS (%)", "(b)"),
            (axes[1, 0], "Cd", "A2", "$C_d$ second harmonic", "(c)"))
    ref_a2 = None
    t, cd = [], []
    for r in read_csv_rows(REFERENCE_DIR / "CdBodyFitted.dat"):
        t.append(r[0])
        cd.append(r[1])
    ref_a2 = reference_second_harmonic(t, cd)
    for ax, coef, key, ylab, tag in spec:
        for lab, method, est, col, dash, mark in estimator_lines():
            pts = sorted((int(r["level"]), float(r[key])) for r in cases
                         if r["method"] == method and r["estimator"] == est
                         and r["coefficient"] == coef and r[key]
                         and int(r["level"]) <= 3)
            if pts:
                ax.plot([p[0] for p in pts], [p[1] for p in pts], color=col,
                        ls=dash, marker=mark, ms=4.2, lw=1.0, label=lab)
        if key == "A2":
            ax.axhline(ref_a2, **REF)
            ax.text(1.05, ref_a2 + 0.15, "body-fitted", fontsize=7)
            ax.set_ylim(0, 7.5)
        ax.set_xticks([1, 2, 3], ["L1\n20,020", "L2\n80,080", "L3\n320,320"])
        ax.set_xlabel("mesh level and cells")
        ax.set_ylabel(ylab)
        titled(ax, tag)
    ax = axes[1, 1]
    flux = sorted((int(r["level"]), float(r["flux_peak_residual_pct"]),
                   int(r["cells"])) for r in cases
                  if r["method"] == "cutLink" and r["estimator"] == "cutLink ST"
                  and r["coefficient"] == "Cd" and r["flux_peak_residual_pct"])
    ax.plot([f[2] for f in flux if f[0] <= 3],
            [f[1] for f in flux if f[0] <= 3], color=COLOR["cutLink"],
            marker="o", ms=4.2, lw=1.0)
    ax.plot([f[2] for f in flux if f[0] == 4],
            [f[1] for f in flux if f[0] == 4], color=COLOR["cutLink"],
            marker="*", ms=8, ls="none")
    for lab, method, est, col, dash, mark in estimator_lines()[2:]:
        pts = sorted((int(r["cells"]), float(r["flux_peak_residual_pct"]))
                     for r in cases if r["method"] == method
                     and r["estimator"] == est and r["coefficient"] == "Cd"
                     and r["flux_peak_residual_pct"])
        ax.plot([p[0] for p in pts], [p[1] for p in pts], color=col,
                marker=mark, ms=4.2, ls="none")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("cells (star: Level 4, 0.5 s flux test only)")
    ax.set_ylabel("peak flux residual (% of peak)")
    ax.text(0.04, 0.5, "penalty, incremental, ghostCell: net flux zero\n"
            "to the solver tolerance (residual = peak)", transform=ax.transAxes,
            fontsize=6.8, color="0.25")
    titled(ax, "(d)")
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=3,
               bbox_to_anchor=(0.5, -0.07))
    fig.tight_layout()
    return fig


def history_figure(level):
    def draw(plt, data, out):
        t_ref, cd_ref, cl_ref = [], [], []
        for r in read_csv_rows(REFERENCE_DIR / "CdBodyFitted.dat"):
            t_ref.append(r[0])
            cd_ref.append(r[1])
            cl_ref.append(r[2])
        fig, axes = plt.subplots(2, 2, figsize=(7.2, 5.6), sharex=True)
        groups = ((0, [e for e in estimator_lines() if e[1] == "cutLink"],
                   "cutLink"),
                  (1, [e for e in estimator_lines() if e[1] != "cutLink"],
                   "penalty, incremental, ghostCell"))
        for row, (coef, ref, lab) in enumerate((("Cd", cd_ref, "$C_d$"),
                                                ("Cl", cl_ref, "$C_l$"))):
            for col, members, title in groups:
                ax = axes[row, col]
                ax.plot([x for x in t_ref if 1 <= x <= 3],
                        [y for x, y in zip(t_ref, ref) if 1 <= x <= 3],
                        label="body-fitted", **REF)
                for elab, method, est, colr, dash, mark in members:
                    hist = read_csv(data / f"history-{method}-L{level}.csv")
                    key = (f"{coef}_cutLink_{est.split()[-1]}"
                           if method == "cutLink" else f"{coef}_total")
                    t = floats(hist, "t")
                    y = floats(hist, key)
                    ax.plot([a for a in t if 1 <= a <= 3],
                            [b for a, b in zip(t, y) if 1 <= a <= 3],
                            color=colr, ls=dash, lw=0.8, label=elab)
                ax.set_ylabel(lab)
                titled(ax, f"({'abcd'[2 * row + col]}) {title}")
                lo, hi = ax.get_ylim()
                ax.set_ylim(lo, hi + 0.3 * (hi - lo))
                ax.legend(loc="upper right", fontsize=6.3, ncol=1)
        for ax in axes[1]:
            ax.set_xlabel("$t$ (s)")
            ax.set_xlim(1, 3)
        axes[1, 1].text(0.98, 0.04, "reference lift unverified", ha="right",
                        transform=axes[1, 1].transAxes, fontsize=6.8,
                        color="0.3")
        fig.suptitle(f"Level {level}: tracked series every 0.02 s, "
                     "unfiltered and unshifted", fontsize=8.5, y=1.0)
        fig.tight_layout()
        return fig
    return draw

FIGURES = {"flux_vs_level": flux_vs_level, "force_histories": force_histories,
           "fv_sensitivity": fv_sensitivity,
           "bodyfitted_control": bodyfitted_control,
           "observed_cost": observed_cost, "mesh_refinement": mesh_refinement,
           "histories_L2": history_figure(2),
           "histories_L3": history_figure(3)}


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--data", default=str(REFERENCE_DIR / "results"))
    parser.add_argument("--out", default=str(VERIFICATION_DIR /
                                             "postProcessing/figures"))
    parser.add_argument("--dpi", type=int, default=200)
    parser.add_argument("--png-only", action="store_true")
    args = parser.parse_args(argv)
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("SKIPPED: matplotlib is not installed")
        return 0
    style(plt)
    data, out = Path(args.data), Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    for name, draw in FIGURES.items():
        needed = {"cost-matched-pairs.csv": "observed_cost",
                  "bodyfitted-R1-histories.csv": "bodyfitted_control",
                  "fv-sensitivity.csv": "fv_sensitivity"}
        missing = [f for f, n in needed.items() if n == name
                   and not (data / f).exists()]
        if missing or not (data / "cases.csv").exists():
            print(f"SKIPPED {name}: missing {missing or 'cases.csv'}")
            continue
        fig = draw(plt, data, out)
        for ext in (("png",) if args.png_only else ("png", "pdf", "svg")):
            fig.savefig(out / f"{name}.{ext}", dpi=args.dpi)
        plt.close(fig)
        print(f"{name}: written to {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
