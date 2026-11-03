#!/usr/bin/env python3
"""Figure 2: rotational verification and the combined-member error floor. Reads data/fig2_convergence.csv and data/fig2_timestep_strict.csv only.
Order lines are drawn ONLY where the observed order is resolved (torsion alpha_B, sharp interval); never as a fit elsewhere.
Usage: python3 fig2_rotational_convergence.py [DATA_DIR] [OUT_PREFIX]"""
import csv, os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams.update({"font.size": 9, "font.family": "DejaVu Sans"})
import matplotlib.pyplot as plt
HERE = os.path.dirname(os.path.abspath(__file__)); DATA = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "data"); OUT = sys.argv[2] if len(sys.argv) > 2 else "fig2_rotational_convergence"
OI = dict(blue="#0072B2", verm="#D55E00", green="#009E73", grey="#555555")
rows = list(csv.DictReader(open(os.path.join(DATA, "fig2_convergence.csv"))))
def get(f, mem, q):
    d = {r["mesh"]: r for r in rows if r["formulation"] == f and r["member"] == mem and r["qoi"] == q}
    return [float(d[m]["E_max_over_peak"]) for m in ("M1", "M2", "M3")], [float(d["M3"]["noise_over_peak"]) if d["M3"]["noise_over_peak"] else 0.0][0]
ntheta = np.array([64, 128, 256])
SER = [("A", "torsion", OI["blue"], "o", "-", "torsion, A (tutorial: leastSquares $\\nabla p$)"), ("A", "combined", OI["verm"], "s", "-", "combined, A (tutorial: leastSquares $\\nabla p$)"),
       ("B", "combined", OI["green"], "^", "--", "combined, B (optional: Gauss linear $\\nabla p$ only)")]
_p = os.path.join(DATA, "fig2_timestep_strict.csv")
has_dt = os.path.exists(_p) and {"P0", "P1", "P2"} <= {r["run"] for r in csv.DictReader(open(_p))}
fig, axs = plt.subplots(1, 3 if has_dt else 2, figsize=(12.6 if has_dt else 8.6, 4.5))
for ax, q, title in ((axs[0], "alpha_B", r"(a) wall angular displacement $\alpha_B$"), (axs[1], "T_B", r"(b) wall torque $T_B$")):
    for f, mem, col, mk, ls, lab in SER:
        e, noise = get(f, mem, q)
        ax.loglog(ntheta, e, color=col, marker=mk, ls=ls, lw=1.6, ms=6, label=lab)
        if noise > 0:
            lim = [i for i, v in enumerate(e) if v < 2 * noise]
            ax.errorbar(ntheta, e, yerr=noise, fmt="none", ecolor=col, elinewidth=1.0, capsize=2.5, alpha=0.8)
            for i in lim: ax.plot(ntheta[i], e[i], marker=mk, ms=6.5, mfc="white", mec=col, mew=1.6, ls="none", zorder=5)
    ax.set_title(title, fontsize=10, loc="left"); ax.set_xlabel(r"azimuthal cells $n_\theta$ (M1, M2, M3; $\Delta t=2.5$ ms)"); ax.set_xticks(ntheta); ax.set_xticklabels(["64\n(M1)", "128\n(M2)", "256\n(M3)"])
    ax.minorticks_off(); ax.grid(True, which="major", color="#dddddd", lw=0.6); ax.set_ylabel("max-norm error over 0–2 s / max|reference|")
e0, _ = get("A", "torsion", "alpha_B"); xs = np.array([64, 256]); axs[0].loglog(xs, 0.45 * e0[0] * (xs / 64.0) ** -2, color="#888888", lw=1.0, ls=":")
axs[0].text(0.97, 0.30, "dotted: slope −2 guide (offset);\nobserved torsion-A order 1.98 (sharp)", transform=axs[0].transAxes, fontsize=7, color="#777777", ha="right", va="top")
axs[1].text(0.03, 0.04, "bars: documented noise\n(torsion ±2.4e-4, combined ±1.3e-3)\nopen: error < 2× noise;\nno order drawn here", transform=axs[1].transAxes, fontsize=7, color=OI["grey"], va="bottom")
axs[0].text(0.03, 0.04, "noise bars (±5.7e-6) smaller than markers", transform=axs[0].transAxes, fontsize=7, color=OI["grey"])
for ax in axs[:2]:
    ax.set_xticks(ntheta); ax.set_xticklabels(["64\n(M1)", "128\n(M2)", "256\n(M3)"]); ax.minorticks_off(); ax.set_xlim(52, 310)
h, l = axs[0].get_legend_handles_labels(); fig.legend(h, l, loc="lower center", ncol=3, fontsize=7.8, frameon=False, bbox_to_anchor=(0.5, 0.075))
if has_dt:
    ax = axs[2]; d = list(csv.DictReader(open(os.path.join(DATA, "fig2_timestep_strict.csv"))))
    def v(run, q):
        r = [r for r in d if r["run"] == run and r["qoi"] == q][0]; return float(r["E_max_over_peak"]) / float(r["E_P0_same_window"])
    x = [0.625, 1.25, 2.5]
    for q, mk, col, lab in (("alpha_B", "o", OI["blue"], r"$\alpha_B$"), ("T_B", "s", OI["verm"], r"$T_B$"), ("H_f", "D", OI["green"], r"$H_f$")):
        y = [v("P2", q), v("P1", q), v("P0", q)]; ax.semilogx(x, np.array(y), marker=mk, color=col, lw=1.6, ms=6, label=lab + " (P0–P2)")
        if any(r["run"] == "P3" for r in d): ax.plot([2.5], [v("P3", q)], marker=mk, mfc="white", mec=col, mew=1.6, ls="none", ms=8)
    ax.axhline(1, color="#999999", lw=0.7); ax.set_ylim(0.5, 1.5); ax.set_xticks(x); ax.set_xticklabels(["0.625", "1.25", "2.5"]); ax.minorticks_off()
    ax.set_xlabel(r"time step $\Delta t$ (ms), combined M2, tutorial formulation A"); ax.set_ylabel(r"error / error at $\Delta t=2.5$ ms (P0)")
    ax.set_title("(c) time step and strict coupling (new runs, M2)", fontsize=10, loc="left"); ax.grid(True, color="#dddddd", lw=0.6)
    ax.plot([], [], marker="o", mfc="white", mec="k", ls="none", label="open: strict coupling P3 (aborted at 1.825 s; window 0–1.8225 s)")
    ax.legend(fontsize=7.2, frameon=False, loc="lower left")
fig.suptitle("Rotational verification of radialTorsionalRing: errors against the semi-analytical reference (provisional category A)", fontsize=10.5, x=0.01, ha="left", y=0.995)
fig.text(0.01, 0.005, "Formulation A = tutorial as stored (default leastSquares for p and U, stock newMovingWallVelocity); B = A with grad(p) Gauss linear only (optional; Layer-A acceptance not clean at M3). Runs: A: F2_M1T/M2T, F1_M3T, F3_M1C/M2C/M3C; B: C1–C3; "
         "P-runs: 2 ranks. Data: data/fig2_*.csv. Single time step per mesh: temporal order not established. Order lines only where resolved.", fontsize=6.6, color=OI["grey"], va="bottom", wrap=True)
fig.tight_layout(rect=[0, 0.12, 1, 0.965]); fig.savefig(OUT + ".png", dpi=150)
print("written", OUT + ".png")
