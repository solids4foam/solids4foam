#!/usr/bin/env python3
"""Figure 3: why the viscous normal wall stress is NOT verified, and why the total normal traction still looks accurate.
Sign convention (transmitted normal traction, outward normal n of the solid): sigma_rr = -p_b + sigma_nv  (closes to 1e-10 in every run).
Matched data only, combined member, M2, dt = 2.5 ms: formulation A (tutorial) and B (A with Gauss-linear grad(p)). Reads data/fig3_traction_timeseries.csv.
No continuity-imposed (Cont/ContX) data are used anywhere in this figure."""
import csv, os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams.update({"font.size": 9, "font.family": "DejaVu Sans"})
import matplotlib.pyplot as plt
HERE = os.path.dirname(os.path.abspath(__file__)); DATA = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "data"); OUT = sys.argv[2] if len(sys.argv) > 2 else "fig3_normal_traction"
OI = dict(blue="#0072B2", verm="#D55E00", green="#009E73", grey="#555555", black="#222222")
R = list(csv.DictReader(open(os.path.join(DATA, "fig3_traction_timeseries.csv"))))
def col(f, k): return np.array([float(r[k]) for r in R if r["formulation"] == f])
t = col("A", "t")
fig, axs = plt.subplots(1, 3, figsize=(12.6, 4.7))
# (a) viscous normal stress, reference vs CFD (A)
ax = axs[0]; ref = col("A", "ref_sig_nv"); cfd = col("A", "cfd_sig_nv")
m = np.abs(ref) > 0.2 * np.max(np.abs(ref)); frac = float(np.median(cfd[m] / ref[m]))
ax.plot(t, ref, color=OI["black"], lw=1.8, label="reference  $\\sigma_{nv}=-2\\mu\\dot b/b$"); ax.plot(t, cfd, color=OI["verm"], lw=1.6, label="CFD, formulation A (tutorial)")
ax.set_title("(a) viscous normal wall stress $\\sigma_{nv}$ (combined, M2)", fontsize=9.5, loc="left"); ax.set_xlabel("t (s)"); ax.set_ylabel("Pa"); ax.grid(True, color="#dddddd", lw=0.6)
ax.text(0.97, 0.04, "CFD transmits ≈ %.2f of the\nreference (0.29–0.32 on M1–M3,\nnot converging)" % frac, transform=ax.transAxes, fontsize=7.8, color=OI["verm"], va="bottom", ha="right")
ax.legend(fontsize=7.3, frameon=False, loc="upper left", bbox_to_anchor=(0.0, 1.0)); ax.set_ylim(-125, 105)
def panel(ax, f, title):
    pr, pc = col(f, "ref_p_b"), col(f, "cfd_p_b"); nr, nc = col(f, "ref_sig_nv"), col(f, "cfd_sig_nv"); rr, rc = col(f, "ref_sig_rr"), col(f, "cfd_sig_rr")
    e_p = -(pc - pr); e_nv = nc - nr; e_rr = rc - rr
    l2 = lambda x: float(np.sqrt(np.mean(x ** 2)))
    ax.plot(t, e_p, color=OI["blue"], lw=1.4, label="pressure term: $-(p_b^{CFD}-p_b^{ref})$")
    ax.plot(t, e_nv, color=OI["verm"], lw=1.4, label="viscous term: $\\sigma_{nv}^{CFD}-\\sigma_{nv}^{ref}$")
    ax.plot(t, e_rr, color=OI["black"], lw=2.0, label="total $\\sigma_{rr}$ error (sum of the two)")
    ax.axhline(0, color="#999999", lw=0.6); ax.set_title(title, fontsize=9.5, loc="left"); ax.set_xlabel("t (s)"); ax.set_ylabel("error of the contribution to $\\sigma_{rr}$ (Pa)"); ax.grid(True, color="#dddddd", lw=0.6)
    ax.text(0.03, 0.03, "L2 errors: pressure %.0f, viscous %.0f,\ntotal %.1f Pa" % (l2(e_p), l2(e_nv), l2(e_rr)), transform=ax.transAxes, fontsize=7.5, ha="left", va="bottom", color=OI["grey"]); return l2(e_rr) / max(l2(e_p), l2(e_nv)), float(np.max(np.abs(rr)))
ra, pk = panel(axs[1], "A", "(b) error budget of $\\sigma_{rr}$, formulation A")
rb, _ = panel(axs[2], "B", "(c) same, formulation B (Gauss $\\nabla p$)")
for ax in axs[1:]: ax.set_ylim(-100, 118)
h, l = axs[1].get_legend_handles_labels(); fig.legend(h, l, loc="lower center", ncol=3, fontsize=7.8, frameon=False, bbox_to_anchor=(0.66, 0.085))
axs[1].text(0.97, 0.955, "compensation: total error / larger part = %.2f" % ra, transform=axs[1].transAxes, fontsize=7.8, ha="right", va="top", color=OI["grey"])
axs[2].text(0.97, 0.955, "compensation: total error / larger part = %.2f" % rb, transform=axs[2].transAxes, fontsize=7.8, ha="right", va="top", color=OI["grey"])
fig.suptitle("Normal wall traction of radialTorsionalRing: the pressure and viscous errors are large but compensate; $\\sigma_{nv}$ and its interaction are not verified", fontsize=10.2, x=0.01, ha="left", y=0.995)
fig.text(0.01, 0.005, "Sign convention: total normal traction σ_rr = −p_b + σ_nv (verified to 1e-10 in every run). Peak |σ_rr| = %.0f Pa, so the total error is ≈ %.1f %% (A) of peak although the viscous term is wrong by a factor ≈ 3. "
         "Matched data: combined member, mesh M2 (128×16×32), Δt = 2.5 ms; A = tutorial (leastSquares ∇p and ∇U, stock wall BC), run F3_M2C; B = A with grad(p) Gauss linear only, run C2_M2C. "
         "Mechanism (source analysis + controlled research run, see README): stock newMovingWallVelocity removes the wall-normal part of the wall velocity gradient. Data: data/fig3_traction_timeseries.csv." % (pk, 100 * float(np.sqrt(np.mean((col("A","cfd_sig_rr")-col("A","ref_sig_rr"))**2))) / pk), fontsize=6.6, color=OI["grey"], va="bottom", wrap=True)
fig.tight_layout(rect=[0, 0.15, 1, 0.955]); fig.savefig(OUT + ".png", dpi=150); print("written", OUT + ".png", "transmitted fraction %.3f; cancellation ratios A %.2f B %.2f" % (frac, ra, rb))
