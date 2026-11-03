#!/usr/bin/env python3
"""Figure 1 (schematic, vector): geometry and loading of radialTorsionalRing.
Left: SCHEMATIC of the annulus (drawn to scale in radius: A = 0.5, B = 1, c = 2 m; the deformation of the interface is exaggerated).
Right: the two drives, plotted from the tutorial's own tables constant/timeVsAngle (alpha_A, tabulated in degrees, converted to rad) and constant/pcTable (p_c/rho in m^2 s^-2, with rho = 1000 kg m^-3 numerically equal to p_c in kPa).
Usage: python3 fig1_geometry.py [TUTORIAL_DIR] [OUT_PREFIX]   (defaults: .. and fig1_geometry)"""
import os, re, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams.update({"svg.fonttype": "path", "font.size": 9, "font.family": "DejaVu Sans"})
import matplotlib.pyplot as plt
from matplotlib.patches import Wedge, Circle, FancyArrowPatch, Arc
TUT = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
OUT = sys.argv[2] if len(sys.argv) > 2 else "fig1_geometry"
OI = dict(blue="#0072B2", orange="#E69F00", green="#009E73", verm="#D55E00", purple="#CC79A7", sky="#56B4E9", grey="#555555")
def table(name):
    a = np.array([[float(x) for x in m] for m in re.findall(r"\(\s*([-\d.eE+]+)\s+([-\d.eE+]+)\s*\)", open(os.path.join(TUT, "constant", name)).read())])
    return a[:, 0], a[:, 1]
t, ang = table("timeVsAngle"); tp, pc = table("pcTable")
ang = np.radians(ang)      # constant/timeVsAngle is tabulated in degrees (fixedRotation); peak 5.7296 deg = 0.1 rad
A, B, c = 0.5, 1.0, 2.0
fig = plt.figure(figsize=(11.0, 4.8))
ax = fig.add_axes([0.0, 0.02, 0.56, 0.96]); ax.set_aspect("equal"); ax.axis("off")
ax.add_patch(Wedge((0, 0), c, 0, 360, width=c - B, facecolor="#cfe5f3", edgecolor="none"))
ax.add_patch(Wedge((0, 0), B, 0, 360, width=B - A, facecolor="#e8d9b5", edgecolor="none", hatch="////", linewidth=0))
for r, lw, ls in ((A, 1.4, "-"), (B, 2.2, "-"), (c, 1.4, "-")):
    ax.add_patch(Circle((0, 0), r, fill=False, lw=lw, ls=ls, edgecolor="k"))
th = np.linspace(0, 2 * np.pi, 400); bd = B * (1 + 0.06 * np.cos(2 * th) * 0)           # deformed interface, schematic (breathing: uniform)
ax.add_patch(Circle((0, 0), B * 1.07, fill=False, lw=1.2, ls="--", edgecolor=OI["verm"]))
ax.text(0.0, -0.74, "solid ring", ha="center", va="center", fontsize=8, color="#5a4a20")
ax.text(0.0, -1.5, "viscous fluid (Newtonian)\nopen outer boundary", ha="center", va="center", fontsize=8, color="#1f4e79")
# radial dimension line along +x: A, B, c
ax.annotate("", xy=(c, 0.0), xytext=(0, 0.0), arrowprops=dict(arrowstyle="<->", lw=0.8, color="k"))
for r, lab in ((A, "A=0.5"), (B, "B=1"), (c, "c=2 m")):
    ax.plot([r, r], [-0.05, 0.05], color="k", lw=0.8); ax.text(r, -0.14, lab, ha="center", va="top", fontsize=8, bbox=dict(fc="white", ec="none", pad=0.6, alpha=0.8))
# imposed inner torsion (arc at r = 0.36 inside the ring hole)
ax.add_patch(Arc((0, 0), 2 * 0.34, 2 * 0.34, theta1=20, theta2=140, lw=2, color=OI["blue"]))
ax.add_patch(FancyArrowPatch((0.34 * np.cos(np.radians(138)), 0.34 * np.sin(np.radians(138))), (0.34 * np.cos(np.radians(150)), 0.34 * np.sin(np.radians(150))), arrowstyle="-|>", mutation_scale=12, color=OI["blue"], lw=2))
ax.text(0.0, 0.47, r"$\alpha_A(t)$", color=OI["blue"], fontsize=10, ha="center")
ax.annotate("imposed inner-wall torsion (r = A: clamped radially, rotated)\n" + r"$\alpha_A=\theta_0 S(t)\sin\Omega t$, $\theta_0=0.1$ rad, $\Omega=3$ rad/s", xy=(-0.30, 0.19), xytext=(-3.45, 2.45), fontsize=8, color=OI["blue"], arrowprops=dict(arrowstyle="-", color=OI["blue"], lw=0.8), ha="left", va="center")
# outer pressure (negative: directed outward)
for a in np.radians(np.arange(15, 360, 30)):
    ax.add_patch(FancyArrowPatch((c * np.cos(a), c * np.sin(a)), ((c + 0.22) * np.cos(a), (c + 0.22) * np.sin(a)), arrowstyle="-|>", mutation_scale=9, color=OI["green"], lw=1.5))
ax.text(0, -2.62, r"outer boundary $p_c(t)=P_0S(t)$, $P_0=-3$ kPa (negative: suction, arrows outward)", ha="center", fontsize=8, color=OI["green"])
# interface label and directions
ax.annotate("fluid/solid interface\nradius $b(t)=B+u_B(t)$\n(dashed: deformed, exaggerated)", xy=(B * 1.07 * np.cos(np.radians(-25)), B * 1.07 * np.sin(np.radians(-25))), xytext=(2.45, -1.0), fontsize=8, color=OI["verm"], arrowprops=dict(arrowstyle="-", color=OI["verm"], lw=0.8), ha="left", va="center")
a0 = np.radians(60); P0 = ((B + 0.05) * np.cos(a0), (B + 0.05) * np.sin(a0))
er = np.array([np.cos(a0), np.sin(a0)]); et = np.array([-np.sin(a0), np.cos(a0)])
ax.add_patch(FancyArrowPatch(P0, tuple(np.array(P0) + 0.38 * er), arrowstyle="-|>", mutation_scale=11, color="k", lw=1.4))
ax.add_patch(FancyArrowPatch(P0, tuple(np.array(P0) + 0.38 * et), arrowstyle="-|>", mutation_scale=11, color="k", lw=1.4))
ax.text(*(np.array(P0) + 0.52 * er), r"$\mathbf{e}_r$", fontsize=9, ha="center", va="center"); ax.text(*(np.array(P0) + 0.55 * et), r"$\mathbf{e}_\theta$", fontsize=9, ha="center", va="center")
ax.text(4.5, 2.55, "SCHEMATIC\n(radii to scale)", fontsize=8, fontweight="bold", color=OI["grey"], ha="right", va="top")
ax.set_xlim(-3.5, 4.6); ax.set_ylim(-2.8, 2.7)
# drives
a1 = fig.add_axes([0.64, 0.58, 0.33, 0.33]); a1.plot(t, ang, color=OI["blue"], lw=1.6); a1.axhline(0, color="#999999", lw=0.6)
a1.set_ylabel(r"$\alpha_A$ (rad)"); a1.set_xlabel("t (s)"); a1.set_title("inner-wall rotation drive (torsion member and combined)", fontsize=8.5, loc="left")
a1.axvspan(0, 1, color="#eeeeee", zorder=0); a1.text(0.5, -0.095, r"ramp $S(t)$ active", ha="center", va="bottom", fontsize=8, color=OI["grey"])
a2 = fig.add_axes([0.64, 0.12, 0.33, 0.33]); a2.plot(tp, pc, color=OI["green"], lw=1.6); a2.axhline(0, color="#999999", lw=0.6)
a2.set_ylabel(r"$p_c$ (kPa)"); a2.set_xlabel("t (s)"); a2.set_title(r"outer pressure drive (radial and combined); table stores $p_c/\rho$", fontsize=8.5, loc="left")
for a in (a1, a2):
    for s in ("top", "right"): a.spines[s].set_visible(False)
fig.savefig(OUT + ".svg"); fig.savefig(OUT + ".png", dpi=170)
with open(OUT + ".svg") as fh: _t = [ln.rstrip() for ln in fh.read().split("\n")]       # no trailing whitespace (repository whitespace check)
with open(OUT + ".svg", "w") as fh: fh.write("\n".join(_t))
print("written", OUT, "alpha_A max %.4f rad, p_c/rho min %.3f" % (ang.max(), pc.min()))
