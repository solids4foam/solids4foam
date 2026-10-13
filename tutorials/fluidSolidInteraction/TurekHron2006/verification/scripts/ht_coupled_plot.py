#!/usr/bin/env python3
"""Coupled FSI3 production sequence: QoIs against level, with the time-step runs.

    python3 ht_coupled_plot.py coupled.json out.png
"""
import json
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

runs = {r["name"]: r for r in json.load(open(sys.argv[1]))}
PANELS = [("uy_amplitude", 1e3, "u_y amplitude (mm)"), ("ux_midrange", 1e3, "u_x (max+min)/2 (mm)"),
          ("ux_timemean", 1e3, "u_x time mean (mm)"), ("frequency", 1, "frequency (Hz)"),
          ("drag_timemean", 1, "drag mean (N/m)"), ("lift_amplitude_smoothed", 1, "lift amplitude (N/m)")]
C1, C2, C3, MUTED, INK = "#2a78d6", "#eb6834", "#1baf7a", "#6b6a63", "#1a1a19"
fig, axes = plt.subplots(2, 3, figsize=(12, 6.5))
for ax, (q, s, lab) in zip(axes.flat, PANELS):
    x = [1, 2, 4]
    y = [runs[f"{n}x_dt5e-4"][q] * s for n in (1, 2, 4)]
    ax.plot(x, y, "-o", color=C1, lw=2, ms=7, mec="white", mew=1.5, label="dt 5e-4 (spatial)")
    for n, off in ((1, -0.06), (2, 0.06)):
        for dt, col, mk in (("1e-3", C2, "s"), ("2.5e-4", C3, "^")):
            k = f"{n}x_dt{dt}"
            if k in runs:
                ax.plot(n * (1 + off), runs[k][q] * s, mk, color=col, ms=7, mec="white", mew=1,
                        label=f"dt {dt}" if (n == 2 and q == PANELS[0][0]) else None)
    ax.set_xscale("log", base=2)
    ax.set_xticks(x, ["1x'", "2x'", "4x'"])
    ax.set_title(lab, fontsize=10, color=INK)
    ax.grid(True, color="#e6e5e0", lw=0.8)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    ax.tick_params(colors=MUTED)
handles, labels = [], []
for a in axes.flat:
    h, l = a.get_legend_handles_labels()
    for hh, ll in zip(h, l):
        if ll not in labels:
            handles.append(hh)
            labels.append(ll)
fig.legend(handles, labels, loc="lower center", ncol=3, frameon=False)
fig.suptitle("Coupled FSI3, ring+box fluid with solid 2x/4x/8x, five periods ending at 6.0 s", fontsize=11, color=INK)
fig.tight_layout(rect=(0, 0.06, 1, 0.96))
fig.savefig(sys.argv[2], dpi=150)
print(sys.argv[2])
