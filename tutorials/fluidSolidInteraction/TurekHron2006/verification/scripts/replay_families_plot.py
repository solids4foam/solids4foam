#!/usr/bin/env python3
"""Q_in and plate lift amplitude of the two replay mesh families against cell count."""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# common window ending 30.2 s (reference/second_order/families.csv, reference/production/replay_8x_runs.csv)
STD = {"cells": [5336, 21344, 48024, 85376], "Q_in": [43.66, 76.17, 81.22, 81.96],
       "lift": [124.26, 142.38, 147.71, 150.49], "lab": ["1x", "2x", "3x", "4x"]}
RB = {"cells": [11764, 47056, 105876, 188224, 752896], "Q_in": [67.18, 79.63, 79.99, 79.27, 77.16],
      "lift": [155.30, 155.59, 156.02, 157.22, 160.85], "lab": ["1x'", "2x'", "3x'", "4x'", "8x'"]}
C1, C2, INK, MUTED = "#2a78d6", "#eb6834", "#1a1a19", "#6b6a63"

fig, axes = plt.subplots(1, 2, figsize=(10, 4))
for ax, key, ylab in ((axes[0], "Q_in", "Q_in (N/m)"), (axes[1], "lift", "plate lift fundamental (N/m)")):
    for fam, col, name in ((STD, C1, "standard"), (RB, C2, "ring+box")):
        ax.plot(fam["cells"], fam[key], "-o", color=col, lw=2, ms=6, mec="white", mew=1.5, label=name)
        ax.annotate(fam["lab"][-1], (fam["cells"][-1], fam[key][-1]), textcoords="offset points",
                    xytext=(6, -3), color=INK, fontsize=9)
    ax.set_xscale("log")
    ax.set_xlabel("fluid cells", color=MUTED)
    ax.set_ylabel(ylab, color=MUTED)
    ax.grid(True, color="#e6e5e0", lw=0.8)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.tick_params(colors=MUTED)
axes[0].set_ylim(40, 86)
axes[0].legend(frameon=False)
fig.suptitle("Replay of the coupled 2x motion: central, displacementLaplacian, movingWallPressure",
             fontsize=10, color=INK)
fig.tight_layout()
out = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("replay_families.png")
fig.savefig(out, dpi=150)
print(out)
