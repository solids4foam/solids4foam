#!/usr/bin/env python3
"""
Compare the solids4foam displacement histories at p0 and p1 with the
solutions submitted to the Arostica et al. (2025) benchmark.

The reference curves (verification/reference/<step>/<group>.csv) were
converted from the pickles distributed with the cardiac-benchmark-toolkit
(https://github.com/Reidmen/cardiac_benchmark_toolkit, MIT licence).

Usage (from the case directory, after ./Allrun <step>):
    python3 scripts/compareWithReference.py step1     # benchmark 1
    python3 scripts/compareWithReference.py coarse    # benchmark 2
writes verification/<step>.png and prints, per component, the largest
difference from the median of the reference solutions over t in [0, 1] s.
"""
import glob
import os
import sys

import numpy as np


def read_s4f(name):
    # Columns: time Dx Dy Dz |D|, from the solidPointDisplacement function object
    files = sorted(glob.glob(f"postProcessing/*/solidPointDisplacement_{name}.dat"))
    if not files:
        sys.exit(f"No postProcessing/*/solidPointDisplacement_{name}.dat found")
    data = np.vstack([np.loadtxt(f, comments="#") for f in files])
    _, idx = np.unique(data[:, 0], return_index=True)
    return data[idx]


def main():
    steps = sorted(os.path.basename(d) for d in glob.glob("verification/reference/*"))
    step = sys.argv[1] if len(sys.argv) > 1 else steps[0]
    if step not in steps:
        sys.exit(f"Unknown reference set {step}; choose from {steps}")

    files = sorted(glob.glob(f"verification/reference/{step}/*.csv"))
    header = open(files[0]).readline().strip().split(",")[1:]
    points = sorted({h.split("_")[0] for h in header})
    refs = {os.path.splitext(os.path.basename(f))[0]:
            np.loadtxt(f, delimiter=",", skiprows=1) for f in files}
    nc = len(header)

    ours, t_ours, labels = [], [], []
    for p in points:
        d = read_s4f(p)
        for c, comp in enumerate(["ux", "uy", "uz"]):
            ours.append(d[:, 1 + c])
            t_ours.append(d[:, 0])
            labels.append(f"{p} {comp}")

    t_ref = next(iter(refs.values()))[:, 0]
    stack = np.stack([[np.interp(t_ref, r[:, 0], r[:, 1 + c]) for c in range(nc)]
                      for r in refs.values()])
    median = np.median(stack, axis=0)
    spread = stack.max(axis=0) - stack.min(axis=0)

    print(f"{step}: solids4foam vs median of {len(refs)} reference solutions")
    print(f"{'component':10s} {'max|diff| [mm]':>15s} {'ref spread [mm]':>16s}")
    for k in range(nc):
        mask = t_ref <= t_ours[k][-1] + 1e-12
        o = np.interp(t_ref[mask], t_ours[k], ours[k])
        diff = np.abs(o - median[k, mask]).max() * 1e3
        print(f"{labels[k]:10s} {diff:15.3f} {spread[k].max()*1e3:16.3f}")

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib not available: no plot written")
        return

    nP = len(points)
    fig, axes = plt.subplots(3, nP, figsize=(5.5 * nP, 10), sharex=True,
                             squeeze=False)
    for k in range(nc):
        ax = axes[k % 3, k // 3]
        for g, r in refs.items():
            ax.plot(r[:, 0], r[:, 1 + k] * 1e3, color="0.75", lw=1)
        ax.plot(t_ours[k], ours[k] * 1e3, color="C3", lw=2, label="solids4foam")
        ax.set_title(labels[k])
        ax.set_ylabel("displacement [mm]")
        ax.grid(alpha=0.3)
    axes[0, 0].plot([], [], color="0.75", label="benchmark participants")
    axes[0, 0].legend(loc="best")
    for ax in axes[-1]:
        ax.set_xlabel("time [s]")
    fig.suptitle(f"Arostica et al. (2025), {os.path.basename(os.getcwd())}, {step}")
    fig.tight_layout()
    out = f"verification/{step}.png"
    fig.savefig(out, dpi=120)
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
