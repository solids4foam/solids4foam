#!/usr/bin/env python3
"""Plot the Richter campaign evidence from reference/campaign/*.csv.

Writes PNG figures next to the data:
- coupling_iterations.png: FSI iterations per step, PIMPLE 1 vs 3 (L1) and
  the coupling-tolerance sweep (L1),
- solid_static.png: static solid u_x(A) against cells through the thickness,
- convergence.png: production-sequence u_x(A) and F_x, with rigid-beam F_x,
  against the refinement factor.
"""

from __future__ import annotations

import csv
import gzip
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

DATA = Path(__file__).resolve().parents[1] / "reference" / "campaign"
COLOURS = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#4a3aa7"]
INK, MUTED, GRID = "#0b0b0b", "#52514e", "#e4e3df"


def style(ax):
    ax.grid(True, color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(MUTED)
    ax.tick_params(colors=MUTED)
    ax.xaxis.label.set_color(INK)
    ax.yaxis.label.set_color(INK)


def read_steps(path: Path):
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as handle:
        rows = list(csv.DictReader(handle))
    return [float(r["time"]) for r in rows], [int(r["iterations"]) for r in rows]


def binned(times, values, width=0.25):
    out_t, out_v, bucket, edge = [], [], [], width
    for t, v in zip(times, values):
        while t > edge:
            if bucket:
                out_t.append(edge - width / 2)
                out_v.append(sum(bucket) / len(bucket))
            bucket, edge = [], edge + width
        bucket.append(v)
    return out_t, out_v


def coupling_plot():
    fig, axes = plt.subplots(1, 2, figsize=(10, 3.6))
    diag = DATA / "coupling_diagnostic_L1m4"
    for colour, (name, label) in zip(
        COLOURS,
        [("L1m4_p1", "PIMPLE 1 corrector"), ("L1m4_p3", "PIMPLE 3 correctors")],
    ):
        t, n = read_steps(diag / f"{name}_steps.csv")
        bt, bn = binned(t, n)
        axes[0].plot(bt, bn, color=colour, linewidth=2, label=label)
    axes[0].set_xlabel("time (s)")
    axes[0].set_ylabel("FSI iterations per step (0.25 s mean)")
    axes[0].set_title("L1, strict test 1e-6", color=INK, fontsize=10)
    axes[0].legend(frameon=False)
    matrix = DATA / "L1_coupling_matrix"
    for colour, tol in zip(COLOURS, ["1e-4", "1e-5", "1e-6", "1e-7"]):
        path = matrix / f"L1_tol{tol}_steps.csv.gz"
        if not path.exists():
            continue
        t, n = read_steps(path)
        bt, bn = binned(t, n, 0.5)
        axes[1].plot(bt, bn, color=colour, linewidth=2, label=f"strict {tol}")
    axes[1].set_xlabel("time (s)")
    axes[1].set_ylabel("FSI iterations per step (0.5 s mean)")
    axes[1].set_title("L1, PIMPLE 3, tolerance sweep (cap 50)", color=INK, fontsize=10)
    axes[1].legend(frameon=False)
    for ax in axes:
        style(ax)
    fig.tight_layout()
    fig.savefig(DATA / "coupling_iterations.png", dpi=150)


def solid_plot():
    with open(DATA / "solid_static.csv") as handle:
        rows = list(csv.DictReader(handle))
    series = {}
    for r in rows:
        key = (r["family"], r["discretisation"])
        series.setdefault(key, []).append(r)
    labels = {
        ("m1", "standard"): "standard, N_t = 4f (original family)",
        ("m2", "standard"): "standard, N_t = 8f (production family)",
        ("thickness x4 at fixed in-plane factor", "standard"): "standard, N_t = 16f",
        ("in-plane x2", "standard"): "standard, N_t = 4f, in-plane x2",
        ("m1", "high-order"): "high-order, N_t = 4f",
    }
    fig, ax = plt.subplots(figsize=(6, 4))
    for colour, (key, label) in zip(COLOURS, labels.items()):
        data = sorted(series.get(key, []), key=lambda r: int(r["cells_through_thickness"]))
        if not data:
            continue
        ax.plot(
            [int(r["cells_through_thickness"]) for r in data],
            [1e4 * float(r["ux_A"]) for r in data],
            marker="o", markersize=5, color=colour, linewidth=2, label=label,
        )
    ax.set_xscale("log", base=2)
    ax.set_xlabel("cells through the beam thickness N_t")
    ax.set_ylabel("u_x(A) (1e-4 m), 100 Pa on upstream face")
    ax.legend(frameon=False, fontsize=8)
    style(ax)
    fig.tight_layout()
    fig.savefig(DATA / "solid_static.png", dpi=150)


def convergence_plot():
    production = DATA / "production_sequence.csv"
    if not production.exists():
        return
    with open(production) as handle:
        prod = list(csv.DictReader(handle))
    with open(DATA / "rigid_fluid.csv") as handle:
        rigid = [r for r in csv.DictReader(handle) if abs(float(r["dt"]) - 0.00625) < 1e-9]
    fig, axes = plt.subplots(1, 2, figsize=(10, 3.8))
    f = [int(r["factor"]) for r in prod]
    axes[0].plot(f, [1e5 * float(r["ux"]) for r in prod], marker="o", color=COLOURS[0],
                 linewidth=2, markersize=6, label="FSI u_x(A)")
    axes[0].set_ylabel("u_x(A) (1e-5 m)")
    axes[1].plot(f, [float(r["fx"]) for r in prod], marker="o", color=COLOURS[0],
                 linewidth=2, markersize=6, label="FSI F_x")
    axes[1].plot([int(r["factor"]) for r in rigid], [float(r["fx"]) for r in rigid],
                 marker="s", color=COLOURS[1], linewidth=2, markersize=6,
                 label="rigid beam F_x")
    axes[1].set_ylabel("F_x (N)")
    for ax in axes:
        ax.set_xscale("log", base=2)
        ax.set_xlabel("refinement factor f (L0 = 1)")
        ax.legend(frameon=False)
        style(ax)
    fig.tight_layout()
    fig.savefig(DATA / "convergence.png", dpi=150)


if __name__ == "__main__":
    coupling_plot()
    solid_plot()
    convergence_plot()
