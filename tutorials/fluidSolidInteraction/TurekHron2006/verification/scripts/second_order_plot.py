#!/usr/bin/env python3
"""Plot Q_in and |Q_in - limit| against the level for the second-order families.

    python3 second_order_plot.py reference/second_order/families.json reference/second_order/q_in_families.png
"""
import json
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

LABELS = {
    "baseline_velLap": "baseline (linearUpwind, velocityLaplacian; earlier windows)",
    "LU_dL": "linearUpwind + dispLap",
    "central_dL": "central + dispLap",
    "central_dL_pwall": "central + dispLap + pwall",
    "rb_central_dL_pwall": "ring+box mesh family, central + dispLap + pwall",
}


def main() -> int:
    data = json.load(open(sys.argv[1]))["families"]
    fig, (a, b) = plt.subplots(1, 2, figsize=(11, 4.2))
    h = [1.0, 0.5, 0.25]
    for name, f in data.items():
        q = f["qoi"]["Q_in"]
        v = [q["v1"], q["v2"], q["v3"]]
        a.plot([1, 2, 4], v, "o-", label=LABELS.get(name, name))
        d = [abs(q["d21"]), abs(q["d32"])]
        b.loglog([0.75, 0.375], d, "o-", label=LABELS.get(name, name))
    b.loglog([0.75, 0.375], [30, 7.5], "k--", lw=0.8, label="slope 2")
    b.loglog([0.75, 0.375], [30, 15], "k:", lw=0.8, label="slope 1")
    a.set_xscale("log", base=2)
    a.set_xticks([1, 2, 4], ["1x", "2x", "4x"])
    a.set_xlabel("level")
    a.set_ylabel("Q_in (N/m)")
    a.grid(alpha=0.3)
    b.set_xlabel("mean relative cell size of the pair")
    b.set_ylabel("|change of Q_in| between levels (N/m)")
    b.grid(alpha=0.3, which="both")
    b.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(sys.argv[2], dpi=130)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
