#!/usr/bin/env python3
"""Time order of the two fluxConsistentPatches implementations on the
fluid-only static womersleyTube problem (section 15 of
womersley_temporal_investigation.md).

For runs iqnils_<mesh>_n<N>_p6_fluidOnly_transpiration_<impl>[_<extra>]
already analysed into verification/temporal/runs, prints the
successive-difference orders of the complex flow-rate coefficient and of
selected QoIs (periods 5-6) for the exact-end-flux form (fluxConsistent) and
the final form (fluxConsistent2), and D = fix2 - fix1, the effect of the
remaining boundary-flux term, with its halving ratios.

    womersley_boundary_order.py m1 200,400,...,12800 [extra]
"""

import cmath
import json
import math
import sys
from pathlib import Path

RUNS = Path(__file__).resolve().parents[1] / "temporal" / "runs"


def window(mesh, steps, impl, extra):
    name = f"iqnils_{mesh}_n{steps}_p6_fluidOnly_transpiration_{impl}"
    name += f"_{extra}" if extra else ""
    return json.loads((RUNS / f"{name}.json").read_text())["windows"]["5-6"]


def coefficient(w, q):
    return (1 + w[f"{q}_amp"]) * cmath.exp(1j * w[f"{q}_phase"])


def orders(values):
    d = [values[i] - values[i + 1] for i in range(len(values) - 1)]
    return d, [math.log2(abs(d[i]) / abs(d[i + 1])) for i in range(len(d) - 1)]


def main():
    mesh, steps = sys.argv[1], [int(n) for n in sys.argv[2].split(",")]
    extra = sys.argv[3] if len(sys.argv) > 3 else ""
    runs = {}
    for impl, label in (("fluxConsistent", "exact end flux"),
                        ("fluxConsistent2", "HbyA_b = U_b + rAtU grad(p)_P")):
        ws = [window(mesh, n, impl, extra) for n in steps]
        runs[impl] = ws
        print(f"== {mesh} {impl} ({label}), steps {steps}")
        d, o = orders([coefficient(w, "flow") for w in ws])
        print("  complex flow coefficient: |d| "
              + " ".join(f"{abs(x):.2e}" for x in d)
              + " | p " + " ".join(f"{x:.2f}" for x in o))
        for q in ("flow_amp", "flow_phase", "speed_ux", "attenuation",
                  "profile"):
            d, o = orders([w[q] for w in ws])
            print(f"  {q:12s} d " + " ".join(f"{x:+.1e}" for x in d)
                  + " | p " + " ".join(f"{x:.2f}" for x in o))
    D = [coefficient(b, "flow") - coefficient(a, "flow")
         for a, b in zip(runs["fluxConsistent"], runs["fluxConsistent2"])]
    print("D = fix2 - fix1 (complex flow): "
          + " ".join(f"{abs(x):.2e}" for x in D) + " | halving ratios "
          + " ".join(f"{abs(D[i]) / abs(D[i + 1]):.2f}"
                     for i in range(len(D) - 1)))


if __name__ == "__main__":
    main()
