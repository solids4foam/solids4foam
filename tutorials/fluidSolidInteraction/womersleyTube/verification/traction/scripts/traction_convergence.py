#!/usr/bin/env python3
"""Convergence analysis and plots of the interface-traction diagnostic.

usage: traction_convergence.py OUT_DIR STUDY_DIR[:COUPLING] [...]   (COUPLING robin|iqnils, default robin)
Each STUDY_DIR must already contain verification/postProcessing/
interfaceTraction_summary.json (from womersley_interface_traction.py).

For each study and each complex error measure E (relative to the exact
value), it reports per case |E| and, for the mesh series m1/m2/m4 at 200
steps per period and the time series 100/200/400 on m2:
  - observed order from successive differences of the complex values
    p = log(|E1 - E2| / |E2 - E4|) / log 2   (reference-free),
  - Richardson limit E_inf = E4 + (E4 - E2)/(2^p_used - 1) with p_used = 2
    (the formal order) and with the observed p, as the estimate of the
    continuum (model-form) part of the error; their difference is quoted as
    the uncertainty of that estimate.
Writes OUT_DIR/traction_convergence.json, .csv, a Markdown table and plots.
"""
import csv
import json
import math
import os
import sys
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", os.environ.get("TMPDIR", "/tmp"))
try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
except ImportError:      # plots are optional
    plt = None

MEASURES = ["P_err/|p|", "Vn_ratio", "Vn_err/|p|", "Sn_err/|p|", "Vx_err/|Vx|",
            "VnEst_ratio", "mom_vs_transfer_n/|Vn|", "mom_vs_transfer_x/|Vx|",
            "solid_vs_fluid_n/|p|", "solid_vs_fluid_x/|Vx|"]


def c(v):
    return complex(v[0], v[1])


def series_analysis(vals):
    """vals: complex values on levels (coarse, mid, fine), ratio 2."""
    e1, e2, e4 = vals
    d1, d2 = abs(e1 - e2), abs(e2 - e4)
    p = math.log(d1 / d2) / math.log(2) if d1 > 0 and d2 > 0 else float("nan")
    lim2 = e4 + (e4 - e2) / 3.0
    limp = e4 + (e4 - e2) / (2 ** p - 1) if math.isfinite(p) and p > 0.3 else lim2
    return p, lim2, limp


def main():
    out = Path(sys.argv[1])
    out.mkdir(parents=True, exist_ok=True)
    report = {}
    md = ["# Interface-traction diagnostic: convergence", ""]
    rows = []
    for arg in sys.argv[2:]:
        sd, _, cp = arg.partition(":")
        cp = cp or "robin"
        sd = Path(sd)
        label = (sd.parent.name if sd.name == "womersleyTube" else sd.name) + f"_{cp}"
        summ = json.loads((sd / "verification/postProcessing/interfaceTraction_summary.json").read_text())
        cases = summ["cases"]
        rep = {"exact": summ["exact"], "cases": {}, "mesh": {}, "time": {}, "time_m4": {}}
        for name, cdata in cases.items():
            rep["cases"][name] = {m: cdata["errors"][m] for m in MEASURES if m in cdata["errors"]}
            row = {"study": label, "case": name}
            for m in MEASURES:
                if m in cdata["errors"]:
                    v = c(cdata["errors"][m])
                    row[m + " re"], row[m + " im"], row[m + " abs"] = v.real, v.imag, abs(v)
            rows.append(row)
        for kind, names in (("mesh", [f"{cp}_m1_n200", f"{cp}_m2_n200", f"{cp}_m4_n200"]),
                            ("time", [f"{cp}_m2_n100", f"{cp}_m2_n200", f"{cp}_m2_n400"]),
                            ("time_m4", [f"{cp}_m4_n200", f"{cp}_m4_n400"])):
            if not all(n in cases for n in names):
                continue
            for m in MEASURES:
                vals = [c(cases[n]["errors"][m]) for n in names]
                if len(vals) == 2:       # time pair on m4: q = 2 limit only
                    lim2 = vals[1] + (vals[1] - vals[0]) / 3.0
                    rep[kind][m] = {"values": [[v.real, v.imag] for v in vals], "order": float("nan"),
                                    "limit_p2": [lim2.real, lim2.imag], "limit_pobs": [lim2.real, lim2.imag],
                                    "limit_uncertainty": float("nan")}
                    continue
                p, lim2, limp = series_analysis(vals)
                rep[kind][m] = {"values": [[v.real, v.imag] for v in vals], "order": p,
                                "limit_p2": [lim2.real, lim2.imag],
                                "limit_pobs": [limp.real, limp.imag],
                                "limit_uncertainty": abs(lim2 - limp)}
        report[label] = rep
        md += [f"## {label}", "", "| measure | series | finest value | observed order | Richardson limit (p=2) | |limit p=2 - p_obs| |",
               "|---|---|---|---|---|---|"]
        for kind in ("mesh", "time", "time_m4"):
            for m, r in rep[kind].items():
                fin = complex(*r["values"][-1])
                lim = complex(*r["limit_p2"])
                md.append(f"| {m} | {kind} | {fin.real:+.3e}{fin.imag:+.3e}i | {r['order']:.2f} | "
                          f"{lim.real:+.3e}{lim.imag:+.3e}i | {r['limit_uncertainty']:.1e} |")
        md.append("")
        # plots: |error| vs h for the physical measures
        if rep["mesh"] and plt is not None:
            fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
            h = [1, 0.5, 0.25]
            for m in ("P_err/|p|", "Sn_err/|p|", "Vn_err/|p|", "Vx_err/|Vx|"):
                vals = [abs(complex(*v)) for v in rep["mesh"][m]["values"]]
                ax[0].loglog(h, vals, "o-", label=m)
            ax[0].loglog(h, [1e-3 * x**2 for x in h], "k--", lw=0.8, label="slope 2")
            ax[0].loglog(h, [1e-3 * x for x in h], "k:", lw=0.8, label="slope 1")
            ax[0].set_xlabel("relative cell size h/h(m1)")
            ax[0].set_ylabel("|relative error| of wave-projected first harmonic")
            ax[0].set_title(f"{label}: mesh series (200 steps/period)")
            ax[0].legend(fontsize=7)
            for m, mk in (("Vn_ratio", "o-"), ("VnEst_ratio", "s--")):
                vals = [complex(*v) for v in rep["mesh"][m]["values"]]
                ax[1].plot(h, [v.real for v in vals], mk, label=m + " (Re)")
                ax[1].plot(h, [v.imag for v in vals], mk, alpha=0.5, label=m + " (Im)")
            ax[1].axhline(1 / 3, color="gray", lw=0.8, ls=":")
            ax[1].axhline(1.0, color="gray", lw=0.8)
            ax[1].set_xscale("log")
            ax[1].set_xlabel("relative cell size h/h(m1)")
            ax[1].set_title("normal viscous stress: transferred / exact")
            ax[1].legend(fontsize=7)
            fig.tight_layout()
            fig.savefig(out / f"traction_mesh_{label}.png", dpi=130)
            plt.close(fig)
    (out / "traction_convergence.json").write_text(json.dumps(report, indent=1))
    keys = sorted({k for r in rows for k in r})
    keys = ["study", "case"] + [k for k in keys if k not in ("study", "case")]
    with (out / "traction_cases_all.csv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)
    (out / "traction_convergence.md").write_text("\n".join(md) + "\n")
    print("\n".join(md))


if __name__ == "__main__":
    main()
