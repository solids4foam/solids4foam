#!/usr/bin/env python3
"""Energy-balance test of the FSI3 amplitude change between meshes.

Reads the per-case output of energy_analysis.py (--json) and the coupled FSI3
results (fsi3_iqnils_refinement_study_levels.csv), and evaluates:

  dW/dA, dW/df   at the 2x mesh (dt 5e-4) from the amplitude-scaled (x0.95,
                 x1.05) and frequency-scaled (x0.99) replays
  dW_m           W_m - W_2x at the same prescribed (2x coupled) motion
  dA_pred        -(dW_m + dW/df * df_m) / (dW/dA), with df_m the coupled
                 frequency change, against the observed coupled dA_m
  dA/dW implied  dA_obs / dW_m
  dA_pred_inphase  alternative in-phase (added-mass/stiffness) balance at fixed
                 frequency: -dQ_in,m / (dQ_in/dA) (fluid only), and
                 dQ_in,m / (S_A - dQ_in/dA) with a linear structure S_A = Q_in/A

for both definitions of the discrete work (power-based W and the
displacement-increment W_modal_sum).

    python3 energy_summary.py energy_runs.json levels.csv --out energy_balance
"""
import argparse
import csv
import json
from pathlib import Path


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("runs", type=Path)
    ap.add_argument("levels", type=Path)
    ap.add_argument("--out", type=Path, default=Path("energy_balance"))
    a = ap.parse_args()
    runs = {r["case"]: r for r in json.loads(a.runs.read_text())}
    coupled = {}
    for r in csv.DictReader(a.levels.open()):
        coupled[r["case"]] = {"A": float(r["uy_amplitude"]), "f": float(r["uy_frequency"]),
                              "ux_mean": float(r["ux_mean"]), "dt": float(r["delta_t"])}
    c2 = coupled["iqnils_mesh_2x"]
    base = runs["e_2x_dt0.0005"]
    out = {"coupled": coupled, "runs": runs, "definitions": {}}
    rows = []
    for wkey in ("W", "W_modal_sum"):
        d = {}
        if "e_2x_dt0.0005_A1.05" in runs and "e_2x_dt0.0005_A0.95" in runs:
            hi, lo = runs["e_2x_dt0.0005_A1.05"], runs["e_2x_dt0.0005_A0.95"]
            d["dW_dA"] = (hi[wkey] - lo[wkey]) / (hi["A_tip_m"] - lo["A_tip_m"])
            d["dW_dA_upper"] = (hi[wkey] - base[wkey]) / (hi["A_tip_m"] - base["A_tip_m"])
            d["dW_dA_lower"] = (base[wkey] - lo[wkey]) / (base["A_tip_m"] - lo["A_tip_m"])
            d["dQin_dA"] = (hi["Q_in"] - lo["Q_in"]) / (hi["A_tip_m"] - lo["A_tip_m"])
            d["S_A_linear_structure"] = base["Q_in"] / base["A_tip_m"]
        if "e_2x_dt0.0005_F0.99" in runs:
            fr = runs["e_2x_dt0.0005_F0.99"]
            d["dW_df"] = (fr[wkey] - base[wkey]) / (fr["frequency_hz"] - base["frequency_hz"])
        paths = {"matched": [("e_1x_dt0.001", "iqnils_mesh_1x"), ("e_2x_dt0.0005", "iqnils_mesh_2x"),
                             ("e_4x_dt0.00025", "iqnils_mesh_4x")],
                 "fixed_dt_2.5e-4": [("e_1x_dt0.00025", "iqnils_mesh_1x"),
                                     ("e_2x_dt0.00025", "iqnils_mesh_2x"),
                                     ("e_4x_dt0.00025", "iqnils_mesh_4x")]}
        for pname, spec in paths.items():
            ref = runs.get(spec[1][0])
            if ref is None:
                continue
            for case, cname in spec:
                if case not in runs or cname not in coupled:
                    continue
                r = runs[case]
                dW = r[wkey] - ref[wkey]
                dA = coupled[cname]["A"] - c2["A"]
                df = coupled[cname]["f"] - c2["f"]
                row = {"definition": wkey, "path": pname, "case": case, "W": r[wkey], "dW": dW,
                       "Q_in": r["Q_in"], "Q_quad": r["Q_quad"], "W_p": r["W_p"], "W_v": r["W_v"],
                       "gross_abs_work": r["gross_abs_work"], "W_over_gross": r[wkey] / r["gross_abs_work"],
                       "coupled_A_mm": 1e3 * coupled[cname]["A"], "coupled_dA_mm": 1e3 * dA,
                       "coupled_df_hz": df}
                if "dW_dA" in d and abs(dW) > 0:
                    row["dA_pred_amp_only_mm"] = -1e3 * dW / d["dW_dA"]
                    if "dW_df" in d:
                        row["dA_pred_amp_freq_mm"] = -1e3 * (dW + d["dW_df"] * df) / d["dW_dA"]
                    row["implied_dA_dW_mm_per_J"] = 1e3 * dA / dW
                if "dQin_dA" in d:
                    dQ = r["Q_in"] - ref["Q_in"]
                    row["dQ_in"] = dQ
                    if abs(dQ) > 0:
                        row["dA_pred_inphase_fluid_only_mm"] = -1e3 * dQ / d["dQin_dA"]
                        row["dA_pred_inphase_with_structure_mm"] = 1e3 * dQ / (d["S_A_linear_structure"] - d["dQin_dA"])
                rows.append(row)
        out["definitions"][wkey] = d
    out["table"] = rows
    a.out.with_suffix(".json").write_text(json.dumps(out, indent=2))
    keys = sorted({k for r in rows for k in r}, key=lambda k: list(rows[0]).index(k) if k in rows[0] else 99)
    with a.out.with_suffix(".csv").open("w", newline="") as h:
        w = csv.DictWriter(h, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)
    for wkey, d in out["definitions"].items():
        print(wkey, {k: round(v, 4) for k, v in d.items()})
    for r in rows:
        print(" ".join(f"{k}={v:.4g}" if isinstance(v, float) else f"{k}={v}" for k, v in r.items()))


if __name__ == "__main__":
    main()
