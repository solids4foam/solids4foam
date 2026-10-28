#!/usr/bin/env python3
"""Where on the flag do the level-to-level changes of Q_in and the lift sit?

Reads faces_<level>.csv written by `tip_localisation.py extract` (windowed
first-harmonic face forces, five periods ending at the last step) for any
levels, sums per region (x bins on the upper and lower faces, face x-extent
overlap weights as in tip_localisation.py, plus the end face):

  Q_in      Re[G1y_f (a_f + i b_f)] / |D_tip|   (sums to energy_analysis Q_in)
  lift      complex first harmonic of the y force (sums to the plate lift phasor)
  lift_par  lift projected on the direction of the finest level's total lift
            phasor: the region's contribution to a change of the lift amplitude

and prints the differences between successive levels per region, with each
region's share of the total difference.

    python3 traction_regions.py --dir <dir> --levels 2x,3x,4x --json traction_regions.json
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from tip_localisation import XA, XB, classify, face_q, load_faces  # noqa: E402

EDGES = [XA, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.58, XB]


def regions():
    out = [("end face", "end", None)]
    for side in ("upper", "lower"):
        for a, b in zip(EDGES[:-1], EDGES[1:]):
            out.append((f"{side} {a:.2f}-{b:.2f}", side, (a, b)))
    return out


def sums(c, d_tip):
    xc, yc, kind = classify(c)
    qc, _, g = face_q(c)
    res = {}
    for name, k, iv in regions():
        if k == "end":
            w = (kind == "end").astype(float)
        else:
            ov = np.clip(np.minimum(c["x1"], iv[1]) - np.maximum(c["x0"], iv[0]), 0, None)
            w = np.where(kind == k, ov / np.maximum(c["x1"] - c["x0"], 1e-12), 0.0)
        res[name] = {"Q_in": float(np.sum(w * qc.real) / d_tip), "lift": complex(np.sum(w * g))}
    return res


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--dir", type=Path, required=True)
    ap.add_argument("--levels", default="2x,3x,4x")
    ap.add_argument("--json", type=Path)
    a = ap.parse_args()
    lv = a.levels.split(",")
    data = {}
    for l in lv:
        meta, c = load_faces(a.dir / f"faces_{l}.csv")
        data[l] = sums(c, meta["d_tip"])
    names = [r[0] for r in regions()]
    tot = {l: {"Q_in": sum(data[l][n]["Q_in"] for n in names), "lift": sum(data[l][n]["lift"] for n in names)} for l in lv}
    u = tot[lv[-1]]["lift"] / abs(tot[lv[-1]]["lift"])
    print("totals: " + "  ".join(f"{l}: Q_in {tot[l]['Q_in']:.3f} lift ({tot[l]['lift'].real:.2f}, {tot[l]['lift'].imag:.2f}) "
                                 f"|lift| {abs(tot[l]['lift']):.3f}" for l in lv))
    out = {"levels": lv, "totals": {l: {"Q_in": tot[l]["Q_in"], "lift_re": tot[l]["lift"].real,
                                        "lift_im": tot[l]["lift"].imag} for l in lv}, "regions": {}}
    pairs = list(zip(lv[:-1], lv[1:]))
    hdr = "region".ljust(18) + "".join(f"  Q_in {l:>3s}" for l in lv) + "".join(
        f"  dQin {p[0]}-{p[1]} (share)" for p in pairs) + "".join(f"  dlift_par {p[0]}-{p[1]} (share)" for p in pairs)
    print(hdr)
    for n in names:
        row = {"Q_in": [data[l][n]["Q_in"] for l in lv],
               "lift_re": [data[l][n]["lift"].real for l in lv], "lift_im": [data[l][n]["lift"].imag for l in lv]}
        s = n.ljust(18) + "".join(f"  {v:8.3f}" for v in row["Q_in"])
        for p in pairs:
            d = data[p[1]][n]["Q_in"] - data[p[0]][n]["Q_in"]
            D = tot[p[1]]["Q_in"] - tot[p[0]]["Q_in"]
            row[f"dQ_in_{p[0]}_{p[1]}"] = d
            s += f"  {d:8.3f} ({d / D:6.2f})      "
        for p in pairs:
            d = ((data[p[1]][n]["lift"] - data[p[0]][n]["lift"]) * u.conjugate()).real
            D = ((tot[p[1]]["lift"] - tot[p[0]]["lift"]) * u.conjugate()).real
            row[f"dlift_par_{p[0]}_{p[1]}"] = d
            s += f"  {d:8.3f} ({d / D:6.2f})          "
        out["regions"][n] = row
        print(s)
    for p in pairs:
        D = ((tot[p[1]]["lift"] - tot[p[0]]["lift"]) * u.conjugate()).real
        print(f"total {p[0]}-{p[1]}: dQ_in {tot[p[1]]['Q_in'] - tot[p[0]]['Q_in']:.3f}, dlift_par {D:.3f}, "
              f"d|lift| {abs(tot[p[1]]['lift']) - abs(tot[p[0]]['lift']):.3f}")
    if a.json:
        a.json.write_text(json.dumps(out, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
