#!/usr/bin/env python3
"""Spatial decomposition of the replay flag traction (tip localisation).

extract  reads the traction_<rank>.bin files of a hron_turek_tip.py case and
         writes one CSV row per flag face with the windowed first-harmonic
         phasors of the pressure and viscous traction (same window and phase
         convention as energy_analysis.py: last --periods periods, ending at
         the last step, phase exp(-i omega T), T = t - t0).
analyse  sums the face contributions to Q_in over geometric regions (face
         x-extent overlap weights, so regions are comparable across meshes),
         checks the sum against energy_analysis.py, and writes CSV/JSON/PNG.

Q_in contribution of a face: Re[ G1y_f (a_f + i b_f) ] / |D_tip|, with G1y_f the
complex first harmonic of the face y force (N/m) and (a_f, b_f) the face's
h = 1, y replay coefficients; the sum over faces equals energy_analysis Q_in.
"""
from __future__ import annotations
import argparse, csv, json, math, sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
XA, XB, THK = 0.24899, 0.6, 0.02   # flag attachment x (cylinder surface), tip x, thickness


def read_bins(case: Path):
    faces, recs = [], []
    for f in sorted(case.glob("traction_*.bin"), key=lambda p: int(p.stem.split("_")[1])):
        b = f.read_bytes()
        nf, nc = (int(v) for v in np.frombuffer(b[:16], dtype=np.int64))
        nh = 5 + 2 * nc
        hd = np.frombuffer(b[16:16 + 8 * nf * nh], dtype=np.float64).reshape(nf, nh)
        r = np.frombuffer(b[16 + 8 * nf * nh:], dtype=np.float64)
        nrec = len(r) // (1 + 4 * nf)
        r = r[:nrec * (1 + 4 * nf)].reshape(nrec, 1 + 4 * nf)
        faces.append(hd)
        recs.append(r)
    nrec = min(len(r) for r in recs)
    t = recs[0][:nrec, 0]
    for r in recs:
        assert np.allclose(r[:nrec, 0], t)
    return np.vstack(faces), t, np.hstack([r[:nrec, 1:] for r in recs]), nc


def extract(a):
    from energy_analysis import setup, plate_kinematics
    case = a.case
    s = setup(case)
    t0 = float(s["t0"])
    omega, d_tip, *_ = plate_kinematics(case, 1.0)
    hd, t, rec, nc = read_bins(case)
    _, idx = np.unique(t, return_index=True)
    idx = np.sort(idx)
    t, rec = t[idx], rec[idx]
    period = 2 * math.pi / omega
    sel = (t > t[-1] - a.periods * period) & (t <= t[-1])
    n = int(sel.sum())
    ph = np.exp(-1j * omega * (t[sel] - t0))
    r = rec[sel].reshape(n, -1, 4)
    out = []
    for i in range(hd.shape[0]):
        fpx, fpy, fvx, fvy = (r[:, i, k] for k in range(4))
        row = dict(x0=hd[i, 0], x1=hd[i, 1], y0=hd[i, 2], y1=hd[i, 3], len=hd[i, 4],
                   a1y=hd[i, 5 + nc + 1], b1y=hd[i, 5 + nc + 2], a1x=hd[i, 5 + 1], b1x=hd[i, 5 + 2],
                   mean_fpx=fpx.mean(), mean_fvx=fvx.mean(), mean_fpy=fpy.mean(), mean_fvy=fvy.mean())
        for tag, v in (("fpx", fpx), ("fvx", fvx), ("fpy", fpy), ("fvy", fvy)):
            g = 2.0 / n * np.sum(v * ph)
            row["G1_" + tag + "_re"], row["G1_" + tag + "_im"] = g.real, g.imag
        out.append(row)
    order = sorted(range(len(out)), key=lambda i: (out[i]["y0"] + out[i]["y1"] < 0.4 and 0, out[i]["x0"], out[i]["y0"]))
    out = [out[i] for i in order]
    meta = dict(case=case.name, omega=omega, d_tip=abs(d_tip), t_window=[float(t[sel][0]), float(t[-1])],
                steps=n, period=period, nfaces=len(out))
    with open(a.out, "w", newline="") as fh:
        fh.write("# " + json.dumps(meta) + "\n")
        w = csv.DictWriter(fh, fieldnames=list(out[0]))
        w.writeheader()
        w.writerows(out)
    print(meta)


def load_faces(path: Path):
    lines = path.read_text().splitlines()
    meta = json.loads(lines[0][2:])
    rows = list(csv.DictReader(lines[1:]))
    cols = {k: np.array([float(r[k]) for r in rows]) for k in rows[0]}
    return meta, cols


def classify(c):
    xc = 0.5 * (c["x0"] + c["x1"])
    yc = 0.5 * (c["y0"] + c["y1"])
    end = (c["x1"] - c["x0"]) < 1e-9
    kind = np.where(end, "end", np.where(yc > 0.2, "upper", "lower"))
    return xc, yc, kind


def face_q(c):
    """Per-face Q_in, Q_quad contributions (unnormalised by |D_tip|)."""
    g = (c["G1_fpy_re"] + c["G1_fvy_re"]) + 1j * (c["G1_fpy_im"] + c["G1_fvy_im"])
    gp = c["G1_fpy_re"] + 1j * c["G1_fpy_im"]
    ab = c["a1y"] + 1j * c["b1y"]
    return g * ab, gp * ab, g


def regions(L=THK):
    """Ordered list of (name, kind, x-interval) partitions of the flag surface."""
    e = [XB - 0.1 * L, XB - L]   # tip band edges
    ca = [XA + 0.1 * L, XA + L]  # cylinder band edges
    mid = np.linspace(ca[1], e[1], 6)
    R = [("tip_end", "end", None),
         ("tip_band_0.1t", "side", (e[0], XB)),
         ("tip_band_0.1t-1t", "side", (e[1], e[0])),
         ("cyl_band_0.1t", "side", (XA, ca[0])),
         ("cyl_band_0.1t-1t", "side", (ca[0], ca[1]))]
    for i in range(5):
        R.append((f"mid_{i + 1}", "side", (mid[i], mid[i + 1])))
    return R


def region_sums(c, d_tip):
    xc, yc, kind = classify(c)
    qc, qcp, g = face_q(c)
    side = kind != "end"
    out = {}
    for name, k, iv in regions():
        if k == "end":
            w = (kind == "end").astype(float)
        else:
            ov = np.clip(np.minimum(c["x1"], iv[1]) - np.maximum(c["x0"], iv[0]), 0, None)
            w = np.where(side, ov / np.maximum(c["x1"] - c["x0"], 1e-12), 0.0)
        length = float(np.sum(w * c["len"]))
        out[name] = dict(
            length=length, Q_in=float(np.sum(w * qc.real) / d_tip), Q_quad=float(np.sum(w * qc.imag) / d_tip),
            Q_in_p=float(np.sum(w * qcp.real) / d_tip),
            lift_re=float(np.sum(w * g.real)), lift_im=float(np.sum(w * g.imag)),
            drag_mean=float(np.sum(w * (c["mean_fpx"] + c["mean_fvx"]))),
            nfaces=int(np.sum(w > 0)))
    return out


def obs_order(q1, q2, q4):
    d12, d24 = q2 - q1, q4 - q2
    mono = d12 * d24 > 0
    if not mono or abs(d12) < 1e-12 or abs(d24) < 1e-12:
        return d12, d24, None, "non-monotone" if not mono else "no change"
    r = d12 / d24
    return d12, d24, (math.log2(r) if r > 1 else None), ("ok" if r > 1 else "diverging (ratio<1)")


def analyse(a):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    lev = ["1x", "2x", "4x"]
    data = {}
    for l in lev:
        meta, c = load_faces(a.dir / f"faces_{l}.csv")
        data[l] = (meta, c, region_sums(c, meta["d_tip"]))
    names = [r[0] for r in regions()]
    groups = {"tip_corner_0.1t": ["tip_end", "tip_band_0.1t"],
              "tip_corner_1t": ["tip_end", "tip_band_0.1t", "tip_band_0.1t-1t"],
              "cyl_corner_1t": ["cyl_band_0.1t", "cyl_band_0.1t-1t"],
              "smooth_mid": [f"mid_{i}" for i in range(1, 6)]}
    res = {"levels": lev, "regions": {}}
    for n, mem in list(((n, [n]) for n in names)) + list(groups.items()):
        row = {}
        for q in ("Q_in", "Q_quad", "lift_re", "lift_im", "drag_mean", "length"):
            row[q] = [sum(data[l][2][m][q] for m in mem) for l in lev]
        row["Q_in_per_len"] = [qq / ll if ll > 0 else float("nan") for qq, ll in zip(row["Q_in"], row["length"])]
        row["lift_amp"] = [math.hypot(x, y) for x, y in zip(row["lift_re"], row["lift_im"])]
        d12, d24, p, st = obs_order(*row["Q_in"])
        row.update(dQ_1x_2x=d12, dQ_2x_4x=d24, order=p, status=st)
        d12, d24, p, st = obs_order(*row["drag_mean"])
        row.update(d_drag_1x_2x=d12, d_drag_2x_4x=d24, drag_order=p, drag_status=st)
        res["regions"][n] = row
    tot = {l: sum(data[l][2][n]["Q_in"] for n in names) for l in lev}
    res["totals"] = {"Q_in_sum_regions": tot,
                     "lift_amp_sum_regions": {l: math.hypot(sum(data[l][2][n]["lift_re"] for n in names),
                                                            sum(data[l][2][n]["lift_im"] for n in names)) for l in lev},
                     "drag_mean_sum_regions": {l: sum(data[l][2][n]["drag_mean"] for n in names) for l in lev},
                     "meta": {l: data[l][0] for l in lev}}
    if a.baseline:
        res["totals"]["baseline_energy_analysis"] = json.loads(a.baseline.read_text())
    (a.dir / "tip_localisation.json").write_text(json.dumps(res, indent=1))
    with open(a.dir / "tip_localisation_regions.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["region", "length_m"] + [f"Q_in_{l}" for l in lev] + [f"Q_in_per_m_{l}" for l in lev] +
                   ["d_1x_2x", "d_2x_4x", "order", "status"] + [f"lift_amp_{l}" for l in lev] +
                   [f"drag_mean_{l}" for l in lev] + ["d_drag_2x_4x"])
        for n, r in res["regions"].items():
            w.writerow([n, f"{r['length'][2]:.5g}"] + [f"{x:.4f}" for x in r["Q_in"]] +
                       [f"{x:.2f}" for x in r["Q_in_per_len"]] + [f"{r['dQ_1x_2x']:.4f}", f"{r['dQ_2x_4x']:.4f}",
                       "" if r["order"] is None else f"{r['order']:.2f}", r["status"]] +
                       [f"{x:.3f}" for x in r["lift_amp"]] + [f"{x:.4f}" for x in r["drag_mean"]] +
                       [f"{r['d_drag_2x_4x']:.4f}"])
    # plots
    col = {"1x": "#1b9e77", "2x": "#d95f02", "4x": "#7570b3"}
    fig, ax = plt.subplots(3, 2, figsize=(11, 9), sharex=True)
    for l in lev:
        meta, c, _ = data[l]
        xc, yc, kind = classify(c)
        gy = (c["G1_fpy_re"] + c["G1_fvy_re"]) + 1j * (c["G1_fpy_im"] + c["G1_fvy_im"])
        qc = face_q(c)[0].real / meta["d_tip"]
        for j, side in enumerate(("upper", "lower")):
            m = kind == side
            o = np.argsort(xc[m])
            x = xc[m][o]
            tr = (gy[m] / c["len"][m])[o]       # first-harmonic y traction, N/m per m = Pa
            ax[0, j].plot(x, np.abs(tr), "-", color=col[l], label=l, lw=1.2)
            ax[1, j].plot(x, np.degrees(np.unwrap(np.angle(tr))), "-", color=col[l], lw=1.2)
            ax[2, j].plot(x, (qc[m] / c["len"][m])[o], "-", color=col[l], lw=1.2)
            ax[0, j].set_title(f"{side} face")
    ax[0, 0].set_ylabel("|first harmonic of t_y| (Pa)")
    ax[1, 0].set_ylabel("phase of t_y (deg)")
    ax[2, 0].set_ylabel("Q_in density (N/m per m)")
    for j in range(2):
        ax[2, j].set_xlabel("x (m)")
    ax[0, 0].legend()
    for r in ax:
        for x in r:
            x.grid(alpha=0.3)
    fig.suptitle("Replay flag traction, first harmonic (last five periods)")
    fig.tight_layout()
    fig.savefig(a.dir / "traction_distribution.png", dpi=130)
    plt.close(fig)
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    for l in lev:
        meta, c, _ = data[l]
        xc, yc, kind = classify(c)
        qc = face_q(c)[0].real / meta["d_tip"]
        # cumulative Q_in from the tip: end face plus side faces at x > xs, face-overlap weights
        xs = np.linspace(XB, XA, 400)
        cum = []
        for x in xs:
            side = kind != "end"
            ov = np.clip(c["x1"] - np.maximum(c["x0"], x), 0, None) / np.maximum(c["x1"] - c["x0"], 1e-12)
            cum.append(np.sum(np.where(side, ov, 0.0) * qc) + np.sum(qc[kind == "end"]))
        ax[0].plot((XB - xs) / THK, cum, color=col[l], label=l)
        ax[1].plot((XB - xs) / THK, np.array(cum) / meta["d_tip"] * 0 + np.array(cum) / sum(qc), color=col[l], label=l)
    ax[0].set_ylabel("cumulative Q_in from the tip (N/m)")
    ax[1].set_ylabel("fraction of total Q_in")
    for x in ax:
        x.set_xlabel("distance from the tip / thickness")
        x.set_xlim(0, 4)
        x.grid(alpha=0.3)
    ax[0].legend()
    fig.tight_layout()
    fig.savefig(a.dir / "tip_cumulative.png", dpi=130)
    print(json.dumps({"tot": res["totals"]["Q_in_sum_regions"]}))


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sp = ap.add_subparsers(dest="cmd", required=True)
    e = sp.add_parser("extract")
    e.add_argument("case", type=Path)
    e.add_argument("out", type=Path)
    e.add_argument("--periods", type=int, default=5)
    n = sp.add_parser("analyse")
    n.add_argument("dir", type=Path)
    n.add_argument("--baseline", type=Path)
    a = ap.parse_args()
    extract(a) if a.cmd == "extract" else analyse(a)


if __name__ == "__main__":
    main()
