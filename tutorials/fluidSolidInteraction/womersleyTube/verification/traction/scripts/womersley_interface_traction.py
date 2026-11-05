#!/usr/bin/env python3
"""womersleyTube interface-traction diagnostic: analysis.

Reads postProcessing/interfaceTraction/{fluid,solid}.dat written by the
optional interfaceTractionDiagnostic function object in each run of a
womersleyTube verification study, and compares the first-harmonic,
wave-projected interface tractions with the exact linear solution.

Definitions (fluid stress at the wall, outward fluid normal n = e_r,
axial tangent e_x; fields = Re[a exp(i w t)], axial dependence exp(-i k x)):
  P      = rho p_b                          (pressure)
  Vn     = rho n.(n & -devReff_b)           (transferred normal viscous)
  Vx     = rho e_x.(n & -devReff_b)         (transferred shear)
  Sn     = -P + Vn                          (total fluid sigma_rr)
  VnMom, VxMom: the same from the momentum equation's own wall flux,
         rho nuEff (laplacian boundary coefficients + n & dev2(T(grad U)))
  Ss_n   = n_s.t_applied, Ss_x = -e_x.t_applied: the traction the
           solidTraction condition applies (traction - n pressure), as
           sigma_rr and sigma_rx on the solid side
  SigN   = n_s.sigma_b.n_s: boundary value of the solid stress (with the
           PETSc SNES solid this is NOT the enforced traction; reference only)
  VnEst  = -2 rho nuEff tr(grad U)_b, the continuum normal viscous stress
           2 mu du_n/dn estimated from the fields (incompressibility)
First Fourier coefficient over the last `window` periods (all samples),
then the wave projection over the central half of the tube,
  a_hat = sum_f a(x_f) exp(+i k x_f) w_f / sum_f w_f,
which for the exact field returns its amplitude at x = 0. The weights w_f are
the recorded face areas magSf and x_f the undeformed position, both of the
solid patch (linear geometry: the solid mesh does not move), matched one-to-one
to the fluid faces on the full face-centre coordinates (tolerance 0.25 of the
smallest distance between face centres), so several faces per axial station are
handled. The moving fluid mesh is not used for the positions: it oscillates
axially by ~7e-4 m (a phase k u_x ~ 1e-4, linear in the amplitude), and its
time-mean position is displaced by the static offset D(0) of the initial solid
displacement (the fluid mesh moves by D(t) - D(0)), so the window mean of the
fluid face centres is not the undeformed position either. --fluid-mean-positions
uses it anyway (to quantify the effect).

Exact values: womersley_exact.Tube with the run's reference parameters,
using the reference A, B, k: p, 2 mu du_r/dr, sigma_rx at r = R_i.

usage: womersley_interface_traction.py STUDY_DIR [--window 2] [--out DIR]
                                        [--fluid-mean-positions]
  STUDY_DIR is a womersleyTube copy whose verification/work/<case> runs
  include the diagnostic.
"""
from __future__ import annotations

import argparse
import cmath
import csv
import json
import math
import sys
from pathlib import Path

import numpy as np


def load_exact(study: Path):
    ref = json.loads((study / "verification/reference/"
                      "womersleyTube_verification_references.json").read_text())
    sys.path.insert(0, str(study / "verification/scripts"))
    import womersley_exact as W  # the study's own (possibly regenerated) copy
    prm = ref["reference"]["parameters"]
    tube = W.Tube(prm)
    exact = ref["reference"]["exact"]
    k = complex(*exact["k"])
    A, B = complex(*exact["A"]), complex(*exact["B"])
    rows = tube.fluid_rows(k, tube.Ri)
    amp = np.array([A, B])
    p = rows["p"] @ amp
    srr = rows["srr"] @ amp
    srx = rows["srx"] @ amp
    return {
        "omega": ref["reference"]["derived"]["omega"],
        "period": ref["reference"]["derived"]["period"],
        "periods": ref["study"]["periods"], "L": prm["length"], "R": prm["innerRadius"],
        "rho": prm["fluidDensity"], "mu": prm["dynamicViscosity"], "k": k,
        "P": p, "Vn": srr + p, "Sn": srr, "Vx": srx,
    }


def read(path: Path, ncol: int):
    data = np.loadtxt(path, comments="#", ndmin=2)
    if data.shape[1] != ncol:
        raise RuntimeError(f"{path}: expected {ncol} columns, found {data.shape[1]}")
    return data


def harmonic(t, values, omega, t0, t1):
    """First Fourier coefficient a of Re[a exp(i omega t)] from all samples
    with t0 < t <= t1 (uniform), per face."""
    mask = (t > t0 + 1e-9) & (t <= t1 + 1e-9)
    tt = t[mask]
    n = tt.size
    if n < 8:
        raise RuntimeError("too few samples in the analysis window")
    dt = np.diff(np.unique(tt))
    if dt.size and np.max(np.abs(dt - dt.mean())) > 1e-6 * dt.mean():
        raise RuntimeError("non-uniform samples")
    e = np.exp(-1j * omega * tt)
    return 2.0 / n * (values[mask] * e[:, None]).sum(axis=0), n


def per_time(data, ncol_face):
    """Reshape rows (t, face, ...) into arrays [time, face, column]."""
    times = np.unique(data[:, 0])
    nface = int(data[:, 1].max()) + 1
    if data.shape[0] != times.size * nface:
        raise RuntimeError("incomplete diagnostic output")
    order = np.lexsort((data[:, 1], data[:, 0]))
    return times, data[order].reshape(times.size, nface, ncol_face)


def mean_positions(F, t, t0, t1):
    """Face centres (x, y, z) averaged over the samples of the analysis
    window (t0, t1]: the undeformed positions to first order."""
    mask = (t > t0 + 1e-9) & (t <= t1 + 1e-9)
    if not mask.any():
        raise RuntimeError("no samples in the analysis window")
    return F[mask][:, :, 2:5].mean(axis=0)


def solid_positions(S, t, t0, t1, L, require_static=True):
    """Undeformed (x, y, z) face centres and areas of the solid patch, time
    mean over the window; they must not move (linear geometry)."""
    mask = (t > t0 + 1e-9) & (t <= t1 + 1e-9)
    if not mask.any():
        raise RuntimeError("no samples in the analysis window")
    X = S[mask][:, :, 2:5]
    if require_static and np.max(np.abs(X - X[0])) > 1e-6 * L:
        raise RuntimeError("the solid interface face centres move (not linear geometry); "
                           "use --fluid-mean-positions")
    return X.mean(axis=0), S[mask][:, :, 5].mean(axis=0)


def match_faces(Xf, Xs, fraction=0.25):
    """Index of the solid face of each fluid face, from the full face-centre
    coordinates. The match must be one-to-one and closer than `fraction` times
    the smallest distance between two fluid face centres (and than the same
    for the solid); anything else raises."""
    Xf, Xs = np.asarray(Xf, float), np.asarray(Xs, float)
    if len(Xf) != len(Xs):
        raise RuntimeError(f"{len(Xf)} fluid and {len(Xs)} solid interface faces")
    if len(Xf) < 2:
        raise RuntimeError("need at least two interface faces")

    def min_spacing(X):
        d = np.linalg.norm(X[:, None, :] - X[None, :, :], axis=2)
        d[np.diag_indices(len(X))] = np.inf
        return d.min()

    h = min(min_spacing(Xf), min_spacing(Xs))
    if not h > 0:
        raise RuntimeError("coincident interface face centres")
    dist = np.linalg.norm(Xf[:, None, :] - Xs[None, :, :], axis=2)
    match = dist.argmin(axis=1)
    if np.unique(match).size != len(Xf):
        raise RuntimeError("solid interface faces are not matched one-to-one to fluid faces")
    worst = dist[np.arange(len(Xf)), match].max()
    if worst > fraction * h:
        raise RuntimeError(f"solid and fluid interface faces do not coincide: worst distance "
                           f"{worst:.3g} m > {fraction} x smallest face spacing {h:.3g} m")
    return match


def wave_projection(a, x, w, k):
    """sum_f a_f exp(+i k x_f) w_f / sum_f w_f over faces with weight w_f."""
    w = np.asarray(w, float)
    if not np.sum(w) > 0:
        raise RuntimeError("no interface faces in the projection range")
    return complex(np.sum(a * np.exp(1j * k * x) * w) / np.sum(w))


def analyse_case(case: Path, ex: dict, window: int, fluid_mean_positions: bool = False) -> dict:
    fl = read(case / "postProcessing/interfaceTraction/fluid.dat", 24)
    so = read(case / "postProcessing/interfaceTraction/solid.dat", 15)
    tf, F = per_time(fl, 24)
    ts, S = per_time(so, 15)
    if not np.allclose(tf, ts):
        raise RuntimeError("fluid and solid samples differ in time")
    rho = ex["rho"]
    t_end = tf[-1]
    t0 = t_end - window * ex["period"]
    # fluid face centres (window mean) only to match the faces; the projection
    # uses the undeformed solid positions and areas
    Xfm = mean_positions(F, tf, t0, t_end)
    Xsm, areas_s = solid_positions(S, ts, t0, t_end, ex["L"], not fluid_mean_positions)
    match = match_faces(Xfm, Xsm)          # solid face of each fluid face
    if fluid_mean_positions:
        Xf, magSf = Xfm, F[(tf > t0 + 1e-9) & (tf <= t_end + 1e-9)][:, :, 5].mean(axis=0)
    else:
        Xf, magSf = Xsm[match], areas_s[match]
    xf = Xf[:, 0]
    n = F[:, :, 6:9]
    ex_axis = np.array([1.0, 0.0, 0.0])
    # tangent: e_x projected on the face plane, normalised
    tang = ex_axis[None, None, :] - n * (n @ ex_axis)[:, :, None]
    tang /= np.linalg.norm(tang, axis=2)[:, :, None]
    P = rho * F[:, :, 9]
    tTr = F[:, :, 10:13]
    tMom = F[:, :, 13:16] + F[:, :, 16:19]
    nu = F[:, :, 23]
    Vn = rho * np.einsum("tfi,tfi->tf", tTr, n)
    Vx = rho * np.einsum("tfi,tfi->tf", tTr, tang)
    VnMom = rho * nu * np.einsum("tfi,tfi->tf", tMom, n)
    VxMom = rho * nu * np.einsum("tfi,tfi->tf", tMom, tang)
    VnEst = -2 * rho * nu * F[:, :, 22]
    # solid side, matched to fluid faces on the full face-centre coordinates
    ns = S[:, match, 6:9]
    tsol = S[:, match, 9:12]          # applied traction (solidTraction)
    tsig = S[:, match, 12:15]         # n & sigma_b (not the enforced traction)
    SsN = np.einsum("tfi,tfi->tf", tsol, ns)
    SsX = -np.einsum("tfi,tfi->tf", tsol, tang)
    SigN = np.einsum("tfi,tfi->tf", tsig, ns)
    series = {"P": P, "Vn": Vn, "Vx": Vx, "Sn": -P + Vn, "VnMom": VnMom,
              "VxMom": VxMom, "VnEst": VnEst, "SolidSn": SsN, "SolidSx": SsX,
              "SolidSigmaN": SigN}
    # central half of the tube, weighted by the face areas
    sel = (xf >= 0.25 * ex["L"]) & (xf <= 0.75 * ex["L"])
    w = np.where(sel, magSf, 0.0)
    out = {"nSamples": None, "nFaces": int(sel.sum())}
    for name, val in series.items():
        a, nsamp = harmonic(tf, val, ex["omega"], t0, t_end)
        out[name] = wave_projection(a, xf, w, ex["k"])
        out["nSamples"] = nsamp
        # max over faces of the deviation from the projected travelling wave
        model = out[name] * np.exp(-1j * ex["k"] * xf)
        scale = max(abs(out[name]), 1e-300)
        out[name + "_waveDev"] = float(np.max(np.abs((a - model)[sel])) / scale)
    return out


def errors(res: dict, ex: dict) -> dict:
    Pabs = abs(ex["P"])
    e = {
        "P_err/|p|": (res["P"] - ex["P"]) / Pabs,
        "Vn_ratio": res["Vn"] / ex["Vn"],
        "Vn_err/|p|": (res["Vn"] - ex["Vn"]) / Pabs,
        "Sn_err/|p|": (res["Sn"] - ex["Sn"]) / Pabs,
        "Vx_err/|Vx|": (res["Vx"] - ex["Vx"]) / abs(ex["Vx"]),
        "VnEst_ratio": res["VnEst"] / ex["Vn"],
        "mom_vs_transfer_n/|Vn|": (res["VnMom"] - res["Vn"]) / abs(ex["Vn"]),
        "mom_vs_transfer_x/|Vx|": (res["VxMom"] - res["Vx"]) / abs(ex["Vx"]),
        "solid_vs_fluid_n/|p|": (res["SolidSn"] - res["Sn"]) / Pabs,
        "solid_vs_fluid_x/|Vx|": (res["SolidSx"] - res["Vx"]) / abs(ex["Vx"]),
        "solidSigma_vs_fluid_n/|p|": (res["SolidSigmaN"] - res["Sn"]) / Pabs,
    }
    return e


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("study", type=Path)
    ap.add_argument("--window", type=int, default=2)
    ap.add_argument("--out", type=Path, default=None)
    ap.add_argument("--fluid-mean-positions", action="store_true",
                    help="project on the window-mean fluid face centres instead of the "
                         "undeformed solid ones")
    args = ap.parse_args()
    ex = load_exact(args.study)
    out_dir = args.out or (args.study / "verification/postProcessing")
    out_dir.mkdir(parents=True, exist_ok=True)
    print(f"exact (Pa): p {ex['P']:.6e}, 2mu u_r,r {ex['Vn']:.6e} "
          f"(|.|/|p| {abs(ex['Vn'])/abs(ex['P']):.3e}), sigma_rx {ex['Vx']:.6e}; "
          f"mu = {ex['mu']}")
    rows, summary = [], {"exact": {k: [v.real, v.imag] for k, v in ex.items()
                                   if isinstance(v, complex)}, "cases": {}}
    for case in sorted((args.study / "verification/work").iterdir()):
        if not (case / "postProcessing/interfaceTraction/fluid.dat").is_file():
            continue
        res = analyse_case(case, ex, args.window, args.fluid_mean_positions)
        err = errors(res, ex)
        summary["cases"][case.name] = {
            "values": {k: [v.real, v.imag] for k, v in res.items() if isinstance(v, complex)},
            "errors": {k: [v.real, v.imag] for k, v in err.items()},
            "waveDeviation": {k: v for k, v in res.items() if k.endswith("_waveDev")},
            "nSamples": res["nSamples"], "nFaces": res["nFaces"]}
        row = {"case": case.name}
        for k, v in err.items():
            row[k + "_re"], row[k + "_im"], row[k + "_abs"] = v.real, v.imag, abs(v)
        rows.append(row)
        print(f"{case.name:16s} " + "  ".join(
            f"{k} {v.real:+.3e}{v.imag:+.3e}i" for k, v in err.items()), flush=True)
    with (out_dir / "interfaceTraction_cases.csv").open("w", newline="") as fh:
        wr = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        wr.writeheader()
        wr.writerows(rows)
    (out_dir / "interfaceTraction_summary.json").write_text(json.dumps(summary, indent=1))
    print(f"wrote {out_dir/'interfaceTraction_cases.csv'} and interfaceTraction_summary.json")
    return 0


if __name__ == "__main__":
    sys.exit(main())
