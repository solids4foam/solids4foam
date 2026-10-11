#!/usr/bin/env python3
"""Per-cycle fluid work on the FSI3 flag from the energy restarts.

For every case written by hron_turek_energy.py (energy.dat) this evaluates,
over the last --periods full periods of the prescribed motion (cumulative work
interpolated at the window ends):

  W      net work of the fluid on the flag per cycle, sum_n P^n dt (J/m),
         split into pressure and viscous parts; positive = into the structure
  W_hc   the same split by harmonic h and direction c of the replayed motion
         (discrete work increments), whose sum must reproduce W
  G1     complex first harmonic of the generalised force on the h = 1, y shape,
         G1 = sum_f F_f conj(D_f), normalised by the tip amplitude |D_tip|:
         Q_in = Re G1/|D_tip| (in phase with the displacement: added stiffness),
         Q_quad = Im G1/|D_tip| (in phase with the velocity: negative damping),
         W_1y ~= pi A_tip Q_quad
  E1     energy of the h = 1, y motion of the plate as a thin beam,
         0.5 rho_s t_p omega^2 int |D_1y|^2 dx, used to normalise W

and a restart check against the pressure-power trace (POWR) of the original
replay run.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np

RHO_S, T_P, Y_LO, Y_HI = 1000.0, 0.02, 0.19, 0.21
COLS = ["t", "T", "g", "Pp", "Pv", "Fpx", "Fpy", "Fvx", "Fvy",
        "dW1x", "dW1y", "dW2x", "dW2y", "dW3x", "dW3y", "dW4x", "dW4y", "dW0x", "dW0y",
        "Q1xa", "Q1xb", "Q1ya", "Q1yb", "Q2xa", "Q2xb", "Q2ya", "Q2yb", "Qp1ya", "Qp1yb"]


def setup(case: Path) -> dict:
    d = {}
    for line in (case / "energy_setup.txt").read_text().splitlines():
        k, v = line.split(maxsplit=1)
        d[k] = v
    return d


def table(case: Path):
    lines = (case / "constant/plateMotion.tab").read_text().split("\n")
    omega, _, ncoef = lines[0].split()
    rows = {}
    for line in lines[1:]:
        if line.strip():
            v = line.split()
            rows[(int(v[0]), int(v[1]))] = np.array([float(x) for x in v[2:]])
    return float(omega), int(ncoef), rows


def plate_kinematics(case: Path, k: float):
    """Tip h=1 y phasor, E1, inextensible tip shortening, fitted tip x mean."""
    omega, nc, rows = table(case)
    key = min(rows, key=lambda q: math.hypot(q[0] / 1e5 - 0.6, q[1] / 1e5 - 0.19))
    c = rows[key]
    d_tip = complex(c[nc + 1], -c[nc + 2])
    lo = {q[0]: rows[q] for q in rows if q[1] == round(Y_LO * 1e5)}
    hi = {q[0]: rows[q] for q in rows if q[1] == round(Y_HI * 1e5)}
    xs = sorted(set(lo) & set(hi))
    x = np.array(xs) / 1e5

    def centre(h):
        return np.array([complex(0.5 * (lo[q][nc + 2 * h - 1] + hi[q][nc + 2 * h - 1]),
                                 -0.5 * (lo[q][nc + 2 * h] + hi[q][nc + 2 * h])) for q in xs])

    e1 = 0.5 * RHO_S * T_P * omega ** 2 * np.trapezoid(np.abs(k * centre(1)) ** 2, x)
    short = 0.0
    for h in range(1, nc // 2 + 1):
        slope = np.gradient(centre(h), x)
        short += 0.25 * np.trapezoid(np.abs(k * slope) ** 2, x)
    xmean = k * 0.5 * (lo[xs[-1]][0] + hi[xs[-1]][0])
    return omega, d_tip, e1, short, xmean, float(x[0]), float(x[-1])


def analyse(case: Path, periods: int, source_log: Path | None) -> dict:
    s = setup(case)
    k, r = float(s["amplitude_factor"]), float(s["frequency_factor"])
    data = np.loadtxt(case / "energy.dat")
    _, idx = np.unique(data[:, 0], return_index=True)
    data = data[np.sort(idx)]
    col = {n: data[:, i] for i, n in enumerate(COLS)}
    t = col["t"]
    dt = float(np.median(np.diff(t)))
    omega, d_tip, e1, short, xmean, x0, x1 = plate_kinematics(case, k)
    f = omega * r / (2 * math.pi)
    period = 1.0 / f
    out = {"case": case.name, "amplitude_factor": k, "frequency_factor": r, "dt": dt,
           "frequency_hz": f, "A_tip_m": k * abs(d_tip), "E1_J_per_m": e1,
           "inextensible_ux_mean_tip_m": -short, "fitted_ux_mean_tip_centre_m": xmean,
           "centreline_x_range": [x0, x1], "steps": int(len(t)),
           "t_first": float(t[0]), "t_last": float(t[-1])}
    tc = np.concatenate([[t[0] - dt], t])

    def cum(series):
        return np.concatenate([[0.0], np.cumsum(series)])

    def window(t1):
        t0 = t1 - periods * period
        res = {"window": [t0, t1]}
        for name, series in (("W_p", col["Pp"] * dt), ("W_v", col["Pv"] * dt)):
            c = cum(series)
            res[name] = (np.interp(t1, tc, c) - np.interp(t0, tc, c)) / periods
        res["W"] = res["W_p"] + res["W_v"]
        for name in COLS[9:19]:
            c = cum(col[name])
            res[name.replace("dW", "W_")] = (np.interp(t1, tc, c) - np.interp(t0, tc, c)) / periods
        res["W_modal_sum"] = sum(res[n.replace("dW", "W_")] for n in COLS[9:19])
        sel = (t > t0) & (t <= t1)
        res["gross_abs_work"] = float(np.sum(np.abs(col["Pp"][sel] + col["Pv"][sel])) * dt / periods)
        ph = np.exp(-1j * omega * col["T"][sel])
        n = int(sel.sum())
        for tag, a, b in (("", "Q1ya", "Q1yb"), ("_p", "Qp1ya", "Qp1yb")):
            g = 2.0 / n * np.sum((col[a][sel] + 1j * col[b][sel]) * ph)
            res["Q_in" + tag] = g.real / abs(d_tip)
            res["Q_quad" + tag] = g.imag / abs(d_tip)
        res["W1y_from_Qquad"] = math.pi * k * abs(d_tip) * res["Q_quad"]
        fy = col["Fpy"][sel] + col["Fvy"][sel]
        fx = col["Fpx"][sel] + col["Fvx"][sel]
        res["plate_lift_amp1"] = float(abs(2.0 / n * np.sum(fy * ph)))
        res["plate_drag_mean"] = float(fx.mean())
        return res

    main = window(t[-1])
    alt = window(t[-1] - 0.5 * period)
    out.update(main)
    out["W_window_shift_change"] = main["W"] - alt["W"]
    out["W_over_E1"] = main["W"] / e1
    if source_log and source_log.is_file():
        pw = {}
        for line in source_log.read_text(errors="replace").splitlines():
            if line.startswith("POWR "):
                _, tt, v = line.split()
                pw[round(float(tt), 7)] = float(v)
        common = [i for i, tt in enumerate(t) if round(tt, 7) in pw and tt > t[0] + 0.2]
        if common:
            a = np.array([col["Pp"][i] for i in common])
            b = np.array([pw[round(t[i], 7)] for i in common])
            out["restart_check_interval"] = [float(t[common[0]]), float(t[common[-1]])]
            out["restart_check_max_abs_dP"] = float(np.max(np.abs(a - b)))
            out["restart_check_rms_P"] = float(np.sqrt(np.mean(b ** 2)))
            out["restart_check_work_diff_per_cycle"] = float(np.sum(a - b) * dt / ((t[common[-1]] - t[common[0]]) / period))
    return out


KEYS = ["W", "W_p", "W_v", "W_1y", "W_1x", "W_2x", "W_2y", "W_0x", "W_0y", "W_modal_sum",
        "W_window_shift_change", "Q_in", "Q_quad", "Q_in_p", "Q_quad_p", "W1y_from_Qquad", "A_tip_m",
        "E1_J_per_m", "W_over_E1", "gross_abs_work", "plate_lift_amp1", "plate_drag_mean",
        "inextensible_ux_mean_tip_m", "fitted_ux_mean_tip_centre_m",
        "restart_check_max_abs_dP", "restart_check_rms_P", "restart_check_work_diff_per_cycle"]


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("cases", nargs="+", type=Path)
    ap.add_argument("--periods", type=int, default=5)
    ap.add_argument("--json", type=Path)
    a = ap.parse_args()
    res = []
    for c in a.cases:
        src = Path(setup(c)["source"]) / "log.solids4Foam"
        res.append(analyse(c, a.periods, src))
    for r in res:
        print("==", r["case"], "window", [round(x, 4) for x in r["window"]])
        for kk in KEYS:
            if kk in r:
                print(f"   {kk:34s} {r[kk]: .6g}")
    if a.json:
        a.json.write_text(json.dumps(res, indent=2))


if __name__ == "__main__":
    main()
