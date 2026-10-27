#!/usr/bin/env python3
"""Consolidated findings of the womersleyTube interface-traction diagnostic
(mu = 5 Pa s, tutorial configuration, Robin-Neumann coupling).

usage: decomposition_summary.py WORK_DIR OUT_DIR
Inputs in WORK_DIR (produced by the other diagnostic scripts):
  attenuation_limit_mu5.json, amplitude_series.json,
  {mu5,amp2,tol}/womersleyTube/verification/postProcessing/
  interfaceTraction_summary.json
Outputs in OUT_DIR: findings_summary.json, and (if matplotlib is
available) amplitude_linearity.png, attenuation_budget.png.

Attenuation/wave-QoI budget at (m4, 200 steps/period):
  temporal part   E(m4,n200) - E(m4, dt->0), from m4 at 200/400 with q = 2;
  spatial part    -(E(m4) - E(m2))/(2^p - 1) at n200 (p = 2 and observed p);
  amplitude part  slope of the linear fit over a = 0.5, 1, 2, 4 at m2/n200,
                  times a = 1 (assumed independent of h and dt; the m1 slope,
                  from a = 1 and 4, differs from the m2 slope by 2.0 %);
  residual        the remainder.
Uncertainty of the residual: spatial part bounded by |E(m4) - E(m2)| when the
m1/m2/m4 series is non-monotone (Richardson not applicable), otherwise by the
spread between p = 2 and the observed p; plus 2 % of the amplitude part.
"""
import json
import math
import os
import sys
from pathlib import Path

import numpy as np  # noqa: E402

os.environ.setdefault("MPLCONFIGDIR", os.environ.get("TMPDIR", "/tmp"))
try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
except ImportError:      # plots are optional
    plt = None

D = Path(sys.argv[1])
R = Path(sys.argv[2])
R.mkdir(parents=True, exist_ok=True)
c = lambda v: complex(*v)  # noqa: E731

lim = json.loads((D / "attenuation_limit_mu5.json").read_text())
_amp = json.loads((D / "amplitude_series.json").read_text())
amp = {float(k): v for k, v in _amp.get("values", _amp).items()}
summ = lambda s: json.loads((D / s / "womersleyTube/verification/postProcessing/"  # noqa: E731
                             "interfaceTraction_summary.json").read_text())
S5, A2, TOL = summ("mu5"), summ("amp2"), summ("tol")
ex = S5["exact"]
out = {"configuration": "womersleyTube tutorial, mu = 5 Pa s, Robin-Neumann, BDF2",
       "waveQoIs": {}, "traction": {}, "amplitude": {}, "coupling": {}}

# --- wave QoIs: budget of the error at (m4, n200) -------------------------
A = np.array(sorted(amp))
for q in ("attenuation", "speed", "wallMid_amp", "wallMid_phase"):
    E = np.array([amp[a][q] for a in A])
    X = np.vstack([np.ones_like(A), A]).T
    (e0m2, slope), *_ = np.linalg.lstsq(X, E, rcond=None)
    L = lim["limits"][q]
    e4, e4b = L["E_m4_n200"], L["E_m4_n400"]
    t_err_n200 = e4 - (e4b + (e4b - e4) / 3)           # E(n200) - E(dt->0) on m4
    s_err = -L["space_increment_m4_p2"]                # E(m4) - E(h->0) at n200
    s_err_pobs = -L["space_increment_m4_pobs"]
    amp_part = slope * 1.0                             # tutorial amplitude 1 Pa
    resid = e4 - t_err_n200 - s_err - amp_part
    e1m, e2m = lim["cases"]["robin_m1_n200"][q], lim["cases"]["robin_m2_n200"][q]
    monotone = (e1m - e2m) * (e2m - e4) > 0
    u_space = abs(s_err - s_err_pobs) if monotone else abs(e4 - e2m)
    u_amp = 0.02 * abs(amp_part)
    out["waveQoIs"][q] = {
        "E_m4_n200": e4, "temporal_error_at_n200": t_err_n200,
        "spatial_error_m4_p2": s_err, "spatial_error_m4_pobs": s_err_pobs,
        "finite_amplitude_part": amp_part, "amplitude_slope_per_Pa_m2": slope,
        "amplitude_fit_misfit": float(np.max(np.abs(X @ [e0m2, slope] - E))),
        "residual_linear_continuum": resid,
        "mesh_series_monotone": bool(monotone),
        "residual_uncertainty_space": u_space, "residual_uncertainty_amplitude": u_amp,
        "residual_uncertainty_estimate": u_space + u_amp,
        "prediction_f1/3_infiniteTube": lim["prediction_modified_normal_viscous"]["f=0.333"][q],
        "prediction_f1/3_finiteTube": lim["prediction_modified_normal_viscous"]["f=0.333"]["finiteTube"][q],
    }

# --- traction components: space-time limits and amplitude correction -------
cs = S5["cases"]
for q in ("Vn_ratio", "VnEst_ratio", "Vn_err/|p|", "P_err/|p|", "Sn_err/|p|", "Vx_err/|Vx|"):
    e1, e2, e4, e4b = (c(cs[n]["errors"][q]) for n in
                       ("robin_m1_n200", "robin_m2_n200", "robin_m4_n200", "robin_m4_n400"))
    p = math.log(abs((e1 - e2) / (e2 - e4))) / math.log(2)
    dt = (e4b - e4) / 3
    lim2 = e4b + dt + (e4 - e2) / 3
    limp = e4b + dt + (e4 - e2) / (2 ** p - 1)
    ampl = c(A2["cases"]["robin_m2_n200"]["errors"][q]) - e2
    zero = limp - ampl
    out["traction"][q] = {"m1": [e1.real, e1.imag], "m2": [e2.real, e2.imag],
                          "m4": [e4.real, e4.imag], "m4_n400": [e4b.real, e4b.imag],
                          "mesh_order": p, "limit_a1_p2": [lim2.real, lim2.imag],
                          "limit_a1_pobs": [limp.real, limp.imag],
                          "amplitude_slope": [ampl.real, ampl.imag],
                          "limit_a0": [zero.real, zero.imag]}
Vn, P = c(ex["Vn"]), c(ex["P"])
deficit = (1 / 3 - 1) * Vn / abs(P)
out["traction"]["uncompensated_deficit_f1/3_rel_p"] = [deficit.real, deficit.imag]
# total normal traction residual (a -> 0) resolved along / across the deficit:
# along = 1 means no compensation, 0 full compensation; across is the part of
# the residual that the viscous split cannot explain
res = c(out["traction"]["Sn_err/|p|"]["limit_a0"])
proj = res * deficit.conjugate() / abs(deficit) ** 2
out["traction"]["Sn_residual_along_deficit_fraction"] = proj.real
out["traction"]["Sn_residual_across_deficit_rel_p"] = abs(proj.imag) * abs(deficit)
out["exact"] = {k: v for k, v in ex.items()}

# --- coupling tolerance --------------------------------------------------
for q in ("P_err/|p|", "Vn_err/|p|", "Sn_err/|p|", "Vx_err/|Vx|", "solid_vs_fluid_n/|p|"):
    d = c(TOL["cases"]["robin_m2_n200"]["errors"][q]) - c(cs["robin_m2_n200"]["errors"][q])
    out["coupling"][q + "_change_tol_x0.01"] = abs(d)
out["amplitude"] = {str(a): amp[a] for a in A}
(R / "findings_summary.json").write_text(json.dumps(out, indent=1, default=float))

# --- plots -----------------------------------------------------------------
if plt is None:
    print("matplotlib not available: plots skipped")
    sys.exit(0)
fig, ax = plt.subplots(figsize=(6, 4))
for q, mk in (("attenuation", "o-"), ("wallMid_phase", "s-"), ("wallMid_amp", "^-"), ("speed", "d-")):
    ax.plot(A, [amp[a][q] for a in A], mk, label=q)
ax.set_xlabel("pressure amplitude a [Pa] (tutorial: 1)")
ax.set_ylabel("relative error (phase in rad), m2, 200 steps/period")
ax.set_title("womersleyTube: error vs amplitude (linear, mesh-independent slope)")
ax.axvline(1.0, color="gray", lw=0.6, ls=":")
ax.legend(fontsize=8)
fig.tight_layout()
fig.savefig(R / "amplitude_linearity.png", dpi=130)
plt.close(fig)

b = out["waveQoIs"]["attenuation"]
parts = [("temporal (n200)", b["temporal_error_at_n200"]), ("spatial (m4)", b["spatial_error_m4_p2"]),
         ("amplitude-dependent", b["finite_amplitude_part"]), ("residual (linear continuum)", b["residual_linear_continuum"])]
fig, ax = plt.subplots(figsize=(6.5, 3.6))
ax.barh([p[0] for p in parts], [p[1] * 1e3 for p in parts], color=["#4c72b0", "#55a868", "#c44e52", "#8172b2"])
ax.axvline(b["prediction_f1/3_infiniteTube"] * 1e3, color="k", ls="--", lw=0.8, label="f=1/3, infinite tube")
ax.axvline(b["prediction_f1/3_finiteTube"] * 1e3, color="k", ls=":", lw=0.8, label="f=1/3, finite tube")
ax.set_xlabel("contribution to the Im(k) error at m4, 200 steps/period  [x 1e-3]")
ax.set_title(f"attenuation error {b['E_m4_n200']*1e3:.2f}e-3: budget")
ax.legend(fontsize=7)
fig.tight_layout()
fig.savefig(R / "attenuation_budget.png", dpi=130)
plt.close(fig)
print(json.dumps({"Sn_residual_along_deficit_fraction": out["traction"]["Sn_residual_along_deficit_fraction"],
                  "Sn_residual_across_deficit_rel_p": out["traction"]["Sn_residual_across_deficit_rel_p"],
                  "waveQoIs": out["waveQoIs"], "traction_limits_a0": {k: v["limit_a0"] for k, v in out["traction"].items() if isinstance(v, dict)},
                  "deficit": out["traction"]["uncompensated_deficit_f1/3_rel_p"], "coupling": out["coupling"]}, indent=1, default=float))
