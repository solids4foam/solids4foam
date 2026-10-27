#!/usr/bin/env python3
"""Space-time limits of the womersleyTube wave QoIs and comparison with the
exact solution with a modified interface normal viscous stress.

usage: attenuation_limit.py STUDY_DIR [OUT_JSON] [COUPLING]   (COUPLING: robin | iqnils)

1. Re-analyses every completed run in STUDY_DIR/verification/work with the
   study's own verification driver (womersley_tube_verification.analyse).
2. Separates the errors E(h, dt) of attenuation, wave speed and wall
   amplitude/phase into spatial and temporal parts, assuming
       E(h, dt) = E0 + a h^p + b dt^q:
     - temporal part on m4 from (m4,n200) and (m4,n400) with q = 2
       (the observed temporal orders are 1.9-2.1), and on m2 from n100/200/400
       with the observed q;
     - spatial part on m4 at n200 from m1/m2/m4 with the observed p
       (p = 2 if the observed p is not usable);
     - E0 = E(m4,n400) + [E(m4,n400) - E(m4,n200)]/3
            + [E(m4,n200) - E(m2,n200)]/(2^p - 1).
   E0 is the estimate of the continuum (model-form) discrepancy; the spread
   between the p = 2 and observed-p variants is quoted as its uncertainty.
3. Computes, with the study's own exact-solution script, the change of the
   exact wave number and wall response when the interface normal viscous
   stress 2 mu du_r/dr is multiplied by f = 1/3 (code-path reading) and
   f = 0, for comparison with E0.
"""
import cmath
import json
import math
import sys
from pathlib import Path

import numpy as np

study = Path(sys.argv[1]).resolve()
sys.path.insert(0, str(study / "verification/scripts"))
import womersley_exact as W  # noqa: E402
import womersley_tube_verification as V  # noqa: E402

ref = json.loads(V.REFERENCE_FILE.read_text())
CP = sys.argv[3] if len(sys.argv) > 3 else "robin"
QOIS = ["attenuation", "speed", "wallMid_amp", "wallMid_phase", "flow_amp", "flow_phase"]

rows = {}
for case in sorted((study / "verification/work").iterdir()):
    try:
        coupling, m, n = case.name.split("_")
        spec = {"coupling": coupling, "factor": int(m[1:]), "steps": int(n[1:]),
                "name": case.name, "groups": []}
        r = V.analyse(spec, ref)
        rows[case.name] = {q: r[q] for q in QOIS if q in r}
    except Exception as e:  # noqa: BLE001
        print(f"skip {case.name}: {e}")


def get(name, q):
    return rows[name][q] if name in rows and q in rows[name] else None


out = {"study": str(study), "coupling": CP, "cases": rows, "limits": {}}
for q in QOIS:
    e1, e2, e4 = (get(f"{CP}_m{m}_n200", q) for m in (1, 2, 4))
    e4b = get(f"{CP}_m4_n400", q)
    t1, t2, t4 = (get(f"{CP}_m2_n{s}", q) for s in (100, 200, 400))
    res = {}
    if None not in (e1, e2, e4):
        d1, d2 = e1 - e2, e2 - e4
        p = math.log(abs(d1 / d2)) / math.log(2) if d2 != 0 and d1 / d2 > 0 else float("nan")
        res["mesh_order"] = p
        space2 = (e4 - e2) / 3.0
        spacep = (e4 - e2) / (2**p - 1) if math.isfinite(p) and 0.5 < p < 4 else space2
    if None not in (t1, t2, t4):
        d1, d2 = t1 - t2, t2 - t4
        res["time_order_m2"] = (math.log(abs(d1 / d2)) / math.log(2)
                                if d2 != 0 and d1 / d2 > 0 else float("nan"))
    if None not in (e1, e2, e4, e4b):
        # Richardson increments (added to the finer value): time on m4 with
        # q = 2, space at n200 with p = 2 or the observed p
        dtime = (e4b - e4) / 3.0
        E0_2 = e4b + dtime + space2
        E0_p = e4b + dtime + spacep
        res.update({"E_m4_n200": e4, "E_m4_n400": e4b,
                    "time_increment_m4": dtime,
                    "space_increment_m4_p2": space2, "space_increment_m4_pobs": spacep,
                    "E0_p2": E0_2, "E0_pobs": E0_p, "E0_uncertainty": abs(E0_2 - E0_p)})
    out["limits"][q] = res

# exact sensitivity to the interface normal viscous stress
prm = ref["reference"]["parameters"]
base = W.Tube(prm)
k1 = base.solve()
xi1 = base.U[base.N]
pred = {}
for f in (1 / 3, 0.0):
    class T(W.Tube):
        def fluid_rows(self, k, r, f=f):
            rr = super().fluid_rows(k, r)
            rr["srr"] = -rr["p"] + f * (rr["srr"] + rr["p"])
            return rr
    t = T(prm)
    k = t.solve()
    # wall radial displacement at x = L/2 (the driver's wallMid station):
    # amplitude U(Ri) exp(-i k x) for the same end pressure
    ratio_mid = t.U[t.N] / xi1 * cmath.exp(-1j * (k - k1) * prm["length"] / 2)
    pred[f"f={f:.3f}"] = {
        "attenuation": (k.imag - k1.imag) / abs(k1.imag),
        "speed": k1.real / k.real - 1,
        "wallMid_amp": abs(ratio_mid) - 1,
        "wallMid_phase": cmath.phase(ratio_mid),
    }
out["prediction_modified_normal_viscous"] = pred


# Finite-tube estimate: in the CFD both ends carry the exact (f = 1) data, so
# a modified interior cannot carry a single free wave. Approximation: the
# modified forward and backward modes (wave numbers k', -k', same radial
# structure) matched to the exact end pressures at r = R/2, then the
# driver's own least-squares fit of the axial pressure (61 points), and the
# wall displacement at L/2. Pressure-only matching; the end velocity-gradient
# and solid-displacement data are not represented, so this is a bracket, not
# a reference.
L = prm["length"]
xs = [L * i / (V.AXIAL_POINTS - 1) for i in range(V.AXIAL_POINTS)]
r_s = 0.5 * prm["innerRadius"]
A1 = base.A
Pexact = [A1 * W.bessel(0, k1 * r_s, True) * cmath.exp(-1j * k1 * x) for x in xs]
k_fit1 = V.fit_wave_number(xs, Pexact)
for f in (1 / 3, 0.0):
    class T2(W.Tube):
        def fluid_rows(self, k, r, f=f):
            rr = super().fluid_rows(k, r)
            rr["srr"] = -rr["p"] + f * (rr["srr"] + rr["p"])
            return rr
    t = T2(prm)
    kp = t.solve()
    g = W.bessel(0, kp * r_s, True)
    # p(x) = g (a exp(-i k' x) + b exp(i k' x)) matching both ends
    M = np.array([[g, g], [g * cmath.exp(-1j * kp * L), g * cmath.exp(1j * kp * L)]])
    a, b = np.linalg.solve(M, np.array([Pexact[0], Pexact[-1]]))
    pf = [g * (a * cmath.exp(-1j * kp * x) + b * cmath.exp(1j * kp * x)) for x in xs]
    kf = V.fit_wave_number(xs, pf)
    # wall radial displacement per unit pressure-mode coefficient (same for
    # both directions by symmetry of the radial component)
    u_per_A = t.U[t.N] / t.A
    u_mid = u_per_A * (a * cmath.exp(-1j * kp * L / 2) + b * cmath.exp(1j * kp * L / 2))
    u_mid_exact = base.U[base.N] * cmath.exp(-1j * k1 * L / 2)
    pred[f"f={f:.3f}"]["finiteTube"] = {
        "attenuation": kf.imag / k_fit1.imag - 1,
        "speed": k_fit1.real / kf.real - 1,
        "wallMid_amp": abs(u_mid) / abs(u_mid_exact) - 1,
        "wallMid_phase": cmath.phase(u_mid / u_mid_exact),
    }
out["mu"] = prm["dynamicViscosity"]
print(json.dumps({"mu": out["mu"], "limits": out["limits"], "prediction": pred}, indent=1, default=float))
if len(sys.argv) > 2:
    Path(sys.argv[2]).write_text(json.dumps(out, indent=1, default=float))
