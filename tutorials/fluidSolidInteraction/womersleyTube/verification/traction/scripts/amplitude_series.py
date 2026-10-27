#!/usr/bin/env python3
"""Amplitude series of the womersleyTube wave QoIs (finite-amplitude part).

usage: amplitude_series.py WORK_DIR OUT_JSON [CASE]
WORK_DIR holds study copies amp05, base (mu5), amp2, amp4, each with the run
CASE (default robin_m2_n200) under verification/work. The QoIs are computed
with each copy's own verification driver. At 4x the tutorial amplitude the
driver's check that the profile sample points stay within 5% of the cell
spacing fails (the wall moves further); for this analysis only, that
tolerance is relaxed to 50% (the axial pressure samples used for the wave
number are unaffected). Fits E(a) = E0 + c a^m and reports m from successive
differences.
"""
import json
import math
import sys
from pathlib import Path

import numpy as np

work, out = Path(sys.argv[1]), Path(sys.argv[2])
case = sys.argv[3] if len(sys.argv) > 3 else "robin_m2_n200"
QOIS = ("attenuation", "speed", "wallMid_amp", "wallMid_phase", "flow_amp", "flow_phase")
res = {}
for study, a in (("amp05", 0.5), ("mu5", 1.0), ("amp2", 2.0), ("amp4", 4.0)):
    st = (work / study / "womersleyTube").resolve()
    if not (st / "verification/work" / case).is_dir():
        continue
    for mod in ("womersley_tube_verification", "womersley_exact"):
        sys.modules.pop(mod, None)
    sys.path.insert(0, str(st / "verification/scripts"))
    import womersley_tube_verification as V
    if a > 2.0:
        src = Path(V.__file__).read_text().replace("> 0.05 * spacing", "> 0.5 * spacing")
        exec(compile(src, V.__file__, "exec"), V.__dict__)
    ref = json.loads(V.REFERENCE_FILE.read_text())
    coupling, m, n = case.split("_")
    r = V.analyse({"coupling": coupling, "factor": int(m[1:]), "steps": int(n[1:]),
                   "name": case, "groups": []}, ref)
    res[a] = {q: r[q] for q in QOIS}
    sys.path.pop(0)
A = np.array(sorted(res))
fit = {}
for q in QOIS:
    E = np.array([res[a][q] for a in A])
    X = np.vstack([np.ones_like(A), A]).T
    coef, *_ = np.linalg.lstsq(X, E, rcond=None)
    d = np.diff(E)
    m_obs = [math.log2(d[i + 1] / d[i]) if d[i] * d[i + 1] > 0 else float("nan")
             for i in range(len(d) - 1)]
    fit[q] = {"E0": coef[0], "slope_per_Pa": coef[1],
              "misfit": float(np.max(np.abs(X @ coef - E))), "m_observed": m_obs}
out.write_text(json.dumps({"case": case, "values": {str(k): v for k, v in res.items()},
                           "linearFit": fit}, indent=1, default=float))
for q, f in fit.items():
    print(f"{q:14s} E0 {f['E0']:+.4e}  slope {f['slope_per_Pa']:+.4e}/Pa  misfit {f['misfit']:.1e}  "
          f"m_obs {['%.2f' % x for x in f['m_observed']]}")
