#!/usr/bin/env python3
"""Space-time limits of the womersleyTube wave QoIs.

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
   E0 is the estimate of the resolution-independent discrepancy; the spread
   between the p = 2 and observed-p variants is quoted as its uncertainty.
"""
import json
import math
import sys
from pathlib import Path

study = Path(sys.argv[1]).resolve()
sys.path.insert(0, str(study / "verification/scripts"))
import womersley_tube_verification as V  # noqa: E402

ref = json.loads(V.REFERENCE_FILE.read_text())
CP = sys.argv[3] if len(sys.argv) > 3 else "robin"
QOIS = ["attenuation", "speed", "wallMid_amp", "wallMid_phase", "flow_amp", "flow_phase"]

rows, skipped = {}, {}
for case in sorted((study / "verification/work").iterdir()):
    try:
        coupling, m, n = case.name.split("_")
        spec = {"coupling": coupling, "factor": int(m[1:]), "steps": int(n[1:]),
                "name": case.name, "groups": []}
        r = V.analyse(spec, ref)
        rows[case.name] = {q: r[q] for q in QOIS if q in r}
    except Exception as e:  # noqa: BLE001
        skipped[case.name] = str(e)
        print(f"skip {case.name}: {e}")

# the runs the limits below are built from
required = [f"{CP}_m{m}_n{s}" for m, s in ((1, 200), (2, 200), (4, 200), (4, 400),
                                           (2, 100), (2, 400))]
absent = [n for n in required if n not in rows]
if absent:
    sys.exit("attenuation_limit.py: required run(s) missing or not analysable in "
             f"{study / 'verification/work'}:\n  "
             + "\n  ".join(f"{n}: {skipped.get(n, 'not found')}" for n in absent))


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

out["mu"] = ref["reference"]["parameters"]["dynamicViscosity"]
print(json.dumps({"mu": out["mu"], "limits": out["limits"]}, indent=1, default=float))
if len(sys.argv) > 2:
    Path(sys.argv[2]).write_text(json.dumps(out, indent=1, default=float))
