#!/usr/bin/env python3
"""Amplitude slope of the wave QoIs after 6 and after 10 periods (m2/200).

usage: periods_slope.py WORK_DIR OUT_JSON [CASE]
WORK_DIR holds the study copies mu5 (a = 1) and amp4 (a = 4) of AllverifyTraction
(6 periods, analysed over periods 5-6) and p10 (a = 1) and p10a4 (a = 4) made by
prepare_study.py --periods 10 [--amplitude 4] (analysed over periods 9-10), each
with the run CASE (default robin_m2_n200). Each QoI is computed with its copy's
own verification driver (profile tolerance 0.5 at a = 4, as in
amplitude_series.py). Reported per QoI: E(1) and E(4) at both ends of the run and
the two-point slope (E(4) - E(1))/3, with the ratio of the 10-period slope to
the 6-period slope. A free mode excited at start-up and decaying over the first
periods would make that ratio fall below 1; a persistent contribution leaves it
at 1.
"""
import json
import sys
from pathlib import Path

QOIS = ("attenuation", "speed", "wallMid_amp", "wallMid_phase", "flow_amp", "flow_phase")
PAIRS = {"6": (("mu5", 1.0), ("amp4", 4.0)), "10": (("p10", 1.0), ("p10a4", 4.0))}


def analyse_copy(work: Path, study: str, a: float, case: str):
    st = (work / study / "womersleyTube").resolve()
    run = st / "verification/work" / case
    if not run.is_dir():
        sys.exit(f"periods_slope.py: required run missing: {run}")
    for mod in ("womersley_tube_verification", "womersley_exact"):
        sys.modules.pop(mod, None)
    sys.path.insert(0, str(st / "verification/scripts"))
    try:
        import womersley_tube_verification as V
        ref = json.loads(V.REFERENCE_FILE.read_text())
        coupling, m, n = case.split("_")
        r = V.analyse({"coupling": coupling, "factor": int(m[1:]), "steps": int(n[1:]),
                       "name": case, "groups": []}, ref,
                      profile_tolerance=0.5 if a >= 4.0 else 0.05)
        return {q: r[q] for q in QOIS}, ref["study"]["periods"]
    finally:
        sys.path.pop(0)


def main(argv):
    work, out = Path(argv[1]), Path(argv[2])
    case = argv[3] if len(argv) > 3 else "robin_m2_n200"
    res = {}
    for label, pair in PAIRS.items():
        values, periods = {}, None
        for study, a in pair:
            values[a], p = analyse_copy(work, study, a, case)
            if periods not in (None, p):
                sys.exit(f"{study}: {p} periods, expected {periods}")
            periods = p
        if str(periods) != label:
            sys.exit(f"copies {pair}: {periods} periods, expected {label}")
        (a1, a4) = (pair[0][1], pair[1][1])
        res[label] = {"periods": periods, "E_a1": values[a1], "E_a4": values[a4],
                      "slope_per_Pa": {q: (values[a4][q] - values[a1][q]) / (a4 - a1)
                                       for q in QOIS}}
    ratio = {q: res["10"]["slope_per_Pa"][q] / res["6"]["slope_per_Pa"][q] for q in QOIS}
    out.write_text(json.dumps({"case": case, "results": res, "slope_ratio_10_over_6": ratio},
                              indent=1))
    for q in QOIS:
        print(f"{q:14s} E(a=1) 6p {res['6']['E_a1'][q]:+.4e} 10p {res['10']['E_a1'][q]:+.4e}  "
              f"slope 6p {res['6']['slope_per_Pa'][q]:+.4e} 10p {res['10']['slope_per_Pa'][q]:+.4e}"
              f"  ratio {ratio[q]:.3f}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
