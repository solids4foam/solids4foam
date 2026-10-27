#!/usr/bin/env python3
"""Pass/fail checks of the interface-traction diagnostic.

usage: traction_checks.py STUDY_DIR [--quick]
Checks only implementation consistency and the convergence of the shear
traction. The normal-viscous ratio, Richardson limits and amplitude-series
estimates are reported by the other scripts and are deliberately NOT
pass/fail criteria.
  1. viscous traction transferred to the solid = the momentum equation's own
     wall flux (normal and shear), relative difference <= 1e-10;
  2. traction applied by the solid condition = transferred fluid traction:
     normal <= 1e-4 of |p| (coupling-tolerance level), shear <= 1e-8 of |Vx|;
  3. (full study) wave-projected shear traction error converges on m1/m2/m4
     at an observed order >= 1.5.
"""
import json
import math
import sys
from pathlib import Path

study = Path(sys.argv[1])
quick = "--quick" in sys.argv
S = json.loads((study / "verification/postProcessing/interfaceTraction_summary.json").read_text())
fails = []
for name, cdata in S["cases"].items():
    e = {k: complex(*v) for k, v in cdata["errors"].items()}
    for key, tol in (("mom_vs_transfer_n/|Vn|", 1e-10), ("mom_vs_transfer_x/|Vx|", 1e-10),
                     ("solid_vs_fluid_n/|p|", 1e-4), ("solid_vs_fluid_x/|Vx|", 1e-8)):
        ok = abs(e[key]) <= tol
        print(f"{'PASS' if ok else 'FAIL'}: {name}: |{key}| = {abs(e[key]):.2e} (<= {tol:g})")
        if not ok:
            fails.append(f"{name} {key}")
if not quick:
    names = ["robin_m1_n200", "robin_m2_n200", "robin_m4_n200"]
    if all(n in S["cases"] for n in names):
        v = [complex(*S["cases"][n]["errors"]["Vx_err/|Vx|"]) for n in names]
        p = math.log(abs((v[0] - v[1]) / (v[1] - v[2]))) / math.log(2)
        ok = p >= 1.5
        print(f"{'PASS' if ok else 'FAIL'}: shear traction mesh order {p:.2f} (>= 1.5)")
        if not ok:
            fails.append("shear traction order")
    else:
        fails.append("mesh series incomplete")
        print("FAIL: mesh series incomplete")
print("PASSED" if not fails else f"FAILED ({len(fails)})")
sys.exit(1 if fails else 0)
