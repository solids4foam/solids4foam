#!/usr/bin/env python3
"""Prepare a womersleyTube study copy with the interface-traction diagnostic.

usage: prepare_study.py SRC_TUTORIAL DEST [--mu MU] [--amplitude A]
                        [--coupling-tol-factor F]

- copies the tutorial (without results) to DEST;
- includes system/interfaceTractionDiagnostic in the controlDict functions,
  sampling every time step over the last two periods (timeStart = 4 periods);
- --mu / --amplitude: change the dynamic viscosity (Pa s) or the pressure
  amplitude (Pa) and REGENERATE, consistently, constant/womersleyProperties
  (end data and initial fields) and the reference JSON with the copy's own
  womersley_exact.py; constant/fluid/transportProperties nu is set to mu/rho;
- --coupling-tol-factor: multiply the Robin and IQN-ILS coupling tolerances.
Nothing in SRC is modified.
"""
import argparse
import json
import re
import shutil
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parents[3] / "system"   # tutorial system/


def sub_one(path: Path, pattern: str, repl: str, count_expected=1):
    text = path.read_text()
    new, n = re.subn(pattern, repl, text, flags=re.MULTILINE)
    if n != count_expected:
        raise SystemExit(f"{path}: expected {count_expected} match(es) of {pattern!r}, found {n}")
    path.write_text(new)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("src", type=Path)
    ap.add_argument("dest", type=Path)
    ap.add_argument("--mu", type=float)
    ap.add_argument("--amplitude", type=float)
    ap.add_argument("--coupling-tol-factor", type=float)
    a = ap.parse_args()
    if a.dest.exists():
        raise SystemExit(f"{a.dest} exists")
    shutil.copytree(a.src, a.dest, symlinks=True,
                    ignore=shutil.ignore_patterns("postProcessing", "work", "dynamicCode",
                                                  "log.*", "processor*"))
    d = a.dest
    # diagnostic function object
    shutil.copy(HERE / "interfaceTractionDiagnostic", d / "system")
    fo = d / "system/interfaceTractionDiagnostic"
    sub_one(fo, r"^(\s*executeInterval\s+)\d+;", r"\g<1>1;")
    sub_one(fo, r"^(\s*writeInterval\s+)\d+;", r"\g<1>1;")
    ref_file = d / "verification/reference/womersleyTube_verification_references.json"
    ref = json.loads(ref_file.read_text())
    period = ref["reference"]["derived"]["period"]
    periods = ref["study"]["periods"]
    window = ref["study"]["analysisPeriods"]
    sub_one(fo, r"^(\s*timeStart\s+)[^;]+;", rf"\g<1>{(periods - window) * period:.12g};")
    sub_one(d / "system/controlDict", r"^functions\s*\n\{",
            'functions\n{\n    #include "interfaceTractionDiagnostic"\n')
    changes = {}
    if a.mu is not None or a.amplitude is not None:
        sys.path.insert(0, str(d / "verification/scripts"))
        import womersley_exact as W
        if a.mu is not None:
            W.PARAMETERS["dynamicViscosity"] = a.mu
            changes["dynamicViscosity"] = a.mu
        if a.amplitude is not None:
            W.PARAMETERS["pressureAmplitude"] = a.amplitude
            changes["pressureAmplitude"] = a.amplitude
        tube = W.Tube(W.PARAMETERS)
        tube.solve()
        (d / "constant/womersleyProperties").write_text(tube.properties())
        ref["reference"] = W.reference()
        ref["description"] = (ref.get("description", "") +
                              f" [regenerated for the diagnostic with {changes}]")
        ref_file.write_text(json.dumps(ref, indent=2))
        # make the copy's exact script self-consistent with the new parameters
        script = d / "verification/scripts/womersley_exact.py"
        for key, val in changes.items():
            sub_one(script, rf'^(\s*"{key}":\s*)[^,]+,', rf"\g<1>{val!r},")
        if a.mu is not None:
            nu = a.mu / W.PARAMETERS["fluidDensity"]
            sub_one(d / "constant/fluid/transportProperties",
                    r"^(nu\s+nu\s+\[0 2 -1 0 0 0 0\]\s+)[^;]+;", rf"\g<1>{nu!r};")
    if a.coupling_tol_factor is not None:
        f = a.coupling_tol_factor
        for name in ("fsiProperties.robin", "fsiProperties.iqnils"):
            p = d / "constant" / name
            text = p.read_text()
            for key in ("outerCorrTolerance", "robinPressureTolerance", "robinFluxTolerance"):
                text = re.sub(rf"^(\s*{key}\s+)([0-9.eE+-]+);",
                              lambda m: f"{m.group(1)}{float(m.group(2)) * f:.6g};",
                              text, flags=re.MULTILINE)
            p.write_text(text)
        changes["couplingTolFactor"] = f
    (d / "diagnostic_settings.json").write_text(json.dumps(
        {"source": str(a.src), "changes": changes}, indent=2))
    print(f"prepared {d} with {changes or 'tutorial parameters'}")


if __name__ == "__main__":
    main()
