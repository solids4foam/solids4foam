#!/usr/bin/env python3
"""Prepare a womersleyTube study copy with the interface-traction diagnostic.

usage: prepare_study.py SRC_TUTORIAL DEST [--mu MU] [--amplitude A]
                        [--coupling-tol-factor F]

- copies the tutorial (without results) to DEST;
- installs scripts/interfaceTractionDiagnostic as system/interfaceTractionDiagnostic
  of the copy and includes it in the controlDict functions,
  sampling every time step over the last two periods (timeStart = 4 periods);
- --mu / --amplitude: change the dynamic viscosity (Pa s) or the pressure
  amplitude (Pa) and REGENERATE, consistently, constant/womersleyProperties
  (end data and initial fields) and the reference JSON with the copy's own
  womersley_exact.py; constant/fluid/transportProperties nu is set to mu/rho;
- --coupling-tol-factor: multiply the Robin and IQN-ILS coupling tolerances.
Nothing in SRC is modified.

The copy is built under DEST.partial and renamed to DEST only when complete, and
diagnostic_settings.json (written last) records a hash of everything copied
and the requested settings. With --ensure, an existing DEST is kept only if
that record matches; a stale or unrecorded DEST is re-prepared if it has no
runs, and otherwise refused (exit status 3) unless --force, so valid runs are
never deleted silently. --overwrite-diagnostic allows replacing a different
system/interfaceTractionDiagnostic in SRC.
"""
import argparse
import hashlib
import json
import re
import shutil
import sys
from pathlib import Path

DIAGNOSTIC = Path(__file__).resolve().parent / "interfaceTractionDiagnostic"
IGNORE = ("postProcessing", "work", "dynamicCode", "log.*", "processor*", "__pycache__")
MARKER = "diagnostic_settings.json"
# Entries scaled by --coupling-tol-factor: each must be found exactly once
TOLERANCE_KEYS = {
    "fsiProperties.robin": ("outerCorrTolerance", "robinPressureTolerance",
                            "robinFluxTolerance"),
    "fsiProperties.iqnils": ("outerCorrTolerance",),
}
NUMBER = r"[0-9]*\.?[0-9]+(?:[eE][+-]?[0-9]+)?"


def sub_one(path: Path, pattern: str, repl: str, count_expected=1):
    text = path.read_text()
    new, n = re.subn(pattern, repl, text, flags=re.MULTILINE)
    if n != count_expected:
        raise SystemExit(f"{path}: expected {count_expected} match(es) of {pattern!r}, found {n}")
    path.write_text(new)


def scale_tolerances(text: str, keys, factor: float, label: str = "file") -> str:
    """Multiply the numeric value of each key by factor. Every key must be
    matched exactly once, otherwise the study would silently run with the
    baseline tolerance."""
    for key in keys:
        text, n = re.subn(
            rf"^(\s*{key}\s+)({NUMBER})(\s*;)",
            lambda m: f"{m.group(1)}{float(m.group(2)) * factor:.6g}{m.group(3)}",
            text, flags=re.MULTILINE)
        if n != 1:
            raise SystemExit(f"{label}: expected exactly one entry '{key} <number>;' "
                             f"to scale, found {n}")
    return text


def source_hash(src: Path) -> str:
    """Hash of everything a copy is made from: the files copied from SRC (path
    and content, symlink targets) and the diagnostic and this script."""
    h = hashlib.sha256()
    ignore = shutil.ignore_patterns(*IGNORE)

    def walk(directory: Path, rel: Path):
        names = sorted(p.name for p in directory.iterdir())
        skipped = ignore(str(directory), names)
        for name in names:
            if name in skipped:
                continue
            path = directory / name
            r = rel / name
            if path.is_symlink():
                h.update(f"L {r} {path.readlink()}\n".encode())
            elif path.is_dir():
                h.update(f"D {r}\n".encode())
                walk(path, r)
            else:
                h.update(f"F {r} {hashlib.sha256(path.read_bytes()).hexdigest()}\n".encode())

    walk(src, Path("."))
    for extra in (DIAGNOSTIC, Path(__file__).resolve()):
        h.update(f"X {extra.name} {hashlib.sha256(extra.read_bytes()).hexdigest()}\n".encode())
    return h.hexdigest()


def has_runs(dest: Path) -> bool:
    work = dest / "verification/work"
    return work.is_dir() and any(work.iterdir())


def ensure_state(dest: Path, wanted: dict) -> str:
    """'missing', 'current' (complete and matching the requested hash and
    settings), 'stale-empty' (to be re-prepared; no runs would be lost) or
    'stale-runs' (has runs: only with --force)."""
    if not dest.exists():
        return "missing"
    marker = dest / MARKER
    try:
        record = json.loads(marker.read_text())
    except (OSError, json.JSONDecodeError):
        record = {}
    if record.get("complete") is True and all(record.get(k) == v for k, v in wanted.items()):
        return "current"
    return "stale-runs" if has_runs(dest) else "stale-empty"


def prepare(src: Path, dest: Path, mu=None, amplitude=None, coupling_tol_factor=None,
            overwrite_diagnostic=False):
    src, dest = src.resolve(), dest.resolve()
    wanted = {"sourceHash": source_hash(src),
              "requested": {"mu": mu, "amplitude": amplitude,
                            "couplingTolFactor": coupling_tol_factor}}
    work = dest.with_name(dest.name + ".partial")
    if work.exists():
        shutil.rmtree(work)     # an unfinished earlier preparation, holds no runs
    shutil.copytree(src, work, symlinks=True, ignore=shutil.ignore_patterns(*IGNORE))
    d = work
    # diagnostic function object: always the one next to this script. A
    # different diagnostic already in SRC is not replaced silently.
    installed = d / "system/interfaceTractionDiagnostic"
    if installed.exists() and installed.read_bytes() != DIAGNOSTIC.read_bytes() \
            and not overwrite_diagnostic:
        shutil.rmtree(work)
        raise SystemExit(f"{src}/system/interfaceTractionDiagnostic differs from "
                         f"{DIAGNOSTIC}; use --overwrite-diagnostic to replace it")
    shutil.copy(DIAGNOSTIC, d / "system")
    fo = installed
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
    if mu is not None or amplitude is not None:
        sys.path.insert(0, str(d / "verification/scripts"))
        try:
            import womersley_exact as W
            if mu is not None:
                W.PARAMETERS["dynamicViscosity"] = mu
                changes["dynamicViscosity"] = mu
            if amplitude is not None:
                W.PARAMETERS["pressureAmplitude"] = amplitude
                changes["pressureAmplitude"] = amplitude
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
            if mu is not None:
                nu = mu / W.PARAMETERS["fluidDensity"]
                sub_one(d / "constant/fluid/transportProperties",
                        r"^(nu\s+nu\s+\[0 2 -1 0 0 0 0\]\s+)[^;]+;", rf"\g<1>{nu!r};")
        finally:
            sys.path.remove(str(d / "verification/scripts"))
            sys.modules.pop("womersley_exact", None)
    if coupling_tol_factor is not None:
        for name, keys in TOLERANCE_KEYS.items():
            p = d / "constant" / name
            p.write_text(scale_tolerances(p.read_text(), keys, coupling_tol_factor,
                                          str(p.relative_to(d))))
        changes["couplingTolFactor"] = coupling_tol_factor
    (d / MARKER).write_text(json.dumps(
        dict(wanted, complete=True, source=str(src), changes=changes), indent=2))
    if dest.exists():
        shutil.rmtree(dest)
    work.rename(dest)
    print(f"prepared {dest} with {changes or 'tutorial parameters'}")


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("src", type=Path)
    ap.add_argument("dest", type=Path)
    ap.add_argument("--mu", type=float)
    ap.add_argument("--amplitude", type=float)
    ap.add_argument("--coupling-tol-factor", type=float)
    ap.add_argument("--ensure", action="store_true",
                    help="keep DEST if complete and matching SRC and the settings")
    ap.add_argument("--force", action="store_true",
                    help="with --ensure, re-prepare a stale DEST even if it holds runs")
    ap.add_argument("--overwrite-diagnostic", action="store_true")
    a = ap.parse_args(argv)
    kw = dict(mu=a.mu, amplitude=a.amplitude, coupling_tol_factor=a.coupling_tol_factor,
              overwrite_diagnostic=a.overwrite_diagnostic)
    if not a.ensure:
        if a.dest.exists():
            raise SystemExit(f"{a.dest} exists")
        return prepare(a.src, a.dest, **kw)
    wanted = {"sourceHash": source_hash(a.src.resolve()),
              "requested": {"mu": a.mu, "amplitude": a.amplitude,
                            "couplingTolFactor": a.coupling_tol_factor}}
    state = ensure_state(a.dest, wanted)
    if state == "current":
        print(f"{a.dest}: up to date")
        return
    if state == "stale-runs" and not a.force:
        print(f"ERROR: {a.dest} is incomplete or differs from the tutorial/scripts/settings "
              "and holds runs; they are kept. Re-run with --force to discard it.",
              file=sys.stderr)
        raise SystemExit(3)
    if state.startswith("stale"):
        print(f"{a.dest}: stale or incomplete, preparing again")
    prepare(a.src, a.dest, **kw)


if __name__ == "__main__":
    main()
