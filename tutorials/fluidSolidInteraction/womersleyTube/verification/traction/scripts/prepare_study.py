#!/usr/bin/env python3
"""Prepare a womersleyTube study copy with the interface-traction diagnostic.

usage: prepare_study.py SRC_TUTORIAL DEST [--mu MU] [--amplitude A]
                        [--coupling-tol-factor F] [--periods N]

- copies the tutorial (without results) to DEST;
- installs scripts/interfaceTractionDiagnostic as system/interfaceTractionDiagnostic
  of the copy and includes it in the controlDict functions,
  sampling every time step over the last two periods (timeStart = 4 periods);
- --mu / --amplitude: change the dynamic viscosity (Pa s) or the pressure
  amplitude (Pa) and REGENERATE, consistently, constant/womersleyProperties
  (end data and initial fields) and the reference JSON with the copy's own
  womersley_exact.py; constant/fluid/transportProperties nu is set to mu/rho;
- --coupling-tol-factor: multiply the Robin and IQN-ILS coupling tolerances;
- --periods N: run N periods instead of the study's 6 (set in the copy's
  reference JSON, which the verification driver reads for the end time); the
  analysis window stays the last analysisPeriods (2), and the diagnostic
  samples them (timeStart = (N - 2) periods). A negative --amplitude is valid
  (use --amplitude=-1).
Nothing in SRC is modified.

The copy is built in a staging directory DEST.partial (any leftover one is
discarded), every copied file and every substitution is checked, the marker
diagnostic_settings.json is written last, and the staging directory is then
published with renames only: an existing DEST is first renamed aside to
DEST.old-<time> and deleted only after the new one is in place. A crash can
therefore leave DEST complete (old or new) or, for an instant between the two
renames, absent with DEST.old-* intact, which --ensure puts back. The marker
records two hashes: the solver hash (everything the solver reads, the
diagnostic and this script) and the analysis hash (verification/scripts/*.py
and the reference JSON). With --ensure, an existing DEST is kept if both match;
if only the analysis hash differs it is rebuilt and its runs moved into the new
copy (whether a run is still valid is decided by its solver-run key, see
study_runs.py); a copy that differs in the solver hash or lacks a complete
marker is re-prepared if it has no runs, and otherwise refused (exit status 3)
unless --force, so valid runs are never deleted silently. --overwrite-diagnostic
allows replacing a different system/interfaceTractionDiagnostic in SRC.
"""
import argparse
import hashlib
import json
import os
import re
import shutil
import sys
import time
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


def is_analysis_file(rel: Path) -> bool:
    """Files that only determine how results are analysed or checked."""
    return rel.parts[:2] == ("verification", "scripts") and rel.suffix == ".py" \
        or rel.parts[:2] == ("verification", "reference")


def source_hashes(src: Path):
    """(solver hash, analysis hash). The solver hash covers the files copied
    from SRC (path and content, symlink targets; not documentation, images,
    verification/traction, or the analysis files) and the diagnostic and this
    script. The analysis hash covers verification/scripts/*.py and the
    reference JSON."""
    solver, analysis = hashlib.sha256(), hashlib.sha256()
    ignore = shutil.ignore_patterns(*IGNORE)

    def walk(directory: Path, rel: Path):
        names = sorted(p.name for p in directory.iterdir())
        skipped = ignore(str(directory), names)
        for name in names:
            if name in skipped:
                continue
            path = directory / name
            r = rel / name
            if r.suffix == ".md" or r.parts[:1] == ("images",) \
                    or r.parts[:2] == ("verification", "traction"):
                continue    # documentation and this diagnostic's own files
            h = analysis if is_analysis_file(r) else solver
            if path.is_symlink():
                h.update(f"L {r} {path.readlink()}\n".encode())
            elif path.is_dir():
                h.update(f"D {r}\n".encode())
                walk(path, r)
            else:
                h.update(f"F {r} {hashlib.sha256(path.read_bytes()).hexdigest()}\n".encode())

    walk(src, Path("."))
    for extra in (DIAGNOSTIC, Path(__file__).resolve()):
        solver.update(f"X {extra.name} {hashlib.sha256(extra.read_bytes()).hexdigest()}\n".encode())
    return solver.hexdigest(), analysis.hexdigest()


def expected_paths(src: Path):
    """Relative paths of the files and symlinks copytree must have made."""
    ignore = shutil.ignore_patterns(*IGNORE)
    found = set()

    def walk(directory: Path, rel: Path):
        names = [p.name for p in directory.iterdir()]
        skipped = ignore(str(directory), names)
        for name in names:
            if name in skipped:
                continue
            path = directory / name
            if path.is_dir() and not path.is_symlink():
                walk(path, rel / name)
            else:
                found.add(rel / name)

    walk(src, Path("."))
    return found


def validate_copy(src: Path, d: Path, requested: dict):
    """Every file of SRC is in the copy, and the preparation's edits are in
    place; raises SystemExit otherwise."""
    missing = sorted(str(r) for r in expected_paths(src)
                     if not (d / r).exists() and not (d / r).is_symlink())
    if missing:
        raise SystemExit(f"incomplete copy, missing: {', '.join(missing[:5])}"
                         + (" ..." if len(missing) > 5 else ""))
    control = (d / "system/controlDict").read_text()
    if control.count('#include "interfaceTractionDiagnostic"') != 1:
        raise SystemExit("system/controlDict does not include the diagnostic exactly once")
    diag = (d / "system/interfaceTractionDiagnostic").read_text()
    for key in ("executeInterval", "writeInterval"):
        if len(re.findall(rf"^\s*{key}\s+1;", diag, re.MULTILINE)) != 1:
            raise SystemExit(f"diagnostic: {key} is not 1")
    if requested.get("amplitude") is not None or requested.get("mu") is not None:
        if not (d / "constant/womersleyProperties").is_file():
            raise SystemExit("constant/womersleyProperties was not regenerated")


def has_runs(dest: Path) -> bool:
    work = dest / "verification/work"
    return work.is_dir() and any(work.iterdir())


def ensure_state(dest: Path, wanted: dict, analysis_hash: str) -> str:
    """'missing', 'current' (complete and matching the requested solver hash
    and settings, and the analysis hash), 'refresh-analysis' (complete, solver
    side matching, analysis files changed: rebuilt, runs moved over),
    'stale-empty' (to be re-prepared; no runs would be lost) or 'stale-runs'
    (has runs: only with --force)."""
    if not dest.exists():
        return "missing"
    marker = dest / MARKER
    try:
        record = json.loads(marker.read_text())
    except (OSError, json.JSONDecodeError):
        record = {}
    if record.get("complete") is True and all(record.get(k) == v for k, v in wanted.items()):
        return "current" if record.get("analysisHash") == analysis_hash else "refresh-analysis"
    return "stale-runs" if has_runs(dest) else "stale-empty"


def recover_interrupted(dest: Path):
    """If a preparation was interrupted between moving DEST aside and
    publishing the new copy, DEST is absent and DEST.old-* holds the previous
    copy: put the newest back."""
    if dest.exists():
        return
    olds = sorted(dest.parent.glob(dest.name + ".old-*"))
    if olds:
        os.rename(olds[-1], dest)
        print(f"{dest}: restored the previous copy from {olds[-1].name} "
              "(an earlier preparation was interrupted)")


def publish(staging: Path, dest: Path, keep_runs: bool):
    """Replace DEST by the complete staging directory using renames only. The
    previous DEST is renamed aside and deleted after success. With keep_runs
    its verification/work is moved into the new copy first."""
    old = None
    if dest.exists():
        old = dest.with_name(f"{dest.name}.old-{time.strftime('%Y%m%dT%H%M%S')}-{os.getpid()}")
        os.rename(dest, old)
    try:
        os.replace(staging, dest)
    except BaseException:
        if old is not None and not dest.exists():
            os.rename(old, dest)
        raise
    if old is not None:
        runs = old / "verification/work"
        if keep_runs and runs.is_dir():
            os.rename(runs, dest / "verification/work")
        shutil.rmtree(old)


def prepare(src: Path, dest: Path, mu=None, amplitude=None, coupling_tol_factor=None,
            overwrite_diagnostic=False, keep_runs=False, periods=None):
    src, dest = src.resolve(), dest.resolve()
    staging = dest.with_name(dest.name + ".partial")
    if staging.exists():
        shutil.rmtree(staging)  # an unfinished earlier preparation, holds no runs
    try:
        build(src, staging, mu, amplitude, coupling_tol_factor, overwrite_diagnostic, periods)
        publish(staging, dest, keep_runs)
    except BaseException:
        shutil.rmtree(staging, ignore_errors=True)
        raise
    changes = json.loads((dest / MARKER).read_text())["changes"]
    print(f"prepared {dest} with {changes or 'tutorial parameters'}")


def build(src: Path, work: Path, mu, amplitude, coupling_tol_factor, overwrite_diagnostic,
          periods=None):
    """Build the study copy in the staging directory `work`; the marker is
    written last, after the copy and the edits have been validated."""
    solver_hash, analysis_hash = source_hashes(src)
    wanted = {"solverHash": solver_hash,
              "requested": requested_settings(mu, amplitude, coupling_tol_factor, periods)}
    shutil.copytree(src, work, symlinks=True, ignore=shutil.ignore_patterns(*IGNORE))
    d = work
    # diagnostic function object: always the one next to this script. A
    # different diagnostic already in SRC is not replaced silently.
    installed = d / "system/interfaceTractionDiagnostic"
    if installed.exists() and installed.read_bytes() != DIAGNOSTIC.read_bytes() \
            and not overwrite_diagnostic:
        raise SystemExit(f"{src}/system/interfaceTractionDiagnostic differs from "
                         f"{DIAGNOSTIC}; use --overwrite-diagnostic to replace it")
    shutil.copy(DIAGNOSTIC, d / "system")
    fo = installed
    sub_one(fo, r"^(\s*executeInterval\s+)\d+;", r"\g<1>1;")
    sub_one(fo, r"^(\s*writeInterval\s+)\d+;", r"\g<1>1;")
    ref_file = d / "verification/reference/womersleyTube_verification_references.json"
    ref = json.loads(ref_file.read_text())
    period = ref["reference"]["derived"]["period"]
    window = ref["study"]["analysisPeriods"]
    changes = {}
    if periods is not None:
        if periods <= window:
            raise SystemExit(f"--periods {periods} must exceed the {window} analysed periods")
        ref["study"]["periods"] = periods
        ref_file.write_text(json.dumps(ref, indent=2))
        changes["periods"] = periods
    periods = ref["study"]["periods"]
    sub_one(fo, r"^(\s*timeStart\s+)[^;]+;", rf"\g<1>{(periods - window) * period:.12g};")
    sub_one(d / "system/controlDict", r"^functions\s*\n\{",
            'functions\n{\n    #include "interfaceTractionDiagnostic"\n')
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
                if key not in ("dynamicViscosity", "pressureAmplitude"):
                    continue    # e.g. periods: not a parameter of the exact script
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
    validate_copy(src, d, wanted["requested"])
    if periods is not None and json.loads(ref_file.read_text())["study"]["periods"] != periods:
        raise SystemExit("the reference JSON of the copy does not hold the requested periods")
    (d / MARKER).write_text(json.dumps(
        dict(wanted, analysisHash=analysis_hash, complete=True, source=str(src),
             changes=changes), indent=2))


def requested_settings(mu, amplitude, coupling_tol_factor, periods=None):
    """The requested settings recorded in the marker. periods is recorded only
    when requested, so that copies prepared without it keep their marker."""
    requested = {"mu": mu, "amplitude": amplitude, "couplingTolFactor": coupling_tol_factor}
    if periods is not None:
        requested["periods"] = periods
    return requested


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("src", type=Path)
    ap.add_argument("dest", type=Path)
    ap.add_argument("--mu", type=float)
    ap.add_argument("--amplitude", type=float)
    ap.add_argument("--coupling-tol-factor", type=float)
    ap.add_argument("--periods", type=int)
    ap.add_argument("--ensure", action="store_true",
                    help="keep DEST if complete and matching SRC and the settings")
    ap.add_argument("--force", action="store_true",
                    help="with --ensure, re-prepare a stale DEST even if it holds runs")
    ap.add_argument("--overwrite-diagnostic", action="store_true")
    a = ap.parse_args(argv)
    kw = dict(mu=a.mu, amplitude=a.amplitude, coupling_tol_factor=a.coupling_tol_factor,
              overwrite_diagnostic=a.overwrite_diagnostic, periods=a.periods)
    if not a.ensure:
        if a.dest.exists():
            raise SystemExit(f"{a.dest} exists")
        return prepare(a.src, a.dest, **kw)
    recover_interrupted(a.dest.resolve())
    solver_hash, analysis_hash = source_hashes(a.src.resolve())
    wanted = {"solverHash": solver_hash,
              "requested": requested_settings(a.mu, a.amplitude, a.coupling_tol_factor,
                                              a.periods)}
    state = ensure_state(a.dest, wanted, analysis_hash)
    if state == "current":
        print(f"{a.dest}: up to date")
        return
    if state == "stale-runs" and not a.force:
        print(f"ERROR: {a.dest} is incomplete or differs from the tutorial/scripts/settings "
              "and holds runs; they are kept. Re-run with --force to discard it.",
              file=sys.stderr)
        raise SystemExit(3)
    if state == "refresh-analysis":
        print(f"{a.dest}: the analysis files changed, preparing again and keeping the runs "
              "(a run is reused only if its solver-run key still matches)")
    elif state.startswith("stale"):
        print(f"{a.dest}: stale or incomplete, preparing again")
    prepare(a.src, a.dest, keep_runs=(state == "refresh-analysis"), **kw)


if __name__ == "__main__":
    main()
