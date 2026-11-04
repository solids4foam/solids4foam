"""Run and reuse single womersleyTube cases of a study copy.

Used by run_single.py. A run is a solver run only: it is not analysed here
(analysis is done afterwards by the analysis scripts, which can use their own
checks), so an analysis failure of one case cannot stop later runs.

A completed run is reused only if everything that defines it matches what is
stored in its verification_settings.json: the case specification, mesh,
reference parameters, content hashes of the tutorial files of the study copy
(this includes the amplitude, viscosity, coupling tolerances and the diagnostic)
and the solver build signature. It must also have a complete solver log and a
non-empty diagnostic output. A run that does not match is not deleted unless it
is about to be rerun.
"""
import json
from pathlib import Path

DIAGNOSTIC_FILES = ("fluid.dat", "solid.dat")


def reuse_problem(V, case: Path, wanted: dict):
    """Why the existing case cannot be reused, or None if it can."""
    if not case.is_dir():
        return "no run"
    stored = case / V.SETTINGS_FILE
    if not stored.is_file():
        return "no stored settings"
    try:
        previous = json.loads(stored.read_text())
    except json.JSONDecodeError:
        return "unreadable stored settings"
    if previous != wanted:
        keys = sorted(k for k in set(previous) | set(wanted)
                      if previous.get(k) != wanted.get(k))
        return "settings differ (" + ", ".join(keys) + ")"
    problem = V.solver_log_problem(case)
    if problem:
        return problem
    for name in DIAGNOSTIC_FILES:
        f = case / "postProcessing/interfaceTraction" / name
        if not f.is_file() or f.stat().st_size == 0:
            return f"missing or empty diagnostic output {name}"
    return None


def run_one(V, reference: dict, spec: dict) -> str:
    """Run the case unless a valid completed run exists."""
    case = V.WORK_ROOT / spec["name"]
    problem = reuse_problem(V, case, V.settings(spec, reference))
    if problem is None:
        return f"{spec['name']}: reusing the completed run"
    if case.exists():
        print(f"{spec['name']}: running again ({problem})", flush=True)
    return V.run_case(spec, reference, reuse=False)


def run_all(V, reference: dict, specs: list):
    """Run every spec; a failure of one does not stop the others. Returns
    the list of (name, error) of the failed ones."""
    failed = []
    for spec in specs:
        try:
            print(run_one(V, reference, spec), flush=True)
        except Exception as error:   # noqa: BLE001 - reported, not hidden
            print(f"{spec['name']}: FAILED: {error}", flush=True)
            failed.append((spec["name"], str(error)))
    return failed
