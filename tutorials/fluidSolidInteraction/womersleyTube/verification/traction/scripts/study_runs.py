"""Run and reuse single womersleyTube cases of a study copy.

Used by run_single.py. A run is a solver run only: it is not analysed here
(analysis is done afterwards by the analysis scripts, which can use their own
checks), so an analysis failure of one case cannot stop later runs.

Two keys are kept apart:

- the SOLVER-RUN key is the driver's settings(): the case specification, mesh,
  physical parameters, time schedule, content hashes of the tutorial files of
  the study copy (amplitude, viscosity, coupling tolerances, the diagnostic),
  a digest of the configured case inputs and the content hashes of the solver
  executable and libraries. A completed run is reused only if it equals what is
  stored in the run's verification_settings.json, the solver log shows a
  completed run and the diagnostic output is non-empty. A run that does not
  match is not deleted unless it is about to be rerun. Editing a README, the
  reference's exact values or any analysis code does not change this key.
- the ANALYSIS key (analysis_key) covers the python files that determine the
  numbers (the verification driver and exact solution, this directory's
  analysis scripts) and the reference JSON. It is recorded with the results
  (check_runs.py --keys-out) so that results can be told stale; it never
  triggers a solver rerun.
"""
import hashlib
import json
import re
from pathlib import Path

DIAGNOSTIC_FILES = ("fluid.dat", "solid.dat")


def reuse_problem(V, case: Path, wanted: dict, ignore_keys=()):
    """Why the existing case cannot be reused, or None if it can. Keys named
    in ignore_keys (e.g. "build") are not compared."""
    if not case.is_dir():
        return "no run"
    stored = case / V.SETTINGS_FILE
    if not stored.is_file():
        return "no stored settings"
    try:
        previous = json.loads(stored.read_text())
    except json.JSONDecodeError:
        return "unreadable stored settings"
    keys = sorted(k for k in set(previous) | set(wanted)
                  if k not in ignore_keys and previous.get(k) != wanted.get(k))
    if keys:
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


def spec_from_name(V, name: str) -> dict:
    """The spec of a case name such as robin_m2_n200."""
    match = re.fullmatch(r"([A-Za-z0-9]+)_m(\d+)_n(\d+)", name)
    if not match:
        raise ValueError(f"not a case name: {name}")
    spec = {"coupling": match.group(1), "factor": int(match.group(2)),
            "steps": int(match.group(3))}
    spec["name"] = V.case_name(spec)
    return spec


def check_required(V, reference: dict, names, ignore_keys=()):
    """Problems of the required runs of one study copy, as (name, problem)."""
    problems = []
    for name in names:
        wanted = V.settings(spec_from_name(V, name), reference)
        problem = reuse_problem(V, V.WORK_ROOT / name, wanted, ignore_keys)
        if problem:
            problems.append((name, problem))
    return problems


def solver_run_key(V, spec: dict, reference: dict) -> str:
    """Digest of the solver-run key (the driver's settings)."""
    return hashlib.sha1(json.dumps(V.settings(spec, reference), sort_keys=True)
                        .encode()).hexdigest()


def analysis_key(files: dict) -> dict:
    """Content hash of every file that determines the analysed numbers;
    files maps a label (independent of where the study is) to a path."""
    return {label: hashlib.sha1(Path(path).read_bytes()).hexdigest()
            for label, path in sorted(files.items())}
