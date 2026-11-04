#!/usr/bin/env python3
"""Check that the solver runs a study needs are valid: each has stored
settings equal to what the study copy's driver computes now (case, mesh,
parameters, tutorial files, configured inputs, solver build), a completed
solver log and non-empty diagnostic output. A non-empty fluid.dat alone is not
enough.

usage: check_runs.py WORK_DIR [--skip-build-check] [--keys-out FILE]
                     study:case [study:case ...]

  WORK_DIR/<study>/womersleyTube is a prepared study copy.
  --skip-build-check  do not compare the solver build (analysing stored runs
                      on a machine without the solver environment).
  --keys-out FILE     write the analysis key and the digests of the solver-run
                      keys of the checked runs (provenance of the results).
Exit status 1 if any run is not valid.
"""
import argparse
import importlib.util
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import study_runs  # noqa: E402


def load_driver(study: Path):
    path = study / "verification/scripts/womersley_tube_verification.py"
    spec = importlib.util.spec_from_file_location(f"wtv_{study.parent.name}", path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("work", type=Path)
    ap.add_argument("runs", nargs="+", help="study:case")
    ap.add_argument("--skip-build-check", action="store_true")
    ap.add_argument("--keys-out", type=Path)
    args = ap.parse_args(argv)
    ignore = ("build",) if args.skip_build_check else ()
    by_study = {}
    for item in args.runs:
        study, case = item.split(":")
        by_study.setdefault(study, []).append(case)
    problems, keys, files = [], {}, {}
    for study, names in by_study.items():
        root = args.work / study / "womersleyTube"
        V = load_driver(root)
        reference = json.loads(V.REFERENCE_FILE.read_text())
        for name, problem in study_runs.check_required(V, reference, names, ignore):
            problems.append(f"{study}/{name}: {problem}")
        for name in names:
            keys[f"{study}/{name}"] = study_runs.solver_run_key(
                V, study_runs.spec_from_name(V, name), reference)
        for f in (root / "verification/scripts").glob("*.py"):
            files[f"{study}/verification/scripts/{f.name}"] = f
        files[f"{study}/reference/{V.REFERENCE_FILE.name}"] = V.REFERENCE_FILE
    for f in HERE.glob("*.py"):
        files[f"traction/scripts/{f.name}"] = f
    if problems:
        print(f"ERROR: {len(problems)} required run(s) missing or not valid:",
              file=sys.stderr)
        for line in problems:
            print("  " + line, file=sys.stderr)
        return 1
    if args.keys_out:
        record = {"solverRunKeys": keys,
                  "analysisKey": study_runs.analysis_key(files),
                  "buildChecked": not args.skip_build_check}
        args.keys_out.write_text(json.dumps(record, indent=1, sort_keys=True))
    return 0


if __name__ == "__main__":
    sys.exit(main())
