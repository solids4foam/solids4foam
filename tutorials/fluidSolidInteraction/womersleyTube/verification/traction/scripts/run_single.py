#!/usr/bin/env python3
"""Run single womersleyTube cases of a prepared study copy through the
study's own verification driver (same case configuration, mesh scaling,
time step and solver checks), e.g. for targeted variants.

usage: run_single.py STUDY_DIR coupling:factor:steps [...]

Only the solver run is made here; the analysis scripts analyse the results.
A completed run whose settings, tutorial files and solver build match is
reused, never deleted. A failed case does not prevent the remaining cases from
running; the exit status is 1 if any case failed.
"""
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import study_runs  # noqa: E402


def parse_specs(V, items):
    specs = []
    for item in items:
        coupling, factor, steps = item.split(":")
        spec = {"coupling": coupling, "factor": int(factor), "steps": int(steps)}
        spec["name"] = V.case_name(spec)
        specs.append(spec)
    return specs


def main(argv) -> int:
    study = Path(argv[1]).resolve()
    sys.path.insert(0, str(study / "verification/scripts"))
    import womersley_tube_verification as V
    reference = json.loads(V.REFERENCE_FILE.read_text())
    V.WORK_ROOT.mkdir(parents=True, exist_ok=True)
    failed = study_runs.run_all(V, reference, parse_specs(V, argv[2:]))
    if failed:
        print(f"ERROR: {len(failed)} case(s) failed: "
              + ", ".join(n for n, _ in failed), file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
