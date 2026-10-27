#!/usr/bin/env python3
"""Run single womersleyTube cases of a prepared study copy through the
study's own verification driver (same case configuration, mesh scaling,
time step and checks), e.g. for targeted variants.

usage: run_single.py STUDY_DIR coupling:factor:steps [...]
"""
import sys
from pathlib import Path

study = Path(sys.argv[1]).resolve()
sys.path.insert(0, str(study / "verification/scripts"))
import json  # noqa: E402
import womersley_tube_verification as V  # noqa: E402

reference = json.loads(V.REFERENCE_FILE.read_text())
V.WORK_ROOT.mkdir(parents=True, exist_ok=True)
for item in sys.argv[2:]:
    coupling, factor, steps = item.split(":")
    spec = {"coupling": coupling, "factor": int(factor), "steps": int(steps)}
    spec["name"] = V.case_name(spec)
    print(V.run_case(spec, reference, reuse=False), flush=True)
    row = V.analyse(dict(spec, groups=["single"]), reference)
    print({k: v for k, v in row.items() if not isinstance(v, (list, dict))}, flush=True)
