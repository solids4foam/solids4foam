#!/usr/bin/env python3

"""Summarise imposed-displacement fluid mesh checks."""

import argparse
import re
from pathlib import Path

import numpy as np


def internal_scalars(path):
    text = path.read_text()
    match = re.search(
        r"internalField\s+nonuniform\s+List<scalar>\s+\d+\s*\((.*?)\)\s*;",
        text,
        re.DOTALL,
    )
    if not match:
        raise ValueError(f"{path}: cannot find nonuniform internal field")
    return np.fromstring(match.group(1), sep=" ")


def log_value(text, pattern):
    match = re.search(pattern, text)
    return float(match.group(1)) if match else np.nan


def analyse(label, case):
    log_path = case/"log.checkMesh.moved"
    initial_volume = case/"0"/"fluid"/"V"
    moved_volume = case/"1"/"fluid"/"V"
    if not log_path.exists() or not moved_volume.exists():
        return {"label": label, "status": "motion failed"}

    log = log_path.read_text()
    volumes_0 = internal_scalars(initial_volume)
    volumes_1 = internal_scalars(moved_volume)
    if len(volumes_0) != len(volumes_1):
        raise ValueError(f"{case}: volume fields have different lengths")
    return {
        "label": label,
        "status": "OK" if "Mesh OK." in log else "checkMesh failed",
        "ratio": np.min(volumes_1/volumes_0),
        "volume": log_value(
            log,
            r"Min volume = ([-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?)",
        ),
        "non_orthogonality": log_value(
            log, r"Mesh non-orthogonality Max: ([^ ]+)"
        ),
        "skewness": log_value(log, r"Max skewness = ([^ ]+)"),
    }


def markdown(results):
    lines = [
        "| Case | Motion/check | Min volume ratio | Min volume (m^3) | "
        "Max non-orthogonality (deg) | Max skewness |",
        "| :-- | :--: | --: | --: | --: | --: |",
    ]
    for result in results:
        if "ratio" not in result:
            values = ["-", "-", "-", "-"]
        else:
            values = [
                f"{result['ratio']:.6g}",
                f"{result['volume']:.6g}",
                f"{result['non_orthogonality']:.4g}",
                f"{result['skewness']:.4g}",
            ]
        lines.append(
            f"| {result['label']} | {result['status']} | "
            + " | ".join(values) + " |"
        )
    return "\n".join(lines) + "\n"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--case",
        action="append",
        nargs=2,
        metavar=("LABEL", "DIRECTORY"),
        required=True,
    )
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    table = markdown([
        analyse(label, Path(directory)) for label, directory in args.case
    ])
    print(table, end="")
    if args.output:
        args.output.write_text(table)


if __name__ == "__main__":
    main()
