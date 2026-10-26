#!/usr/bin/env python3

"""Summarise incrementally imposed fluid mesh-motion checks."""

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


def metadata(path):
    result = {}
    for line in path.read_text().splitlines():
        if "=" in line:
            key, value = line.split("=", 1)
            result[key] = value
    return result


def analyse_checkpoint(case_label, path):
    information = metadata(path/"metadata")
    log_path = path/"log.checkMesh"
    initial_volume = path/"V0"
    moved_volume = path/"V"
    result = {
        "label": case_label,
        "checkpoint": information["checkpoint"],
        "status": information["status"],
        "reached": float(information["reached_mm"]),
    }
    if not log_path.exists() or not moved_volume.exists():
        return result

    log = log_path.read_text()
    volumes_0 = internal_scalars(initial_volume)
    volumes_1 = internal_scalars(moved_volume)
    if len(volumes_0) != len(volumes_1):
        raise ValueError(f"{path}: volume fields have different lengths")
    result.update({
        "status": (
            information["status"]
            if "Mesh OK." in log
            else f"{information['status']}; checkMesh failed"
        ),
        "ratio": np.min(volumes_1/volumes_0),
        "volume": log_value(
            log,
            r"Min volume = ([-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?)",
        ),
        "non_orthogonality": log_value(
            log, r"Mesh non-orthogonality Max: ([^ ]+)"
        ),
        "skewness": log_value(log, r"Max skewness = ([^ ]+)"),
    })
    return result


def analyse_case(label, case):
    checkpoints = []
    metrics = case/"metrics"
    if metrics.exists():
        for path in sorted(metrics.iterdir()):
            if path.is_dir() and (path/"metadata").exists():
                checkpoints.append(analyse_checkpoint(label, path))
    status_path = case/"run-status"
    status = metadata(status_path) if status_path.exists() else {
        "status": "not run",
        "largest_positive_mm": "0",
        "largest_negative_mm": "0",
    }
    status["label"] = label
    return checkpoints, status


def markdown(results, summaries):
    lines = [
        "| Case | Checkpoint | Reached (mm) | Motion/check | "
        "Min volume ratio | Min volume (m^3) | Max non-orthogonality "
        "(deg) | Max skewness |",
        "| :-- | :--: | --: | :--: | --: | --: | --: | --: |",
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
            f"| {result['label']} | {result['checkpoint']} | "
            f"{result['reached']:.3f} | {result['status']} | "
            + " | ".join(values) + " |"
        )
    lines.extend([
        "",
        "| Case | Run status | Largest positive (mm) | "
        "Largest negative magnitude (mm) |",
        "| :-- | :-- | --: | --: |",
    ])
    for summary in summaries:
        lines.append(
            f"| {summary['label']} | {summary['status']} | "
            f"{float(summary['largest_positive_mm']):.3f} | "
            f"{abs(float(summary['largest_negative_mm'])):.3f} |"
        )

    old_summaries = [
        summary for summary in summaries if summary["label"].startswith("old")
    ]
    representative = len(old_summaries) == 2 and all(
        float(summary["largest_positive_mm"]) >= 6.2
        for summary in old_summaries
    )
    lines.extend([
        "",
        "Old-mesh representativeness check: "
        + (
            "**PASS** (both old-mesh mode ramps reached at least 6.2 mm)."
            if representative
            else "**FAIL** (an old-mesh mode ramp did not reach 6.2 mm)."
        ),
    ])
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

    analyses = [
        analyse_case(label, Path(directory)) for label, directory in args.case
    ]
    table = markdown(
        [result for checkpoints, _ in analyses for result in checkpoints],
        [summary for _, summary in analyses],
    )
    print(table, end="")
    if args.output:
        args.output.write_text(table)


if __name__ == "__main__":
    main()
