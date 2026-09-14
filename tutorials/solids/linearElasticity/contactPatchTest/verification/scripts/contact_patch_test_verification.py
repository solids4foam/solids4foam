#!/usr/bin/env python3
"""Run the contact patch test with its two legacy solid formulations."""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
import shutil
import subprocess
import sys
from pathlib import Path


VARIANTS = (
    "linearGeometryTotalDisplacement",
    "unsLinearGeometry",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--verification-dir", required=True, type=Path)
    parser.add_argument(
        "--variant",
        choices=("all", *VARIANTS),
        default="all",
        help="formulation to run (default: both)",
    )
    parser.add_argument(
        "--reuse",
        action="store_true",
        help="reuse a completed case instead of rerunning it",
    )
    return parser.parse_args()


def require_commands(commands: tuple[str, ...]) -> None:
    missing = [command for command in commands if shutil.which(command) is None]
    if missing:
        raise RuntimeError("missing required command(s): " + ", ".join(missing))


def prepare_case(parent: Path, case_dir: Path, variant: str) -> None:
    if case_dir.exists():
        shutil.rmtree(case_dir)

    shutil.copytree(
        parent,
        case_dir,
        ignore=shutil.ignore_patterns(
            "verification",
            "regressionTests",
            "[1-9]*",
            "log.*",
            "case.foam",
            "lnInclude",
        ),
    )

    properties_path = case_dir / "constant" / "solidProperties"
    properties = properties_path.read_text(encoding="utf-8")
    properties = re.sub(
        r"(?m)^solidModel\s+\w+\s*;",
        f"solidModel     {variant};",
        properties,
        count=1,
    )

    if variant == "unsLinearGeometry":
        properties += """

unsLinearGeometryCoeffs
{
    nCorrectors          10000;
    solutionTolerance    1e-06;
    alternativeTolerance 1e-06;
    materialTolerance    1e-04;
    infoFrequency        100;
}
"""

    properties_path.write_text(properties, encoding="utf-8")


def run_case(case_dir: Path) -> None:
    output_path = case_dir / "log.Allverify"
    with output_path.open("w", encoding="utf-8") as output:
        completed = subprocess.run(
            ["./Allrun"],
            cwd=case_dir,
            stdout=output,
            stderr=subprocess.STDOUT,
            check=False,
        )
    if completed.returncode:
        raise RuntimeError(
            f"case failed with exit code {completed.returncode}; see {output_path}"
        )


def final_match(pattern: str, text: str, description: str) -> float:
    matches = re.findall(pattern, text, flags=re.MULTILINE)
    if not matches:
        raise RuntimeError(f"could not extract {description} from solver log")
    value = float(matches[-1])
    if not math.isfinite(value):
        raise RuntimeError(f"non-finite {description} in solver log")
    return value


def read_result(case_dir: Path, variant: str, references: dict[str, float]) -> dict:
    log_path = case_dir / "log.solids4Foam"
    if not log_path.is_file():
        raise RuntimeError(f"missing solver log: {log_path}")
    log = log_path.read_text(encoding="utf-8", errors="replace")

    epsilon_eq = final_match(
        r"Max epsilonEq\s*=\s*([^\s]+)", log, "maximum equivalent strain"
    )
    sigma_eq = final_match(
        r"Max sigmaEq \(von Mises stress\)\s*=\s*([^\s]+)",
        log,
        "maximum equivalent stress",
    )
    sigma_y_error = final_match(
        r"Average relative error in sigma_y field:\s*([^%\s]+)%",
        log,
        "average vertical-stress error",
    )

    analytical_sigma = references["analytical_sigma_y_magnitude_pa"]
    sigma_eq_error = abs(sigma_eq - analytical_sigma) / analytical_sigma * 100.0
    converged = "The momentum equation converged in all time-steps" in log
    passed = (
        converged
        and sigma_y_error
        <= references["max_average_sigma_y_error_percent"]
        and sigma_eq_error
        <= references["max_sigma_eq_relative_error_percent"]
    )

    return {
        "variant": variant,
        "converged": converged,
        "max_epsilon_eq": epsilon_eq,
        "max_sigma_eq_pa": sigma_eq,
        "average_sigma_y_error_percent": sigma_y_error,
        "sigma_eq_relative_error_percent": sigma_eq_error,
        "status": "PASS" if passed else "FAIL",
    }


def write_results(output_dir: Path, rows: list[dict], references: dict) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    csv_path = output_dir / "results.csv"
    fields = list(rows[0])
    with csv_path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)

    summary_path = output_dir / "verification_summary.md"
    lines = [
        "# Contact patch test verification summary",
        "",
        "| Formulation | Converged | Max epsilonEq | Max sigmaEq (Pa) | "
        "Average sigma_y error (%) | SigmaEq error (%) | Status |",
        "| --- | --- | ---: | ---: | ---: | ---: | --- |",
    ]
    for row in rows:
        lines.append(
            f"| `{row['variant']}` | {row['converged']} | "
            f"{row['max_epsilon_eq']:.8g} | {row['max_sigma_eq_pa']:.8g} | "
            f"{row['average_sigma_y_error_percent']:.6g} | "
            f"{row['sigma_eq_relative_error_percent']:.6g} | "
            f"{row['status']} |"
        )
    lines.extend(
        [
            "",
            "Acceptance limits: average sigma_y error <= "
            f"{references['max_average_sigma_y_error_percent']}%; max sigmaEq "
            "relative error <= "
            f"{references['max_sigma_eq_relative_error_percent']}%.",
            "",
        ]
    )
    summary_path.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    args = parse_args()
    verification_dir = args.verification_dir.resolve()
    parent = verification_dir.parent
    references = json.loads(
        (
            verification_dir
            / "reference"
            / "contact_patch_test_verification_references.json"
        ).read_text(encoding="utf-8")
    )
    variants = VARIANTS if args.variant == "all" else (args.variant,)

    require_commands(("blockMesh", "solids4Foam", "wmake"))
    rows = []
    for variant in variants:
        case_dir = verification_dir / "work" / variant
        log_path = case_dir / "log.solids4Foam"
        if args.reuse and log_path.is_file():
            print(f"Reusing {variant}")
        else:
            print(f"Running {variant}")
            prepare_case(parent, case_dir, variant)
            run_case(case_dir)
        row = read_result(case_dir, variant, references)
        rows.append(row)
        print(
            f"  sigma_y error={row['average_sigma_y_error_percent']:.6g}%, "
            f"sigmaEq error={row['sigma_eq_relative_error_percent']:.6g}%: "
            f"{row['status']}"
        )

    write_results(verification_dir / "postProcessing", rows, references)
    print(f"Wrote {verification_dir / 'postProcessing' / 'results.csv'}")
    return 0 if all(row["status"] == "PASS" for row in rows) else 1


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (OSError, RuntimeError, ValueError, KeyError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        sys.exit(2)
