#!/usr/bin/env python3
"""Run opt-in beamInCrossFlow verification studies in isolated case copies."""

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

SCRIPT = Path(__file__).resolve()
VERIFICATION = SCRIPT.parents[1]
TUTORIAL = VERIFICATION.parent
REFERENCE_FILE = VERIFICATION / "reference" / "beamInCrossFlow_verification_references.json"
WORK_ROOT = VERIFICATION / "work"
OUTPUT_ROOT = VERIFICATION / "postProcessing"


def fail(message: str) -> None:
    raise RuntimeError(message)


def command_exists(name: str) -> bool:
    return shutil.which(name) is not None


def replace_entry(path: Path, key: str, value: str) -> None:
    text = path.read_text()
    pattern = rf"^(\s*{re.escape(key)}\s+)[^;]+;"
    text, count = re.subn(pattern, rf"\g<1>{value};", text, flags=re.MULTILINE)
    if count != 1:
        fail(f"Expected one '{key}' entry in {path}, found {count}")
    path.write_text(text)


def set_case_form(case: Path, form: str, delta_t: float | None = None,
                  end_time: float | None = None) -> None:
    u_file = case / "0/fluid/U"
    mechanical = case / "constant/solid/mechanicalProperties"
    control = case / "system/controlDict.iqnils"
    if form == "original":
        replace_entry(u_file, "maxVelocity", "0.2")
        replace_entry(u_file, "timeVaryingEndTime", "4.0")
        replace_entry(u_file, "timeAtMaxVelocity", "4.0")
        replace_entry(mechanical, "E", "E [1 -1 -2 0 0 0 0] 1.4e6")
    elif form == "modified":
        replace_entry(u_file, "maxVelocity", "0.3")
        replace_entry(u_file, "timeVaryingEndTime", "1.0")
        replace_entry(u_file, "timeAtMaxVelocity", "1.0")
        replace_entry(mechanical, "E", "E [1 -1 -2 0 0 0 0] 1e4")
    else:
        fail(f"Unknown case form: {form}")
    if delta_t is not None:
        replace_entry(control, "deltaT", f"{delta_t:.8g}")
    if end_time is not None:
        replace_entry(control, "endTime", f"{end_time:.8g}")


def configure_output(case: Path, delta_t: float, end_time: float,
                     requested_interval: int | None) -> None:
    """Write fields and function-object data only when verification needs them."""
    interval = requested_interval or round(end_time / delta_t)
    if interval < 1:
        fail("Verification write interval must be at least one time step")
    replace_entry(case / "system/controlDict.iqnils", "writeInterval", str(interval))
    replace_entry(case / "system/functions", "writeInterval", str(interval))


def configure_time_scheme(case: Path, scheme: str) -> None:
    """Apply one scheme consistently to fluid and solid time derivatives."""
    fluid_schemes = case / "system/fluid/fvSchemes"
    fluid_text = fluid_schemes.read_text()
    fluid_text, count = re.subn(
        r"(ddtSchemes\s*\{\s*default\s+)[^;]+;",
        rf"\g<1>{scheme};",
        fluid_text,
        count=1,
        flags=re.DOTALL,
    )
    if count != 1:
        fail(f"Expected one ddtSchemes/default entry in {fluid_schemes}")
    fluid_schemes.write_text(fluid_text)

    # The solid has both first- and second-time-derivative operators. The
    # latter is used for the inertia term and must use the same setting.
    solid_schemes = case / "system/solid/fvSchemes"
    text = solid_schemes.read_text()
    for section in ("d2dt2Schemes", "ddtSchemes"):
        text, count = re.subn(
            rf"({section}\s*\{{\s*default\s+)[^;]+;",
            rf"\g<1>{scheme};",
            text,
            count=1,
            flags=re.DOTALL,
        )
        if count != 1:
            fail(f"Expected one {section}/default entry in {solid_schemes}")
    solid_schemes.write_text(text)

def refine_mesh(path: Path, factor: int) -> None:
    """Multiply every block cell-count tuple, retaining the tutorial topology."""
    text = path.read_text()

    def scale(match: re.Match[str]) -> str:
        counts = [int(value) * factor for value in match.group(1).split()]
        return "(" + " ".join(str(value) for value in counts) + ")"

    text, count = re.subn(r"hex\s+\([^)]*\)\s+\(([^()]+)\)\s+simpleGrading", lambda m: m.group(0).replace("(" + m.group(1) + ")", scale(m)), text)
    if count == 0:
        fail(f"No block cell counts found in {path}")
    path.write_text(text)


def copy_case(name: str) -> Path:
    destination = WORK_ROOT / name
    if destination.exists():
        shutil.rmtree(destination)
    ignore = shutil.ignore_patterns("verification", "regressionTests", "postProcessing", "processor*", "log.*")
    shutil.copytree(TUTORIAL, destination, symlinks=True, ignore=ignore)
    return destination


def run_case(case: Path, label: str, cores: int) -> None:
    allrun = case / "Allrun"
    if not allrun.is_file():
        fail(f"Missing {allrun}")
    log = case / "log.Allverify"
    command = [str(allrun), "iqnils"]
    if cores > 1:
        for dictionary in (
            case / "system/decomposeParDict",
            case / "system/fluid/decomposeParDict",
            case / "system/solid/decomposeParDict",
        ):
            replace_entry(dictionary, "numberOfSubdomains", str(cores))
        command.append("parallel")
    with log.open("w") as handle:
        result = subprocess.run(command, cwd=case, stdout=handle,
                                stderr=subprocess.STDOUT, text=True)
    if result.returncode:
        fail(f"{label} failed on {cores} core(s); see {log}")
    solver_log = case / "log.solids4Foam"
    if not solver_log.is_file():
        fail(f"{label} did not create {solver_log}")
    solver_text = solver_log.read_text(errors="replace")
    if re.search(r"FOAM FATAL|FOAM aborting|^ERROR$", solver_text, re.MULTILINE):
        fail(f"{label} failed inside Allrun; see {solver_log}")


def numeric_rows(path: Path):
    for line in path.read_text(errors="replace").splitlines():
        fields = line.split()
        if fields and re.fullmatch(r"[-+0-9.eE]+", fields[0]):
            yield line, fields


def find_displacement(case: Path) -> Path:
    candidates = list(case.glob("postProcessing/**/solidPointDisplacement*_displacement.dat"))
    if not candidates:
        fail(f"Point-A displacement data not found in {case}")
    return candidates[0]


def vector_at_or_before(path: Path, target: float | None) -> tuple[float, tuple[float, float, float]]:
    header = ""
    for line in path.read_text(errors="replace").splitlines():
        if line.startswith("#"):
            header = line
    names = header.replace("#", " ").split()
    indices = {name.lower(): i for i, name in enumerate(names)}
    rows = list(numeric_rows(path))
    if not rows:
        fail(f"No numeric displacement rows in {path}")
    selected = rows[-1]
    if target is not None:
        eligible = [row for row in rows if float(row[1][0]) <= target + 1e-10]
        if not eligible:
            fail(f"No displacement data at or before t={target} in {path}")
        selected = eligible[-1]
    line, fields = selected
    vectors = re.findall(r"\(([^()]*)\)", line)
    if vectors:
        try:
            vector = tuple(float(value) for value in vectors[0].split())
            if len(vector) == 3:
                return float(fields[0]), vector
        except ValueError:
            pass
    # OpenFOAM.com names columns Dx/Dy/Dz; other supported variants use 2..4.
    if all(key in indices for key in ("dx", "dy", "dz")):
        values = tuple(float(fields[indices[key]]) for key in ("dx", "dy", "dz"))
    elif len(fields) >= 4:
        values = tuple(float(fields[i]) for i in (1, 2, 3))
    else:
        fail(f"Could not parse displacement vector from {path}: {line}")
    return float(fields[0]), values


def find_force(case: Path) -> Path:
    candidates = list(case.glob("postProcessing/**/force*.dat")) + list(case.glob("postProcessing/**/forces.dat"))
    if not candidates:
        fail(f"Force data not found in {case}")
    return candidates[0]


def force_at_or_before(path: Path, target: float | None) -> tuple[float, tuple[float, float, float]]:
    rows = list(numeric_rows(path))
    if not rows:
        fail(f"No numeric force rows in {path}")
    if target is not None:
        rows = [row for row in rows if float(row[1][0]) <= target + 1e-10]
    if not rows:
        fail(f"No force data at or before t={target} in {path}")
    line, fields = rows[-1]
    vectors = re.findall(r"\(([^()]*)\)", line)
    if not vectors:
        # OpenFOAM.com writes total, pressure, and viscous vectors as nine
        # flat columns. The first vector is already the total force.
        if len(fields) >= 4:
            try:
                return float(fields[0]), tuple(float(fields[i]) for i in (1, 2, 3))
            except ValueError:
                pass
        fail(f"Could not parse force vector from {path}: {line}")
    parsed = [tuple(float(value) for value in vector.split()) for vector in vectors]
    if len(parsed[0]) != 3:
        fail(f"Invalid force vector in {path}: {line}")
    # The forces function object writes pressure and viscous contributions first.
    total = tuple(sum(vector[i] for vector in parsed[:2]) for i in range(3)) if len(parsed) >= 2 else parsed[0]
    return float(fields[0]), total


def cell_count(case: Path) -> int:
    total = 0
    for mesh in (case / "system/fluid/blockMeshDict", case / "system/solid/blockMeshDict"):
        for match in re.finditer(r"hex\s+\([^)]*\)\s+\(([^()]+)\)\s+simpleGrading", mesh.read_text()):
            values = [int(value) for value in match.group(1).split()]
            total += math.prod(values)
    return total


def relative_error(value: float, reference: float) -> float:
    return abs(value - reference) / abs(reference) if reference else abs(value - reference)


def extract(case: Path, evaluation_time: float | None) -> dict[str, float]:
    time, displacement = vector_at_or_before(find_displacement(case), evaluation_time)
    force_time, force = force_at_or_before(find_force(case), evaluation_time)
    return {"evaluation_time": time, "force_time": force_time, "ux": displacement[0],
            "uy": displacement[1], "uz": displacement[2], "fx": force[0],
            "fy": force[1], "fz": force[2], "cell_count": cell_count(case)}


def observed_order(coarse: float, medium: float, fine: float) -> float | None:
    numerator, denominator = coarse - medium, medium - fine
    if numerator == 0 or denominator == 0 or numerator * denominator <= 0:
        return None
    return math.log(abs(numerator / denominator), 2.0)


def reference_error_order(coarse_error: float, fine_error: float,
                          refinement_ratio: float) -> float | None:
    """Return the net order from reference-error reduction over a mesh sweep."""
    if refinement_ratio <= 1.0 or coarse_error < 0.0 or fine_error < 0.0:
        return None
    if fine_error == 0.0:
        return math.inf
    if coarse_error == 0.0:
        return None
    return math.log(coarse_error / fine_error) / math.log(refinement_ratio)


def study_cores(requested: str, study: str, refinement: int = 1) -> int:
    if requested != "auto":
        return int(requested)
    if study == "mesh":
        return {1: 1, 2: 4, 4: 8, 8: 64}[refinement]
    return 4


def create_mesh_plot(name: str) -> None:
    """Render the PNG convergence plot for a completed mesh sweep."""
    plot_script = VERIFICATION / "scripts" / "plotMeshConvergence.gnuplot"
    data = OUTPUT_ROOT / f"{name}.csv"
    output = OUTPUT_ROOT / f"{name}.png"
    result = subprocess.run(
        [
            "gnuplot",
            "-e",
            f"data='{data}'; output='{output}'",
            str(plot_script),
        ],
        cwd=VERIFICATION,
        text=True,
    )
    if result.returncode:
        fail(f"Could not create mesh plot {output}")
    print(f"Mesh plot: {output}")


def write_results(name: str, rows: list[dict], references: dict, order: float | None,
                  reference_orders: dict[str, float | None] | None = None,
                  reference_order_limit: float | None = None,
                  reference_row: int = -1,
                  published_references: dict | None = None) -> bool:
    OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)
    csv_path = OUTPUT_ROOT / f"{name}.csv"
    quantities = ("ux", "uy", "uz", "fx", "fy", "fz")
    columns = ["case", "study", "time_scheme", "mesh_level", "cell_count", "delta_t", "cores", "evaluation_time"]
    columns += list(quantities)
    columns += [item for quantity in quantities for item in
                (f"reference_{quantity}", f"{quantity}_relative_error", f"{quantity}_pass")]
    for source, source_references in (published_references or {}).items():
        for quantity in source_references:
            columns += [f"{source}_{quantity}", f"{source}_{quantity}_relative_error"]
    columns += ["observed_order", "pass"]
    columns += [f"{quantity}_reference_error_order" for quantity in (reference_orders or {})]
    passed = True
    for row_index, row in enumerate(rows):
        for quantity, spec in references.items():
            row[f"reference_{quantity}"] = spec["value"]
            row[f"{quantity}_relative_error"] = relative_error(row[quantity], spec["value"])
            row[f"{quantity}_pass"] = row[f"{quantity}_relative_error"] <= spec["tolerance"]
            if row_index == reference_row % len(rows) and spec.get("primary", True):
                passed = passed and row[f"{quantity}_pass"]
        for source, source_references in (published_references or {}).items():
            for quantity, value in source_references.items():
                row[f"{source}_{quantity}"] = value
                row[f"{source}_{quantity}_relative_error"] = relative_error(
                    row[quantity], value
                )
        row["observed_order"] = order if order is not None else ""
    if reference_order_limit is not None:
        passed = passed and all(
            reference_order is not None and reference_order > reference_order_limit
            for reference_order in (reference_orders or {}).values()
        )
    for row in rows:
        row["pass"] = passed
        for quantity, reference_order in (reference_orders or {}).items():
            row[f"{quantity}_reference_error_order"] = reference_order
    with csv_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    summary = OUTPUT_ROOT / "verification_summary.md"
    with summary.open("a") as handle:
        handle.write(f"## {name}\\n\\n")
        handle.write(f"- Result: {'PASS' if passed else 'FAIL'}\\n")
        if order is not None:
            handle.write(f"- Solution-change order for u_x(A): {order:.3f}\\n")
        for quantity, reference_order in (reference_orders or {}).items():
            displayed_order = (
                "infinite" if reference_order == math.inf
                else f"{reference_order:.3f}" if reference_order is not None
                else "undefined"
            )
            handle.write(f"- Reference-error order for {quantity}: {displayed_order}\\n")
        handle.write(f"- Data: `{csv_path.name}`\\n\\n")
    print(f"{name}: {'PASS' if passed else 'FAIL'}; results: {csv_path}")
    return passed


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", choices=("original",), required=True)
    parser.add_argument("--study", choices=("mesh",), required=True)
    parser.add_argument("--cores", default="auto", help="MPI ranks per case: positive integer or auto (default)")
    parser.add_argument("--write-interval", type=int, help="write every N time steps (default: final time only)")
    parser.add_argument(
        "--time-scheme",
        choices=("backward", "Euler"),
        default="backward",
        help="time scheme for the mesh diagnostic (default: backward)",
    )
    args = parser.parse_args()
    if args.cores != "auto" and (not args.cores.isdecimal() or int(args.cores) < 1):
        parser.error("--cores must be a positive integer or auto")
    if args.write_interval is not None and args.write_interval < 1:
        parser.error("--write-interval must be a positive integer")
    required_executables = ["python3", "blockMesh", "solids4Foam"]
    if args.study == "mesh":
        required_executables.append("gnuplot")
    for executable in required_executables:
        if not command_exists(executable):
            fail(f"Required executable '{executable}' is unavailable. Source OpenFOAM and build solids4foam first.")
    if not (TUTORIAL / "Allrun").is_file():
        fail(f"Tutorial not found at {TUTORIAL}")
    references = json.loads(REFERENCE_FILE.read_text())
    OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)
    (OUTPUT_ROOT / "verification_summary.md").write_text("# beamInCrossFlow verification summary\\n\\n")
    rows: list[dict] = []
    if args.study == "mesh":
        mesh_end_time = references["mesh"]["endTime"]
        mesh_factors = references["mesh"]["refinementFactors"]
        mesh_delta_ts = references["mesh"]["deltaTs"]
        if len(mesh_factors) != len(mesh_delta_ts):
            fail("mesh refinementFactors and deltaTs must have the same length")
        scheme_suffix = "" if args.time_scheme == "backward" else f"_{args.time_scheme}"
        for level, (factor, mesh_delta_t) in enumerate(zip(mesh_factors, mesh_delta_ts)):
            case = copy_case(f"original_mesh{scheme_suffix}_{factor}x")
            set_case_form(case, "original", mesh_delta_t, mesh_end_time)
            configure_time_scheme(case, args.time_scheme)
            configure_output(case, mesh_delta_t, mesh_end_time, args.write_interval)
            refine_mesh(case / "system/fluid/blockMeshDict", factor)
            refine_mesh(case / "system/solid/blockMeshDict", factor)
            cores = study_cores(args.cores, "mesh", factor)
            run_case(case, f"mesh level {factor}x", cores)
            row = extract(case, None)
            row.update({"case": "original", "study": "mesh", "time_scheme": args.time_scheme, "mesh_level": level, "delta_t": mesh_delta_t, "cores": cores})
            rows.append(row)
        name = f"original_mesh_sweep{scheme_suffix}"
        order = observed_order(rows[-3]["ux"], rows[-2]["ux"], rows[-1]["ux"])
        refinement_ratio = mesh_factors[-1] / mesh_factors[0]
        reference_orders = {
            quantity: reference_error_order(
                relative_error(rows[0][quantity], spec["value"]),
                relative_error(rows[-1][quantity], spec["value"]),
                refinement_ratio,
            )
            for quantity, spec in references["original"]["references"].items()
            if spec.get("primary", True)
        }
        passed = write_results(
            name,
            rows,
            references["original"]["references"],
            order,
            reference_orders,
            references["mesh"]["minimumReferenceErrorOrder"],
            -1,
            references["original"]["publishedReferences"],
        )
        create_mesh_plot(name)
        return 0 if passed else 1
    fail(f"Unsupported study: {args.study}")


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (RuntimeError, subprocess.SubprocessError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        sys.exit(2)
