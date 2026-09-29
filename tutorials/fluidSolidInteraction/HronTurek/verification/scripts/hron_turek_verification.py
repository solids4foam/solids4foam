#!/usr/bin/env python3
"""Run opt-in HronTurek verification studies in isolated case copies.

The studies compare the periodic response of the Turek-Hron FSI3 benchmark
(mean, amplitude and frequency of the point-A displacement and of the drag and
lift on the cylinder and plate) with the published reference values.
"""

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
REFERENCE_FILE = VERIFICATION / "reference" / "HronTurek_verification_references.json"
WORK_ROOT = VERIFICATION / "work"
OUTPUT_ROOT = VERIFICATION / "postProcessing"

QUANTITIES = ("ux", "uy", "drag", "lift")
STATISTICS = ("mean", "amplitude", "frequency")


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


def dictionary_scalar(path: Path, key: str) -> float:
    match = re.search(
        rf"^\s*{re.escape(key)}\s+([-+0-9.eE]+)\s*;",
        path.read_text(),
        flags=re.MULTILINE,
    )
    if not match:
        fail(f"Could not read '{key}' from {path}")
    return float(match.group(1))


# ---------------------------------------------------------------------------
# Case preparation
# ---------------------------------------------------------------------------

def copy_case(name: str) -> Path:
    destination = WORK_ROOT / name
    if destination.exists():
        shutil.rmtree(destination)
    ignore = shutil.ignore_patterns(
        "verification", "regressionTests", "postProcessing", "processor*",
        "log.*", "*.pdf", "case.foam",
    )
    shutil.copytree(TUTORIAL, destination, symlinks=True, ignore=ignore)
    return destination


def refine_mesh(path: Path, factor: int) -> None:
    """Multiply the in-plane block divisions, keeping the single spanwise cell."""
    text = path.read_text()

    def scale(match: re.Match[str]) -> str:
        counts = [int(value) for value in match.group(2).split()]
        if len(counts) != 3:
            fail(f"Unexpected block division tuple in {path}: {match.group(2)}")
        counts = [counts[0] * factor, counts[1] * factor, counts[2]]
        return match.group(1) + "(" + " ".join(str(value) for value in counts) + ")"

    text, count = re.subn(
        r"(hex\s+\([^)]*\)\s+)\(([^()]+)\)", scale, text
    )
    if count == 0:
        fail(f"No block divisions found in {path}")
    path.write_text(text)


def cell_count(case: Path) -> int:
    total = 0
    for mesh in (case / "system/fluid/blockMeshDict", case / "system/solid/blockMeshDict"):
        for match in re.finditer(r"hex\s+\([^)]*\)\s+\(([^()]+)\)", mesh.read_text()):
            total += math.prod(int(value) for value in match.group(1).split())
    return total


def configure_case(case: Path, spec: dict, factor: int, delta_t: float,
                   end_time: float, write_interval: int | None,
                   cores: int) -> None:
    refine_mesh(case / "system/fluid/blockMeshDict", factor)
    refine_mesh(case / "system/solid/blockMeshDict", factor)

    control = case / "system/controlDict"
    replace_entry(control, "deltaT", f"{delta_t:.8g}")
    replace_entry(control, "endTime", f"{end_time:.8g}")
    # Write fields only when needed; the point displacement and force
    # histories are written every time step regardless.
    interval = write_interval or round(end_time / delta_t)
    if interval < 1:
        fail("Verification write interval must be at least one time step")
    replace_entry(control, "writeInterval", str(interval))

    # The published drag and lift act on the cylinder and the plate together
    # and are given per unit depth for a fluid of density 1000 kg/m^3.
    functions = case / "system/functions"
    replace_entry(functions, "patches", "(plate cylinder)")
    replace_entry(functions, "rhoInf", f"{spec['fluid_density']:.8g}")

    # The benchmark specifies a St. Venant-Kirchhoff material.
    replace_entry(
        case / "constant/solid/mechanicalProperties",
        "type",
        spec["mechanical_law"],
    )

    # The tutorial's interface tolerance of 1e-6 sits at the floor below
    # which the IQN-ILS secant stalls on this case (about one step in five
    # thousand), so a long run would abort. Both couplings use the same,
    # still stringent, tolerance so that the coupling study compares like with
    # like.
    for coupling in ("iqnils", "robin"):
        replace_entry(
            case / f"constant/fsiProperties.{coupling}",
            "outerCorrTolerance",
            f"{spec['outer_corr_tolerance']:.8g}",
        )

    if cores > 1:
        # The tutorial decomposes with the simple method on a fixed 2x2
        # arrangement; scotch accepts any rank count.
        replace_entry(case / "system/decomposeParDict", "numberOfSubdomains", str(cores))
        for dictionary in (
            case / "system/fluid/decomposeParDict",
            case / "system/solid/decomposeParDict",
        ):
            replace_entry(dictionary, "numberOfSubdomains", str(cores))
            replace_entry(dictionary, "method", "scotch")


def run_case(case: Path, label: str, cores: int, coupling: str,
             reuse: bool) -> None:
    solver_log = case / "log.solids4Foam"
    if reuse and solver_log.is_file() and re.search(
        r"^End\s*$", solver_log.read_text(errors="replace"), re.MULTILINE
    ):
        print(f"{label}: reusing the completed run in {case}")
        return
    allrun = case / "Allrun"
    if not allrun.is_file():
        fail(f"Missing {allrun}")
    command = [str(allrun), coupling]
    if cores > 1:
        command.append("parallel")
    log = case / "log.Allverify"
    print(f"{label}: running {' '.join(command[1:])} in {case}")
    with log.open("w") as handle:
        result = subprocess.run(command, cwd=case, stdout=handle,
                                stderr=subprocess.STDOUT, text=True)
    if result.returncode:
        fail(f"{label} failed on {cores} core(s); see {log}")
    if not solver_log.is_file():
        fail(f"{label} did not create {solver_log}")
    solver_text = solver_log.read_text(errors="replace")
    if re.search(
        r"FOAM FATAL|FOAM aborting|PRTE ERROR|No sockets were able|^ERROR$",
        solver_text,
        re.MULTILINE,
    ):
        fail(f"{label} failed inside Allrun; see {solver_log}")
    if not re.search(r"^End\s*$", solver_text, re.MULTILINE):
        fail(f"{label} did not run to completion; see {solver_log}")


# ---------------------------------------------------------------------------
# History extraction
# ---------------------------------------------------------------------------

def numeric_rows(path: Path):
    for line in path.read_text(errors="replace").splitlines():
        fields = line.replace("(", " ").replace(")", " ").split()
        if fields and re.fullmatch(r"[-+0-9.eE]+", fields[0]):
            yield fields


def find_file(case: Path, patterns: tuple[str, ...], description: str) -> Path:
    for pattern in patterns:
        candidates = sorted(case.glob(pattern))
        if candidates:
            return candidates[0]
    fail(f"{description} not found in {case}")


def displacement_history(case: Path) -> tuple[list[float], list[float], list[float]]:
    path = find_file(
        case,
        ("postProcessing/**/solidPointDisplacement*_displacement.dat",
         "postProcessing/**/solidPointDisplacement*.dat"),
        "Point-A displacement history",
    )
    time, ux, uy = [], [], []
    for fields in numeric_rows(path):
        if len(fields) < 3:
            continue
        time.append(float(fields[0]))
        ux.append(float(fields[1]))
        uy.append(float(fields[2]))
    if not time:
        fail(f"No displacement rows in {path}")
    return time, ux, uy


def force_history(case: Path, thickness: float) -> tuple[list[float], list[float], list[float]]:
    path = find_file(
        case,
        ("postProcessing/**/force.dat", "postProcessing/**/forces.dat",
         "forces/**/forces.dat"),
        "Force history",
    )
    time, drag, lift = [], [], []
    for fields in numeric_rows(path):
        values = [float(value) for value in fields]
        if len(values) >= 13:
            # OpenFOAM.org and foam-extend: pressure, viscous, then moments
            fx, fy = values[1] + values[4], values[2] + values[5]
        elif len(values) >= 4:
            # OpenFOAM.com: total, pressure, viscous
            fx, fy = values[1], values[2]
        else:
            continue
        time.append(values[0])
        drag.append(fx / thickness)
        lift.append(fy / thickness)
    if not time:
        fail(f"No force rows in {path}")
    return time, drag, lift


def reference_history(path: Path) -> dict[str, list[float]]:
    columns: dict[str, list[float]] = {name: [] for name in ("time", *QUANTITIES)}
    with path.open() as handle:
        reader = csv.DictReader(
            line for line in handle if not line.startswith("#")
        )
        for row in reader:
            for name in columns:
                columns[name].append(float(row[name]))
    return columns


# ---------------------------------------------------------------------------
# Periodic statistics
# ---------------------------------------------------------------------------

def periodic_statistics(time: list[float], values: list[float],
                        window: float,
                        fundamental: list[float] | None = None) -> dict[str, float]:
    """Mean, amplitude and frequency of a periodic signal.

    The benchmark reports the mean as (max + min)/2 and the amplitude as
    (max - min)/2 over the last full period. Periods are delimited by upward
    crossings of the window mean with a hysteresis band, so that the harmonics
    in the force signals do not register as additional periods. The frequency
    is the average over all full periods in the window.
    """
    end = time[-1]
    start_index = next((i for i, t in enumerate(time) if t >= end - window), None)
    if start_index is None or len(time) - start_index < 8:
        fail(f"The analysis window of {window} s contains too few samples")
    t = time[start_index:]
    y = values[start_index:]
    mean = sum(y) / len(y)
    band = 0.25 * (max(y) - min(y))
    if band <= 0.0:
        fail("The signal is constant over the analysis window")

    crossings: list[float] = []
    armed = y[0] < mean - band
    for k in range(1, len(y)):
        if y[k] < mean - band:
            armed = True
        elif armed and y[k] >= mean:
            a, b = y[k - 1] - mean, y[k] - mean
            crossings.append(t[k - 1] + a / (a - b) * (t[k] - t[k - 1]) if a != b else t[k])
            armed = False
    if len(crossings) < 3:
        fail(f"Fewer than two full periods found in the analysis window of {window} s")

    def extrema(t_start: float, t_end: float) -> tuple[float, float]:
        segment = [v for tv, v in zip(t, y) if t_start <= tv <= t_end]
        return max(segment), min(segment)

    # The extrema are taken over the fundamental period (that of u_y) when it
    # is given: u_x and the drag oscillate at twice that frequency with
    # alternating troughs, and the benchmark reports their mean and amplitude
    # over the full period of the plate motion
    bounds = fundamental if fundamental is not None else crossings
    last_max, last_min = extrema(bounds[-2], bounds[-1])
    previous_max, previous_min = extrema(bounds[-3], bounds[-2])
    # Mean amplitude over the last two and the two preceding periods, to
    # detect a trend without mistaking cycle-to-cycle scatter for one
    amplitudes = [0.5 * (high - low) for high, low in
                  (extrema(a, b) for a, b in zip(bounds[:-1], bounds[1:]))]
    if len(amplitudes) >= 4:
        late_amplitude = 0.5 * (amplitudes[-1] + amplitudes[-2])
        early_amplitude = 0.5 * (amplitudes[-3] + amplitudes[-4])
    else:
        late_amplitude = amplitudes[-1]
        early_amplitude = amplitudes[-2]
    return {
        "crossings": crossings,
        "mean": 0.5 * (last_max + last_min),
        "amplitude": 0.5 * (last_max - last_min),
        "frequency": (len(crossings) - 1) / (crossings[-1] - crossings[0]),
        "previous_amplitude": 0.5 * (previous_max - previous_min),
        "late_amplitude": late_amplitude,
        "early_amplitude": early_amplitude,
        "periods": len(crossings) - 1,
        "last_maximum_time": max(
            ((v, tv) for tv, v in zip(t, y) if crossings[-2] <= tv <= crossings[-1])
        )[1],
    }


def extract(case: Path, spec: dict, window: float,
            end_time: float) -> dict[str, float]:
    thickness = spec["thickness_m"]
    time_d, ux, uy = displacement_history(case)
    time_f, drag, lift = force_history(case, thickness)
    for path_time, name in ((time_d[-1], "displacement"), (time_f[-1], "force")):
        if not math.isclose(path_time, end_time, rel_tol=1e-6, abs_tol=1e-9):
            fail(f"The {name} history ends at t={path_time:g}, not at the end time t={end_time:g}")
    row: dict[str, float] = {"cell_count": cell_count(case)}
    histories = {"ux": (time_d, ux), "uy": (time_d, uy),
                 "drag": (time_f, drag), "lift": (time_f, lift)}
    fundamental = periodic_statistics(time_d, uy, window)["crossings"]
    for quantity, (time, values) in histories.items():
        statistics = periodic_statistics(time, values, window, fundamental)
        for name in STATISTICS:
            row[f"{quantity}_{name}"] = statistics[name]
        row[f"{quantity}_previous_amplitude"] = statistics["previous_amplitude"]
        row[f"{quantity}_late_amplitude"] = statistics["late_amplitude"]
        row[f"{quantity}_early_amplitude"] = statistics["early_amplitude"]
        row[f"{quantity}_periods"] = statistics["periods"]
        row[f"{quantity}_last_maximum_time"] = statistics["last_maximum_time"]
    return row


def write_aligned_history(case: Path, name: str, spec: dict, window: float,
                          row: dict[str, float]) -> None:
    """Write the last window of the run, with time relative to the last uy peak."""
    time_d, ux, uy = displacement_history(case)
    time_f, drag, lift = force_history(case, spec["thickness_m"])
    shift = row["uy_last_maximum_time"]
    path = OUTPUT_ROOT / f"{name}_history.csv"
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["time", "drag", "lift", "ux", "uy"])
        forces = {round(t, 9): (d, l) for t, d, l in zip(time_f, drag, lift)}
        for t, x, y in zip(time_d, ux, uy):
            if t < time_d[-1] - window:
                continue
            force = forces.get(round(t, 9))
            if force is None:
                continue
            writer.writerow([f"{t - shift:.6f}", f"{force[0]:.6e}",
                             f"{force[1]:.6e}", f"{x:.6e}", f"{y:.6e}"])


def write_reference_history(spec: dict, window: float) -> None:
    reference = reference_history(VERIFICATION / "reference" / spec["history_reference"])
    statistics = periodic_statistics(reference["time"], reference["uy"], window)
    shift = statistics["last_maximum_time"]
    path = OUTPUT_ROOT / "reference_history.csv"
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["time", "drag", "lift", "ux", "uy"])
        for i, t in enumerate(reference["time"]):
            if t < reference["time"][-1] - window:
                continue
            writer.writerow([f"{t - shift:.6f}"] + [
                f"{reference[name][i]:.6e}" for name in ("drag", "lift", "ux", "uy")
            ])


def reference_history_statistics(spec: dict, window: float) -> dict[str, float]:
    reference = reference_history(VERIFICATION / "reference" / spec["history_reference"])
    values: dict[str, float] = {}
    fundamental = periodic_statistics(reference["time"], reference["uy"], window)["crossings"]
    for quantity in QUANTITIES:
        statistics = periodic_statistics(reference["time"], reference[quantity], window, fundamental)
        for name in STATISTICS:
            values[f"{quantity}_{name}"] = statistics[name]
    return values


# ---------------------------------------------------------------------------
# Checks and reporting
# ---------------------------------------------------------------------------

def relative_error(value: float, reference: float) -> float:
    return abs(value - reference) / abs(reference) if reference else abs(value - reference)


def periodicity_failures(row: dict[str, float], tolerance: float,
                         report_only: tuple[str, ...] = ("drag",)) -> list[str]:
    failures = []
    for quantity in QUANTITIES:
        change = relative_error(row[f"{quantity}_late_amplitude"],
                                row[f"{quantity}_early_amplitude"])
        row[f"{quantity}_amplitude_change"] = change
        # The FSI3 drag amplitude scatters by about 2% from cycle to cycle with
        # no trend, so its change is reported but not tested
        if quantity not in report_only and change > tolerance:
            failures.append(
                f"{quantity} mean amplitude over the last two periods differs by "
                f"{100 * change:.2f}% from the two before (tolerance {100 * tolerance:g}%)"
            )
    return failures


def create_plot(script_name: str, run_name: str, options: str = "") -> None:
    plot_script = VERIFICATION / "scripts" / script_name
    output = OUTPUT_ROOT / f"{run_name}_history.png"
    result = subprocess.run(
        [
            "gnuplot",
            "-e",
            f"run='{OUTPUT_ROOT / (run_name + '_history.csv')}'; "
            f"reference='{OUTPUT_ROOT / 'reference_history.csv'}'; "
            f"output='{output}'; label='{run_name}'" + options,
            str(plot_script),
        ],
        cwd=VERIFICATION,
        text=True,
    )
    if result.returncode:
        fail(f"Could not create the history plot {output}")
    print(f"History plot: {output}")


def study_cores(requested: str, factor: int, spec: dict) -> int:
    if requested != "auto":
        return int(requested)
    return int(spec["cores"].get(str(factor), spec["cores"]["default"]))


def format_value(quantity: str, name: str, value: float) -> str:
    if name == "frequency":
        return f"{value:.3f}"
    if quantity in ("ux", "uy"):
        return f"{1000 * value:.3f}"
    return f"{value:.2f}"


def summary_table(rows: list[dict], references: dict, published: dict,
                  key: str) -> str:
    lines = [
        "| Quantity | " + " | ".join(str(row[key]) for row in rows)
        + " | Featflow L4 | " + " | ".join(published) + " |",
        "|---|" + "---:|" * (len(rows) + 1 + len(published)),
    ]
    units = {"ux": "mm", "uy": "mm", "drag": "N/m", "lift": "N/m"}
    for quantity in QUANTITIES:
        for name in STATISTICS:
            column = f"{quantity}_{name}"
            unit = "Hz" if name == "frequency" else units[quantity]
            cells = [format_value(quantity, name, row[column]) for row in rows]
            cells.append(format_value(quantity, name, references[column]["value"]))
            cells += [format_value(quantity, name, source[column]) for source in published.values()]
            lines.append(f"| {quantity} {name} ({unit}) | " + " | ".join(cells) + " |")
    return "\n".join(lines)


def write_summary(title: str, passed: bool, notes: list[str], table: str,
                  csv_name: str) -> None:
    summary = OUTPUT_ROOT / "verification_summary.md"
    with summary.open("a") as handle:
        handle.write(f"## {title}\n\n")
        handle.write(f"- Result: {'PASS' if passed else 'FAIL'}\n")
        for note in notes:
            handle.write(f"- {note}\n")
        handle.write(f"- Data: `{csv_name}`\n\n")
        handle.write(table + "\n\n")


def check_references(row: dict, references: dict) -> list[str]:
    failures = []
    for column, spec in references.items():
        error = relative_error(row[column], spec["value"])
        row[f"{column}_reference"] = spec["value"]
        row[f"{column}_error"] = error
        row[f"{column}_pass"] = error <= spec["tolerance"]
        if spec.get("primary", True) and error > spec["tolerance"]:
            failures.append(
                f"{column} = {row[column]:.6g} differs from the reference "
                f"{spec['value']:.6g} by {100 * error:.2f}% (tolerance {100 * spec['tolerance']:g}%)"
            )
    return failures


def write_csv(path: Path, rows: list[dict]) -> None:
    columns: list[str] = []
    for row in rows:
        for key in row:
            if key not in columns:
                columns.append(key)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)


# ---------------------------------------------------------------------------
# Studies
# ---------------------------------------------------------------------------

def run_level(args: argparse.Namespace, spec: dict, factor: int, coupling: str,
              end_time: float, window: float, name: str) -> dict:
    index = spec["mesh"]["refinementFactors"].index(factor)
    delta_t = spec["mesh"]["deltaTs"][index]
    cores = study_cores(args.cores, factor, spec)
    case = WORK_ROOT / name
    if not (args.reuse and case.is_dir()):
        case = copy_case(name)
        configure_case(case, spec, factor, delta_t, end_time,
                       args.write_interval, cores)
        if spec.get("benchmark") == "fsi2":
            configure_fsi2(case, spec["fsi2"])
    run_case(case, name, cores, coupling, args.reuse)
    row: dict = {"case": name, "coupling": coupling, "mesh_level": index + 1,
                 "refinement": factor, "delta_t": delta_t, "cores": cores,
                 "end_time": end_time, "window": window}
    if spec.get("benchmark") == "fsi2":
        # A reused run keeps the rank count it was run with
        row["cores"] = len(list(case.glob("processor*"))) or 1
        row["execution_time"] = execution_time(case)
        row.update(iqnils_iteration_summary(case))
        if "coupled_time_steps" in row:
            row["mean_outer_iterations"] = (row["total_outer_iterations"]
                                            / row["coupled_time_steps"])
    if args.quick:
        return row
    row.update(extract(case, spec, window, end_time))
    write_aligned_history(case, name, spec, window, row)
    return row


def mesh_study(args: argparse.Namespace, spec: dict) -> bool:
    end_time = args.end_time or spec["mesh"]["endTime"]
    window = args.window or spec["mesh"]["analysisWindow"]
    factors = [int(level) for level in args.levels.split(",")] if args.levels else spec["mesh"]["defaultLevels"]
    for factor in factors:
        if factor not in spec["mesh"]["refinementFactors"]:
            fail(f"Unsupported refinement factor {factor}; choose from {spec['mesh']['refinementFactors']}")
    if args.quick:
        end_time = spec["quick"]["endTime"]
    # FSI2 runs, results and plots carry a prefix so that they do not
    # overwrite the FSI3 ones
    prefix = "fsi2_" if spec.get("benchmark") == "fsi2" else ""
    rows = []
    notes = []
    failures = []
    for factor in factors:
        name = f"{prefix}{args.coupling}_mesh_{factor}x"
        row = run_level(args, spec, factor, args.coupling, end_time, window, name)
        rows.append(row)
        if args.quick:
            continue
        failures += [f"{name}: {message}" for message in
                     periodicity_failures(row, spec["periodicity"]["amplitudeTolerance"],
                                          tuple(spec["periodicity"].get("reportOnly", ("drag",))))]
        level_failures = check_references(row, spec["references"])
        if factor == factors[-1]:
            failures += [f"{name}: {message}" for message in level_failures]
    csv_name = f"{prefix}{args.coupling}_mesh_sweep.csv"
    OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)
    if args.quick:
        write_csv(OUTPUT_ROOT / csv_name, rows)
        print(f"Quick run completed for levels {factors}; results: {OUTPUT_ROOT / csv_name}")
        return True

    # Convergence: a primary reference error may not grow from the coarsest
    # to the finest level by more than the allowed margin.
    if len(rows) > 1:
        growth_limit = spec["mesh"]["maximumErrorGrowth"]
        for column, reference in spec["references"].items():
            if not reference.get("primary", True):
                continue
            growth = rows[-1][f"{column}_error"] - rows[0][f"{column}_error"]
            # The FSI2 lift amplitude does not converge monotonically, so its
            # error growth is reported but not tested
            if not reference.get("errorGrowthChecked", True):
                rows[-1][f"{column}_error_growth"] = growth
                continue
            rows[-1][f"{column}_error_growth"] = growth
            if growth > growth_limit:
                failures.append(
                    f"{column} reference error grew from {100 * rows[0][f'{column}_error']:.2f}% "
                    f"to {100 * rows[-1][f'{column}_error']:.2f}% between the coarsest and finest levels"
                )
    passed = not failures
    for row in rows:
        row["pass"] = passed
    write_csv(OUTPUT_ROOT / csv_name, rows)

    published = dict(spec.get("publishedReferences", {}))
    published["Turek-Hron history"] = reference_history_statistics(spec, window)
    table = summary_table(rows, spec["references"], published, "refinement")
    notes += failures
    for row in rows:
        notes.append(
            f"{row['case']}: {row['cell_count']} cells, dt = {row['delta_t']:g} s, "
            f"{row['uy_periods']} uy periods in the analysis window"
        )
        if "mean_outer_iterations" in row:
            notes[-1] += (
                f", {row['cores']} core(s), {row['execution_time']:.0f} s, mean "
                f"{row['mean_outer_iterations']:.2f} (maximum "
                f"{row['maximum_outer_iterations']:.0f}) FSI iterations per step"
            )
    title = f"{prefix.upper().replace('_', ' ')}{args.coupling} mesh study"
    write_summary(title, passed, notes, table, csv_name)
    write_reference_history(spec, window)
    if command_exists("gnuplot"):
        # FSI2 periods are about 0.52 s long; two are shown
        options = "; benchmark='FSI2'; xmin=-1.1" if prefix else ""
        create_plot("plotPeriodicHistory.gnuplot", rows[-1]["case"], options)
    print(f"{title}: {'PASS' if passed else 'FAIL'}; results: {OUTPUT_ROOT / csv_name}")
    for message in failures:
        print(f"  - {message}")
    return passed


def robin_residual_summary(case: Path, end_time: float) -> dict[str, float]:
    path = case / "postProcessing/fsiResiduals.dat"
    if not path.is_file():
        fail(f"Robin residual data not found in {case}")
    header = path.read_text(errors="replace").splitlines()[0].split()
    rows = [fields for fields in numeric_rows(path)]
    if not rows or any(len(fields) < 5 for fields in rows):
        fail(f"Robin residual columns are missing from {path}")
    # Columns: time, iteration, displacement, pressure and leakage-flux
    # residuals and (when written) the Robin convergence state (0 not
    # converged, 1 converged, 2 stalled within the stall tolerance).
    state_index = header.index("robinConvergenceState") if "robinConvergenceState" in header else None
    final_by_time: dict[float, list[float]] = {}
    states: dict[float, float] = {}
    for fields in rows:
        values = [float(value) for value in fields[:5]]
        final_by_time[values[0]] = values
        if state_index is not None and len(fields) > state_index:
            states[values[0]] = float(fields[state_index])
    if not math.isclose(max(final_by_time), end_time, rel_tol=1e-6, abs_tol=1e-9):
        fail(f"{path} ends at t={max(final_by_time):g}, not at the end time t={end_time:g}")
    properties = case / "constant/fsiProperties.robin"
    tolerances = (
        dictionary_scalar(properties, "outerCorrTolerance"),
        dictionary_scalar(properties, "robinPressureTolerance"),
        dictionary_scalar(properties, "robinFluxTolerance"),
    )
    coupling_start = dictionary_scalar(properties, "couplingStartTime")

    def converged(values: list[float]) -> bool:
        if values[0] in states:
            return states[values[0]] > 0
        return all(values[2 + i] <= tolerances[i] for i in range(3))

    coupled = [values for time, values in final_by_time.items() if time > coupling_start]
    unconverged = [values for values in coupled if not converged(values)]
    if unconverged:
        fail(
            f"{len(unconverged)} Robin time step(s) reached nOuterCorr without "
            "satisfying all convergence criteria"
        )
    return {
        "coupled_time_steps": len(coupled),
        "total_outer_iterations": sum(values[1] for values in coupled),
        "maximum_outer_iterations": max(values[1] for values in coupled),
        "maximum_final_displacement_residual": max(values[2] for values in coupled),
        "maximum_final_pressure_residual": max(values[3] for values in coupled),
        "maximum_final_flux_residual": max(values[4] for values in coupled),
    }


def iqnils_iteration_summary(case: Path) -> dict[str, float]:
    path = case / "postProcessing/fsiResiduals.dat"
    if not path.is_file():
        return {}
    final_by_time: dict[float, float] = {}
    for fields in numeric_rows(path):
        final_by_time[float(fields[0])] = float(fields[1])
    coupling_start = dictionary_scalar(case / "constant/fsiProperties.iqnils", "couplingStartTime")
    coupled = [n for time, n in final_by_time.items() if time > coupling_start]
    if not coupled:
        fail(f"No coupled time steps recorded in {path}")
    return {
        "coupled_time_steps": len(coupled),
        "total_outer_iterations": sum(coupled),
        "maximum_outer_iterations": max(coupled),
    }


def coupled_histories(case: Path, spec: dict,
                      coupling_start: float) -> dict[str, tuple[list[float], list[float]]]:
    """Histories of the four quantities after the coupling start."""
    time_d, ux, uy = displacement_history(case)
    time_f, drag, lift = force_history(case, spec["thickness_m"])
    histories = {}
    for quantity, (time, values) in (("ux", (time_d, ux)), ("uy", (time_d, uy)),
                                     ("drag", (time_f, drag)), ("lift", (time_f, lift))):
        selected = [(t, v) for t, v in zip(time, values) if t > coupling_start]
        histories[quantity] = ([t for t, _ in selected], [v for _, v in selected])
    return histories


def enable_restart(case: Path) -> None:
    """Ask the solid model to write, or read, the full restart state."""
    path = case / "constant/solid/solidProperties"
    text, count = re.subn(r"(Coeffs\"?\s*\{)", r"\g<1>\n    restart yes;", path.read_text(), count=1)
    if count != 1:
        fail(f"Could not find the solid model coefficients in {path}")
    path.write_text(text)


def coupling_level(args: argparse.Namespace, spec: dict, factor: int,
                   end_time: float) -> dict:
    """Run both couplings on one mesh level from a common start state."""
    index = spec["mesh"]["refinementFactors"].index(factor)
    delta_t = spec["mesh"]["deltaTs"][index]
    # The restarts below are run in serial so that the common start state
    # does not have to be decomposed
    cores = 1
    coupling_start = dictionary_scalar(
        TUTORIAL / "constant/fsiProperties.iqnils", "couplingStartTime"
    )

    # The fluid boundary conditions of the two variants differ even before
    # the coupling starts (elasticWallPressure is not a zero-gradient wall
    # while the plate is at rest), so separate runs from t = 0 would enter
    # the coupling from different flow states. Both variants therefore
    # restart from the state of one uncoupled Dirichlet-Neumann run.
    base_name = f"precoupling_{factor}x"
    base = WORK_ROOT / base_name
    if not (args.reuse and base.is_dir()):
        base = copy_case(base_name)
        configure_case(base, spec, factor, delta_t, coupling_start, None, cores)
        replace_entry(base / "system/controlDict", "writePrecision", "12")
        enable_restart(base)
    run_case(base, base_name, cores, "iqnils", args.reuse)
    start_directory = base / f"{coupling_start:g}"
    if not start_directory.is_dir():
        fail(f"The uncoupled run did not write the start state {start_directory}")

    rows = []
    for coupling in ("iqnils", "robin"):
        name = f"{coupling}_coupling_{factor}x"
        case = WORK_ROOT / name
        if not (args.reuse and case.is_dir()):
            case = copy_case(name)
            configure_case(case, spec, factor, delta_t, end_time,
                           args.write_interval, cores)
            enable_restart(case)
            shutil.copytree(start_directory, case / start_directory.name)
            if coupling == "robin":
                for field, condition in (("p", "elasticWallPressure"),
                                         ("U", "elasticWallVelocity")):
                    path = case / start_directory.name / "fluid" / field
                    text, count = re.subn(
                        r"(\bplate\s*\{\s*type\s+)\w+;",
                        rf"\g<1>{condition};",
                        path.read_text(),
                    )
                    if count != 1:
                        fail(f"Could not set the Robin condition on plate in {path}")
                    path.write_text(text)
        run_case(case, name, cores, coupling, args.reuse)
        row: dict = {"case": name, "coupling": coupling, "refinement": factor,
                     "delta_t": delta_t, "cores": cores, "end_time": end_time,
                     "cell_count": cell_count(case)}
        if coupling == "robin":
            row.update(robin_residual_summary(case, end_time))
        else:
            row.update(iqnils_iteration_summary(case))
        rows.append(row)
    histories = {
        row["coupling"]: coupled_histories(WORK_ROOT / row["case"], spec, coupling_start)
        for row in rows
    }
    iqnils, robin = rows
    for quantity in QUANTITIES:
        time_i, values_i = histories["iqnils"][quantity]
        time_r, values_r = histories["robin"][quantity]
        if len(time_i) != len(time_r) or any(
            not math.isclose(a, b, rel_tol=1e-9, abs_tol=1e-12) for a, b in zip(time_i, time_r)
        ):
            fail(f"The {quantity} histories of the two couplings are sampled at different times")
        scale = max(abs(v) for v in values_i)
        robin[f"{quantity}_vs_iqnils"] = max(
            abs(a - b) for a, b in zip(values_i, values_r)
        ) / scale
        robin[f"{quantity}_scale"] = scale

    OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)
    history_name = f"coupling_comparison_{factor}x_history.csv"
    with (OUTPUT_ROOT / history_name).open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["time"] + [f"{c}_{q}" for c in ("iqnils", "robin") for q in QUANTITIES])
        time = histories["iqnils"]["ux"][0]
        for i, t in enumerate(time):
            writer.writerow([f"{t:.6f}"] + [
                f"{histories[c][q][1][i]:.6e}" for c in ("iqnils", "robin") for q in QUANTITIES
            ])
    if command_exists("gnuplot"):
        output = OUTPUT_ROOT / f"coupling_comparison_{factor}x_history.png"
        result = subprocess.run(
            [
                "gnuplot",
                "-e",
                f"data='{OUTPUT_ROOT / history_name}'; output='{output}'",
                str(VERIFICATION / "scripts" / "plotCouplingHistory.gnuplot"),
            ],
            cwd=VERIFICATION,
            text=True,
        )
        if result.returncode:
            fail(f"Could not create the coupling history plot {output}")
        print(f"History plot: {output}")
    return {"rows": rows, "histories": histories}


def coupling_study(args: argparse.Namespace, spec: dict) -> bool:
    """Compare the Robin-Neumann and IQN-ILS transients after the coupling start.

    Both couplings are converged at every time step, but the two interface
    formulations are different spatial discretisations of the interface
    conditions: the converged Robin condition prescribes the pressure gradient
    from the interpolated solid acceleration, whereas the Dirichlet-Neumann
    wall is kinematically exact in the discrete sense. Their difference is
    therefore a discretisation error, not a coupling error. The study checks
    that it decreases under mesh refinement and that on the finest level it
    is smaller than the discretisation error of either solution, estimated
    from the change in the IQN-ILS history between the two finest levels.
    """
    end_time = args.end_time or spec["coupling"]["endTime"]
    factors = (
        [int(level) for level in args.levels.split(",")]
        if args.levels else spec["coupling"]["levels"]
    )
    if len(factors) < 2:
        fail("The coupling study needs at least two mesh levels")
    results = {factor: coupling_level(args, spec, factor, end_time) for factor in factors}
    rows = [row for factor in factors for row in results[factor]["rows"]]
    csv_name = "coupling_comparison.csv"
    if args.quick:
        write_csv(OUTPUT_ROOT / csv_name, rows)
        print(f"Quick coupling run completed; results: {OUTPUT_ROOT / csv_name}")
        return True

    coarse, fine = factors[-2], factors[-1]
    fine_robin = results[fine]["rows"][1]
    coarse_robin = results[coarse]["rows"][1]
    failures = []
    for quantity in QUANTITIES:
        # Discretisation error estimate: the change in the IQN-ILS history
        # between the two finest levels, at the times both levels share
        time_c, values_c = results[coarse]["histories"]["iqnils"][quantity]
        time_f, values_f = results[fine]["histories"]["iqnils"][quantity]
        fine_values = {round(t, 9): v for t, v in zip(time_f, values_f)}
        pairs = [(v, fine_values[round(t, 9)]) for t, v in zip(time_c, values_c)
                 if round(t, 9) in fine_values]
        if not pairs:
            fail(f"The {quantity} histories of the two mesh levels share no sample times")
        scale = max(abs(b) for _, b in pairs)
        mesh_change = max(abs(a - b) for a, b in pairs) / scale
        fine_robin[f"{quantity}_iqnils_mesh_change"] = mesh_change
        gap_coarse = coarse_robin[f"{quantity}_vs_iqnils"]
        gap_fine = fine_robin[f"{quantity}_vs_iqnils"]
        if gap_fine >= gap_coarse:
            failures.append(
                f"{quantity}: the Robin vs IQN-ILS difference does not decrease with "
                f"refinement ({100 * gap_coarse:.3f}% on {coarse}x, {100 * gap_fine:.3f}% on {fine}x)"
            )
        if gap_fine >= mesh_change:
            failures.append(
                f"{quantity}: the Robin vs IQN-ILS difference on the {fine}x mesh "
                f"({100 * gap_fine:.3f}%) is not below the IQN-ILS change between the "
                f"{coarse}x and {fine}x meshes ({100 * mesh_change:.3f}%)"
            )
    passed = not failures
    for row in rows:
        row["pass"] = passed
    write_csv(OUTPUT_ROOT / csv_name, rows)

    notes = list(failures)
    for row in rows:
        notes.append(
            f"{row['case']}: {row['total_outer_iterations']:.0f} FSI iterations over "
            f"{row['coupled_time_steps']:.0f} coupled time steps "
            f"(mean {row['total_outer_iterations'] / row['coupled_time_steps']:.1f}, "
            f"maximum {row['maximum_outer_iterations']:.0f} per step)"
        )
    for factor in factors:
        robin = results[factor]["rows"][1]
        notes.append(
            f"Worst converged Robin residuals on the {factor}x mesh: pressure "
            f"{robin['maximum_final_pressure_residual']:.3g}, leakage flux "
            f"{robin['maximum_final_flux_residual']:.3g}"
        )
    header = "| Quantity | " + " | ".join(
        f"Robin vs IQN-ILS, {factor}x (%)" for factor in factors
    ) + f" | IQN-ILS change {coarse}x to {fine}x (%) |"
    table = header + "\n|---|" + "---:|" * (len(factors) + 1) + "\n" + "\n".join(
        f"| {quantity} | " + " | ".join(
            f"{100 * results[factor]['rows'][1][f'{quantity}_vs_iqnils']:.3f}" for factor in factors
        ) + f" | {100 * fine_robin[f'{quantity}_iqnils_mesh_change']:.3f} |"
        for quantity in QUANTITIES
    )
    write_summary(
        f"coupling comparison to t = {end_time:g} s on the "
        + ", ".join(f"{factor}x" for factor in factors) + " meshes",
        passed, notes, table, csv_name,
    )
    print(f"coupling comparison: {'PASS' if passed else 'FAIL'}; results: {OUTPUT_ROOT / csv_name}")
    for message in failures:
        print(f"  - {message}")
    return passed


# ---------------------------------------------------------------------------
# FSI1: steady benchmark
# ---------------------------------------------------------------------------

def configure_fsi1(case: Path, fsi1: dict) -> None:
    """Turn a configured FSI3 copy into the steady FSI1 benchmark.

    FSI1 has the same geometry and fluid as FSI3 but a mean inflow of
    0.2 m/s and a softer plate. The inflow is ramped smoothly, the coupling
    starts early in the ramp, and the run is a pseudo-transient route to the
    steady state. The full solid state is written at the end time so that the
    steady coupling comparison can restart from it.
    """
    for condition in ("dirichletNeumann", "robin"):
        path = case / f"0/fluid/U.{condition}"
        replace_entry(path, "maxValue", f"{fsi1['inletMaxValue']:.8g}")
        replace_entry(path, "transitionPeriod", f"{fsi1['transitionPeriod']:.8g}")
    replace_entry(
        case / "constant/solid/mechanicalProperties",
        "E",
        f"E [1 -1 -2 0 0 0 0] {fsi1['youngsModulus']:.8g}",
    )
    for coupling in ("iqnils", "robin"):
        replace_entry(
            case / f"constant/fsiProperties.{coupling}",
            "couplingStartTime",
            f"{fsi1['couplingStartTime']:.8g}",
        )
    # The FSI1 loads and displacements are one to two orders of magnitude
    # smaller than those of FSI3, so the absolute fluid solver tolerances are
    # tightened accordingly; otherwise the solver residuals set a floor on
    # the interface residual above outerCorrTolerance on the finer meshes
    path = case / "system/fluid/fvSolution"
    text, count = re.subn(r"^(\s*tolerance\s+)[^;]+;",
                          rf"\g<1>{fsi1['fluidSolverTolerance']:.8g};",
                          path.read_text(), flags=re.MULTILINE)
    if count == 0:
        fail(f"No solver tolerances found in {path}")
    path.write_text(text)
    replace_entry(case / "system/controlDict", "writePrecision", "12")
    enable_restart(case)


def set_robin_conditions(directory: Path) -> None:
    """Switch the fluid plate conditions of a time directory to Robin."""
    for field, condition in (("p", "elasticWallPressure"),
                             ("U", "elasticWallVelocity")):
        path = directory / "fluid" / field
        text, count = re.subn(
            r"(\bplate\s*\{\s*type\s+)\w+;",
            rf"\g<1>{condition};",
            path.read_text(),
        )
        if count != 1:
            fail(f"Could not set the Robin condition on plate in {path}")
        path.write_text(text)


def steady_value(time: list[float], values: list[float],
                 window_fraction: float) -> tuple[float, float]:
    """Final value and its relative spread over the closing part of the run."""
    start = time[-1] * (1.0 - window_fraction)
    window = [v for t, v in zip(time, values) if t >= start]
    if len(window) < 2:
        fail("The steady-state window contains too few samples")
    final = values[-1]
    return final, (max(window) - min(window)) / max(abs(final), 1e-30)


def execution_time(case: Path) -> float:
    matches = re.findall(r"ExecutionTime\s*=\s*([0-9.eE+-]+)\s*s",
                         (case / "log.solids4Foam").read_text(errors="replace"))
    return float(matches[-1]) if matches else float("nan")


def fsi1_run(args: argparse.Namespace, spec: dict, factor: int,
             coupling: str) -> dict:
    """Run, or reuse, one FSI1 level and extract its steady values."""
    fsi1 = spec["fsi1"]
    mesh = fsi1["mesh"]
    if factor not in mesh["refinementFactors"]:
        fail(f"Unsupported refinement factor {factor}; choose from {mesh['refinementFactors']}")
    index = mesh["refinementFactors"].index(factor)
    delta_t = mesh["deltaTs"][index]
    end_time = args.end_time or (fsi1["quick"]["endTime"] if args.quick else mesh["endTime"])
    if args.cores == "auto":
        cores = int(fsi1["cores"].get(str(factor), fsi1["cores"]["default"]))
    else:
        cores = int(args.cores)
    name = f"fsi1_{coupling}_mesh_{factor}x"
    case = WORK_ROOT / name
    if not (args.reuse and case.is_dir()):
        case = copy_case(name)
        configure_case(case, spec, factor, delta_t, end_time,
                       args.write_interval, cores)
        configure_fsi1(case, fsi1)
    run_case(case, name, cores, coupling, args.reuse)
    # A reused run keeps the rank count it was run with
    cores = len(list(case.glob("processor*"))) or 1
    row: dict = {"case": name, "coupling": coupling, "mesh_level": index + 1,
                 "refinement": factor, "delta_t": delta_t, "cores": cores,
                 "end_time": end_time, "cell_count": cell_count(case),
                 "execution_time": execution_time(case)}
    if coupling == "robin":
        row.update(robin_residual_summary(case, end_time))
    else:
        row.update(iqnils_iteration_summary(case))
    row["mean_outer_iterations"] = row["total_outer_iterations"] / row["coupled_time_steps"]
    if args.quick:
        return row
    fraction = fsi1["steadyState"]["windowFraction"]
    for quantity, (time, values) in fsi1_histories(case, spec).items():
        if not math.isclose(time[-1], end_time, rel_tol=1e-6, abs_tol=1e-9):
            fail(f"The {quantity} history of {name} ends at t={time[-1]:g}, not at t={end_time:g}")
        row[quantity], row[f"{quantity}_spread"] = steady_value(time, values, fraction)
    return row


def fsi1_histories(case: Path, spec: dict) -> dict[str, tuple[list[float], list[float]]]:
    time_d, ux, uy = displacement_history(case)
    time_f, drag, lift = force_history(case, spec["thickness_m"])
    return {"ux": (time_d, ux), "uy": (time_d, uy),
            "drag": (time_f, drag), "lift": (time_f, lift)}


def steady_failures(row: dict, tolerance: float) -> list[str]:
    return [
        f"{row['case']}: {quantity} has not settled: relative spread "
        f"{row[f'{quantity}_spread']:.2e} over the closing part of the run "
        f"(tolerance {tolerance:g})"
        for quantity in QUANTITIES if row[f"{quantity}_spread"] > tolerance
    ]


def fsi1_format(quantity: str, value: float) -> str:
    return f"{1000 * value:.6f}" if quantity in ("ux", "uy") else f"{value:.4f}"


def fsi1_mesh_study(args: argparse.Namespace, spec: dict) -> bool:
    fsi1 = spec["fsi1"]
    factors = ([int(level) for level in args.levels.split(",")] if args.levels
               else fsi1["mesh"]["defaultLevels"])
    rows = [fsi1_run(args, spec, factor, args.coupling) for factor in factors]
    csv_name = f"fsi1_{args.coupling}_mesh_sweep.csv"
    if args.quick:
        write_csv(OUTPUT_ROOT / csv_name, rows)
        print(f"Quick FSI1 run completed for levels {factors}; results: {OUTPUT_ROOT / csv_name}")
        return True

    failures = []
    tolerance = fsi1["steadyState"]["relativeTolerance"]
    for row in rows:
        failures += steady_failures(row, tolerance)
        for quantity, reference in fsi1["references"].items():
            error = relative_error(row[quantity], reference["value"])
            row[f"{quantity}_error"] = error
            # Tolerances are given for the finest level of the sweep
            limit = reference["tolerance"].get(str(row["refinement"]))
            if (row is rows[-1] and reference["primary"] and limit is not None
                    and error > limit):
                failures.append(
                    f"{row['case']}: {quantity} = {row[quantity]:.6g} differs from the "
                    f"Featflow value {reference['value']:.6g} by {100 * error:.3f}% "
                    f"(tolerance {100 * limit:g}%)"
                )
        if row is rows[-1] and str(row["refinement"]) not in fsi1["references"]["uy"]["tolerance"]:
            failures.append(f"{row['case']}: no reference tolerances for the finest level; "
                            "include level 2 or 4 in the sweep")
    if len(rows) > 1:
        for quantity, reference in fsi1["references"].items():
            growth = rows[-1][f"{quantity}_error"] - rows[0][f"{quantity}_error"]
            if reference["primary"] and growth > fsi1["mesh"]["maximumErrorGrowth"]:
                failures.append(
                    f"{quantity} reference error grew from {100 * rows[0][f'{quantity}_error']:.3f}% "
                    f"to {100 * rows[-1][f'{quantity}_error']:.3f}% between the coarsest and finest levels"
                )
    passed = not failures
    for row in rows:
        row["pass"] = passed
    write_csv(OUTPUT_ROOT / csv_name, rows)

    finest = fsi1["featflowLevels"][-1]
    lines = [
        "| Quantity | " + " | ".join(f"{row['refinement']}x" for row in rows)
        + f" | Featflow {finest['level']} |",
        "|---|" + "---:|" * (len(rows) + 1),
    ]
    units = {"ux": "mm", "uy": "mm", "drag": "N/m", "lift": "N/m"}
    for quantity in QUANTITIES:
        cells = [f"{fsi1_format(quantity, row[quantity])} "
                 f"({100 * row[f'{quantity}_error']:.2f}%)" for row in rows]
        cells.append(fsi1_format(quantity, finest[quantity]))
        lines.append(f"| {quantity} ({units[quantity]}) | " + " | ".join(cells) + " |")
    notes = list(failures)
    for row in rows:
        notes.append(
            f"{row['case']}: {row['cell_count']} cells, dt = {row['delta_t']:g} s, "
            f"{row['cores']} core(s), {row['execution_time']:.0f} s, mean "
            f"{row['mean_outer_iterations']:.2f} FSI iterations per step, largest "
            f"steady spread {max(row[f'{q}_spread'] for q in QUANTITIES):.2e}"
        )
    write_summary(f"FSI1 {args.coupling} mesh study", passed, notes,
                  "\n".join(lines), csv_name)
    print(f"FSI1 {args.coupling} mesh study: {'PASS' if passed else 'FAIL'}; "
          f"results: {OUTPUT_ROOT / csv_name}")
    for message in failures:
        print(f"  - {message}")
    return passed


def fsi1_continuation(args: argparse.Namespace, spec: dict, base: dict,
                      coupling: str) -> dict:
    """Continue the steady IQN-ILS state with one coupling at a small dt."""
    fsi1 = spec["fsi1"]
    factor = base["refinement"]
    start = base["end_time"]
    start_directory = WORK_ROOT / base["case"] / f"{start:g}"
    if not start_directory.is_dir():
        fail(f"The IQN-ILS run did not write its end state {start_directory}")
    delta_t = fsi1["coupling"]["deltaT"]
    end_time = start + (fsi1["quick"]["continuation"] if args.quick
                        else fsi1["coupling"]["continuation"])
    name = f"fsi1_{coupling}_continuation_{factor}x"
    case = WORK_ROOT / name
    if not (args.reuse and case.is_dir()):
        case = copy_case(name)
        configure_case(case, spec, factor, delta_t, end_time, None, 1)
        configure_fsi1(case, fsi1)
        shutil.copytree(start_directory, case / start_directory.name)
        if coupling == "robin":
            set_robin_conditions(case / start_directory.name)
    run_case(case, name, 1, coupling, args.reuse)
    row: dict = {"case": name, "coupling": coupling, "refinement": factor,
                 "delta_t": delta_t, "cores": 1, "start_time": start,
                 "end_time": end_time, "cell_count": cell_count(case),
                 "execution_time": execution_time(case)}
    if coupling == "robin":
        row.update(robin_residual_summary(case, end_time))
    else:
        row.update(iqnils_iteration_summary(case))
    row["mean_outer_iterations"] = row["total_outer_iterations"] / row["coupled_time_steps"]
    return row


def fsi1_coupling_study(args: argparse.Namespace, spec: dict) -> bool:
    """Compare the couplings from one converged FSI1 steady state.

    The steady solution does not depend on the path taken to it, so both
    couplings restart from the steady state of the IQN-ILS mesh-study run.
    The Robin-Neumann fixed-point iterations converge far too slowly at the
    pseudo-transient time step, so both continuations use a smaller time
    step; the change of time step perturbs the state slightly, identically
    for both couplings. The differences are reported only: the study fails
    only if the IQN-ILS start state has not settled or a run does not
    converge.
    """
    fsi1 = spec["fsi1"]
    factor = int(args.levels) if args.levels else fsi1["coupling"]["refinement"]
    # The restarts are run in serial so that the start state does not have to
    # be decomposed
    if args.cores != "auto" and int(args.cores) != 1:
        fail("The FSI1 coupling study runs in serial")
    args.cores = "1"
    base = fsi1_run(args, spec, factor, "iqnils")
    rows = [base] + [fsi1_continuation(args, spec, base, coupling)
                     for coupling in ("iqnils", "robin")]
    csv_name = f"fsi1_coupling_comparison_{factor}x.csv"
    if args.quick:
        write_csv(OUTPUT_ROOT / csv_name, rows)
        print(f"Quick FSI1 coupling run completed; results: {OUTPUT_ROOT / csv_name}")
        return True

    _, iqnils, robin = rows
    failures = steady_failures(base, fsi1["steadyState"]["relativeTolerance"])
    histories = {row["coupling"]: fsi1_histories(WORK_ROOT / row["case"], spec)
                 for row in (iqnils, robin)}
    for quantity in QUANTITIES:
        selected = {}
        for coupling, history in histories.items():
            time, values = history[quantity]
            selected[coupling] = [(t, v) for t, v in zip(time, values)
                                  if t > base["end_time"] + 1e-9]
        time_i = [t for t, _ in selected["iqnils"]]
        time_r = [t for t, _ in selected["robin"]]
        if not time_i or len(time_i) != len(time_r) or any(
            not math.isclose(a, b, rel_tol=1e-9, abs_tol=1e-12) for a, b in zip(time_i, time_r)
        ):
            fail(f"The {quantity} continuations are sampled at different times")
        values_i = [v for _, v in selected["iqnils"]]
        values_r = [v for _, v in selected["robin"]]
        # Each coupling converges the interface to outerCorrTolerance at
        # every step, which leaves step-to-step noise in the small resultants
        # (lift, ux), so the continued states are compared through their
        # means over the closing half of the continuation. Differences are
        # relative to the steady value of the quantity.
        scale = abs(base[quantity])
        closing = time_i[0] + 0.5 * (time_i[-1] - time_i[0])
        means = {}
        for coupling, values in (("iqnils", values_i), ("robin", values_r)):
            window = [v for t, v in zip(time_i, values) if t >= closing]
            means[coupling] = sum(window) / len(window)
            row = iqnils if coupling == "iqnils" else robin
            row[quantity] = means[coupling]
            row[f"{quantity}_noise"] = (max(window) - min(window)) / scale
        difference = abs(means["robin"] - means["iqnils"]) / scale
        robin[f"{quantity}_vs_iqnils"] = difference
        robin[f"{quantity}_max_vs_iqnils"] = max(
            abs(a - b) for a, b in zip(values_i, values_r)) / scale
        robin[f"{quantity}_first_vs_iqnils"] = abs(values_r[0] - values_i[0]) / scale
    passed = not failures
    for row in rows:
        row["pass"] = passed
    write_csv(OUTPUT_ROOT / csv_name, rows)

    table = [f"| Quantity | Steady (dt = {base['delta_t']:g} s) | IQN-ILS mean | "
             "Robin mean | Mean difference | First-step difference | "
             "IQN-ILS noise | Robin noise |",
             "|---|---:|---:|---:|---:|---:|---:|---:|"]
    for quantity in QUANTITIES:
        table.append(
            f"| {quantity} | {base[quantity]:.7g} | {iqnils[quantity]:.7g} | "
            f"{robin[quantity]:.7g} | {100 * robin[f'{quantity}_vs_iqnils']:.4f}% | "
            f"{100 * robin[f'{quantity}_first_vs_iqnils']:.4f}% | "
            f"{100 * iqnils[f'{quantity}_noise']:.3f}% | "
            f"{100 * robin[f'{quantity}_noise']:.3f}% |"
        )
    notes = list(failures)
    notes.append(
        f"Both couplings restarted from the IQN-ILS steady state at t = {base['end_time']:g} s "
        f"and run to t = {robin['end_time']:g} s with dt = {robin['delta_t']:g} s; "
        "differences are relative to the steady value"
    )
    for row in rows:
        notes.append(
            f"{row['case']}: {row['total_outer_iterations']:.0f} FSI iterations over "
            f"{row['coupled_time_steps']:.0f} coupled time steps (mean "
            f"{row['mean_outer_iterations']:.2f}, maximum "
            f"{row['maximum_outer_iterations']:.0f} per step), {row['execution_time']:.0f} s"
        )
    notes.append(
        "Worst converged Robin pressure residual: "
        f"{robin['maximum_final_pressure_residual']:.3g}; leakage-flux residual: "
        f"{robin['maximum_final_flux_residual']:.3g}"
    )
    notes.append("The Robin vs IQN-ILS differences are reported only, not checked")
    write_summary(f"FSI1 steady coupling comparison on the {factor}x mesh",
                  passed, notes, "\n".join(table), csv_name)
    print(f"FSI1 coupling comparison (differences reported only): "
          f"{'PASS' if passed else 'FAIL'}; results: {OUTPUT_ROOT / csv_name}")
    for message in failures:
        print(f"  - {message}")
    return passed


# ---------------------------------------------------------------------------
# FSI2: periodic benchmark with large deformation
# ---------------------------------------------------------------------------

def fsi2_spec(spec: dict) -> dict:
    """The FSI3 specification with the FSI2 entries in place of the FSI3 ones."""
    fsi2 = spec["fsi2"]
    merged = dict(spec)
    for key in ("history_reference", "outer_corr_tolerance", "mesh", "quick",
                "cores", "periodicity", "references", "publishedReferences",
                "references_source"):
        merged[key] = fsi2[key]
    merged["benchmark"] = "fsi2"
    return merged


def configure_fsi2(case: Path, fsi2: dict) -> None:
    """Turn a configured FSI3 copy into the periodic FSI2 benchmark.

    FSI2 has the same geometry and fluid as FSI3, half its mean inflow
    (1 m/s) and a heavier, softer plate. The inflow is ramped smoothly over
    the benchmark's 2 s ramp and the coupling starts at the end of the ramp.
    """
    for condition in ("dirichletNeumann", "robin"):
        path = case / f"0/fluid/U.{condition}"
        replace_entry(path, "maxValue", f"{fsi2['inletMaxValue']:.8g}")
        replace_entry(path, "transitionPeriod", f"{fsi2['transitionPeriod']:.8g}")
    properties = case / "constant/solid/mechanicalProperties"
    replace_entry(properties, "rho", f"rho [1 -3 0 0 0 0 0] {fsi2['solidDensity']:.8g}")
    replace_entry(properties, "E", f"E [1 -1 -2 0 0 0 0] {fsi2['youngsModulus']:.8g}")
    for coupling in ("iqnils", "robin"):
        replace_entry(
            case / f"constant/fsiProperties.{coupling}",
            "couplingStartTime",
            f"{fsi2['couplingStartTime']:.8g}",
        )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--benchmark", choices=("fsi3", "fsi1", "fsi2"), default="fsi3",
                        help="periodic FSI3 (default), steady FSI1 or periodic FSI2 benchmark")
    parser.add_argument("--study", choices=("mesh", "coupling"), default="mesh")
    parser.add_argument("--coupling", choices=("iqnils", "robin"), default="iqnils",
                        help="coupling variant for the mesh study (default: iqnils)")
    parser.add_argument("--levels", help="comma-separated refinement factors, e.g. 1,2,4 "
                        "(mesh study) or a single factor (coupling study)")
    parser.add_argument("--cores", default="auto",
                        help="MPI ranks per case: positive integer or auto (default)")
    parser.add_argument("--end-time", type=float,
                        help="override the end time in seconds (mesh study default 7, coupling study 2.3)")
    parser.add_argument("--window", type=float,
                        help="length of the closing analysis window in seconds")
    parser.add_argument("--write-interval", type=int,
                        help="write fields every N time steps (default: final time only)")
    parser.add_argument("--reuse", action="store_true",
                        help="reuse completed runs found under verification/work")
    parser.add_argument("--quick", action="store_true",
                        help="short smoke run that only checks the cases complete")
    args = parser.parse_args()
    if args.cores != "auto" and (not args.cores.isdecimal() or int(args.cores) < 1):
        parser.error("--cores must be a positive integer or auto")
    if args.write_interval is not None and args.write_interval < 1:
        parser.error("--write-interval must be a positive integer")
    for executable in ("blockMesh", "solids4Foam"):
        if not command_exists(executable):
            fail(f"Required executable '{executable}' is unavailable. "
                 "Source OpenFOAM and build solids4foam first.")
    if not (TUTORIAL / "Allrun").is_file():
        fail(f"Tutorial not found at {TUTORIAL}")
    spec = json.loads(REFERENCE_FILE.read_text())
    OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)
    (OUTPUT_ROOT / "verification_summary.md").write_text(
        "# HronTurek verification summary\n\n"
    )
    if args.benchmark == "fsi1":
        study = fsi1_coupling_study if args.study == "coupling" else fsi1_mesh_study
        return 0 if study(args, spec) else 1
    if args.benchmark == "fsi2":
        if args.study == "coupling":
            fail("The coupling study is not available for FSI2")
        return 0 if mesh_study(args, fsi2_spec(spec)) else 1
    if args.study == "coupling":
        return 0 if coupling_study(args, spec) else 1
    return 0 if mesh_study(args, spec) else 1


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (RuntimeError, subprocess.SubprocessError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        sys.exit(2)
