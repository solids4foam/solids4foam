#!/usr/bin/env python3
"""Run the mokLidDrivenCavity verification studies.

The lid-driven cavity with a flexible bottom (Mok 2001, as defined by Valdes
2007) is run on copies of the tutorial under verification/work. Three studies
are available:

- mesh:     the fluid and solid meshes are refined together;
- timestep: the time step is refined on a fixed mesh;
- solid:    the standard and the high-order solid discretisations are compared
            on the same meshes, including a through-thickness refinement.

The midpoint displacement history of each run is compared with the group B
references, whose boundary conditions are fully specified: Valdes (2007),
Kratos (Zorrilla) and the scalar values of Tiba et al. (2026). The late-time
peak, trough and mean of the periodic response must lie within a tolerance of
the envelope spanned by these references, and the history must stay close to
the band between the Valdes and Kratos curves. The curves of Mok (2001), Wall
(1999), Gerbeau and Vidrascu (2003) and Kassiotis et al. (2011) are reported
for context only: their openings are defined differently.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import re
import shutil
import subprocess
import sys
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent
VERIFY_DIR = SCRIPT_DIR.parent
CASE_DIR = VERIFY_DIR.parent
WORK_DIR = VERIFY_DIR / "work"
POST_DIR = VERIFY_DIR / "postProcessing"
REFERENCE_DIR = VERIFY_DIR / "reference"
REFERENCE_FILE = REFERENCE_DIR / "mokLidDrivenCavity_verification_references.json"
STUDIES = ("mesh", "timestep", "solid")


# --------------------------------------------------------------------------- #
# Command line
# --------------------------------------------------------------------------- #

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run the mokLidDrivenCavity verification studies"
    )
    parser.add_argument(
        "--study",
        default="all",
        help="comma-separated studies from: mesh, timestep, solid, all",
    )
    parser.add_argument(
        "--quick",
        action="store_true",
        help="smoke test: the two coarsest members of each study, to t = 10 s",
    )
    parser.add_argument(
        "--reuse",
        action="store_true",
        help="reuse completed cases under verification/work",
    )
    parser.add_argument(
        "--cores",
        type=int,
        default=1,
        help="MPI ranks per case (default 1, serial)",
    )
    parser.add_argument(
        "--end-time",
        type=float,
        help="override the end time (s); the checks need at least 30 s",
    )
    parser.add_argument(
        "--keep-going",
        action="store_true",
        help="continue after an individual case fails",
    )
    return parser.parse_args()


# --------------------------------------------------------------------------- #
# Case set-up
# --------------------------------------------------------------------------- #

def ignored(directory: str, names: list[str]) -> set[str]:
    """Copy only the inputs a fresh run needs, never a previous result."""
    ignored_names = {
        "verification", "regressionTests", "postProcessing", "images",
        "dynamicCode", "README.md",
    }
    ignored_names.update(name for name in names if name.startswith("processor"))
    ignored_names.update(name for name in names if name.startswith("log."))
    directory_path = Path(directory)
    if directory_path == CASE_DIR:
        ignored_names.update(
            name
            for name in names
            if name != "0" and re.fullmatch(r"[0-9]+(?:\.[0-9]+)?(?:e[-+]?\d+)?", name)
        )
    if directory_path.name in {"fluid", "solid"} and directory_path.parent.name == "constant":
        ignored_names.add("polyMesh")
    return ignored_names.intersection(names)


def set_blocks(path: Path, divisions: list[tuple[int, int]]) -> None:
    """Replace the in-plane divisions of each hex block, in order."""
    pattern = re.compile(r"(hex\s*\([^)]*\)\s*)\(\s*\d+\s+\d+\s+(\d+)\s*\)")
    text = path.read_text()
    matches = list(pattern.finditer(text))
    if len(matches) != len(divisions):
        raise RuntimeError(
            f"expected {len(divisions)} blocks in {path}, found {len(matches)}"
        )
    for match, (nx, ny) in reversed(list(zip(matches, divisions))):
        text = (
            text[: match.start()]
            + f"{match.group(1)}({nx} {ny} {match.group(2)})"
            + text[match.end():]
        )
    path.write_text(text)


def set_mesh(run_dir: Path, nx: int, solid_ny: int) -> None:
    """nx cells across the cavity and along the membrane (directMap needs
    matching interface faces), 7/8 of them below the openings, and solid_ny
    cells through the membrane thickness."""
    if nx % 8:
        raise RuntimeError("nx must be a multiple of 8 to resolve the openings")
    set_blocks(
        run_dir / "system" / "fluid" / "blockMeshDict",
        [(nx, 7 * nx // 8), (nx, nx // 8)],
    )
    set_blocks(run_dir / "system" / "solid" / "blockMeshDict", [(nx, solid_ny)])


def set_entry(path: Path, key: str, value: str) -> None:
    text, count = re.subn(
        rf"(^\s*{re.escape(key)}\s+)[^;]+;",
        rf"\g<1>{value};",
        path.read_text(),
        count=1,
        flags=re.MULTILINE,
    )
    if count != 1:
        raise RuntimeError(f"could not set {key} in {path}")
    path.write_text(text)


def set_decomposition(run_dir: Path, cores: int) -> None:
    for path in (
        run_dir / "system" / "decomposeParDict",
        run_dir / "system" / "fluid" / "decomposeParDict",
        run_dir / "system" / "solid" / "decomposeParDict",
    ):
        # The region dictionaries use scotch, which needs no further input
        set_entry(path, "numberOfSubdomains", str(cores))


def prepare_case(run_dir: Path, member: dict, end_time: float, cores: int) -> None:
    if run_dir.exists():
        shutil.rmtree(run_dir)
    shutil.copytree(CASE_DIR, run_dir, ignore=ignored, symlinks=True)
    set_mesh(run_dir, int(member["nx"]), int(member["solid_ny"]))
    control = run_dir / "system" / "controlDict"
    delta_t = float(member["delta_t"])
    set_entry(control, "deltaT", f"{delta_t:.10g}")
    set_entry(control, "endTime", f"{end_time:.10g}")
    # Write fields only at the end: the studies use function-object data
    set_entry(control, "writeInterval", str(int(round(end_time / delta_t))))
    if cores > 1:
        set_decomposition(run_dir, cores)


# --------------------------------------------------------------------------- #
# Running and reading results
# --------------------------------------------------------------------------- #

def numeric_rows(path: Path) -> list[list[float]]:
    rows = []
    for line in path.read_text(errors="replace").splitlines():
        fields = line.replace("(", " ").replace(")", " ").split()
        if not fields or fields[0].startswith("#"):
            continue
        try:
            rows.append([float(field) for field in fields])
        except ValueError:
            continue
    return rows


def read_history(run_dir: Path) -> list[tuple[float, float]]:
    candidates = sorted(
        run_dir.glob("postProcessing/**/solidPointDisplacement_midpoint.dat")
    )
    if not candidates:
        raise RuntimeError(f"no midpoint displacement data in {run_dir}")
    history = [
        (row[0], row[2]) for row in numeric_rows(candidates[-1]) if len(row) >= 3
    ]
    if not all(math.isfinite(time) and math.isfinite(value)
               for time, value in history):
        raise RuntimeError(f"non-finite midpoint displacement in {candidates[-1]}")
    return history


def check_history(
    history: list[tuple[float, float]], end: float, delta_t: float, label: str
) -> None:
    """The history must be increasing in time, without gaps, and reach the end
    of the comparison window."""
    if not history:
        raise RuntimeError(f"{label}: empty midpoint displacement history")
    times = [time for time, _ in history]
    if any(right <= left for left, right in zip(times, times[1:])):
        raise RuntimeError(f"{label}: midpoint history times are not increasing")
    gap = max(
        [right - left for left, right in zip(times, times[1:])] + [times[0]]
    )
    if gap > 1.5 * delta_t:
        raise RuntimeError(
            f"{label}: midpoint history has a gap of {gap:g} s "
            f"(time step {delta_t:g} s)"
        )
    if times[-1] < end - 0.5 * delta_t:
        raise RuntimeError(
            f"{label}: midpoint history stops at t = {times[-1]:g} s, before "
            f"the end of the comparison window at {end:g} s"
        )


def read_iterations(run_dir: Path) -> list[int]:
    """FSI iterations per time step, from the coupling residual file."""
    path = run_dir / "postProcessing" / "fsiResiduals.dat"
    if not path.is_file():
        return []
    last: dict[float, int] = {}
    for row in numeric_rows(path):
        if len(row) >= 3:
            last[row[0]] = int(row[1])
    return list(last.values())


def check_solver_log(solver_log: Path, end_time: float, label: str) -> None:
    text = solver_log.read_text(errors="replace")
    if re.search(r"FOAM FATAL|FOAM aborting|^ERROR$|\[stack trace\]", text,
                 re.MULTILINE):
        raise RuntimeError(f"{label} diverged or aborted; see {solver_log}")
    if not re.search(r"^End\s*$", text, re.MULTILINE):
        raise RuntimeError(
            f"{label} did not terminate normally (no End); see {solver_log}"
        )
    times = [float(value) for value in
             re.findall(r"^Time = ([0-9.eE+-]+)\s*$", text, re.MULTILINE)]
    if not times or times[-1] < end_time - 1.0e-8:
        reached = times[-1] if times else 0.0
        raise RuntimeError(
            f"{label} stopped at t = {reached:g} before the end time "
            f"{end_time:g}; see {solver_log}"
        )


def tutorial_fingerprint() -> str:
    """Hash of the tutorial inputs that a member copies, so that a member run
    from different inputs is not reused."""
    digest = hashlib.sha256()
    for path in sorted(CASE_DIR.rglob("*")):
        relative = path.relative_to(CASE_DIR)
        if (
            relative.parts[0] not in {"0", "constant", "system", "Allrun"}
            or "polyMesh" in relative.parts
            or not path.is_file()
        ):
            continue
        digest.update(str(relative).encode())
        digest.update(path.read_bytes())
    return digest.hexdigest()


def member_metadata(member: dict, end_time: float) -> dict:
    return {
        "member": {key: member[key] for key in sorted(member)},
        "end_time_s": end_time,
        "tutorial": tutorial_fingerprint(),
    }


def run_member(
    study: str, member: dict, end_time: float, cores: int, reuse: bool
) -> dict:
    run_dir = WORK_DIR / study / member["name"]
    solver_log = run_dir / "log.solids4Foam"
    metadata_file = run_dir / "verification_member.json"
    label = f"{study}/{member['name']}"
    metadata = member_metadata(member, end_time)
    completed_case = False
    if reuse and solver_log.exists() and metadata_file.is_file():
        stored = json.loads(metadata_file.read_text())
        completed_case = (
            all(stored.get(key) == value for key, value in metadata.items())
            and re.search(r"^End\s*$", solver_log.read_text(errors="replace"),
                          re.MULTILINE) is not None
        )
        if completed_case:
            cores = int(stored.get("cores", cores))
    if completed_case:
        print(f"Reusing {label} in {run_dir}", flush=True)
    else:
        prepare_case(run_dir, member, end_time, cores)
        arguments = ["./Allrun", member.get("coupling", "iqnils"), member["solid"]]
        if cores > 1:
            arguments.append("parallel")
        print(
            f"Running {label}: nx = {member['nx']}, solid_ny = "
            f"{member['solid_ny']}, dt = {member['delta_t']}, "
            f"{member['solid']} solid, {cores} rank(s)",
            flush=True,
        )
        with (run_dir / "log.Allverify").open("w") as log:
            completed = subprocess.run(
                arguments, cwd=run_dir, stdout=log, stderr=subprocess.STDOUT,
                check=False,
            )
        if completed.returncode:
            raise RuntimeError(f"{label} failed; see {run_dir / 'log.Allverify'}")
        if not solver_log.exists():
            if "Skipping this case as PETSc is not installed" in (
                run_dir / "log.Allverify"
            ).read_text(errors="replace"):
                raise RuntimeError("this case requires solids4foam built with PETSc")
            raise RuntimeError(f"solver log was not created in {run_dir}")
        check_solver_log(solver_log, end_time, label)
        metadata_file.write_text(
            json.dumps(dict(metadata, cores=cores), indent=2) + "\n"
        )
    check_solver_log(solver_log, end_time, label)
    log_text = solver_log.read_text(errors="replace")
    clock = re.findall(r"ClockTime\s*=\s*([0-9.eE+-]+)", log_text)
    iterations = read_iterations(run_dir)
    return {
        "study": study,
        "name": member["name"],
        "nx": int(member["nx"]),
        "solid_ny": int(member["solid_ny"]),
        "delta_t_s": float(member["delta_t"]),
        "solid": member["solid"],
        "coupling": member.get("coupling", "iqnils"),
        "cores": cores,
        "history": read_history(run_dir),
        "fsi_iterations_mean": (
            sum(iterations) / len(iterations) if iterations else math.nan
        ),
        "fsi_iterations_max": max(iterations) if iterations else 0,
        "clock_time_s": float(clock[-1]) if clock else math.nan,
    }


# --------------------------------------------------------------------------- #
# Metrics
# --------------------------------------------------------------------------- #

def read_curve(path: Path) -> list[tuple[float, float]]:
    points = []
    for line in path.read_text().splitlines():
        if not line.strip() or line.startswith("#") or line.startswith("time"):
            continue
        time, value = (float(field) for field in line.split(",")[:2])
        points.append((time, value))
    return sorted(points)


def interpolate(curve: list[tuple[float, float]], time: float) -> float:
    if time <= curve[0][0]:
        return curve[0][1]
    if time >= curve[-1][0]:
        return curve[-1][1]
    low, high = 0, len(curve) - 1
    while high - low > 1:
        mid = (low + high) // 2
        if curve[mid][0] <= time:
            low = mid
        else:
            high = mid
    (t0, u0), (t1, u1) = curve[low], curve[high]
    return u0 + (u1 - u0) * (time - t0) / (t1 - t0)


def periodic_stats(
    curve: list[tuple[float, float]], start: float, end: float, period: float
) -> dict:
    """Mean peak and trough of the periodic response, one per lid period,
    the time-mean, and the mean phase of the peaks within the period."""
    peaks, troughs, peak_phases = [], [], []
    cycle_start = start
    while cycle_start + period <= end + 1.0e-9:
        window = [
            point for point in curve
            if cycle_start <= point[0] < cycle_start + period
        ]
        if len(window) >= 10:
            peak = max(window, key=lambda point: point[1])
            peaks.append(peak[1])
            peak_phases.append((peak[0] - cycle_start) % period)
            troughs.append(min(point[1] for point in window))
        cycle_start += period
    if not peaks:
        raise RuntimeError("no complete lid period in the comparison window")
    samples = [
        interpolate(curve, start + (end - start) * index / 2000.0)
        for index in range(2001)
    ]
    return {
        "peak": sum(peaks) / len(peaks),
        "trough": sum(troughs) / len(troughs),
        "mean": sum(samples) / len(samples),
        "peak_phase": sum(peak_phases) / len(peak_phases),
        "cycles": len(peaks),
    }


def history_difference(
    curve: list[tuple[float, float]],
    reference: list[tuple[float, float]],
    start: float,
    end: float,
) -> float:
    end = min(end, reference[-1][0])
    times = [time for time, _ in curve if start <= time <= end]
    if not times:
        return math.nan
    return max(
        abs(interpolate(curve, time) - interpolate(reference, time))
        for time in times
    )


def band_distance(
    curve: list[tuple[float, float]],
    references: list[list[tuple[float, float]]],
    start: float,
    end: float,
) -> float:
    """Largest distance of the curve from the band between the references,
    which at each time runs from the lowest to the highest reference value."""
    end = min([end] + [reference[-1][0] for reference in references])
    distances = []
    for time, value in curve:
        if not start <= time <= end:
            continue
        values = [interpolate(reference, time) for reference in references]
        distances.append(max(value - max(values), min(values) - value, 0.0))
    return max(distances) if distances else math.nan


QUANTITIES = ("peak", "trough", "mean")


def envelope(references: dict) -> dict:
    """Range of each quantity over the group B references that give it."""
    result = {}
    for quantity in QUANTITIES:
        values = [
            reference["stats"][quantity]
            for reference in references.values()
            if reference["group_b"] and quantity in reference["stats"]
        ]
        result[quantity] = (min(values), max(values))
    return result


def envelope_deviation(value: float, low: float, high: float) -> float:
    """Relative distance outside [low, high]: positive above, negative below."""
    if value > high:
        return (value - high) / abs(high)
    if value < low:
        return (value - low) / abs(low)
    return 0.0


def evaluate(row: dict, references: dict, config: dict) -> None:
    window = config["comparison_window_s"]
    period = float(config["lid_period_s"])
    start, end = float(window[0]), float(window[1])
    label = f"{row['study']}/{row['name']}"
    check_history(row["history"], end, row["delta_t_s"], label)
    stats = periodic_stats(row["history"], start, end, period)
    row.update({f"{key}_m" if key in QUANTITIES else key: value
                for key, value in stats.items()})

    ranges = envelope(references)
    for quantity in QUANTITIES:
        low, high = ranges[quantity]
        row[f"{quantity}_envelope"] = (low, high)
        row[f"{quantity}_envelope_deviation"] = envelope_deviation(
            stats[quantity], low, high
        )

    # The band between the group B curves, normalised by the smaller of their
    # peaks
    curves = [
        reference for reference in references.values()
        if reference["group_b"] and reference["curve"]
    ]
    scale = min(abs(reference["stats"]["peak"]) for reference in curves)
    row["band_distance"] = band_distance(
        row["history"], [reference["curve"] for reference in curves], start, end
    ) / scale

    for name, reference in references.items():
        if not reference["curve"]:
            continue
        ref_scale = abs(reference["stats"]["peak"])
        for quantity in QUANTITIES:
            row[f"{name}_{quantity}_error"] = (
                stats[quantity] - reference["stats"][quantity]
            ) / ref_scale
        row[f"{name}_history_difference"] = history_difference(
            row["history"], reference["curve"], start, end
        ) / ref_scale

    checked = [f"{quantity}_m" for quantity in QUANTITIES] + [
        f"{quantity}_envelope_deviation" for quantity in QUANTITIES
    ] + ["band_distance"]
    if not all(math.isfinite(row[key]) for key in checked):
        raise RuntimeError(f"{label}: non-finite verification metric")


# --------------------------------------------------------------------------- #
# Studies
# --------------------------------------------------------------------------- #

def study_members(config: dict, study: str, quick: bool) -> list[dict]:
    members = [dict(member) for member in config["studies"][study]["members"]]
    if quick:
        members = members[:2]
    return members


def check_study(study: str, rows: list[dict], config: dict,
                quick: bool) -> list[tuple[str, bool, str]]:
    """Return (description, passed, detail) for every check of a study."""
    acceptance = config["acceptance"]
    checks: list[tuple[str, bool, str]] = []
    if quick:
        checks.append((f"{study}: all members completed", True, ""))
        return checks

    expected = len(study_members(config, study, quick=False))
    if len(rows) < expected:
        checks.append((
            f"{study}: all {expected} members completed",
            False,
            f"{len(rows)} completed",
        ))
        if len(rows) < 2:
            return checks

    if study in {"mesh", "timestep"}:
        # The change in the periodic response (the larger of the peak and the
        # trough changes) between successive members must shrink, unless it is
        # already below the noise floor set by the FSI and solver tolerances,
        # and the last change must be small compared with the reference
        # tolerances
        changes = [
            max(
                abs(right["peak_m"] - left["peak_m"]),
                abs(right["trough_m"] - left["trough_m"]),
            )
            for left, right in zip(rows, rows[1:])
        ]
        scale = abs(rows[-1]["peak_m"])
        floor = float(acceptance["noise_floor_fraction"]) * scale
        detail = ", ".join(f"{100.0 * change / scale:.2f}%" for change in changes)
        checks.append((
            f"{study}: response change decreases with refinement",
            all(b < a or b < floor for a, b in zip(changes, changes[1:])),
            detail,
        ))
        tolerance = float(acceptance["converged_change_fraction"])
        checks.append((
            f"{study}: finest response change below {100.0 * tolerance:g}% of "
            "the peak",
            changes[-1] / scale < tolerance,
            f"{100.0 * changes[-1] / scale:.2f}%",
        ))
        compared = [rows[-1]]
    else:
        # The solid discretisations and through-thickness resolutions must
        # agree with the last member (the high-order solid) to within the
        # converged-change tolerance
        last = rows[-1]
        scale = abs(last["peak_m"])
        tolerance = float(acceptance["converged_change_fraction"])
        for row in rows[:-1]:
            change = max(
                abs(row["peak_m"] - last["peak_m"]),
                abs(row["trough_m"] - last["trough_m"]),
            )
            checks.append((
                f"{study}: {row['name']} agrees with {last['name']} within "
                f"{100.0 * tolerance:g}% of the peak",
                change / scale < tolerance,
                f"{100.0 * change / scale:.2f}%",
            ))
        compared = [row for row in rows if row.get("compare", True)]

    envelope_tolerance = float(acceptance["envelope_tolerance"])
    band_tolerance = float(acceptance["band_tolerance"])
    for row in compared:
        for quantity in QUANTITIES:
            deviation = row[f"{quantity}_envelope_deviation"]
            low, high = row[f"{quantity}_envelope"]
            value = row[f"{quantity}_m"]
            # Distance to the nearer acceptance limit, relative to its bound
            margin = min(
                ((1.0 + envelope_tolerance) * high - value) / abs(high),
                (value - (1.0 - envelope_tolerance) * low) / abs(low),
            )
            checks.append((
                f"{study}/{row['name']}: {quantity} within "
                f"{100.0 * envelope_tolerance:g}% of the group B envelope "
                f"[{low:.4f}, {high:.4f}] m",
                abs(deviation) <= envelope_tolerance,
                f"{value:.4f} m, outside by {100.0 * deviation:+.2f}%, "
                f"margin to the limit {100.0 * margin:.2f}%",
            ))
        checks.append((
            f"{study}/{row['name']}: history within "
            f"{100.0 * band_tolerance:g}% of the peak of the group B band",
            row["band_distance"] <= band_tolerance,
            f"{100.0 * row['band_distance']:.2f}%, margin "
            f"{100.0 * (band_tolerance - row['band_distance']):.2f}%",
        ))
    return checks


# --------------------------------------------------------------------------- #
# Output
# --------------------------------------------------------------------------- #

SUMMARY_COLUMNS = [
    "study", "name", "solid", "coupling", "nx", "solid_ny", "delta_t_s",
    "cores", "peak_m", "trough_m", "mean_m", "peak_phase", "cycles",
    "fsi_iterations_mean", "fsi_iterations_max", "clock_time_s",
    "peak_envelope_deviation", "trough_envelope_deviation",
    "mean_envelope_deviation", "band_distance",
]


def write_study_csv(study: str, rows: list[dict], curve_names: list[str]) -> None:
    columns = list(SUMMARY_COLUMNS)
    for name in curve_names:
        columns += [f"{name}_{quantity}_error" for quantity in QUANTITIES]
        columns.append(f"{name}_history_difference")
    with (POST_DIR / f"{study}_study.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=columns, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    # One history file per study, sampled at the coarsest time step of the
    # study, for the plots and for inspection
    step = max(row["delta_t_s"] for row in rows)
    end = min(row["history"][-1][0] for row in rows)
    count = int(round(end / step))
    with (POST_DIR / f"{study}_histories.csv").open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["time_s"] + [row["name"] for row in rows])
        for index in range(count + 1):
            time = index * step
            writer.writerow(
                [f"{time:.6g}"]
                + [f"{interpolate(row['history'], time):.6g}" for row in rows]
            )


def create_plot(study: str, rows: list[dict]) -> None:
    plot_script = SCRIPT_DIR / "plotHistories.gnuplot"
    if shutil.which("gnuplot") is None or not plot_script.is_file():
        return
    output = POST_DIR / f"{study}_histories.png"
    completed = subprocess.run(
        [
            "gnuplot",
            "-e",
            f"histories='{POST_DIR / (study + '_histories.csv')}'; "
            f"ncurves={len(rows)}; "
            f"refdir='{REFERENCE_DIR}'; "
            f"output='{output}'",
            str(plot_script),
        ],
        cwd=VERIFY_DIR,
        text=True,
        check=False,
    )
    if completed.returncode:
        print("WARNING: gnuplot failed; the CSV results are still valid.")
    else:
        print(f"Plot: {output}")


def summary_table(rows: list[dict]) -> list[str]:
    header = (
        "| Member | Solid | nx | Solid ny | dt (s) | Peak (m) | Trough (m) | "
        "Mean (m) | Peak dev. | Trough dev. | Mean dev. | Band dist. | "
        "FSI its (mean/max) | Clock (s) |"
    )
    lines = [header, "|" + "---|" * (header.count("|") - 1)]
    for row in rows:
        lines.append(
            f"| {row['name']} | {row['solid']} | {row['nx']} | "
            f"{row['solid_ny']} | {row['delta_t_s']:g} | {row['peak_m']:.4f} | "
            f"{row['trough_m']:.4f} | {row['mean_m']:.4f} | "
            + " | ".join(
                f"{100.0 * row[f'{quantity}_envelope_deviation']:+.2f}%"
                for quantity in QUANTITIES
            )
            + f" | {100.0 * row['band_distance']:.2f}% | "
            f"{row['fsi_iterations_mean']:.1f}/{row['fsi_iterations_max']} | "
            f"{row['clock_time_s']:.0f} |"
        )
    return lines


# --------------------------------------------------------------------------- #
# Main
# --------------------------------------------------------------------------- #

def main() -> int:
    args = parse_args()
    config = json.loads(REFERENCE_FILE.read_text())

    studies = (
        list(STUDIES) if args.study == "all"
        else [value.strip() for value in args.study.split(",")]
    )
    invalid = [study for study in studies if study not in STUDIES]
    if invalid or not studies:
        raise SystemExit(f"--study must be a selection from: {', '.join(STUDIES)}, all")
    if args.cores < 1:
        raise SystemExit("--cores must be positive")

    end_time = (
        args.end_time if args.end_time is not None
        else float(config["quick_end_time_s"] if args.quick else config["end_time_s"])
    )
    window = config["comparison_window_s"]
    if end_time < float(window[0]) + float(config["lid_period_s"]) and not args.quick:
        raise SystemExit("--end-time must cover at least one lid period after "
                         f"t = {window[0]} s")

    required = ["blockMesh", "solids4Foam"]
    if args.cores > 1:
        required += ["decomposePar", "reconstructPar", "mpirun"]
    missing = [command for command in required if shutil.which(command) is None]
    if missing:
        raise SystemExit(f"required command(s) not found: {', '.join(missing)}")

    period = float(config["lid_period_s"])
    quick_window = [0.0, end_time] if args.quick else window
    config["comparison_window_s"] = quick_window
    references = {}
    for name, entry in config["references"].items():
        if "file" in entry:
            curve = read_curve(REFERENCE_DIR / entry["file"])
            # Some published curves stop before the end of the window
            stats_end = min(float(quick_window[1]), curve[-1][0])
            stats = periodic_stats(
                curve, float(quick_window[0]), stats_end, period
            )
        else:
            # A reference published as scalar values only
            curve = []
            stats = {key: float(value) for key, value in entry["values"].items()}
        references[name] = {
            "curve": curve,
            "group_b": bool(entry["group_b"]),
            "stats": stats,
        }
    reference_names = list(references)
    group_b = [name for name in reference_names if references[name]["group_b"]]
    curve_names = [name for name in reference_names if references[name]["curve"]]
    ranges = envelope(references)

    WORK_DIR.mkdir(parents=True, exist_ok=True)
    POST_DIR.mkdir(parents=True, exist_ok=True)

    lines = [
        "# mokLidDrivenCavity verification summary",
        "",
        f"- Mode: {'quick smoke test' if args.quick else 'full verification'}",
        f"- End time: {end_time:g} s; comparison window: "
        f"{quick_window[0]:g}-{min(quick_window[1], end_time):g} s",
        f"- Group B references: {', '.join(group_b)}; context only: "
        + ", ".join(name for name in reference_names if name not in group_b),
        "",
        "## References",
        "",
        "| Reference | Peak (m) | Trough (m) | Mean (m) |",
        "|---|---|---|---|",
    ]
    for name, reference in references.items():
        stats = reference["stats"]
        lines.append(
            f"| {name} | "
            + " | ".join(
                f"{stats[quantity]:.4f}" if quantity in stats else "-"
                for quantity in QUANTITIES
            )
            + " |"
        )
    lines += [
        "",
        "Group B envelope: "
        + ", ".join(
            f"{quantity} {ranges[quantity][0]:.4f}-{ranges[quantity][1]:.4f} m"
            for quantity in QUANTITIES
        ),
        "",
        "Deviations are relative to the nearer envelope bound (zero inside",
        "it). The band distance is the largest distance of the history from",
        "the band between the group B curves in the comparison window, as a",
        "fraction of the smaller group B curve peak.",
    ]

    all_checks: list[tuple[str, bool, str]] = []
    failures = 0
    for study in studies:
        rows = []
        for member in study_members(config, study, args.quick):
            try:
                row = run_member(study, member, end_time, args.cores, args.reuse)
                row["compare"] = bool(member.get("compare", True))
                evaluate(row, references, config)
            except RuntimeError as error:
                print(f"ERROR: {error}", file=sys.stderr)
                failures += 1
                if not args.keep_going:
                    return 1
                continue
            rows.append(row)
        if not rows:
            all_checks.append((f"{study}: no member completed", False, ""))
            lines += ["", f"## {study} study", "", f"- FAIL: {study}: no member completed"]
            continue
        write_study_csv(study, rows, curve_names)
        create_plot(study, rows)
        checks = check_study(study, rows, config, args.quick)
        all_checks += checks
        lines += ["", f"## {study} study", ""]
        lines += summary_table(rows)
        lines += [""]
        lines += [
            f"- {'PASS' if passed else 'FAIL'}: {description}"
            + (f" ({detail})" if detail else "")
            for description, passed, detail in checks
        ]

    passed = failures == 0 and all(passed for _, passed, _ in all_checks)
    lines += ["", f"- Result: {'PASS' if passed else 'FAIL'}"]
    summary = "\n".join(lines) + "\n"
    (POST_DIR / "verification_summary.md").write_text(summary)
    print(summary)
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
