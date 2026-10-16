#!/usr/bin/env python3
"""Run the opt-in womersleyTube verification study in isolated case copies.

The tube starts from the exact solution and is driven at both ends by it. The
driver takes the first Fourier coefficient, at the forcing frequency, of each
sampled history over the last period, and compares it with the complex
amplitude of the exact linear solution in
reference/womersleyTube_verification_references.json: the velocity profile at
x = L/2, the flow rate there, the wall displacement at x = L/4, L/2 and 3L/4,
and the wave number fitted to the pressure along the tube.
"""

from __future__ import annotations

import argparse
import cmath
import concurrent.futures
import csv
import hashlib
import json
import math
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

SCRIPT = Path(__file__).resolve()
VERIFICATION = SCRIPT.parents[1]
TUTORIAL = VERIFICATION.parent
REFERENCE_FILE = (
    VERIFICATION / "reference" / "womersleyTube_verification_references.json"
)
WORK_ROOT = VERIFICATION / "work"
OUTPUT_ROOT = VERIFICATION / "postProcessing"
SETTINGS_FILE = "verification_settings.json"
WALL = "postProcessing/0/solidPointDisplacement_wall{}.dat"
STATIONS = {"Quarter": 0.25, "Mid": 0.5, "ThreeQuarter": 0.75}
SAMPLES = "postProcessing/{}Samples/fluid"
RESIDUALS = "postProcessing/fsiResiduals.dat"
# Sampled fields per period: the Fourier coefficient at the forcing
# frequency is exact for a sinusoid, and the harmonics that could alias onto
# it (9 and 11 times the frequency) are negligible
SAMPLES_PER_PERIOD = 10
# Points of the axial pressure set in system/controlDict
AXIAL_POINTS = 61


def fail(message: str) -> None:
    raise RuntimeError(message)


# Errors from analysing damaged or degenerate results, reported as a failure
# of that case rather than aborting the study
ANALYSIS_ERRORS = (RuntimeError, ValueError, ZeroDivisionError)


def replace_entry(path: Path, key: str, value: str) -> None:
    """Replace the one top-level (unindented) entry key."""
    text = path.read_text()
    pattern = rf"^({re.escape(key)}\s+)[^;]+;"
    text, found = re.subn(pattern, rf"\g<1>{value};", text, flags=re.MULTILINE)
    if found != 1:
        fail(f"Expected one '{key}' entry in {path}, found {found}")
    path.write_text(text)


def is_time_directory(name: str) -> bool:
    try:
        float(name)
    except ValueError:
        return False
    return name != "0"


# ---------------------------------------------------------------------------
# Case preparation
# ---------------------------------------------------------------------------

def ignore_results(directory: str, names: list[str]) -> set[str]:
    """Skip run output and anything that is not part of the tutorial set-up,
    in particular time directories left by an earlier tutorial run."""
    ignored = set()
    for name in names:
        if (
            name in ("verification", "regressionTests", "postProcessing",
                     "dynamicCode", "case.foam")
            or name.startswith("processor")
            or name.startswith("log.")
            or name.endswith(".pdf")
            or is_time_directory(name)
        ):
            ignored.add(name)
        if name == "polyMesh" and Path(directory).name in ("fluid", "solid") \
                and Path(directory).parent.name == "constant":
            ignored.add(name)
    return ignored


def case_name(spec: dict) -> str:
    return "_".join([spec["coupling"], f"m{spec['factor']}",
                     f"n{spec['steps']}"])


def settings(spec: dict, reference: dict) -> dict:
    """Everything that defines a run, used to validate a reused case."""
    data = {k: v for k, v in spec.items() if k != "groups"}
    data["mesh"] = mesh_divisions(spec["factor"], reference)
    data["parameters"] = reference["reference"]["parameters"]
    data["tutorialFiles"] = tutorial_signature()
    data["build"] = build_signature()
    return data


def build_signature() -> dict:
    """OpenFOAM version and the solids4foam executable and library, so that a
    run is not reused with a different or rebuilt solver."""
    signature = {key: os.environ.get(key, "")
                 for key in ("WM_PROJECT", "WM_PROJECT_VERSION", "WM_OPTIONS")}
    candidates = [shutil.which("solids4Foam")]
    for variable in ("FOAM_USER_LIBBIN", "FOAM_MODULE_LIBBIN", "FOAM_SITE_LIBBIN"):
        libbin = os.environ.get(variable, "")
        if libbin:
            candidates += sorted(str(p) for p in Path(libbin).glob(
                "libsolids4FoamModels.*"))
    for candidate in candidates:
        if candidate:
            status = Path(candidate).stat()
            signature[candidate] = f"{status.st_size}:{status.st_mtime_ns}"
    return signature


_SIGNATURE: dict | None = None


def tutorial_signature() -> dict:
    """Content hashes of the tutorial files that define the runs."""
    global _SIGNATURE
    if _SIGNATURE is None:
        signature = {}
        for path in sorted(TUTORIAL.rglob("*")):
            relative = path.relative_to(TUTORIAL)
            if relative.parts[0] not in ("0", "constant", "system", "Allrun") \
                    or "polyMesh" in relative.parts or path.is_dir():
                continue
            if path.is_symlink():
                signature[str(relative)] = "->" + str(path.readlink())
            else:
                signature[str(relative)] = hashlib.sha1(
                    path.read_bytes()).hexdigest()
        _SIGNATURE = signature
    return _SIGNATURE


def mesh_divisions(factor: int, reference: dict) -> dict:
    base = reference["study"]["baseMesh"]
    return {key: value * factor for key, value in base.items()}


def configure(case: Path, spec: dict, reference: dict) -> None:
    shutil.copytree(TUTORIAL, case, symlinks=True, ignore=ignore_results)
    mesh = mesh_divisions(spec["factor"], reference)
    for region in ("fluid", "solid"):
        dict_path = case / "system" / region / "blockMeshDict"
        replace_entry(dict_path, "nAxial", str(mesh["axial"]))
        replace_entry(dict_path, "nRadial", str(mesh[f"{region}Radial"]))
    period = reference["reference"]["derived"]["period"]
    periods = reference["study"]["periods"]
    steps = spec["steps"]
    if steps % SAMPLES_PER_PERIOD:
        fail(f"{steps} steps per period is not a multiple of "
             f"{SAMPLES_PER_PERIOD}")
    control = case / "system/controlDict"
    replace_entry(control, "deltaT", f"{period / steps:.12g}")
    replace_entry(control, "endTime", f"{periods * period:.12g}")
    replace_entry(control, "writeInterval", str(periods * steps))
    replace_entry(control, "timePrecision", "12")
    # The fluid samples over the analysed periods
    window = reference["study"]["analysisPeriods"]
    text = control.read_text()
    text, found = re.subn(
        r"((?:profile|axial)Samples\s*\{[^}]*?timeStart\s+)[^;]+;([^}]*?writeInterval\s+)"
        r"\d+;",
        rf"\g<1>{(periods - window) * period:.12g};\g<2>"
        rf"{steps // SAMPLES_PER_PERIOD};", text, flags=re.DOTALL)
    if found != 2:
        fail("Could not set the fluid sampling times")
    control.write_text(text)


def solver_log_problem(case: Path) -> str | None:
    """Why the solver log does not show a completed run, or None. A
    completed log has an End line and no fatal error; PETSc and MPI can
    print warnings after End on shut-down."""
    log = case / "log.solids4Foam"
    if not log.is_file():
        return f"did not run the solver; see {case / 'log.Allverify'}"
    text = log.read_text(errors="replace")
    if re.search(r"FOAM FATAL|FOAM aborting", text):
        return f"failed; see {log}"
    if not re.search(r"^End\s*$", text, re.MULTILINE):
        return f"did not run to completion; see {log}"
    return None


def run_case(spec: dict, reference: dict, reuse: bool) -> str:
    case = WORK_ROOT / spec["name"]
    wanted = settings(spec, reference)
    stored = case / SETTINGS_FILE
    if reuse and case.is_dir():
        previous = None
        if stored.is_file():
            try:
                previous = json.loads(stored.read_text())
            except json.JSONDecodeError:
                previous = None
        if previous == wanted and solver_log_problem(case) is None:
            try:
                analyse(spec, reference)
                return f"{spec['name']}: reusing the completed run"
            except ANALYSIS_ERRORS as error:
                print(f"{spec['name']}: not reused ({error})")
        else:
            print(f"{spec['name']}: cached run does not match the requested "
                  "settings or did not finish; running it again")
    if case.exists():
        shutil.rmtree(case)
    configure(case, spec, reference)
    stored.write_text(json.dumps(wanted, indent=2, sort_keys=True))

    command = ["./Allrun", spec["coupling"]]
    log = case / "log.Allverify"
    with log.open("w") as handle:
        result = subprocess.run(command, cwd=case, stdout=handle,
                                stderr=subprocess.STDOUT, text=True)
    problem = solver_log_problem(case)
    if problem:
        fail(f"{spec['name']} {problem}")
    if result.returncode:
        # The tutorial Allrun only fails after the solver when plotting
        print(f"{spec['name']}: Allrun returned {result.returncode} after the "
              f"solver completed; see {log}")
    return f"{spec['name']}: completed"


# ---------------------------------------------------------------------------
# Exact solution
# ---------------------------------------------------------------------------

def besselI(n: int, z: complex) -> complex:
    term = (z / 2) ** n
    total = term
    for m in range(1, 500):
        term *= (z / 2) ** 2 / (m * (m + n))
        total += term
        if abs(term) <= 1e-17 * abs(total):
            return total
    fail(f"Bessel series did not converge for z = {z}")


def cplx(pair: list[float]) -> complex:
    return complex(pair[0], pair[1])


class Exact:
    """Complex amplitudes of the exact solution, fields = Re[a exp(i w t)]."""

    def __init__(self, reference: dict):
        ref = reference["reference"]
        prm = ref["parameters"]
        self.rho = prm["fluidDensity"]
        self.R = prm["innerRadius"]
        self.L = prm["length"]
        self.omega = ref["derived"]["omega"]
        self.period = ref["derived"]["period"]
        exact = ref["exact"]
        self.k = cplx(exact["k"])
        self.s = cplx(exact["s"])
        self.A = cplx(exact["A"])
        self.B = cplx(exact["B"])
        self.eta = cplx(exact["wallRadialDisplacement"])
        self.xi = cplx(exact["wallAxialDisplacement"])
        self.Q = cplx(exact["flowRate"])

    def wave(self, x: float) -> complex:
        return cmath.exp(-1j * self.k * x)

    def pressure(self, r: float, x: float) -> complex:
        return self.A * besselI(0, self.k * r) * self.wave(x)

    def velocity(self, r: float, x: float) -> complex:
        return (self.k * self.A / (self.rho * self.omega) * besselI(0, self.k * r)
                + self.B * besselI(0, self.s * r)) * self.wave(x)


# ---------------------------------------------------------------------------
# History analysis
# ---------------------------------------------------------------------------

def read_rows(path: Path, columns: int) -> list[list[float]]:
    if not path.is_file():
        fail(f"Missing history {path}")
    rows = []
    for line in path.read_text(errors="replace").splitlines():
        fields = line.split()
        if not fields or fields[0].startswith("#"):
            continue
        try:
            values = [float(v) for v in fields[:columns]]
        except ValueError:
            fail(f"Malformed line in {path}: {line!r}")
        if len(values) < columns:
            fail(f"Truncated line in {path}: {line!r}")
        if not all(math.isfinite(v) for v in values):
            fail(f"Non-finite value in {path}: {line!r}")
        rows.append(values)
    return rows


def fourier(times: list[float], values: list[complex | float],
            omega: float) -> complex:
    """First Fourier coefficient a of f = Re[a exp(i omega t)] from uniform
    samples covering one period."""
    n = len(times)
    return 2 / n * sum(v * cmath.exp(-1j * omega * t)
                       for t, v in zip(times, values))


def period_times(start: float, period: float, count: int) -> list[float]:
    return [start + period * i / count for i in range(count)]


def wall_history(case: Path, station: str, exact: Exact, steps: int,
                 periods: int) -> tuple[list[float], list[float], list[float]]:
    """Radial and axial inner-wall displacement at a station, checking the
    history is complete."""
    rows = read_rows(case / WALL.format(station), 4)
    times = [r[0] for r in rows]
    if any(b <= a for a, b in zip(times, times[1:])):
        fail(f"Non-increasing time in {case / WALL.format(station)}")
    n_steps = steps * periods
    end = periods * exact.period
    if len(rows) < n_steps or abs(times[-1] - end) > 0.25 * exact.period / steps:
        fail(f"Truncated history in {case}: {len(rows)} samples up to "
             f"t = {times[-1] if times else 'none'}, expected {n_steps} up to "
             f"t = {end:.6g}")
    # The point lies on the front wedge plane, at -0.5 degrees
    theta = -math.radians(0.5)
    radial = [r[2] * math.cos(theta) + r[3] * math.sin(theta) for r in rows]
    return times, radial, [r[1] for r in rows]


def last_period(times, values, period, steps, before=0, count=1):
    """The samples of the last `count` periods, t in
    (T_end - count period, T_end], or of the period `before` periods
    earlier."""
    end = times[-1] - before * period
    margin = 0.25 * period / steps
    selected = [(t, v) for t, v in zip(times, values)
                if end - count * period + margin < t < end + margin]
    if len(selected) != steps * count:
        fail(f"Found {len(selected)} samples in the last {count} periods, "
             f"expected {steps * count}")
    # The Fourier sum assumes uniform samples
    n = len(selected)
    for i, (t, _) in enumerate(selected):
        if abs(t - (end - (n - 1 - i) * period / steps)) > 1e-4 * period / steps:
            fail(f"Samples are not uniformly spaced in time near t = {t}")
    return [t for t, _ in selected], [v for _, v in selected]


def read_samples(case: Path, exact: Exact, periods: int, window: int,
                 expected: dict) -> dict:
    """Fourier coefficients of the profile and axial samples over the last
    `window` periods; `expected` holds the number of points of each set."""
    start = (periods - window) * exact.period
    wanted = period_times(start, window * exact.period,
                          window * SAMPLES_PER_PERIOD)
    data = {"profile": None, "axial": None}
    for set_name in data:
        root = case / SAMPLES.format(set_name)
        if not root.is_dir():
            fail(f"Missing fluid samples in {root}")
        available = {}
        for directory in root.iterdir():
            try:
                available[float(directory.name)] = directory
            except ValueError:
                continue
        series = []
        for t in wanted:
            match = [d for tt, d in available.items()
                     if abs(tt - t) < 1e-6 * exact.period]
            if len(match) != 1:
                fail(f"Missing fluid sample at t = {t:.6g} in {root}")
            rows = read_rows(match[0] / f"{set_name}_p_U.xy", 5)
            series.append(rows)
        n_points = len(series[0])
        if n_points != expected[set_name] \
                or any(len(s) != n_points for s in series):
            fail(f"Expected {expected[set_name]} {set_name} samples at every "
                 f"time in {root}, found {[len(s) for s in series]}")
        coord = [row[0] for row in series[0]]
        if any(b <= a for a, b in zip(coord, coord[1:])):
            fail(f"The {set_name} sample coordinates are not increasing in "
                 f"{root}")
        # The samples must be at the same places at every time, apart from
        # the motion of the mesh (the wall moves by 1% of the tutorial radial
        # cell size)
        spacing = min(b - a for a, b in zip(coord, coord[1:]))
        for s in series[1:]:
            if any(abs(row[0] - c) > 0.05 * spacing for row, c in zip(s, coord)):
                fail(f"The {set_name} sample coordinates change between "
                     f"times in {root}")
        p = [fourier(wanted, [s[i][1] for s in series], exact.omega)
             for i in range(n_points)]
        ux = [fourier(wanted, [s[i][2] for s in series], exact.omega)
              for i in range(n_points)]
        data[set_name] = {"coord": coord, "p": p, "ux": ux}
    return data


def flow_rate_history(case: Path, exact: Exact, steps: int,
                      periods: int) -> tuple[list[float], list[float]]:
    """Flow rate through the faceZone at x = L/2 (surfaceFieldValue)."""
    found = sorted(case.glob("postProcessing/**/flowRate/0/surfaceFieldValue.dat"))
    if len(found) != 1:
        fail(f"Expected one flow-rate history in {case}, found {len(found)}")
    rows = read_rows(found[0], 2)
    times = [r[0] for r in rows]
    if any(b <= a for a, b in zip(times, times[1:])):
        fail(f"Non-increasing time in {found[0]}")
    end = periods * exact.period
    if len(rows) < steps * periods \
            or abs(times[-1] - end) > 0.25 * exact.period / steps:
        fail(f"Truncated flow-rate history in {case}")
    return times, [r[1] for r in rows]


def fit_wave_number(x: list[float], p: list[complex]) -> complex:
    """Least-squares fit of log p = log C - i k x, with the phase unwrapped."""
    logs, previous = [], None
    for value in p:
        if value == 0:
            fail("Zero pressure harmonic in the axial samples")
        z = cmath.log(value)
        if previous is not None:
            while z.imag - previous.imag > math.pi:
                z -= 2j * math.pi
            while z.imag - previous.imag < -math.pi:
                z += 2j * math.pi
        logs.append(z)
        previous = z
    n = len(x)
    mean_x, mean_y = sum(x) / n, sum(logs) / n
    slope = sum((a - mean_x) * (b - mean_y) for a, b in zip(x, logs)) \
        / sum((a - mean_x) ** 2 for a in x)
    return 1j * slope


def phase_error(value: complex, exact: complex) -> float:
    if value == 0 or exact == 0:
        fail("Zero harmonic amplitude")
    return cmath.phase(value / exact)


def coupling_iterations(case: Path) -> tuple[float, int]:
    path = case / RESIDUALS
    if not path.is_file():
        return math.nan, 0
    per_step: dict[str, int] = {}
    for line in path.read_text(errors="replace").splitlines()[1:]:
        fields = line.split()
        if len(fields) >= 2:
            per_step[fields[0]] = max(per_step.get(fields[0], 0), int(fields[1]))
    if not per_step:
        return math.nan, 0
    counts = list(per_step.values())
    return sum(counts) / len(counts), max(counts)


def execution_time(case: Path) -> float:
    matches = re.findall(r"^ExecutionTime = ([0-9.eE+-]+) s",
                         (case / "log.solids4Foam").read_text(errors="replace"),
                         re.MULTILINE)
    return float(matches[-1]) if matches else math.nan


def analyse(spec: dict, reference: dict) -> dict:
    case = WORK_ROOT / spec["name"]
    exact = Exact(reference)
    periods = reference["study"]["periods"]
    window = reference["study"]["analysisPeriods"]
    steps = spec["steps"]
    row = dict(spec)

    # Wall displacement at the three stations
    for station, fraction in STATIONS.items():
        times, radial, axial = wall_history(case, station, exact, steps, periods)
        t, v = last_period(times, radial, exact.period, steps, count=window)
        eta = fourier(t, v, exact.omega)
        eta_exact = exact.eta * exact.wave(fraction * exact.L)
        row[f"wall{station}_amp"] = abs(eta) / abs(eta_exact) - 1
        row[f"wall{station}_phase"] = phase_error(eta, eta_exact)
        t, v = last_period(times, axial, exact.period, steps, count=window)
        xi = fourier(t, v, exact.omega)
        xi_exact = exact.xi * exact.wave(fraction * exact.L)
        row[f"axial{station}_amp"] = abs(xi) / abs(xi_exact) - 1
        row[f"axial{station}_phase"] = phase_error(xi, xi_exact)
        if station == "Mid":
            # Periodicity: the change from the last but one period to the
            # last, left by the start-up transient
            t0, v0 = last_period(times, radial, exact.period, steps, before=1)
            t1, v1 = last_period(times, radial, exact.period, steps)
            last = fourier(t1, v1, exact.omega)
            if last == 0:
                fail("Zero wall displacement harmonic")
            row["periodicity"] = abs(fourier(t0, v0, exact.omega) / last - 1)

    mesh = mesh_divisions(spec["factor"], reference)
    samples = read_samples(case, exact, periods, window,
                           {"profile": mesh["fluidRadial"],
                            "axial": AXIAL_POINTS})
    # Velocity profile and flow rate at x = L/2
    prof = samples["profile"]
    x_mid = 0.5 * exact.L
    # The profile holds the values of the cells just downstream of x = L/2,
    # centred half a cell further on, and at the centroid radius of each
    # wedge cell, (2/3)(r2^3 - r1^3)/(r2^2 - r1^2), rather than at the
    # midpoint between its faces, (r1 + r2)/2, that the set reports; both are
    # on the wedge mid-plane, where the faces are at r cos(0.5 degrees)
    x_cells = x_mid + 0.5 * exact.L / mesh["axial"]
    dr = exact.R / mesh["fluidRadial"]
    chord = math.cos(math.radians(0.5))
    radii = []
    for j, r in enumerate(prof["coord"]):
        if abs(r - (j + 0.5) * dr) > 0.1 * dr:
            fail(f"Unexpected profile sample radius {r} for cell {j}")
        r1, r2 = j * dr, (j + 1) * dr
        radii.append(chord * 2 / 3 * (r2**3 - r1**3) / (r2**2 - r1**2))
    u_exact = [exact.velocity(r, x_cells) for r in radii]
    scale = max(abs(u) for u in u_exact)
    row["profile"] = max(abs(u - ue) for u, ue in zip(prof["ux"], u_exact)) / scale
    # The profile itself, normalised by the exact peak, to compare runs;
    # entries starting with an underscore are not written to the CSV file
    row["_profile"] = [u / scale for u in prof["ux"]]
    for label, phi in (("0", 0.0), ("90", 0.5 * math.pi), ("180", math.pi),
                       ("270", 1.5 * math.pi)):
        row[f"profile{label}"] = max(
            abs(((u - ue) * cmath.exp(1j * phi)).real)
            for u, ue in zip(prof["ux"], u_exact)) / scale
    times, flux = flow_rate_history(case, exact, steps, periods)
    t, v = last_period(times, flux, exact.period, steps, count=window)
    # The flow rate is through the 1 degree wedge
    q = 360 * fourier(t, v, exact.omega)
    q_exact = exact.Q * exact.wave(x_mid)
    row["flow_amp"] = abs(q) / abs(q_exact) - 1
    row["flow_phase"] = phase_error(q, q_exact)

    # Pressure along the tube at r = R/2: wave number
    ax = samples["axial"]
    k = fit_wave_number(ax["coord"], ax["p"])
    row["k_real"] = k.real
    row["k_imag"] = k.imag
    row["speed"] = (exact.omega / k.real) / (exact.omega / exact.k.real) - 1
    row["attenuation"] = k.imag / exact.k.imag - 1
    # The pimpleFluid pressure is kinematic
    p_exact = [exact.pressure(0.5 * exact.R, x) for x in ax["coord"]]
    row["pressure"] = max(abs(exact.rho * p - pe)
                          for p, pe in zip(ax["p"], p_exact)) \
        / max(abs(pe) for pe in p_exact)

    row["mean_iterations"], row["max_iterations"] = coupling_iterations(case)
    row["execution_time"] = execution_time(case)
    return row


# ---------------------------------------------------------------------------
# Study definition
# ---------------------------------------------------------------------------

def build_specs(args: argparse.Namespace, reference: dict) -> list[dict]:
    study = reference["study"]
    if args.quick:
        studies = {"quick"}
    else:
        studies = {"main", "coupling"} if args.study == "all" else {args.study}

    specs: dict[str, dict] = {}

    def add(factor: int, steps: int, coupling: str, group: str) -> None:
        spec = {"coupling": coupling, "factor": factor, "steps": steps}
        spec["name"] = case_name(spec)
        entry = specs.setdefault(spec["name"], spec)
        entry.setdefault("groups", set()).add(group)

    if "quick" in studies:
        quick = study["quick"]
        add(quick["meshFactor"], quick["steps"], study["studyCoupling"],
            "quick")
    if "main" in studies:
        for factor in study["meshFactors"]:
            add(factor, study["meshStudySteps"], study["studyCoupling"], "mesh")
        for steps in study["timeSteps"]:
            add(study["timeStudyFactor"], steps, study["studyCoupling"], "time")
    if "coupling" in studies:
        for coupling in ("iqnils", "robin"):
            add(study["couplingFactor"], study["couplingSteps"], coupling,
                "coupling")
    for spec in specs.values():
        spec["groups"] = sorted(spec["groups"])
    return list(specs.values())


# ---------------------------------------------------------------------------
# Checks and reporting
# ---------------------------------------------------------------------------

METRICS = [
    ("profile", "velocity profile at L/2, max |error| / max |u|"),
    ("flow_amp", "flow rate amplitude at L/2"),
    ("flow_phase", "flow rate phase at L/2 (rad)"),
    ("wallMid_amp", "wall radial displacement amplitude at L/2"),
    ("wallMid_phase", "wall radial displacement phase at L/2 (rad)"),
    ("speed", "wave speed from the pressure along the tube"),
    ("attenuation", "attenuation rate, Im(k)"),
]

# Signed errors, whose differences between runs are differences of the
# computed values
SIGNED = ("flow_amp", "flow_phase", "wallMid_amp", "wallMid_phase", "speed",
          "attenuation")


def error_order(errors: list[float]) -> float:
    """Order from the last two errors at a refinement ratio of 2."""
    a, b = abs(errors[-2]), abs(errors[-1])
    if a == 0 or b == 0:
        return math.nan
    return math.log(a / b) / math.log(2)


def difference_order(values: list[float]) -> float:
    """Order from three values at a refinement ratio of 2."""
    d1, d2 = values[0] - values[1], values[1] - values[2]
    if d1 == 0 or d2 == 0 or d1 * d2 < 0:
        return math.nan
    return math.log(abs(d1 / d2)) / math.log(2)


def find(rows: list[dict], **keys) -> dict | None:
    for row in rows:
        if all(row.get(k) == v for k, v in keys.items()):
            return row
    return None


def fmt(value: float, fmt_spec: str = ".2e") -> str:
    return "–" if value is None or not math.isfinite(value) \
        else format(value, fmt_spec)


def check_rows(rows: list[dict], reference: dict, args: argparse.Namespace,
               lines: list[str]) -> list[str]:
    tol = reference["tolerances"]
    study = reference["study"]
    failures: list[str] = []
    orders = ["", "## Observed orders", "",
              "Unsigned errors (profile): from the errors of the two finest "
              "runs. Signed errors: from the differences between the three "
              "runs, which cancel the error shared by the series, such as "
              "the time error of the mesh study. Time-step series of more "
              "than three runs give one order per three successive steps.",
              "",
              "| Series | Quantity | Errors | Order |", "|---|---|---|---:|"]

    def check(condition: bool, message: str) -> None:
        lines.append(f"- {'PASS' if condition else 'FAIL'}: {message}")
        if not condition:
            failures.append(message)

    def check_errors(row: dict, limits: dict, label: str) -> None:
        for key, description in METRICS:
            if key in limits:
                check(abs(row[key]) <= limits[key],
                      f"{label}: {description} {fmt(row[key])} within "
                      f"{limits[key]:g}")

    for row in rows:
        check(row["periodicity"] <= tol["periodicity"],
              f"{row['name']}: last two periods differ by "
              f"{fmt(row['periodicity'])} (within {tol['periodicity']:g})")

    if args.quick:
        for row in rows:
            check_errors(row, tol["quick"], row["name"])
        return failures

    def series(group: str, key: str, values: list, **keys):
        found = [find(rows, **{key: v}, **keys) for v in values]
        found = [f for f in found if f is not None and group in f["groups"]]
        return found if len(found) == len(values) else None

    coupling = study["studyCoupling"]
    mesh = series("mesh", "factor", study["meshFactors"], coupling=coupling,
                  steps=study["meshStudySteps"])
    if mesh:
        check_errors(mesh[-1], tol["finest"], f"finest mesh ({mesh[-1]['name']})")
        for key, _ in METRICS:
            errors = [r[key] for r in mesh]
            # The signed errors share the time error of the common step, which
            # their differences cancel
            order = difference_order(errors) if key in SIGNED \
                else error_order(errors)
            orders.append(f"| mesh | {key} | "
                          f"{', '.join(fmt(e) for e in errors)} | "
                          f"{fmt(order, '.2f')} |")
            if key in tol["meshOrderQuantities"]:
                check(math.isfinite(order) and order >= tol["minMeshOrder"],
                      f"mesh: observed order {fmt(order, '.2f')} of the {key} "
                      f"error at least {tol['minMeshOrder']}")
    time = series("time", "steps", study["timeSteps"], coupling=coupling,
                  factor=study["timeStudyFactor"])
    if time:
        for key in SIGNED:
            values = [r[key] for r in time]
            # One order per three successive step sizes; the check is on the
            # finest three, where an O(dt) term would show first
            triplet_orders = [difference_order(values[i:i + 3])
                              for i in range(len(values) - 2)]
            order = triplet_orders[-1]
            orders.append(f"| time step | {key} | "
                          f"{', '.join(fmt(e) for e in values)} | "
                          f"{', '.join(fmt(o, '.2f') for o in triplet_orders)} |")
            if key in tol["timeOrderQuantities"]:
                check(math.isfinite(order) and order >= tol["minTimeOrder"],
                      f"time step: observed order {fmt(order, '.2f')} of the "
                      f"{key} (finest three steps) at least "
                      f"{tol['minTimeOrder']}")
        # The time error of the mesh-study step size, estimated from the
        # difference with the next smaller step
        finest = time[-1]
        coarser = time[-2]
        lines.append(f"  - time error at {coarser['steps']} steps per period: "
                     "wall amplitude "
                     f"{fmt(coarser['wallMid_amp'] - finest['wallMid_amp'])}, "
                     "wave speed "
                     f"{fmt(coarser['speed'] - finest['speed'])}")
    iq = find(rows, coupling="iqnils", factor=study["couplingFactor"],
              steps=study["couplingSteps"])
    rb = find(rows, coupling="robin", factor=study["couplingFactor"],
              steps=study["couplingSteps"])
    if iq is not None and rb is not None \
            and "coupling" in iq["groups"] and "coupling" in rb["groups"]:
        check_errors(iq, tol["coupling"], f"IQN-ILS ({iq['name']})")
        check_errors(rb, tol["coupling"], f"Robin ({rb['name']})")
        for key, description in METRICS:
            if key == "profile":
                # The profiles themselves, not their error norms
                description = "velocity profile at L/2, max |difference| / " \
                    "max |u|"
                difference = max(abs(a - b) for a, b in
                                 zip(rb["_profile"], iq["_profile"]))
            else:
                difference = abs(rb[key] - iq[key])
            check(difference <= tol["couplingAgreement"],
                  f"Robin and IQN-ILS agree on the {description} to "
                  f"{fmt(difference)} (within {tol['couplingAgreement']:g})")
        lines.append(f"  - mean coupling iterations per step: "
                     f"{fmt(iq['mean_iterations'], '.2f')} (IQN-ILS), "
                     f"{fmt(rb['mean_iterations'], '.2f')} (Robin); run time "
                     f"{fmt(iq['execution_time'], '.0f')} s and "
                     f"{fmt(rb['execution_time'], '.0f')} s")
    if len(orders) > 7:
        lines += orders
    return failures


def table(rows: list[dict]) -> list[str]:
    head = ["Case"] + [key for key, _ in METRICS] + [
        "wallQuarter_amp", "wallThreeQuarter_amp", "axialMid_amp",
        "periodicity", "Iter./step", "Time (s)"]
    out = ["| " + " | ".join(head) + " |",
           "|---|" + "---:|" * (len(head) - 1)]
    for r in rows:
        values = [fmt(r[key]) for key, _ in METRICS] + [
            fmt(r["wallQuarter_amp"]), fmt(r["wallThreeQuarter_amp"]),
            fmt(r["axialMid_amp"]), fmt(r["periodicity"]),
            fmt(r["mean_iterations"], ".2f"), fmt(r["execution_time"], ".0f")]
        out.append(f"| {r['name']} | " + " | ".join(values) + " |")
    return out


def write_csv(path: Path, rows: list[dict]) -> None:
    fields = ["name", "coupling", "factor", "steps"] + sorted(
        {k for r in rows for k in r if not k.startswith("_")}
        - {"name", "coupling", "factor", "steps", "groups"})
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--study", choices=("all", "main", "coupling"),
                        default="all",
                        help="main: mesh and time-step studies; coupling: "
                        "IQN-ILS against Robin-Neumann (default: all)")
    parser.add_argument("--cores", type=int, default=1,
                        help="number of serial cases run at the same time")
    parser.add_argument("--reuse", action="store_true",
                        help="reuse completed runs under verification/work whose "
                        "settings match")
    parser.add_argument("--quick", action="store_true",
                        help="short smoke run on the coarsest mesh")
    args = parser.parse_args()
    if args.cores < 1:
        parser.error("--cores must be a positive integer")
    for executable in ("blockMesh", "solids4Foam"):
        if shutil.which(executable) is None:
            fail(f"Required executable '{executable}' is unavailable. "
                 "Source OpenFOAM and build solids4foam first.")
    reference = json.loads(REFERENCE_FILE.read_text())
    specs = build_specs(args, reference)
    WORK_ROOT.mkdir(parents=True, exist_ok=True)
    OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)

    print(f"Running {len(specs)} cases, {args.cores} at a time")
    errors = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.cores) as pool:
        futures = {pool.submit(run_case, spec, reference, args.reuse): spec
                   for spec in specs}
        for future in concurrent.futures.as_completed(futures):
            try:
                print(future.result(), flush=True)
            except RuntimeError as error:
                errors.append(str(error))
                print(f"ERROR: {error}", flush=True)

    rows = []
    for spec in specs:
        try:
            rows.append(analyse(spec, reference))
        except ANALYSIS_ERRORS as error:
            if not any(spec["name"] in e for e in errors):
                errors.append(f"{spec['name']}: {error}")
    rows.sort(key=lambda r: (r["coupling"], r["factor"], r["steps"]))

    lines = ["# womersleyTube verification summary", "",
             "Errors are relative to the exact linear solution unless stated.",
             ""]
    lines += table(rows)
    lines += ["", "## Checks", ""]
    failures = check_rows(rows, reference, args, lines)
    for error in errors:
        lines.append(f"- FAIL: {error}")
    failures += errors
    lines += ["", "PASSED" if not failures else f"FAILED ({len(failures)} checks)"]
    summary = OUTPUT_ROOT / "verification_summary.md"
    summary.write_text("\n".join(lines) + "\n")
    write_csv(OUTPUT_ROOT / "womersleyTube_results.csv", rows)
    print("\n".join(lines))
    print(f"\nWrote {summary}")
    return 0 if not failures else 1


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (RuntimeError, subprocess.SubprocessError, OSError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        sys.exit(2)
