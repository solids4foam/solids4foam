#!/usr/bin/env python3
"""Run the opt-in Hessenthaler et al. (2017) FSI validation studies.

calibration
    The flap alone under its net buoyancy (no fluid), relaxed with damping to
    its static deflection at three buoyancy load factors. The static
    deflection of a neo-Hookean solid with a fixed Poisson's ratio depends
    only on rho*g/mu, so the load factors map to shear moduli, and
    interpolation to the measured zero-flow tip deflection of 29.50 mm gives
    the calibrated shear modulus. The study covers three solid meshes, two
    Poisson's ratios, the stabilisation, the high-order solid, and one
    parallel run.

phaseI
    The coupled steady-inflow Phase I case on the coarse and medium fluid
    meshes, run to a steady state and compared with the measured centreline
    of the flap and with the phase-contrast MRI velocity, averaged over each
    MRI voxel, on the planes z = 10 and 30 mm.

Every run is a copy of the tutorial under validation/work; results are
written to validation/postProcessing.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import re
import shutil
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

SCRIPT = Path(__file__).resolve()
VALIDATION = SCRIPT.parents[1]
TUTORIAL = VALIDATION.parent
REFERENCE_DIR = VALIDATION / "reference"
REFERENCE_FILE = REFERENCE_DIR / "hessenthalerFsi_validation_references.json"
WORK_ROOT = VALIDATION / "work"
OUTPUT_ROOT = VALIDATION / "postProcessing"
SETTINGS_FILE = "validation_settings.json"

# Hexahedral solid meshes (cells in x, y, z) of the calibration study
SOLID_MESHES = {
    "coarse": (6, 4, 33),
    "medium": (12, 8, 65),
    "fine": (18, 12, 98),
}

# Voxel sub-sampling of the MRI comparison (points per voxel in x, y, z)
VOXEL_SAMPLES = (4, 4, 12)


def fail(message: str) -> None:
    raise RuntimeError(message)


# ---------------------------------------------------------------------------
# File helpers
# ---------------------------------------------------------------------------

def read_csv(path: Path) -> list[dict[str, float]]:
    """Read a CSV file with '#' comment lines into a list of float rows."""
    lines = [line for line in path.read_text().splitlines()
             if line.strip() and not line.startswith("#")]
    rows = []
    for row in csv.DictReader(lines):
        rows.append({key: float(value) for key, value in row.items()})
    return rows


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        return
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def replace_entry(path: Path, key: str, value: str,
                  top_level: bool = False, count: int = 1) -> None:
    """Set dictionary entries without GNU sed; top_level skips indented ones."""
    text = path.read_text()
    indent = "" if top_level else r"[ \t]*"
    pattern = rf"^({indent}{re.escape(key)}\s+)[^;\n]+;"
    text, found = re.subn(pattern, rf"\g<1>{value};", text, flags=re.MULTILINE)
    if found != count:
        fail(f"Expected {count} '{key}' entries in {path}, found {found}")
    path.write_text(text)


def replace_text(path: Path, old: str, new: str) -> None:
    text = path.read_text()
    if text.count(old) != 1:
        fail(f"Expected one '{old.strip()}' in {path}")
    path.write_text(text.replace(old, new))


def numeric_rows(path: Path) -> list[list[float]]:
    """All data rows of a whitespace-separated file.

    Comment lines and one leading header line are skipped. Any other row
    that is not entirely finite numbers, such as one cut short by an aborted
    run, is an error rather than being dropped.
    """
    rows = []
    header_allowed = True
    for number, line in enumerate(
        path.read_text(errors="replace").splitlines(), 1
    ):
        fields = line.replace("(", " ").replace(")", " ").split()
        if not fields or fields[0].startswith("#"):
            continue
        try:
            values = [float(field) for field in fields]
        except ValueError:
            if header_allowed:
                header_allowed = False
                continue
            fail(f"{path}:{number} is not a numeric data row: {line.strip()}")
        header_allowed = False
        if not all(math.isfinite(value) for value in values):
            fail(f"{path}:{number} contains non-finite values")
        rows.append(values)
    return rows


def last_per_time(rows: list[list[float]]) -> list[list[float]]:
    """Keep the final row written for each time, in time order."""
    by_time: dict[float, list[float]] = {}
    for row in rows:
        by_time[row[0]] = row
    return [by_time[time] for time in sorted(by_time)]


def require_history(path: Path, rows: list[list[float]], end_time: float,
                    delta_t: float, columns: int) -> None:
    """Require complete, finite rows every time step up to the end time."""
    if not rows:
        fail(f"{path} has no data rows")
    short = [row for row in rows if len(row) < columns]
    if short:
        fail(f"{path} has {len(short)} row(s) with fewer than {columns} "
             f"columns, starting at t = {short[0][0]:g}")
    times = [row[0] for row in rows]
    if not math.isclose(times[-1], end_time, rel_tol=1e-6, abs_tol=1e-9):
        fail(f"{path} ends at t = {times[-1]:g}, not at the end time "
             f"t = {end_time:g}")
    if times[0] > 1.5 * delta_t:
        fail(f"{path} starts at t = {times[0]:g}, not at the first time "
             f"step t = {delta_t:g}")
    gaps = [b - a for a, b in zip(times, times[1:])]
    if gaps and max(gaps) > 1.5 * delta_t:
        fail(f"{path} has a gap of {max(gaps):g} s in its history")
    expected = round(end_time / delta_t)
    if len(times) < expected:
        fail(f"{path} has {len(times)} time levels, expected {expected}")


# ---------------------------------------------------------------------------
# Fingerprints and reuse
# ---------------------------------------------------------------------------

def tutorial_fingerprint() -> str:
    """Hash of the tutorial inputs that a run copies."""
    digest = hashlib.sha256()
    for directory in ("0", "constant", "system", "geometry"):
        for path in sorted((TUTORIAL / directory).rglob("*")):
            relative = path.relative_to(TUTORIAL)
            if "polyMesh" in path.parts or "triSurface" in path.parts:
                continue
            if path.is_symlink():
                digest.update(f"{relative}->{os.readlink(path)}".encode())
            elif path.is_file():
                digest.update(str(relative).encode())
                digest.update(path.read_bytes())
    for name in ("Allrun", "makeFluidSurface.py"):
        digest.update((TUTORIAL / name).read_bytes())
    return digest.hexdigest()


def driver_fingerprint() -> str:
    """Hash of this script, which generates or edits the run inputs."""
    return hashlib.sha256(SCRIPT.read_bytes()).hexdigest()


def build_fingerprint() -> str:
    """Identify the solids4foam build: solver and library size and time."""
    parts = [os.environ.get("WM_PROJECT", ""),
             os.environ.get("WM_PROJECT_VERSION", ""),
             os.environ.get("PETSC_DIR", "")]
    solver = shutil.which("solids4Foam")
    candidates = [Path(solver)] if solver else []
    for variable in ("FOAM_USER_LIBBIN", "FOAM_SITE_LIBBIN", "FOAM_LIBBIN"):
        directory = os.environ.get(variable)
        if directory:
            candidates.extend(Path(directory).glob("libsolids4FoamModels.*"))
    for path in candidates:
        if path.is_file():
            stat = path.stat()
            parts.append(f"{path}:{stat.st_size}:{int(stat.st_mtime)}")
    return hashlib.sha256("\n".join(parts).encode()).hexdigest()


def solver_completed(case: Path) -> bool:
    log = case / "log.solids4Foam"
    return log.is_file() and bool(
        re.search(r"^End\s*$", log.read_text(errors="replace"), re.MULTILINE)
    )


def check_solver_log(case: Path, label: str) -> None:
    """Reject a run that aborted, even if its script returned zero."""
    log = case / "log.solids4Foam"
    if not log.is_file():
        fail(f"{label} did not create {log}")
    text = log.read_text(errors="replace")
    if re.search(r"FOAM FATAL|FOAM aborting|PETSC ERROR|PRTE ERROR|"
                 r"MPI_ABORT|^ERROR$", text, re.MULTILINE):
        fail(f"{label} failed; see {log}")
    if not re.search(r"^End\s*$", text, re.MULTILINE):
        fail(f"{label} did not run to completion; see {log}")


def reusable(case: Path, settings: dict, args: argparse.Namespace,
             label: str) -> bool:
    if not args.reuse or not solver_completed(case):
        return False
    settings_file = case / SETTINGS_FILE
    stored = (json.loads(settings_file.read_text())
              if settings_file.is_file() else None)
    if stored != settings:
        print(f"Not reusing {case}: it was run with other settings or "
              "another build")
        return False
    check_solver_log(case, label)
    print(f"Reusing {label} in {case}")
    return True


def run_environment() -> dict[str, str]:
    """The environment of the OpenFOAM commands.

    macOS drops DYLD_LIBRARY_PATH when a protected interpreter such as
    /usr/bin/env starts this script; restore it as OpenFOAM does.
    """
    env = dict(os.environ)
    if sys.platform == "darwin" and not env.get("DYLD_LIBRARY_PATH") \
            and env.get("FOAM_LD_LIBRARY_PATH"):
        env["DYLD_LIBRARY_PATH"] = env["FOAM_LD_LIBRARY_PATH"]
    return env


def run(command: list[str], case: Path, log_name: str) -> None:
    log = case / log_name
    with log.open("w") as handle:
        result = subprocess.run(command, cwd=case, stdout=handle,
                                stderr=subprocess.STDOUT, text=True,
                                env=run_environment())
    if result.returncode:
        fail(f"'{' '.join(command)}' failed in {case}; see {log}")


def clock_time(case: Path) -> float:
    matches = re.findall(
        r"ClockTime\s*=\s*([0-9.eE+-]+)",
        (case / "log.solids4Foam").read_text(errors="replace"),
    )
    return float(matches[-1]) if matches else math.nan


def set_subdomains(case: Path, cores: int, regions: tuple[str, ...]) -> None:
    for region in regions:
        path = case / "system" / region / "decomposeParDict"
        replace_entry(path, "numberOfSubdomains", str(cores), True)


def use_parallel_lu(fv_solution: Path) -> None:
    """Replace the block-Jacobi LU by an exact parallel LU (MUMPS).

    Block Jacobi drops the coupling between subdomains, so its Krylov
    iterations depend on the decomposition; it stalls for the high-order
    solid of this thin flap, although the tutorial's updated Lagrangian solid
    also converged with it on 24 ranks. The exact parallel LU makes the solid
    solve independent of the decomposition.
    """
    text = fv_solution.read_text()
    text, found = re.subn(
        r"pc_type\s+bjacobi\s*;\s*sub_pc_type\s+lu\s*;",
        "pc_type lu;\n            pc_factor_mat_solver_type mumps;", text)
    if found != 1:
        fail(f"Expected one block-Jacobi LU preconditioner in {fv_solution}")
    fv_solution.write_text(text)


# ---------------------------------------------------------------------------
# Calibration study: the flap alone under buoyancy
# ---------------------------------------------------------------------------

# Buoyancy load factors, relative to the Phase I value, of the plateaus of a
# calibration run. The static deflection of a neo-Hookean solid with a fixed
# Poisson's ratio depends only on rho*g/mu, so load factor f at shear modulus
# MU_REFERENCE is equivalent to shear modulus MU_REFERENCE/f at load 1.
LOAD_LEVELS = (0.9, 1.05, 1.2)
# Each plateau: a smooth 0.25 s transition, then damped settling
PLATEAU_TIME = 2.0
TRANSITION_TIME = 0.25
CAL_DELTA_T = 0.005
CAL_DAMPING = 15.0
MU_REFERENCE = 61000.0

LOAD_STEPS = """
loadSteps
{
    type            coded;
    libs            (utilityFunctionObjects);
    name            loadSteps;

    codeInclude
    #{
        #include "uniformDimensionedFields.H"
    #};

    codeData
    #{
        vector gFinal_{Zero};
    #};

    localCode
    #{
        // Piecewise-constant buoyancy load factor with smooth transitions
        static void setLoad(const fvMesh& mesh, const vector& gFinal)
        {
            const scalarList levels({LEVELS});
            const scalar plateau = PLATEAU;
            const scalar transition = TRANSITION;
            const scalar t = mesh.time().value() + mesh.time().deltaTValue();
            const label k = min(label(t/plateau), levels.size() - 1);
            const scalar previous = k > 0 ? levels[k - 1] : 0.0;
            const scalar s = min(max((t - k*plateau)/transition, 0.0), 1.0);
            const scalar f = previous + (levels[k] - previous)*s*s*(3 - 2*s);
            uniformDimensionedVectorField& g =
                const_cast<uniformDimensionedVectorField&>
                (
                    mesh.lookupObject<uniformDimensionedVectorField>("g")
                );
            g.value() = f*gFinal;
        }
    #};

    codeRead
    #{
        // The unscaled value is read from constant/g on disk, not from the
        // registered field, which holds the scaled value
        const uniformDimensionedVectorField gFile
        (
            IOobject
            (
                "g",
                mesh().time().constant(),
                mesh(),
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false
            )
        );
        gFinal_ = gFile.value();
        setLoad(mesh(), gFinal_);
    #};

    codeExecute
    #{
        setLoad(mesh(), gFinal_);
    #};
}
"""

# High-order total Lagrangian solid of the comparison: cubic residual; the
# Jacobian (compact or the linear high-order one) and the extra stencil cells
# are set per run
HIGH_ORDER_PROPERTIES = """
solidModel     nonLinearGeometryTotalLagrangianTotalDisplacement;

nonLinearGeometryTotalLagrangianTotalDisplacementCoeffs
{
    solutionAlgorithm PETScSNES;

    predictor yes;

    highOrderCoeffs
    {
        highOrderResidual true;
        highOrderJacobian true;
        displacement
        {
            polynomialOrder 3;
            faceStencilExtraCells 65;
            weightFunctionCoeffs
            {
                type radiallySymmetricExponential;
                k 6;
            }
            calcConditionNumber false;
        }
    }

    stabilisation
    {
        momentum
        {
            type        alpha;
            scaleFactor 0.1;
        }
    }
}
"""


# Momentum stabilisation scale factor of the tutorial's solid
TUTORIAL_STABILISATION = 0.01

# Standard solid models of the calibration study; "tutorial" keeps the
# tutorial's own
SOLID_MODELS = {
    "tutorial": None,
    "ul": "nonLinearGeometryUpdatedLagrangian",
    "tl": "nonLinearGeometryTotalLagrangianTotalDisplacement",
}


def set_displacement_field(directory: Path, name: str) -> None:
    """Name the initial displacement field D (total) or DD (increment)."""
    other = "DD" if name == "D" else "D"
    if (directory / other).is_file():
        text = (directory / other).read_text().replace(
            f"object      {other};", f"object      {name};")
        (directory / other).unlink()
        (directory / name).write_text(text)


def solid_case(name: str, spec: dict, end_time: float) -> Path:
    """Build a solid-only copy of the tutorial's solid region."""
    case = WORK_ROOT / name
    if case.exists():
        shutil.rmtree(case)
    (case / "system").mkdir(parents=True)
    shutil.copytree(TUTORIAL / "0" / "solid", case / "0")
    shutil.copytree(TUTORIAL / "constant" / "solid", case / "constant")
    for path in (TUTORIAL / "system" / "solid").iterdir():
        shutil.copy(path, case / "system" / path.name)
    (case / "constant" / "physicsProperties").write_text(
        (TUTORIAL / "constant" / "physicsProperties").read_text()
        .replace("fluidSolidInteraction;", "solid;"))

    # Material: the reference shear modulus with the requested nu
    mechanical = case / "constant" / "mechanicalProperties"
    nu = spec["nu"]
    bulk = 2.0 * MU_REFERENCE * (1.0 + nu) / (3.0 * (1.0 - 2.0 * nu))
    text = mechanical.read_text()
    text, found = re.subn(r"(mu\s+mu\s+\[[^]]*\]\s+)[^;]+;",
                          rf"\g<1>{MU_REFERENCE:.10g};", text)
    text, found2 = re.subn(r"(K\s+K\s+\[[^]]*\]\s+)[^;]+;",
                           rf"\g<1>{bulk:.10g};", text)
    if found != 1 or found2 != 1:
        fail(f"Could not set mu and K in {mechanical}")
    mechanical.write_text(text)

    solid_properties = case / "constant" / "solidProperties"
    if spec["solid"] == "highOrder":
        text = solid_properties.read_text()
        header = text[:text.index("solidModel")]
        body = HIGH_ORDER_PROPERTIES.replace(
            "highOrderJacobian true;",
            f"highOrderJacobian {'true' if spec.get('hoj') else 'false'};")
        body = re.sub(r"faceStencilExtraCells\s+\d+;",
                      f"faceStencilExtraCells {spec.get('extra', 65)};", body)
        solid_properties.write_text(header + body.lstrip())
        set_displacement_field(case / "0", "D")
        solution = case / "system" / "fvSolution"
        # Accept an inexact linear solve rather than stopping: the Krylov
        # solver does not always reach its tolerance for this thin flap
        replace_entry(solution, "ksp_max_it", '"200"')
        replace_text(solution, '            snes_monitor;\n',
                     '            snes_monitor;\n'
                     '            snes_max_linear_solve_fail "1000";\n')
    else:
        text = solid_properties.read_text()
        model = SOLID_MODELS[spec.get("model", "tutorial")]
        if model:
            text = re.sub(r"nonLinearGeometry(UpdatedLagrangian|"
                          r"TotalLagrangianTotalDisplacement)\b", model, text)
            set_displacement_field(
                case / "0", "DD" if "Updated" in model else "D")
        solid_properties.write_text(text)
        replace_entry(solid_properties, "scaleFactor", f"{spec['sf']:g}")
    replace_text(solid_properties, "    solutionAlgorithm PETScSNES;\n",
                 "    solutionAlgorithm PETScSNES;\n\n"
                 f"    dampingCoeff    [0 0 -1 0 0 0 0] {CAL_DAMPING:g};\n")

    # First-order implicit time scheme: only the static limit matters
    schemes = case / "system" / "fvSchemes"
    text = schemes.read_text()
    text, found = re.subn(r"(d2dt2Schemes|ddtSchemes)(\s*\{\s*default\s+)"
                          r"[^;]+;", r"\1\2Euler;", text)
    if found != 2:
        fail(f"Could not set the time schemes in {schemes}")
    schemes.write_text(text)

    if spec["cores"] > 1:
        use_parallel_lu(case / "system" / "fvSolution")
    replace_entry(case / "system" / "decomposeParDict", "numberOfSubdomains",
                  str(spec["cores"]), True)

    block = case / "system" / "blockMeshDict"
    for key, value in zip(("nx", "ny", "nz"), spec["mesh"]):
        replace_entry(block, key, str(value), True)

    steps = (LOAD_STEPS
             .replace("LEVELS", ", ".join(f"{v:g}" for v in LOAD_LEVELS))
             .replace("PLATEAU", f"{PLATEAU_TIME:g}")
             .replace("TRANSITION", f"{TRANSITION_TIME:g}"))
    control = (TUTORIAL / "system" / "controlDict").read_text()
    control = re.sub(r"^functions\s*\{.*\Z", "", control,
                     flags=re.MULTILINE | re.DOTALL)
    control += ("functions\n{\n" + steps + """
    flapTip
    {
        type            solidPointDisplacement;
        point           (0 0 0.065);
    }
}
""")
    control_path = case / "system" / "controlDict"
    control_path.write_text(control)
    replace_entry(control_path, "startFrom", "startTime", True)
    replace_entry(control_path, "endTime", f"{end_time:.10g}", True)
    replace_entry(control_path, "deltaT", f"{CAL_DELTA_T:.10g}", True)
    replace_entry(control_path, "writeControl", "timeStep", True)
    replace_entry(control_path, "writeInterval",
                  str(round(end_time / CAL_DELTA_T)), True)
    return case


def run_solid_case(case: Path, cores: int) -> None:
    run(["blockMesh"], case, "log.blockMesh")
    if cores > 1:
        run(["decomposePar"], case, "log.decomposePar")
        run(["mpirun", "-np", str(cores), "solids4Foam", "-parallel"], case,
            "log.solids4Foam")
    else:
        run(["solids4Foam"], case, "log.solids4Foam")


def monitor_history(case: Path, name: str, end_time: float,
                    delta_t: float) -> list[list[float]]:
    """History of one solidPointDisplacement monitor.

    A restarted run writes one file per start time; the segments are merged
    in start-time order, a later segment replacing overlapping times.
    """
    candidates = sorted(
        case.glob(f"postProcessing/*/solidPointDisplacement_{name}.dat"),
        key=lambda path: float(path.parent.name))
    if not candidates:
        fail(f"No '{name}' history in {case}")
    rows = []
    for path in candidates:
        # A restart supersedes everything the earlier segments wrote after
        # its start time
        start = float(path.parent.name)
        rows = [row for row in rows if row[0] <= start + 1e-12]
        rows.extend(numeric_rows(path))
    rows = last_per_time(rows)
    require_history(candidates[-1], rows, end_time, delta_t, 5)
    return rows


def tip_history(case: Path, end_time: float, delta_t: float
                ) -> list[list[float]]:
    return monitor_history(case, "flapTip", end_time, delta_t)


def calibration_runs(args: argparse.Namespace) -> list[dict]:
    """Every run of the calibration study."""
    meshes = [SOLID_MESHES[level] for level in
              (args.levels or ["coarse", "medium", "fine"])]
    if args.quick:
        meshes = meshes[:1]

    def run(**changes) -> dict:
        spec = {"solid": "standard", "model": "tutorial",
                "mesh": SOLID_MESHES["medium"], "nu": 0.45,
                "sf": TUTORIAL_STABILISATION, "cores": 1}
        spec.update(changes)
        return spec

    runs = [run(mesh=m) for m in meshes]
    if args.quick:
        return runs
    runs += [run(sf=0.05), run(sf=0.001), run(nu=0.49),
             run(model="ul", sf=0.05),
             # Cubic high-order solid with the compact Jacobian and a larger
             # least-squares stencil
             run(solid="highOrder", mesh=(6, 4, 33), sf=0.1, extra=60)]
    cores = int(args.cores) if args.cores != "auto" else 8
    if cores > 1:
        runs.append(run(cores=cores))
    return runs


def run_name(spec: dict) -> str:
    solid = spec["solid"] if spec["solid"] == "highOrder" else \
        f"standard{spec['model'].upper() if spec['model'] != 'tutorial' else ''}"
    extra = f"_x{spec['extra']}" if spec["solid"] == "highOrder" else ""
    return (f"calibration_{solid}{extra}_{'x'.join(map(str, spec['mesh']))}"
            f"_nu{spec['nu']:g}_sf{spec['sf']:g}_np{spec['cores']}")


def plateau_tips(history: list[list[float]]) -> list[tuple[float, float]]:
    """(load factor, settled tip y in mm) at the end of every plateau."""
    points = []
    for k, level in enumerate(LOAD_LEVELS):
        end = (k + 1) * PLATEAU_TIME
        rows = [r for r in history if end - 0.1 <= r[0] <= end + 1e-9]
        if not rows:
            fail(f"No tip history at the end of plateau {k + 1}")
        drift = 1.0e3 * (max(r[2] for r in rows) - min(r[2] for r in rows))
        if drift > 0.05:
            fail(f"The tip has not settled on plateau {k + 1}: it moves by "
                 f"{drift:.3f} mm over the last 0.1 s")
        points.append((level, 1.0e3 * rows[-1][2]))
    return points


def run_calibration(args: argparse.Namespace, reference: dict) -> bool:
    global LOAD_LEVELS, PLATEAU_TIME
    if args.quick:
        # Smoke test: one short plateau at the Phase I load, no checks
        LOAD_LEVELS = (1.0,)
        PLATEAU_TIME = 0.5
    target = float(reference["targets"]["zero_flow_tip_y_mm"])
    uniaxial = float(reference["materials_phase_I"]
                     ["uniaxial_neo_hookean_c1_Pa"])
    cheart = float(reference["materials_phase_I"]["cheart_neo_hookean_mu_Pa"])
    end_time = len(LOAD_LEVELS) * PLATEAU_TIME
    rows = []
    passed = True
    specs = calibration_runs(args)

    def settings_of(spec: dict) -> dict:
        settings = {
            "study": "calibration", **spec,
            "load_levels": list(LOAD_LEVELS), "plateau_time": PLATEAU_TIME,
            "transition_time": TRANSITION_TIME,
            "driver": driver_fingerprint(),
            "delta_t": CAL_DELTA_T, "damping": CAL_DAMPING,
            "mu_reference": MU_REFERENCE,
            "tutorial_inputs": tutorial_fingerprint(),
            "build": build_fingerprint(),
        }
        settings["mesh"] = list(spec["mesh"])
        return settings

    def execute(spec: dict) -> None:
        name = run_name(spec)
        label = f"calibration {name}"
        settings = settings_of(spec)
        if reusable(WORK_ROOT / name, settings, args, label):
            return
        case = solid_case(name, spec, end_time)
        print(f"Running {label}", flush=True)
        run_solid_case(case, spec["cores"])
        check_solver_log(case, label)
        (case / SETTINGS_FILE).write_text(json.dumps(settings, indent=2)
                                          + "\n")

    # The runs are independent: run up to --jobs of them at once
    with ThreadPoolExecutor(max_workers=max(1, args.jobs)) as pool:
        for future in [pool.submit(execute, spec) for spec in specs]:
            future.result()

    for spec in specs:
        name = run_name(spec)
        case = WORK_ROOT / name
        history = tip_history(case, end_time, CAL_DELTA_T)
        write_csv(OUTPUT_ROOT / f"{name}_tip_history.csv",
                  [{"time_s": f"{r[0]:.6g}", "tip_dy_mm": f"{1e3 * r[2]:.6g}",
                    "tip_dz_mm": f"{1e3 * r[3]:.6g}"} for r in history])
        if args.quick:
            print(f"  {name}: tip y = {1e3 * history[-1][2]:.3f} mm at "
                  f"t = {history[-1][0]:g} s (quick run, no checks)")
            continue
        points = plateau_tips(history)
        # Equivalent (mu, tip) pairs, in increasing tip order
        curve = sorted((MU_REFERENCE / load, tip) for load, tip in points)
        curve.sort(key=lambda item: item[1])
        mu_cal = interpolate([(tip, mu) for mu, tip in curve], target)
        row = {
            "solid": run_name(spec).split("_")[1],
            "mesh": "x".join(map(str, spec["mesh"])),
            "cells": math.prod(spec["mesh"]), "nu": spec["nu"],
            "stabilisation": spec["sf"], "ranks": spec["cores"],
            "mu_calibrated_Pa": mu_cal,
            "tip_mm_at_61kPa": interpolate(
                [(mu, tip) for mu, tip in sorted(curve)], MU_REFERENCE),
            "tip_mm_at_uniaxial": extrapolate_tip(curve, uniaxial),
            "clock_time_s": clock_time(case),
        }
        for load, tip in points:
            row[f"tip_mm_load{load:g}"] = tip
        rows.append(row)
        print(f"  {name}: mu = {mu_cal:.5g} Pa", flush=True)
        if not math.isfinite(mu_cal):
            passed = False
    if args.quick:
        return True
    write_csv(OUTPUT_ROOT / "calibration.csv",
              [{k: (f"{v:.5g}" if isinstance(v, float) else v)
                for k, v in r.items()} for r in rows])
    plot_calibration(rows, target, cheart)
    lines = ["# Calibration (flap alone under buoyancy)", "",
             f"Target zero-flow tip deflection: {target:.2f} mm; CHeart "
             f"used {cheart / 1e3:.0f} kPa, the uniaxial neo-Hookean fit "
             f"gives {uniaxial / 1e3:.1f} kPa.", "",
             "| solid | mesh | cells | nu | stab. | ranks | calibrated mu "
             "(kPa) | tip at 61 kPa (mm) | time (s) |",
             "| --- | --- | --- | --- | --- | --- | --- | --- | --- |"]
    for row in rows:
        lines.append(
            f"| {row['solid']} | {row['mesh']} | {row['cells']} | "
            f"{row['nu']:g} | {row['stabilisation']:g} | {row['ranks']} | "
            f"{row['mu_calibrated_Pa'] / 1e3:.2f} | "
            f"{row['tip_mm_at_61kPa']:.2f} | {row['clock_time_s']:.0f} |")
    for row in [r for r in rows if r["ranks"] > 1]:
        serial = [r for r in rows if r["ranks"] == 1 and r["mesh"] ==
                  row["mesh"] and r["nu"] == row["nu"]
                  and r["solid"] == row["solid"]
                  and r["stabilisation"] == row["stabilisation"]]
        if serial:
            a = serial[0]["mu_calibrated_Pa"]
            b = row["mu_calibrated_Pa"]
            same = math.isfinite(a) and abs(a - b) <= 1e-3 * a
            passed &= same
            lines += ["", f"Parallel check ({row['ranks']} ranks): mu = "
                      f"{b:.5g} Pa against {a:.5g} Pa in serial: "
                      f"{'PASS' if same else 'FAIL'}"]
    write_summary("calibration", lines)
    return passed


def interpolate(curve: list[tuple[float, float]], x: float) -> float:
    """Linear interpolation of an increasing (x, y) curve."""
    for (x0, y0), (x1, y1) in zip(curve, curve[1:]):
        if x0 <= x <= x1 and x1 > x0:
            return y0 + (y1 - y0) * (x - x0) / (x1 - x0)
    return math.nan


def extrapolate_tip(curve: list[tuple[float, float]], mu: float) -> float:
    """Tip at a shear modulus outside the plateaus, from the last two."""
    ordered = sorted(curve)
    (m0, t0), (m1, t1) = ordered[-2], ordered[-1]
    # The tip is close to linear in 1/mu over this range
    return t1 + (t0 - t1) * (1 / mu - 1 / m1) / (1 / m0 - 1 / m1)


def plot_calibration(rows: list[dict], target: float, cheart: float) -> None:
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib is unavailable: skipping the calibration plot")
        return
    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    ax.set_prop_cycle(color=plt.cm.tab20.colors)
    for row in rows:
        mus = [MU_REFERENCE / load / 1e3 for load in LOAD_LEVELS]
        tips = [row[f"tip_mm_load{load:g}"] for load in LOAD_LEVELS]
        style = "--" if row["ranks"] > 1 else "-o"
        ax.plot(mus, tips, style, ms=3,
                label=f"{row['solid']} {row['mesh']}, nu = {row['nu']:g}, "
                f"stab. {row['stabilisation']:g}"
                + (f", {row['ranks']} ranks" if row["ranks"] > 1 else ""))
    ax.axhline(target, color="k", lw=0.8)
    ax.axvline(cheart / 1e3, color="0.5", lw=0.8, ls=":")
    ax.text(cheart / 1e3, ax.get_ylim()[0], " CHeart", fontsize=7,
            va="bottom")
    ax.set_xlabel("shear modulus mu (kPa)")
    ax.set_ylabel("static tip position y (mm)")
    ax.set_title("zero-flow tip deflection; measured 29.50 mm", fontsize=9)
    ax.legend(fontsize=6)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(OUTPUT_ROOT / "calibration_tip_vs_mu.png", dpi=150)
    plt.close(fig)


def read_csv_plain(path: Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle))


# ---------------------------------------------------------------------------
# Phase I study: steady inflow
# ---------------------------------------------------------------------------

def ignored(directory: str, names: list[str]) -> set[str]:
    """Copy only the inputs a fresh run needs, never a previous result."""
    skip = {"validation", "regressionTests", "postProcessing", "dynamicCode",
            "case.foam"}
    skip.update(name for name in names if name.startswith("processor"))
    skip.update(name for name in names if name.startswith("log."))
    directory_path = Path(directory)
    if directory_path == TUTORIAL:
        skip.update(name for name in names
                    if name != "0" and re.fullmatch(r"[0-9.eE+-]+", name))
    if directory_path.name in ("fluid", "solid") \
            and directory_path.parent.name == "constant":
        skip.add("polyMesh")
    if directory_path == TUTORIAL / "constant":
        skip.update({"polyMesh", "triSurface"})
    return skip.intersection(names)


def prepare_phase_i(name: str, cores: int, end_time: float,
                    delta_t: float) -> Path:
    case = WORK_ROOT / name
    if case.exists():
        shutil.rmtree(case)
    shutil.copytree(TUTORIAL, case, symlinks=True, ignore=ignored)
    control = case / "system" / "controlDict"
    replace_entry(control, "startFrom", "startTime", True)
    replace_entry(control, "endTime", f"{end_time:.10g}", True)
    replace_entry(control, "deltaT", f"{delta_t:.10g}", True)
    replace_entry(control, "writeInterval", f"{end_time:.10g}", True)
    set_subdomains(case, cores, ("", "fluid", "solid"))
    if cores > 1:
        use_parallel_lu(case / "system" / "solid" / "fvSolution")
    return case


def voxel_points(velocity: list[dict[str, float]]) -> list[tuple]:
    """Sub-sampling points (m) of every MRI voxel, voxel index first."""
    sizes = (1.302e-3, 1.302e-3, 6.0e-3)
    points = []
    for index, row in enumerate(velocity):
        centre = (1e-3 * row["x_mm"], 1e-3 * row["y_mm"], 1e-3 * row["z_mm"])
        for i in range(VOXEL_SAMPLES[0]):
            for j in range(VOXEL_SAMPLES[1]):
                for k in range(VOXEL_SAMPLES[2]):
                    offsets = [((n + 0.5) / VOXEL_SAMPLES[d] - 0.5) * sizes[d]
                               for d, n in enumerate((i, j, k))]
                    points.append((index, tuple(c + o for c, o in
                                                zip(centre, offsets))))
    return points


def write_sampling_dict(case: Path, points: list[tuple]) -> Path:
    lines = "\n".join(f"                ({p[0]:.8g} {p[1]:.8g} {p[2]:.8g})"
                      for _, p in points)
    path = case / "system" / "sampleVoxels"
    path.write_text(f"""FoamFile
{{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      sampleVoxels;
}}

functions
{{
    sampleVoxels
    {{
        type            sets;
        libs            (sampling);
        region          fluid;
        interpolationScheme cellPoint;
        setFormat       raw;
        fields          (U);
        sets
        {{
            voxels
            {{
                type    cloud;
                axis    xyz;
                points
                (
{lines}
                );
            }}
        }}
    }}
}}
""")
    return path


def latest_time(case: Path, cores: int) -> str:
    base = case / ("processor0" if cores > 1 else "")
    times = [p.name for p in base.iterdir()
             if p.is_dir() and re.fullmatch(r"[0-9.eE+-]+", p.name)
             and p.name != "0"]
    if not times:
        fail(f"No result time in {base}")
    return max(times, key=float)


def sample_voxels(case: Path, cores: int, velocity: list[dict],
                  end_time: float) -> list[dict]:
    """Average the computed velocity over every MRI voxel."""
    points = voxel_points(velocity)
    write_sampling_dict(case, points)
    command = ["postProcess", "-region", "fluid", "-dict",
               "system/sampleVoxels", "-fields", "(U)", "-latestTime"]
    if cores > 1:
        command = ["mpirun", "-np", str(cores)] + command + ["-parallel"]
    run(command, case, "log.sampleVoxels")
    time = latest_time(case, cores)
    if not math.isclose(float(time), end_time, rel_tol=1e-9, abs_tol=1e-9):
        fail(f"The latest field time of {case} is {time}, not the end "
             f"time {end_time:g}")
    files = sorted(path for path in
                   (case / "postProcessing").rglob("voxels*U*")
                   if "sampleVoxels" in path.parts
                   and path.parent.name == time)
    if not files:
        fail(f"No voxel samples written in {case}")
    samples = numeric_rows(files[0])
    # Match the written samples to the requested points on a 10 um grid:
    # the sample points are at least 0.3 mm apart, and the written
    # coordinates are rounded
    def key(point) -> tuple:
        return tuple(round(v * 1.0e5) for v in point)

    lookup = {}
    for row in samples:
        lookup[key(row[:3])] = row[3:6]
    sums = [[0.0, 0.0, 0.0, 0] for _ in velocity]
    for index, point in points:
        value = lookup.get(key(point))
        if value is not None:
            for d in range(3):
                sums[index][d] += value[d]
            sums[index][3] += 1
    per_voxel = math.prod(VOXEL_SAMPLES)
    averaged = []
    for row, total in zip(velocity, sums):
        found = total[3]
        entry = dict(row)
        entry["fluid_fraction"] = found / per_voxel
        for key in ("vx", "vy", "vz"):
            entry[f"{key}_sim_mm_s"] = None
        # Compare only voxels that lie mostly in the computed fluid
        if found >= 0.5 * per_voxel:
            for d, key in enumerate(("vx", "vy", "vz")):
                entry[f"{key}_sim_mm_s"] = 1.0e3 * total[d] / found
        averaged.append(entry)
    return averaged


def centreline(case: Path, end_time: float, delta_t: float
               ) -> list[tuple[float, float]]:
    """Deformed centreline (z, y) in mm of the flap at the final time.

    Built from the displacement monitors on the undeformed centreline
    (system/flapCentreline) and the tip monitor.
    """
    points = {0.0: (0.0, 0.0)}
    names = sorted({
        path.name[len("solidPointDisplacement_"):-len(".dat")]
        for path in case.glob(
            "postProcessing/*/solidPointDisplacement_flapCentrelineZ*.dat")
    })
    if len(names) < 10:
        fail(f"Centreline monitors missing in {case}")
    for name in names + ["flapTip"]:
        match = re.search(r"Z(\d+)$", name)
        z0 = 0.1 * int(match.group(1)) if match else 65.0
        rows = monitor_history(case, name, end_time, delta_t)
        points[round(z0, 3)] = (z0 + 1e3 * rows[-1][3], 1e3 * rows[-1][2])
    return sorted(points.values())


def centreline_errors(sim: list[tuple[float, float]],
                      measured: list[dict[str, float]]) -> dict:
    errors = []
    for row in measured:
        y_sim = interpolate(sim, row["z_mm"])
        if math.isfinite(y_sim):
            errors.append(y_sim - row["y_mm"])
    if not errors:
        fail("The computed centreline does not overlap the measurement")
    return {
        "points": len(errors),
        "rms_mm": math.sqrt(sum(e * e for e in errors) / len(errors)),
        "max_mm": max(abs(e) for e in errors),
        "mean_mm": sum(errors) / len(errors),
    }


def velocity_errors(voxels: list[dict], scale: float) -> dict:
    valid = [v for v in voxels if v["vx_sim_mm_s"] is not None]
    if not valid:
        fail("No voxel lies in the computed fluid domain")
    distances = [math.sqrt(sum((v[f"{k}_sim_mm_s"] - v[f"{k}_mm_s"]) ** 2
                               for k in ("vx", "vy", "vz"))) for v in valid]
    result = {"voxels": len(valid), "voxels_total": len(voxels)}
    for plane in (10.0, 30.0):
        d = [dist for dist, v in zip(distances, valid)
             if abs(v["z_mm"] - plane) < 1e-6]
        if d:
            result[f"d_mean_z{plane:g}"] = sum(d) / len(d) / scale
            result[f"d_max_z{plane:g}"] = max(d) / scale
    result["d_mean"] = sum(distances) / len(distances) / scale
    result["d_max"] = max(distances) / scale
    for key in ("vx", "vy", "vz"):
        diff = [v[f"{key}_sim_mm_s"] - v[f"{key}_mm_s"] for v in valid]
        result[f"rms_{key}_mm_s"] = math.sqrt(sum(x * x for x in diff)
                                              / len(diff))
    return result


def fsi_iterations(case: Path) -> tuple[float, int]:
    """Mean and maximum FSI iterations per time step."""
    text = (case / "log.solids4Foam").read_text(errors="replace")
    counts: dict[str, int] = {}
    for time, iteration in re.findall(
            r"^Time = ([0-9.eE+-]+), iteration: (\d+)", text, re.MULTILINE):
        counts[time] = max(counts.get(time, 0), int(iteration))
    if not counts:
        return math.nan, 0
    values = list(counts.values())
    return sum(values) / len(values), max(values)


def run_phase_i(args: argparse.Namespace, reference: dict) -> bool:
    levels = args.levels or ["coarse", "medium"]
    if args.quick:
        levels = levels[:1]
    end_time = 0.02 if args.quick else args.end_time
    delta_t = args.delta_t
    cores = int(args.cores) if args.cores != "auto" else 16
    measured_line = read_csv(REFERENCE_DIR / "phaseI_centreline.csv")
    velocity = read_csv(REFERENCE_DIR / "phaseI_velocity.csv")
    scale = float(reference["measurement"]["velocity_normalisation_mm_s"])
    tip_target = float(reference["targets"]["phase_I_tip_y_mm"])
    rows = []
    passed = True
    def execute(level: str) -> None:
        name = f"phaseI_{level}"
        label = f"Phase I {level}"
        settings = {
            "study": "phaseI", "level": level, "cores": cores,
            "end_time": end_time, "delta_t": delta_t,
            "driver": driver_fingerprint(),
            "tutorial_inputs": tutorial_fingerprint(),
            "build": build_fingerprint(),
        }
        if reusable(WORK_ROOT / name, settings, args, label):
            return
        case = prepare_phase_i(name, cores, end_time, delta_t)
        print(f"Running {label} on {cores} ranks", flush=True)
        run(["./Allrun", level] + (["parallel"] if cores > 1 else []),
            case, "log.Allrun")
        check_solver_log(case, label)
        (case / SETTINGS_FILE).write_text(json.dumps(settings, indent=2)
                                          + "\n")

    external = {}
    if args.case:
        # Evaluate a completed run made elsewhere from a copy of this
        # tutorial; its own controlDict and decomposition apply
        case = Path(args.case).resolve()
        check_solver_log(case, f"Phase I run in {case}")
        control = (case / "system" / "controlDict").read_text()
        end_time = float(re.search(r"^endTime\s+([^;]+);", control,
                                   re.MULTILINE).group(1))
        delta_t = float(re.search(r"^deltaT\s+([^;]+);", control,
                                  re.MULTILINE).group(1))
        cores = sum(1 for p in case.iterdir()
                    if re.fullmatch(r"processor\d+", p.name)) or 1
        levels = levels[:1]
        external[levels[0]] = case
    else:
        with ThreadPoolExecutor(max_workers=max(1, args.jobs)) as pool:
            for future in [pool.submit(execute, level) for level in levels]:
                future.result()

    for level in levels:
        name = f"phaseI_{level}"
        case = external.get(level, WORK_ROOT / name)
        history = tip_history(case, end_time, delta_t)
        write_csv(OUTPUT_ROOT / f"{name}_tip_history.csv",
                  [{"time_s": f"{r[0]:.6g}", "tip_dx_mm": f"{1e3 * r[1]:.6g}",
                    "tip_dy_mm": f"{1e3 * r[2]:.6g}",
                    "tip_dz_mm": f"{1e3 * r[3]:.6g}"} for r in history])
        mean_it, max_it = fsi_iterations(case)
        row = {"level": level, "ranks": cores,
               "tip_y_mm": 1e3 * history[-1][2],
               "tip_z_mm": 65.0 + 1e3 * history[-1][3],
               "fsi_iterations_mean": mean_it, "fsi_iterations_max": max_it,
               "clock_time_h": clock_time(case) / 3600.0}
        # Change of the tip over the last tenth of the run
        tail = [r[2] for r in history if r[0] >= 0.9 * end_time]
        row["tip_drift_last_tenth_mm"] = 1e3 * (max(tail) - min(tail))
        if not args.quick:
            line = centreline(case, end_time, delta_t)
            write_csv(OUTPUT_ROOT / f"{name}_centreline.csv",
                      [{"z_mm": f"{z:.5g}", "y_mm": f"{y:.5g}"}
                       for z, y in line])
            row.update({f"centreline_{k}": v for k, v in
                        centreline_errors(line, measured_line).items()})
            voxels = sample_voxels(case, cores, velocity, end_time)
            write_csv(OUTPUT_ROOT / f"{name}_voxels.csv", voxels)
            row.update(velocity_errors(voxels, scale))
            tip_error = abs(row["tip_y_mm"] - tip_target)
            # Steady: the tip moves by less than 0.1 mm over the last tenth
            ok = (tip_error <= 1.0 and row["centreline_rms_mm"] <= 1.0
                  and row["tip_drift_last_tenth_mm"] <= 0.1)
            row["pass"] = ok
            passed &= ok
        rows.append(row)
    write_csv(OUTPUT_ROOT / "phaseI_summary.csv",
              [{k: (f"{v:.5g}" if isinstance(v, float) else v)
                for k, v in r.items()} for r in rows])
    if not args.quick:
        plot_phase_i(rows, measured_line, reference)
    lines = ["# Phase I (steady inflow)", ""]
    for row in rows:
        lines.append(f"## {row['level']} fluid mesh, {row['ranks']} ranks")
        lines.append("")
        for key, value in row.items():
            if key in ("level", "ranks"):
                continue
            text = f"{value:.4g}" if isinstance(value, float) else str(value)
            lines.append(f"- {key}: {text}")
        lines.append("")
    write_summary("phaseI", lines)
    return passed


def plot_phase_i(rows: list[dict], measured_line: list[dict],
                 reference: dict) -> None:
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib is unavailable: skipping the Phase I plots")
        return
    fit = reference["phase_I_centreline_fit"]
    zero = reference["zero_flow_centreline_fit"]

    def poly(fitting: dict, z: float) -> float:
        s = z / fitting["zhat_mm"]
        return fitting["yhat_mm"] * sum(p * s ** (i + 1)
                                        for i, p in enumerate(fitting["p"]))

    # Centreline
    fig, ax = plt.subplots(figsize=(6.4, 3.6))
    ax.plot([r["z_mm"] for r in measured_line],
            [r["y_mm"] for r in measured_line], "ko", ms=3,
            label="experiment (MRI)")
    zs = [zero["zhat_mm"] * i / 100 for i in range(101)]
    ax.plot(zs, [poly(zero, z) for z in zs], "k:", lw=0.8,
            label="experiment, zero flow (fit)")
    loz = read_csv(REFERENCE_DIR / "Lozovskiy2019_fig9_upper_surface.csv")
    ax.plot([r["z_mm"] for r in loz], [r["y_upper_mm"] - 1.0 for r in loz],
            "c-.", lw=1, label="Lozovskiy et al. (2019), digitised")
    published = reference["published_phase_I_tip_y_mm"]
    ax.plot([62.5], [published["cheart_inf_sup_stable_rl3"]], "m^", ms=5,
            label="CHeart inf-sup stable, tip (read from figure)")
    ax.plot([62.5], [published["cheart_cG1cG1"]], "mv", ms=5,
            label="CHeart cG(1)cG(1), tip (read from figure)")
    for row in rows:
        data = read_csv(OUTPUT_ROOT / f"phaseI_{row['level']}_centreline.csv")
        ax.plot([d["z_mm"] for d in data], [d["y_mm"] for d in data],
                label=f"solids4foam, {row['level']} mesh")
    ax.set_xlabel("z (mm)")
    ax.set_ylabel("centreline y (mm)")
    ax.legend(fontsize=7)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(OUTPUT_ROOT / "phaseI_centreline.png", dpi=150)
    plt.close(fig)

    # Tip history
    fig, ax = plt.subplots(figsize=(6.4, 3.2))
    for row in rows:
        data = read_csv(OUTPUT_ROOT / f"phaseI_{row['level']}_tip_history.csv")
        ax.plot([d["time_s"] for d in data], [d["tip_dy_mm"] for d in data],
                label=f"solids4foam, {row['level']} mesh")
    ax.axhline(reference["targets"]["phase_I_tip_y_mm"], color="k", lw=0.8,
               label="experiment, 16.41 mm")
    ax.set_xlabel("time (s)")
    ax.set_ylabel("tip displacement y (mm)")
    ax.legend(fontsize=7)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(OUTPUT_ROOT / "phaseI_tip_history.png", dpi=150)
    plt.close(fig)

    # Velocity profiles along y at x = 0 on both planes
    fig, axes = plt.subplots(2, 3, figsize=(10, 5.6), sharex=True)
    for row in rows:
        voxels = read_csv_plain(OUTPUT_ROOT / f"phaseI_{row['level']}_voxels.csv")
        for i, plane in enumerate((10.0, 30.0)):
            line = sorted((float(v["y_mm"]), v) for v in voxels
                          if abs(float(v["z_mm"]) - plane) < 1e-6
                          and abs(float(v["x_mm"])) < 1e-6)
            for j, key in enumerate(("vx", "vy", "vz")):
                ax = axes[i][j]
                if row is rows[0]:
                    ax.errorbar([y for y, _ in line],
                                [float(v[f"{key}_mm_s"]) for _, v in line],
                                yerr=0.05 * reference["measurement"]
                                ["venc_mm_s"][key], fmt="ko", ms=3,
                                capsize=2, label="experiment")
                sim = [(y, float(v[f"{key}_sim_mm_s"])) for y, v in line
                       if v.get(f"{key}_sim_mm_s")]
                ax.plot([y for y, _ in sim], [s for _, s in sim], "-",
                        label=f"{row['level']} mesh")
                ax.set_title(f"{key}, z = {plane:g} mm, x = 0", fontsize=9)
                ax.grid(alpha=0.3)
    for ax in axes[1]:
        ax.set_xlabel("y (mm)")
    for ax in axes[:, 0]:
        ax.set_ylabel("velocity (mm/s)")
    axes[0][0].legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(OUTPUT_ROOT / "phaseI_velocity_profiles.png", dpi=150)
    plt.close(fig)


def write_summary(stem: str, lines: list[str]) -> None:
    path = OUTPUT_ROOT / f"{stem}_summary.md"
    path.write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--study", choices=("calibration", "phaseI"),
                        default="calibration")
    parser.add_argument(
        "--levels",
        help="comma-separated mesh levels: solid meshes coarse,medium,fine "
             "for calibration (default coarse,medium), fluid meshes "
             "coarse,medium for phaseI (default coarse,medium)")
    parser.add_argument("--cores", default="auto",
                        help="MPI ranks: the parallel calibration run and "
                             "every Phase I run (default 8 and 16)")
    parser.add_argument("--end-time", type=float, default=10.0,
                        help="Phase I end time in s (default 10)")
    parser.add_argument("--delta-t", type=float, default=0.002,
                        help="Phase I time step in s (default 0.002)")
    parser.add_argument("--jobs", type=int, default=1,
                        help="runs executed at once (default 1)")
    parser.add_argument("--case",
                        help="phaseI: evaluate this completed run of the "
                             "first --levels mesh instead of running it")
    parser.add_argument("--quick", action="store_true",
                        help="smoke test: first level only, short runs, "
                             "no accuracy checks")
    parser.add_argument("--reuse", action="store_true",
                        help="re-evaluate completed runs under "
                             "validation/work")
    args = parser.parse_args()
    if args.cores != "auto" and (not args.cores.isdecimal()
                                 or int(args.cores) < 1):
        parser.error("--cores must be a positive integer or auto")
    if args.levels:
        args.levels = [value.strip() for value in args.levels.split(",")]
        allowed = (SOLID_MESHES if args.study == "calibration"
                   else ("coarse", "medium"))
        if any(level not in allowed for level in args.levels):
            parser.error(f"--levels must be from {', '.join(allowed)}")
    for executable in ("blockMesh", "solids4Foam"):
        if shutil.which(executable) is None:
            fail(f"Required executable '{executable}' is unavailable. "
                 "Source OpenFOAM and build solids4foam first.")
    reference = json.loads(REFERENCE_FILE.read_text())
    WORK_ROOT.mkdir(parents=True, exist_ok=True)
    OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)
    if args.study == "calibration":
        passed = run_calibration(args, reference)
    else:
        passed = run_phase_i(args, reference)
    return 0 if passed else 1


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (RuntimeError, subprocess.SubprocessError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        sys.exit(2)
