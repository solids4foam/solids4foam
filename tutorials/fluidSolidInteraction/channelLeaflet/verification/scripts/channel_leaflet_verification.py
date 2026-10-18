#!/usr/bin/env python3
"""Run the channelLeaflet code-to-code comparison with oomph-lib.

Each case is a copy of the tutorial under verification/work/ with a different
time step, solid discretisation, fluid mesh or leaflet thickness. The history
of the leaflet tip displacement is compared with oomph-lib solutions of the
same problem (reference/channelLeaflet_oomph_*.csv). The default studies
(static, comparison, time) complete and have acceptance checks; the optional
studies (mesh, thickness) are exploratory, are run only when named with
--study, and do not currently complete (see ../README.md).
"""

from __future__ import annotations

import argparse
import bisect
import csv
import hashlib
import inspect
import json
import math
import os
import re
import shutil
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent
VERIFY_DIR = SCRIPT_DIR.parent
CASE_DIR = VERIFY_DIR.parent
WORK_DIR = VERIFY_DIR / "work"
POST_DIR = VERIFY_DIR / "postProcessing"
REFERENCE_DIR = VERIFY_DIR / "reference"
REFERENCE_JSON = REFERENCE_DIR / "channelLeaflet_verification_references.json"

MONITOR = "tip"

# Geometry: leaflet root x0, leaflet height, channel height, downstream length
X0 = 1.0
LEAFLET = 0.5
CHANNEL = 1.0
LDOWN = 7.0

# Tutorial leaflet: thickness, plane-strain modulus E/(1 - nu^2), Poisson's
# ratio and density
H_TUTORIAL = 0.05
E_EFF = 5000.0
NU = 0.3
RHO = 1.0

FINGERPRINT_FILE = "fingerprint.dat"
HEADER = (
    "FoamFile\n{\n    version     2.0;\n    format      ascii;\n"
    "    class       dictionary;\n    object      blockMeshDict;\n}\n\n"
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run the channelLeaflet verification studies"
    )
    parser.add_argument(
        "--study",
        default="default",
        help="comma-separated studies: static, comparison, time, mesh,"
        " thickness; default (the default) runs the non-optional ones and"
        " all runs every study",
    )
    parser.add_argument(
        "--quick",
        action="store_true",
        help="smoke test: the coarsest case of each study to a short end time",
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
        help="number of cases run at the same time, each in serial (default 1)",
    )
    parser.add_argument(
        "--cases",
        help="comma-separated case names to run, from the selected studies",
    )
    parser.add_argument(
        "--keep-going",
        action="store_true",
        help="continue after an individual case fails",
    )
    parser.add_argument(
        "--list",
        action="store_true",
        help="list the cases of the selected studies and exit",
    )
    args = parser.parse_args()
    if args.cores < 1:
        parser.error("--cores must be at least 1")
    return args


def load_references() -> dict:
    return json.loads(REFERENCE_JSON.read_text())


# ---------------------------------------------------------------------------
# Cases
# ---------------------------------------------------------------------------

def fluid_divisions(level: float) -> dict[str, int]:
    """Block divisions of the fluid mesh at a refinement level: level 1 is
    the tutorial mesh, and level n divides every block n times as finely"""
    return {
        "thick": int(round(4*level)),     # across the leaflet top
        "height": int(round(40*level)),   # along the leaflet, and above it
        "up": int(round(40*level)),       # upstream of the leaflet
        "down1": int(round(100*level)),   # fine downstream section
        "down2": int(round(60*level)),    # coarse downstream section
    }


def case_name(case: dict) -> str:
    name = (
        f"F{case['fluid']:g}_S{case['solid'][0]}x{case['solid'][1]}"
        f"_dt{1.0/case['dt']:g}_{case['solidType']}"
    )
    if abs(case["h"] - H_TUTORIAL) > 1e-12:
        name += f"_h{case['h']:g}"
    return name


def make_case(spec: dict, **overrides) -> dict:
    case = {
        "fluid": spec.get("fluidLevel", 1),
        "solid": tuple(spec.get("solidMesh", (4, 40))),
        "dt": spec.get("deltaT", 0.01),
        "solidType": spec.get("solidType", "highOrder"),
        "h": spec.get("h", H_TUTORIAL),
    }
    case.update(overrides)
    case["solid"] = tuple(case["solid"])
    case["name"] = case_name(case)
    return case


def study_cases(study: str, refs: dict, quick: bool) -> list[dict]:
    spec = refs["studies"][study]
    cases = []
    if study == "static":
        meshes = spec["quickSolidMeshes"] if quick else spec["solidMeshes"]
        for h, solid_type, solid in (
            (h, solid_type, solid)
            for h in spec.get("thickness", [H_TUTORIAL])
            for solid_type in spec["solidTypes"]
            for solid in meshes.get(f"{h:g}", [])
        ):
            name = f"static_S{solid[0]}x{solid[1]}_{solid_type}"
            if abs(h - H_TUTORIAL) > 1e-12:
                name += f"_h{h:g}"
            cases.append({
                "name": name,
                "solid": tuple(solid),
                "solidType": solid_type,
                "h": h,
                "static": True,
            })
    elif study == "mesh":
        levels = spec["quickLevels"] if quick else spec["levels"]
        for level in levels:
            cases.append(make_case(spec, fluid=level))
    elif study == "time":
        steps = spec["quickDeltaT"] if quick else spec["deltaT"]
        for dt in steps:
            cases.append(make_case(spec, dt=dt))
    elif study == "comparison":
        meshes = spec["quickSolidMeshes"] if quick else spec["solidMeshes"]
        for solid_type in spec["solidTypes"]:
            for solid in meshes:
                cases.append(
                    make_case(spec, solid=tuple(solid), solidType=solid_type)
                )
    elif study == "thickness":
        values = spec["quickThickness"] if quick else spec["thickness"]
        for h in values:
            cases.append(make_case(spec, h=h))
    else:
        raise RuntimeError(f"unknown study {study}")
    return cases


# ---------------------------------------------------------------------------
# Meshes and case set-up
# ---------------------------------------------------------------------------

def expansion_ratio(length: float, n: int, first: float) -> float:
    """blockMesh expansion ratio (last/first cell) of n geometric cells over
    `length` whose first cell is `first`"""
    if abs(n*first - length) <= 1e-12*length:
        return 1.0
    growing = n*first < length
    lo, hi = (1.0 + 1e-14, 2.0) if growing else (0.5, 1.0 - 1e-14)
    for _ in range(200):
        r = 0.5*(lo + hi)
        total = first*(r**n - 1.0)/(r - 1.0)
        if (total > length) == growing:
            hi = r
        else:
            lo = r
    return r**(n - 1)


def fluid_block_mesh(h: float, level: float) -> str:
    """The fluid domain around a leaflet of thickness h; `level` multiplies
    every block division. Level 1 with h = 0.05 is the tutorial mesh."""
    a, b, x3 = X0 - h/2, X0 + h/2, X0 + LDOWN
    n = fluid_divisions(level)
    n_thick, n_height, n_up = n["thick"], n["height"], n["up"]
    # Cell width next to the leaflet faces: the geometric mean of the cells
    # across the leaflet top and the cells along it
    d_near = math.sqrt((h/4)*(LEAFLET/40))/level
    r_up = 1.0/expansion_ratio(a, n_up, d_near)
    # Downstream: a fine section over two channel heights, then a coarse one
    n_d1, n_d2 = n["down1"], n["down2"]
    l_d1 = 2.0
    r_d1 = expansion_ratio(l_d1, n_d1, d_near)
    d_end1 = d_near*r_d1
    r_d2 = expansion_ratio(
        x3 - b - l_d1, n_d2, d_end1*r_d1**(1.0/(n_d1 - 1))
    )
    xs = (0.0, a, b, x3)
    ys = (0.0, LEAFLET, CHANNEL)

    def v(i, j, k=0):
        return i + 4*j + 12*k

    lines = ["// Channel of unit height, 0 < x < 8, with the leaflet of height 0.5",
             f"// and thickness {h:g} centred on x = 1 cut out of the fluid domain",
             "scale 1;", "", "vertices", "("]
    for k, z in enumerate((0.0, 0.1)):
        for j, y in enumerate(ys):
            for i, x in enumerate(xs):
                lines.append(f"    ({x:.10g} {y:.10g} {z:.10g}) // {v(i, j, k)}")
    lines += [");", "", "blocks", "("]
    down = (
        f"((({l_d1:.10g} {n_d1} {r_d1:.8g}) "
        f"({x3 - b - l_d1:.10g} {n_d2} {r_d2:.8g})) 1 1)"
    )
    for (i, j), cells, grading in (
        ((0, 0), (n_up, n_height), f"({r_up:.8g} 1 1)"),
        ((0, 1), (n_up, n_height), f"({r_up:.8g} 1 1)"),
        ((1, 1), (n_thick, n_height), "(1 1 1)"),
        ((2, 0), (n_d1 + n_d2, n_height), down),
        ((2, 1), (n_d1 + n_d2, n_height), down),
    ):
        ids = (v(i, j), v(i + 1, j), v(i + 1, j + 1), v(i, j + 1),
               v(i, j, 1), v(i + 1, j, 1), v(i + 1, j + 1, 1), v(i, j + 1, 1))
        lines.append(
            f"    hex ({' '.join(map(str, ids))}) ({cells[0]} {cells[1]} 1)"
            f" simpleGrading {grading}"
        )
    lines += [");", "", "edges", "(", ");", "", "boundary", "("]

    def quad(*ids):
        return "            (" + " ".join(map(str, ids)) + ")"

    empty = []
    for i, j in ((0, 0), (0, 1), (1, 1), (2, 0), (2, 1)):
        empty.append(quad(v(i, j), v(i, j + 1), v(i + 1, j + 1), v(i + 1, j)))
        empty.append(
            quad(v(i, j, 1), v(i + 1, j, 1), v(i + 1, j + 1, 1), v(i, j + 1, 1))
        )
    patches = (
        ("inlet", "patch", [
            quad(v(0, 0), v(0, 0, 1), v(0, 1, 1), v(0, 1)),
            quad(v(0, 1), v(0, 1, 1), v(0, 2, 1), v(0, 2))]),
        ("outlet", "patch", [
            quad(v(3, 0), v(3, 1), v(3, 1, 1), v(3, 0, 1)),
            quad(v(3, 1), v(3, 2), v(3, 2, 1), v(3, 1, 1))]),
        ("interface", "wall", [
            quad(v(1, 0), v(1, 1), v(1, 1, 1), v(1, 0, 1)),
            quad(v(1, 1), v(2, 1), v(2, 1, 1), v(1, 1, 1)),
            quad(v(2, 0), v(2, 0, 1), v(2, 1, 1), v(2, 1))]),
        ("walls", "wall", [
            quad(v(0, 0), v(1, 0), v(1, 0, 1), v(0, 0, 1)),
            quad(v(2, 0), v(3, 0), v(3, 0, 1), v(2, 0, 1)),
            quad(v(0, 2), v(0, 2, 1), v(1, 2, 1), v(1, 2)),
            quad(v(1, 2), v(1, 2, 1), v(2, 2, 1), v(2, 2)),
            quad(v(2, 2), v(2, 2, 1), v(3, 2, 1), v(3, 2))]),
        ("frontAndBack", "empty", empty),
    )
    for name, kind, faces in patches:
        lines += [f"    {name}", "    {", f"        type {kind};",
                  "        faces", "        ("] + faces + ["        );", "    }"]
    lines += [");", "", "mergePatchPairs", "(", ");"]
    return "\n".join(lines) + "\n"


def solid_block_mesh(h: float, cells: tuple[int, int], patches: str) -> str:
    """The leaflet, `cells` = (across, along); `patches` is "fsi" (one
    interface patch) or "static" (separate faces)"""
    a, b = X0 - h/2, X0 + h/2
    lines = [f"// Leaflet of thickness {h:g} and height 0.5 centred on x = 1,",
             "// clamped at its base, y = 0", "scale 1;", "", "vertices", "("]
    for z in (0.0, 0.1):
        for x, y in ((a, 0.0), (b, 0.0), (b, LEAFLET), (a, LEAFLET)):
            lines.append(f"    ({x:.10g} {y:.10g} {z:.10g})")
    lines += [");", "", "blocks", "(",
              f"    hex (0 1 2 3 4 5 6 7) ({cells[0]} {cells[1]} 1)"
              " simpleGrading (1 1 1)", ");", "", "edges", "(", ");", "",
              "boundary", "("]
    if patches == "fsi":
        groups = (("interface", ["(0 4 7 3)", "(3 7 6 2)", "(1 2 6 5)"]),)
    else:
        groups = (("upstream", ["(0 4 7 3)"]), ("top", ["(3 7 6 2)"]),
                  ("downstream", ["(1 2 6 5)"]))
    for name, faces in groups + (("base", ["(0 1 5 4)"]),):
        lines += [f"    {name}", "    {", "        type patch;", "        faces",
                  "        ("] + [f"            {f}" for f in faces] + [
                  "        );", "    }"]
    lines += ["    frontAndBack", "    {", "        type empty;", "        faces",
              "        (", "            (0 3 2 1)", "            (4 5 6 7)",
              "        );", "    }", ");", "", "mergePatchPairs", "(", ");"]
    return "\n".join(lines) + "\n"


def leaflet_properties(h: float) -> tuple[float, float]:
    """Young's modulus and density of a leaflet of thickness h with the
    bending stiffness and mass per unit length of the tutorial leaflet"""
    scale = H_TUTORIAL/h
    return E_EFF*(1.0 - NU*NU)*scale**3, RHO*scale


def ignored(directory: str, names: list[str]) -> set[str]:
    """Copy only the inputs a fresh run needs, never a previous result."""
    ignored_names = {
        "verification", "regressionTests", "postProcessing", "images",
        "dynamicCode",
    }
    ignored_names.update(n for n in names if n.startswith("processor"))
    ignored_names.update(n for n in names if n.startswith("log."))
    directory_path = Path(directory)
    if directory_path == CASE_DIR:
        ignored_names.update(n for n in names if is_result_time(n))
        ignored_names.update({"README.md", "regressionTest.sh"})
    if directory_path.name in {"fluid", "solid"} and directory_path.parent.name == "constant":
        ignored_names.add("polyMesh")
    return ignored_names.intersection(names)


def is_result_time(name: str) -> bool:
    """A time directory other than 0, in any OpenFOAM time format"""
    try:
        value = float(name)
    except ValueError:
        return False
    return math.isfinite(value) and name != "0"


def solver_fingerprint() -> str:
    """The solver on the path and the solids4foam libraries it loads: the
    models and the radial-basis-function mesh motion, wherever the build put
    them"""
    parts = []
    solver = shutil.which("solids4Foam")
    directories = [
        os.environ.get(name, "")
        for name in ("FOAM_MODULE_LIBBIN", "FOAM_USER_LIBBIN", "FOAM_SITE_LIBBIN")
    ]
    libraries = [
        str(Path(d) / f"lib{name}{suffix}")
        for d in directories if d
        for name in ("solids4FoamModels", "RBFMeshMotionSolver")
        for suffix in (".so", ".dylib")
    ]
    for path in [solver] + libraries:
        if path and Path(path).is_file():
            stat = Path(path).stat()
            parts.append(f"{path}:{stat.st_size}:{stat.st_mtime_ns}")
    parts.append(os.environ.get("WM_PROJECT_VERSION", ""))
    return ";".join(parts)


def fingerprint(case: dict, end_time: float) -> str:
    """Hash of the tutorial inputs, the case settings and the build.

    A reused case must have been run from the same inputs: the files a fresh
    copy takes from the tutorial (not generated plots, backups or the
    variant links Allrun creates), the functions of this script that set up
    and run the cases, the solver and solids4foam library found in the
    environment, the case parameters and the end time. The evaluation of the
    results is not part of it, so that it can change without re-running.
    """
    digest = hashlib.sha256()
    for directory, dirs, files in os.walk(CASE_DIR):
        skipped = ignored(directory, dirs + files)
        dirs[:] = sorted(d for d in dirs if d not in skipped)
        for name in sorted(f for f in files if f not in skipped):
            path = Path(directory) / name
            if (
                path.is_symlink() or name.endswith((".bak", ".pdf"))
                or name.endswith(".withDefaultValues")
            ):
                continue
            digest.update(str(path.relative_to(CASE_DIR)).encode())
            digest.update(path.read_bytes())
    for function in SETUP_FUNCTIONS:
        digest.update(inspect.getsource(function).encode())
    digest.update(
        repr((X0, LEAFLET, CHANNEL, LDOWN, H_TUTORIAL, E_EFF, NU, RHO,
              HEADER)).encode()
    )
    digest.update(solver_fingerprint().encode())
    digest.update(json.dumps(case, sort_keys=True, default=str).encode())
    digest.update(f"{end_time:.10g}".encode())
    return digest.hexdigest()


def replace_once(path: Path, pattern: str, replacement: str) -> None:
    text, count = re.subn(pattern, replacement, path.read_text(), count=1,
                          flags=re.MULTILINE)
    if count != 1:
        raise RuntimeError(f"pattern {pattern!r} not found in {path}")
    path.write_text(text)


def set_leaflet_material(path: Path, h: float) -> None:
    modulus, density = leaflet_properties(h)
    replace_once(path, r"^(\s*E\s+E\s+\[[^]]*\]\s+)[^;]+;", rf"\g<1>{modulus:.10g};")
    replace_once(path, r"^(\s*rho\s+rho\s+\[[^]]*\]\s+)[^;]+;", rf"\g<1>{density:.10g};")


def configure_case(run_dir: Path, case: dict, end_time: float) -> None:
    (run_dir / "system" / "fluid" / "blockMeshDict").write_text(
        HEADER + fluid_block_mesh(case["h"], case["fluid"])
    )
    (run_dir / "system" / "solid" / "blockMeshDict").write_text(
        HEADER + solid_block_mesh(case["h"], case["solid"], "fsi")
    )
    set_leaflet_material(run_dir / "constant" / "solid" / "mechanicalProperties",
                         case["h"])
    control = run_dir / "system" / "controlDict"
    replace_once(control, r"^(deltaT\s+)[^;]+;", rf"\g<1>{case['dt']:.10g};")
    replace_once(control, r"^(endTime\s+)[^;]+;", rf"\g<1>{end_time:.10g};")
    replace_once(control, r"^(writeInterval\s+)[^;]+;", rf"\g<1>{end_time:.10g};")
    # Matching interface faces can be mapped directly; otherwise use AMI
    n = fluid_divisions(case["fluid"])
    if tuple(case["solid"]) != (n["thick"], n["height"]):
        replace_once(
            run_dir / "constant" / "fsiProperties",
            r"^(\s*interfaceTransferMethod\s+)\w+;",
            r"\g<1>AMI;",
        )


# ---------------------------------------------------------------------------
# Running and checking
# ---------------------------------------------------------------------------

def numeric_rows(path: Path) -> list[list[float]]:
    rows = []
    for line in path.read_text(errors="replace").splitlines():
        fields = line.replace("(", " ").replace(")", " ").split()
        if not fields or fields[0].startswith("#"):
            continue
        try:
            rows.append([float(field) for field in fields])
        except ValueError:
            # A malformed line is kept so that the finiteness check fails
            rows.append([math.nan])
    return rows


def read_history(run_dir: Path) -> list[tuple[float, float, float]]:
    """(time, tip x displacement, tip y displacement)"""
    candidates = sorted(
        run_dir.glob(f"postProcessing/**/solidPointDisplacement_{MONITOR}.dat")
    )
    if not candidates:
        raise RuntimeError(f"no {MONITOR} displacement data in {run_dir}")
    # A restart writes a second file; keep the latest value at each time
    values: dict[float, tuple[float, float]] = {}
    for path in candidates:
        for row in numeric_rows(path):
            if len(row) < 3:
                raise RuntimeError(f"{path}: truncated line")
            values[round(row[0], 10)] = (row[1], row[2])
    return [(t, dx, dy) for t, (dx, dy) in sorted(values.items())]


def check_run(run_dir: Path, name: str, end_time: float, step: float) -> None:
    """Raise unless the run ended cleanly with a complete, finite history."""
    log = run_dir / "log.solids4Foam"
    text = log.read_text(errors="replace") if log.is_file() else ""
    if not re.search(r"^End\s*$", text, re.MULTILINE):
        raise RuntimeError(f"{name} did not finish; see {run_dir}")
    if not (math.isfinite(step) and step > 0 and math.isfinite(end_time)):
        raise RuntimeError(f"{name}: invalid time step or end time")
    expected = int(round(end_time/step))
    if abs(expected*step - end_time) > 1e-6*end_time:
        raise RuntimeError(f"{name}: the end time is not a multiple of the step")
    history = read_history(run_dir)
    if not all(all(math.isfinite(v) for v in point) for point in history):
        raise RuntimeError(f"{name}: the {MONITOR} history is not finite")
    steps = {
        int(round(t/step)) for t, _, _ in history
        if t > 0 and abs(t/step - round(t/step)) < 1e-4
    }
    missing = set(range(1, expected + 1)) - steps
    if missing:
        raise RuntimeError(
            f"{name}: the {MONITOR} history is incomplete: {len(missing)} of"
            f" {expected} steps missing, from t = {min(missing)*step:g}"
        )


def completed(run_dir: Path, name: str, end_time: float, step: float,
              expected: str) -> bool:
    stamp = run_dir / FINGERPRINT_FILE
    if not stamp.is_file() or stamp.read_text().strip() != expected:
        return False
    try:
        check_run(run_dir, name, end_time, step)
    except RuntimeError:
        return False
    return True


def run_case(case: dict, end_time: float, reuse: bool) -> Path:
    run_dir = WORK_DIR / case["name"]
    stamp = fingerprint(case, end_time)
    if reuse and completed(run_dir, case["name"], end_time, case["dt"], stamp):
        print(f"  reusing {case['name']}")
        return run_dir
    if run_dir.exists():
        shutil.rmtree(run_dir)
    shutil.copytree(CASE_DIR, run_dir, ignore=ignored, symlinks=True)
    configure_case(run_dir, case, end_time)

    command = ["./Allrun"]
    if case["solidType"] == "linear":
        command.append("linear")
    print(f"  running {case['name']}", flush=True)
    with (run_dir / "log.Allverify").open("w") as handle:
        result = subprocess.run(
            command, cwd=run_dir, stdout=handle, stderr=subprocess.STDOUT
        )
    if result.returncode:
        raise RuntimeError(f"{case['name']} failed; see {run_dir}")
    check_run(run_dir, case["name"], end_time, case["dt"])
    (run_dir / FINGERPRINT_FILE).write_text(stamp + "\n")
    return run_dir


def run_static_case(case: dict, refs: dict, reuse: bool) -> Path:
    """Solve the leaflet alone under a small uniform traction, raised in
    steps, on its upstream face."""
    spec = refs["studies"]["static"]
    steps = spec["loadSteps"]
    run_dir = WORK_DIR / case["name"]
    # The applied traction is part of the fingerprint: the reference
    # deflection is computed from it
    stamp = fingerprint({**case, "traction": spec["traction"]}, steps)
    if reuse and completed(run_dir, case["name"], steps, 1.0, stamp):
        print(f"  reusing {case['name']}")
        return run_dir
    if run_dir.exists():
        shutil.rmtree(run_dir)
    for sub in ("0", "constant", "system"):
        (run_dir / sub).mkdir(parents=True)
    for name in ("dynamicMeshDict", "g", "mechanicalProperties"):
        shutil.copy(CASE_DIR / "constant" / "solid" / name, run_dir / "constant")
    set_leaflet_material(run_dir / "constant" / "mechanicalProperties", case["h"])
    shutil.copy(
        CASE_DIR / "constant" / "solid" / f"solidProperties.{case['solidType']}",
        run_dir / "constant" / "solidProperties",
    )
    shutil.copy(CASE_DIR / "system" / "solid" / "fvSchemes", run_dir / "system")
    shutil.copy(
        CASE_DIR / "system" / "solid" / f"fvSolution.{case['solidType']}",
        run_dir / "system" / "fvSolution",
    )
    (run_dir / "system" / "blockMeshDict").write_text(
        HEADER + solid_block_mesh(case["h"], case["solid"], "static")
    )
    (run_dir / "constant" / "physicsProperties").write_text(
        "FoamFile { version 2.0; format ascii; class dictionary;"
        " object physicsProperties; }\ntype solid;\n"
    )
    load = spec["traction"]
    (run_dir / "constant" / "load").write_text(
        f"(\n    (0 (0 0 0))\n    ({steps} ({load:.10g} 0 0))\n)\n"
    )
    free = "type solidTraction; traction uniform (0 0 0); pressure uniform 0;" \
        " value uniform (0 0 0);"
    (run_dir / "0" / "D").write_text(
        "FoamFile { version 2.0; format ascii; class volVectorField; object D; }\n"
        "dimensions [0 1 0 0 0 0 0];\ninternalField uniform (0 0 0);\n"
        "boundaryField\n{\n"
        "    upstream { type solidTraction; tractionSeries { fileName"
        ' "$FOAM_CASE/constant/load"; outOfBounds clamp; } pressure uniform 0;'
        " value uniform (0 0 0); }\n"
        f"    top {{ {free} }}\n    downstream {{ {free} }}\n"
        "    base { type fixedDisplacement; value uniform (0 0 0); }\n"
        "    frontAndBack { type empty; }\n}\n"
    )
    (run_dir / "0" / "pointD").write_text(
        "FoamFile { version 2.0; format ascii; class pointVectorField;"
        " object pointD; }\n"
        "dimensions [0 1 0 0 0 0 0];\ninternalField uniform (0 0 0);\n"
        "boundaryField\n{\n"
        "    upstream { type calculated; value uniform (0 0 0); }\n"
        "    top { type calculated; value uniform (0 0 0); }\n"
        "    downstream { type calculated; value uniform (0 0 0); }\n"
        "    base { type fixedValue; value uniform (0 0 0); }\n"
        "    frontAndBack { type empty; }\n}\n"
    )
    (run_dir / "system" / "controlDict").write_text(
        "FoamFile { version 2.0; format ascii; class dictionary;"
        " object controlDict; }\n"
        "application solids4Foam;\nstartFrom startTime;\nstartTime 0;\n"
        f"stopAt endTime;\nendTime {steps};\ndeltaT 1;\n"
        f"writeControl timeStep;\nwriteInterval {steps};\n"
        "writeFormat ascii;\nwritePrecision 10;\ntimeFormat general;\n"
        "timePrecision 6;\nrunTimeModifiable no;\n"
        f"functions\n{{\n    {MONITOR} {{ type solidPointDisplacement;"
        f" point ({X0:g} {LEAFLET:g} 0.05); }}\n}}\n"
    )
    # Run through the tutorial run functions, which also restore the library
    # path that macOS strips from child processes
    allrun = run_dir / "Allrun"
    allrun.write_text(
        "#!/bin/bash\n"
        ". $WM_PROJECT_DIR/bin/tools/RunFunctions\n"
        "source solids4FoamScripts.sh\n"
        "solids4Foam::runApplication blockMesh\n"
        "solids4Foam::runApplication solids4Foam\n"
    )
    allrun.chmod(0o755)
    print(f"  running {case['name']}", flush=True)
    with (run_dir / "log.Allverify").open("w") as handle:
        subprocess.run(
            ["./Allrun"], cwd=run_dir, stdout=handle, stderr=subprocess.STDOUT
        )
    check_run(run_dir, case["name"], steps, 1.0)
    (run_dir / FINGERPRINT_FILE).write_text(stamp + "\n")
    return run_dir


# ---------------------------------------------------------------------------
# Evaluation
# ---------------------------------------------------------------------------

def read_reference(path: Path) -> list[tuple[float, float, float]]:
    with path.open() as handle:
        rows = list(csv.DictReader(
            line for line in handle if not line.startswith("#")
        ))
    curve = [
        (float(row["time"]), float(row["tipDx"]), float(row["tipDy"]))
        for row in rows
    ]
    if not curve or not all(all(math.isfinite(v) for v in p) for p in curve):
        raise RuntimeError(f"{path}: empty or non-finite reference")
    return curve


def interpolate(curve: list[tuple[float, float, float]], time: float
                ) -> tuple[float, float]:
    times = [point[0] for point in curve]
    if not times[0] - 1e-9 <= time <= times[-1] + 1e-9:
        raise RuntimeError(
            f"t = {time:g} is outside the reference, {times[0]:g} to {times[-1]:g}"
        )
    index = bisect.bisect_left(times, time - 1e-12)
    if index <= 0:
        return curve[0][1], curve[0][2]
    if index >= len(curve):
        return curve[-1][1], curve[-1][2]
    (t0, x0, y0), (t1, x1, y1) = curve[index - 1], curve[index]
    w = (time - t0)/(t1 - t0)
    return x0 + w*(x1 - x0), y0 + w*(y1 - y0)


def amplitude(curve: list[tuple[float, float, float]]) -> float:
    """The largest tip displacement of a reference"""
    return max(math.hypot(dx, dy) for _, dx, dy in curve)


def period_statistics(history, start: float, end: float) -> dict:
    """Mean, minimum and maximum of the tip displacements over [start, end]"""
    inside = [p for p in history if start - 1e-9 <= p[0] <= end + 1e-9]
    if len(inside) < 3:
        raise RuntimeError(f"too few samples between t = {start:g} and {end:g}")
    stats = {}
    for index, label in ((1, "dx"), (2, "dy")):
        values = [p[index] for p in inside]
        # Trapezoidal mean over the period
        mean = sum(
            0.5*(a[index] + b[index])*(b[0] - a[0])
            for a, b in zip(inside, inside[1:])
        )/(inside[-1][0] - inside[0][0])
        stats[f"{label}Mean"] = mean
        stats[f"{label}Min"] = min(values)
        stats[f"{label}Max"] = max(values)
    return stats


def history_difference(history, reference, start: float, end: float,
                       scale: float) -> tuple[float, float]:
    """Largest and RMS difference of the tip displacement vector over
    [start, end], relative to `scale`"""
    errors = []
    for t, dx, dy in history:
        if start - 1e-9 <= t <= end + 1e-9:
            rx, ry = interpolate(reference, t)
            errors.append(math.hypot(dx - rx, dy - ry))
    if not errors:
        raise RuntimeError("no samples to compare")
    return (
        max(errors)/scale,
        math.sqrt(sum(e*e for e in errors)/len(errors))/scale,
    )


def solver_statistics(run_dir: Path) -> dict:
    text = (run_dir / "log.solids4Foam").read_text(errors="replace")
    iterations: dict[str, int] = {}
    for time, iteration in re.findall(
        r"^Time = (\S+), iteration: (\d+)", text, re.MULTILINE
    ):
        iterations[time] = max(iterations.get(time, 0), int(iteration))
    clock = re.findall(r"ClockTime = ([0-9.eE+-]+) s", text)
    counts = list(iterations.values())
    return {
        "steps": len(counts),
        "meanFsiIterations": sum(counts)/len(counts) if counts else math.nan,
        "maxFsiIterations": max(counts) if counts else 0,
        "clockTime": float(clock[-1]) if clock else math.nan,
    }


def cell_counts(run_dir: Path) -> tuple[int, int]:
    counts = []
    for region in ("fluid", "solid"):
        log = run_dir / f"log.blockMesh.{region}"
        match = re.search(r"nCells:\s*(\d+)", log.read_text(errors="replace"))
        counts.append(int(match.group(1)) if match else 0)
    return counts[0], counts[1]


def reference_path(refs: dict, spec: dict, case: dict) -> Path:
    """The oomph-lib reference a case is compared with: the converged one,
    the one at the case's own time step, or the one for its leaflet
    thickness"""
    kind = spec.get("reference", "converged")
    if kind == "converged":
        return REFERENCE_DIR / refs["reference"]["file"]
    if kind == "thickness":
        return REFERENCE_DIR / refs["reference"]["thickness"][f"{case['h']:g}"]
    return REFERENCE_DIR / refs["reference"]["sameDeltaT"].format(
        steps=f"{1.0/case['dt']:g}"
    )


def evaluate(case: dict, run_dir: Path, refs: dict, spec: dict,
             end_time: float) -> dict:
    reference = read_reference(reference_path(refs, spec, case))
    scale = amplitude(read_reference(REFERENCE_DIR / refs["reference"]["file"]))
    history = [p for p in read_history(run_dir) if 0 < p[0] <= end_time + 1e-9]
    period = refs["period"]
    row = {
        "case": case["name"],
        "fluidLevel": case["fluid"],
        "solidCells": f"{case['solid'][0]}x{case['solid'][1]}",
        "deltaT": case["dt"],
        "solidType": case["solidType"],
        "h": case["h"],
    }
    row["fluidCells"], row["solidCellCount"] = cell_counts(run_dir)
    row["maxError"], row["rmsError"] = history_difference(
        history, reference, 0.0, end_time, scale
    )
    start = end_time - period
    row["periodMaxError"], _ = history_difference(
        history, reference, start, end_time, scale
    )
    stats = period_statistics(history, start, end_time)
    ref_stats = period_statistics(reference, start, end_time)
    for key, value in stats.items():
        row[key] = value
        row[f"{key}Error"] = (value - ref_stats[key])/scale
    # Periodicity: the change of the tip history over the last period
    if end_time >= 2*period:
        previous = [(t + period, dx, dy) for t, dx, dy in history
                    if start - period - 1e-9 <= t <= start + 1e-9]
        row["periodicity"], _ = history_difference(
            [p for p in history if p[0] >= start - 1e-9], previous, start,
            end_time, scale
        )
    row.update(solver_statistics(run_dir))
    require_finite(row)
    return row


def evaluate_static(case: dict, run_dir: Path, refs: dict) -> dict:
    spec = refs["studies"]["static"]
    # Kirchhoff-Love (Euler-Bernoulli) cantilever under a uniform load, whose
    # bending stiffness per unit width, E_eff h^3/12, is that of the tutorial
    # leaflet for every thickness of the study
    stiffness = E_EFF*H_TUTORIAL**3/12.0
    reference = spec["traction"]*LEAFLET**4/(8.0*stiffness)
    value = read_history(run_dir)[-1][1]
    text = (run_dir / "log.solids4Foam").read_text(errors="replace")
    clock = re.findall(r"ClockTime = ([0-9.eE+-]+) s", text)
    row = {
        "case": case["name"],
        "solidCells": f"{case['solid'][0]}x{case['solid'][1]}",
        "solidType": case["solidType"],
        "h": case["h"],
        "tipDx": value,
        # Positive: the leaflet deflects more than the Kirchhoff-Love beam
        "relError": (value - reference)/reference,
        "clockTime": float(clock[-1]) if clock else math.nan,
    }
    require_finite(row)
    return row


def require_finite(row: dict) -> None:
    for key, value in row.items():
        if isinstance(value, float) and not math.isfinite(value):
            if key == "clockTime":
                continue
            raise RuntimeError(f"{row['case']}: {key} is not finite")


def self_convergence(results, spacing, scale: float, start: float = 0.0):
    """Differences between successive cases of a refinement study over
    [start, end]; `spacing` gives the mesh spacing or time step of a case,
    by which the cases are ordered from coarse to fine. No order of
    convergence is computed: the studies have too few completed levels."""
    ordered = sorted(results, key=lambda item: -spacing(item[0]))
    lines, diffs = [], []
    for (case_a, dir_a, _), (case_b, dir_b, _) in zip(ordered, ordered[1:]):
        coarse = read_history(dir_a)
        fine = read_history(dir_b)
        end = min(coarse[-1][0], fine[-1][0])
        diff, _ = history_difference(coarse, fine, start, end, scale)
        diffs.append(diff)
        lines.append(
            f"- {case_a['name']} vs {case_b['name']}: max difference {diff:.3%}"
        )
    return lines, diffs


def write_rows(study: str, rows: list[dict]) -> Path:
    POST_DIR.mkdir(parents=True, exist_ok=True)
    path = POST_DIR / f"{study}_study.csv"
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
    return path


def describe(row: dict) -> str:
    if "relError" in row:
        return (
            f"{row['case']}: tip {row['tipDx']:.7f}, error"
            f" {row['relError']:+.3%}, {row['clockTime']:.0f} s"
        )
    return (
        f"{row['case']}: max error {row['maxError']:.3%}, last period"
        f" {row['periodMaxError']:.3%}, {row['meanFsiIterations']:.1f} FSI"
        f" iterations/step, {row['clockTime']:.0f} s"
    )


def summary_table(study: str, rows: list[dict]) -> list[str]:
    if study == "static":
        lines = ["| Case | Tip (m) | Error | Clock (s) |", "|---|---:|---:|---:|"]
        for row in rows:
            lines.append(
                f"| {row['case']} | {row['tipDx']:.7f} | {row['relError']:+.3%}"
                f" | {row['clockTime']:.0f} |"
            )
        return lines
    lines = [
        "| Case | Max error | RMS error | Last period | Mean dx error |"
        " Max dx error | FSI it./step | Clock (s) |",
        "|---|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in rows:
        lines.append(
            f"| {row['case']} | {row['maxError']:.3%} | {row['rmsError']:.3%} |"
            f" {row['periodMaxError']:.3%} | {row['dxMeanError']:+.3%} |"
            f" {row['dxMaxError']:+.3%} | {row['meanFsiIterations']:.1f} |"
            f" {row['clockTime']:.0f} |"
        )
    return lines


def check_study(study: str, results, refs: dict, quick: bool,
                scale: float) -> tuple[list[str], list[str]]:
    """Return the summary lines and the failed acceptance checks of a study"""
    rows = [row for _, _, row in results]
    criteria = refs["studies"][study].get("acceptance", {})
    lines, failures = [], []
    if quick:
        return lines, failures
    if study == "static":
        for row in rows:
            tolerance = criteria["maxRelError"].get(row["solidType"])
            if tolerance is not None and abs(row["relError"]) > tolerance:
                failures.append(
                    f"static: {row['case']} error {row['relError']:.3g}"
                    f" > {tolerance}"
                )
    elif study == "comparison":
        # Code-to-code comparison over the last period of the periodic state
        for case, _, row in results:
            lines.append(
                f"- {case['name']}: largest difference over the last period"
                f" {row['periodMaxError']:.3%}, over the whole run"
                f" {row['maxError']:.3%}; change of the last period from the"
                f" one before {row.get('periodicity', math.nan):.3%}"
            )
            tolerance = criteria["maxPeriodError"].get(case["solidType"])
            if tolerance is not None and row["periodMaxError"] > tolerance:
                failures.append(
                    f"comparison: {case['name']} differs from oomph-lib by"
                    f" {row['periodMaxError']:.3g} over the last period"
                    f" > {tolerance}"
                )
            if row.get("periodicity", math.inf) > criteria["maxPeriodicity"]:
                failures.append(
                    f"comparison: {case['name']} is not periodic: its last"
                    f" period changed by {row.get('periodicity', math.nan):.3g}"
                    f" > {criteria['maxPeriodicity']}"
                )
    elif study == "time":
        # Consistency of two time steps over the periodic state; with only
        # two steps that complete, no order of convergence is claimed
        start = refs["studies"]["time"]["comparisonStart"]
        lines, diffs = self_convergence(
            results, lambda case: case["dt"], scale, start
        )
        lines = [f"{line} (t >= {start:g} s)" for line in lines]
        if diffs and diffs[-1] > criteria["maxSelfDifference"]:
            failures.append(
                f"time: time-step difference {diffs[-1]:.3g}"
                f" > {criteria['maxSelfDifference']}"
            )
    elif study == "mesh":
        # Exploratory: the refined fluid meshes do not currently complete
        lines, _ = self_convergence(
            results, lambda case: 1.0/case["fluid"], scale
        )
    elif study == "thickness":
        # Exploratory: the difference from the beam model does not fall as
        # the leaflet is thinned, so no limit is extrapolated
        for case, _, row in results:
            lines.append(
                f"- h = {case['h']:g}: largest difference from oomph-lib at the"
                f" same bending stiffness {row['maxError']:.3%}, over the last"
                f" period {row['periodMaxError']:.3%}"
            )
    return lines, failures


def write_plot(study: str, cases, reference: Path, end_time: float) -> None:
    if not shutil.which("gnuplot"):
        return
    lines = [
        'set terminal pngcairo enhanced size 1000,700 font ",11"',
        f'set output "{POST_DIR / (study + "_tip.png")}"',
        'set multiplot layout 2,1',
        'set grid',
        'set key outside right',
        f'set xrange [0:{end_time:g}]',
    ]
    ref = f'"< grep -v ^# {reference} | tail -n +2 | tr , \' \'"'
    for column, label in ((2, "x"), (3, "y")):
        lines.append(f'set ylabel "tip {label} displacement (m)"')
        if column == 3:
            lines.append('set xlabel "t (s)"')
        plot = [f'plot {ref} using 1:{column} with lines lw 3 lc rgb "black"'
                ' title "oomph-lib"']
        for case, run_dir in cases:
            data = sorted(run_dir.glob(
                f"postProcessing/**/solidPointDisplacement_{MONITOR}.dat"
            ))
            if data:
                plot.append(
                    f'"{data[-1]}" using 1:{column} with lines lw 1.5'
                    f' title "{case["name"]}" noenhanced'
                )
        lines.append(", \\\n     ".join(plot))
    lines.append("unset multiplot")
    path = POST_DIR / f"{study}_tip.gnuplot"
    path.write_text("\n".join(lines) + "\n")
    subprocess.run(["gnuplot", str(path)], check=False)


# The functions whose source is part of a case's fingerprint: everything that
# sets up or runs a case
SETUP_FUNCTIONS = (
    fluid_divisions, expansion_ratio, fluid_block_mesh, solid_block_mesh,
    leaflet_properties, set_leaflet_material, replace_once, ignored,
    is_result_time, configure_case, run_case, run_static_case,
)


def main() -> int:
    args = parse_args()
    refs = load_references()
    all_studies = list(refs["studies"])
    if args.study == "all":
        studies = all_studies
    elif args.study == "default":
        studies = [
            name for name in all_studies
            if not refs["studies"][name].get("optional", False)
        ]
    else:
        studies = [s.strip() for s in args.study.split(",")]
    unknown = [s for s in studies if s not in refs["studies"]]
    if unknown:
        print(f"ERROR: unknown studies: {', '.join(unknown)}", file=sys.stderr)
        return 2
    end_time = refs["quickEndTime"] if args.quick else refs["endTime"]
    selected = set(args.cases.split(",")) if args.cases else None

    if args.list:
        for study in studies:
            for case in study_cases(study, refs, args.quick):
                print(f"{study}: {case['name']}")
        return 0

    if selected:
        known = {
            case["name"] for study in studies
            for case in study_cases(study, refs, args.quick)
        }
        unknown = selected - known
        if unknown:
            print(f"ERROR: unknown cases: {', '.join(sorted(unknown))}",
                  file=sys.stderr)
            return 2

    for tool in ("blockMesh", "solids4Foam"):
        if not shutil.which(tool):
            print(f"ERROR: {tool} not found; source OpenFOAM and solids4foam",
                  file=sys.stderr)
            return 2
    version = os.environ.get("WM_PROJECT_VERSION", "")
    if not version.startswith("v"):
        print(f"ERROR: the channelLeaflet study needs OpenFOAM.com"
              f" (tested with v2412); found {version or 'none'}",
              file=sys.stderr)
        return 2
    if not os.environ.get("PETSC_DIR"):
        print("ERROR: the channelLeaflet study needs solids4foam built"
              " with PETSc (PETSC_DIR is not set)", file=sys.stderr)
        return 2

    scale = amplitude(read_reference(REFERENCE_DIR / refs["reference"]["file"]))
    all_failures = []
    partial = False
    summary = ["# channelLeaflet verification summary", "",
               f"Errors are relative to the largest tip displacement of the"
               f" reference, {scale:.5f} m.", ""]

    # Run every case of the selected studies first, a case shared between
    # studies once, so that the cases of all studies share the --cores slots
    selections = {}
    unique: dict[str, dict] = {}
    for study in studies:
        all_cases = study_cases(study, refs, args.quick)
        cases = [
            case for case in all_cases
            if not selected or case["name"] in selected
        ]
        partial = partial or len(cases) < len(all_cases)
        selections[study] = (all_cases, cases)
        for case in cases:
            unique.setdefault(case["name"], case)

    def run_one(case: dict):
        try:
            if case.get("static"):
                return case["name"], run_static_case(case, refs, args.reuse)
            return case["name"], run_case(case, end_time, args.reuse)
        except RuntimeError as error:
            return case["name"], str(error)

    with ThreadPoolExecutor(max_workers=args.cores) as pool:
        finished = dict(pool.map(run_one, unique.values()))

    for study in studies:
        print(f"Study: {study}", flush=True)
        spec = refs["studies"][study]
        all_cases, cases = selections[study]
        results = []
        for case in cases:
            run_dir = finished[case["name"]]
            try:
                if isinstance(run_dir, str):
                    raise RuntimeError(run_dir)
                if case.get("static"):
                    row = evaluate_static(case, run_dir, refs)
                else:
                    row = evaluate(case, run_dir, refs, spec, end_time)
            except RuntimeError as error:
                print(f"  FAILED: {error}", flush=True)
                all_failures.append(str(error))
                continue
            results.append((case, run_dir, row))
            print(f"  {describe(row)}", flush=True)
        if all_failures and not args.keep_going:
            return 1
        if not results:
            continue
        if not args.quick and len(results) < len(all_cases) and not selected:
            all_failures.append(f"{study}: not every case completed")
        lines, failures = check_study(study, results, refs, args.quick, scale)
        for line in lines:
            print(f"  {line[2:]}", flush=True)
        rows = [row for _, _, row in results]
        path = write_rows(study, rows)
        if study != "static":
            write_plot(
                study,
                [(case, run_dir) for case, run_dir, _ in results],
                reference_path(refs, spec, results[-1][0]),
                end_time,
            )
        summary += [f"## {study}", "", f"Results: `{path.name}`", ""]
        summary += summary_table(study, rows)
        summary += [""] + lines + [""] + [f"- FAIL: {f}" for f in failures] + [""]
        all_failures.extend(failures)

    POST_DIR.mkdir(parents=True, exist_ok=True)
    (POST_DIR / "verification_summary.md").write_text("\n".join(summary) + "\n")
    if all_failures:
        print("Verification FAILED:")
        for failure in all_failures:
            print(f"  {failure}")
        return 1
    if partial:
        print("Verification INCOMPLETE: --cases ran only part of the studies")
        return 3
    print("Verification PASSED")
    return 0


if __name__ == "__main__":
    sys.exit(main())
