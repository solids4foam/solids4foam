#!/usr/bin/env python3
"""Prescribed-motion fluid-only ALE test on the HronTurek FSI3 fluid meshes.

The flag of the rigid CFD3 case is given a prescribed periodic bending: the
point displacement of the plate patch is eta(x, t) = A s(t) phi(xi) sin(2 pi f
(t - t0)), in the transverse (y) direction, with phi the normalised first
cantilever mode (zero displacement and slope at the cylinder), xi the
distance along the flag from the cylinder, s a smooth ramp from 0 to 1 over
`ramp` seconds from the start time t0, and A = 35 mm, f = 5.5 Hz (the
FSI3 tip amplitude and frequency). There is no solid and no coupling. The
mesh motion is the FSI3 one (velocityLaplacian, quadratic inverseDistance
diffusivity on plate and cylinder), the plate boundary velocity is
newMovingWallVelocity, and the point velocity is the finite difference
(eta(t) - eta(t - dt))/dt, so that the Euler-integrated plate boundary
follows eta exactly.

The run restarts from the developed rigid CFD3 state of the same level (the
cfd3_<level>x_dt<dt> run, time 25 s), so it needs that run to exist.

    python3 scripts/hron_turek_ale.py --level 2 --delta-t 0.0005 --duration 8 --cores 8
"""

from __future__ import annotations

import argparse
import re
import shutil
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import hron_turek_verification as driver  # noqa: E402
from hron_turek_cfd3 import case_name as cfd3_case_name  # noqa: E402

WORK_ROOT = driver.WORK_ROOT
TUTORIAL = driver.TUTORIAL
AMPLITUDE = 0.035
FREQUENCY = 5.5
RAMP = 1.0
START = 25.0
# Root of the flag (vertex 12 of the blockMeshDict) and tip, in metres
X_ROOT, X_TIP = 0.2489898, 0.6
# Clamped-free beam modes: beta_n L and sigma_n. FSI3 flaps at 5.5 Hz, about
# 2.5 times the first in-vacuum frequency, i.e. close to the second mode
MODES = {1: (1.875104068711961, 0.734095514), 2: (4.694091132974175, 1.018467319)}
MODE = 2

BC_CODE = r"""
        const scalar t = this->db().time().value();
        const scalar dt = this->db().time().deltaTValue();
        const scalar A = %(A)g, f = %(f)g, t0 = %(t0)g, ramp = %(ramp)g;
        const scalar xr = %(xr)g, L = %(xt)g - %(xr)g;
        const scalar beta = %(beta)s;
        const scalar sigma = %(sigma)s;
        auto phi = [&](scalar x) -> scalar
        {
            scalar xi = min(max((x - xr)/L, 0.0), 1.0);
            scalar b = beta*xi;
            return
            (
                cosh(b) - cos(b) - sigma*(sinh(b) - sin(b))
            )/(cosh(beta) - cos(beta) - sigma*(sinh(beta) - sin(beta)));
        };
        auto s = [&](scalar tt) -> scalar
        {
            scalar r = min(max((tt - t0)/ramp, 0.0), 1.0);
            return r*r*(3.0 - 2.0*r);
        };
        auto eta = [&](scalar x, scalar tt) -> scalar
        {
            return A*s(tt)*phi(x)*sin(2.0*constant::mathematical::pi*f*max(tt - t0, 0.0));
        };
        const pointField& pts = this->patch().localPoints();
        vectorField v(pts.size(), vector::zero);
        forAll(pts, i)
        {
            v[i] = vector(0, (eta(pts[i].x(), t) - eta(pts[i].x(), t - dt))/dt, 0);
        }
        operator==(v);
"""


GENERALISED_CODE = r"""
            const fvMesh& m = mesh();
            const label pid = m.boundaryMesh().findPatchID("plate");
            const volScalarField& p = m.lookupObject<volScalarField>("p");
            const scalarField& pp = p.boundaryField()[pid];
            const vectorField& Sf = m.Sf().boundaryField()[pid];
            const vectorField& Cf = m.Cf().boundaryField()[pid];
            const scalar xr = %(xr)g, L = %(xt)g - %(xr)g;
            const scalar beta = %(beta)s, sigma = %(sigma)s;
            auto phi = [&](scalar x) -> scalar
            {
                scalar xi = min(max((x - xr)/L, 0.0), 1.0);
                scalar b = beta*xi;
                return (cosh(b) - cos(b) - sigma*(sinh(b) - sin(b)))
                    /(cosh(beta) - cos(beta) - sigma*(sinh(beta) - sin(beta)));
            };
            scalar q = 0, u = 0;
            forAll(pp, i)
            {
                const scalar f = 1000.0*pp[i]*Sf[i].y()/0.015;
                q += f*phi(Cf[i].x());
                u += f;
            }
            reduce(q, sumOp<scalar>());
            reduce(u, sumOp<scalar>());
            Info<< "GENF " << m.time().value() << " " << q << " " << u << endl;
"""


def case_name(level: int, delta_t: float, tag: str = "") -> str:
    return f"ale_{level}x_dt{delta_t:g}{tag}"


def build_case(level: int, delta_t: float, duration: float, cores: int,
               tag: str = "", plate_bc: str | None = None,
               extra_functions: str = "") -> Path:
    source = WORK_ROOT / cfd3_case_name(level, delta_t)
    if not (source / "25").is_dir():
        driver.fail(f"{source} has no developed state at 25 s")
    case = WORK_ROOT / case_name(level, delta_t, tag)
    if case.exists():
        shutil.rmtree(case)
    for sub in ("constant", "system"):
        shutil.copytree(source / sub, case / sub, symlinks=True)
    shutil.copytree(source / "25", case / "25")
    # Mesh motion: the FSI3 fluid-mesh motion
    shutil.copyfile(TUTORIAL / "constant/fluid/dynamicMeshDict",
                    case / "constant/dynamicMeshDict")
    # Plate is a moving wall, as in FSI3
    u = case / "25/U"
    text = u.read_text()
    text, count = re.subn(r"(plate\s*\{\s*type\s+)fixedValue", r"\g<1>newMovingWallVelocity", text)
    if count != 1:
        driver.fail("Could not restore the moving-wall velocity on the plate")
    u.write_text(text)
    # Point velocity of the prescribed motion
    pm = (TUTORIAL / "0/fluid/pointMotionU").read_text()
    beta, sigma = MODES[MODE]
    code = BC_CODE % dict(A=AMPLITUDE, f=FREQUENCY, t0=START, ramp=RAMP,
                          xr=X_ROOT, xt=X_TIP, beta=repr(beta), sigma=repr(sigma))
    plate = plate_bc or ("plate\n    {\n        type            codedFixedValue;\n"
             "        value           uniform (0 0 0);\n"
             "        name            plateBendingMode%d;" % MODE + "\n"
             "        code\n        #{" + code + "        #};\n    }")
    pm, count = re.subn(r"plate\s*\{[^}]*\}", lambda m: plate, pm, count=1)
    if count != 1:
        driver.fail("Could not set the plate point motion")
    (case / "25/pointMotionU").write_text(pm)
    # Run control: restart from 25 s; whole-plate and plate-only forces
    control = case / "system/controlDict"
    driver.replace_entry(control, "startFrom", "startTime")
    driver.replace_entry(control, "startTime", f"{START:g}")
    driver.replace_entry(control, "endTime", f"{START + duration:.8g}")
    driver.replace_entry(control, "deltaT", f"{delta_t:.8g}")
    driver.replace_entry(control, "writeInterval", str(round(0.5 / delta_t)))
    template = (
        "    %(name)s\n    {\n"
        "        type                forces;\n"
        "        libs                ( \"libforces.so\" );\n"
        "        writeControl        timeStep;\n"
        "        writeInterval       1;\n"
        "        patches             %(patches)s;\n"
        "        \"pName|p\"           p;\n"
        "        \"UName|U\"           U;\n"
        "        \"rhoName|rho\"       rhoInf;\n"
        "        log                 false;\n"
        "        rhoInf              1000;\n"
        "        CofR                (0.5 0.1 0);\n"
        "    }\n")
    # Pressure force on the flag projected on the prescribed mode shape
    # (generalised force), printed to the log as "GENF t Q_mode Q_uniform"
    generalised = (
        "    generalisedForce\n    {\n"
        "        type            coded;\n"
        "        libs            ( utilityFunctionObjects );\n"
        "        name            generalisedForce;\n"
        "        codeExecute\n        #{\n"
        + GENERALISED_CODE % dict(xr=X_ROOT, xt=X_TIP, beta=repr(MODES[MODE][0]),
                                  sigma=repr(MODES[MODE][1]))
        + "        #};\n    }\n")
    (case / "system/functions").write_text(
        "functions\n{\n"
        + template % dict(name="forces", patches="(plate cylinder)")
        + template % dict(name="forcesPlate", patches="(plate)")
        + generalised + extra_functions + "}\n")
    driver.replace_entry(case / "system/decomposeParDict", "numberOfSubdomains", str(cores))
    return case


def run(case: Path, cores: int) -> None:
    def execute(command: list[str], log: str) -> None:
        with (case / log).open("w") as handle:
            result = subprocess.run(command, cwd=case, stdout=handle,
                                    stderr=subprocess.STDOUT, text=True)
        if result.returncode:
            driver.fail(f"{' '.join(command)} failed; see {case / log}")

    if cores > 1:
        driver.replace_entry(case / "system/decomposeParDict", "method", "scotch")
        execute(["decomposePar", "-force", "-latestTime"], "log.decomposePar")
        execute(["mpirun", "-np", str(cores), "solids4Foam", "-parallel"], "log.solids4Foam")
        execute(["reconstructPar", "-latestTime"], "log.reconstructPar")
    else:
        execute(["solids4Foam"], "log.solids4Foam")
    text = (case / "log.solids4Foam").read_text(errors="replace")
    if re.search(r"FOAM FATAL|FOAM aborting", text) or not re.search(
            r"^End\s*$", text, re.MULTILINE):
        driver.fail(f"{case.name} did not complete; see log.solids4Foam")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--level", type=int, required=True, choices=(1, 2, 4))
    parser.add_argument("--delta-t", type=float, required=True)
    parser.add_argument("--duration", type=float, default=8.0)
    parser.add_argument("--cores", type=int, default=1)
    parser.add_argument("--tag", default="")
    parser.add_argument("--setup-only", action="store_true")
    args = parser.parse_args()
    case = build_case(args.level, args.delta_t, args.duration, args.cores, args.tag)
    print(f"{case.name}: case written to {case}")
    if not args.setup_only:
        run(case, args.cores)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
