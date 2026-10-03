#!/usr/bin/env python3
"""Diagnostics for the temporal order of the womersleyTube verification study.

Reuses the case set-up, history readers and Fourier extraction of
womersley_tube_verification.py, and adds:

- longer runs, with the fluid sampled over every period after the first;
- the quantities of interest for every period and for every two-period
  window, so that the start-up transient and the choice of the analysis
  window can be examined;
- case variants that change only the initial data or the tolerances
  (see VARIANTS), never the numerical method of the production case;
- submission of each run as a one-core Slurm job (or a local run).

Usage (from any directory, with OpenFOAM and solids4foam loaded):

    womersley_temporal_investigation.py run  SPEC [SPEC ...] [--local]
    womersley_temporal_investigation.py analyse [SPEC ...]
    womersley_temporal_investigation.py orders

A SPEC is coupling_m<factor>_n<steps>_p<periods>[_<variant>...], for example
robin_m2_n100_p6 or robin_m2_n200_p12_exactOldBoundary. Results go to
verification/temporal/ (summaries, committed) and the runs to
verification/work/temporal/ (ignored by git).
"""

from __future__ import annotations

import argparse
import cmath
import copy
import csv
import gzip
import json
import math
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import womersley_tube_verification as wtv  # noqa: E402

RUN_ROOT = wtv.VERIFICATION / "work" / "temporal"
RESULTS = wtv.VERIFICATION / "temporal"
QOIS = ("flow_amp", "flow_phase", "wallMid_amp", "wallMid_phase", "speed",
        "attenuation", "speed_ux", "attenuation_ux", "profile", "wallQuarter_amp", "wallThreeQuarter_amp",
        "axialMid_amp", "axialMid_phase", "pressure")
SIGNED = ("flow_amp", "flow_phase", "wallMid_amp", "wallMid_phase", "speed",
          "attenuation", "speed_ux", "attenuation_ux", "wallQuarter_amp", "wallThreeQuarter_amp",
          "axialMid_amp", "axialMid_phase")


# ---------------------------------------------------------------------------
# Specs and variants
# ---------------------------------------------------------------------------

def parse_spec(text: str) -> dict:
    match = re.fullmatch(r"(robin|iqnils)_m(\d+)_n(\d+)_p(\d+)((?:_[A-Za-z0-9]+)*)",
                         text)
    if not match:
        wtv.fail(f"Bad spec {text!r}")
    variants = [v for v in match.group(5).split("_") if v]
    for v in variants:
        if v not in VARIANTS:
            wtv.fail(f"Unknown variant {v!r}; known: {sorted(VARIANTS)}")
    return {"name": text, "coupling": match.group(1),
            "factor": int(match.group(2)), "steps": int(match.group(3)),
            "periods": int(match.group(4)), "variants": variants}


def set_option(path: Path, key: str, value: str) -> None:
    """Replace every `key value;` entry (any indentation) in a dictionary."""
    text = path.read_text()
    text, found = re.subn(rf"^(\s*{re.escape(key)}\s+)[^;\n]+;",
                          rf"\g<1>{value};", text, flags=re.MULTILINE)
    if not found:
        wtv.fail(f"No '{key}' entry in {path}")
    path.write_text(text)


def set_petsc_option(path: Path, key: str, value: str) -> None:
    text = path.read_text()
    text, found = re.subn(rf'^(\s*{re.escape(key)}\s+)"[^"]*";',
                          rf'\g<1>"{value}";', text, flags=re.MULTILINE)
    if found != 1:
        wtv.fail(f"Expected one '{key}' PETSc option in {path}")
    path.write_text(text)


# The exact displacement on a boundary patch, at the time level given by the
# womersleyTimeLevel entry of the field file, for the boundary values of the
# initial and old-time solid displacement fields
PATCH_CODE = r"""
womersleyPatchDisplacement
#{
    const dictionary& top = dict.topDict();
    const fvMesh& mesh =
        refCast<const fvMesh>(static_cast<const IOdictionary&>(top).db());
    const womersley::mode& m = womersley::lookup(mesh.time());
    const scalar level =
        top.found("womersleyTimeLevel")
      ? readScalar(top.lookup("womersleyTimeLevel"))
      : 0;
    const scalar t = mesh.time().value() - level*mesh.time().deltaTValue();
    const label patchI = mesh.boundaryMesh().findPatchID(dict.dictName());
    if (patchI < 0)
    {
        FatalErrorInFunction << "No patch " << dict.dictName()
            << exit(FatalError);
    }
    const vectorField& Cf = mesh.boundaryMesh()[patchI].faceCentres();
    vectorField value(Cf.size());
    forAll(Cf, faceI)
    {
        value[faceI] = m.displacement(Cf[faceI], t);
    }
    os  << word("nonuniform") << token::SPACE << value;
#};
"""

EXACT_VALUE = """value #codeStream
        {
            codeInclude     $womersleyCodeInclude;
            codeOptions     $womersleyCodeOptions;
            codeLibs        $womersleyCodeLibs;
            code            $womersleyPatchDisplacement;
        };"""


def exact_old_boundary(case: Path) -> None:
    """Initial-data variant: the boundary values of D and of the old-time
    levels D_0, D_0_0 and D_0_0_0 on the inlet, outlet, interface and outer
    patches are the exact solution at their time level, instead of the
    placeholder zero (old-time fields are never re-evaluated, and the
    boundary values enter the interface acceleration of the Robin
    condition). The cell values, and every other field, are unchanged."""
    code = case / "system" / "womersleyCode"
    text = code.read_text()
    marker = "// ************"
    text = text.replace(marker, PATCH_CODE + "\n" + marker, 1) \
        if marker in text else text + PATCH_CODE
    code.write_text(text)
    for name in ("D", "D_0", "D_0_0", "D_0_0_0"):
        path = case / "0" / "solid" / name
        text = path.read_text()
        # Split the regex entry so that each patch knows its name
        head, sep, rest = text.partition('"inlet|outlet"')
        if not sep:
            wtv.fail(f"No inlet|outlet entry in {path}")
        body_end = rest.index("}") + 1
        entry = rest[:body_end]
        rest = rest[body_end:]
        text = head + "inlet" + entry + "\n    outlet" + entry + rest
        # Old-time levels: the inlet/outlet are calculated; D: codedFixedValue,
        # whose value is re-evaluated anyway
        for patch in ("inlet", "outlet", "interface", "outer"):
            pattern = (rf"(\n    {patch}\n    \{{[^}}]*?)"
                       r"value\s+uniform \(0 0 0\);")
            text, found = re.subn(pattern, rf"\g<1>{EXACT_VALUE}", text,
                                  flags=re.DOTALL)
            if found != 1:
                wtv.fail(f"Could not set the {patch} value in {path}")
        path.write_text(text)


def tight_tolerances(case: Path) -> None:
    """Coupling and solid solver tolerances tightened by about 100 times."""
    for name in ("fsiProperties.robin", "fsiProperties.iqnils"):
        path = case / "constant" / name
        set_option(path, "outerCorrTolerance", "1e-8")
        set_option(path, "nOuterCorr", "1000")
        if name.endswith("robin"):
            set_option(path, "robinPressureTolerance", "1e-7")
            set_option(path, "robinFluxTolerance", "5e-5")
    solid = case / "system" / "solid" / "fvSolution"
    set_petsc_option(solid, "snes_rtol", "1e-10")
    set_petsc_option(solid, "snes_atol", "1e-11")
    set_petsc_option(solid, "snes_stol", "1e-8")
    fluid = case / "system" / "fluid" / "fvSolution"
    set_option(fluid, "tolerance", "1e-14")


def pimple4(case: Path) -> None:
    """Four PIMPLE outer correctors per coupling iteration instead of two."""
    set_option(case / "system" / "fluid" / "fvSolution", "nOuterCorrectors",
               "4")


FLUID_ONLY_ALLRUN = """#!/bin/bash
# Fluid-only diagnostic: the fluid region of womersleyTube with the exact
# wall motion imposed on the moving mesh (no solid, no coupling)
. $WM_PROJECT_DIR/bin/tools/RunFunctions
source solids4FoamScripts.sh
ln -nsf U.dirichletNeumann 0/U
ln -nsf p.dirichletNeumann 0/p
runApplication blockMesh
runApplication topoSet
runApplication solids4Foam
mkdir -p postProcessing/profileSamples postProcessing/axialSamples
ln -nsf . postProcessing/profileSamples/fluid
ln -nsf . postProcessing/axialSamples/fluid
"""

# Exact displacement increment over the time-step at the interface points,
# divided by the time-step: the fluid mesh then follows the exact wall from
# its undeformed start, as it follows the solid increments in the coupled case
MESH_VELOCITY = """
    interface
    {
        type            codedFixedValue;
        value           uniform (0 0 0);
        name            womersleyMeshVelocity;
        codeInclude     $womersleyCodeInclude;
        code
        #{
            const womersley::mode& m = womersley::lookup(this->db().time());
            const scalar t = this->db().time().value();
            const scalar dt = this->db().time().deltaTValue();
            const pointField& x = this->patch().localPoints();
            vectorField v(x.size());
            forAll(x, i)
            {
                // Undeformed position, from x = X + d(X, t - dt) - d(X, 0)
                vector X = x[i];
                for (int it = 0; it < 3; ++it)
                {
                    X = x[i] - (m.displacement(X, t - dt) - m.displacement(X, 0));
                }
                v[i] = (m.displacement(X, t) - m.displacement(X, t - dt))/dt;
            }
            operator==(v);
        #};
    }
"""


def fluid_only(case: Path) -> None:
    """Fluid sub-problem alone: the fluid region with the exact wall motion
    imposed, newMovingWallVelocity and zero-gradient pressure on the wall
    (the Dirichlet-Neumann fluid conditions), and nothing else changed."""
    for sub in ("0", "constant", "system"):
        region = case / sub / "fluid"
        for item in region.iterdir():
            target = case / sub / item.name
            if target.exists() or target.is_symlink():
                if target.is_dir() and not target.is_symlink():
                    shutil.rmtree(target)
                else:
                    target.unlink()
            item.rename(target)
        region.rmdir()
        solid = case / sub / "solid"
        if solid.exists():
            shutil.rmtree(solid)
    for name in ("U", "p"):
        link = case / "0" / name
        if link.is_symlink() or link.exists():
            link.unlink()
    physics = case / "constant" / "physicsProperties"
    physics.write_text(re.sub(r"^type\s+\w+;", "type    fluid;",
                              physics.read_text(), flags=re.MULTILINE))
    motion = case / "0" / "pointMotionU"
    text = motion.read_text()
    text = text.replace("dimensions",
                        '#include "$FOAM_CASE/system/womersleyCode"\n\n'
                        "dimensions", 1)
    text, found = re.subn(r"\n    interface\n    \{[^}]*\}\n", MESH_VELOCITY,
                          text)
    if found != 1:
        wtv.fail("Could not set the interface mesh velocity")
    motion.write_text(text)
    control = case / "system" / "controlDict"
    text = control.read_text()
    text = re.sub(r"\n    wall(Quarter|Mid|ThreeQuarter)\n    \{[^}]*\}\n", "\n",
                  text)
    text = re.sub(r"\n\s*region\s+fluid;", "", text)
    control.write_text(text)
    allrun = case / "Allrun"
    allrun.write_text(FLUID_ONLY_ALLRUN)
    allrun.chmod(0o755)


def pimple10(case: Path) -> None:
    """Ten PIMPLE outer correctors per time-step (or coupling iteration)."""
    path = case / "system" / "fluid" / "fvSolution"
    if not path.exists():
        path = case / "system" / "fvSolution"
    set_option(path, "nOuterCorrectors", "10")


WALL_VELOCITY = """
    interface
    {
        type            codedFixedValue;
        value           uniform (0 0 0);
        name            womersleyWallVelocity;
        codeInclude     $womersleyCodeInclude;
        code
        #{
            const womersley::mode& m = womersley::lookup(this->db().time());
            const scalar t = this->db().time().value();
            const vectorField& Cf = this->patch().Cf();
            vectorField value(Cf.size());
            forAll(Cf, faceI)
            {
                value[faceI] = m.velocity(Cf[faceI], t);
            }
            operator==(value);
        #};
    }
"""


def transpiration(case: Path) -> None:
    """Fluid-only variant on a static mesh: the exact fluid (= wall)
    velocity is imposed at the undeformed wall, so there is no mesh motion,
    no ALE flux and no moving-wall velocity condition."""
    if not (case / "constant" / "dynamicMeshDict").exists():
        wtv.fail("transpiration must follow fluidOnly")
    (case / "constant" / "dynamicMeshDict").write_text(
        (case / "constant" / "dynamicMeshDict").read_text().split(
            "dynamicFvMesh")[0] + "dynamicFvMesh   staticFvMesh;\n")
    (case / "0" / "pointMotionU").unlink()
    path = case / "0" / "U.dirichletNeumann"
    text, found = re.subn(r"\n    interface\n    \{[^}]*\}\n", WALL_VELOCITY,
                          path.read_text())
    if found != 1:
        wtv.fail("Could not set the interface velocity")
    path.write_text(text)


def fluid_dict(case: Path, name: str) -> Path:
    path = case / "system" / "fluid" / name
    return path if path.exists() else case / "system" / name


def ddt_coeff(value: str):
    def apply(case: Path) -> None:
        """Diagnostic only: a constant Rhie-Chow ddtCorr coupling
        coefficient (backward <value>) instead of OpenFOAM's flux-ratio
        limiter, 1 - min(|phi - U_f.S|/|phi|, 1)."""
        path = fluid_dict(case, "fvSchemes")
        text, found = re.subn(r"(ddtSchemes\s*\{\s*default\s+)backward;",
                              rf"\g<1>backward {value};", path.read_text())
        if found != 1:
            wtv.fail(f"Could not set the ddt scheme in {path}")
        path.write_text(text)
    return apply


def no_ddt_corr(case: Path) -> None:
    """Diagnostic only: no Rhie-Chow ddtCorr term (PIMPLE ddtCorr false)."""
    path = fluid_dict(case, "fvSolution")
    text, found = re.subn(r"(PIMPLE\s*\{)", r"\g<1>\n    ddtCorr             false;",
                          path.read_text())
    if found != 1:
        wtv.fail(f"No PIMPLE dictionary in {path}")
    path.write_text(text)


def rigid(case: Path) -> None:
    """Diagnostic only (self-convergence, no exact solution): the static-
    mesh fluid problem with a no-slip wall instead of the exact wall
    velocity; the tube ends keep the exact travelling-wave data."""
    path = case / "0" / "U.dirichletNeumann"
    text, found = re.subn(r"\n    interface\n    \{\n.*?\n    \}\n",
                          "\n    interface\n    {\n        type            noSlip;\n    }\n",
                          path.read_text(), flags=re.DOTALL)
    if found != 1:
        wtv.fail("rigid must follow transpiration")
    path.write_text(text)


def ends_zero_gradient(case: Path) -> None:
    """Diagnostic only (self-convergence): zero-gradient velocity at the
    tube ends instead of the exact normal gradient (codedMixed)."""
    path = fluid_dict(case, "fvSchemes").parent.parent / "0"
    path = path / "fluid" / "U.dirichletNeumann" \
        if (path / "fluid").is_dir() else path / "U.dirichletNeumann"
    text, found = re.subn(r'\n    "inlet\|outlet"\n    \{\n.*?\n    \}\n',
                          '\n    "inlet|outlet"\n    {\n        type            zeroGradient;\n    }\n',
                          path.read_text(), flags=re.DOTALL)
    if found != 1:
        wtv.fail(f"Could not set the end velocity condition in {path}")
    path.write_text(text)


# The same exact normal velocity gradient at the tube ends, imposed through a
# fixedGradient condition (assignable, so constrainHbyA does not replace HbyA
# by the boundary velocity) instead of codedMixed with valueFraction 0
# (which OpenFOAM treats as fixing the value). A coded function object sets
# the gradient for the coming time level at the start of the run and at the
# end of every time-step
SET_END_GRADIENT = """
            const Foam::fvMesh& mesh = this->mesh();
            const Foam::Time& runTime = mesh.time();
            const womersley::mode& m = womersley::lookup(runTime);
            const Foam::scalar t = runTime.value() + runTime.deltaTValue();
            Foam::volVectorField& U = const_cast<Foam::volVectorField&>
            (
                mesh.lookupObject<Foam::volVectorField>("U")
            );
            for (const Foam::word patch : {"inlet", "outlet"})
            {
                const Foam::label patchI =
                    mesh.boundaryMesh().findPatchID(patch);
                auto& pf = Foam::refCast<Foam::fixedGradientFvPatchVectorField>
                (
                    U.boundaryFieldRef()[patchI]
                );
                const Foam::vectorField& Cf = mesh.boundary()[patchI].Cf();
                const Foam::vectorField nf(mesh.boundary()[patchI].nf());
                forAll(Cf, faceI)
                {
                    pf.gradient()[faceI] = m.velocity(Cf[faceI], t, &nf[faceI]);
                }
            }
"""

END_GRADIENT_FO = """
    womersleyEndGradient
    {
        type            coded;
        libs            (utilityFunctionObjects);
        name            womersleyEndGradient;
        REGION
        codeInclude     $womersleyCodeInclude;
        codeOptions
        #{
            -I$(LIB_SRC)/finiteVolume/lnInclude \\
            -I$(LIB_SRC)/meshTools/lnInclude \\
            -include fixedGradientFvPatchFields.H
        #};
        codeLibs
        #{
            -lfiniteVolume \\
            -lmeshTools
        #};
        codeRead
        #{""" + SET_END_GRADIENT + """        #};
        codeExecute
        #{""" + SET_END_GRADIENT + """        #};
    }
"""


def ends_fixed_gradient(case: Path) -> None:
    """The exact end velocity gradient through fixedGradient (assignable)
    instead of codedMixed (non-assignable). Same continuum problem."""
    zero = case / "0" / "fluid"
    zero = zero if zero.is_dir() else case / "0"
    for name in ("U.dirichletNeumann", "U.robin"):
        path = zero / name
        if not path.exists():
            continue
        text, found = re.subn(
            r'\n    "inlet\|outlet"\n    \{\n.*?\n    \}\n',
            '\n    "inlet|outlet"\n    {\n        type            fixedGradient;\n'
            '        gradient        uniform (0 0 0);\n    }\n',
            path.read_text(), flags=re.DOTALL)
        if found != 1:
            wtv.fail(f"Could not set the end velocity condition in {path}")
        path.write_text(text)
    control = case / "system" / "controlDict"
    text = control.read_text()
    region = "region          fluid;" if (case / "0" / "fluid").is_dir() else ""
    text = text.replace("functions\n{", "functions\n{\n"
                        + END_GRADIENT_FO.replace("REGION", region), 1)
    if "womersleyEndGradient" not in text:
        wtv.fail("Could not add the end-gradient function object")
    text = text.replace("adjustTimeStep", '#include "$FOAM_CASE/system/womersleyCode"\n\nadjustTimeStep', 1)
    control.write_text(text)


def least_squares_s4f(case: Path) -> None:
    """Fluid gradients by leastSquaresS4f (full boundary delta vectors)
    instead of OpenFOAM's leastSquares."""
    path = fluid_dict(case, "fvSchemes")
    text, found = re.subn(r"(gradSchemes\s*\{\s*default\s+)leastSquares;",
                          r"\g<1>leastSquaresS4f;", path.read_text())
    if found != 1:
        wtv.fail(f"Could not set the gradient scheme in {path}")
    path.write_text(text)


END_VELOCITY = """
    "inlet|outlet"
    {
        type            codedFixedValue;
        value           uniform (0 0 0);
        name            womersleyEndVelocity;
        codeInclude     $womersleyCodeInclude;
        code
        #{
            const womersley::mode& m = womersley::lookup(this->db().time());
            const scalar t = this->db().time().value();
            const vectorField& Cf = this->patch().Cf();
            vectorField value(Cf.size());
            forAll(Cf, faceI)
            {
                value[faceI] = m.velocity(Cf[faceI], t);
            }
            operator==(value);
        #};
    }
"""


def ends_dirichlet_velocity(case: Path) -> None:
    """Fluid-only diagnostic: the exact velocity at the tube ends with
    fixedFluxPressure, the pair for which OpenFOAM makes the end flux equal
    to the boundary velocity flux; the pressure level is then set by a
    reference cell, so only gauge-free quantities are meaningful."""
    zero = case / "0"
    path = zero / "U.dirichletNeumann"
    text, found = re.subn(r'\n    "inlet\|outlet"\n    \{\n.*?\n    \}\n',
                          END_VELOCITY, path.read_text(), flags=re.DOTALL)
    if found != 1:
        wtv.fail("Could not set the end velocity")
    path.write_text(text)
    path = zero / "p.dirichletNeumann"
    text, found = re.subn(r'\n    "inlet\|outlet"\n    \{\n.*?\n    \}\n',
                          '\n    "inlet|outlet"\n    {\n        type            fixedFluxPressure;\n'
                          '        value           uniform 0;\n    }\n',
                          path.read_text(), flags=re.DOTALL)
    if found != 1:
        wtv.fail("Could not set the end pressure")
    path.write_text(text)
    fv = fluid_dict(case, "fvSolution")
    text, found = re.subn(r"(PIMPLE\s*\{)", r"\g<1>\n    pRefCell            0;\n    pRefValue           0;",
                          fv.read_text())
    fv.write_text(text)


def flux_consistent(case: Path) -> None:
    """pimpleFluid fluxConsistentPatches on the tube ends: the same end
    conditions, with the end flux made equal to the boundary velocity flux
    (needs a solids4foam build with that option)."""
    path = fluid_dict(case, "fvSolution")
    text, found = re.subn(r"(PIMPLE\s*\{)",
                          r"\g<1>\n    fluxConsistentPatches (inlet outlet);",
                          path.read_text())
    if found != 1:
        wtv.fail(f"No PIMPLE dictionary in {path}")
    path.write_text(text)


def no_change(case: Path) -> None:
    """A tag only: of2412 runs the same case with the OpenFOAM-v2412 build."""


VARIANTS = {
    "exactOldBoundary": exact_old_boundary,
    "of2412": no_change,
    "fluidOnly": fluid_only,
    "transpiration": transpiration,
    "pimple10": pimple10,
    "ddtCoeff1": ddt_coeff("1"),
    "ddtCoeff0": ddt_coeff("0"),
    "noDdtCorr": no_ddt_corr,
    "rigid": rigid,
    "endsZeroGrad": ends_zero_gradient,
    "endsFixedGradient": ends_fixed_gradient,
    "lsS4f": least_squares_s4f,
    "endsDirichletU": ends_dirichlet_velocity,
    "fluxConsistent": flux_consistent,
    # The same option, second implementation (HbyA_b = U_b + rAtU grad(p)_P)
    "fluxConsistent2": flux_consistent,
    "tight": tight_tolerances,
    "pimple4": pimple4,
}


def reference_for(spec: dict) -> dict:
    reference = json.loads(wtv.REFERENCE_FILE.read_text())
    reference = copy.deepcopy(reference)
    reference["study"]["periods"] = spec["periods"]
    # The fluid is sampled over every period after the first
    reference["study"]["analysisPeriods"] = spec["periods"] - 1
    return reference


def configure(spec: dict) -> Path:
    case = RUN_ROOT / spec["name"]
    if case.exists():
        shutil.rmtree(case)
    case.parent.mkdir(parents=True, exist_ok=True)
    wtv.configure(case, {"factor": spec["factor"], "steps": spec["steps"]},
                  reference_for(spec))
    for variant in spec["variants"]:
        VARIANTS[variant](case)
    (case / "temporal_spec.json").write_text(json.dumps(
        {**spec, "build": wtv.build_signature()}, indent=2, sort_keys=True))
    return case


# ---------------------------------------------------------------------------
# Running
# ---------------------------------------------------------------------------

def run(specs: list[dict], local: bool, env_script: str, version: str,
        hours: int, nprocs: int = 1) -> None:
    for spec in specs:
        case = configure(spec)
        if "of2412" in spec["variants"]:
            version = "v2412"
        mode = ""
        if nprocs > 1:
            # Axial slabs, with no processor boundary near a sampling
            # station: the flow-rate faceZone (L/2), the profile set
            # (L/2 + 0.03) and the wall points (L/4, L/2, 3L/4); a faceZone
            # on a processor boundary gives a wrong flow rate
            length = 15.0
            cuts = [length * i / nprocs for i in range(1, nprocs)]
            stations = (3.75, 7.5, 7.53, 11.25)
            if any(abs(c - x) < 0.5 for c in cuts for x in stations):
                wtv.fail(f"{nprocs} axial slabs put a processor boundary "
                         "near a sampling station; choose another --np")
            for path in case.glob("system/**/decomposeParDict"):
                set_option(path, "numberOfSubdomains", str(nprocs))
                text = path.read_text()
                text = re.sub(r"^method\s+\w+;\n", "", text,
                              flags=re.MULTILINE)
                text = text.replace(
                    f"numberOfSubdomains {nprocs};",
                    f"numberOfSubdomains {nprocs};\n\nmethod          simple;"
                    f"\n\nsimpleCoeffs\n{{\n    n               "
                    f"({nprocs} 1 1);\n}}", 1)
                path.write_text(text)
            mode = " parallel"
        command = (f"source {env_script} {version} && cd {case} && "
                   f"./Allrun {spec['coupling']}{mode} > log.Allverify 2>&1")
        if local:
            subprocess.run(["bash", "-c", command], check=False)
            print(f"{spec['name']}: {wtv.solver_log_problem(case) or 'done'}")
        else:
            job = subprocess.run(
                ["sbatch", "--parsable", f"--job-name=wt_{spec['name']}",
                 "--partition=main", f"--ntasks={nprocs}", "--cpus-per-task=1",
                 f"--time={hours}:00:00", f"--output={case}/log.slurm",
                 f"--chdir={case}", "--wrap", f"bash -c '{command}'"],
                check=True, capture_output=True, text=True).stdout.strip()
            print(f"{spec['name']}: Slurm job {job}")


# ---------------------------------------------------------------------------
# Analysis
# ---------------------------------------------------------------------------

def coefficient(times, values, omega, period, steps, first, count):
    """Fourier coefficient over periods first..first+count-1 (1-based)."""
    t0 = (first - 1) * period
    t1 = t0 + count * period
    margin = 0.25 * period / steps
    selected = [(t, v) for t, v in zip(times, values)
                if t0 + margin < t < t1 + margin]
    if len(selected) != count * steps:
        wtv.fail(f"{len(selected)} samples in periods {first}.."
                 f"{first + count - 1}, expected {count * steps}")
    return wtv.fourier([t for t, _ in selected], [v for _, v in selected],
                       omega)


def relative(c: complex, exact: complex) -> tuple[float, float]:
    return abs(c) / abs(exact) - 1, cmath.phase(c / exact)


def fluid_series(case: Path, exact: wtv.Exact, periods: int, expected: dict):
    """Fluid samples at every sampled time: {set: (times, rows)}."""
    out = {}
    for set_name in ("profile", "axial"):
        root = case / wtv.SAMPLES.format(set_name)
        times = sorted((float(d.name), d) for d in root.iterdir()
                       if re.fullmatch(r"[0-9.eE+-]+", d.name))
        rows = [wtv.read_rows(d / f"{set_name}_p_U.xy", 5) for _, d in times]
        # In parallel, a set point on a processor boundary is written by
        # both processors, with coordinates that differ by round-off: merge
        # points closer than 1e-2 m (the set spacings are 0.05 m or more;
        # the copies differ by up to about 1e-4 m on the moving mesh)
        def unique(rs):
            out = []
            for r in sorted(rs, key=lambda r: r[0]):
                if out and abs(r[0] - out[-1][0][0]) < 1e-2:
                    out[-1].append(r)
                else:
                    out.append([r])
            return [[sum(c) / len(group) for c in zip(*group)]
                    for group in out]
        rows = [unique(rs) for rs in rows]
        if any(len(r) != expected[set_name] for r in rows):
            wtv.fail(f"Unexpected {set_name} sample count in {root}")
        out[set_name] = ([t for t, _ in times], rows)
    return out


def fluid_window(series, exact: wtv.Exact, mesh: dict, first: int,
                 count: int) -> dict:
    """Profile, wave number and pressure QoIs over a window of periods."""
    per = wtv.SAMPLES_PER_PERIOD
    t0 = (first - 1) * exact.period
    wanted = wtv.period_times(t0, count * exact.period, count * per)
    result = {}
    data = {}
    for set_name, (times, rows) in series.items():
        index = []
        for t in wanted:
            match = [i for i, tt in enumerate(times)
                     if abs(tt - t) < 1e-6 * exact.period]
            if len(match) != 1:
                wtv.fail(f"Missing fluid sample at t = {t}")
            index.append(match[0])
        n_points = len(rows[index[0]])
        data[set_name] = {
            "coord": [row[0] for row in rows[index[0]]],
            "p": [wtv.fourier(wanted, [rows[i][j][1] for i in index],
                              exact.omega) for j in range(n_points)],
            "ux": [wtv.fourier(wanted, [rows[i][j][2] for i in index],
                               exact.omega) for j in range(n_points)],
        }
    prof = data["profile"]
    x_cells = 0.5 * exact.L + 0.5 * exact.L / mesh["axial"]
    dr = exact.R / mesh["fluidRadial"]
    chord = math.cos(math.radians(0.5))
    radii = [chord * 2 / 3 * ((j + 1)**3 - j**3) * dr / ((j + 1)**2 - j**2)
             for j in range(len(prof["coord"]))]
    u_exact = [exact.velocity(r, x_cells) for r in radii]
    scale = max(abs(u) for u in u_exact)
    result["profile"] = max(abs(u - ue) for u, ue in
                            zip(prof["ux"], u_exact)) / scale
    ax = data["axial"]
    k = wtv.fit_wave_number(ax["coord"], ax["p"])
    result["k_real"], result["k_imag"] = k.real, k.imag
    result["speed"] = exact.k.real / k.real - 1
    result["attenuation"] = k.imag / exact.k.imag - 1
    p_exact = [exact.pressure(0.5 * exact.R, x) for x in ax["coord"]]
    result["pressure"] = max(abs(exact.rho * p - pe) for p, pe in
                             zip(ax["p"], p_exact)) / max(abs(pe) for pe in
                                                          p_exact)
    # Alternative wave-number extraction (H8): a fit of the complex
    # pressure itself, p = C exp(-i k x), by Gauss-Newton from the log fit,
    # which weights the points by the pressure rather than equally in log
    result["k_alt_real"], result["k_alt_imag"] = \
        complex_exponential_fit(ax["coord"], ax["p"], k)
    # Gauge-free wave number from the axial velocity at r = R/2
    ku = wtv.fit_wave_number(ax["coord"], ax["ux"])
    result["speed_ux"] = exact.k.real / ku.real - 1
    result["attenuation_ux"] = ku.imag / exact.k.imag - 1
    result["speed_alt"] = exact.k.real / result["k_alt_real"] - 1
    result["attenuation_alt"] = result["k_alt_imag"] / exact.k.imag - 1
    return result


def complex_exponential_fit(x, p, k0):
    """Least-squares fit of p(x) = C exp(-i k x), complex C and k."""
    k = k0
    for _ in range(50):
        e = [cmath.exp(-1j * k * xi) for xi in x]
        c = sum(pi * ei.conjugate() for pi, ei in zip(p, e)) \
            / sum(abs(ei) ** 2 for ei in e)
        # Gauss-Newton on k with C eliminated (variable projection, simple)
        r = [pi - c * ei for pi, ei in zip(p, e)]
        j = [-1j * xi * c * ei for xi, ei in zip(x, e)]
        dk = sum(ji.conjugate() * ri for ji, ri in zip(j, r)) \
            / sum(abs(ji) ** 2 for ji in j)
        k += dk
        if abs(dk) < 1e-15 * abs(k):
            break
    return k.real, k.imag


def analyse_case(name: str) -> dict:
    spec = parse_spec(name)
    case = RUN_ROOT / name
    problem = wtv.solver_log_problem(case)
    if problem:
        wtv.fail(f"{name} {problem}")
    reference = reference_for(spec)
    exact = wtv.Exact(reference)
    P, n = spec["periods"], spec["steps"]
    mesh = wtv.mesh_divisions(spec["factor"], reference)
    out = {"spec": spec, "per_period": {}, "windows": {}, "transient": {}}

    histories = {}
    fluid_only_run = "fluidOnly" in spec["variants"]
    for station, fraction in ({} if fluid_only_run else wtv.STATIONS).items():
        times, radial, axial = wtv.wall_history(case, station, exact, n, P)
        histories[station] = (times, radial, axial)
    times_q, flux = wtv.flow_rate_history(case, exact, n, P)
    series = fluid_series(case, exact, P, {"profile": mesh["fluidRadial"],
                                           "axial": wtv.AXIAL_POINTS})

    def wall_and_flow(first: int, count: int) -> dict:
        row = {}
        for station, fraction in ({} if fluid_only_run else
                                  wtv.STATIONS).items():
            times, radial, axial = histories[station]
            c = coefficient(times, radial, exact.omega, exact.period, n,
                            first, count)
            row[f"wall{station}_amp"], row[f"wall{station}_phase"] = \
                relative(c, exact.eta * exact.wave(fraction * exact.L))
            c = coefficient(times, axial, exact.omega, exact.period, n,
                            first, count)
            row[f"axial{station}_amp"], row[f"axial{station}_phase"] = \
                relative(c, exact.xi * exact.wave(fraction * exact.L))
        c = 360 * coefficient(times_q, flux, exact.omega, exact.period, n,
                              first, count)
        row["flow_amp"], row["flow_phase"] = \
            relative(c, exact.Q * exact.wave(0.5 * exact.L))
        return row

    for p in range(1, P + 1):
        row = wall_and_flow(p, 1)
        if p >= 2:
            row.update(fluid_window(series, exact, mesh, p, 1))
        out["per_period"][p] = row
    for p in range(2, P):
        row = wall_and_flow(p, 2)
        row.update(fluid_window(series, exact, mesh, p, 2))
        out["windows"][f"{p}-{p + 1}"] = row

    if fluid_only_run:
        return finish(case, name, out, times_q, None, None, flux)

    # Transient: the wall displacement at L/2 minus the exact periodic
    # solution, relative to the exact amplitude: its largest value in each
    # period, and its largest value and dominant frequency after removing
    # the fundamental fitted to the last two periods
    times, radial, _ = histories["Mid"]
    eta = exact.eta * exact.wave(0.5 * exact.L)
    err = [(v - (eta * cmath.exp(1j * exact.omega * t)).real) / abs(eta)
           for t, v in zip(times, radial)]
    out["transient"]["max_error_per_period"] = {
        p: max(abs(e) for t, e in zip(times, err)
               if (p - 1) * exact.period < t <= p * exact.period + 1e-9)
        for p in range(1, P + 1)}
    out["transient"]["first_steps_error"] = err[:10]
    c_last = coefficient(times, radial, exact.omega, exact.period, n, P - 1, 2)
    resid = [v - (c_last * cmath.exp(1j * exact.omega * t)).real
             for t, v in zip(times, radial)]
    out["transient"]["max_residual_per_period"] = {
        p: max(abs(r) for t, r in zip(times, resid)
               if (p - 1) * exact.period < t <= p * exact.period + 1e-9)
        / abs(eta) for p in range(1, P + 1)}
    out["transient"]["residual_spectrum"] = dominant_frequencies(
        times, resid, exact.period, n, P)
    _, _, axial = histories["Mid"]
    return finish(case, name, out, times, radial, axial, flux)


def finish(case, name, out, times, radial, axial, flux):
    log = (case / "log.solids4Foam").read_text(errors="replace")
    build = re.search(r"^Build\s*:\s*(.*)$", log, re.MULTILINE)
    out["openfoam_build"] = build.group(1).strip() if build else ""
    out["iterations"] = wtv.coupling_iterations(case)
    out["execution_time"] = wtv.execution_time(case)

    RESULTS.mkdir(parents=True, exist_ok=True)
    (RESULTS / "runs").mkdir(exist_ok=True)
    (RESULTS / "runs" / f"{name}.json").write_text(
        json.dumps(out, indent=1, sort_keys=True))
    # Compact histories (wall displacement at L/2 and flow rate), for
    # re-analysis without the run directories
    with gzip.open(RESULTS / "runs" / f"{name}_history.csv.gz", "wt") as f:
        f.write("t,wallMid_radial,wallMid_axial,flowRate_wedge\n")
        if radial is None:
            radial = axial = [math.nan] * len(times)
        for t, r, a, q in zip(times, radial, axial, flux):
            f.write(f"{t:.10g},{r:.12e},{a:.12e},{q:.12e}\n")
    return out


def dominant_frequencies(times, values, period, steps, periods):
    """Amplitudes of the residual at multiples of 1/(2 T) up to 5/T over the
    first two periods, to identify free modes excited at the start."""
    count = 2
    sel = [(t, v) for t, v in zip(times, values) if t <= count * period + 1e-9]
    out = {}
    for m in range(1, 21):
        w = 2 * math.pi * m / (count * period)
        c = 2 / len(sel) * sum(v * cmath.exp(-1j * w * t) for t, v in sel)
        out[f"{m / count:g}"] = abs(c)
    return out


# ---------------------------------------------------------------------------
# Orders
# ---------------------------------------------------------------------------

def order3(values):
    d1, d2 = values[0] - values[1], values[1] - values[2]
    if d1 == 0 or d2 == 0 or d1 * d2 < 0:
        return math.nan
    return math.log(abs(d1 / d2)) / math.log(2)


def orders() -> None:
    runs = {}
    for path in sorted((RESULTS / "runs").glob("*.json")):
        data = json.loads(path.read_text())
        runs[data["spec"]["name"]] = data
    groups: dict[tuple, list] = {}
    for name, data in runs.items():
        s = data["spec"]
        key = (s["coupling"], s["factor"], s["periods"],
               "_".join(s["variants"]))
        groups.setdefault(key, []).append(data)
    rows = []
    for key, members in sorted(groups.items()):
        members.sort(key=lambda d: d["spec"]["steps"])
        steps = [d["spec"]["steps"] for d in members]
        windows = sorted(members[0]["windows"])
        for window in windows:
            for qoi in QOIS:
                vals = [d["windows"].get(window, {}).get(qoi) for d in members]
                if any(v is None for v in vals):
                    continue
                row = {"coupling": key[0], "factor": key[1],
                       "periods": key[2], "variant": key[3] or "base",
                       "window": window, "qoi": qoi}
                for s, v in zip(steps, vals):
                    row[f"n{s}"] = v
                for i in range(len(vals) - 1):
                    if qoi in SIGNED:
                        row[f"d{steps[i]}-{steps[i + 1]}"] = \
                            vals[i] - vals[i + 1]
                for i in range(len(vals) - 2):
                    if qoi in SIGNED:
                        row[f"order{steps[i]}-{steps[i + 2]}"] = \
                            order3(vals[i:i + 3])
                    a, b = abs(vals[i + 1]), abs(vals[i + 2])
                    row[f"errorOrder{steps[i + 1]}-{steps[i + 2]}"] = \
                        math.log(a / b) / math.log(2) if a and b else math.nan
                rows.append(row)
    fields = []
    for row in rows:
        for k in row:
            if k not in fields:
                fields.append(k)
    with (RESULTS / "temporal_orders.csv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)
    print(f"Wrote {RESULTS / 'temporal_orders.csv'} ({len(rows)} rows)")


def summary() -> None:
    """One row per run and analysis window: the set-up and every QoI error
    against the exact solution."""
    rows = []
    period = 50.0
    for path in sorted((RESULTS / "runs").glob("*.json")):
        data = json.loads(path.read_text())
        spec = data["spec"]
        base = {
            "name": spec["name"], "coupling": spec["coupling"],
            "mesh_factor": spec["factor"], "steps_per_period": spec["steps"],
            "dt": period / spec["steps"], "periods": spec["periods"],
            "variant": "_".join(spec["variants"]) or "base",
            "start": "exact initial fields"
            + (", exact old-time boundary values"
               if "exactOldBoundary" in spec["variants"] else ""),
            "tolerances": "tight (x0.01)" if "tight" in spec["variants"]
            else "tutorial",
            "openfoam_build": data.get("openfoam_build", ""),
            "mean_coupling_iterations": data["iterations"][0]
            if data.get("iterations") else "",
            "execution_time_s": data.get("execution_time", ""),
        }
        for window, values in sorted(data["windows"].items()):
            row = dict(base, window=window)
            row.update({k: v for k, v in values.items()})
            rows.append(row)
    fields = []
    for row in rows:
        for k in row:
            if k not in fields:
                fields.append(k)
    with (RESULTS / "temporal_runs.csv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)
    print(f"Wrote {RESULTS / 'temporal_runs.csv'} ({len(rows)} rows)")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = parser.add_subparsers(dest="command", required=True)
    p_run = sub.add_parser("run")
    p_run.add_argument("specs", nargs="+")
    p_run.add_argument("--local", action="store_true")
    p_run.add_argument("--version", default="v2512")
    p_run.add_argument("--hours", type=int, default=24)
    p_run.add_argument("--np", type=int, default=1,
                       help="MPI ranks per run (decomposes both regions)")
    p_run.add_argument("--env", default=os.environ.get("S4F_ENV_SCRIPT", ""),
                       help="script sourced with the version as argument to "
                       "load OpenFOAM and solids4foam in the job")
    p_an = sub.add_parser("analyse")
    p_an.add_argument("specs", nargs="*")
    sub.add_parser("orders")
    sub.add_parser("summary")
    args = parser.parse_args()
    if args.command == "run":
        if not args.env:
            parser.error("--env (or S4F_ENV_SCRIPT) is required")
        run([parse_spec(s) for s in args.specs], args.local, args.env,
            args.version, args.hours, args.np)
    elif args.command == "analyse":
        names = args.specs or sorted(p.name for p in RUN_ROOT.iterdir()
                                     if p.is_dir())
        for name in names:
            try:
                data = analyse_case(name)
                print(f"{name}: analysed ({data['execution_time']:.0f} s)")
            except wtv.ANALYSIS_ERRORS as error:
                print(f"{name}: {error}")
    elif args.command == "summary":
        summary()
    else:
        orders()
    return 0


if __name__ == "__main__":
    sys.exit(main())
