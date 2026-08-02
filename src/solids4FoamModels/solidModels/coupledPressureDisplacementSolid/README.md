---
sort: 9
---

# coupledPressureDisplacementSolid

This page documents the block-coupled pressure-displacement solid model. The
runtime type is:

```text
coupledPressureDisplacementSolid
```

The model is aimed at **incompressible and nearly incompressible** solids,
where a displacement-only formulation locks. It solves the displacement
increment `DD` and the hydrostatic pressure increment `Dp` together in one
block system of four unknowns per cell, in an updated-Lagrangian frame with
the mesh moved to the current configuration.

```warning
This model is **only available with foam-extend**. It is listed in
`Make/files.foamextend` but not in `Make/files.openfoam`, so it does not exist
at all in OpenFOAM.com and OpenFOAM.org builds. Its tutorials call
`solids4Foam::caseOnlyRunsWithFoamExtend`.
```

---

## User Guide

### What this model solves

`coupledPressureDisplacementSolid`:

- assembles a `volVector4Field` unknown, `DDDp`, holding the three components
  of the displacement increment `DD` plus the pressure increment `Dp`;
- inserts the momentum equation at rows 0-2 and the pressure equation at row
  3, with the two coupling blocks — the pressure gradient in the momentum
  equation, and the displacement divergence in the pressure equation —
  inserted implicitly;
- uses a Rhie-Chow-style flux `phi` to avoid pressure-velocity decoupling on
  a collocated mesh, with an optional consistent variant;
- moves the mesh to the current configuration each outer iteration when
  `nonLinear` is enabled;
- supports both first-order Euler and backward temporal schemes, and a
  steady-state mode.

The pressure unknown is what distinguishes this from the other coupled model,
[coupledUnsLinGeomLinearElasticSolid](https://www.solids4foam.com/documentation/solid-models/coupledUnsLinGeomLinearElasticSolid.html):
that one is displacement-only and linear-elastic, this one is a general
finite-strain mixed formulation.

### Supported solution algorithms

There is a single block-coupled path. This model does not read
`solutionAlgorithm`.

### Model options

| Entry | Default | Description |
| --- | --- | --- |
| `nonLinear` | `true` | Nonlinear geometry and mesh motion |
| `K` | `0` | Mass-proportional damping coefficient |
| `consistentRhieChow` | `false` | Consistent Rhie-Chow correction |
| `composite` | `false` | Composite-material mode |
| `stdDispGrad` | `true` | Standard displacement-gradient evaluation |
| `debug` | `false` | Extra per-iteration solver output |
| `writeConvergenceData` | `false` | Write `convergData.dat` in the case |
| `solutionTolerance` | inherited | Used here as a _relative_ tolerance |

`K` has dimensions `[0 0 -1 0 0 0 0]`. `consistentRhieChow` allocates extra
interpolated face fields when enabled.

```warning
As in `unsNonLinGeomTotalLagSolid`, `solutionTolerance` is used as a relative
tolerance in this model, not an absolute one.
```

`K` is unusual in that it is re-read from the dictionary inside
`updateTotalFields()`, so it can be changed while a case is running — useful
for ramping damping down as a quasi-static solution settles.

The relevant inherited `solidModel` entries are:

| Entry | Default | Relevance |
| --- | --- | --- |
| `nCorrectors` | `10000` | Maximum number of outer correctors |
| `alternativeTolerance` | `1e-07` | Secondary convergence tolerance |
| `materialTolerance` | `1e-05` | Mechanical-law convergence tolerance |
| `infoFrequency` | `100` | Frequency for solver progress output |
| `restart` | `false` | Writes extra fields needed for a consistent restart |
| `stabilisation` | auto-created | Both sub-dictionaries are used |
| `cellDisplacements` | optional | Internal-cell displacement constraints |

The stabilisation sub-models accept a `scaleFactorJacobian` alongside
`scaleFactor`, so the stabilisation used when forming the block matrix can
differ from the one used in the residual. Setting both to zero, as the
`heartTissueBeam` tutorial does, disables stabilisation entirely; the comment
there notes that momentum stabilisation contributes to `sigmaHyd`
oscillations.

### Required input files

Because the model is incremental and mixed, the fields differ from the other
solid models:

- `DD`, the displacement increment, with the boundary conditions;
- `Dp`, the pressure increment (`MUST_READ`);
- `pointDD`, usually `calculated`;
- any material-specific fields, such as a fibre direction `f0`;
- `constant/solidProperties`, `constant/mechanicalProperties`, `constant/g`.

Note that `D` and `p` are outputs, accumulated from `DD` and `Dp`.

### Recommended dictionary setup

Example for `constant/solidProperties`:

```text
solidModel coupledPressureDisplacementSolid;

coupledPressureDisplacementSolidCoeffs
{
    nCorrectors             2000;
    solutionTolerance       1e-7;
    alternativeTolerance    1e-06;
    materialTolerance       1e-05;
    infoFrequency           1;

    stabilisation
    {
        momentum
        {
            type                laplacian;
            scaleFactor         0;
            scaleFactorJacobian 0;
        }

        pressure
        {
            type                laplacian;
            scaleFactor         0;
            scaleFactorJacobian 0;
        }
    }

    nonLinear             true;
    debug                 false;

    consistentRhieChow    false;
    composite             false;
    stdDispGrad           true;

    restart               true;

    K                     K [0 0 -1 0 0 0 0] 150;

    writeConvergenceData  true;
}
```

The block system is solved with a `DDDp` solver in `fvSolution`, not with
separate `DD` and `Dp` solvers:

```text
solvers
{
    DDDp
    {
        solver BiCGStab;

        preconditioner
        {
            preconditioner Cholesky;
        }
    }
}
```

### Boundary conditions

Displacement boundaries are set on `DD`. The type designed for this model is
`tractionPressureDisplacement`, which applies a traction consistently with the
pressure unknown; `symmetryPlane` and the standard fixed-displacement types
also work. Pressure boundaries on `Dp` are usually `zeroGradient`.

### Field glossary

- `DDDp`: the block unknown; components 0-2 are `DD`, component 3 is `Dp`.
- `DD`, `D`: displacement increment and accumulated total displacement.
- `Dp`, `p`: pressure increment and accumulated pressure.
- `DU`: velocity increment.
- `phi`: the Rhie-Chow flux increment.
- `delta`, `weights`: geometric quantities cached across the mesh motion,
  written with `AUTO_WRITE` so restarts are consistent.
- `force`: the assembled surface force used in the momentum equation.
- `pointDD`, `pointD`, `U`, `sigma`: as for the other incremental models.

---

## Developer Notes

### Class role

`coupledPressureDisplacementSolid` inherits directly from `solidModel`. The
key design choices are:

- `DD` is the primary solution variable (`DDisRequired()` is called), with
  `Dp` as a coupled second unknown;
- the block matrix is a foam-extend `fvBlockMatrix<vector4>`, which is why the
  class is foam-extend-only;
- the mesh is moved to the current configuration inside the iteration when
  `nonLinear` is set, so the formulation is updated-Lagrangian in effect while
  the equations are assembled in total-Lagrangian form on the reference
  points, kept in `refPoints_`.

The implementation is split across a set of `.H` include files —
`calcPhi.H`, `calcTraction.H`, `calcForce.H`, `completeMomentumEqAssembly.H`,
`addBlockCoupledBC.H`, `finalizeMomentumEqn.H`, `updateAD.H`,
`solidContinuityErrs.H` and others — which are `#include`d directly into
`evolve()`. This keeps `evolve()` readable but means the control flow is
spread over a dozen files.

### `evolve()`

Each outer iteration:

1. computes the tangential and normal parts of the displacement-increment
   gradient at faces, and the traction (`calcTraction.H`);
2. assembles the momentum equation as a Laplacian in `DD` with explicit
   corrections for the tangential gradient, the nonlinear force correction and
   the old-time force, plus the temporal terms — Euler, backward or
   steady-state — and the optional `K` damping;
3. inserts it into the block matrix at row 0, completes the assembly with the
   normal-derivative term, applies the block-coupled boundary corrections, and
   inserts the pressure-gradient coupling block at (0, 3);
4. moves the mesh to the current configuration if `nonLinear`;
5. computes the Rhie-Chow flux (`calcPhi.H`), assembles the pressure equation
   as `Sp(rKappa, Dp) - laplacian(rADf, Dp) + div(phi)`, inserts it at row 3,
   and inserts the displacement-divergence coupling block at (3, 0);
6. solves the block system, retrieves `DD` and `Dp`, corrects their boundary
   conditions, updates the flux and reports the continuity errors;
7. accumulates `D = D.oldTime() + DD` and `p = p.oldTime() + Dp`.

### `updateTotalFields()`

Refreshes `impKf_` from the mechanical law, re-reads `K` from the dictionary,
and — when `nonLinear` is set — updates the accumulated configuration. This is
the hook that makes `K` adjustable mid-run.

### Extension points

- the pressure equation and `calcPhi.H` if a different
  pressure-velocity coupling is wanted;
- `completeMomentumEqAssembly.H` for changes to the implicit momentum terms;
- `addBlockCoupledBC.H` for a new block-aware boundary condition.

The `composite` and `consistentRhieChow` switches guard code paths that are
less exercised than the defaults; treat them as experimental.

---

## Tutorials

Cases that select `coupledPressureDisplacementSolid`, all foam-extend only:

- `solids/hyperelasticity/heartTissueBeam`
- `solids/hyperelasticity/ratCarotid`
- `solids/hyperelasticity/idealisedVentricle`
- `solids/hyperelasticity/cylindricalPressureVessel`, in its
  `caseOptions/pressureDisplacement` variants
- `solids/linearElasticity/plateHole`, in its
  `caseOptions/pressureDisplacement` variant
