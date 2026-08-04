---
sort: 7
---

# vertexCentredLinGeomSolid

This page documents the vertex-centred linear-geometry solid model. The
runtime type is:

```text
vertexCentredLinearGeometry
```

The model assumes small strains and small rotations and solves for the total
displacement at the mesh **points**, `pointD`, rather than at cell centres.
The governing equations are integrated over a _dual_ mesh, which is
constructed automatically from the primary mesh by `meshDualiser`. The linear
system couples all displacement components implicitly and is solved with
PETSc SNES.

The method is described in
[10.13140/RG.2.2.22896.33283](https://doi.org/10.13140/RG.2.2.22896.33283).

```note
This model works best on triangular and tetrahedral grids. It requires PETSc
for the implicit path.
```

---

## User Guide

### What this model solves

`vertexCentredLinGeomSolid`:

- solves the linear momentum equation in linear-geometry form over the dual
  mesh, with the primary unknown at the mesh points;
- computes the displacement gradient and the stress at the **dual mesh faces**
  (`dualGradDf`, `dualSigmaf`), through a `dualMechanicalModel`;
- assembles the exact material tangent into the Jacobian by default, giving
  Newton-Raphson convergence;
- enforces displacement boundary conditions as fixed degrees of freedom rather
  than through patch-field manipulation;
- also offers an explicit central-difference path.

Being vertex-centred, it avoids the checkerboarding and the stabilisation
terms that the cell-centred segregated models need, at the cost of a larger
implicit system and a dependency on PETSc.

### Supported solution algorithms

The `solutionAlgorithm` entry selects the path:

| Value | Meaning | Notes |
| --- | --- | --- |
| `PETScSNES` | Newton-Raphson via PETSc SNES | The intended path |
| `explicit` | Explicit time integration | Time-step limited by stability |
| `implicitSegregated` | Not implemented | Fatal error; use `PETScSNES` |
| `implicitCoupled` | Not implemented | Fatal error; use `PETScSNES` |

```warning
The base-class default is `implicitSegregated`, which this model rejects with
a fatal error. Every case **must** set `solutionAlgorithm` explicitly.
```

### Model options

| Entry | Default | Description |
| --- | --- | --- |
| `zeta` | `1.0` | Compact edge gradient fraction |
| `zetaImplicit` | `zeta` | The `zeta` used to form the Jacobian |
| `approximateJacobian` | `false` | Compact Laplacian instead of the tangent |
| `compactImplicitStencil` | none | Compact stencil for that Jacobian |
| `fixedDofScale` | see below | Scaling of the constraint equations |
| `predictor` | `false` | Predict `pointD` at a new time step |
| `predictorMethod` | `linear` | `linear` or `quadratic` predictor |
| `optionsFile` | `petscOptions` | PETSc options file name |
| `stopOnPetscError` | `true` | Stop if PETSc reports an error |

`zeta` is the compact edge gradient fraction: `0` is more accurate but more
prone to oscillations, `1` is less accurate but more robust. It is the entry
most worth tuning, and the shipped tutorial uses `0.1`. `zetaImplicit`
defaults to `zeta`, and lets the Jacobian use a different value from the
residual.

`compactImplicitStencil` has no default and is only read when
`approximateJacobian` is enabled, in which case it must be present.

`fixedDofScale` defaults to `average(impK)*sqrt(gAverage(magSf))`, so that the
constraint equations are conditioned like the momentum equations.

The `linear` predictor assumes constant velocity over the step; `quadratic`
assumes constant acceleration.

The relevant inherited `solidModel` entries are:

| Entry | Default | Relevance |
| --- | --- | --- |
| `solutionAlgorithm` | `implicitSegregated` | Must be overridden, see above |
| `solvePressure` | `false` | Adds pressure as an extra unknown per point |
| `infoFrequency` | `100` | Frequency for solver progress output |

The block size of the PETSc system is 2 in 2-D and 3 in 3-D, or 3 and 4
respectively when `solvePressure` is enabled.

### Recommended dictionary setup

Minimal example for `constant/solidProperties`:

```text
solidModel     vertexCentredLinearGeometry;

vertexCentredLinearGeometryCoeffs
{
    // Solution algorithm (PETScSNES, explicit)
    solutionAlgorithm PETScSNES;

    // Use the exact Jacobian or a small-stencil approximation
    approximateJacobian no;

    // Compact edge discretisation fraction
    // 0 -> more accurate but oscillations more likely
    // 1 -> less accurate but oscillations less likely
    zeta            0.1;
}
```

A `petscOptions` file must be present in the case directory, unless the name is
overridden with `optionsFile`.

### Required fields and boundary conditions

The boundary conditions are set on `pointD`, a `pointVectorField`, not on the
cell-centred `D`. Displacement constraints become fixed degrees of freedom in
the linear system.

The point patch types provided by solids4foam are:

| Patch type | Purpose |
| --- | --- |
| `pointSolidTraction` | Prescribed traction and pressure |
| `pointNormalDisplacement` | Prescribed normal displacement |
| `pointFixedDisplacementZeroShear` | Normal displacement fixed, zero shear |
| `pointFixedRotation` | Prescribed rotation about an axis |
| `pointLinearSpatialDisplacement` | Displacement varying linearly in space |
| `pointSolidForce` | Prescribed force |
| `pointSolidContact` | Contact |
| `componentMixed` | Mixed per-component condition |

Ordinary `fixedValue` and `zeroGradient` point patch types work as usual.

```warning
`tractionBoundarySnGrad()` is `notImplemented` for this model, as is
`solutionD()`; both are cell-centred concepts. Do not expect cell-centred
traction boundary conditions to work.
```

### Field glossary

- `pointD`: total displacement at the mesh points; the primary unknown.
- `pointDD`: point displacement increment, `pointD - pointD.oldTime()`.
- `pointU`, `pointA`: point velocity and acceleration.
- `pointRho`: point density.
- `pointVol`: dual cell volume associated with each point. On a processor
  boundary this is the local contribution only.
- `pointGlobalVol`: the same, with processor-boundary contributions summed.
- `pointDivSigma`: divergence of stress at each point; the residual of the
  momentum equation.
- `dualGradDf`, `dualSigmaf`: displacement gradient and stress at the dual
  mesh faces.
- `D`, `sigma`: cell-centred fields, interpolated from the point fields for
  post-processing.

### Explicit path

With `solutionAlgorithm explicit`, `setDeltaT()` computes the stable time step
from the elastic wave speed and the dual mesh `deltaCoeffs`, scaled by `maxCo`
from `controlDict` (default `0.1`).

---

## Developer Notes

### Class role

`vertexCentredLinGeomSolid` inherits from `solidModel` and, when PETSc is
available, from `foamPetscSnesHelper`. The key design choices are:

- the unknowns live at points, so `solutionD()` is `notImplemented`;
- the dual mesh comes from the base class's `dualMesh()` and
  `dualMeshMap()`, built by `meshDualiser`;
- constitutive evaluation goes through a `dualMechanicalModel`, which maps
  dual faces back to primary cells so that ordinary mechanical laws can be
  reused unchanged;
- displacement constraints are held as `fixedDofs_`, `fixedDofValues_` and
  `fixedDofDirections_`, all sized by the number of points, rather than being
  applied through patch fields.

### Construction

The constructor initialises the PETSc SNES helper with `pointD` as the
solution field and `solutionLocation::POINTS`, builds the
`dualMechanicalModel`, computes `blockSize_`, allocates the fixed-degree-of-
freedom lists, and reads `fixedDofScale` — whose default is derived from the
average implicit stiffness and the mesh size, so that the constraint equations
are conditioned like the momentum equations.

### PETSc SNES path

`evolveSnes()` drives `foamPetscSnesHelper`, which calls back into:

- `initialiseSolution()` and `initialiseJacobian()` to size the vector and set
  the matrix's non-zero structure;
- `formResidual()`, which extracts `pointD` from the solution vector, corrects
  the boundary conditions, computes `dualGradDf` via `vfvc::fGrad` with
  `zeta`, evaluates `dualSigmaf`, and assembles `pointDivSigma`;
- `formJacobian()`, which either adds the exact linearisation of `div(sigma)`
  through `vfvm::divSigma` using the material tangent from
  `materialTangentFaceField()`, or — when `approximateJacobian` is set — adds
  a compact Laplacian through `vfvm::laplacian` with `zetaImplicit` and
  `dualImpKf()`.

Fixed degrees of freedom are removed from the system through
`fixedDofRowsIS()`, a PETSc index set built from `fixedDofs_`.

The exact Jacobian gives quadratic convergence but is more expensive to form;
the approximate one is cheaper per iteration and needs more of them.

### Explicit path implementation

`evolveExplicit()` advances `pointU` and `pointD` with a central-difference
update, using `pointDivSigma` and `pointGlobalVol`. `setDeltaT()` handles
stability.

### Extension points

- `formResidual()` and `formJacobian()` for a different discretisation;
- `setFixedDofs()` for a new kind of displacement constraint;
- `enforceTractionBoundaries()` for a new traction-style point patch.

Note the distinction between `pointVol_` and `pointGlobalVol_`: use the global
one wherever a physical volume is needed, and the local one only where
processor-local partial sums are correct.

---

## Tutorials

No tutorial selects `vertexCentredLinearGeometry` by default. Two cases are
set up for it:

- `solids/linearElasticity/cantilever2d`, which ships
  `constant/solidProperties.vertexCentred` alongside its other variants;
- `solids/linearElasticity/wobblyNewton`, where it is a commented alternative.
