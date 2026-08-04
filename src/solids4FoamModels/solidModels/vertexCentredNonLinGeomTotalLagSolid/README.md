---
sort: 8
---

# vertexCentredNonLinGeomTotalLagSolid

This page documents the vertex-centred nonlinear total-Lagrangian solid model.
The runtime type is:

```text
vertexCentredNonLinTotalLagGeometry
```

The model assumes finite strains and rotations and solves for the total
displacement at the mesh **points**, `pointD`, in the total-Lagrangian frame,
so the mesh is not moved. As with
[vertexCentredLinGeomSolid](https://www.solids4foam.com/documentation/solid-models/vertexCentredLinGeomSolid.html),
the governing equations are integrated over a dual mesh built automatically by
`meshDualiser`, and the implicit system is solved with PETSc SNES.

```note
This model works best on triangular and tetrahedral grids. It requires PETSc
for the implicit path.
```

---

## User Guide

### What this model solves

`vertexCentredNonLinGeomTotalLagSolid`:

- solves the momentum equation in total-Lagrangian form over the dual mesh,
  with the primary unknown at the mesh points;
- carries the kinematics at the **dual mesh faces** — `dualFf`, `dualFinvf`,
  `dualJf` — alongside `dualGradDf` and `dualSigmaf`;
- assembles the material tangent _and_, optionally, the geometric stiffness
  into the Jacobian;
- can add pressure as an extra unknown per point, for near-incompressible
  materials;
- enforces displacement boundary conditions as fixed degrees of freedom;
- also offers an explicit central-difference path.

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
| `useGeometricStiffness` | none | Geometric stiffness in the Jacobian |
| `zeta` | `1.0` | Compact edge gradient fraction |
| `zetaImplicit` | `zeta` | The `zeta` used to form the Jacobian |
| `approximateJacobian` | `false` | Compact Laplacian instead of the tangent |
| `compactImplicitStencil` | none | Compact stencil for that Jacobian |
| `tangentEps` | `1e-10` | Geometric stiffness perturbation size |
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

`useGeometricStiffness` is read with `lookup()`, so it has no default and the
model fails at construction if it is missing. Including the geometric
stiffness costs one Jacobian assembly pass but usually pays for itself in
Newton iterations on strongly nonlinear problems.

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
solidModel     vertexCentredNonLinTotalLagGeometry;

vertexCentredNonLinTotalLagGeometryCoeffs
{
    // Solution algorithm (PETScSNES, explicit)
    solutionAlgorithm     PETScSNES;

    // Include the geometric stiffness in the Jacobian
    useGeometricStiffness yes;

    // Use the exact Jacobian or a small-stencil approximation
    approximateJacobian   no;

    // Compact edge discretisation fraction
    // 0 -> more accurate but oscillations more likely
    // 1 -> less accurate but oscillations less likely
    zeta                  0.1;
}
```

A `petscOptions` file must be present in the case directory, unless the name is
overridden with `optionsFile`.

### Required fields and boundary conditions

The boundary conditions are set on `pointD`, a `pointVectorField`, not on the
cell-centred `D`. The available point patch types are the same as for the
linear-geometry vertex-centred model:

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

```warning
`tractionBoundarySnGrad()` is `notImplemented` for this model, as is
`solutionD()`; both are cell-centred concepts.
```

### Field glossary

- `pointD`: total displacement at the mesh points; the primary unknown.
- `pointP`: point pressure, when `solvePressure` is enabled.
- `pointDD`, `pointU`, `pointA`, `pointRho`: point increment, velocity,
  acceleration and density.
- `pointVol`, `pointGlobalVol`: dual cell volume per point, local and with
  processor-boundary contributions summed.
- `pointDivSigma`: divergence of stress at each point; the momentum residual.
- `dualGradDf`, `dualSigmaf`: displacement gradient and stress at the dual
  faces.
- `dualFf`, `dualFinvf`, `dualJf`: total deformation gradient, its inverse and
  Jacobian at the dual faces.
- `dualPf`, `volP`: pressure at the dual faces and at cell centres.
- `D`, `sigma`: cell-centred fields, interpolated for post-processing.

### Explicit path

With `solutionAlgorithm explicit`, `setDeltaT()` computes the stable time step
from the elastic wave speed and the dual mesh `deltaCoeffs`, scaled by `maxCo`
from `controlDict` (default `0.1`).

---

## Developer Notes

### Class role

`vertexCentredNonLinGeomTotalLagSolid` inherits from `solidModel` and, when
PETSc is available, from `foamPetscSnesHelper`. It is the finite-strain
counterpart of `vertexCentredLinGeomSolid` and shares its structure: dual
mesh, `dualMechanicalModel`, fixed degrees of freedom, and a PETSc SNES
driver.

- `nonLinGeom()` returns `TOTAL_LAGRANGIAN`, so the mesh is never moved;
- the unknowns live at points, so `solutionD()` is `notImplemented`;
- the extra state relative to the linear model is the dual-face kinematics and
  the optional pressure unknown.

### PETSc SNES path

`evolveSnes()` drives `foamPetscSnesHelper`, which calls back into:

- `formResidual()`, which extracts `pointD` (and `pointP` when solving
  pressure) from the solution vector, corrects the boundary conditions,
  computes `dualGradDf` via `vfvc::fGrad` with `zeta`, forms `dualFf`,
  `dualFinvf` and `dualJf`, evaluates `dualSigmaf`, and assembles
  `pointDivSigma` on the deformed dual areas `J*Finv.T() & Sf0`;
- `formJacobian()`, which adds either the exact linearisation of `div(sigma)`
  through `vfvm::divSigma` using `materialTangentFaceField()`, or a compact
  Laplacian through `vfvm::laplacian` when `approximateJacobian` is set, and
  then — if `useGeometricStiffness` is enabled — the geometric stiffness.

### Geometric stiffness

`geometricStiffnessField()` builds the contribution to the tangent that comes
from the deformed area vector `gamma = J*F^-T & Sf0` changing with the
displacement. It is formed **numerically**: for each of the nine components of
`gradD` in turn, the component is perturbed by `tangentEps`, `F`, `Finv`, `J`
and the deformed `Sf` are recomputed, and the tangent column is the finite
difference `(SfPerturb - SfRef)/eps`.

Because it is a finite difference, `tangentEps` trades truncation error
against round-off. The default `1e-10` suits displacement gradients of order
one; on problems with very small or very large strains it may need adjusting.

### Extension points

- `formResidual()` and `formJacobian()` for a different discretisation;
- `geometricStiffnessField()` if an analytical geometric tangent is derived,
  which would remove the `tangentEps` sensitivity;
- `setFixedDofs()` and `enforceTractionBoundaries()` for new point patch
  types.

---

## Tutorials

No tutorial currently selects `vertexCentredNonLinTotalLagGeometry`. The
finite-strain cases under `solids/hyperelasticity` are the natural place to
try it, but note that they use cell-centred `D` boundary conditions, so the
`0/pointD` field and its point patch types have to be set up first.
