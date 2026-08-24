---
sort: 6
---

# unsNonLinGeomUpdatedLagSolid

This page documents the unstructured nonlinear updated-Lagrangian
incremental-displacement solid model. The runtime type is:

```text
unsNonLinearGeometryUpdatedLagrangian
```

The model assumes finite strains and rotations and solves for the increment of
displacement `DD`. The mesh is moved at the end of each time step, so the
governing equations are written in the updated configuration. It is the
updated-Lagrangian counterpart of
[unsLinGeomSolid](https://www.solids4foam.com/documentation/solid-models/unsLinGeomSolid.html):
the "uns" prefix indicates that the face tangential gradient is calculated
with a face-Gauss approach rather than from the cell-centred gradient.

The stress is computed by the run-time selectable mechanical law.

---

## User Guide

### What this model solves

`unsNonLinGeomUpdatedLagSolid`:

- solves the momentum equation in updated-Lagrangian form;
- uses the displacement increment `DD` as the primary unknown, and recovers
  `D` as `D.oldTime() + DD`;
- carries the kinematics as _surface_ fields — `relFf`, `relFinvf`, `relJf`,
  `Ff`, `Jf` — as well as their cell-centred counterparts;
- uses the surface stress `sigmaf` in the momentum equation;
- maintains its own density field `rho`, updated as `rho.oldTime()/relJ`;
- moves the mesh at the end of each time step;
- solves with a segregated implicit algorithm only.

The model is appropriate when deformation is large enough that linear geometry
is not valid and mesh motion between time steps is part of the formulation.

### Supported solution algorithms

This model implements a single, segregated implicit algorithm. It does not
read `solutionAlgorithm`, and there is no PETSc SNES or explicit path. Use
[nonLinGeomUpdatedLagSolid](https://www.solids4foam.com/documentation/solid-models/nonLinGeomUpdatedLagSolid.html)
if you need those.

### Model options

`unsNonLinGeomUpdatedLagSolid` adds no options of its own. Unlike its
total-Lagrangian sibling it has no `nonLinear` switch and no automatic linear
fallback, and `solutionTolerance` keeps its usual absolute meaning. The
relevant entries in `unsNonLinearGeometryUpdatedLagrangianCoeffs` are the
inherited `solidModel` ones:

| Entry | Default | Relevance |
| --- | --- | --- |
| `nCorrectors` | `10000` | Maximum number of outer correctors |
| `solutionTolerance` | `1e-06` | Primary convergence tolerance |
| `alternativeTolerance` | `1e-07` | Secondary convergence tolerance |
| `materialTolerance` | `1e-05` | Mechanical-law convergence tolerance |
| `relaxationMethod` | `fixed` | Under-relaxation method (`fixed`, `aitken`) |
| `infoFrequency` | `100` | Frequency for solver progress output |
| `writeResidualField` | `false` | Writes a residual field during output |
| `residualFile` | `false` | Writes `residual.dat` in the case directory |
| `cellDisplacements` | optional | Internal-cell displacement constraints |

```note
`dampingCoeff` is read by the base class but is not applied by this model's
momentum equation.
```

### Recommended dictionary setup

Minimal example for `constant/solidProperties`:

```text
solidModel        unsNonLinearGeometryUpdatedLagrangian;

unsNonLinearGeometryUpdatedLagrangianCoeffs
{
    nCorrectors          10000;
    solutionTolerance    1e-06;
    alternativeTolerance 1e-07;
    materialTolerance    1e-05;
    infoFrequency        100;
}
```

The momentum equation is assembled with the name `laplacian(DDD,DD)`, so
`fvSchemes` must provide a `laplacian(DDD,DD)` entry, and `fvSolution` must
provide a solver for `DD`.

### Required fields and boundary conditions

Because this model is incremental, the boundary conditions are specified on
`DD`, not on `D`:

- `DD` in the time directory;
- `solidProperties` in `constant`;
- `g` in `constant`;
- a `mechanicalProperties` configuration for the chosen material law.

The base class enforces this through `DDisRequired()`. Traction boundaries are
handled through `tractionBoundarySnGrad()`, which uses the deformed normal
`relJ*relFinv.T() & n` obtained from Nanson's relation between the updated and
deformed configurations.

### Field glossary

- `DD`: displacement increment at cell centres; the primary unknown.
- `D`: total displacement, recovered as `D.oldTime() + DD`.
- `pointDD`, `pointD`: incremental and total displacement at mesh points.
- `grad(DD)`, `grad(D)`: cell-centre gradients.
- `grad(DD)f`: face gradient of `DD`, computed with the face-Gauss approach.
- `relF`, `relFinv`, `relJ`: relative deformation gradient, its inverse and
  relative Jacobian, mapping the updated to the deformed configuration.
- `F`, `J`: total deformation gradient and Jacobian, accumulated as
  `relF & F.oldTime()` and `relJ*J.oldTime()`.
- `relFf`, `relFinvf`, `relJf`, `Ff`, `Jf`: the same quantities at faces;
  these are the ones the momentum equation uses.
- `rho`: density in the current configuration, updated as `rho.oldTime()/relJ`.
- `sigma`, `sigmaf`: cell-centre and face stress.
- `U`: velocity field, computed as `ddt(D)`.

---

## Developer Notes

### Class role

`unsNonLinGeomUpdatedLagSolid` is the updated-Lagrangian, face-gradient member
of the family. The key design choices are:

- `DD` is the primary solution variable (`solutionD()` returns `DD()`, and
  `incremental()` returns `true`);
- `nonLinGeom()` returns `UPDATED_LAGRANGIAN`;
- `movingMesh()` returns `true`, so the base class expects the mesh to be
  moved once per time step;
- the momentum equation is written on the updated mesh, with
  `relJf*relFinvf.T() & Sf` mapping the updated face areas to the deformed
  configuration;
- stress is delegated to `mechanicalModel`.

The class name in the header's `Class` tag is `nonLinGeomUpdatedLagSolid`,
which is a historical leftover; the actual class is
`unsNonLinGeomUpdatedLagSolid`.

### Construction

The constructor:

1. calls the base `solidModel` constructor and `DDisRequired()`;
2. creates `sigmaf_`, `gradDDf_`, the relative and total kinematic fields, and
   its own `rho_`, all `READ_IF_PRESENT` so that a restart picks them up;
3. takes `impK_`, `impKf_` from the mechanical law and stores `rImpK_`;
4. on restart, recomputes `grad(DD)`, `grad(DD)f`, `relFf`, `relFinvf`, `Ff`,
   `relJf` and `Jf`, stores the old times of `Ff` and `Jf`, and calls
   `mechanical().setRestart()`.

### `evolve()`

A single outer loop:

- store `DD` for under-relaxation and residual calculation;
- assemble the momentum equation

  ```text
  d2dt2(rho, DD) + fvc::d2dt2(rho, D.oldTime())
      == laplacian(impKf, DD) - fvc::laplacian(impKf, DD)
       + fvc::div((relJf*relFinvf.T() & Sf) & sigmaf) + rho*g
  ```

  i.e. an implicit Laplacian with a deferred correction that recovers the
  updated-Lagrangian divergence of stress at convergence;
- relax the linear system, apply `solidModel::setCellDisps()`, and solve;
- under-relax `DD`, then update `pointDD`, `grad(DD)`, `grad(DD)f`, `relFf`,
  `relFinvf`, `Ff`, `relJf` and `Jf`;
- update `sigmaf` from the mechanical law.

The loop exits through the base-class `converged()` helper, or when `iCorr`
reaches `nCorrectors`. After the loop, the cell-centred `relF`, `relFinv`,
`F`, `relJ`, `J` and `sigma` are updated once, then `D`, `grad(D)`, `pointD`
and `U`.

`evolve()` always returns `true`; there is no divergence fallback in this
model.

### `updateTotalFields()`

Called once per time step by the base class after `evolve()`. It:

1. updates the density as `rho = rho.oldTime()/relJ`;
2. moves the mesh by `pointDD` through `moveMesh()`, so that the next time
   step starts in the new updated configuration;
3. calls `solidModel::updateTotalFields()`.

This is the step that makes the formulation updated- rather than
total-Lagrangian, and it is why `movingMesh()` must return `true`.

### Traction boundary treatment

`tractionBoundarySnGrad()` applies the traction and pressure on the deformed
normal `nCurrent = relJ*relFinv.T() & n`, where `n` is the normal in the
updated configuration, and subtracts the `impK*gradDD` term that the implicit
Laplacian already accounts for.

### Extension points

- `evolve()` for a different outer-iteration strategy;
- `updateTotalFields()` if the mesh-motion or density update changes;
- `tractionBoundarySnGrad()` if the deferred-correction split changes.

As with the other `uns` models, the face fields are the primary state: any new
algorithm must keep `sigmaf_`, `gradDDf_`, `relFf_`, `relFinvf_`, `Ff_`,
`relJf_` and `Jf_` consistent, not just their cell-centred counterparts.

---

## Tutorials

Cases that select `unsNonLinearGeometryUpdatedLagrangian`:

- `solids/hyperelasticity/cantileverBeam`

It is offered as a commented alternative in several further finite-strain
cases, which are the natural places to try it:

- `solids/hyperelasticity/cylindricalPressureVessel`
- `solids/hyperelasticity/longWall`
- `solids/elastoplasticity/cooksMembrane`
- `solids/elastoplasticity/pipeCrush`
- `solids/elastoplasticity/cylinderCrush`
