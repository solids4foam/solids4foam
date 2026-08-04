---
sort: 5
---

# unsNonLinGeomTotalLagSolid

This page documents the unstructured nonlinear total-Lagrangian
total-displacement solid model. The runtime type is:

```text
unsNonLinearGeometryTotalLagrangian
```

The model assumes finite strains and rotations and solves for the total
displacement field `D` in the total-Lagrangian frame, so the mesh is not
moved. It is the nonlinear-geometry counterpart of
[unsLinGeomSolid](https://www.solids4foam.com/documentation/solid-models/unsLinGeomSolid.html):
the "uns" prefix indicates that the face tangential gradient is calculated
with a face-Gauss approach rather than from the cell-centred gradient.

The stress is computed by the run-time selectable mechanical law.

---

## User Guide

### What this model solves

`unsNonLinGeomTotalLagSolid`:

- solves the momentum equation in total-Lagrangian form;
- uses total displacement `D` as the primary unknown;
- carries the kinematics as _surface_ fields — `Ff`, `Finvf`, `Jf` — as well
  as the cell-centred `F`, `Finv`, `J`;
- uses the surface stress `sigmaf` in the momentum equation;
- solves with a segregated implicit algorithm only;
- can fall back to a linear-geometry equation when the outer iterations
  diverge (see below).

The model is appropriate when deformation is large enough that linear geometry
is not valid, and the mesh is unstructured enough to warrant the face-Gauss
gradient.

### Supported solution algorithms

This model implements a single, segregated implicit algorithm. It does not
read `solutionAlgorithm`, and there is no PETSc SNES or explicit path. Use
[nonLinGeomTotalLagTotalDispSolid](https://www.solids4foam.com/documentation/solid-models/nonLinGeomTotalLagTotalDispSolid.html)
if you need those.

### Model options

| Entry | Default | Description |
| --- | --- | --- |
| `nonLinear` | `true` | Allow the automatic linear fallback |
| `debug` | `false` | Extra per-iteration solver output |
| `solutionTolerance` | inherited | Used here as a _relative_ tolerance |

```warning
`solutionTolerance` has a different meaning in this model. It is used as a
relative tolerance: the convergence target is the largest relative momentum
residual seen so far, multiplied by `solutionTolerance`, floored at the base
class `solutionTolerance`. A value that works for the other solid models may
converge sooner or later here.
```

The relevant inherited `solidModel` entries are:

| Entry | Default | Relevance |
| --- | --- | --- |
| `nCorrectors` | `10000` | Maximum number of outer correctors |
| `materialTolerance` | `1e-05` | Mechanical-law convergence tolerance |
| `dampingCoeff` | `0` | Adds `dampingCoeff*rho*ddt(D)` to the equation |
| `relaxationMethod` | `fixed` | Under-relaxation method (`fixed`, `aitken`) |
| `infoFrequency` | `100` | Frequency for solver progress output |
| `cellDisplacements` | optional | Internal-cell displacement constraints |

### The linear fallback

Nonlinear geometry can drive the outer iterations into a state where the
deformation gradient becomes unphysical. After each update the model calls
`checkEnforceLinear(Jf_)`; if the relative Jacobian indicates divergence, the
`enforceLinear` flag is set and the nonlinear part of the divergence-of-stress
term is replaced by its linear-geometry equivalent for the remainder of the
time step. The traction boundary condition switches to its linear form at the
same time.

When this happens, `evolve()` returns `false`, which signals the calling solver
to reduce the time step and try again. If `nonLinear` is set to `false`, the
check is skipped and the model runs without the fallback.

### Recommended dictionary setup

Minimal example for `constant/solidProperties`:

```text
solidModel        unsNonLinearGeometryTotalLagrangian;

unsNonLinearGeometryTotalLagrangianCoeffs
{
    nonLinear            yes;

    nCorrectors          10000;
    solutionTolerance    1e-06;
    materialTolerance    1e-05;
    infoFrequency        100;
}
```

The momentum equation is assembled with the name `laplacian(DD,D)`, so
`fvSchemes` must provide a `laplacian(DD,D)` entry, and `fvSolution` must
provide a solver for `D`.

### Required fields and boundary conditions

At minimum, the model expects:

- `D` in the time directory;
- `solidProperties` in `constant`;
- `g` in `constant`;
- a `mechanicalProperties` configuration for the chosen material law.

Traction boundaries are handled through `tractionBoundarySnGrad()`, which uses
the deformed-configuration normal `J*Finv.T() & n` from Nanson's relation, or
the reference normal `n` when the linear fallback is active.

### Field glossary

- `D`, `DD`: total and incremental displacement at cell centres.
- `pointD`, `pointDD`: total and incremental displacement at mesh points.
- `grad(D)`, `grad(DD)`: cell-centre gradients.
- `grad(D)f`: face gradient of `D`, computed with the face-Gauss approach.
- `F`, `Finv`, `J`: total deformation gradient, its inverse and Jacobian at
  cell centres.
- `Ff`, `Finvf`, `Jf`: the same quantities at faces; these are the ones the
  momentum equation uses.
- `sigma`, `sigmaf`: cell-centre and face stress.
- `U`: velocity field, computed as `ddt(D)`.

The cell-centred `F`, `Finv`, `J` and `sigma` are updated once at the end of
`evolve()`, for post-processing and for the mechanical law's cell-centred
interface; the loop itself runs entirely on the face fields.

---

## Developer Notes

### Class role

`unsNonLinGeomTotalLagSolid` is the total-Lagrangian, face-gradient member of
the family. The key design choices are:

- `D` is the primary solution variable (`solutionD()` returns `D()`);
- `nonLinGeom()` returns `TOTAL_LAGRANGIAN`, so the mesh is never moved;
- the momentum equation is written on the reference mesh, with
  `Jf*Finvf.T() & Sf` mapping the reference face areas to the deformed
  configuration;
- stress is delegated to `mechanicalModel`.

The class name in the header's `Class` tag is `nonLinGeomTotalLagSolid`, which
is a historical leftover; the actual class is `unsNonLinGeomTotalLagSolid`.

### Construction

The constructor:

1. calls the base `solidModel` constructor and `DisRequired()`;
2. creates `sigmaf_`, `gradDf_` and the kinematic fields, all of which are
   `READ_IF_PRESENT` so that a restart picks them up;
3. takes `impK_`, `impKf_` from the mechanical law and stores `rImpK_`;
4. reads `nonLinear`, `debug` and the relative `solutionTolerance`;
5. on restart, recomputes `grad(D)`, `grad(D)f`, `Ff`, `Finvf` and `Jf` from
   `D` and calls `mechanical().setRestart()`.

### `evolve()`

A single outer loop:

- reset `enforceLinear` and store `D` for under-relaxation;
- assemble the momentum equation

  ```text
  rho*d2dt2(D) == laplacian(impKf, D) - fvc::laplacian(impKf, D)
                + fvc::div((Jf*Finvf.T() & Sf) & sigmaf) + rho*g
  ```

  i.e. an implicit Laplacian with a deferred correction that recovers the
  total-Lagrangian divergence of stress at convergence;
- add damping if `dampingCoeff` is non-zero;
- if `enforceLinear` is set, add the difference between the nonlinear and
  linear divergence terms, which cancels the nonlinear contribution;
- relax the linear system, apply `solidModel::setCellDisps()`, and solve;
- under-relax `D`, then update `pointD`, `grad(D)`, `grad(D)f`, `grad(DD)`,
  `Ff`, `Finvf` and `Jf`;
- call `checkEnforceLinear(Jf_)` unless the fallback is already active or
  `nonLinear` is `false`;
- update `sigmaf` and evaluate the relative residual.

Convergence is checked against `maxRes*relativeTol_`, floored at the base-class
`solutionTolerance`; the first iteration is always followed by a second. After
the loop, `U`, `F`, `Finv`, `J`, `sigma`, `DD` and `pointDD` are updated once.

`evolve()` returns `false` if the linear fallback was triggered, and `true`
otherwise.

### Traction boundary treatment

`tractionBoundarySnGrad()` has two branches. In the nonlinear branch the
traction and pressure are applied on the deformed normal:

```text
((traction - nCurrent*pressure) - (nCurrent & sigma) + (n & (impK*gradD)))*rImpK
```

with `nCurrent = J*Finv.T() & n`. In the `enforceLinear` branch, `nCurrent` is
replaced by the reference normal `n`. In both cases the `impK*gradD` term
removes the part of the traction that the implicit Laplacian already accounts
for.

### Extension points

- `evolve()` for a different outer-iteration or fallback strategy;
- `tractionBoundarySnGrad()` if the deferred-correction split changes.

As with `unsLinGeomSolid`, the face fields are the primary state: any new
algorithm must keep `sigmaf_`, `gradDf_`, `Ff_`, `Finvf_` and `Jf_`
consistent, not just their cell-centred counterparts.

---

## Tutorials

No tutorial selects `unsNonLinearGeometryTotalLagrangian` by default. The
finite-strain cases that exercise the sibling models —
`solids/hyperelasticity/cantileverBeam`,
`solids/hyperelasticity/cylindricalPressureVessel` and
`solids/elastoplasticity/cooksMembrane` — are the natural places to try it, by
changing the `solidModel` entry in `constant/solidProperties`.
