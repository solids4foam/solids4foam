---
sort: 4
---

# unsLinGeomSolid

This page documents the unstructured linear-geometry total-displacement solid
model. The runtime type is:

```text
unsLinearGeometry
```

The model assumes small strains and small rotations, and solves for the total
displacement field `D`. It differs from
[linGeomTotalDispSolid](https://www.solids4foam.com/documentation/solid-models/linGeomTotalDispSolid.html)
in how the gradient is calculated: the "uns" prefix stands for
_unstructured_, and indicates that the face tangential gradient is calculated
using a face-Gauss approach rather than from the cell-centred gradient. This
is more accurate on unstructured meshes, at extra cost per iteration.

The stress is computed by the run-time selectable mechanical law.

---

## User Guide

### What this model solves

`unsLinGeomSolid`:

- solves the linear momentum equation in linear-geometry form;
- uses total displacement `D` as the primary unknown;
- stores the stress and the displacement gradient as _surface_ fields
  (`sigmaf` and `grad(D)f`), evaluated directly at the faces;
- solves the momentum equation with a segregated implicit algorithm only;
- uses the selected mechanical law to calculate stress.

The model is suitable when changes in geometry are small enough that the mesh
can be treated as undeformed, and when the mesh is unstructured enough that
the face-Gauss gradient is worth its extra cost.

### Supported solution algorithms

This model implements a single, segregated implicit algorithm. It does not
read `solutionAlgorithm`, and there is no PETSc SNES or explicit path. Use
`linGeomTotalDispSolid` if you need those.

### Model options

`unsLinGeomSolid` adds no options of its own. The relevant entries in
`unsLinearGeometryCoeffs` are the inherited `solidModel` ones:

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

The under-relaxation factor for `D` itself is set in `fvSolution` under
`relaxationFactors/equations`, as for the other segregated solid models.

```note
`dampingCoeff` is read by the base class but is not applied by this model's
momentum equation.
```

### Recommended dictionary setup

Minimal example for `constant/solidProperties`:

```text
solidModel        unsLinearGeometry;

unsLinearGeometryCoeffs
{
    nCorrectors          10000;
    solutionTolerance    1e-06;
    alternativeTolerance 1e-07;
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

Traction boundaries are handled through `tractionBoundarySnGrad()`, which
returns the surface-normal gradient consistent with the face stress `sigmaf`
and the face implicit stiffness. All the standard solids4foam displacement and
traction patch types apply.

### Field glossary

- `D`: total displacement field at cell centres.
- `pointD`, `pointDD`: total and incremental displacement at mesh points.
- `grad(D)`: cell-centre gradient of `D`.
- `grad(D)f`: face gradient of `D`, computed with the face-Gauss approach.
- `sigma`: cell-centre symmetric stress tensor.
- `sigmaf`: face symmetric stress tensor; this is the field the momentum
  equation actually uses.
- `U`: velocity field, computed as `ddt(D)`.

`sigmaf` is written with `AUTO_WRITE`, so it appears in the time directories.

---

## Developer Notes

### Class role

`unsLinGeomSolid` is the face-gradient counterpart of
`linGeomTotalDispSolid`. The key design choices are:

- `D` is the primary solution variable (`solutionD()` returns `D()`);
- `nonLinGeom()` returns `LINEAR_GEOMETRY`, so the geometry is fixed;
- the stress used in the momentum equation is the _surface_ field `sigmaf_`,
  not an interpolation of the cell-centred `sigma`;
- stress is delegated to `mechanicalModel`.

The class inherits directly from `solidModel`.

### Construction

The constructor:

1. calls the base `solidModel` constructor and `DisRequired()`;
2. creates `sigmaf_` (`AUTO_WRITE`) and `gradDf_` (`NO_WRITE`);
3. takes `impK_`, `impKf_` from the mechanical law and stores `rImpK_`, the
   reciprocal used on every traction boundary evaluation;
4. calls `fvm::d2dt2(D())` purely to force creation of the old-time fields;
5. evaluates the boundary conditions and the gradient once, so that a restart
   starts from a consistent `grad(D)f`;
6. stores the old times of `gradDf_` and `sigmaf_`.

### `evolve()`

A single outer loop:

- store `D` for under-relaxation and residual calculation;
- assemble the momentum equation

  ```text
  rho*d2dt2(D) == laplacian(impKf, D) - fvc::laplacian(impKf, D)
                + fvc::div(Sf & sigmaf) + rho*g
  ```

  i.e. an implicit Laplacian with a deferred correction that recovers the
  original divergence of stress at convergence;
- relax the linear system, apply `solidModel::setCellDisps()`, and solve;
- under-relax `D` with `relaxField()`;
- interpolate `D` to `pointD` and update `grad(D)` and `grad(D)f` through
  `mechanical().grad()`;
- update `sigmaf` and `sigma` from the mechanical law.

The loop exits through the base-class `converged()` helper, or when `iCorr`
reaches `nCorrectors`. Afterwards `pointDD` and `U` are updated once.

Note that `DD` and `grad(DD)` are deliberately not updated inside the loop;
the corresponding lines are commented out in the source.

### Traction boundary treatment

`tractionBoundarySnGrad()` returns

```text
((traction - n*pressure) - (n & (sigma - impK*gradD)))*rImpK
```

evaluated with the _face_ stress and gradient. The `impK*gradD` term removes
the part of the traction that the implicit Laplacian already accounts for, so
that the boundary condition is consistent with the deferred-correction split
used in `evolve()`.

### Extension points

- `evolve()` if a different outer-iteration strategy is needed;
- `tractionBoundarySnGrad()` if the deferred-correction split changes.

Because the class holds face fields as its primary state, any new algorithm
must keep `sigmaf_` and `gradDf_` up to date, not just their cell-centred
counterparts.

---

## Tutorials

No tutorial selects `unsLinearGeometry` by default. It is offered as a
commented alternative in several linear-elasticity and contact cases, which
are the natural places to try it:

- `solids/linearElasticity/patchTest`
- `solids/linearElasticity/plateHole`
- `solids/linearElasticity/ellipticPlate`
- `solids/linearElasticity/narrowTmember`
- `solids/linearElasticity/punch`
- `solids/linearElasticity/flatEndedRigidIndenter`
- `solids/linearElasticity/rigidCylinderContactBrick`
- `solids/linearElasticity/slidingFrictionBall`
- `solids/multiMaterial/layeredPipe`
- `solids/abaqusUMATs/plateHoleTotalDispUMAT`

Swap the `solidModel` entry in `constant/solidProperties` to switch a case
over; no other case settings need to change.
