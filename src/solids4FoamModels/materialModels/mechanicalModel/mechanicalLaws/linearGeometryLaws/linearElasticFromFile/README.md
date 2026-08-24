---
sort: 13
---

# linearElasticFromFile

This page documents the small-strain Hookean law that reads a cell-wise
Young's modulus field from disk. The runtime type is:

```text
linearElasticFromFile
```

```warning
The shear term of this law appears to be twice the Hookean value. It sets
`mu = E/(1 + nu)`, whereas the shear modulus is `E/(2*(1 + nu))`, and then
forms the stress as `mu*twoSymm(grad(D))`, which is `2*mu*epsilon`. The
resulting deviatoric stress is therefore `2*E/(1 + nu)*epsilon` rather than
`E/(1 + nu)*epsilon`. For comparison, `linearElastic` sets
`mu = E/(2*(1 + nu))` for the same stress expression. The volumetric term and
the `planeStress` handling are correct. Treat results from this law with
caution until this is resolved; `linearElasticCt` contains the same
expression. This is tracked as
[issue #335](https://github.com/solids4foam/solids4foam/issues/335).
```

---

## User Guide

### What it computes

The law reads the `E` field and evaluates the elastic coefficients as
implemented:

```text
mu     = E/(1 + nu)
lambda = E*nu/((1 + nu)*(1 - 2*nu))        // plane strain / 3-D
lambda = E*nu/((1 + nu)*(1 - nu))          // planeStress yes

epsilon = symm(grad(D))
sigma   = 2*mu*epsilon + lambda*tr(grad(D))*I
```

The face-centred calculation uses interpolated `mu` and `lambda` fields with
`grad(D)f` in the same relation.

The law implements `correct(volSymmTensorField&)` and
`correct(surfaceSymmTensorField&)`. Both abort with a fatal error for an
incremental solid solver. The inherited point-tensor overload aborts with
`notImplemented`. On non-foam-extend forks, the inherited
`CompactListList` quadrature-point overload also aborts with
`notImplemented`.

### Model options

| Entry | Required | Description |
| --- | --- | --- |
| `nu` | yes | Poisson's ratio, dimensionless |
| `rho` | yes | Density, `[1 -3 0 0 0 0 0]` |
| `solvePressureEqn` | no | Base switch; default `false` |
| `pressureSmoothingScaleFactor` | no | Base scalar; default `100` |
| `regionName` | no | Base solid-mesh region name |

Young's modulus is not a dictionary entry. A volume scalar field named `E`
must exist in the current time directory when the law is constructed. It is
read with `MUST_READ` and marked `AUTO_WRITE`. Its dimensions must be pressure
dimensions, `[1 -1 -2 0 0 0 0]`, for the field operations to be dimensionally
consistent.

The base class reads `solvePressureEqn` and
`pressureSmoothingScaleFactor`, but this law does not call
`updateSigmaHyd()`, so they do not change its stress calculation.
`regionName` selects the solid model queried to determine whether the solver
is incremental. If omitted, the base class searches for `solid` and then
`region0`.

The top-level `planeStress` entry selects the alternative expression for
`lambda`. The constructor does not validate the range of `nu`.

### Recommended dictionary setup

This example describes an aluminium component whose spatially varying
Young's modulus is supplied in the `E` volume field:

```text
planeStress     no;

mechanical
(
    gradedAluminium
    {
        type            linearElasticFromFile;
        rho             rho [1 -3 0 0 0 0 0] 2700;
        nu              nu [0 0 0 0 0 0 0] 0.33;

        // Optional base-class entries
        // solvePressureEqn             no;
        // pressureSmoothingScaleFactor 100;
        // regionName                   region0;
    }
);
```

Place the pressure-dimensioned `E` field in the initial time directory before
starting the solver. Its internal and boundary values define the spatial
distribution used by the law.

### Field glossary

- `E`: cell-wise Young's modulus, read from the current time directory and
  marked for automatic writing.
- `mu`: cell-wise coefficient `E/(1 + nu)`.
- `lambda`: cell-wise first Lame coefficient, with the plane-stress form used
  when requested.
- `muf_`, `lambdaf_`: face interpolations of `mu` and `lambda`, calculated
  during construction.
- `grad(D)`, `grad(D)f`: solver-provided displacement-gradient fields looked
  up by the cell-centred and face-centred stress corrections.
- `impK`: temporary cell field containing `2*mu + lambda`.

---

## Developer Notes

### Class role

`linearElasticFromFile` derives directly from `mechanicalLaw` and is
registered in the linear-geometry mechanical-law selection table. It stores
`E_`, `mu_` and `lambda_` as cell fields, `muf_` and `lambdaf_` as face
fields, and `nu_` as a `dimensionedScalar`.

The header uses `#ifdef OPENFOAM_NOT_EXTEND` only to include the surface-field
API for OpenFOAM.com and OpenFOAM.org. The inherited quadrature-point
interface is guarded by `#ifndef FOAMEXTEND`. The source is listed in both
`Make/files.openfoam` and `Make/files.foamextend`.

### Construction

The base constructor first reads `solvePressureEqn` and
`pressureSmoothingScaleFactor`, selects `regionName`, and aborts if no solid
region can be found. The derived constructor then reads `E` from the current
time directory, reads `nu`, calculates `mu` and `lambda`, and interpolates the
two coefficients to faces. If `planeStress` is enabled, it replaces `lambda`
with the plane-stress expression and interpolates that field again.

A missing or unreadable `E` field causes the mandatory field read to abort.
Both `correct` methods emit `Not implemented for incremental solid solver` as
a fatal error when the selected solid model is incremental. The constructor
does not issue a warning or validate `nu`.

### Key methods

- `impK()`: returns the cell-wise implicit stiffness `2*mu + lambda`. This is
  the diffusivity of the solid model's Laplacian term and affects the
  outer-iteration convergence rate rather than the converged answer.
- `correct(volSymmTensorField&)`: looks up `grad(D)` and evaluates the
  cell-centred stress. It rejects incremental solvers.
- `correct(surfaceSymmTensorField&)`: looks up `grad(D)f` and evaluates the
  face-centred stress from the interpolated coefficients. It rejects
  incremental solvers.
- `materialTangent()`: is not overridden; the inherited implementation aborts
  with `notImplemented`.

### Extension points

A related spatially varying elastic law can copy this class and replace the
`E` field read or the coefficient relations. A law that must support
incremental, point-centred or quadrature-point solid models must implement the
corresponding `correct` overloads. Newton-based models also require a
`materialTangent()` override. The new source must remain in the appropriate
fork build lists and its runtime registration must use the intended geometry
selection table.

The source is at
[linearElasticFromFile.C](https://github.com/solids4foam/solids4foam/blob/master/src/solids4FoamModels/materialModels/mechanicalModel/mechanicalLaws/linearGeometryLaws/linearElasticFromFile/linearElasticFromFile.C).

---

## Tutorials

No tutorial currently uses `linearElasticFromFile`.
