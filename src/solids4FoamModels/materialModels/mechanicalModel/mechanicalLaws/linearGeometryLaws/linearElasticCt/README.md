---
sort: 12
---

# linearElasticCt

This page documents the small-strain elastic law that assigns a cell-wise
Young's modulus from CT image values. The runtime type is:

```text
linearElasticCt
```

```warning
`linearElasticCt.C` is listed in neither `Make/files.openfoam` nor
`Make/files.foamextend`, so it is **not compiled** into the standard
`solids4FoamModels` library. Its run-time registration therefore never
happens, and selecting `linearElasticCt` in `mechanicalProperties` fails with
an unknown-type error. The source is present and this page describes it, but
the law has to be added to a build list before it can be used. This is tracked
as [issue #336](https://github.com/solids4foam/solids4foam/issues/336); note
that [issue #335](https://github.com/solids4foam/solids4foam/issues/335)
applies to this law too once it is built.
```

---

## User Guide

### What it computes

The law converts each selected CT value, `Hu`, to a relative density and then
to Young's modulus:

```text
relRho = 0.000769*Hu + 1.028                 for Hu > 816
relRho = max(0.0019*Hu + 0.105, 0.2)        otherwise

E = 0.09*relRho^7.4*1e9                     for relRho > 1.54
E = (0.06 + 0.9*relRho^2)*1e9               otherwise
```

Cells whose selected slice, row or column is on an outer CT-array index are
instead assigned `Hu = 0`, `relRho = 0` and `E = 0`. The lookup uses one CT
sample without interpolation. The row and column indices are reversed after
they are found.

The elastic coefficients and stress are then evaluated as implemented:

```text
mu     = E/(1 + nu)
lambda = E*nu/((1 + nu)*(1 - 2*nu))         // plane strain / 3-D
lambda = E*nu/((1 + nu)*(1 - nu))           // planeStress yes

epsilon = symm(grad(D))
sigma   = 2*mu*epsilon + lambda*tr(grad(D))*I
```

The law implements `correct(volSymmTensorField&)` and
`correct(surfaceSymmTensorField&)`. Both abort with a fatal error for an
incremental solid solver. The inherited point-tensor overload aborts with
`notImplemented`. On non-foam-extend forks, the inherited
`CompactListList` quadrature-point overload also aborts with
`notImplemented`.

### Model options

| Entry | Required | Description |
| --- | --- | --- |
| `rho` | yes | Density, `[1 -3 0 0 0 0 0]` |
| `nu` | yes | Poisson's ratio, dimensionless |
| `constantE` | no | Use a uniform `E`; default `false` |
| `value` | if `constantE` | Uniform `E`, read as a scalar with pressure units |
| `nSlices` | for CT | Number of CT slices |
| `nRows` | for CT | Number of rows per slice |
| `nColumns` | for CT | Number of columns per row |
| `sliceSpacing` | for CT | Slice spacing in mesh coordinate units |
| `rowSpacing` | for CT | Row spacing in mesh coordinate units |
| `columnSpacing` | for CT | Column spacing in mesh coordinate units |
| `sliceOffset` | for CT | Offset in the mesh z direction |
| `rowOffset` | for CT | Offset in the mesh y direction |
| `columnOffset` | for CT | Offset in the mesh x direction |
| `ctImagesFilePath` | for CT | Parent path of the `ctImages` directory |
| `useRotationMatrix` | no | Rotate lookup coordinates; default `false` |
| `vectorBeforeRotation` | if rotating | First direction vector |
| `vectorAfterRotation` | if rotating | Second direction vector |
| `centreOfRotation` | if rotating | Rotation centre in mesh coordinates |
| `smooth` | no | Its presence enables smoothing |
| `weight` | if smoothing | Weight retained from the current cell |
| `nSmoothIters` | if smoothing | Number of smoothing iterations |
| `bound` | no | Its presence enables a lower bound on `E` |
| `boundLowerValue` | if bounding | Lower `E`, `[1 -1 -2 0 0 0 0]` |
| `solvePressureEqn` | no | Base switch; default `false` |
| `pressureSmoothingScaleFactor` | no | Base scalar; default `100` |
| `regionName` | no | Base solid-mesh region name |

`value`, the spacings and the offsets are read as plain scalars, so the code
does not check their dimensions. The spacings and offsets are compared
directly with the mesh cell-centre coordinates. The image files must be named
`ct_scan_0.txt` through `ct_scan_<nSlices - 1>.txt` under
`<ctImagesFilePath>/ctImages`. Each file is read in row-major order with
`nRows*nColumns` scalar values.

The presence of `smooth` or `bound` activates that operation; the value of
either entry is not read. Smoothing is in-place and uses inverse-distance
weighted neighbouring-cell values. `weight` blends the current value with
the neighbour average. Bounding is applied after smoothing.

The base class reads `solvePressureEqn` and
`pressureSmoothingScaleFactor`, but this law does not call
`updateSigmaHyd()`, so they do not change its stress calculation.
`regionName` selects the solid model queried to determine whether the solver
is incremental. If omitted, the base class searches for `solid` and then
`region0`.

### Recommended dictionary setup

This example describes a femur CT volume with mesh coordinates in metres:

```text
planeStress     no;

mechanical
(
    femur
    {
        type            linearElasticCt;
        rho             rho [1 -3 0 0 0 0 0] 1900;
        nu              nu [0 0 0 0 0 0 0] 0.3;

        constantE       no;
        nSlices         240;
        nRows           512;
        nColumns        512;
        sliceSpacing    0.001;
        rowSpacing      0.0005;
        columnSpacing   0.0005;
        sliceOffset     0;
        rowOffset       -0.128;
        columnOffset    -0.128;
        ctImagesFilePath constant;

        // Optional coordinate transformation
        // useRotationMatrix    yes;
        // vectorBeforeRotation (0 0 1);
        // vectorAfterRotation  (0 1 0);
        // centreOfRotation     (0 0 0);

        // Optional smoothing
        // smooth          yes;
        // weight          0.5;
        // nSmoothIters    5;

        // Optional lower bound
        // bound           yes;
        // boundLowerValue Emin [1 -1 -2 0 0 0 0] 1e6;
    }
);
```

For a uniform pressure-valued modulus instead of CT data, set `constantE yes`
and supply a scalar `value`, for example `value 17e9;`. The CT size, spacing,
offset and path entries are then not read.

### Field glossary

- `E`: cell-wise Young's modulus. It has zero-gradient boundaries and is
  explicitly written during construction.
- `mu`: cell-wise coefficient `E/(1 + nu)` with its face interpolation
  maintained in `muf_`.
- `lambda`: cell-wise first Lame coefficient with its face interpolation
  maintained in `lambdaf_`.
- `Hu`: temporary cell-wise selected CT value, created only while `E` is
  assembled from CT data.
- `relRho`: temporary cell-wise relative density derived from `Hu`, created
  only while `E` is assembled from CT data.

---

## Developer Notes

### Class role

`linearElasticCt` derives directly from `mechanicalLaw` and is registered in
the linear-geometry mechanical-law selection table. It stores `E_`, `mu_` and
`lambda_` as cell fields, `muf_` and `lambdaf_` as face fields, and `nu_` as a
`dimensionedScalar`. Coordinate transformation is controlled by
`useRotationMatrix_`, `rotationMatrix_` and `centreOfRotation_`.

The header and source use `OPENFOAM_NOT_EXTEND` guards to select field-access
APIs for OpenFOAM.com and OpenFOAM.org versus foam-extend. The inherited
quadrature-point interface is guarded by `#ifndef FOAMEXTEND`. The source is
currently listed in neither `Make/files.openfoam` nor
`Make/files.foamextend`, so the law is not compiled into the standard library
for any fork.

### Construction

The base constructor first reads `solvePressureEqn` and
`pressureSmoothingScaleFactor`, selects `regionName`, and aborts if no solid
region can be found. The derived constructor then creates zero-valued `E`,
reads `nu`, creates the cell and face elastic-coefficient fields, and reads
`useRotationMatrix`.

When rotation is enabled, both direction vectors are read, normalised and used
to construct the rotation tensor that maps the post-rotation direction back
to the pre-rotation direction. A zero-magnitude direction causes a fatal
error. The centre of rotation is then read.

`setYoungsModulusFromCt()` next either assigns the uniform `value`, or reads
the CT geometry and files, maps cell centres to image indices, evaluates the
density and modulus correlations, optionally smooths and bounds `E`, and
writes `E`. Failure to open an image file prints `Cannot open <file>.` but
does not abort. Finally, the constructor recalculates all cell and face
coefficients and substitutes the plane-stress expression for `lambda` when
requested. It does not validate the range of `nu`.

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
- `setYoungsModulusFromCt()`: builds `E` from a uniform value or the CT lookup,
  then performs the optional smoothing and lower-bound operations.

### Extension points

A related image-driven law can copy this class and replace the CT-to-density
and density-to-property relations in `setYoungsModulusFromCt()`. A different
image layout requires changing the file-reading and coordinate-index mapping.
New interfaces should override the corresponding base `correct` overload or
`materialTangent()` rather than inheriting its `notImplemented` behaviour.
The new source must also be added to the appropriate fork build lists and its
runtime registration must remain in the linear-geometry selection table.

The source is at
[linearElasticCt.C](https://github.com/solids4foam/solids4foam/blob/master/src/solids4FoamModels/materialModels/mechanicalModel/mechanicalLaws/linearGeometryLaws/linearElasticCt/linearElasticCt.C).

---

## Tutorials

No tutorial currently uses `linearElasticCt`.
