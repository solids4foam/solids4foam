---
sort: 19
---

# YeohElastic

This page documents the compressible three-parameter Yeoh hyperelastic law.
The runtime type is:

```text
YeohElastic
```

---

## User Guide

### What it computes

The law uses the following isochoric strain-energy function:

```text
Psi_iso = c1*(isoI1 - 3) + c2*(isoI1 - 3)^2 + c3*(isoI1 - 3)^3
```

The implemented cell-centred stress update is:

```text
J       = det(F)
isoB    = J^(-2/3)*symm(F & F.T)
isoI1   = tr(isoB)
q       = c1 + 2*c2*(isoI1 - 3) + 3*c3*(isoI1 - 3)^2
sigmaHydExplicit = 0.5*K*(J^2 - 1)
J*sigma = 2*q*dev(isoB) + sigmaHyd*I + symm(F & sigma0 & F.T)
```

For the cell-centred overload, `sigmaHyd` is passed through
`updateSigmaHyd()`. It equals `sigmaHydExplicit` unless `solvePressureEqn` is
enabled, in which case the base class obtains it from a smoothed Laplacian
equation. The face-centred overload uses `Ff`, interpolated material fields and
the explicit hydrostatic term directly.

The law implements `correct(volSymmTensorField&)` and
`correct(surfaceSymmTensorField&)`. The inherited point-centred overload aborts
with `notImplemented`. On OpenFOAM.com and OpenFOAM.org, the inherited
`CompactListList` quadrature-point overload also aborts with `notImplemented`;
that overload is excluded from foam-extend by `#ifndef FOAMEXTEND` in the base
class.

### Model options

The three Yeoh coefficients are required. The bulk response is specified by
either `K` or `nu`; when both are present, `K` takes precedence.

| Entry | Required | Description |
| --- | --- | --- |
| `c1` | yes | First Yeoh coefficient, `[1 -1 -2 0 0 0 0]` |
| `c2` | yes | Second Yeoh coefficient, `[1 -1 -2 0 0 0 0]` |
| `c3` | yes | Third Yeoh coefficient, `[1 -1 -2 0 0 0 0]` |
| `K` | `K` or `nu` | Bulk modulus, `[1 -1 -2 0 0 0 0]` |
| `nu` | `K` or `nu` | Poisson's ratio, `[0 0 0 0 0 0 0]` |
| `rho` | yes | Density, `[1 -3 0 0 0 0 0]` |
| `solvePressureEqn` | no | Solve hydrostatic stress; default `false` |
| `pressureSmoothingScaleFactor` | no | Smoothing scale; default `100.0` |
| `regionName` | no | Base mesh region name; detected when absent |

Each coefficient dictionary value supplies the uniform fallback for a volume
field of the same name. A `c1`, `c2` or `c3` field in the current time
directory is read instead when present, and all three fields are written
automatically.

The constructor sets the linearised shear modulus to `mu = 2*c1`. When `nu`
is used, it forms `E = 6*c1` and computes:

```text
K = E/(3*(1 - 2*nu)) = 2*c1/(1 - 2*nu)
```

The constructor does not validate the signs or ranges of `c1`, `c2`, `c3`,
`K` or `nu`. The inherited `regionName` selection first uses an explicit
entry, then looks for `solid`, and finally `region0`; construction aborts if
none exists.

### Recommended dictionary setup

The following example represents a generic nearly incompressible silicone
rubber:

```text
mechanical
(
    siliconeRubber
    {
        type            YeohElastic;
        rho             rho [1 -3 0 0 0 0 0] 1100;
        c1              c1 [1 -1 -2 0 0 0 0] 100e3;
        c2              c2 [1 -1 -2 0 0 0 0] 10e3;
        c3              c3 [1 -1 -2 0 0 0 0] 1e3;
        nu              nu [0 0 0 0 0 0 0] 0.49;

        // Alternative to nu
        // K             K [1 -1 -2 0 0 0 0] 10e6;

        // Optional
        // solvePressureEqn false;
        // pressureSmoothingScaleFactor 100.0;
        // regionName    region0;
    }
);
```

### Field glossary

- `c1`, `c2`, `c3`: cell-centred coefficient fields, read if present and
  written automatically.
- `c1f`, `c2f`, `c3f`: face interpolations of the coefficient fields, formed
  during construction.
- `mu`, `muf`: linearised cell and face shear-modulus fields; `mu = 2*c1`.
- `K`, `Kf`: cell and face bulk-modulus fields; `Kf` is interpolated from `K`.
- `F`, `Ff`: total deformation gradients at cell centres and faces, created by
  the base class and written for restart.
- `relF`, `relFf`: relative deformation gradients maintained by the base
  class.
- `sigmaHyd`: cell-centred hydrostatic stress, written automatically.
- `sigma0`: initial stress read from its field file when present; the law
  pushes it forward as `symm(F & sigma0 & F.T)/J`.
- `impK`: implicit stiffness field, read if present and not written.

---

## Developer Notes

### Class role

`YeohElastic` derives directly from `mechanicalLaw` and is registered in the
`nonLinGeomMechLaw` runtime-selection table. It stores `c1_`, `c2_` and `c3_`
as `volScalarField`s, together with the face fields `c1f_`, `c2f_` and `c3f_`.

The class has no fork-specific preprocessor guards. Its source appears in both
`Make/files.openfoam` and `Make/files.foamextend`. The inherited
quadrature-point interface is guarded by `#ifndef FOAMEXTEND` in
`mechanicalLaw`.

### Construction

The base constructor first reads or inserts `solvePressureEqn` and
`pressureSmoothingScaleFactor`, then selects the base mesh region. The derived
constructor proceeds in this order:

1. It constructs `c1_`, `c2_` and `c3_` from field files when present, using
   the required dictionary values as fallbacks, then interpolates their face
   fields.
2. It sets the base-class linearised shear modulus to `2*c1_`.
3. It reads `K` when present. Otherwise it reads `nu`, forms `E = 6*c1_`, and
   sets `K = E/(3*(1 - 2*nu))`.
4. It prints the maxima of the three coefficients, `mu` and `K`.
5. It interpolates `mu` and `K` to their face fields.

If neither `K` nor `nu` is present, construction terminates with a fatal error
stating that one must be specified. The base constructor terminates with a
fatal error if it cannot identify the solid region. There are no law-specific
parameter-range checks or warnings. When `sigma0` is first created, the base
class warns if a single material uses a material mesh different from the base
mesh because the field may not be correct.

### Key methods

- `impK()`: returns `(4/3)*mu + K`. This is the diffusivity of the solid
  model's implicit Laplacian term and affects outer-iteration convergence rate
  rather than the converged answer.
- `correct(volSymmTensorField&)`: updates `F`, calculates the isochoric tensor
  and Yeoh response, passes the volumetric stress to `updateSigmaHyd()`, pushes
  forward `sigma0`, and assigns the Cauchy stress.
- `correct(surfaceSymmTensorField&)`: performs the corresponding calculation
  on faces, but uses the explicit hydrostatic term instead of
  `updateSigmaHyd()`.
- `materialTangent()`: is not overridden. The inherited implementation aborts
  with `notImplemented`, so this law does not provide a 6x6 Newton tangent.

Both `correct()` implementations first call the base-class `updateF()`. That
method can return early after applying linear elasticity when the solid model
enforces linearity. It also rejects a non-incremental updated-Lagrangian call
with a fatal error.

### Extension points

A related invariant-based law can copy this class, replace the coefficient
`q` and stress expressions in both `correct()` overloads, and add any required
cell properties and face interpolations. It should retain registration in the
`nonLinGeomMechLaw` table and entries in both build lists. Implement the point,
quadrature-point and material-tangent interfaces if the target solid models
require them.

The source is at
[YeohElastic.C](https://github.com/solids4foam/solids4foam/blob/master/src/solids4FoamModels/materialModels/mechanicalModel/mechanicalLaws/nonLinearGeometryLaws/YeohElastic/YeohElastic.C).

---

## Tutorials

No tutorial currently uses `YeohElastic`.
