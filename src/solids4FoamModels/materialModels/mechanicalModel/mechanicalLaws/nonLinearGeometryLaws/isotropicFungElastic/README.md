---
sort: 22
---

# isotropicFungElastic

This page documents the compressible isotropic Fung-like hyperelastic law. The
runtime type is:

```text
isotropicFungElastic
```

---

## User Guide

### What it computes

The law forms an isochoric left Cauchy-Green tensor and its first invariant,
then updates the Cauchy stress as follows:

```text
J       = det(F)
isoB    = J^(-2/3)*symm(F & F.T)
isoI1   = tr(isoB)
psi1    = c1*exp(0.5*c2*(isoI1 - 3))
sigmaHydExplicit = 0.5*K*(J^2 - 1)
sigma   = (1/J)*(psi1*dev(isoB) + sigmaHyd*I
          + symm(F & sigma0 & F.T))
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

`c1` and `c2` are required. The bulk response is specified by either `K` or
`nu`; when both are present, `K` takes precedence.

| Entry | Required | Description |
| --- | --- | --- |
| `c1` | yes | Shear parameter, `[1 -1 -2 0 0 0 0]` |
| `c2` | yes | Exponential parameter, `[0 0 0 0 0 0 0]` |
| `K` | `K` or `nu` | Bulk modulus, `[1 -1 -2 0 0 0 0]` |
| `nu` | `K` or `nu` | Poisson's ratio, `[0 0 0 0 0 0 0]` |
| `rho` | yes | Density, `[1 -3 0 0 0 0 0]` |
| `solvePressureEqn` | no | Solve for hydrostatic stress; default `no` |
| `pressureSmoothingScaleFactor` | no | Pressure smoothing; default `100` |
| `regionName` | no | Base mesh region name |

The `c1` dictionary value supplies the uniform fallback for a `c1` volume
field. A `c1` field in the current time directory is read instead when
present. The field is written automatically.

When `nu` is used, the constructor sets `mu = c1`, forms `E = 3*c1`, and
computes:

```text
K = E/(3*(1 - 2*nu)) = c1/(1 - 2*nu)
```

The constructor does not validate the signs or ranges of `c1`, `c2`, `K` or
`nu`. The inherited `regionName` selection first uses an explicit entry, then
looks for `solid`, and finally `region0`; construction aborts if none exists.

### Recommended dictionary setup

The following example represents a generic nearly incompressible soft tissue:

```text
mechanical
(
    softTissue
    {
        type            isotropicFungElastic;
        rho             rho [1 -3 0 0 0 0 0] 1000;
        c1              c1 [1 -1 -2 0 0 0 0] 10e3;
        c2              c2 [0 0 0 0 0 0 0] 5;
        nu              nu [0 0 0 0 0 0 0] 0.49;

        // Alternative to nu
        // K             K [1 -1 -2 0 0 0 0] 1e6;

        // Optional
        // solvePressureEqn no;
        // pressureSmoothingScaleFactor 100;
        // regionName    region0;
    }
);
```

### Field glossary

- `c1`: cell-centred shear-parameter field, read if present and written
  automatically.
- `c1f`: face interpolation of `c1`, formed during construction.
- `mu`, `muf`: linearised cell and face shear-modulus fields; `mu` is set to
  `c1` and `muf` is interpolated from it.
- `K`, `Kf`: cell and face bulk-modulus fields; `Kf` is interpolated from `K`.
- `F`, `Ff`: total deformation gradients at cell centres and faces, created by
  the base class and written for restart.
- `relF`, `relFf`: relative deformation gradients maintained by the base
  class.
- `sigmaHyd`: cell-centred hydrostatic stress, written automatically.
- `sigma0`: initial stress read from the field file when present; the law
  pushes it forward as `symm(F & sigma0 & F.T)/J`.
- `impK`: implicit stiffness field, read if present and not written.

---

## Developer Notes

### Class role

`isotropicFungElastic` derives directly from `mechanicalLaw` and is registered
in the `nonLinGeomMechLaw` runtime-selection table. It stores `c1_` as a
`volScalarField`, `c1f_` as a `surfaceScalarField`, and `c2_` as a constant
`dimensionedScalar`.

The class has no fork-specific preprocessor guards. Its source appears in both
`Make/files.openfoam` and `Make/files.foamextend`. The inherited
quadrature-point interface is guarded by `#ifndef FOAMEXTEND` in
`mechanicalLaw`.

### Construction

The base constructor first reads or inserts `solvePressureEqn` and
`pressureSmoothingScaleFactor`, then selects the base mesh region. The derived
constructor proceeds in this order:

1. It constructs `c1_` from a field file when present, using the required
   dictionary value as its fallback, interpolates `c1f_`, and reads `c2_`.
2. It sets the base-class linearised shear modulus to `c1_`.
3. It reads `K` when present. Otherwise it reads `nu`, forms `E = 3*c1_`, and
   sets `K = E/(3*(1 - 2*nu))`.
4. It prints the maxima of `c1`, `mu` and `K`, together with the value of `c2`.
5. It interpolates `mu` and `K` to their face fields and forces `sigma0` to be
   created or read.

If neither `K` nor `nu` is present, construction terminates with a fatal error
stating that one must be specified. The base constructor terminates with a
fatal error if it cannot identify the solid region. Creating `sigma0` for a
single material whose material mesh differs from the base mesh emits a warning
that the field may not be correct. There are no law-specific parameter-range
checks or warnings.

### Key methods

- `impK()`: returns the field `(4/3)*mu + K`. This is the diffusivity of the
  solid model's implicit Laplacian term and affects outer-iteration convergence
  rate rather than the converged answer.
- `correct(volSymmTensorField&)`: updates `F`, calculates the isochoric tensor
  and exponential response, passes the volumetric stress to
  `updateSigmaHyd()`, pushes forward `sigma0`, and assigns the Cauchy stress.
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

A related invariant-based law can copy this class, replace `psi1` and the
stress expressions in both `correct()` overloads, and add any required cell
properties and face interpolations. It should retain registration in the
`nonLinGeomMechLaw` table and entries in both build lists. Implement the point,
quadrature-point and material-tangent interfaces if the target solid models
require them.

The source is at
[isotropicFungElastic.C](https://github.com/solids4foam/solids4foam/blob/master/src/solids4FoamModels/materialModels/mechanicalModel/mechanicalLaws/nonLinearGeometryLaws/isotropicFungElastic/isotropicFungElastic.C).

---

## Tutorials

No tutorial currently uses `isotropicFungElastic`.
