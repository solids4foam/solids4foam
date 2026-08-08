---
sort: 20
---

# OgdenElastic

This page documents the three-term Ogden hyperelastic mechanical law. The
runtime type is:

```text
OgdenElastic
```

---

## User Guide

### What it computes

The cell-centred `correct` forms the right Cauchy-Green tensor and its
principal stretches. With `c_i` denoting an eigenvalue of `C`, the implemented
relation is:

```text
C        = F^T F
J        = det(F)
lambda_i = max(sqrt(c_i), VSMALL)
p_i      = mu1*lambda_i^alpha1
         + mu2*lambda_i^alpha2
         + mu3*lambda_i^alpha3
s        = rotation of diag(p_1, p_2, p_3) to the spatial basis
sigmaHyd = (K/2)*(J^2 - 1)
sigma    = (1/J)*[dev(s - (mu1 + mu2 + mu3)*I)
           + sigmaHyd*I + symm(F sigma0 F^T)]
```

The subtraction of `(mu1 + mu2 + mu3)*I` occurs inside `dev`, so it does not
alter the deviatoric result. `updateSigmaHyd()` either assigns the expression
above directly or obtains a smoothed field by solving the pressure equation.
If the solid model enforces material linearity, `updateF()` returns a Hookean
stress based on `mu1 + mu2 + mu3` and `K` before the Ogden calculation.

The law implements `correct(volSymmTensorField&)`. Its
`correct(surfaceSymmTensorField&)` overload aborts with `notImplemented`. The
inherited point-tensor overload also aborts. On OpenFOAM.com and OpenFOAM.org,
the inherited `CompactListList` quadrature-point overload is present and
aborts; that interface is not compiled on foam-extend.

### Model options

| Entry | Required | Description |
| --- | --- | --- |
| `mu1` | yes | First modulus parameter, `[1 -1 -2 0 0 0 0]` |
| `mu2` | yes | Second modulus parameter, `[1 -1 -2 0 0 0 0]` |
| `mu3` | yes | Third modulus parameter, `[1 -1 -2 0 0 0 0]` |
| `alpha1` | yes | First exponent, `[0 0 0 0 0 0 0]` |
| `alpha2` | yes | Second exponent, `[0 0 0 0 0 0 0]` |
| `alpha3` | yes | Third exponent, `[0 0 0 0 0 0 0]` |
| `K` | yes | Bulk modulus, `[1 -1 -2 0 0 0 0]` |
| `rho` | yes | Density, `[1 -3 0 0 0 0 0]` |
| `solvePressureEqn` | no | Solve for pressure; default `no` |
| `pressureSmoothingScaleFactor` | no | Pressure scale; default `100` |
| `regionName` | no | Base mesh region; otherwise detected |

All seven law-specific entries are read as `dimensionedScalar` values. The
constructor does not apply explicit range or sign checks. `rho` is read lazily
by the base-class density accessor. The remaining base-class entries are
described on the
[mechanical-law page](https://www.solids4foam.com/documentation/material-models/mechanical-law.html).

Initial stress is not a dictionary entry for this law. The base class instead
reads a `sigma0` volume field if it is present and otherwise creates a zero
field.

### Recommended dictionary setup

The following rubber parameters are used by the commented Ogden alternative
in the `cylinderCrush` tutorial:

```text
mechanical
(
    rubber
    {
        type            OgdenElastic;
        rho             rho [1 -3 0 0 0 0 0] 1000;
        K               K [1 -1 -2 0 0 0 0] 1.410e9;
        mu1             mu1 [1 -1 -2 0 0 0 0] 0.746e6;
        mu2             mu2 [1 -1 -2 0 0 0 0] -0.306e6;
        mu3             mu3 [1 -1 -2 0 0 0 0] 6.609e-5;
        alpha1          alpha1 [0 0 0 0 0 0 0] 1.748;
        alpha2          alpha2 [0 0 0 0 0 0 0] -1.656;
        alpha3          alpha3 [0 0 0 0 0 0 0] 7.671;

        // Optional
        // solvePressureEqn             no;
        // pressureSmoothingScaleFactor 100;
        // regionName                   region0;
    }
);
```

### Field glossary

- `mu_<name>`, `K_<name>`: uniform cell fields set to `mu1 + mu2 + mu3` and
  `K` for use by the base class.
- `muf_<name>`, `Kf_<name>`: corresponding interpolated face fields.
- `F_`, `relF_`: total and relative deformation gradients maintained by
  `updateF()`; `F_` is written for incremental total-Lagrangian restarts.
- `sigmaHyd`, `grad(sigmaHyd)`: hydrostatic stress quantity and its gradient,
  created when the stress is first corrected.
- `sigma0`: initial stress field, read if present and otherwise initialised to
  zero. It is pushed forward as `symm(F sigma0 F^T)/J`.
- `impK`: implicit stiffness field returned as `(4/3)*mu + K`.
- `eigenVal(C)`, `eigenVec(C)`, `s`: non-writing temporary fields created
  during each cell-centred stress correction.

---

## Developer Notes

### Class role

`OgdenElastic` derives directly from `mechanicalLaw` and is registered in the
`nonLinGeomMechLaw` runtime selection table. It stores seven uniform
`dimensionedScalar` members: `mu1_`, `mu2_`, `mu3_`, `alpha1_`, `alpha2_`,
`alpha3_` and `K_`.

The `#ifdef FOAMEXTEND` branches in `correct()` select the appropriate mutable
internal- and boundary-field APIs. The quadrature-point base interface is
guarded by `#ifndef FOAMEXTEND`. The law is listed in both `Make/files.openfoam`
and `Make/files.foamextend`.

### Construction

The base constructor first adds `solvePressureEqn` with default `false` and
`pressureSmoothingScaleFactor` with default `100.0` when they are absent. It
then uses `regionName`, or detects `solid` and then `region0`. Failure to find a
base mesh region produces a fatal error.

The derived constructor reads `mu1`, `mu2`, `mu3`, `alpha1`, `alpha2`,
`alpha3` and `K`, in that order, and prints their values. A missing or malformed
required entry produces the dictionary reader's fatal error; the class has no
additional parameter validation. It sets the inherited cell shear modulus to
`mu1 + mu2 + mu3`, sets the cell bulk modulus to `K`, and interpolates both to
the faces. Construction itself emits no law-specific warnings.

When `sigma0` is later created on a distinct single-material mesh, the base
class warns that mapping from the base mesh may be incorrect. `updateF()` also
aborts for a non-incremental updated-Lagrangian solid and warns once per time
step when material linearity is enforced.

### Key methods

- `impK()`: returns `(4/3)*(mu1 + mu2 + mu3) + K`. This is the diffusivity of
  the solid model's implicit Laplacian term and affects outer-iteration
  convergence rate rather than the converged answer.
- `correct(volSymmTensorField&)`: updates `F`, computes the eigenpairs of
  `F^T F`, builds and rotates the three principal stress values, updates the
  hydrostatic quantity, and adds the pushed-forward initial stress.
- `correct(surfaceSymmTensorField&)`: aborts with `notImplemented`.
- `materialTangent()`: is not overridden. The inherited implementation aborts
  with `notImplemented`, so this law does not provide a 6x6 Newton tangent.

### Extension points

A related principal-stretch law can copy this class, retain registration in
the nonlinear-geometry selection table, and replace the principal-stress and
hydrostatic expressions in the cell-centred `correct()`. Implement the face,
point or quadrature-point `correct()` overloads for solid models that use those
interfaces, and override `materialTangent()` for Newton methods that request a
consistent tangent. Add the new source file to each supported fork's build
list.

The source is at
[OgdenElastic.C](https://github.com/solids4foam/solids4foam/blob/master/src/solids4FoamModels/materialModels/mechanicalModel/mechanicalLaws/nonLinearGeometryLaws/OgdenElastic/OgdenElastic.C).

---

## Tutorials

- `solids/hyperelasticity/cylinderCrush` contains an equivalent Ogden setup,
  but it is commented out; the active configuration uses
  `MooneyRivlinElastic`.
