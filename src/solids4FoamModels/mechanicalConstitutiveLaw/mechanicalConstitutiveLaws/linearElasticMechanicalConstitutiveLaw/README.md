---
sort: 3
---

# linearElastic

The isotropic Hookean linear elastic law, for small strains. The runtime type
is:

```text
linearElastic
```

It is the most widely used law in the solids4foam tutorials, and the law most
often nested inside `thermoMechanicalLaw` and `poroMechanicalLaw`. The
dictionary layout common to all laws is described in the
[material models page](../../../materialModels/README.md), and the framework
itself in the [framework README](../../README.md).

---

## User Guide

### What it computes

The stress is computed from the displacement gradient supplied by the solid
model, with `epsilon = symm(grad(D))`:

```text
sigma = 2*mu*epsilon + lambda*tr(epsilon)*I + sigma0
```

where `sigma0` is an optional initial stress (see below). The law is
small-strain only: it implements the small-strain update and a nonlinear
geometry solid model, which asks for the finite-strain update, stops with an
error.

The constants are derived from whichever pair of elastic entries was given:

```text
mu     = E/(2*(1 + nu))
lambda = nu*E/((1 + nu)*(1 - 2*nu))        // plane strain / 3-D
lambda = nu*E/((1 + nu)*(1 - nu))          // planeStress yes
K      = lambda + (2/3)*mu
```

If `mu` and `K` are given, `E` and `nu` are first recovered from them with
`E = 9*K*mu/(3*K + mu)` and `nu = (3*K - 2*mu)/(2*(3*K + mu))`.

### Model options

The elastic constants are given as **either** `E` and `nu` **or** `mu` and
`K`. The `E`/`nu` pair is checked first; if neither pair is complete the law
stops at construction.

| Entry | Required | Default | Description |
| --- | --- | --- | --- |
| `rho` | yes | - | Density, `[1 -3 0 0 0 0 0]` |
| `E` | with `nu` | - | Young's modulus, `[1 -1 -2 0 0 0 0]` |
| `nu` | with `E` | - | Poisson's ratio, dimensionless |
| `mu` | with `K` | - | Shear modulus, `[1 -1 -2 0 0 0 0]` |
| `K` | with `mu` | - | Bulk modulus, `[1 -1 -2 0 0 0 0]` |
| `sigma0` | no | zero | Uniform initial stress, `[1 -1 -2 0 0 0 0]` |

The dimensions of `E`, `mu` and `K` (pressure) and of `nu` (dimensionless) are
checked. `nu` must satisfy `-1 < nu < 0.5`: a value at or above `0.5` (to
within `SMALL`) is a fatal error, because the bulk modulus would be
ill-conditioned.

`planeStress` is a top-level entry of `mechanicalProperties`; it changes
`lambda` as shown above and must not be given inside the law's dictionary.

The optional `solvePressureEqn` and `pressureSmoothingScaleFactor` entries are
read by the solid model, not by this law. This law separates its deviatoric
and volumetric responses, so, when selected directly, it can be used with
hydrostatic stress smoothing and with the mixed displacement-pressure
formulations. The `thermoMechanicalLaw` and `poroMechanicalLaw` wrappers do
not provide the split, so neither is available through them.

### Initial stress

An initial (residual) stress can be supplied in two ways:

- a uniform value, as the `sigma0` entry in the law's dictionary; or
- a spatially varying `volSymmTensorField` named `sigma0` in the start time
  directory (or in `0`), for example written by `setFields`.

If both are present, the dictionary value wins and the field is ignored. The
initial stress is added unchanged to every evaluation.

### Tangents and volumetric split

- Scalar tangent: `2*mu + lambda`; the deviatoric scalar tangent is
  `(4/3)*mu`.
- Fourth-order tangent: the analytical isotropic 6x6 elasticity matrix. The
  finite-difference tangent is also available.
- Volumetric split: yes. The volumetric response is `K*tr(epsilon)`, and the
  deviatoric part returned with it keeps any initial stress, including the
  spherical part of `sigma0`.

### State variables

| Name | Type | Role | Description |
| --- | --- | --- | --- |
| `sigma0` | symmTensor | prescribed | Initial stress, read from a field |

`sigma0` is *prescribed* when it comes from a field and *fixed* when the
dictionary gives it. Neither role is history, so this law writes no restart
files and no state fields.

### Example

From the `plateHole` tutorial:

```text
planeStress     no;

mechanical
(
    steel
    {
        type            linearElastic;
        rho             rho [1 -3 0 0 0 0 0] 7854;
        E               E [1 -1 -2 0 0 0 0] 200e+9;
        nu              nu [0 0 0 0 0 0 0] 0.3;
    }
);
```

With a uniform initial stress:

```text
        sigma0          sigma0 [1 -1 -2 0 0 0 0] (10e6 2e6 -3e6 15e6 0 -5e6);
```

### Migrating from the legacy law

- `regionName` is no longer read.
- `nu = 0.5` is now a fatal error; the legacy law accepted it and set `K` and
  `lambda` to `GREAT`. There is no longer a warning for `nu > 0.49`. For a
  nearly incompressible material use a mixed displacement-pressure solid model
  or `solvePressureEqn`.
- `solvePressureEqn` combined with `planeStress yes` is no longer refused.
- One evaluation serves every integration-point location (cells, faces,
  points, quadrature points). The legacy point overload refused a non-zero
  `sigma0`; this law does not.
- `E()`, `nu()`, `mu()` and `lambda()` are still available to solid models:
  `coupledUnsLinGeomLinearElasticSolid` and `kirchhoffPlateSolid` require this
  law for that reason.

---

## Tutorials

`linearElastic` is used by most of the `solids/linearElasticity` cases, for
example:

- [plateHole](../../../../../tutorials/solids/linearElasticity/plateHole/README.md)
- [cantilever2d](../../../../../tutorials/solids/linearElasticity/cantilever2d/README.md),
  whose regression tests cover the dictionary and field forms of `sigma0`
- [cooksMembrane](../../../../../tutorials/solids/linearElasticity/cooksMembrane/README.md)
- [pressurisedCylinder](../../../../../tutorials/solids/linearElasticity/pressurisedCylinder/README.md)
- [contactPatchTest](../../../../../tutorials/solids/linearElasticity/contactPatchTest/README.md)
- [punch](../../../../../tutorials/solids/linearElasticity/punch/README.md)
- [squarePlate](../../../../../tutorials/solids/beamsPlatesShells/squarePlate/README.md)
- [layeredPipe](../../../../../tutorials/solids/multiMaterial/layeredPipe/README.md)
- [3dTube](../../../../../tutorials/fluidSolidInteraction/3dTube/README.md)
- [perpendicularFlap](../../../../../tutorials/fluidSolidInteraction/perpendicularFlap/README.md)

It is also the law nested inside `thermoMechanicalLaw` in the
`solids/thermoelasticity` tutorials, for example
[hotCylinder](../../../../../tutorials/solids/thermoelasticity/hotCylinder/hotCylinder/README.md)
and [hotSphere](../../../../../tutorials/solids/thermoelasticity/hotSphere/README.md).
