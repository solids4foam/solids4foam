---
sort: 3
---

# linearElastic

This page documents the isotropic Hookean linear elastic mechanical law. The
runtime type is:

```text
linearElastic
```

This is the default choice for small-strain elasticity and by far the most
widely used law in the solids4foam tutorial suite. It is also the law most
often nested inside the `poroMechanicalLaw` and `thermoMechanicalLaw`
wrappers.

---

## User Guide

### What it computes

The cell-centred `correct` splits the stress into deviatoric and hydrostatic
parts:

```text
sigmaHyd = K*tr(epsilon)
sigma    = 2*mu*dev(epsilon) + sigmaHyd*I + sigma0
```

The hydrostatic part is passed through `updateSigmaHyd()`, so when
`solvePressureEqn` is enabled it is obtained from a smoothed Laplacian
equation instead of directly from `K*tr(epsilon)`. The face-centred and
point-centred overloads instead use the plain form

```text
sigma = 2*mu*epsilon + lambda*tr(epsilon)*I + sigma0f
```

and both abort with `notImplemented` if `solvePressureEqn` is enabled. The
point overload additionally aborts if `sigma0` is non-zero.

### Model options

Elastic constants are given as **either** `E` and `nu`, **or** `mu` and `K`.
The `E`/`nu` pair is checked first; if neither pair is present the law aborts
at construction.

| Entry | Required | Description |
| --- | --- | --- |
| `E` | with `nu` | Young's modulus, `[1 -1 -2 0 0 0 0]` |
| `nu` | with `E` | Poisson's ratio, dimensionless |
| `mu` | with `K` | Shear modulus, `[1 -1 -2 0 0 0 0]` |
| `K` | with `mu` | Bulk modulus, `[1 -1 -2 0 0 0 0]` |
| `sigma0` | no | Uniform initial stress, symmetric tensor |
| `rho` | yes | Density, `[1 -3 0 0 0 0 0]` |

`sigma0` is optional. If it is absent but a non-uniform `sigma0` field already
carries a non-zero value, that field is used instead and a message is printed.

The remaining `mechanicalLaw` base-class entries — `solvePressureEqn`
(default `no`), `pressureSmoothingScaleFactor` (default `100`) and
`regionName` — apply here as described in the
[subsection index](https://www.solids4foam.com/documentation/material-models/linear-geometry-laws.html).

### Derived quantities and checks

The constructor derives whichever pair was not supplied, and also the first
Lame parameter:

```text
mu     = E/(2*(1 + nu))
lambda = nu*E/((1 + nu)*(1 - 2*nu))        // plane strain / 3-D
lambda = nu*E/((1 + nu)*(1 - nu))          // planeStress yes
```

Three checks are applied:

- `nu` outside `[-1, 0.5]` is a fatal error;
- `nu` above `0.49` with `solvePressureEqn no` produces a warning suggesting
  the pressure equation be enabled;
- `planeStress yes` combined with `solvePressureEqn yes` is a fatal error.

When `nu` is exactly `0.5`, `K` and `lambda` are set to `GREAT` and `impK`
falls back to `2*mu`; use a hybrid or pressure-solving solid model in that
case.

### Recommended dictionary setup

```text
planeStress     no;

mechanical
(
    steel
    {
        type            linearElastic;
        rho             rho [1 -3 0 0 0 0 0] 7854;
        E               E [1 -1 -2 0 0 0 0] 200e9;
        nu              nu [0 0 0 0 0 0 0] 0.3;

        // Optional
        // sigma0        sigma0 [1 -1 -2 0 0 0 0] (0 0 0 0 0 0);
        // solvePressureEqn no;
    }
);
```

### Field glossary

- `epsilon`, `epsilonf`: small-strain tensor at cell centres and faces; both
  have their old-time values stored so that incremental solid models work.
- `sigmaHyd`: hydrostatic stress, `tr(sigma)/3`, maintained by the base class.
- `sigma0`, `sigma0f`: initial (residual) stress at cells and faces.
- `impK`: implicit stiffness handed to the solid model's Laplacian term;
  `2*mu + lambda`, or `2*mu` for `nu == 0.5`.

---

## Developer Notes

### Class role

`linearElastic` derives directly from `mechanicalLaw` and stores five
`dimensionedScalar` members: `mu_`, `K_`, `E_`, `nu_` and `lambda_`. Because
the properties are uniform scalars rather than fields, the law is cheap and is
the natural base case against which other laws are compared.

It implements an unusually complete set of interfaces, which is why it is the
only law usable by every solid model in the toolbox:

- `correct(volSymmTensorField&)` and the `epsilon`-taking overload;
- `correct(surfaceSymmTensorField&)` and its `epsilon`-taking overload, used
  by the `uns` (cell-centred, face-based) solid models;
- `correct(pointSymmTensorField&, const pointTensorField&)`, used by the
  vertex-centred solid models;
- `correct(CompactListList<symmTensor>&, const CompactListList<tensor>&)`, the
  face quadrature-point form used by the higher-order vertex-centred path.

The quadrature-point overload is guarded by `#ifndef FOAMEXTEND` in both the
header and the source, so it is compiled only for the OpenFOAM.com and
OpenFOAM.org forks. Everything else is common to all three forks; the law is
listed in both `Make/files.openfoam` and `Make/files.foamextend`.

### Construction

The constructor stores old-time `epsilon` and `epsilonf`, reads the elastic
pair, derives the remaining constants, validates `nu`, and finally reads
`sigma0`. Note that the error messages emitted for a missing elastic pair name
`linearElasticMisesPlastic` rather than `linearElastic`, a copy-paste artefact
that does not affect behaviour.

### Key methods

- `impK()`: the implicit stiffness, `2*mu + lambda` (or `2*mu` for
  `nu == 0.5`). This is the diffusivity of the Laplacian term in the solid
  model's momentum equation and controls the outer-iteration convergence rate
  rather than the converged answer.
- `materialTangent()`: returns the full 6x6 isotropic elasticity matrix, used
  by the block-coupled and PETSc SNES solid models.
- `bulkModulus()`, `shearModulus()`: uniform fields with zero-gradient
  boundaries, provided so that wrapper laws such as `thermoMechanicalLaw` can
  query them generically.
- `mu()`, `K()`, `E()`, `nu()`, `lambda()`: direct scalar accessors, used by
  the coupled solid models.

### Extension points

A new isotropic elastic law is most cheaply written by copying this class and
replacing the two `correct` bodies. If the new law has spatially varying
properties, the `dimensionedScalar` members must become `volScalarField`
members, and the scalar accessors above have to be dropped or reimplemented —
see `linearElasticFromFile` for that pattern.

The source is at
[linearElastic.C](https://github.com/solids4foam/solids4foam/blob/master/src/solids4FoamModels/materialModels/mechanicalModel/mechanicalLaws/linearGeometryLaws/linearElastic/linearElastic.C).

---

## Tutorials

`linearElastic` is used by most of the `solids/linearElasticity` cases, for
example:

- `solids/linearElasticity/plateHole`
- `solids/linearElasticity/cooksMembrane`
- `solids/linearElasticity/pressurisedCylinder`
- `solids/linearElasticity/contactPatchTest`
- `solids/beamsPlatesShells/squarePlate`
- `solids/multiMaterial/layeredPipe`
- `solids/fracture/crackingPlateHole`
- `fluidSolidInteraction/3dTube`
- `fluidSolidInteraction/perpendicularFlap`

It is also the law nested inside `thermoMechanicalLaw` in every
`solids/thermoelasticity` and `thermoFluidSolidInteraction` tutorial.
