---
sort: 15
---

# StVenantKirchhoffElastic

This page documents the St. Venant-Kirchhoff hyperelastic law. The runtime
type is:

```text
StVenantKirchhoffElastic
```

The St. Venant-Kirchhoff model is the simplest finite-strain extension of
Hooke's law: the second Piola-Kirchhoff stress is a linear function of the
Green-Lagrange strain. It captures large rotations exactly but is only valid
for small to moderate strains, since the model softens unphysically under
large compression.

---

## User Guide

### Strain energy function

```text
Psi = mu*(E && E) + 0.5*lambda*(tr(E))^2

E = 0.5*(C - I),   C = F.T() & F
```

Differentiating gives the second Piola-Kirchhoff stress that the code
evaluates,

```text
S = 2*mu*E + lambda*tr(E)*I
```

which is then pushed forward to the Cauchy stress,

```text
sigma = (1/J)*F & S & F.T(),   J = det(F)
```

The law is registered in the nonlinear-geometry table, so it works with both
total-Lagrangian and updated-Lagrangian solid models. It is not available to
linear-geometry solid models; use `linearElastic` there instead.

```warning
Because `Psi` is quadratic in the Green strain, the model does not penalise
`J` approaching zero. Do not use it for problems with large volumetric
compression; prefer `neoHookeanElastic` or another of the rubber-like laws.
```

### Model options

Elastic constants must be given as **either** `E` and `nu` **or** `mu` and
`K`. Mixing the two pairs is a fatal error.

| Entry | Requirement | Description |
| --- | --- | --- |
| `E` | with `nu` | Young's modulus, `[1 -1 -2 0 0 0 0]` |
| `nu` | with `E` | Poisson's ratio, dimensionless |
| `mu` | with `K` | Shear modulus, `[1 -1 -2 0 0 0 0]` |
| `K` | with `mu` | Bulk modulus, `[1 -1 -2 0 0 0 0]` |
| `rho` | required | Density, `[1 -3 0 0 0 0 0]` |

The first Lame parameter is derived internally:

```text
lambda = nu*E/((1 + nu)*(1 - 2*nu))          plane strain / 3-D
lambda = nu*E/((1 + nu)*(1 - nu))            plane stress
```

The plane-stress form is used when the top-level `planeStress` entry of
`constant/mechanicalProperties` is `yes`. The bulk modulus is then recomputed
as `K = lambda + (2/3)*mu`, so a user-supplied `K` is overwritten under plane
stress.

The base class also accepts the optional entries `solvePressureEqn` (default
`false`), `pressureSmoothingScaleFactor` (default `100`) and `regionName`, but
this law never calls `updateSigmaHyd()`, so the pressure smoothing has no
effect here.

### Recommended dictionary setup

Minimal example for `constant/mechanicalProperties`:

```text
planeStress     no;

mechanical
(
    steel
    {
        type            StVenantKirchhoffElastic;
        rho             rho [1 -3 0 0 0 0 0] 7800;
        E               E [1 -1 -2 0 0 0 0] 200e+09;
        nu              nu [0 0 0 0 0 0 0] 0.3;
    }
);
```

### Field glossary

- `F`, `Ff`: total deformation gradient at cell centres and faces.
- `relF`, `relFf`: relative deformation gradient, used in the updated
  Lagrangian formulation.
- `sigma`: Cauchy stress tensor returned to the solid model.
- `impK`: implicit stiffness used as the Laplacian diffusivity.

---

## Developer Notes

### Class role

`StVenantKirchhoffElastic` inherits from `mechanicalLaw` and is added to the
`nonLinGeomMechLaw` runtime selection table. It holds `lambda_`, `mu_` and
`K_` as uniform `dimensionedScalar` values, so a single material entry is
homogeneous.

### Construction

The constructor resolves `E` and `nu` from whichever pair the user supplied,
sets `lambda_` from the plane-stress or plane-strain relation, recomputes
`K_ = lambda_ + (2/3)*mu_`, and finally calls `F().storeOldTime()` and
`Ff().storeOldTime()` so that the old-time deformation gradients exist from
the first time step.

### `correct()` walkthrough

Both the volume-field and the surface-field overloads follow the same five
steps:

1. call `updateF(sigma, mu_, K_)`; a `true` return means the solid model
   requested enforced linearity for stability, the base class has already
   written a Hooke's law stress, and the function returns;
2. materialise `FT` into its own field before contracting;
3. form `C = symm(FT & F)` and the Green strain `E = 0.5*(C - I)`;
4. evaluate `S = 2*mu*E + lambda*tr(E)*I`;
5. push forward with `sigma = (1/J)*transform(F, S)`.

```warning
This file carries the canonical `NOTE [IMPORTANT]` comment that other laws
refer to: never write `F.T() & F` as a single field expression. The lazy
`tmp<>` pathway can drop shear-squared contributions when the transpose is not
materialised, which was observed with g++ 11.x. Always build `FT` first.
```

### Implicit stiffness

`impK()` returns `2*mu + lambda`, the standard segregated-solver diffusivity
for an isotropic elastic material. It is a uniform field, so the same value is
used everywhere within the material.

### Extension points

- The constitutive relation is duplicated between the two `correct()`
  overloads; any change must be applied to both.
- `setRestart()` sets `AUTO_WRITE` on `F` and `Ff` so a restart reproduces the
  same deformation state.

---

## Tutorials

Cases that select `StVenantKirchhoffElastic`:

- `solids/hyperelasticity/rigidRotation/rotatingBlock`
- `solids/hyperelasticity/rigidRotation/rotatingCylinder`
- `solids/hyperelasticity/rigidRotation/rotatingSphere`
- `fluidSolidInteraction/cavityFlexibleBottom`

The `fluidSolidInteraction/3dTube` and `fluidSolidInteraction-preCICE/3dTube`
cases ship the law as a commented alternative to `linearElastic`.
