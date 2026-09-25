---
sort: 15
---

# StVenantKirchhoffElastic

St. Venant-Kirchhoff hyperelasticity: the second Piola-Kirchhoff stress is a
linear function of the Green-Lagrange strain. The runtime type is:

```text
StVenantKirchhoffElastic
```

It is the simplest finite-strain extension of Hooke's law. It captures large
rotations exactly but is only suited to small or moderate strains. The
dictionary layout common to all laws is described in the
[material models page](../../../materialModels/README.md).

---

## User Guide

### What it computes

At each integration point, from `F` and `J = det(F)`:

```text
E     = 0.5*(symm(F.T() & F) - I)
S     = 2*mu*E + lambda*tr(E)*I
sigma = symm(F & S & F.T())/J
```

which is the Cauchy stress of `Psi = mu*(E && E) + 0.5*lambda*tr(E)^2`.

The law is finite-strain only; a linear geometry solid model stops with an
error, and `linearElastic` should be used there. It holds no history.

```warning
`Psi` is quadratic in the Green strain, so it does not penalise `J`
approaching zero, and the law has no check on `J`. Do not use it for large
volumetric compression; prefer `neoHookeanElastic`.
```

### Model options

| Entry | Required | Default | Description |
| --- | --- | --- | --- |
| `rho` | yes | - | Density, `[1 -3 0 0 0 0 0]` |
| `E` | with `nu` | - | Young's modulus, `[1 -1 -2 0 0 0 0]` |
| `nu` | with `E` | - | Poisson's ratio, dimensionless |
| `mu` | with `K` | - | Shear modulus, `[1 -1 -2 0 0 0 0]` |
| `K` | with `mu` | - | Bulk modulus, `[1 -1 -2 0 0 0 0]` |

Exactly one of the pairs `E`/`nu` and `mu`/`K` must be given; both, or
neither, is a fatal error. Their dimensions are checked, and `nu` must satisfy
`-1 < nu < 0.5`. A `mu`/`K` pair is converted to `E` and `nu` with the 3-D
relations, and then

```text
mu     = E/(2*(1 + nu))
lambda = E*nu/((1 + nu)*(1 - 2*nu))    planeStress no
lambda = E*nu/((1 + nu)*(1 - nu))      planeStress yes
K      = lambda + (2/3)*mu
```

so under `planeStress yes` the bulk modulus used is the plane-stress one, not
a `K` given in the dictionary.

### Tangents and volumetric split

- Scalar tangent: `2*mu + lambda`; the deviatoric scalar tangent is
  `(4/3)*mu`.
- Fourth-order tangent: `fourthOrderFiniteDifference` only. Asking for the
  analytical `fourthOrder` tangent is a fatal error.
- Volumetric split: no. The energy is written on the full Green strain, so
  the law cannot separate its isochoric stress from its volumetric response.
  A mixed displacement-pressure formulation refuses it, naming the material,
  and `solvePressureEqn` is refused at start-up.

### State variables

None. The law writes no restart files and no state fields.

### Example

From the `beamInCrossFlow` fluid-solid interaction tutorial
(`constant/solid/mechanicalProperties`):

```text
planeStress no;

mechanical
(
    rubber
    {
        type        StVenantKirchhoffElastic;
        rho         rho [1 -3 0 0 0 0 0] 1000;
        E           E [1 -1 -2 0 0 0 0] 1e4;
        nu          nu [0 0 0 0 0 0 0] 0.4;
    }
);
```

### Differences from the legacy law

- `solvePressureEqn` is now a fatal error. The legacy law accepted it and
  `pressureSmoothingScaleFactor`, but never used them.
- `regionName` is no longer read.
- `nu` is now validated.
- The legacy fallback to a Hookean stress when the solid model enforced
  linearity is gone.
- A fourth-order tangent is available by finite differences.

---

## Tutorials

- [rigidRotation](../../../../../tutorials/solids/hyperelasticity/rigidRotation/README.md):
  `rotatingBlock`, `rotatingCylinder` and `rotatingSphere`
- [beamInCrossFlow](../../../../../tutorials/fluidSolidInteraction/beamInCrossFlow/README.md)
- [cavityFlexibleBottom](../../../../../tutorials/fluidSolidInteraction/cavityFlexibleBottom/README.md)

The [3dTube](../../../../../tutorials/fluidSolidInteraction/3dTube/README.md)
tutorial keeps the law as a commented alternative to `linearElastic`.
