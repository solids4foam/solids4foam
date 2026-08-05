---
sort: 17
---

# neoHookeanElastic

This page documents the compressible neo-Hookean hyperelastic law. The runtime
type is:

```text
neoHookeanElastic
```

The law uses the volumetric-isochoric split of Simo and Hughes (1998), Eqn
9.2.6. It is the default choice for large-strain rubber-like materials in
solids4foam, and it is also the elastic backbone of several other laws in the
toolbox.

Reference: J. C. Simo and T. J. R. Hughes, _Computational Inelasticity_,
Springer, 1998,
[10.1007/b98904](https://doi.org/10.1007/b98904).

---

## User Guide

### Strain energy function

The stored energy is split into an isochoric and a volumetric part:

```text
Psi = 0.5*mu*(tr(bBar) - 3) + U(J)

bBar = J^(-2/3)*b,   b = F & F.T(),   J = det(F)
```

with the default volumetric part

```text
U(J) = 0.5*K*(0.5*(J^2 - 1) - ln(J))
```

The resulting Cauchy stress, which is what the code evaluates, is

```text
sigma = (1/J)*(0.5*K*(J^2 - 1)*I + mu*dev(bBar))
```

Setting `alternatePressureDefinition` replaces the Kirchhoff hydrostatic term
`0.5*K*(J^2 - 1)` by `K*(J - 1)`, which corresponds to
`U(J) = K*(J - ln(J) - 1)`.

The law is registered in the nonlinear-geometry table, so it works with both
total-Lagrangian and updated-Lagrangian solid models; the base class
`updateF()` selects the right deformation-gradient update. It is not available
to linear-geometry solid models.

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

When `E` and `nu` are given, `mu = E/(2*(1 + nu))` and `K` follows from the
plane-stress or plane-strain relation, as set by the top-level `planeStress`
entry of `constant/mechanicalProperties`.

Optional entries:

| Entry | Default | Description |
| --- | --- | --- |
| `alternatePressureDefinition` | `false` | Use `K*(J - 1)` instead |
| `pressureDisplacement` | `false` | Pressure-displacement coupling |
| `pressureDisplacementCoeff` | `-1` | Negative means auto `3*nu/(1 + nu)` |
| `tangentEps` | none | Perturbation for `materialTangentField()` |
| `solvePressureEqn` | `false` | Smooth `sigmaHyd` with a Laplacian |
| `pressureSmoothingScaleFactor` | `100` | Scale factor for that smoothing |

`tangentEps` is read with `lookup()` inside `materialTangentField()`, so it is
only required by solid models that ask for the material tangent; it is not
read at construction.

### Pressure-displacement mode

Setting `pressureDisplacement true` switches the law into the form expected by
[coupledPressureDisplacementSolid](https://www.solids4foam.com/documentation/solid-models/coupledPressureDisplacementSolid.html).
In that mode:

- the deviatoric stress becomes `mu*(b - I)/J`, i.e. the full left
  Cauchy-Green tensor is used rather than its isochoric part;
- the hydrostatic part is taken from the separately solved `p` field (`pf` on
  faces) and multiplied by `pressureDisplacementCoeff`;
- `impK()` returns `mu` alone rather than `(4/3)*mu + K`;
- if `E` and `nu` are supplied with `nu` at or above `0.5 - SMALL`, `K` is set
  to `GREAT`, which allows a truly incompressible material to be specified.

`pressureDisplacementCoeff` defaults to `3*nu/(1 + nu)`, with `nu` recovered
from `mu` and `K`, or taken as `0.5` when `K` has been set to `GREAT`.

### Recommended dictionary setup

Minimal example for `constant/mechanicalProperties`:

```text
planeStress     no;

mechanical
(
    rubber
    {
        type            neoHookeanElastic;
        rho             rho [1 -3 0 0 0 0 0] 1000;
        E               E [1 -1 -2 0 0 0 0] 2e+06;
        nu              nu [0 0 0 0 0 0 0] 0.45;
    }
);
```

For a nearly incompressible material solved with a pressure equation, add
`solvePressureEqn yes;` and, if needed, `pressureSmoothingScaleFactor`.

### Field glossary

- `F`, `Ff`: total deformation gradient at cell centres and faces.
- `relF`, `relFf`: relative deformation gradient, used in the updated
  Lagrangian formulation.
- `sigmaHyd`: hydrostatic Cauchy stress, `-p`; written when
  `solvePressureEqn` is enabled.
- `sigma`: Cauchy stress tensor returned to the solid model.
- `impK`: implicit stiffness used as the Laplacian diffusivity.

---

## Developer Notes

### Class role

`neoHookeanElastic` inherits from `mechanicalLaw` and is added to the
`nonLinGeomMechLaw` runtime selection table. It stores `mu_` and `K_` as
uniform `dimensionedScalar` values, so the law is homogeneous within a given
material entry; heterogeneity is obtained by using several `mechanical`
entries with the multi-material machinery.

### Construction

The constructor reads the two switches and the coefficient first, then decides
between the `(E, nu)` and `(mu, K)` branches. It finishes by calling
`F().storeOldTime()` and `Ff().storeOldTime()`, which forces the old-time
deformation gradients to exist from the first time step.

### `correct()` walkthrough

`correct(volSymmTensorField&)`:

1. calls `updateF(sigma, mu_, K_)`; if that returns `true` the solid model has
   requested enforced linearity for stability, the base class has already set
   a Hooke's law stress, and the function returns immediately;
2. materialises `FT` explicitly before contracting, to avoid a `tmp<>`
   evaluation problem documented in `StVenantKirchhoffElastic.C`;
3. forms `J`, `b`, and `bBar`, then the deviatoric stress `s`;
4. in pressure-displacement mode looks up `p` and returns
   `-pressureDisplacementCoeff*p*I + s` directly;
5. otherwise calls `updateSigmaHyd()` with either `K*(J - 1)` or
   `0.5*K*(J^2 - 1)`, then assembles `sigma = (1/J)*(sigmaHyd*I + s)`.

The surface overloads route through the private `correctF()` helper, which is
also used by the `correct(sigma, gradD)` overload that builds `F` from a
supplied displacement gradient.

### Implicit stiffness

`impK()` returns `(4/3)*mu + K`, which equals `2*mu + lambda` and is the
standard choice for the segregated Laplacian operator. In
pressure-displacement mode it returns `mu` only, because the volumetric
response is carried by the pressure equation.

`bulkModulus()` and `shearModulus()` return uniform fields built from `K_` and
`mu_`.

### Material tangent

`materialTangentField()` builds a face-based `6x6` tangent by finite
differences: each component of `grad(D)f` is perturbed by `tangentEps` and the
resulting change in `sigmaf` is divided by the perturbation. It looks up
`sigmaf` and `grad(D)f` from the registry, so it only works with solid models
that create those surface fields.

### Extension points

- `correctF()` is the single place where the constitutive relation is
  evaluated for surface fields; a new volumetric-isochoric split can be added
  there and in `correct(volSymmTensorField&)`.
- `setRestart()` sets `AUTO_WRITE` on `F` and `Ff` so a restart reproduces the
  same deformation state.

---

## Tutorials

Cases that select `neoHookeanElastic`:

- `solids/hyperelasticity/blockPunch`
- `solids/hyperelasticity/cantileverBeam`
- `solids/hyperelasticity/cylindricalPressureVessel`, in the three
  `caseOptions/pressureDisplacement` variants (`cylinder`, `cylinderLin`
  and `cylinderUnsteady`)
- `solids/linearElasticity/plateHole`, in both
  `caseOptions/pressureDisplacement` variants (`compressible` and
  `incompressible`)
- `fluidSolidInteraction/beamInCrossFlow`
- `fluidSolidInteraction/fillingElasticContainer`
- `fluidSolidInteraction/flexibleDamBreak`
- `fluidSolidInteraction/HronTurekFsi3`
