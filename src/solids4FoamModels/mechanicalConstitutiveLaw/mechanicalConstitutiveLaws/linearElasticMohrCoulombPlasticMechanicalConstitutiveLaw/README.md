---
sort: 6
---

# linearElasticMohrCoulombPlastic

Small-strain elasto-plasticity with incremental isotropic Hookean elasticity
and a non-associated Mohr-Coulomb yield surface, returned in principal-stress
space. The runtime type is:

```text
linearElasticMohrCoulombPlastic
```

It is normally used as the effective-stress law inside `poroMechanicalLaw`.
The dictionary layout common to all laws is described in the
[material models page](../../../materialModels/README.md).

---

## User Guide

### What it computes

The stress is carried as a variation `deltaSigma` from the initial stress
`sigma0`. At each integration point a trial stress is built from the strain
increment over the time step and the previous step's variation:

```text
DEpsilon    = symm(grad(D)) - symm(grad(D)_0)
DSigmaTrial = 2*mu*DEpsilon + lambda*tr(DEpsilon)*I
sigmaTrial  = deltaSigma_0 + DSigmaTrial + sigma0

mu     = E/(2*(1 + nu))
lambda = nu*E/((1 + nu)*(1 - 2*nu))        // plane strain / 3-D
lambda = nu*E/((1 + nu)*(1 - nu))          // planeStress yes
```

The principal stresses of the trial stress are ordered
`sigma1 >= sigma2 >= sigma3` and the yield function is evaluated:

```text
k = (1 + sin(frictionAngle))/(1 - sin(frictionAngle))
m = (1 + sin(dilationAngle))/(1 - sin(dilationAngle))
f = k*sigma1 - sigma3 - 2*cohesion*sqrt(k)
```

If `f > SMALL`, the principal stresses are returned to a plane, an edge or the
apex of the yield surface, using `k` for the yield surface and `m` for the
plastic flow direction, and transformed back to tensor form. Otherwise the
trial stress is kept. Finally `deltaSigma = sigma - sigma0` is stored.

Because the trial stress is built on the stored `deltaSigma` rather than on
the stress already held by the solid model, the law works under
`poroMechanicalLaw`, which subtracts the pore pressure afterwards.

The law is small-strain only; a nonlinear geometry solid model stops with an
error.

### Model options

| Entry | Required | Default | Description |
| --- | --- | --- | --- |
| `rho` | yes | - | Density, `[1 -3 0 0 0 0 0]` |
| `E` | yes | - | Young's modulus, `[1 -1 -2 0 0 0 0]` |
| `nu` | yes | - | Poisson's ratio, `[0 0 0 0 0 0 0]` |
| `frictionAngle` | yes | - | Friction angle in degrees, dimensionless |
| `cohesion` | yes | - | Cohesion, `[1 -1 -2 0 0 0 0]` |
| `dilationAngle` | yes | - | Dilation angle in degrees, dimensionless |

When the law is nested inside `poroMechanicalLaw`, `rho` is taken from the
enclosing law's dictionary and need not be repeated. `planeStress` is the
top-level entry of `mechanicalProperties` and changes `lambda` as shown above.

There is no range check on `nu`, `frictionAngle` or `dilationAngle`. The apex
stress is `2*cohesion*sqrt(k)/(k - 1)`, which is singular for a zero friction
angle, so `frictionAngle` must be non-zero. For a pressure-independent yield
surface use `linearElasticMisesPlastic` instead.

### Initial stress

The initial stress is read from a `volSymmTensorField` named `sigma0` in the
start time directory (or in `0`), if one exists; otherwise it is zero. Unlike
`linearElastic`, there is no `sigma0` dictionary entry.

### Tangents and volumetric split

- Scalar tangent: the elastic value `2*mu + lambda`, also returned when the
  deviatoric scalar tangent is asked for.
- Fourth-order tangent: `fourthOrderFiniteDifference` only. Asking for the
  analytical `fourthOrder` tangent is a fatal error.
- Volumetric split: no. The law cannot be used with `solvePressureEqn` or with
  the mixed displacement-pressure formulations, which refuse it at start-up.

### State variables

| Name | Type | Role | Description |
| --- | --- | --- | --- |
| `deltaSigma` | symmTensor | persistent | Stress variation from `sigma0` |
| `activeYield` | scalar | persistent | `1` where yielding, else `0` |
| `sigma0` | symmTensor | prescribed | Initial stress, read from a field |

At each write time, for cell-centred solid models, each persistent variable is
written as a volField named after it (`deltaSigma`, `activeYield`) for
viewing. The restart data are written alongside as
`<material>_<topology>_<variable>`; when the law is nested inside
`poroMechanicalLaw` the sub-law's name is included, for example
`soil_cellCentredIntegrationPointTopology_effectiveStressMechanicalLaw_deltaSigma`.
To restart a run from a later time, set `restart yes;` in the solid model's
coefficients dictionary from the start of the run.

### Example

From the `stripFooting` tutorial:

```text
planeStress     no;

mechanical
(
    soil
    {
        type            poroMechanicalLaw;
        rho             rho [1 -3 0 0 0 0 0] 7854;
        biotCoeff       biotCoeff [0 0 0 0 0 0 0] 1.0;
        effectiveStressMechanicalLaw
        {
            type            linearElasticMohrCoulombPlastic;
            E               E [1 -1 -2 0 0 0 0] 20e6;
            nu              nu [0 0 0 0 0 0 0] 0.3;
            frictionAngle   frictionAngle [0 0 0 0 0 0 0] 30;
            dilationAngle   dilationAngle [0 0 0 0 0 0 0] 0;
            cohesion        cohesion [1 -1 -2 0 0 0 0] 1e5;
        }
    }
);
```

### Migrating from the legacy law

- `regionName` is no longer read. The legacy law silently ignored
  `solvePressureEqn`; the solid model now refuses it for this law.
- The legacy check that `frictionAngle` has a magnitude of at least `1e-3`
  degrees is not present, and no warning is given for a zero `sigma0`.
- Only `deltaSigma` and `activeYield` are kept. The accumulated plastic strain
  `epsilonP`, `epsilonPEq` and the increments `DEpsilon` and `DEpsilonP` are
  no longer computed or written, and the end-of-step report of plastic strain
  and yielding cells is gone.
- The strain increment is taken from the displacement gradient and its
  old-time value, rather than from stored `epsilon` fields.
- The law now serves every integration-point location, and a fourth-order
  tangent is available by finite differences.

---

## Tutorials

Both select this law inside the `effectiveStressMechanicalLaw` sub-dictionary
of `poroMechanicalLaw`:

- [stripFooting](../../../../../tutorials/solids/poroelasticity/stripFooting/README.md)
- [suctionCaission](../../../../../tutorials/solids/poroelasticity/suctionCaission/README.md)
