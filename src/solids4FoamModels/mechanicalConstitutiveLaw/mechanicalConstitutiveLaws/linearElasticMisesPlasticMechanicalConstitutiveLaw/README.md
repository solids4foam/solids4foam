---
sort: 5
---

# linearElasticMisesPlastic

Small-strain elasto-plasticity with Hookean elasticity and von Mises (J2)
plasticity with isotropic hardening. The runtime type is:

```text
linearElasticMisesPlastic
```

The return mapping is the radial return of Box 3.2 in Simo and Hughes (1998).
The dictionary layout common to all laws is described in the
[material models page](../../../materialModels/README.md).

---

## User Guide

### What it computes

At each integration point, with `epsilon = symm(grad(D))` from the solid
model and the plastic history from the previous time step (`_0`):

```text
sTrial = 2*mu*(dev(epsilon) - dev(epsilonP_0))
qTrial = sqrt(3/2)*mag(sTrial)
fTrial = qTrial - sigmaY(epsilonPEq_0)
```

If `fTrial <= SMALL` the step is elastic and the history is carried over.
Otherwise the increment of equivalent plastic strain `DLambda` is found from
the consistency condition `qTrial - 3*mu*DLambda = sigmaY(epsilonPEq_0 +
DLambda)` and the state is updated:

```text
n          = sTrial/mag(sTrial)
epsilonP   = epsilonP_0 + sqrt(3/2)*DLambda*n
epsilonPEq = epsilonPEq_0 + DLambda
s          = sTrial - 2*mu*sqrt(3/2)*DLambda*n
sigma      = s + K*tr(epsilon)*I
```

with `mu = E/(2*(1 + nu))` and `K = E/(3*(1 - 2*nu))`. The law is small-strain
only; a nonlinear geometry solid model stops with an error. For finite strains
use `neoHookeanElasticMisesPlastic`.

### Hardening input

The yield stress is a table of yield stress versus equivalent plastic strain,
read by solids4foam's `interpolationTable`. The initial yield stress is the
value at zero plastic strain. The number of rows decides the algorithm:

| Rows | Behaviour | Solution |
| --- | --- | --- |
| 1 | Perfect plasticity | Closed form |
| 2 | Linear hardening | Closed form, `DLambda = fTrial/(3*mu + Hp)` |
| 3 or more | Nonlinear hardening | Local Newton iteration |

`Hp` is the slope between the two rows. The Newton iteration uses a
forward-difference derivative and stops when the correction, relative to the
largest strain magnitude among this law's points, falls below `NewtonLoopTol`;
a warning is printed if `NewtonMaxIter` is reached.

### Model options

| Entry | Required | Default | Description |
| --- | --- | --- | --- |
| `rho` | yes | - | Density, `[1 -3 0 0 0 0 0]` |
| `E` | with `nu` | - | Young's modulus, `[1 -1 -2 0 0 0 0]` |
| `nu` | with `E` | - | Poisson's ratio, dimensionless |
| `mu` | with `K` | - | Shear modulus, `[1 -1 -2 0 0 0 0]` |
| `K` | with `mu` | - | Bulk modulus, `[1 -1 -2 0 0 0 0]` |
| `file` or `fileName` | yes | - | Hardening table file |
| `outOfBounds` | yes | - | `clamp`, `warn`, `error` or `repeat` |
| `readerType` | no | `openFoam` | Table format; `csv` is also accepted |
| `NewtonLoopTol` | no | `1e-8` | Newton relative tolerance |
| `NewtonMaxIter` | no | `200` | Newton iteration limit |
| `NewtonFiniteDiffEps` | no | `0.25e-6` | Newton derivative step |

As for `linearElastic`, the elastic constants are given as the `E`/`nu` pair
or the `mu`/`K` pair, their dimensions are checked, and `nu` must satisfy
`-1 < nu < 0.5`.

The table file holds a list of `(equivalentPlasticStrain yieldStress)` pairs.
Tutorials write the file key as `"file|fileName"` and set `outOfBounds clamp`.

`planeStress yes` is a fatal error. To model a thin plate, solve in 3-D with a
`symmetryPlane` back patch and a traction-free front patch.

### Tangents and volumetric split

- Scalar tangent: `(4/3)*mu + K` at elastic points and
  `theta*(4/3)*mu + K` at yielding points, where
  `theta = 1 - 2*mu*sqrt(3/2)*DLambda/mag(sTrial)`. The deviatoric scalar
  tangent drops the `K`.
- Fourth-order tangent: `fourthOrderFiniteDifference` only. Asking for the
  analytical `fourthOrder` tangent is a fatal error.
- Volumetric split: yes. Plasticity is purely deviatoric, so the volumetric
  response is `K*tr(epsilon)`. The law can be used with `solvePressureEqn` and
  with the mixed displacement-pressure formulations.

### State variables

| Name | Type | Role | Description |
| --- | --- | --- | --- |
| `epsilonP` | symmTensor | persistent | Plastic strain tensor |
| `epsilonPEq` | scalar | persistent | Equivalent plastic strain |
| `sigmaY` | scalar | persistent | Current yield stress |

At each write time, for cell-centred solid models, each persistent variable is
written as a volField named after it (`epsilonP`, `epsilonPEq`, `sigmaY`) for
viewing. The restart data are written alongside as
`<material>_<topology>_<variable>`, for example
`aluminium_cellCentredIntegrationPointTopology_epsilonPEq`, together with a
`<material>_<topology>_integrationPoints` file. To restart a run from a later
time, set `restart yes;` in the solid model's coefficients dictionary from the
start of the run: without the old-time displacement gradient it writes, a
restart of a law with history is refused.

### Diagnostics

At the end of each time step the law prints the maximum increment of
equivalent plastic strain and the number of yielding integration points,
summed over all processors and including boundary points. The report is on
by default and controlled by the `linearElasticMisesPlastic` debug switch.

### Example

From the `perforatedPlate` tutorial:

```text
planeStress     no;

mechanical
(
    aluminium
    {
        type            linearElasticMisesPlastic;
        rho             rho [1 -3 0 0 0 0 0] 3500;
        E               E [1 -1 -2 0 0 0 0] 70e+9;
        nu              nu [0 0 0 0 0 0 0] 0.3;
        "file|fileName" "$FOAM_CASE/constant/plasticStrainVsYieldStress";
        outOfBounds     clamp;
    }
);
```

with `constant/plasticStrainVsYieldStress` holding two rows, so linear
hardening:

```text
(
    ( 0     243e6 )
    ( 1     2414e6 )
)
```

### Migrating from the legacy law

- New optional entries: `NewtonLoopTol`, `NewtonMaxIter` and
  `NewtonFiniteDiffEps`. The legacy Newton loop used fixed values of `1e-8`,
  `100` iterations and a `1e-6` step.
- Removed: `solvePressureEquation` (unused in the legacy law), `tangentEps`,
  and `regionName`. `maxDeltaErr` in `controlDict` and the plastic-strain
  based time-step control are gone. `materialTolerance` is ignored.
- Only `epsilonP`, `epsilonPEq` and `sigmaY` are kept. `DEpsilonP`,
  `DEpsilonPEq`, `activeYield`, the face (`...f`) and point (`p...`) copies,
  and the internal `plasticN`, `DLambda` and `DSigmaY` fields are no longer
  written.
- The end-of-step report no longer prints the maximum `epsilonPEq`, and its
  counts now include boundary integration points.
- A fourth-order tangent is available by finite differences
  (`fourthOrderFiniteDifference`); the legacy law built a per-face numerical
  tangent that needed `tangentEps`.
- The return mapping is equivalent to the legacy one: the legacy plastic
  multiplier is expressed here as the equivalent plastic strain increment.

---

## Tutorials

- [perforatedPlate](../../../../../tutorials/solids/elastoplasticity/perforatedPlate/README.md),
  whose regression tests also cover restart, parallel and reconstructed runs

---

## References

- J.C. Simo and T.J.R. Hughes, _Computational Inelasticity_, Springer, 1998.
- P. Cardiff, Z. Tuković, P. De Jaeger, M. Clancy and A. Ivanković, _A
  Lagrangian cell-centred finite volume method for metal forming simulation_,
  Int. J. Numer. Meth. Eng., 2017,
  [10.1002/nme.5345](https://doi.org/10.1002/nme.5345).
