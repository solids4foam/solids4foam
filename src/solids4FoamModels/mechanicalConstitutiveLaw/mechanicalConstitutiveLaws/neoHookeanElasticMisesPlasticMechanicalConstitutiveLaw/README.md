---
sort: 25
---

# neoHookeanElasticMisesPlastic

Finite-strain elasto-plasticity with compressible neo-Hookean elasticity and
von Mises (J2) plasticity with isotropic hardening. The runtime type is:

```text
neoHookeanElasticMisesPlastic
```

The formulation is the multiplicative one of Simo and Hughes (1998), with the
isochoric elastic left Cauchy-Green tensor `bEbar` carried as history. The
dictionary layout common to all laws is described in the
[material models page](../../../materialModels/README.md).

---

## User Guide

### What it computes

At each integration point, from `F`, `J = det(F)` and their values at the
start of the time step (`_0`), with `bEbar_0` from the previous step:

```text
relJ       = J/J_0
relFbar    = relJ^(-1/3)*(F & inv(F_0))
bEbarTrial = symm(relFbar & bEbar_0 & relFbar.T())
sTrial     = mu*dev(bEbarTrial)
Ibar       = tr(bEbarTrial)/3
muBar      = Ibar*mu
fTrial     = mag(sTrial) - sqrt(2/3)*J*sigmaY(epsilonPEq_0)
```

The step is elastic if `fTrial < SMALL`. Otherwise the plastic multiplier
`DLambda` is found from

```text
0 = mag(sTrial) - 2*muBar*DLambda
  - sqrt(2/3)*J*sigmaY(epsilonPEq_0 + sqrt(2/3)*DLambda)
```

and, with `N = sTrial/mag(sTrial)`,

```text
epsilonPEq = epsilonPEq_0 + sqrt(2/3)*DLambda
s          = sTrial - 2*mu*Ibar*DLambda*N
bEbar      = s/mu + Ibar*I
sigmaHyd   = 0.5*K*(J^2 - 1)
sigma      = (sigmaHyd*I + s)/J
```

The hardening curve is Cauchy yield stress against equivalent true plastic
strain; it is multiplied by `J` to give the Kirchhoff yield stress the return
is written in. With `updateBEbarConsistent yes` the spherical part of `bEbar`
is found from a cubic so that `det(bEbar) = 1` (Rubin and Attia 1996, in the
form of Hollenstein et al. 2013), rather than taken as `Ibar`.

`mu = E/(2*(1 + nu))` and `K = E/(3*(1 - 2*nu))`. The law is finite-strain
only; a linear geometry solid model stops with an error. For small strains use
`linearElasticMisesPlastic`.

### Hardening input

The yield stress is read with OpenFOAM's `interpolationTable`, constructed
from this law's dictionary. The number of rows decides the algorithm:

| Rows | Behaviour | Solution |
| --- | --- | --- |
| 1 | Perfect plasticity | Closed form |
| 2 | Linear hardening | Closed form |
| 3 or more | Nonlinear hardening | Local Newton iteration |

For two rows, `Hp` is the slope between them and
`DLambda = fTrial/(2*muBar*(1 + Hp/(3*muBar)))`. The Newton iteration uses a
forward-difference derivative and stops when the correction, relative to the
largest `mag(bEbarTrial)` over this law's points on all processors, falls
below `NewtonLoopTol`; a warning is printed if `NewtonMaxIter` is reached.

### Model options

| Entry | Required | Default | Description |
| --- | --- | --- | --- |
| `rho` | yes | - | Density, `[1 -3 0 0 0 0 0]` |
| `E` | with `nu` | - | Young's modulus, `[1 -1 -2 0 0 0 0]` |
| `nu` | with `E` | - | Poisson's ratio, dimensionless |
| `mu` | with `K` | - | Shear modulus, `[1 -1 -2 0 0 0 0]` |
| `K` | with `mu` | - | Bulk modulus, `[1 -1 -2 0 0 0 0]` |
| `file` or `fileName` | yes | - | Hardening table file |
| `outOfBounds` | no | `warn` | `clamp`, `warn`, `error` or `repeat` |
| `readerType` | no | `openFoam` | Table format; `csv` is also accepted |
| `updateBEbarConsistent` | no | `yes` | Enforce `det(bEbar) = 1` |
| `NewtonLoopTol` | no | `1e-8` | Newton relative tolerance |
| `NewtonMaxIter` | no | `200` | Newton iteration limit |
| `NewtonFiniteDiffEps` | no | `0.25e-6` | Newton derivative step |
| `solvePressureEqn` | no | `no` | Hydrostatic smoothing, see below |
| `pressureSmoothingScaleFactor` | no | `100` | Scale of that smoothing |

If both `E` and `nu` are present they are used; otherwise `mu` and `K` are
required, and are converted to `E` and `nu`. Their dimensions are checked,
and `nu` must satisfy `-1 < nu < 0.5`. The law has no incompressible limit.

The table file holds a list of `(equivalentPlasticStrain yieldStress)` pairs.
Tutorials write the file key as `"file|fileName"` and set `outOfBounds clamp`;
foam-extend requires `outOfBounds`.

`planeStress yes` is a fatal error. To model a thin plate, solve in 3-D with a
`symmetryPlane` back patch and a traction-free front patch.

### Tangents and volumetric split

- Scalar tangent: `theta*(4/3)*mu + K`, where
  `theta = 1 - 2*muBar*DLambda/mag(sTrial)`, so `theta = 1` at elastic
  points. The deviatoric scalar tangent drops the `K`.
- Fourth-order tangent: `fourthOrderFiniteDifference` only. Asking for the
  analytical `fourthOrder` tangent is a fatal error.
- Volumetric split: yes. J2 plasticity is isochoric, so the return acts on `s`
  alone: the isochoric stress is `s/J` and the volumetric response is
  `sigmaHyd/J = dU/dJ` for `U(J) = 0.25*K*(J^2 - 1 - 2*ln(J))`. The law can be
  used with `solvePressureEqn` and with the mixed displacement-pressure
  formulations, as for
  [`neoHookeanElastic`](../neoHookeanElasticMechanicalConstitutiveLaw/README.md).
- The isochoric stress is not invariant under a superposed dilation once the
  material yields, because the yield surface is scaled by `J`.

### Hydrostatic stress smoothing

`solvePressureEqn yes` asks the solid model to replace the explicit
`sigmaHyd` with a smoothed field of that name, which is written. It is done
by the solid model, not the law, and only with
`nonLinearGeometryUpdatedLagrangian`,
`nonLinearGeometryTotalLagrangianTotalDisplacement` or
`linearGeometryTotalDisplacement`, the `implicitSegregated` algorithm, a
single material, and without `solvePressure`. Any other combination is
refused at start-up.

### State variables

| Name | Type | Role | Description |
| --- | --- | --- | --- |
| `bEbar` | symmTensor | persistent | Isochoric elastic `b`, initially `I` |
| `epsilonPEq` | scalar | persistent | Equivalent plastic strain |
| `sigmaY` | scalar | persistent | Current Cauchy yield stress |

`sigmaY` starts at the table value at zero plastic strain. At each write time,
for cell-centred solid models, each variable is written as a volField named
after it (`bEbar`, `epsilonPEq`, `sigmaY`) for viewing. These fields carry no
dimensions and are never read. The restart data are written alongside as
`<material>_<topology>_<variable>`, for example
`steel_cellCentredIntegrationPointTopology_epsilonPEq`, together with a
`<material>_<topology>_integrationPoints` file. To restart a run from a later
time, set `restart yes;` in the solid model's coefficients dictionary from the
start of the run: without the old-time displacement gradient it writes, a
restart of a law with history is refused. A restart whose state files are
missing is also refused.

### Diagnostics

At the end of each time step the law reports the maximum increment of
equivalent plastic strain and the number of yielding integration points,
summed over all processors and including boundary points. The report is off
by default; enable it with the `neoHookeanElasticMisesPlastic` debug switch.

### Example

From the `neckingBar` tutorial:

```text
planeStress     no;

mechanical
(
    steel
    {
        type            neoHookeanElasticMisesPlastic;
        rho             rho [ 1 -3 0 0 0 0 0 ] 7833;
        E               E [ 1 -1 -2 0 0 0 0 ] 200e9;
        nu              nu [ 0 0 0 0 0 0 0 ] 0.3;
        "file|fileName" "$FOAM_CASE/constant/plasticStrainVsYieldStress";
        outOfBounds     clamp;
    }
);
```

with `constant/plasticStrainVsYieldStress` holding eight rows, so nonlinear
hardening:

```text
(
    (0.000    0.451e9)
    (0.006    0.476e9)
    (0.019    0.525e9)
    (0.038    0.583e9)
    (0.066    0.642e9)
    (0.147    0.710e9)
    (0.500    0.777e9)
    (1.000    0.831e9)
)
```

### Differences from the legacy law

- New optional entries: `NewtonLoopTol`, `NewtonMaxIter` and
  `NewtonFiniteDiffEps`, with the defaults the legacy law had fixed.
- Removed: `regionName`, and `maxDeltaErr` in `controlDict` with the
  plastic-strain based time-step control.
- `planeStress yes` is now a fatal error; the legacy law used it to compute
  `K`. `nu` is now validated.
- Only `bEbar`, `epsilonPEq` and `sigmaY` are kept. `epsilonP`, `DEpsilonP`,
  `DEpsilonPEq`, `DSigmaY`, `DLambda`, `bEbarTrial`, `plasticN`,
  `activeYield`, the Jacobian fields and the face (`...f`) copies are no
  longer written.
- The plastic strain increment is not under-relaxed between outer
  iterations, and the law no longer supplies a plastic-strain residual.
- The hydrostatic smoothing is done by the solid model rather than the law,
  and is refused where it is not supported rather than ignored.
- A fourth-order tangent is available by finite differences; the legacy law
  had none.
- The end-of-step report now counts boundary integration points.

---

## Tutorials

- [neckingBar](../../../../../tutorials/solids/elastoplasticity/neckingBar/README.md),
  whose regression tests also cover restart
- [cylinderExpansion](../../../../../tutorials/solids/elastoplasticity/cylinderExpansion/README.md),
  with `solvePressureEqn`
- [pipeCrush](../../../../../tutorials/solids/elastoplasticity/pipeCrush/README.md)
- [cooksMembrane](../../../../../tutorials/solids/elastoplasticity/cooksMembrane/README.md)
- [cylinderCrush](../../../../../tutorials/solids/elastoplasticity/cylinderCrush/README.md)
- [curvedBeams](../../../../../tutorials/solids/elastoplasticity/curvedBeams/README.md)
- [impactBar](../../../../../tutorials/solids/elastoplasticity/impactBar/README.md)
- [upsetBillet](../../../../../tutorials/solids/elastoplasticity/upsetBillet/README.md)

---

## References

- J.C. Simo and T.J.R. Hughes, _Computational Inelasticity_, Springer, 1998,
  [10.1007/b98904](https://doi.org/10.1007/b98904).
- M.B. Rubin and A. Attia, Int. J. Numer. Meth. Eng., 39, 309-320, 1996.
- M. Hollenstein, M. Jabareen and M.B. Rubin, Comput. Mech., 52, 649-667,
  2013.
