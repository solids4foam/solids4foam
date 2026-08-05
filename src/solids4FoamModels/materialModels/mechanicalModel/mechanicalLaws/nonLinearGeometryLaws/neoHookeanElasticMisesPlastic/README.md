---
sort: 25
---

# neoHookeanElasticMisesPlastic

This law combines compressible neo-Hookean elasticity with finite-strain
Mises plasticity and isotropic hardening. The runtime type is:

```text
neoHookeanElasticMisesPlastic
```

---

## User Guide

### What it computes

The law updates the isochoric elastic left Cauchy-Green tensor incrementally.
For a relative deformation gradient `relF` and Jacobian `J`, the trial state is

```text
relJ         = J/J_old
relFbar      = relJ^(-1/3)*relF
bEbar_trial  = relFbar bEbar_old relFbar^T
s_trial      = mu*dev(bEbar_trial)
Ibar         = tr(bEbar_trial)/3
muBar        = mu*Ibar
f_trial      = |s_trial| - sqrt(2/3)*J*sigmaY
```

An elastic step has `DLambda = 0`. For a plastic step, the return direction is
`N = s_trial/|s_trial|`. A nonlinear hardening table is integrated by solving

```text
0 = |s_trial| - 2*muBar*DLambda
    - sqrt(2/3)*J*sigmaY(epsilonPEq_old + sqrt(2/3)*DLambda)

DEpsilonPEq = sqrt(2/3)*DLambda
DEpsilonP   = Ibar*DLambda*N
s           = s_trial - 2*mu*DEpsilonP
```

The tabulated `sigmaY` values are treated as Cauchy yield stress versus true
plastic strain. The return calculation multiplies the tabulated stress by `J`
to obtain Kirchhoff yield stress. One table point gives perfect plasticity,
two points give linear hardening, and more than two points activate the Newton
iteration for nonlinear hardening.

The final Cauchy stress is

```text
sigmaHyd = 0.5*K*(J^2 - 1)
sigma    = (sigmaHyd*I + s)/J
```

When `updateBEbarConsistent` is enabled, the spherical part of `bEbar` is
obtained from a cubic equation so that `det(bEbar) = 1`. Otherwise, the trial
value `Ibar` is retained.

The law implements `correct(volSymmTensorField&)` and
`correct(surfaceSymmTensorField&)`. It does not override the point-field or
quadrature-point forms. The point-field form therefore aborts with
`notImplemented`; on non-foam-extend forks the `CompactListList` form also
aborts with `notImplemented`, while that interface is not compiled for
foam-extend.

### Model options

Elastic constants are given as either `E` and `nu`, or `mu` and `K`. The
stress-plastic-strain table is stored in a separate OpenFOAM-format file.

| Entry | Required | Description |
| --- | --- | --- |
| `E` | with `nu` | Young's modulus, `[1 -1 -2 0 0 0 0]` |
| `nu` | with `E` | Poisson's ratio, dimensionless |
| `mu` | with `K` | Shear modulus, `[1 -1 -2 0 0 0 0]` |
| `K` | with `mu` | Bulk modulus, `[1 -1 -2 0 0 0 0]` |
| `"fileName\|file"` | yes | Portable key for the hardening-table file |
| `outOfBounds` | yes | Table bounds handling, for example `clamp` |
| `rho` | yes | Density, `[1 -3 0 0 0 0 0]` |
| `updateBEbarConsistent` | no | Enforce unit determinant; default `yes` |
| `solvePressureEqn` | no | Solve for hydrostatic stress; default `no` |
| `pressureSmoothingScaleFactor` | no | Pressure scale; default `100` |
| `regionName` | conditional | Base solid region if it cannot be detected |

The regular-expression key `"fileName|file"` works with forks that request
either `fileName` or `file`. `outOfBounds` is required by foam-extend; the
OpenFOAM implementation otherwise defaults to `warn`.

The stress table contains pairs of equivalent true plastic strain and Cauchy
yield stress in pressure units. Its strain coordinate is dimensionless. With
two points, the code calculates the linear plastic modulus from their slope.

The `planeStress` entry is read from the enclosing `mechanicalProperties`
dictionary. It changes the conversion from `E` and `nu` to `K`. The law does
not validate the supplied Poisson's ratio. The `maxDeltaErr` entry is read from
`controlDict`, defaults to `0.01`, and controls the time-step recommendation
returned by `newDeltaT()`.

### Recommended dictionary setup

For a steel billet with a separate true-plastic-strain/Cauchy-stress table:

```text
planeStress     no;

mechanical
(
    steel
    {
        type            neoHookeanElasticMisesPlastic;
        rho             rho [1 -3 0 0 0 0 0] 7833;
        E               E [1 -1 -2 0 0 0 0] 200e9;
        nu              nu [0 0 0 0 0 0 0] 0.3;
        "fileName|file" "$FOAM_CASE/constant/plasticStrainVsYieldStress";
        outOfBounds     clamp;

        // Optional
        // updateBEbarConsistent yes;
        // solvePressureEqn no;
        // pressureSmoothingScaleFactor 100;
    }
);
```

The referenced table can, for example, contain two points for linear
hardening:

```text
(
    (0 700e6)
    (10 3700e6)
)
```

### Field glossary

- `F_`, `Ff_`: total deformation gradients at cells and faces, maintained by
  the base class and read on restart.
- `lawJ`, `lawJf`: cell and face Jacobians, created lazily and stored with an
  old-time value.
- `sigmaY`, `sigmaYf`: current Cauchy yield stress; restartable and written.
- `DSigmaY`, `DSigmaYf`: current increments of Cauchy yield stress.
- `epsilonP`, `epsilonPf`: accumulated plastic strain tensors.
- `DEpsilonP`, `DEpsilonPf`: current plastic strain increments, relaxed during
  the stress correction.
- `epsilonPEq`, `epsilonPEqf`: accumulated equivalent plastic strain.
- `DEpsilonPEq`, `DEpsilonPEqf`: current equivalent plastic strain increments.
- `DLambda`, `DLambdaf`: plastic multiplier increments.
- `bEbarTrial`, `bEbarTrialf`: trial isochoric elastic left Cauchy-Green
  tensors.
- `bEbar`, `bEbarf`: corrected isochoric elastic left Cauchy-Green tensors.
- `plasticN`, `plasticNf`: return directions; the cell field stores old-time
  values for the time-integration error estimate.
- `activeYield`: cell yielding flag, set to one where `DEpsilonPEq` is
  positive and zero elsewhere at the end of a time step.
- `sigmaHyd`: cell hydrostatic stress maintained by the base class.

---

## Developer Notes

### Class role

`neoHookeanElasticMisesPlastic` derives directly from `mechanicalLaw` and is
registered in its `nonLinGeomMechLaw` selection table. It stores uniform shear
and bulk moduli, the hardening interpolation table, cell and face plastic
history, trial and corrected elastic strain tensors, return directions, and
the plastic integration controls.

The `OPENFOAM_NOT_EXTEND` preprocessor branches select field-access APIs and
surface-field boundary correction behavior for OpenFOAM versus foam-extend.
They do not remove either implemented stress-correction interface. The law is
listed in both `Make/files.openfoam` and `Make/files.foamextend`.

### Construction

The base constructor first reads `solvePressureEqn` and
`pressureSmoothingScaleFactor`, then selects `regionName` or detects `solid` or
`region0`. The derived member initializers then load the hardening table,
construct the cell and face history fields, read `updateBEbarConsistent`, and
read `maxDeltaErr` from `controlDict`. The body stores the old return direction,
forces creation of the cell and face deformation gradients, and reads either
`E` with `nu` or `mu` with `K`.

For `E` and `nu`, the constructor calculates

```text
mu = E/(2*(1 + nu))
K  = nu*E/((1 + nu)*(1 - 2*nu)) + 2*mu/3   // plane strain or 3-D
K  = nu*E/((1 + nu)*(1 - nu)) + 2*mu/3     // planeStress yes
```

It then classifies the table as perfect, linear, or nonlinear plasticity and
calculates `Hp` for a two-point table. Missing both complete elastic-constant
pairs is a fatal error. Failure to find a base solid region is also fatal.
Construction reports the selected plasticity form and consistent-`bEbar`
setting but performs no Poisson-ratio validation.

At run time, `updateF()` gives fatal errors for a non-incremental updated
Lagrangian model or an unknown nonlinear-geometry type. It warns when a solid
model temporarily enforces linear elasticity. `newDeltaT()` warns when its
plastic-strain integration error exceeds 50 times `maxDeltaErr`.

### Key methods

- `impK()`: returns `scaleFactor*(4*mu/3) + K`, where
  `scaleFactor = 1 - 2*muBar*DLambda/|s_trial|`. It is the diffusivity of the
  solid model's implicit Laplacian term and affects the outer-iteration
  convergence rate rather than the converged answer.
- `correct(volSymmTensorField&)`: updates `F`, `J`, the trial elastic state,
  the plastic return, `bEbar`, hydrostatic stress, and cell Cauchy stress. It
  can use the base-class pressure equation for the hydrostatic part.
- `correct(surfaceSymmTensorField&)`: performs the corresponding face update.
  It evaluates the hydrostatic stress directly because the pressure equation
  is not implemented for surface fields.
- `curYieldStress()` and `yieldFunction()`: interpolate the Cauchy hardening
  curve, convert it to Kirchhoff stress, and form the return residual.
- `newtonLoop()`: solves for `DLambda` by a first-order finite-difference
  Newton iteration when the hardening curve has more than two points.
- `Ibar()`: solves a cubic equation for the spherical part of `bEbar` when
  consistent unit-determinant updating is enabled.
- `residual()`: measures the relative change in the current cell or face
  plastic strain increment between outer iterations.
- `updateTotalFields()`: accumulates yield stress and plastic strain and
  updates `activeYield` at the end of the time step.
- `newDeltaT()`: estimates plastic integration error from the change in return
  direction and proposes a scaled time step.
- `materialTangent()`: is not overridden; the base implementation aborts with
  `notImplemented`, so Newton methods requiring this tangent cannot use the
  law.

### Extension points

A related finite-strain plastic law can copy this class and replace the trial
stress, yield function, and return update while retaining the cell and face
history pattern. New history variables must be accumulated in
`updateTotalFields()`, included in `residual()` where needed, and given restart
write settings in `setRestart()`. A law intended for point-based or
quadrature-point solid models must also implement those `correct` overloads;
one intended for Newton solvers must implement `materialTangent()`.

The new class must be registered in the `nonLinGeomMechLaw` table and added to
the appropriate OpenFOAM and foam-extend build lists.

The source is at
[neoHookeanElasticMisesPlastic.C](https://github.com/solids4foam/solids4foam/blob/master/src/solids4FoamModels/materialModels/mechanicalModel/mechanicalLaws/nonLinearGeometryLaws/neoHookeanElasticMisesPlastic/neoHookeanElasticMisesPlastic.C).

---

## Tutorials

- `solids/elastoplasticity/pipeCrush`
- `solids/elastoplasticity/cooksMembrane`
- `solids/elastoplasticity/cylinderCrush`
- `solids/elastoplasticity/curvedBeams`
- `solids/elastoplasticity/neckingBar`
- `solids/elastoplasticity/impactBar`
- `solids/elastoplasticity/cylinderExpansion`
- `solids/elastoplasticity/upsetBillet`
