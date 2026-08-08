---
sort: 6
---

# linearElasticMohrCoulombPlastic

This law combines incremental isotropic Hookean elasticity with a
principal-stress Mohr-Coulomb return calculation. The runtime type is:

```text
linearElasticMohrCoulombPlastic
```

---

## User Guide

### What it computes

The cell-centred and face-centred `correct` functions form a trial stress from
the strain increment and the stress variation stored at the old time:

```text
DEpsilon   = epsilon - epsilon.oldTime()
DSigmaTrial = 2*mu*DEpsilon + lambda*tr(DEpsilon)*I
sigmaTrial  = deltaSigma.oldTime() + DSigmaTrial + sigma0

mu     = E/(2*(1 + nu))
lambda = nu*E/((1 + nu)*(1 - 2*nu))        // plane strain / 3-D
lambda = nu*E/((1 + nu)*(1 - nu))          // planeStress yes
```

The code obtains the principal stresses, orders them as `sigma1 >= sigma2 >=
sigma3`, and evaluates

```text
k = (1 + sin(frictionAngle))/(1 - sin(frictionAngle))
m = (1 + sin(dilationAngle))/(1 - sin(dilationAngle))
f = k*sigma1 - sigma3 - 2*cohesion*sqrt(k)
```

The angles are read in degrees and converted to radians in the sine
expressions. If `f > SMALL`, `calculateStress()` applies a non-associated
return in principal-stress space. It uses `k` for the yield surface and `m`
for the plastic-potential direction, and selects a plane, an edge, or the apex
according to its projection tests. The returned principal stresses are then
transformed back to tensor form. Otherwise, the elastic trial stress is kept.

The law implements `correct(volSymmTensorField&)` and
`correct(surfaceSymmTensorField&)`. The inherited point-centred overload aborts
with `notImplemented`. On OpenFOAM.com and OpenFOAM.org, the inherited
`CompactListList` quadrature-point overload also aborts with `notImplemented`;
that overload is not compiled for foam-extend.

### Model options

| Entry | Required | Description |
| --- | --- | --- |
| `E` | yes | Young's modulus, `[1 -1 -2 0 0 0 0]` |
| `nu` | yes | Poisson's ratio, `[0 0 0 0 0 0 0]` |
| `frictionAngle` | yes | Friction angle in degrees, dimensionless |
| `cohesion` | yes | Cohesion, `[1 -1 -2 0 0 0 0]` |
| `dilationAngle` | yes | Dilation angle in degrees, dimensionless |
| `rho` | yes | Density, `[1 -3 0 0 0 0 0]` |
| `solvePressureEqn` | no | Base switch; default `no` |
| `pressureSmoothingScaleFactor` | no | Base scalar; default `100` |
| `regionName` | no | Base mesh region name override |

`planeStress` is read from the enclosing `mechanicalProperties` dictionary,
not from the law dictionary. It changes the value of `lambda` as shown above.
The constructor does not check the ranges of `nu`, `frictionAngle`, or
`dilationAngle`.

The base class reads `solvePressureEqn` and
`pressureSmoothingScaleFactor`, but this law never calls
`updateSigmaHyd()`. Consequently, those entries do not alter its stress
calculation. `regionName` is used when the base class finds the solid model,
and `rho` is read when the solid model requests the density.

Initial stress is supplied through the `sigma0` field. It is not a dictionary
entry read by this law.

### Recommended dictionary setup

```text
planeStress     no;

mechanical
(
    soil
    {
        type            linearElasticMohrCoulombPlastic;
        rho             rho [1 -3 0 0 0 0 0] 2000;
        E               E [1 -1 -2 0 0 0 0] 20e6;
        nu              nu [0 0 0 0 0 0 0] 0.3;
        frictionAngle   frictionAngle [0 0 0 0 0 0 0] 30;
        cohesion        cohesion [1 -1 -2 0 0 0 0] 1e5;
        dilationAngle   dilationAngle [0 0 0 0 0 0 0] 0;

        // Optional for a non-standard solid region
        // regionName   solid;
    }
);
```

### Field glossary

- `epsilon_`, `epsilonf_`: cell- and face-centred small-strain fields; both
  read if present and retain old-time values.
- `deltaSigma_`, `deltaSigmaf`: stress variation from `sigma0` at cells and
  faces. They are neither read nor written.
- `DEpsilon`, `DEpsilonf`: cell and face strain increments. They are read if
  present and not written.
- `DEpsilonP`, `DEpsilonPf`: cell and face plastic-strain increments.
  `DEpsilonP` is read if present and written; `DEpsilonPf` is neither read nor
  written.
- `epsilonP`: accumulated cell-centred plastic strain, written at output
  times.
- `epsilonPEq`: accumulated equivalent plastic strain, written at output
  times.
- `activeYield`: cell-centred yielding flag, read if present and written; `1`
  denotes active yielding and `0` denotes an elastic state.
- `sigma0`, `sigma0f`: initial stress at cells and its face interpolation,
  maintained by the base class.

---

## Developer Notes

### Class role

`linearElasticMohrCoulombPlastic` derives directly from `mechanicalLaw` and is
registered in the `linGeomMechLaw` runtime-selection table. It stores the
elastic constants `E_`, `nu_`, `lambda_`, `mu_`, and `K_`; the three input
plasticity properties; the derived return-mapping vectors and tensors; and the
stress, strain-increment, plastic-strain, and yielding fields listed above.

The law is listed in both `Make/files.openfoam` and `Make/files.foamextend`.
The implementation uses `OPENFOAM_NOT_EXTEND` guards for mathematical
constants and field access, and an `OPENFOAM_COM` guard for face addressing.
The inherited quadrature-point interface is guarded by `#ifndef FOAMEXTEND`.

### Construction

The base constructor first reads `solvePressureEqn` and
`pressureSmoothingScaleFactor`, then selects `regionName` explicitly or finds
`solid` or `region0`. Failure to find a solid region causes a fatal error.

The derived constructor then reads `E`, `nu`, `frictionAngle`, `cohesion`, and
`dilationAngle`. It derives `lambda`, `mu`, `K`, `k`, `m`, the elastic
principal-stress matrix and its inverse, the plane and edge directions, and
the apex stress. It constructs all history and output fields, stores old-time
values for `epsilon` and `epsilonf`, and creates `sigma0`.

Missing required entries or entries with incompatible dimensions fail during
dictionary lookup. A zero `sigma0` field produces a warning because the return
calculation is stress-state dependent. During stress correction, the custom
eigensolver can warn about a zero root, complex eigenvalues, or eigenvectors
that it cannot determine; the latter two conditions replace the affected
values with zero.

### Key methods

- `impK()`: returns the uniform field `2*mu + lambda`. This is the diffusivity
  of the solid model's implicit Laplacian term and affects the outer-iteration
  convergence rate rather than the converged answer.
- `correct(volSymmTensorField&)`: updates cell strain, forms the incremental
  trial stress, returns every cell and boundary value to the yield surface
  when required, and stores `deltaSigma_` for the residual calculation.
- `correct(surfaceSymmTensorField&)`: performs the same operation on faces and
  transfers each internal-face yielding result to its owner and neighbour
  cells.
- `calculateStress()`: perturbs positive tensor components slightly, computes
  principal values and directions, evaluates the yield function, and applies
  the plane, edge, or apex return.
- `calculateEigens()`: computes eigenvalues analytically and reconstructs the
  associated eigenvectors from sub-determinants.
- `residual()`: returns the global maximum relative change in cell or face
  `deltaSigma`, choosing the face form when a face displacement-gradient field
  exists.
- `updateTotalFields()`: reconstructs the plastic-strain increment in yielding
  cells, accumulates `epsilonP` and `epsilonPEq`, and reports their maxima and
  the number of yielding cells.
- `materialTangent()`: is not overridden. Calling the inherited implementation
  aborts with `notImplemented`.

### Extension points

A related small-strain pressure-dependent plastic law can copy this class and
replace the yield test and the plane, edge, and apex return data in
`calculateStress()`. A new implementation must keep the cell and face paths
consistent, update the plastic history in `updateTotalFields()`, and add a
material tangent if it is to be used by Newton or block-coupled paths that
request one.

The source is at
[linearElasticMohrCoulombPlastic.C](https://github.com/solids4foam/solids4foam/blob/master/src/solids4FoamModels/materialModels/mechanicalModel/mechanicalLaws/linearGeometryLaws/linearElasticMohrCoulombPlastic/linearElasticMohrCoulombPlastic.C).

---

## Tutorials

- `solids/poroelasticity/stripFooting`
- `solids/poroelasticity/suctionCaission`

Both cases select this law inside the `effectiveStressMechanicalLaw`
sub-dictionary of `poroMechanicalLaw`.
