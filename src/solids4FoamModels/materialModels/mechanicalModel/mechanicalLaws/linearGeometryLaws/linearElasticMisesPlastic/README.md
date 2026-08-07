---
sort: 5
---

# linearElasticMisesPlastic

This page documents the small-strain elasto-plastic mechanical law with
Hookean elasticity and von Mises (J2) plasticity. The runtime type is:

```text
linearElasticMisesPlastic
```

The return-mapping algorithm is the radial return of Box 3.2 in Simo and
Hughes, _Computational Inelasticity_ (1998). The implementation is described
in P. Cardiff, Z. Tuković, P. De Jaeger, M. Clancy and A. Ivanković (2017), _A
Lagrangian cell-centred finite volume method for metal forming simulation_,
[10.1002/nme.5345](https://doi.org/10.1002/nme.5345).

---

## User Guide

### What it computes

For each cell and each boundary face the law performs an elastic trial
followed by a radial return:

```text
sTrial   = 2*mu*(dev(epsilon) - dev(epsilonP.oldTime()))
fTrial   = mag(sTrial) - sqrt(2/3)*sigmaY.oldTime()
```

If `fTrial` is below `SMALL` the step is elastic. Otherwise the return
direction is `plasticN = sTrial/mag(sTrial)` and a plastic multiplier
increment `DLambda` is computed. The deviatoric and hydrostatic parts are then
recombined:

```text
DEpsilonP = DLambda*plasticN
s         = sTrial - 2*mu*DEpsilonP
sigmaHyd  = K*tr(epsilon)
sigma     = sigmaHyd*I + s
```

`sigmaHyd` is routed through `updateSigmaHyd()`, so it comes from the smoothed
pressure equation when `solvePressureEqn` is enabled.

### Hardening input

The post-yield behaviour is supplied as a table of equivalent plastic strain
versus yield stress, read by an `interpolationTable`. The number of rows in
that table decides the algorithm:

| Rows | Behaviour | Solution method |
| --- | --- | --- |
| 1 | Perfect plasticity | Closed form |
| 2 | Linear hardening | Closed form, modulus `Hp` |
| 3 or more | Nonlinear hardening | Newton iteration per cell |

With one or two rows, `DLambda = fTrial/(2*mu)`, divided by `1 + Hp/(3*mu)`
when the linear plastic modulus `Hp` is non-zero. With three or more rows a
local Newton loop solves the yield function to a relative tolerance of `1e-8`
in at most 100 iterations, using a first-order finite difference derivative
with a step of `1e-6`; a warning is printed if it fails to converge.

### Model options

| Entry | Required | Description |
| --- | --- | --- |
| `E` | with `nu` | Young's modulus, `[1 -1 -2 0 0 0 0]` |
| `nu` | with `E` | Poisson's ratio, dimensionless |
| `mu` | with `K` | Shear modulus, `[1 -1 -2 0 0 0 0]` |
| `K` | with `mu` | Bulk modulus, `[1 -1 -2 0 0 0 0]` |
| `rho` | yes | Density, `[1 -3 0 0 0 0 0]` |
| `file` or `fileName` | yes | Hardening table file |
| `outOfBounds` | yes | `interpolationTable` bounds handling |

As for `linearElastic`, the elastic constants are given as either the `E`/`nu`
pair or the `mu`/`K` pair, and a missing pair is a fatal error.

The hardening table entries are those of the standard OpenFOAM
`interpolationTable` constructed from a dictionary. Tutorials write the file
key as `"file|fileName"` so that the same dictionary works across forks, and
set `outOfBounds clamp`. The referenced file contains a list of
`(plasticStrain yieldStress)` pairs.

Two further entries are read outside the material sub-dictionary or in a way
worth flagging:

| Entry | Where | Description |
| --- | --- | --- |
| `maxDeltaErr` | `controlDict` | Target plastic-strain integration error |
| `tangentEps` | material dict | Perturbation for the numerical tangent |

`maxDeltaErr` defaults to `0.01` and is only consulted when the solid model
requests adaptive time-stepping through `newDeltaT()`. `tangentEps` is read
with `lookup()` — that is, it is required, with no default — but only inside
`materialTangentField()`, so it is needed only by the block-coupled and PETSc
SNES solid models.

```warning
The constructor also reads a `solvePressureEquation` switch (default `no`),
but that member is never used anywhere else in the class. The switch that
actually takes effect is the base-class `solvePressureEqn`. Do not rely on
`solvePressureEquation`.
```

`planeStress yes` is a fatal error. The header suggests the workaround of
solving in 3-D with a `symmetryPlane` back patch and a traction-free front
patch.

### Recommended dictionary setup

```text
planeStress     no;

mechanical
(
    aluminium
    {
        type            linearElasticMisesPlastic;
        rho             rho [1 -3 0 0 0 0 0] 2700;
        E               E [1 -1 -2 0 0 0 0] 70e9;
        nu              nu [0 0 0 0 0 0 0] 0.3;
        "file|fileName" "$FOAM_CASE/constant/plasticStrainVsYieldStress";
        outOfBounds     clamp;
        solvePressureEqn no;
    }
);
```

with `constant/plasticStrainVsYieldStress` containing, for example:

```text
(
    (0    250e6)
    (0.1  300e6)
    (0.5  400e6)
)
```

### Field glossary

Written to the time directory (`AUTO_WRITE`):

- `sigmaY`, `sigmaYf`, `pSigmaY`: current yield stress at cells, faces and
  points. All three are `READ_IF_PRESENT`, so a restart picks up the hardened
  state.
- `epsilonP`, `epsilonPf`, `pEpsilonP`: accumulated plastic strain tensor.
- `epsilonPEq`, `epsilonPEqf`, `pEpsilonPEq`: accumulated equivalent plastic
  strain.
- `DEpsilonP`, `DEpsilonPf`: plastic strain increment for the current step.
- `DEpsilonPEq`, `DEpsilonPEqf`: equivalent plastic strain increment.
- `activeYield`: `1` in cells that yielded in the current step, `0` elsewhere.

Internal only: `plasticN`, `plasticNf`, `pPlasticN` (return direction),
`DLambda`, `DLambdaf`, `pDLambda` (plastic multiplier increment) and
`DSigmaY`, `DSigmaYf`, `pDSigmaY` (yield-stress increment).

`updateTotalFields()` prints the maximum `DEpsilonPEq`, the maximum
`epsilonPEq` and the number of actively yielding cells at the end of each time
step; those three lines are the quickest health check on a plastic run.

---

## Developer Notes

### Class role

`linearElasticMisesPlastic` derives from `mechanicalLaw` and maintains a
complete triple of state fields — cell, face and point — for every plastic
quantity. This is what allows it to serve the cell-centred, face-based (`uns`)
and vertex-centred solid models from one class. It is compiled for all three
forks, appearing in both `Make/files.openfoam` and `Make/files.foamextend`.

Note that the class is registered with a debug level of `1`, unlike the other
laws in this subsection which use `0`.

### Construction

The initialiser list builds the `interpolationTable` first, because the
initial yield stress `stressPlasticStrainSeries_(0.0)` is needed as the
uniform value of `sigmaY`, `sigmaYf` and `pSigmaY`. The body then forces
old-time storage for every state field, reads the elastic pair, and classifies
the hardening as perfect, linear or nonlinear, computing `Hp_` in the linear
case.

### Key methods

- `updatePlasticity()`: the per-point kernel. It is called once per cell and
  once per boundary face from `correct()`, and decides between the elastic
  branch and the return mapping.
- `newtonLoop()`: local Newton solve for `DLambda` in the nonlinear hardening
  case, calling `yieldFunction()` twice per iteration.
- `correct(volSymmTensorField&)`: the full cell-centred update, ending with
  `updateSigmaHyd()`. The surface overload delegates to a shared
  `epsilon`-taking implementation.
- `residual()`: relative change in `DEpsilonP` (or `DEpsilonPf`, if a face
  gradient field exists) between outer iterations. This feeds the solid
  model's `materialTolerance` check.
- `newDeltaT()`: estimates the plastic-strain time integration error from the
  rotation of the return direction over the step, following Lee and Bathe
  (1994), and scales `deltaT` so that the error matches `maxDeltaErr`. It warns
  if the error exceeds 50 times the target.
- `impK()`, `impKdiagTensor()`: scalar and diagonal-tensor implicit
  stiffnesses.
- `materialTangentField()`: builds a per-face 6x6 tangent by numerical
  perturbation of `grad(D)f`, one symmetric-tensor component at a time. The
  branch is preceded by a commented-out `numericalTangent` lookup; there is
  currently no analytic alternative.

### Extension points

- Replacing the `interpolationTable` with a run-time selectable hardening law
  would remove the row-count-based dispatch in the constructor.
- The Newton loop constants (`LoopTol_`, `MaxNewtonIter_`, `finiteDiff_`) are
  compile-time members rather than dictionary entries; exposing them would
  help for materials with steep hardening curves.
- Kinematic hardening would require a back-stress field alongside `sigmaY` and
  a modified `sTrial`.

The source is at
[linearElasticMisesPlastic.C](https://github.com/solids4foam/solids4foam/blob/master/src/solids4FoamModels/materialModels/mechanicalModel/mechanicalLaws/linearGeometryLaws/linearElasticMisesPlastic/linearElasticMisesPlastic.C).

---

## Tutorials

Cases that select `linearElasticMisesPlastic`:

- `solids/elastoplasticity/perforatedPlate`
