---
sort: 28
---

# diffusionHyperElastic

This law applies a compressible neo-Hookean relation with shear and bulk
stiffnesses scaled by a normalized mesh-motion diffusivity. The runtime type
is:

```text
diffusionHyperElastic
```

---

## User Guide

### What it computes

The law obtains a raw face diffusivity from the selected `motionDiffusivity`,
limits it relative to its face average, and normalizes it:

```text
dRaw_f = motionDiffusivity
dAvg   = average(dRaw_f)
d_f    = min(maxFactor*dAvg, max(minFactor*dAvg, dRaw_f))/dAvg
d      = average(d_f)
```

The cell field `d` receives zero-gradient boundary conditions. The diffusivity
fields are calculated during construction and are not refreshed by
`correct()`, `impK()`, `bulkModulus()`, or `shearModulus()`.

For deformation gradient `F`, Jacobian `J`, and the volume-preserving left
Cauchy-Green tensor `bBar`, the cell-centred Cauchy stress is:

```text
J     = det(F)
bBar  = J^(-2/3)*symm(F & F.T())
s     = mu*dev(bBar)
sigma = d*(s + 0.5*K*(J^2 - 1)*I)/J
```

The base class updates `F` for incremental or non-incremental total
Lagrangian models and incremental updated Lagrangian models. A
non-incremental updated Lagrangian model produces a fatal error. If the solid
model's `enforceLinear` switch is active, the base class returns a linear
elastic stress before the expression above is evaluated.

The law implements `correct(volSymmTensorField&)`. Its explicitly declared
`correct(surfaceSymmTensorField&)` aborts with `NotImplemented`. The inherited
point-tensor and `CompactListList` quadrature-point overloads also abort with
`notImplemented` on forks where those interfaces are compiled.

### Model options

| Entry | Required | Description |
| --- | --- | --- |
| `mu` | yes | Shear stiffness, `[1 -1 -2 0 0 0 0]` |
| `K` | yes | Bulk stiffness, `[1 -1 -2 0 0 0 0]` |
| `diffusivity` | yes | Mesh-motion diffusivity specification |
| `writeStiffScaleFactor` | yes | Write `stiffnessScaleFactor` if `yes` |
| `maxFactor` | no | Upper limiter relative to the average; default `100` |
| `minFactor` | no | Lower limiter relative to the average; default `0.1` |
| `rho` | yes | Density, `[1 -3 0 0 0 0 0]` |
| `solvePressureEqn` | no | Base-class switch; default `no` |
| `pressureSmoothingScaleFactor` | no | Base-class factor; default `100` |
| `regionName` | no | Override the base mesh region name |

`diffusivity` is passed directly to `motionDiffusivity::New`, so its syntax
must match a mesh-motion diffusivity available in the selected OpenFOAM fork.
Only the ordering of the limiter factors is checked: `maxFactor` less than
`minFactor` is fatal.

The base class reads `solvePressureEqn` and
`pressureSmoothingScaleFactor`, but this law calculates the volumetric stress
directly and never calls `updateSigmaHyd()`. Those entries therefore do not
alter its constitutive response. If `regionName` is omitted, the base class
selects `solid` and then `region0`, in that order.

### Recommended dictionary setup

The following values describe a silicone-like mesh pseudo-material with an
inverse-distance diffusivity based on `movingWall`:

```text
mechanical
(
    siliconeMesh
    {
        type            diffusionHyperElastic;
        rho             rho [1 -3 0 0 0 0 0] 1100;
        mu              mu [1 -1 -2 0 0 0 0] 0.6e6;
        K               K [1 -1 -2 0 0 0 0] 60e6;
        diffusivity     quadratic inverseDistance (movingWall);
        writeStiffScaleFactor yes;

        // Optional
        // maxFactor    100;
        // minFactor    0.1;
        // regionName   region0;
        // solvePressureEqn no;
        // pressureSmoothingScaleFactor 100;
    }
);
```

### Field glossary

- `distf`: normalized and limited face diffusivity, `d_f`.
- `stiffnessScaleFactor`: cell-centred `d`, obtained by averaging `distf`;
  written only when `writeStiffScaleFactor yes` is selected.
- `F_`: total deformation gradient. Its old-time value is stored during
  construction, and the field is written for incremental restarts.
- `relF_`: relative deformation gradient maintained by the base class.
- `mu_<name>`, `K_<name>`: cell fields populated from `mu` and `K` when
  `updateF()` is first called.
- `bulkModulus`: temporary cell field returned as `K*d`.
- `shearModulus`: temporary cell field returned as `mu*d`.
- `k`: temporary implicit-stiffness field returned by `impK()`.

---

## Developer Notes

### Class role

`diffusionHyperElastic` derives directly from `mechanicalLaw`. It stores
uniform `dimensionedScalar` values `mu_` and `K_`, scalar limiter values
`maxFactor_` and `minFactor_`, an `autoPtr<motionDiffusivity>`, and the face
and cell scale fields `df_` and `d_`.

The complete implementation is guarded by `#ifdef OPENFOAM_NOT_EXTEND`. Its
source appears in `Make/files.openfoam` but not in `Make/files.foamextend`, so
the law is built for the OpenFOAM.com and OpenFOAM.org path, but not for
foam-extend.

### Construction

The base constructor first reads or inserts `solvePressureEqn` and
`pressureSmoothingScaleFactor`, then resolves `regionName`. The derived
constructor proceeds in this order:

1. It reads `mu` and `K` as dimensioned scalars.
2. It reads `maxFactor` and `minFactor`, using `100` and `0.1` by default.
3. It constructs the selected motion diffusivity from `diffusivity`.
4. It initializes `distf` from that model and initializes
   `stiffnessScaleFactor` by face-to-cell averaging.
5. It creates `F_` and stores its old-time value.
6. It aborts if `maxFactor` is less than `minFactor`.
7. It corrects, limits, and normalizes the diffusivity fields.
8. It reads `writeStiffScaleFactor` and enables automatic writing of the cell
   scale field when requested.

Missing required entries and an unknown motion-diffusivity type produce the
corresponding OpenFOAM dictionary or selection errors. The base class emits a
fatal error if it cannot find `solid` or `region0` and no `regionName` was
given. During stress correction, it emits a fatal error for non-incremental
updated Lagrangian kinematics. It warns once per time step when
`enforceLinear` is active.

### Key methods

- `updateD()`: calls the motion diffusivity's `correct()`, limits and
  normalizes its face field, averages it to cells, and corrects the cell-field
  boundary conditions. It is called by the constructor.
- `impK()`: returns `d*((4/3)*mu + K)`. This is the diffusivity used by the
  solid model's implicit Laplacian term and affects the outer-iteration
  convergence rate rather than the converged answer.
- `correct(volSymmTensorField&)`: asks the base class to update `F`, then
  applies the diffusivity-scaled neo-Hookean relation unless linearity was
  enforced.
- `correct(surfaceSymmTensorField&)`: aborts with `NotImplemented`.
- `bulkModulus()` and `shearModulus()`: return `K*d` and `mu*d` cell fields.
- `materialTangent()`: is not overridden; the base implementation aborts with
  `notImplemented`.

### Extension points

A related diffusivity-scaled finite-strain law can be made by copying this
class, changing the cell-centred constitutive expression, and retaining
`updateD()` when the same mesh-motion diffusivity and limiter are appropriate.
Implement the face, point, and quadrature-point `correct` overloads, and
`materialTangent()`, if the new law must support face-based, vertex-centred,
higher-order, block-coupled, or PETSc SNES solid models. Add the source to each
fork-specific build list that supports its dependencies.

The source is at
[diffusionHyperElastic.C](https://github.com/solids4foam/solids4foam/blob/master/src/solids4FoamModels/materialModels/mechanicalModel/mechanicalLaws/nonLinearGeometryLaws/diffusionHyperElastic/diffusionHyperElastic.C).

---

## Tutorials

No tutorial currently uses `diffusionHyperElastic`.
