---
sort: 11
---

# diffusionElastic

This law applies isotropic small-strain elasticity with shear and bulk
stiffnesses scaled by a normalized mesh-motion diffusivity. The runtime type
is:

```text
diffusionElastic
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
fields are calculated during construction and are not refreshed by either
`correct` overload.

For non-incremental solid models, the strain is calculated from the total
displacement gradient. For incremental solid models, the strain increment is
added to the old-time strain:

```text
epsilon  = symm(grad(D))
epsilon  = epsilon.oldTime() + symm(grad(DD))    // incremental

epsilon_f = symm(grad(D)_f)
epsilon_f = epsilon_f.oldTime() + symm(grad(DD)_f) // incremental
```

The cell-centred and face-centred stresses are then:

```text
sigma   = d*(2*mu*dev(epsilon) + kappa*tr(epsilon)*I)
sigma_f = d_f*(2*mu*dev(epsilon_f) + kappa*tr(epsilon_f)*I)
```

The law implements `correct(volSymmTensorField&)` and
`correct(surfaceSymmTensorField&)`. The point-tensor overload is inherited
from `mechanicalLaw` and aborts with `notImplemented`. The
`CompactListList` quadrature-point overload is also inherited and aborts with
`notImplemented` on forks where that interface is compiled.

### Model options

| Entry | Required | Description |
| --- | --- | --- |
| `mu` | yes | Shear stiffness, `[1 -1 -2 0 0 0 0]` |
| `kappa` | yes | Bulk stiffness, `[1 -1 -2 0 0 0 0]` |
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

The following values use aluminium-like stiffness and density for a mesh
pseudo-material with an inverse-distance diffusivity based on `movingWall`:

```text
mechanical
(
    aluminiumMesh
    {
        type            diffusionElastic;
        rho             rho [1 -3 0 0 0 0 0] 2700;
        mu              mu [1 -1 -2 0 0 0 0] 26e9;
        kappa           kappa [1 -1 -2 0 0 0 0] 69e9;
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
- `epsilon`, `epsilonf`: small-strain tensors at cell centres and faces. Their
  old-time values are stored during construction for incremental solid models.
- `bulkModulus`: temporary cell field returned as `kappa*d`.
- `shearModulus`: temporary cell field returned as `mu*d`.
- `k`: temporary implicit-stiffness field returned by `impK()`.

---

## Developer Notes

### Class role

`diffusionElastic` derives directly from `mechanicalLaw`. It stores uniform
`dimensionedScalar` values `mu_` and `kappa_`, scalar limiter values
`maxFactor_` and `minFactor_`, an `autoPtr<motionDiffusivity>`, and the face
and cell scale fields `df_` and `d_`.

The class contains no fork-specific preprocessor guards. Its source appears in
`Make/files.openfoam` but not in `Make/files.foamextend`, so it is built for
the OpenFOAM.com and OpenFOAM.org path, but not for foam-extend.

### Construction

The base constructor first reads or inserts `solvePressureEqn` and
`pressureSmoothingScaleFactor`, then resolves `regionName`. The derived
constructor proceeds in this order:

1. It reads `mu` and `kappa` as dimensioned scalars.
2. It reads `maxFactor` and `minFactor`, using `100` and `0.1` by default.
3. It constructs the selected motion diffusivity from `diffusivity`.
4. It initializes `distf` from that model and initializes
   `stiffnessScaleFactor` by face-to-cell averaging.
5. It stores old-time values for `epsilon` and `epsilonf`.
6. It aborts if `maxFactor` is less than `minFactor`.
7. It corrects, limits, and normalizes the diffusivity fields.
8. It reads `writeStiffScaleFactor` and enables automatic writing of the cell
   scale field when requested.

Missing required entries and an unknown motion-diffusivity type produce the
corresponding OpenFOAM dictionary or selection errors. The base class emits a
fatal error if it cannot find `solid` or `region0` and no `regionName` was
given. The derived constructor emits no warnings.

### Key methods

- `updateD()`: calls the motion diffusivity's `correct()`, limits and
  normalizes its face field, averages it to cells, and corrects the cell-field
  boundary conditions. It is called by the constructor.
- `impK()`: returns `d*((4/3)*mu + kappa)`. This is the diffusivity used by the
  solid model's implicit Laplacian term and affects the outer-iteration
  convergence rate rather than the converged answer.
- `correct(volSymmTensorField&)`: updates the cell strain and applies the
  diffusivity-scaled constitutive relation.
- `correct(surfaceSymmTensorField&)`: updates the face strain and applies the
  corresponding face-scaled relation.
- `bulkModulus()` and `shearModulus()`: return `kappa*d` and `mu*d` cell
  fields.
- `materialTangent()`: is not overridden; the base implementation aborts with
  `notImplemented`.

### Extension points

A related spatially scaled small-strain law can be made by copying this class,
changing the two constitutive expressions, and retaining `updateD()` when the
same mesh-motion diffusivity and limiter are appropriate. Implement the point
and quadrature-point `correct` overloads, and `materialTangent()`, if the new
law must support vertex-centred, higher-order, block-coupled, or PETSc SNES
solid models. Add the source to each fork-specific build list that supports
its dependencies.

The source is at
[diffusionElastic.C](https://github.com/solids4foam/solids4foam/blob/master/src/solids4FoamModels/materialModels/mechanicalModel/mechanicalLaws/linearGeometryLaws/diffusionElastic/diffusionElastic.C).

---

## Tutorials

No tutorial currently uses `diffusionElastic`.
