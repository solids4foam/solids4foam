---
sort: 11
---

# diffusionElastic

Small-strain linear elasticity whose stiffness is scaled by a mesh motion
diffusivity, for moving a fluid mesh with a solid model. The runtime type is:

```text
diffusionElastic
```

The law is used by the `solids4foamSolidModel` fluid mesh motion solver, which
solves a solid model on a static copy of the fluid mesh; see the
[fluid-solid interface page](../../../fluidSolidInterfaces/README.md). The
dictionary layout common to all laws is described in the
[material models page](../../../materialModels/README.md), and the framework
itself in the [framework README](../../README.md).

---

## User Guide

### What it computes

The stress is that of [`linearElastic`](../linearElasticMechanicalConstitutiveLaw/README.md)
scaled by the stiffness scale factor `d`:

```text
sigma = d*(2*mu*epsilon + lambda*tr(epsilon)*I)
```

`d` is read from the `stiffnessScaleFactor` field, at the law's integration
points (face values are interpolated from the cells). The
`solids4foamSolidModel` mesh motion solver calculates it once, on the
reference mesh, from the `diffusivity` in `constant/fluid/dynamicMeshDict`, as
the face diffusivity limited to between `minFactor` (default `0.1`) and
`maxFactor` (default `100`) times its average, divided by its average, and
averaged to the cells. With `quadratic inverseDistance (interface)`, the mesh
is stiffer, and deforms less, near the interface. Set
`writeStiffScaleFactor yes` in the `dynamicMeshDict` to write the field. The
law can also be used with a `stiffnessScaleFactor` field given in the time
directory.

### Model options

| Entry | Required | Description |
| --- | --- | --- |
| `rho` | yes | Density |
| `E`, `nu` or `mu`, `K` | yes | Elastic constants, as for `linearElastic` |

`kappa` is accepted in place of `K`. An initial stress, `sigma0`, is not accepted.

### Tangents and volumetric split

The scalar tangent and the analytical fourth-order tangent are those of
`linearElastic` scaled by `d`; the finite difference fourth-order tangent is
that of the scaled law. The law does not separate its volumetric response, so
it cannot be used with a mixed displacement-pressure solid model.

### Example

From the `beamInCrossFlow` tutorial (`./Allrun diffusionElastic`), in
`constant/fluidMeshMotion/mechanicalProperties`:

```text
planeStress     no;

mechanical
(
    meshMotion
    {
        type        diffusionElastic;
        rho         rho [1 -3 0 0 0 0 0] 1;
        mu          mu [1 -1 -2 0 0 0 0] 1;
        K           K [1 -1 -2 0 0 0 0] 3;
    }
);
```

### Migrating from the legacy law

- `diffusivity`, `maxFactor`, `minFactor` and `writeStiffScaleFactor` are now
  entries of the `solids4foamSolidModel` mesh motion solver, in
  `constant/fluid/dynamicMeshDict`; a `diffusivity` entry in the law is
  ignored with a warning.
- The face stress uses the face-interpolated cell `d` rather than the limited
  face diffusivity.

---

## Tutorials

- `fluidSolidInteraction/beamInCrossFlow` (`./Allrun diffusionElastic`)
