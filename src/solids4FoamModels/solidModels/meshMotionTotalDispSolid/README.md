---
sort: 16
---

# meshMotionTotalDispSolid

This page documents the solid model for moving a mesh owned by another model,
e.g. the fluid mesh with the `solids4foamSolidModel` mesh motion solver. The
runtime type is:

```text
meshMotionTotalDisplacement
```

The model is available on OpenFOAM.com only, as is the mesh motion solver.

---

## User Guide

### What this model solves

The model solves the quasi-static momentum equation, `div(sigma) = 0`, for
the total displacement `D` from the reference points, but discretised on the
current mesh, `x = X + D`, as `displacementLaplacian` does. So no copy of the
reference mesh is needed, and, because the unknown is the total displacement
rather than an increment, the mesh for a given boundary motion does not depend
on how the boundary got there.

- With `finiteStrain no` (the default), the small-strain law, e.g.
  `diffusionElastic`, is evaluated from the displacement gradient on the
  current mesh.
- With `finiteStrain yes`, the kinematics are those of the current
  configuration, `Finv = I - gradD^T`, `F = inv(Finv)` and `J = det(F)`, and
  the finite-strain law, e.g. `diffusionHyperElastic`, gives the Cauchy
  stress, which acts on the current face areas. Once the mesh is at
  `x = X + D`, this is the spatial form of the total Lagrangian equations.

The mesh is not moved by the model: its owner moves it to the reference
points plus `pointD`. The model has no inertia and no gravity, and the
implicit stiffness is updated before each solution, as the stiffness of the
law can change with the mesh. Only history-independent (elastic) laws are
supported. Traction conditions act on the current face areas;
`useUndeformedArea` is not supported.

The model stops if the displacement inverts the mesh, i.e. if the
determinant of `Finv` is not positive.

The corrector loop stops when the relative residual falls below `rTol`
(default `1e-6`), the residual below `aTol` (default `1e-50`) or the
relative correction below `sTol` (default `1e-10`), and fails after
`nCorrectors` (default `10000`) iterations. The older `solutionTolerance` and
`alternativeTolerance` keys are not read.

### Example

From the `beamInCrossFlow` tutorial (`./Allrun diffusionHyperElastic`), in
`constant/fluid/solidProperties`:

```text
solidModel     meshMotionTotalDisplacement;

meshMotionTotalDisplacementCoeffs
{
    // Finite strain kinematics, for the finite strain law
    finiteStrain    yes;

    nCorrectors     1000;
    rTol            1e-04;
    infoFrequency   100;
}
```

---

## Tutorials

- `fluidSolidInteraction/beamInCrossFlow` (`./Allrun diffusionElastic` and
  `./Allrun diffusionHyperElastic`)
