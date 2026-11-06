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

### Solution algorithms

`solutionAlgorithm` selects how the equations are solved:

- `implicitSegregated` (the default): a segregated corrector loop. It stops
  when the relative residual falls below `rTol` (default `1e-6`), the
  residual below `aTol` (default `1e-50`) or the relative correction below
  `sTol` (default `1e-10`), and fails after `nCorrectors` (default `10000`)
  iterations. The older `solutionTolerance` and `alternativeTolerance` keys
  are not read.
- `PETScSNES`: the PETSc nonlinear solver (SNES). The residual is that of the
  segregated loop, the divergence of the stress plus the stabilisation on
  the current mesh, with the same kinematics and law, so the two algorithms
  converge to the same displacement. The approximate Jacobian, used as the
  preconditioner of a matrix-free Newton-Krylov method (`snes_mf_operator`),
  is the implicit stiffness Laplacian of the segregated loop, assembled on the
  current mesh. Each solution starts from the last displacement with the new
  boundary values. The PETSc options are read from the `D` entry of the
  `solvers` in `system/<region>/fvSolution` of the moving mesh, or else from
  the `optionsFile` (default `petscOptions`); `stopOnPetscError` (default
  `yes`) stops the run if the SNES solve fails. The `rTol`, `nCorrectors`
  and the other corrector loop entries are not used, and the
  `cellDisplacements` and high-order options are not supported.

The quasi-monolithic coupling (`NewtonQuasiMonolithic`) needs `PETScSNES`, and
stops with a message if the segregated algorithm is selected. The partitioned
couplings accept either algorithm. In both cases, the stiffness scale factor
of the `solids4foamSolidModel` mesh motion solver is recalculated on the
current mesh before each solution (`updateStiffnessScaleFactor`, default
`yes` for this model).

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

For the quasi-monolithic coupling, select PETSc SNES:

```text
solidModel     meshMotionTotalDisplacement;

meshMotionTotalDisplacementCoeffs
{
    finiteStrain        yes;
    solutionAlgorithm   PETScSNES;
}
```

with the PETSc options in `system/fluid/fvSolution`, e.g. as in the
`TurekHron2006` tutorial (`system/fluid/fvSolution.monolithic`):

```text
solvers
{
    D // fluid mesh motion
    {
        solver    petsc;

        options
        {
            snes_type newtonls;
            snes_linesearch_type basic;
            snes_rtol "1e-3";
            snes_stol "1e-3";
            snes_max_it "100";
            snes_monitor;
            snes_converged_reason;

            snes_mf;
            snes_mf_operator;

            ksp_type lgmres;
            ksp_gmres_restart "100";
            ksp_converged_reason;

            pc_type bjacobi;
            sub_pc_type lu;

            snes_lag_preconditioner_persists true;
            snes_lag_preconditioner "-2";
        }
    }
}
```

---

## Tutorials

- `fluidSolidInteraction/beamInCrossFlow` (`./Allrun diffusionElastic` and
  `./Allrun diffusionHyperElastic`)
