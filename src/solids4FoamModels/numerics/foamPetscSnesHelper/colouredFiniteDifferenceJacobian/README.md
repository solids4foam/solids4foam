# Coloured Finite-Difference Jacobian

`colouredFiniteDifferenceJacobian` computes the Jacobian of a
`foamPetscSnesHelper` residual by coloured forward differences (PETSc
`MatColoring` and `MatFDColoring`).

## Use in the PETSc SNES Solid Models

The cell-centred solid models with a PETSc SNES path
(`linearGeometryTotalDisplacement`, `nonLinearGeometryTotalLagrangianTotalDisplacement`
and `nonLinearGeometryUpdatedLagrangian`) select it in their coefficients
dictionary:

```text
solutionAlgorithm PETScSNES;
jacobian          colouredFiniteDifference; // default: compact
//jacobianCheck          yes;  // check every build, not only the first
//jacobianCheckTolerance 1e-4;
```

By default (`jacobian compact`) the matrix given to SNES is the compact
approximate Jacobian of the model, a vector Laplacian with the same coefficient
for each displacement component. It omits the coupling between the components,
so it represents rotations, and hence bending, poorly: for slender solids the
Krylov iterations per Newton iteration grow with the slenderness.
`colouredFiniteDifference` replaces it by the Jacobian of the residual:

- The pattern is the square of the compact matrix, i.e. the cells within two
  layers, in full 2x2 or 3x3 component blocks.
- PETSc colours the pattern (distance two) and evaluates the residual once per
  colour (about 40-50 colours in 2-D, about 120 in 3-D).
- The Jacobian is rebuilt at every Jacobian evaluation, i.e. at every Newton
  iteration by default. With `snes_lag_jacobian "-2"` and
  `snes_lag_jacobian_persists true` it is built once per time step, because
  the options are applied again before each solve.
- With `snes_mf_operator` it is the preconditioning matrix; without it, it is
  also the Newton Jacobian.
- After the first build the matrix is compared with a forward difference of
  the residual along a random direction. A relative difference above
  `jacobianCheckTolerance` stops the run.

PETSc options of the colouring and of the differencing take the prefix
`solid_fd_`, e.g. `solid_fd_mat_fd_coloring_err`.

## Limitations

- Cell-centred models only. The vertex-centred models already assemble their
  Jacobian exactly.
- The residual must depend only on cells within two layers. Residuals that
  reach further fail the directional check and stop the run: the high-order
  residual, cells where `leastSquaresS4f` falls back to the point-cell
  stencil, and `tractionPatchSurfaceGradient`. High-order support is
  tracked in issue #789.
- The residual must be a function of the unknowns only, as for
  `snes_mf_operator`.
- In 3-D a hexahedral cell has 25 cells in its two-layer stencil, i.e. 225
  matrix entries per row instead of 63, which increases memory and the cost
  of the preconditioner.
- A build costs one residual evaluation per colour. It pays off when the
  compact Jacobian needs many Krylov iterations (slender, bending solids).
  Otherwise the build time can exceed what the iterations save:

  | case | Krylov/Newton, compact to FD | run time |
  | --- | --- | --- |
  | `cantilever2d` (320 x 16) | 79 to 1 (Newton 3 to 1) | 0.42 to 0.20 s |
  | `cantilever2d`, 4 ranks | 205 to 19 (Newton 3 to 2) | 0.35 to 0.10 s |
  | `cantileverVibration` (3-D, 130 steps) | 27.5 to 1.8 | 21.5 to 24.3 s |
  | `rotatingBlock`, updated Lagrangian | 8.8 to 2.7 | 0.29 to 0.34 s |
