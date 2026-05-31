# `newtonIcoFluid`

`newtonIcoFluid` is an incompressible Newtonian laminar fluid model that solves
the coupled velocity-pressure system with PETSc SNES. It is intended for cases
where the pressure-velocity coupling should be handled as one nonlinear system
instead of through the segregated PIMPLE loop used by `pimpleFluid`.

The model is selected in `constant/fluidProperties`:

```c++
fluidModel newtonIcoFluid;

newtonIcoFluidCoeffs
{
    predictor no;

    forceImplicitFlux yes;
    fluidFluxExtrapolationAlgorithm1 yes;

    addDivPhiUDamping no;
    pressureScaleFactor 1.0;

    stabilisation
    {
        momentum
        {
            type        diffStencilLaplacian;
            scaleFactor 0.0;
        }

        pressure
        {
            type        RhieChow;
            scaleFactor 1.0;
        }
    }
}
```

The default pressure stabilisation is `RhieChow` with scale factor `1.0`. The
default momentum stabilisation is `diffStencilLaplacian` with scale factor
`0.0`, so it is inactive unless explicitly enabled. Older dictionaries that put
the stabilisation entries directly under `stabilisation` are interpreted as
pressure-stabilisation settings.

## User Notes

### PETSc Solver Setup

`newtonIcoFluid` solves the combined `Up` field. The unknown block layout is
`(Ux, Uy, p)` for 2-D and `(Ux, Uy, Uz, p)` for 3-D. PETSc options are normally
specified in the `Up` solver entry in `system/fvSolution`:

```c++
solvers
{
    Up
    {
        solver    petsc;

        options
        {
            snes_type newtonls;
            snes_linesearch_type l2;
            snes_rtol "1e-4";
            snes_stol "1e-4";
            snes_max_it "10";
            snes_monitor;
            snes_converged_reason;
            snes_mf_operator;

            ksp_type lgmres;
            ksp_gmres_restart "100";
            ksp_rtol "1e-2";
            ksp_pc_side right;
            ksp_converged_reason;
            ksp_max_it "400";
        }
    }
}
```

The `snes_mf_operator` option asks PETSc to use a matrix-free operator while
still allowing an assembled preconditioning matrix. This is the usual setup for
this model.

### Schur Fieldsplit Preconditioner

The most robust starting point is a PETSc Schur `fieldsplit` preconditioner.
For a 2-D case, the velocity fields are components `0,1` and pressure is
component `2`:

```c++
options
{
    snes_type newtonls;
    snes_linesearch_type l2;
    snes_rtol "1e-4";
    snes_stol "1e-4";
    snes_max_it "10";
    snes_monitor;
    snes_converged_reason;
    snes_mf_operator;

    ksp_type lgmres;
    ksp_gmres_restart "100";
    ksp_rtol "1e-2";
    ksp_pc_side right;
    ksp_converged_reason;
    ksp_max_it "400";

    pc_type fieldsplit;
    pc_fieldsplit_block_size "3";
    pc_fieldsplit_type schur;
    pc_fieldsplit_schur_factorization_type full;
    pc_fieldsplit_schur_precondition a11;

    pc_fieldsplit_0_fields "0,1";
    pc_fieldsplit_1_fields "2";

    fieldsplit_0_ksp_type preonly;
    fieldsplit_0_pc_type asm;
    fieldsplit_0_sub_pc_type ilu;
    fieldsplit_0_sub_pc_factor_levels "0";

    fieldsplit_1_ksp_type preonly;
    fieldsplit_1_pc_type hypre;
    fieldsplit_1_pc_hypre_type boomeramg;
    fieldsplit_1_pc_hypre_boomeramg_max_iter "1";
    fieldsplit_1_pc_hypre_boomeramg_strong_threshold "0.6";
    fieldsplit_1_pc_hypre_boomeramg_grid_sweeps_up "1";
    fieldsplit_1_pc_hypre_boomeramg_grid_sweeps_down "1";
    fieldsplit_1_pc_hypre_boomeramg_agg_nl "1";
    fieldsplit_1_pc_hypre_boomeramg_agg_num_paths "1";
    fieldsplit_1_pc_hypre_boomeramg_max_levels "25";
    fieldsplit_1_pc_hypre_boomeramg_coarsen_type HMIS;
    fieldsplit_1_pc_hypre_boomeramg_interp_type ext+i;
    fieldsplit_1_pc_hypre_boomeramg_P_max "1";
    fieldsplit_1_pc_hypre_boomeramg_truncfactor "0.3";
}
```

For a 3-D case, use block size `4`, velocity fields `"0,1,2"`, and pressure
field `"3"`.

### SIMPLE-Type Physics Preconditioner

The model also provides a `physics` PETSc preconditioner path that applies one
segregated SIMPLE-type correction through `newtonIcoFluid::precondition(...)`.
This is useful for benchmarking and tuning because it is cheap per application
and can reduce Krylov iterations. It may not always reduce wall time in
parallel, because the inner OpenFOAM pressure correction can dominate.

Example PETSc options:

```c++
options
{
    snes_type newtonls;
    snes_linesearch_type basic;
    snes_rtol "1e-4";
    snes_stol "1e-4";
    snes_max_it "10";
    snes_monitor;
    snes_converged_reason;
    snes_mf_operator;

    ksp_type fgmres;
    ksp_gmres_restart "100";
    ksp_rtol "1e-2";
    ksp_pc_side right;
    ksp_converged_reason;
    ksp_max_it "400";

    pc_type physics;
    pc_physics_type model;
}
```

The preconditioner solves temporary `dU` and `dp` fields, so the corresponding
schemes and solver dictionaries must exist:

```c++
laplacianSchemes
{
    laplacian(nuEff,dU) Gauss linear corrected;
    laplacian(rAUf,dp) Gauss linear corrected;
}
```

```c++
solvers
{
    dU
    {
        solver          PBiCG;
        preconditioner  DILU;
        tolerance       0;
        relTol          0;
        minIter         1;
        maxIter         1;
    }

    dp
    {
        solver          GAMG;
        tolerance       0;
        relTol          0;
        minIter         1;
        maxIter         1;
        smoother        GaussSeidel;
        nPreSweeps      0;
        nPostSweeps     1;
        nFinestSweeps   1;
        scaleCorrection true;
        directSolveCoarsest false;
        cacheAgglomeration true;
        nCellsInCoarsestLevel 20;
        agglomerator    faceAreaPair;
        mergeLevels     1;
    }
}
```

The fixed-iteration settings intentionally do not require the inner `dU` and
`dp` solves to converge. This follows the usual physics-preconditioner pattern:
do a small amount of cheap approximate physics work and let the outer Krylov
solver correct the remaining error.

### Useful Model Controls

The following entries are read from `newtonIcoFluidCoeffs`:

- `predictor`: if enabled, extrapolates the velocity from old time levels before
  the PETSc solve.
- `forceImplicitFlux`: if enabled, updates the flux from the current velocity
  inside the residual and Jacobian paths.
- `fluidFluxExtrapolationAlgorithm1`: selects one of the explicit old-time flux
  extrapolation formulae when `forceImplicitFlux` is disabled.
- `addDivPhiUDamping`: adds a semi-implicit damping term based on
  `div(phiAbs)`. The matrix uses `fvm::SuSp` so only the part that preserves
  the negative-diagonal momentum-block sign is added implicitly.
- `pressureScaleFactor`: scales the pressure-continuity residual and pressure
  Jacobian rows relative to the momentum rows.
- `maxTimeStepRetries`: maximum failed-PETSc-solve retry count when
  `adjustTimeStep` is enabled.
- `stopOnPetscError`: passed to `foamPetscSnesHelper` to control whether PETSc
  errors abort immediately.
- `tolerateSnesNonConvergence`: if `yes`, accept SNES outcomes of
  `SNES_DIVERGED_MAX_IT` and `SNES_DIVERGED_LINE_SEARCH` as usable iterates
  and advance time. Intended for Picard-style runs that perform a small
  fixed number of nonlinear sweeps per step. Genuine numerical failures
  (NaN, function-domain, KSP divergence, dtol) still trigger the
  time-step-reduction retry path. Defaults to `no`.
- `optionsFile`: optional PETSc options file name, defaulting to
  `petscOptions`.

`newtonIcoFluid` honours `relaxationFactors.equations.U` from
`system/fvSolution` by adding an implicit diagonal source to the
velocity `fvMatrix` before PETSc insertion:
`(1/alpha - 1)*fvm::Sp(UEqn.A(), U)`. This scales the segregated
velocity-equation diagonal from `D` to `D/alpha` without calling
`fvMatrix::relax()`, whose diagonal-dominance logic assumes a positive
diagonal and is not valid for this negative-diagonal Jacobian. For
implicit-flux Jacobians, the extra tensor terms from the Newton
linearisation of `div(phi,U)` are still inserted directly into PETSc
and are not modified by the scalar diagonal relaxation. The pressure
stabilisation `rAUf` is still constructed from the upwind
laplacian-ddt-div approximation before optional `addDivPhiUDamping`
terms are added, since the damping term does not guarantee the same
diagonal sign.

`relaxationFactors.equations.p` is not honoured. For step-level
damping (analogous to a global relaxation factor) use
`snes_linesearch_damping` in the PETSc options.

Time-step adaptation uses controls in `system/controlDict`:

- `adjustTimeStep`, `maxDeltaT`, and `minDeltaT`
- `minTargetNIter` and `maxTargetNIter`
- `logTimeStepAdjustments`
- `timeStepLogFile`

## Developer Notes

### Class Responsibilities

`newtonIcoFluid` derives from both `fluidModel` and `foamPetscSnesHelper`.
`fluidModel` supplies the OpenFOAM fields and fluid-model interface.
`foamPetscSnesHelper` owns the PETSc SNES/KSP integration, solution vector,
backup vector, matrix callbacks, and custom `physics` preconditioner hook.

The constructor:

- constructs the base `fluidModel` and `foamPetscSnesHelper` for the `Up`
  field;
- creates the transport and turbulence models;
- reads density, pressure reference data, `pressureScaleFactor`, and PETSc
  behavior flags;
- creates default momentum and pressure stabilisation dictionaries if the user
  did not provide them;
- supports the older stabilisation dictionary layout by treating it as pressure
  stabilisation;
- creates the `momentumStabilisation` and `pressureStabilisation` runtime
  models.

### Main Time-Step Flow

`evolve()` is the top-level solve entry point:

1. Correct `U` boundary conditions and report Courant number.
2. Save the old time, old mesh points, time state, and PETSc solution backup so
   failed solves can be retried.
3. Optionally apply the old-time velocity predictor.
4. Move the mesh and rebuild `phi`, including relative flux correction for
   moving meshes.
5. Call `foamPetscSnesHelper::solve(true)`.
6. Accept a converged solve, or a tolerated max-iteration/line-search outcome
   when `tolerateSnesNonConvergence` is enabled. On other PETSc failures,
   restore the solution and time state, then reduce `deltaT` if
   `adjustTimeStep` is enabled.
7. Extract the accepted PETSc solution vector back into `U` and `p`.
8. Rebuild `phi` and correct transport and turbulence models.

`setDeltaT()` implements the SNES-iteration-based time-step adjustment. It uses
the previous SNES iteration count and convergence reason to reduce or increase
`deltaT`, and can write a small log of the adjustment history.

`restoreOldTimeState()` restores mesh points and field old-time states when a
PETSc solve fails and the time step needs to be retried.

### PETSc Vector Layout

`initialiseSolution()` and `initialiseJacobian()` delegate to
`foamPetscSnesHelper` with a cell-centred block size:

- 2-D: `blockSize_ = 3`, storing `(Ux, Uy, p)`;
- 3-D: `blockSize_ = 4`, storing `(Ux, Uy, Uz, p)`.

The helper functions `ExtractFieldComponents` and `InsertFieldComponents` map
between OpenFOAM fields and the PETSc vector. Momentum components start at
offset `0`; pressure is at `blockSize_ - 1`.

### Residual Assembly

`formResidual(f, x, extrapolatedFlux)` computes the nonlinear residual for a
given PETSc solution vector:

- copies velocity and pressure from `x` into `U` and `p`;
- corrects boundary conditions and updates `gradU`, `gradp`, and `phi`;
- computes the momentum residual from diffusion, pressure gradient, transient,
  advection, and optional momentum stabilisation terms;
- computes the pressure residual from `-div(U)` and pressure stabilisation;
- applies the pressure reference cell and `pressureScaleFactor`;
- inserts the cell residuals into the PETSc residual vector.

The pressure stabilisation uses `rAUf = interpolate(1/A(U))` from an auxiliary
momentum equation. This couples the pressure stabilisation strength to the
current momentum diagonal.

### Jacobian and Preconditioning Matrix

`formJacobian(jac, x, extrapolatedFlux)` assembles the matrix used by PETSc as
the preconditioning matrix:

- inserts the momentum block from an `fvVectorMatrix`;
- inserts implicit advection contributions when the flux is not extrapolated;
- computes `rAUf` from the momentum diagonal using the sign convention expected
  by the pressure block;
- inserts the pressure-stabilisation matrix;
- adds the wider `div(interpolate(grad(p)))` pressure-stabilisation
  contribution for `RhieChow` and `diffStencilLaplacian` pressure
  stabilisation;
- inserts the continuity coupling from velocity to pressure;
- inserts the pressure-gradient coupling from pressure to momentum.

This matrix is an approximation to the nonlinear residual Jacobian. The usual
PETSc setup uses `snes_mf_operator`, so PETSc applies the nonlinear operator
matrix-free and uses this assembled matrix as the preconditioner.

### Model-Side Physics Preconditioner

`precondition(y, x)` implements the custom model-side path used by
`pc_type physics` with `pc_physics_type model`:

1. Create temporary fields `dU` and `dp`.
2. Extract the momentum and pressure components of the incoming PETSc vector
   into OpenFOAM source fields.
3. Assemble and solve an approximate momentum correction equation for `dU`
   using `mesh.solverDict("dU")`.
4. Build `rAUf` from the momentum diagonal.
5. Assemble and solve the pressure correction equation for `dp` using the
   pressure-stabilisation Jacobian and `mesh.solverDict("dp")`.
6. Correct `dU` with `rAU*grad(dp)`.
7. Insert `dU` and `dp` back into the PETSc output vector.

The preconditioner intentionally leaves some higher-order and tensor
contributions to the assembled PETSc matrix path; it is designed to be cheap
and useful as a physics-based approximation rather than a full exact inverse.

For consistency checks, compare this model PC against the assembled operator PC
path in `foamPetscSnesHelper` by using `pc_physics_type operator`. If the
assembled operator is solved exactly, Krylov iterations should approach one;
larger counts usually indicate matrix-free/assembled-operator inconsistency or
boundary/parallel-stencil gaps.

### Related Files

- `src/solids4FoamModels/numerics/foamPetscSnesHelper`
- `src/solids4FoamModels/numerics/stabilisationModels`
- `tutorials/fluids/cylinderInChannel/system/fvSolution.newtonIcoFluid`
- `tutorials/fluids/cylinderInChannel/system/fvSolution.newtonIcoFluid.physicsPC`
