# newtonQuasiMonolithicCouplingInterface

Quasi-monolithic Newton-based fluid-structure interaction coupling, following
Algorithm 2 from Jaiman & Joshi (2022), Chapter 6.

## Overview

This interface assembles a monolithic system of equations coupling the fluid
(velocity U, pressure p) and solid (velocity U) fields into a single PETSc SNES
nonlinear solve per time step. The 2x2 block system takes the form:

```text
| Aff  Afs | | x_fluid |   | f_fluid |
|          | |         | = |         |
| Asf  Ass | | x_solid |   | f_solid |
```

where:

- **Aff**: fluid momentum and pressure equations (assembled by `newtonIcoFluid`)
- **Ass**: solid momentum equations (assembled by
  `nonLinearGeometryTotalLagrangianVelocity` or `linearGeometryVelocity`)
- **Afs**: solid velocity at the FSI interface affecting the fluid
- **Asf**: fluid pressure at the FSI interface affecting the solid traction

The system is stored as a PETSc `MatNest` and solved with SNES Newton line
search. In parallel, the MatNest is converted to a monolithic `MPIAIJ` matrix
for preconditioning.

The solid sets its displacement from its velocity with the backward (BDF2)
scheme, and evaluates the stress and the interface traction at the new time.
Before each solve, the fluid mesh is moved to a predicted interface
displacement; the fluid mesh motion is a solid model
(`motionSolver solids4foamSolidModel` in `constant/fluid/dynamicMeshDict`).

The coupling is available on OpenFOAM.com with PETSc only.

## Key Settings (fsiProperties)

```c++
fluidSolidInterface NewtonQuasiMonolithic;

NewtonQuasiMonolithicCoeffs
{
    solidPatch          interface;
    fluidPatch          interface;
    coupled             yes;

    // directMap (conformal) or AMI (non-conformal)
    interfaceTransferMethod directMap;

    // Pass the viscous traction to the solid, in addition to the pressure
    passViscousStress   yes;

    // Scale the fluid equations to improve the condition number (default 1)
    fluidSystemScaleFactor 1e+08;

    // Optional: scale the solid equations (default 1/(mu*deltaT), averaged
    // over the solid)
    //solidSystemScaleFactor 1;

    // Optional: see Non-conformal Interfaces
    //conservativeTractionTransfer no;
    //interfaceFluxDiagnostics no;

    // Optional: stop when PETSc reports an error (default yes)
    //stopOnPetscError yes;

    // Optional: the maximum number of retries of a failed time step (default
    // 10; see Time-step Adjustment)
    //maxTimeStepRetries 10;

    // Optional: a file of PETSc options, used if system/fvSolution has no
    // solvers/UpU entry (see PETSc Options)
    //optionsFile         petscOptions;
}
```

## Interface Conditions

On the fluid interface patch:

- **Velocity**: `fixedValue`. The coupling sets the value to the solid
  interface velocity.
- **Pressure**: `fixedGradient`. The coupling sets the gradient from the normal
  momentum equation at the moving wall, which makes the interface pressure and
  the added mass second-order accurate in space. A `zeroGradient` interface
  pressure also runs, but is first order, and the coupling gives a warning.
- **Viscous traction** (`passViscousStress yes`): the boundary flux of the
  fluid momentum residual's viscous terms, `rho*(n & -devReff)` including the
  transpose part of the viscous stress, which is nonzero on a deforming wall
  (see the `newtonIcoFluid` README, `viscousStressTranspose`). The coupling
  blocks of the assembled Jacobian approximate it through the boundary-face
  `snGrad` of the adjacent fluid cell velocity; the matrix-free residual
  contains the full form.

On the solid interface patch, use a `solidTraction` condition, as for the
partitioned couplings.

## Non-conformal Interfaces

`interfaceTransferMethod` selects how interface data are transferred between
the fluid and solid interface faces:

- `directMap`: conformal interfaces, where each fluid face matches one solid
  face.
- `AMI`: non-conformal interfaces, where the fluid and solid faces differ. Each
  fluid face value is the face-area-weighted sum of the overlapping solid face
  values, and vice versa.

The same face weights are used for every transfer in the residual (interface
velocity and displacement, traction) and for the coupling blocks of the
Jacobian, so the Newton convergence is the same as for a conformal interface.
Other methods (`RBF`) do not provide face weights and are not supported. The
point displacement that moves the fluid mesh uses the mapping's point transfer.

By default, the solid traction is interpolated from the fluid traction, as in
the partitioned couplings. The total interface force then balances only up to
the discretisation error of the interface geometry. Setting

```c++
    conservativeTractionTransfer yes;
```

instead distributes the fluid face forces to the solid faces with the transpose
of the solid-to-fluid velocity transfer, so the total interface force and power
balance to round-off. The solid applies the resulting nominal traction on its
undeformed face areas, so a nonlinear-geometry solid must set
`useUndeformedArea yes` on its interface `solidTraction` patch.

With `interfaceFluxDiagnostics yes`, the log reports the interface force
balance at each time step (`Interface force balance: ... rel_imbalance=`).

`ringAddedMass/verification` (`./Allverify --study nonconformal`) verifies the
coupling with solid interface faces coarser and finer than the fluid's against
the exact added-mass frequency.

## PETSc Options

The PETSc options for the monolithic `UpU` system are read from the
`solvers/UpU` entry of the case-level `system/fvSolution`. If that entry does
not exist, they are read from the file named by `optionsFile` in
`fsiProperties` (default `petscOptions`), which contains PETSc command-line
options.

The monolithic tutorials all use the same set-up in `system/fvSolution`:

- **Nonlinear solver**: Newton line search (`snes_type newtonls`) with a
  matrix-free Jacobian (`snes_mf_operator`); the assembled Jacobian is used to
  build the preconditioner only.
- **Linear solver**: `ksp_type lgmres`.
- **Preconditioner**: a multiplicative `fieldsplit` between the fluid and
  solid blocks:
  - solid: exact LU on one rank (`redundant` with `lu`);
  - fluid: a nested velocity-pressure Schur `fieldsplit`, with block Jacobi
    and ILU(0) on the velocity block and HYPRE BoomerAMG on the pressure Schur
    complement.

The fluid `fieldsplit` block size and fields depend on the dimension:

```text
// 2-D: velocity components 0,1 and pressure 2
fieldsplit_fluid_pc_fieldsplit_block_size "3";
fieldsplit_fluid_pc_fieldsplit_0_fields "0,1";
fieldsplit_fluid_pc_fieldsplit_1_fields "2";

// 3-D: velocity components 0,1,2 and pressure 3
fieldsplit_fluid_pc_fieldsplit_block_size "4";
fieldsplit_fluid_pc_fieldsplit_0_fields "0,1,2";
fieldsplit_fluid_pc_fieldsplit_1_fields "3";
```

Copy the `solvers/UpU` entry of a tutorial with the same dimension, e.g.
`ringAddedMass` (2-D) or `3dTube` (3-D), as a starting point.

`3dTube` and `beamInCrossFlow` change this set-up as follows, which makes
their serial runs about 2.3 (`beamInCrossFlow`) and 2.9 (`3dTube`) times
faster, with the same results to within the solver tolerances:

- `ksp_rtol "1e-3"`: the first Newton iteration of a time step is limited
  by the nonlinearity of the residual to a reduction of about 1e-3, so a
  tighter linear tolerance only adds Krylov iterations, each of which costs
  a residual evaluation with `snes_mf_operator`;
- `snes_lag_jacobian "8"` and `snes_lag_preconditioner "8"`, with
  `snes_lag_jacobian_persists` and `snes_lag_preconditioner_persists`:
  with `snes_mf_operator`, the assembled Jacobian only builds the
  preconditioner, which can be reused across Newton iterations and time
  steps;
- `fieldsplit_fluid_pc_fieldsplit_schur_fact_type full` and
  `snes_linesearch_type basic`;
- in the fluid mesh motion (`D` in `system/fluid/fvSolution`), HYPRE
  BoomerAMG instead of an LU factorisation, rebuilt every 20 time steps.

Replacing `snes_mf_operator` with the assembled Jacobian is not
recommended: the assembled Jacobian is approximate, and Newton then
converges only linearly.

The physics-based preconditioner (`pc_type physics`) is work in progress: on
`ringAddedMass`, Newton stalls at the first coupled time step with it.

## Time-step Adjustment

With `adjustTimeStep yes` in `system/controlDict`, which then also needs
`minDeltaT` and `maxDeltaT`, the time step is adapted to the number of Newton
iterations: it grows when fewer than `minTargetNIter - 1` are needed and
shrinks when more than `maxTargetNIter + 1` are needed (defaults 3 and 6), by
at most a factor of 1.5 and 2. If the Newton solve fails, the time step is
restored and retried with a quarter of the time step, at most
`maxTimeStepRetries` times (in `fsiProperties`, default 10). Without
`adjustTimeStep`, a failed Newton solve stops the run. The adjustments are
logged to `timeStepLog.dat` (`timeStepLogFile`, `logTimeStepAdjustments`).

## Running

The monolithic tutorial variants are run with `./Allrun monolithic`, and most
of them in parallel with `./Allrun monolithic parallel`.

In parallel, the fluid and solid regions can be decomposed independently with
`decomposePar`, or together with `decomposeParMonolithic`, which can give some
ranks no cells of one region.

## Files

- `newtonQuasiMonolithicCouplingInterface.C` - Main implementation
- `newtonQuasiMonolithicCouplingInterface.H` - Class declaration

## Dependencies

- PETSc (with SNES support)
- `foamPetscSnesHelper` (PETSc SNES framework)
- `newtonIcoFluid` (fluid model with Newton-linearised formJacobian/formResidual)
- A velocity-based PETSc solid model (`nonLinearGeometryTotalLagrangianVelocity`
  or `linearGeometryVelocity`)
