# newtonQuasiMonolithicCouplingInterface

Quasi-monolithic Newton-based fluid-structure interaction coupling, following
Algorithm 2 from Jaiman & Joshi (2022), Chapter 6.

## Overview

This interface assembles a monolithic system of equations coupling the fluid
(velocity U, pressure p) and solid (displacement D) fields into a single PETSc
SNES nonlinear solve per time step. The 2x2 block system takes the form:

```
| Aff  Afs | | x_fluid |   | f_fluid |
|          | |         | = |         |
| Asf  Ass | | x_solid |   | f_solid |
```

where:
- **Aff**: fluid momentum + pressure equations (assembled by `newtonIcoFluid`)
- **Ass**: solid momentum equations (assembled by `nonLinGeomTotalLagVelocitySolid`
  or similar)
- **Afs**: solid velocity at the FSI interface affecting the fluid mesh motion
- **Asf**: fluid pressure at the FSI interface affecting the solid traction

The system is stored as a PETSc `MatNest` and solved with SNES Newton line
search. In parallel, the MatNest is converted to a monolithic `MPIAIJ` matrix
for preconditioning.

## Key Settings (fsiProperties)

```c++
fluidSolidInterface NewtonQuasiMonolithic;

NewtonQuasiMonolithicCoeffs
{
    solidPatch          interface;
    fluidPatch          interface;
    fluidSystemScaleFactor 1e+08;    // Scale fluid equations for conditioning
    coupled             yes;
    interfaceTransferMethod directMap;
    passViscousStress   yes;
    optionsFile         petscOptions.mfjacobi;  // PETSc SNES/KSP/PC options
}
```

## PETSc Options Files

The `optionsFile` entry points to a file containing PETSc command-line options
for the SNES nonlinear solver, KSP linear solver, and PC preconditioner.

### Recommended for parallel runs

**`petscOptions.mfjacobi`** (default): Matrix-free Jacobian with Jacobi
preconditioner. Uses finite-difference Jacobian-vector products for the Newton
direction and the assembled Jacobian only for preconditioning.

**`petscOptions.lupar`**: Matrix-free Jacobian with redundant LU
preconditioner. More robust but gathers the matrix to a single process for
factoring — only practical for small to moderate problem sizes.

### Serial-only options

**`petscOptions.lu`**: Direct LU factorisation of the assembled Jacobian.
The most robust option for serial runs.

## Running in Parallel

```bash
./Allrun parallel
```

This decomposes both fluid and solid regions, runs `mpirun -np N solids4Foam
-parallel`, then reconstructs. The number of processors is set in
`system/decomposeParDict`.

## Files

- `newtonQuasiMonolithicCouplingInterface.C` - Main implementation
- `newtonQuasiMonolithicCouplingInterface.H` - Class declaration

## Dependencies

- PETSc (with SNES support)
- `foamPetscSnesHelper` (base class providing PETSc SNES framework)
- `newtonIcoFluid` (fluid model with Newton-linearised formJacobian/formResidual)
- A PETSc-enabled solid model (e.g. `nonLinGeomTotalLagVelocitySolid`)
