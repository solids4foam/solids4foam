# cavityFlexibleBottom

Lid-driven cavity flow with a flexible bottom wall, solved using both
partitioned and quasi-monolithic FSI approaches.

## Case Variants

- **`base/partitioned/`**: Standard partitioned Dirichlet-Neumann coupling
  (Aitken acceleration).
- **`base/quasiMonolitic/`**: Quasi-monolithic Newton coupling using PETSc SNES.

## Quasi-Monolithic Case

The monolithic case assembles the coupled fluid-solid system into a single
PETSc nonlinear solve per time step. The 2x2 block system couples:
- Fluid: incompressible Navier-Stokes (velocity + pressure)
- Solid: nonlinear geometry total Lagrangian (displacement)

### Running

Serial:
```bash
cd base/quasiMonolitic
./Allrun
```

Parallel (2 processors):
```bash
cd base/quasiMonolitic
./Allrun parallel
```

### PETSc Options

The solver behaviour is controlled via the `optionsFile` entry in
`constant/fsiProperties`. Several configurations are provided:

| File | Description | Notes |
|------|-------------|-------|
| `petscOptions.mf_bjacobi_lu` | Matrix-free J*v, block Jacobi + LU sub-blocks | **Default**. Best parallel performance |
| `petscOptions.mf_asm_ilu` | Matrix-free J*v, additive Schwarz + ILU(1) | Good alternative, ~2x more KSP iters |
| `petscOptions.mf_bjacobi_ilu` | Matrix-free J*v, block Jacobi + ILU(0) | Lighter memory, ~3x more KSP iters |
| `petscOptions.mf_hypre` | Matrix-free J*v, Hypre BoomerAMG | Does not converge on this saddle-point system |
| `petscOptions.lu` | Assembled Jacobian, direct LU | Serial only. Most robust |

To switch, edit `constant/fsiProperties`:
```c++
optionsFile petscOptions.mf_bjacobi_lu;
```

### Mesh Refinement

Alternative `blockMeshDict` files are provided in `system/fluid/` and
`system/solid/` (suffixed `.1` through `.4`) for mesh convergence studies.
Copy the desired file over `blockMeshDict` before running.
