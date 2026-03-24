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

| File | Description | Parallel |
|------|-------------|----------|
| `petscOptions.mfjacobi` | Matrix-free J*v, Jacobi PC, LGMRES | Yes (default) |
| `petscOptions.lupar` | Matrix-free J*v, redundant LU PC | Yes (small cases) |
| `petscOptions.lu` | Assembled Jacobian, direct LU | Serial only |

To switch, edit `constant/fsiProperties`:
```c++
optionsFile petscOptions.mfjacobi;
```

### Mesh Refinement

Alternative `blockMeshDict` files are provided in `system/fluid/` and
`system/solid/` (suffixed `.1` through `.4`) for mesh convergence studies.
Copy the desired file over `blockMeshDict` before running.
