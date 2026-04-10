# foilInWind

## Case description

A flexible airfoil clamped at its base and exposed to a uniform wind flow. The
foil is thin and compliant (rubber-like), so the wind loading causes significant
bending deflection. The case exercises large-displacement nonlinear structural
response coupled to incompressible laminar fluid flow.

The geometry is three-dimensional, extruded in the span direction (8 cell
layers). The fluid and solid meshes share a common interface (`flag`/`interface`
patch pair) across which displacement and traction are exchanged.

### Physical parameters

| Property | Value |
|----------|-------|
| Fluid kinematic viscosity | 1.542e-5 m²/s |
| Fluid density | 1.18 kg/m³ |
| Solid density | 2000 kg/m³ |
| Young's modulus | 2 MPa |
| Poisson's ratio | 0.35 |
| Material model | St. Venant-Kirchhoff nonlinear elastic |

### Simulation parameters

| Parameter | Value |
|-----------|-------|
| End time | 12 s |
| Time step | 0.01 s |
| Total steps | 1200 |

---

## Running the case

```bash
./Allrun                      # partitioned (IQNILS), serial
./Allrun aitken               # partitioned (Aitken), serial
./Allrun monolithic           # monolithic (Newton/PETSc), serial
./Allrun monolithic parallel  # monolithic (Newton/PETSc), parallel
./Allrun aitken parallel      # partitioned (Aitken), parallel
```

The script:

1. Restores the `0` directory from `0.orig`
2. Symlinks the appropriate `*.monolithic` or `*.partitioned` config files
3. Runs `blockMesh` for fluid and solid regions
4. Refines the fluid mesh with `setSet` + `refineMesh`
5. Initialises the fluid velocity with `potentialFoam`
6. For monolithic parallel runs, calls `decomposeParMonolithic -force`
7. For partitioned parallel runs, calls `decomposePar` independently for the
   fluid and solid regions
8. Runs `solids4Foam` (serial or parallel)
9. Generates deflection and force plots with gnuplot (if available)

The monolithic approach requires OpenFOAM.com (ESI) and a PETSc installation
with `PETSC_DIR` set. It does not run with foam-extend or OpenFOAM.org.

If `decomposeParMonolithic` has not been built yet, rebuild the applications
before running the parallel monolithic case, for example:

```bash
cd applications/utilities/decomposeParMonolithic
wmake
```

---

## FSI settings

### Coupling approach selection

The coupling approach is controlled by `constant/fsiProperties`, which is
symlinked by `Allrun` to the appropriate variant:

| File | Approach |
|------|----------|
| `fsiProperties.iqnils` | Partitioned, IQNILS acceleration |
| `fsiProperties.aitken` | Partitioned, Aitken relaxation |
| `fsiProperties.monolithic` | Monolithic Newton/PETSc |

### Partitioned approach (`fsiProperties.iqnils`)

```c++
fluidSolidInterface    IQNILS;

IQNILSCoeffs
{
    solidPatch          interface;
    fluidPatch          flag;

    coupled             no;
    couplingStartTime   2000.0;   // set < endTime to activate coupling

    relaxationFactor    0.1;
    outerCorrTolerance  1e-6;
    nOuterCorr          50;

    interfaceTransferMethod  directMap;
}
```

The `coupled no` / `couplingStartTime 2000` setting keeps the fluid and solid
decoupled throughout a standard run. To activate FSI coupling, set
`couplingStartTime` to a value within the simulation window (e.g. `2.0`).

### Monolithic approach (`fsiProperties.monolithic`)

```c++
fluidSolidInterface    NewtonQuasiMonolithic;

NewtonQuasiMonolithicCoeffs
{
    solidPatch          interface;
    fluidPatch          flag;

    coupled             no;
    couplingStartTime   2.0;       // fluid runs alone for 2 s; FSI from 2 s

    interfaceTransferMethod  directMap;

    passViscousStress   yes;       // include viscous traction in solid loading

    fluidSystemScaleFactor  1e8;   // rescales fluid rows to balance condition number
}
```

The `fluidSystemScaleFactor` multiplies all fluid rows and columns in the
monolithic matrix. The fluid pressure equation has much smaller coefficients
than the solid momentum equation, so scaling brings them to a comparable
magnitude and prevents the fluid block from dominating the preconditioner in
the wrong direction.

The fluid model is `newtonIcoFluid` (see
`constant/fluid/fluidProperties.monolithic`). The solid model is
`nonLinearGeometryTotalLagrangianVelocity` (see
`constant/solid/solidProperties.monolithic`).

---

## Monolithic decomposition

The monolithic parallel branch no longer decomposes the two regions
independently. Instead, `Allrun monolithic parallel` uses the custom
`decomposeParMonolithic` utility.

The active decomposition dictionary is `system/decomposeParDict`, which is
symlinked by `Allrun` to `system/decomposeParDict.monolithic` in monolithic
mode. For this case it contains:

```foam
numberOfSubdomains 4;

method scotch;

monolithicCoeffs
{
    regions (fluid solid);

    regionWeights
    {
        fluid 4;
        solid 3;
    }

    interfaces
    (
        {
            regionA fluid;
            patchA  flag;
            regionB solid;
            patchB  interface;
        }
    );
}
```

### What the utility does

- Loads the `fluid` and `solid` meshes and decomposes them as one combined
  graph.
- Uses the internal mesh connectivity within each region plus extra graph edges
  across the `flag`/`interface` FSI patch pair.
- Applies the optional `regionWeights` as per-cell weights in the combined
  graph.
- Runs the selected decomposition method once, then writes the normal
  `processor*/fluid/...` and `processor*/solid/...` directories by calling the
  standard `decomposePar` utility internally with `method manual`.

On this case, the resulting decomposition typically leaves some ranks with zero
solid cells, which is the intended behaviour for a strongly fluid-dominated
monolithic FSI problem.

### Current v1 assumptions

- The coupled patches must be conformal and have the same number of faces.
- The current implementation matches the two interface patches by nearest face
  centre, which is appropriate for this case's conformal `flag`/`interface`
  pair.
- No separate mapping mode is read from `monolithicCoeffs`; the utility always
  uses this conformal nearest-face-centre matching.
- `Allrun` uses the normal full workflow. The optional
  `decomposeParMonolithic -decompose-only` mode is mainly for debugging.

If you do use `-decompose-only`, make sure that any follow-up manual
`decomposeParDict` uses the same `numberOfSubdomains` as
`system/decomposeParDict.monolithic`.

---

## Monolithic PETSc solver settings (`system/fvSolution.monolithic`)

The monolithic system is a block 5×5 matrix (per cell): fluid [Ux, Uy, Uz, p]
plus solid [U]. PETSc SNES drives the outer nonlinear iteration; PETSc KSP +
preconditioner handle each linear solve.

### Nonlinear solver (SNES)

```c++
snes_type              newtonls;
snes_linesearch_type   l2;       // l2-norm line search
snes_rtol              "1e-4";
snes_stol              "1e-4";
snes_mf_operator;                // matrix-free Jacobian-vector products
```

`snes_mf_operator` means PETSc applies J·v via finite-differenced residual
evaluations rather than assembling the full Jacobian explicitly. The
preconditioner still uses the assembled approximate Jacobian (the MatNest
structure from solids4foam).

### Linear solver (KSP)

```c++
ksp_type          lgmres;
ksp_gmres_restart "200";
ksp_rtol          "1e-4";
ksp_max_it        "1000";
```

### Preconditioner structure

The preconditioner is a two-level block split:

```text
PC = multiplicative fieldsplit( fluid | solid )
        fluid  --> nested fieldsplit Schur( velocity | pressure )
        solid  --> redundant + LU
```

#### Top level: fluid/solid multiplicative split

```c++
pc_type                   fieldsplit;
pc_fieldsplit_type        multiplicative;
```

The fluid and solid sub-problems are solved sequentially, with each pass using
the updated solution from the other. This is more effective than an additive
(Jacobi-style) split for strongly coupled FSI.

#### Solid sub-block: `redundant+lu`

```c++
fieldsplit_solid_ksp_type          preonly;
fieldsplit_solid_pc_type           redundant;
fieldsplit_solid_redundant_pc_type lu;
```

**Why `redundant` rather than `bjacobi`?**

The naive parallel choice for a direct solve is `bjacobi+lu`: each MPI rank
applies LU to its own local rows and ignores coupling to other ranks. For a
small solid mesh (the foilInWind solid has ~14,000 DOFs) this degrades badly
as the rank count increases:

| Ranks | DOFs per rank with `bjacobi+lu` | Max sigmaEq |
|-------|--------------------------------|-------------|
| 1 | 14,000 (full mesh, exact) | 4.95 (correct) |
| 8 | ~1,750 (1/8 of mesh, severe coupling loss) | **2.12 (wrong)** |
| 32 | ~437 (1/32 of mesh, almost no coupling) | **1.08 (wrong)** |

Solid mechanics is globally coupled — stiffness at one node drives response
throughout the mesh. A per-rank LU discards all cross-rank coupling. At 8+
ranks the preconditioner becomes so poor that LGMRES nominally converges to
`ksp_rtol` of the preconditioned residual while the actual Newton step is in
the wrong direction. The SNES outer loop then converges to the wrong nonlinear
fixed point: force (a fluid-dominated quantity) stays correct, but solid
displacement and stress are qualitatively wrong.

`PCREDUNDANT` fixes this by:

1. Gathering the full distributed solid matrix to rank 0
2. Applying exact LU to the complete assembled matrix on rank 0
3. Scattering the solution back to all ranks

The solid block is always factored as a 14K×14K system regardless of how
many ranks are in use. The result is correct at all rank counts.

The trade-off is that rank 0 holds the full solid matrix, so this approach
does not scale to very large solid meshes. For tutorial-sized meshes it is
cheap and exact.

#### Fluid sub-block: nested velocity/pressure Schur

```c++
fieldsplit_fluid_ksp_type                      preonly;
fieldsplit_fluid_pc_type                       fieldsplit;
fieldsplit_fluid_pc_fieldsplit_type            schur;
fieldsplit_fluid_pc_fieldsplit_schur_fact_type upper;     // upper triangular
fieldsplit_fluid_pc_fieldsplit_schur_precondition selfp;  // pressure mass approx
fieldsplit_fluid_pc_fieldsplit_block_size      "4";       // 3-D: [Ux,Uy,Uz,p]
fieldsplit_fluid_pc_fieldsplit_0_fields        "0,1,2";   // velocity
fieldsplit_fluid_pc_fieldsplit_1_fields        "3";       // pressure
```

The `upper` Schur factorization gives 24% fewer KSP iterations than `lower`
on this case. The difference is small on easier 2-D cases but grows with
problem difficulty; `upper` is therefore the better general default.

**Velocity sub-block:**

```c++
fieldsplit_fluid_fieldsplit_0_pc_type          bjacobi;
fieldsplit_fluid_fieldsplit_0_sub_pc_type      ilu;
fieldsplit_fluid_fieldsplit_0_sub_pc_factor_levels "0";  // ILU(0)
```

Per-rank ILU(0) is cheap per iteration. Despite requiring slightly more KSP
iterations than a full AMG velocity solve, the lower cost per iteration makes
it 30% faster in wall-clock time on this case.

**Pressure Schur complement:**

```c++
fieldsplit_fluid_fieldsplit_1_pc_type                           hypre;
fieldsplit_fluid_fieldsplit_1_pc_hypre_type                     boomeramg;
fieldsplit_fluid_fieldsplit_1_pc_hypre_boomeramg_strong_threshold "0.6";
fieldsplit_fluid_fieldsplit_1_pc_hypre_boomeramg_coarsen_type   HMIS;
fieldsplit_fluid_fieldsplit_1_pc_hypre_boomeramg_interp_type    ext+i;
```

HYPRE BoomerAMG with aggressive coarsening (`strong_threshold 0.6`, HMIS
coarsening, extended+i interpolation) handles the pressure Schur complement.
This configuration was validated across all five monolithic tutorial cases.

### Benchmark summary (one coupled step, local serial)

| Config | Newton | KSP total | sigmaEq | Wall time |
|--------|--------|-----------|---------|-----------|
| `bjacobi+lu` solid, `lower` Schur | 5 | 378 | 4.95 | ~200 s |
| `redundant+lu` solid, `upper` Schur | 5 | 321 | 4.998 | **85 s** |

---

## Output

After a successful run, `gnuplot` generates:

- `deflection.pdf` — foil tip deflection vs time
- `force.pdf` — aerodynamic force vs time

These can be used to compare partitioned and monolithic results.
