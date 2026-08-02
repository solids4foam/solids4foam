---
sort: 3
---

# nonLinGeomUpdatedLagSolid

This page documents the nonlinear updated-Lagrangian incremental-displacement
solid model. The runtime type is:

```text
nonLinearGeometryUpdatedLagrangian
```

The model assumes finite strains and rotations and solves for the increment of
displacement `DD`. The mesh is moved each time step so that the governing
equations are written in the updated configuration.

The stress is computed by the run-time selectable mechanical law.

---

## User Guide

### What this model solves

`nonLinGeomUpdatedLagSolid` is the nonlinear-geometry analogue of the updated
Lagrangian solid solver. It:

- solves the momentum equation in updated-Lagrangian form;
- uses the displacement increment `DD` as the primary unknown;
- updates the total displacement `D` from `D.oldTime() + DD`;
- computes the total deformation gradient `F`, relative deformation gradient
  `relF`, and their inverses/Jacobians;
- moves the mesh each time step;
- supports segregated implicit and PETSc SNES solution paths;
- can also solve pressure when requested, but only through the PETSc SNES
  path;
- uses the selected mechanical law to compute stress.

The model is appropriate when deformation is large enough that linear geometry
is not valid and mesh motion between time steps is part of the formulation.

### Supported solution algorithms

The `solutionAlgorithm` entry in the coefficients sub-dictionary selects how
the equations are solved:

| Value | Meaning | Notes |
| --- | --- | --- |
| `implicitSegregated` | Default segregated solve | Recommended |
| `PETScSNES` | Nonlinear solve via PETSc SNES | Requires PETSc |
| `implicitCoupled` | Unsupported | Do not select this value |
| `explicit` | Unsupported | Do not select this value |

The runtime type name for this class is `nonLinearGeometryUpdatedLagrangian`.

### Important inherited options

The following `solidModel` options are particularly relevant here:

| Entry | Default | Relevance |
| --- | --- | --- |
| `solutionAlgorithm` | `implicitSegregated` | Selects the solver path |
| `dampingCoeff` | `0` | Damping in the momentum equation |
| `predictor` | `false` | Linear predictor for `DD` at new time steps |
| `optionsFile` | `petscOptions` | PETSc options file |
| `stopOnPetscError` | `true` | Stop if PETSc reports an error |
| `relaxationMethod` | `fixed` | Under-relaxation method |
| `nCorrectors` | `10000` | Maximum number of outer correctors |
| `solutionTolerance` | `1e-06` | Primary convergence tolerance |
| `alternativeTolerance` | `1e-07` | Secondary convergence tolerance |
| `materialTolerance` | `1e-05` | Mechanical-law convergence tolerance |
| `infoFrequency` | `100` | Frequency for solver progress output |
| `restart` | `false` | Writes extra fields needed for a consistent restart |
| `writeResidualField` | `false` | Writes a residual field during output |
| `residualFile` | `false` | Writes `residual.dat` in the case directory |
| `stabilisation` | auto-created if absent | Has `momentum` and `pressure` |
| `cellDisplacements` | optional | Internal-cell displacement constraints |
| `solvePressure` | `false` | Pressure unknown, only with `PETScSNES` |

For `stabilisation`, the code creates a default dictionary if none is
provided. The default sub-model is `diffStencilLaplacian` with
`scaleFactor 0.1` for both momentum and pressure.

### Pressure coupling special case

This model can also solve for pressure by enabling:

```text
solvePressure true;
```

This is a special case and is only supported when `solutionAlgorithm` is
`PETScSNES`. When pressure solving is enabled:

- the PETSc solution vector includes `p` as an extra unknown;
- the pressure stabilisation sub-model is used;
- the code enforces the `PETScSNES` path at construction time.

### Recommended dictionary setup

Minimal example for `constant/solidProperties`:

```text
solidModel        nonLinearGeometryUpdatedLagrangian;

nonLinearGeometryUpdatedLagrangianCoeffs
{
    solutionAlgorithm    implicitSegregated;

    predictor            false;
    dampingCoeff         [0 0 -1 0 0 0 0] 0;
    relaxationMethod     fixed;

    nCorrectors          10000;
    solutionTolerance    1e-06;
    alternativeTolerance 1e-07;
    materialTolerance    1e-05;
    infoFrequency        100;

    restart              false;
    residualFile         false;
    writeResidualField   false;

    solvePressure        false;

    stabilisation
    {
        momentum
        {
            type        diffStencilLaplacian;
            scaleFactor 0.1;
        }

        pressure
        {
            type        diffStencilLaplacian;
            scaleFactor 0.1;
        }
    }
}
```

If `solvePressure` is enabled, add:

```text
    solvePressure         true;
```

If `solutionAlgorithm` is `PETScSNES`, also ensure the case provides the PETSc
options file expected by `optionsFile` (default name: `petscOptions`), unless
you override that name in the coefficients sub-dictionary.

### Required fields and boundary conditions

At minimum, the model expects:

- `DD` in the time directory;
- `solidProperties` in `constant`;
- `g` in `constant`;
- a compatible mechanical-law configuration for the chosen material model.

The base `solidModel` infrastructure also manages `D`, `pointD`, `pointDD`,
`grad(D)`, `grad(DD)`, `F`, `relF`, `relFinv`, `relJ`, `J`, `sigma`, and
related restart fields as needed.

### Field glossary

The main fields you will typically see in the output are:

- `DD`: displacement increment, the primary unknown in this model.
- `D`: total displacement, reconstructed as `D = D.oldTime() + DD`.
- `pointDD`: point displacement increment, interpolated from `DD`.
- `pointD`: total point displacement, updated from the old-time value plus
  `pointDD`.
- `relF`: relative deformation gradient for the current time step.
- `relFinv`: inverse of `relF`.
- `relJ`: Jacobian of `relF`.
- `F`: total deformation gradient.
- `J`: Jacobian of `F`.
- `grad(DD)`: cell-centre gradient of `DD`.
- `grad(D)`: cell-centre gradient of `D`.
- `U`: velocity field, computed from the time derivative of `D`.
- `A`: acceleration field.
- `sigma`: symmetric stress tensor field.
- `rho`: updated-configuration density field.
- `p`: pressure field, only when `solvePressure` is enabled.

Here `oldTime()` means the value stored at the previous time step.

### Boundary conditions

The solver treats the following patch field types specially:

- `solidTraction`
- `fixedDisplacementZeroShear`
- `symmetry`
- `slip`

On traction patches, the force is formed using the current deformed normals and
area magnitudes. On zero-shear, symmetry, and slip patches, the tangential
force is projected away.

`solidTraction::useUndeformedArea()` is not implemented for the updated-
Lagrangian model.

### Output fields

Typical written fields include:

- `D`, `DD`
- `pointD`, `pointDD`
- `grad(D)`, `grad(DD)`
- `F`, `relF`, `relFinv`, `relJ`, `J`
- `U`
- `A`
- `sigma`
- `rho`
- linear-geometry-style strain and von Mises stress fields via
  `solidModel::writeFields()`

---

## Developer Notes

### Class role

`nonLinGeomUpdatedLagSolid` is the updated-Lagrangian nonlinear geometry
solver with incremental displacement as the primary unknown. The key design
choices are:

- geometry is updated each time step;
- `DD` is the primary variable;
- `D` is reconstructed from `DD`;
- `F`, `relF`, `relFinv`, `relJ`, and `J` are updated from `grad(DD)`;
- the model is marked as incremental and moving-mesh;
- stress is delegated to `mechanicalModel`;
- the solver supports segregated implicit and PETSc SNES paths only.

The class inherits from `solidModel` and `foamPetscSnesHelper`.

### Main control flow

#### Construction

The constructor performs the following setup:

1. Calls the base `solidModel` constructor.
2. Creates the PETSc helper if PETSc is available and SNES is selected.
3. Initialises `F`, `J`, `relF`, `relFinv`, `relJ`, `rho`, `A`, `impK`,
   `impKf`, and `rImpK`.
4. Reads optional settings such as `predictor`, `solvePressure`,
   `dampingCoeff`, and the PETSc options file name.
5. Forces creation of old-time fields for consistent restart behaviour.
6. Creates `p` and validates that `solvePressure` is only used with
   `PETScSNES`.
7. Rebuilds `relF`, `relFinv`, `relJ`, `F`, and `J` when `restart` is enabled.
8. Enforces `leastSquaresS4f` for `grad(DD)` when PETSc SNES is used.
9. Enables `extrapolateValue` on `solidTraction` patch fields for the PETSc
   path.

#### `evolve()`

`evolve()` dispatches to one of:

- `evolveImplicitSegregated()`
- `evolveSnes()`

Any other algorithm value is treated as unsupported for this class.

#### Segregated implicit path

`evolveImplicitSegregated()`:

- corrects `DD` boundary conditions;
- optionally applies the linear predictor on the first call in a new time step;
- computes deformed face normals and face areas from the relative kinematics;
- assembles the updated-Lagrangian momentum equation for `DD`;
- adds stabilisation and optional damping;
- enforces traction boundary conditions with `enforceTractionBoundaries()`;
- applies any cell-displacement constraints through `solidModel::setCellDisps()`;
- solves, relaxes, and checks convergence;
- updates `D`, `DD`, `gradDD`, `pointDD`, `pointD`, `U`, `A`, `F`, `relF`,
  `relFinv`, `relJ`, and `J`.

#### PETSc SNES path

`evolveSnes()`:

- optionally predicts `DD` and inserts the predicted values into the SNES
  solution vector;
- solves the nonlinear system through `foamPetscSnesHelper`;
- extracts the solution back into `DD`;
- optionally extracts `p` when `solvePressure` is enabled;
- recomputes `D`, `gradDD`, `pointDD`, `pointD`, `U`, `A`, `F`, `relF`,
  `relFinv`, `relJ`, and `J`.

The PETSc path uses the same approximate Jacobian structure as the segregated
solver, which makes Jacobian-free Krylov-Newton methods a natural fit.

### Traction boundary enforcement

`enforceTractionBoundaries()` maps traction conditions into face forces on the
updated configuration. It:

- uses the current deformed normals and area magnitudes;
- honours `solidTraction` pressure and traction components;
- projects tangential force away for zero-shear, symmetry, and slip patches.

This is one of the key places where the updated-Lagrangian formulation differs
from the total-Lagrangian model.

### Extension points

If you extend this model, the main places to consider are:

- `evolve()` for additional solution algorithms;
- `enforceTractionBoundaries()` for new traction-type patches;
- the constructor if new validation or default setup is required;
- `updateTotalFields()` if the accumulated output fields need to be extended.

Keep the updated-configuration interpretation in mind when making changes:
the solver updates kinematic quantities relative to the last time step and
moves the mesh each step.
