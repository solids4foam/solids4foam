---
sort: 2
---

# nonLinGeomTotalLagTotalDispSolid

This page documents the nonlinear total-Lagrangian total-displacement solid
model. The runtime type is:

```text
nonLinearGeometryTotalLagrangianTotalDisplacement
```

The model assumes finite strains and rotations and solves for the total
displacement field `D`. The governing equations are written in the total
Lagrangian frame, so the mesh is not moved during the solution procedure.

The stress is computed by the run-time selectable mechanical law.

---

## User Guide

### What this model solves

`nonLinGeomTotalLagTotalDispSolid` is the nonlinear-geometry analogue of the
total-displacement linear solver. It:

- solves the momentum equation in total-Lagrangian form;
- uses total displacement `D` as the primary unknown;
- computes the total deformation gradient `F`, its inverse `Finv`, and the
  Jacobian `J`;
- supports segregated implicit and PETSc SNES solution paths;
- can also solve pressure when requested, but only through the PETSc SNES
  path;
- uses the selected mechanical law to compute stress.

The model is appropriate when deformation is large enough that linear geometry
is not valid.

### Supported solution algorithms

The `solutionAlgorithm` entry in the coefficients sub-dictionary selects how
the equations are solved:

| Value | Meaning | Notes |
| --- | --- | --- |
| `implicitSegregated` | Default segregated solve | Recommended |
| `PETScSNES` | Nonlinear solve via PETSc SNES | Requires PETSc |
| `implicitCoupled` | Unsupported | Do not select this value |
| `explicit` | Unsupported | Do not select this value |

The runtime type name for this class is
`nonLinearGeometryTotalLagrangianTotalDisplacement`.

### Important inherited options

The following `solidModel` options are particularly relevant here:

| Entry | Default | Relevance |
| --- | --- | --- |
| `solutionAlgorithm` | `implicitSegregated` | Selects the solver path |
| `dampingCoeff` | `0` | Damping in the momentum equation |
| `predictor` | `false` | Linear predictor for `D` at new time steps |
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
solidModel        nonLinearGeometryTotalLagrangianTotalDisplacement;

nonLinearGeometryTotalLagrangianTotalDisplacementCoeffs
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

- `D` in the time directory;
- `solidProperties` in `constant`;
- `g` in `constant`;
- a compatible mechanical-law configuration for the chosen material model.

The base `solidModel` infrastructure also manages `DD`, `pointD`, `pointDD`,
`grad(D)`, `grad(DD)`, `F`, `Finv`, `J`, `sigma`, and related restart fields as
needed.

### Field glossary

The main fields you will typically see in the output are:

- `D`: total displacement field at cell centres.
- `DD`: displacement increment, defined as `DD = D - D.oldTime()`.
- `pointD`: total displacement interpolated to the mesh points.
- `pointDD`: point displacement increment, defined as
  `pointDD = pointD - pointD.oldTime()`.
- `F`: total deformation gradient.
- `Finv`: inverse of the total deformation gradient.
- `J`: Jacobian of `F`.
- `grad(D)`: cell-centre gradient of `D`.
- `grad(DD)`: cell-centre gradient of `DD`.
- `U`: velocity field, computed from the time derivative of `D`.
- `A`: acceleration field used by the predictor and SNES paths.
- `sigma`: symmetric stress tensor field.
- `p`: pressure field, only when `solvePressure` is enabled.

Here `oldTime()` means the value stored at the previous time step.

### Boundary conditions

The solver treats the following patch field types specially:

- `solidTraction`
- `fixedDisplacementZeroShear`
- `symmetry`
- `slip`

On traction patches, the force is formed using the current deformed normals and
area measures. On zero-shear, symmetry, and slip patches, the tangential force
is projected away.

The `solidTraction` patch also has an optional `useUndeformedArea()` switch.
This is not implemented for the total-Lagrangian model when it is asked to use
the undeformed area.

### Output fields

Typical written fields include:

- `D`, `DD`
- `pointD`, `pointDD`
- `grad(D)`, `grad(DD)`
- `F`, `Finv`, `J`
- `U`
- `sigma`
- linear-geometry-style strain and von Mises stress fields via
  `solidModel::writeFields()`

---

## Developer Notes

### Class role

`nonLinGeomTotalLagTotalDispSolid` is the total-Lagrangian nonlinear geometry
solver with total displacement as the primary unknown. The key design choices
are:

- geometry is evaluated in the reference configuration;
- `D` is the primary variable;
- `F`, `Finv`, and `J` are updated from `grad(D)`;
- stress is delegated to `mechanicalModel`;
- the solver supports segregated implicit and PETSc SNES paths only.

The class inherits from `solidModel` and `foamPetscSnesHelper`.

### Main control flow

#### Construction

The constructor performs the following setup:

1. Calls the base `solidModel` constructor.
2. Creates the PETSc helper if PETSc is available and SNES is selected.
3. Initialises `F`, `Finv`, `J`, `A`, `impK`, `impKf`, and `rImpK`.
4. Reads optional settings such as `predictor`, `solvePressure`,
   `dampingCoeff`, and the PETSc options file name.
5. Forces creation of old-time fields for consistent restart behaviour.
6. Creates `p` and validates that `solvePressure` is only used with
   `PETScSNES`.
7. Recomputes the deformation gradient and stress if PETSc SNES is active.
8. Enforces `leastSquaresS4f` for `grad(D)` when PETSc SNES is used.
9. Enables `extrapolateValue` on `solidTraction` patch fields for the PETSc
   path.

#### `evolve()`

`evolve()` dispatches to one of:

- `evolveImplicitSegregated()`
- `evolveSnes()`

Any other algorithm value is treated as unsupported for this class.

#### Segregated implicit path

`evolveImplicitSegregated()`:

- corrects `D` boundary conditions;
- optionally applies the linear predictor on the first call in a new time step;
- computes deformed face normals and face areas from `F` and `Finv`;
- assembles the nonlinear momentum equation in the reference configuration;
- adds stabilisation and optional damping;
- enforces traction boundary conditions with `enforceTractionBoundaries()`;
- applies any cell-displacement constraints through `solidModel::setCellDisps()`;
- solves, relaxes, and checks convergence;
- updates `DD`, `gradD`, `gradDD`, `pointD`, `pointDD`, `U`, `A`, `F`, `Finv`,
  and `J`.

#### PETSc SNES path

`evolveSnes()`:

- optionally predicts `D` and inserts the predicted values into the SNES
  solution vector;
- solves the nonlinear system through `foamPetscSnesHelper`;
- extracts the solution back into `D`;
- optionally extracts `p` when `solvePressure` is enabled;
- recomputes `gradD`, `pointD`, `DD`, `pointDD`, `U`, `A`, `F`, `Finv`, and
  `J`.

The PETSc path uses the same approximate Jacobian structure as the segregated
solver, which makes Jacobian-free Krylov-Newton methods a natural fit.

### Traction boundary enforcement

`enforceTractionBoundaries()` is responsible for mapping patch traction
conditions into face forces on the deformed surface. It:

- uses the current deformed normals and area magnitudes;
- honours `solidTraction` pressure and traction components;
- projects tangential force away for zero-shear, symmetry, and slip patches.

This is one of the key places where the total-Lagrangian formulation differs
from the linear-geometry model.

### Extension points

If you extend this model, the main places to consider are:

- `evolve()` for additional solution algorithms;
- `enforceTractionBoundaries()` for new traction-type patches;
- the constructor if new validation or default setup is required;

Keep the reference-configuration interpretation in mind when making changes:
the solver updates kinematic quantities relative to the undeformed state, not
the previous time step.
