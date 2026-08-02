---
sort: 1
---

# linGeomTotalDispSolid

This page documents the linear-geometry total-displacement solid model.
The runtime type is:

```text
linearGeometryTotalDisplacement
```

The model assumes small strains and small rotations, and solves for the total
displacement field `D`. The stress is computed by the run-time selectable
mechanical law, so the model is not tied to one specific constitutive law.

---

## User Guide

### What this model solves

`linGeomTotalDispSolid` is the standard small-deformation solid solver in
solids4foam. It:

- solves the linear momentum equation in linear-geometry form;
- uses total displacement `D` as the primary unknown;
- computes displacement increments `DD`, point displacement `pointD`, point
  displacement increment `pointDD`, velocity `U`, gradient fields, and stress
  fields as part of the standard solidModel workflow;
- uses the selected mechanical law to calculate stress;
- supports segregated implicit, PETSc SNES, and explicit time integration.

The model is suitable when changes in geometry are small enough that the mesh
can be treated as undeformed for the governing equations.

### Supported solution algorithms

The `solutionAlgorithm` entry in the coefficients sub-dictionary selects how
the equations are solved:

| Value | Meaning | Notes |
| --- | --- | --- |
| `implicitSegregated` | Default segregated solve | Recommended |
| `PETScSNES` | Nonlinear solve via PETSc SNES | Requires PETSc |
| `explicit` | Explicit time integration | Time-step limited by stability |
| `implicitCoupled` | Not implemented here | Do not select this value |

The runtime type name for this class is `linearGeometryTotalDisplacement`.

### Important inherited options

The following `solidModel` options are particularly relevant here:

| Entry | Default | Relevance |
| --- | --- | --- |
| `solutionAlgorithm` | `implicitSegregated` | Selects the solver path |
| `dampingCoeff` | `0` | Damping for momentum and explicit updates |
| `predictor` | `false` | Linear predictor for `D` at new time steps |
| `optionsFile` | `petscOptions` | PETSc options file name |
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
| `stabilisation` | auto-created if absent | `momentum` always used |
| `cellDisplacements` | optional | Internal-cell displacement constraints |
| `solvePressure` | `false` | Pressure unknown, only with `PETScSNES` |

For `stabilisation`, the code will create a default dictionary if none is
provided. The `momentum` sub-dictionary is always used. The `pressure`
sub-dictionary is only used when `solvePressure` is enabled. The default
sub-model is `diffStencilLaplacian` with `scaleFactor 0.1`. The old flat
format is not accepted.

### Pressure coupling special case

This model can also solve for pressure by enabling:

```text
solvePressure true;
```

This is a special case and is only supported when `solutionAlgorithm` is
`PETScSNES`. When pressure solving is enabled:

- the PETSc solution vector includes `p` as an extra unknown;
- the model uses the pressure stabilisation sub-model;
- the code enforces the `PETScSNES` path at construction time.

If `solvePressure` is `false`, the `stabilisation/pressure` sub-dictionary is
ignored.

In other words, `solvePressure true` is not a general option for the
segregated or explicit paths.

### Explicit-solve time-step control

When `solutionAlgorithm explicit` is selected, the solver updates `deltaT`
automatically from the estimated elastic wave speed and the mesh spacing. The
time-step scaling uses `maxCo` from `controlDict`, with a default of `0.1` if
it is not set.

### Recommended dictionary setup

Minimal example for `constant/solidProperties`:

```text
solidModel        linearGeometryTotalDisplacement;

linearGeometryTotalDisplacementCoeffs
{
    solutionAlgorithm   implicitSegregated;

    predictor           false;
    dampingCoeff        [0 0 -1 0 0 0 0] 0;
    relaxationMethod    fixed;

    nCorrectors         10000;
    solutionTolerance   1e-06;
    alternativeTolerance 1e-07;
    materialTolerance   1e-05;
    infoFrequency       100;

    restart             false;
    residualFile        false;
    writeResidualField  false;

    solvePressure       false;

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

    // Optional:
    // cellDisplacements
    // {
    //     ...
    // }
}
```

If `solvePressure` is enabled, add:

```text
    solvePressure        true;
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
`grad(D)`, `grad(DD)`, `sigma`, and related restart fields as needed.

### Field glossary

The main fields you will typically see in the output are:

- `D`: total displacement field at cell centres.
- `DD`: displacement increment, defined as `DD = D - D.oldTime()`.
- `pointD`: total displacement interpolated to the mesh points.
- `pointDD`: point displacement increment, defined as
  `pointDD = pointD - pointD.oldTime()`.
- `grad(D)`: cell-centre gradient of `D`.
- `grad(DD)`: cell-centre gradient of `DD`.
- `U`: velocity field, computed from the time derivative of `D`.
- `sigma`: symmetric stress tensor field.
- `A`: acceleration field used by the explicit algorithm.

Here `oldTime()` means the value stored at the previous time step.

The following boundary-condition types are treated specially by the solver:

- `solidTraction` patch fields;
- `fixedDisplacementZeroShear` patch fields;
- `symmetry` patch fields.

On traction boundaries, the solver enforces the traction and pressure
contributions consistently. On zero-shear or symmetry boundaries, the shear
traction is projected away from the normal direction.

### Output fields

Typical written fields include:

- `D`, `DD`
- `pointD`, `pointDD`
- `grad(D)`, `grad(DD)`
- `U`
- `sigma`
- linear-geometry strain fields and von Mises stress via `solidModel::writeFields()`

When the explicit algorithm is used, `A` is also maintained internally for the
central-difference update and may be present in the time directory depending on
write settings.

---

## Developer Notes

### Class role

`linGeomTotalDispSolid` is the linear-geometry analogue of the total-displacement
family. The key design choices are:

- `D` is the primary solution variable;
- the geometry is fixed in the governing equations;
- stress is delegated to `mechanicalModel`;
- the solver can run in three distinct modes: segregated implicit, PETSc SNES,
  and explicit.

The class inherits from `solidModel` and `foamPetscSnesHelper`.

### Main control flow

#### Construction

The constructor performs the following setup:

1. Calls the base `solidModel` constructor.
2. Creates the PETSc helper if PETSc is available and SNES is selected.
3. Builds the implicit stiffness fields `impK_`, `impKf_`, `rImpK_`, and the
   bulk-modulus reciprocal `rKappa_`.
4. Reads optional settings such as `predictor`, `solvePressure`,
   `dampingCoeff`, and the PETSc options file name.
5. Forces creation of old-time fields for consistent restart behavior.
6. Creates the default stabilisation dictionaries if they are missing.
7. Validates key consistency constraints:
   - `solvePressure` requires `PETScSNES`;
   - `predictor` cannot be used with a steady-state ddt scheme;
   - `grad(D)` should typically use `leastSquaresS4f` for the implicit paths.

The constructor also enables `extrapolateValue` on `solidTraction` patch fields
for the implicit and PETSc paths.

#### `evolve()`

`evolve()` dispatches to one of:

- `evolveImplicitSegregated()`
- `evolveSnes()`
- `evolveExplicit()`

Any other algorithm value is treated as unsupported for this class.

#### Segregated implicit path

`evolveImplicitSegregated()`:

- corrects `D` boundary conditions;
- optionally applies the linear predictor on the first call in a new time step;
- assembles a linear momentum equation for `D`;
- adds stabilisation and optional damping;
- enforces traction boundary conditions with `enforceTractionBoundaries()`;
- applies any cell-displacement constraints through `solidModel::setCellDisps()`;
- solves, relaxes, and checks convergence;
- updates `DD`, `gradD`, `gradDD`, `pointD`, `pointDD`, and `U`;
- repeats if the mesh update requests another pass.

This path uses a Laplacian-based approximation with the stabilisation term
stored in `impK_` and `impKf_`.

#### PETSc SNES path

`evolveSnes()`:

- optionally predicts `D` and inserts the predicted values into the SNES
  solution vector;
- solves the nonlinear system through `foamPetscSnesHelper`;
- extracts the solution back into `D`;
- optionally extracts `p` when `solvePressure` is enabled;
- recomputes `gradD`, `pointD`, `DD`, `pointDD`, and `U`.

The PETSc path uses a compact approximate Jacobian consistent with the
segregated discretisation, and `precondition()` is available for a physics-based
preconditioner.

The pressure-coupled branch adds one scalar unknown per cell, so the block size
is `2/3` for displacement only and `3/4` when pressure is included, depending on
the dimensionality.

#### Explicit path

`evolveExplicit()`:

- advances `U` and `D` using a central-difference style update;
- zeros out the out-of-plane components in 2-D cases;
- computes `gradD` and `sigma`;
- builds the traction field and applies traction boundary handling;
- calculates the acceleration field `A_` for the next explicit step.

`setDeltaT()` is responsible for stability control in this mode.

### Traction boundary enforcement

`enforceTractionBoundaries()` is a small but important helper. It recognises the
standard solids4foam displacement/traction patch types and applies model-specific
rules:

- `solidTraction`: uses the patch traction and subtracts the pressure term;
- `fixedDisplacementZeroShear` and `symmetry`: remove the tangential traction
  component by projection onto the normal direction.

This means the solver can work with traction BCs without each algorithm having
to reimplement the patch logic.

### Extension points

If you extend this model, the main places to consider are:

- `evolve()` for new solution algorithms;
- `enforceTractionBoundaries()` for new traction-style patches;
- the constructor if new model-specific checks or options are required.

Keep in mind that this class is intended to remain compatible with the base
`solidModel` workflow and with the standard solids4foam boundary conditions.
