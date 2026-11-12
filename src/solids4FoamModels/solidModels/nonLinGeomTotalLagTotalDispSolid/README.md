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
- uses the selected mechanical law to compute stress;
- adds `fvOptions` sources, for example a manufactured body force, to the
  momentum equation of both solution paths on OpenFOAM.com. The source is per
  unit reference volume. On OpenFOAM.org and foam-extend a derived model can
  supply the same source through the protected `fvOptionsSource()` function,
  as the hyperelasticity `manufacturedSolution` tutorial does.

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
| `infoFrequency` | `100` | Frequency for solver progress output |
| `restart` | `false` | Writes extra fields needed for a consistent restart |
| `writeResidualField` | `false` | Writes a residual field during output |
| `residualFile` | `false` | Writes `residual.dat` in the case directory |
| `stabilisation` | auto-created if absent | Has `momentum` and `pressure` |
| `cellDisplacements` | optional | Internal-cell displacement constraints |
| `solvePressure` | `false` | Pressure unknown, only with `PETScSNES` |
| `mixedImplicitStiffness` | `deviatoric` | `deviatoric` or `twoMu`, only with `solvePressure` |
| `scaleMixedPETScFields` | `false` | Store the pressure unknown as `p/s` in PETSc, only with `solvePressure` |
| `pressureUnknownScale` | `twoMu` | `s`: `twoMu`, `user` (with `pressureUnknownScaleValue`) or `none` |
| `pressureScaleByTwoMu` | as `scaleMixedPETScFields` | Multiply the pressure rows by 2μ |
| `pressureScaleFactor` | `1` | Further factor on the pressure rows |

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

The mechanical law is asked for its stress without its volumetric response,
which the solved pressure replaces (`sigma = sigma_law - p I`, with no
deviatoric projection of the total stress). The stress is therefore whatever
the law defines apart from that response. For the cardiac laws:

- `GuccioneElastic` evaluates its energy on the isochoric deformation
  `Fbar = J^(-1/3) F`, so its stress is the isochoric (deviatoric) Cauchy stress
  `dev(Fbar S Fbar^T)/J`, and its volumetric response `0.5 K (J^2 - 1)/J` is
  the one the pressure equation uses;
- `electroMechanicalLaw` adds the active stress `F (Ta f0 f0) F^T / J` in
  full, with the physical `F`, including its spherical part, which is not part
  of any energy and so is not replaced by the pressure.

The pressure variable `p` is therefore not the mean stress where an active
stress acts: the mean Cauchy stress is `-p + tr(sigma_active)/3`.

Parallel runs are supported across `processor` patches. The assembled mixed
Jacobian blocks recognise only that patch type: a mesh with any other coupled
patch (`cyclic`, `cyclicAMI`, `cyclicACMI`, or the `processorCyclic` patches
that decomposing a cyclic mesh creates) stops in the first Jacobian assembly,
in serial or parallel.

### Mixed implicit stiffness

With `solvePressure true`, the implicit stiffness `impK` of the mixed
formulation is selected by

```text
mixedImplicitStiffness deviatoric;   // default
mixedImplicitStiffness twoMu;
```

- `deviatoric` (default, unchanged behaviour): the law's deviatoric scalar
  tangent, the Laplacian surrogate for div(dev(sigma)). For `linearElastic`,
  `neoHookeanElastic` and `GuccioneElastic` it is (4/3) mu; other laws define
  their own.
- `twoMu`: 2 mu, as declared by the mechanical constitutive law
  (`tangentRequest::scalarShear`). Supported by `linearElastic`,
  `neoHookeanElastic`, `GuccioneElastic` and an `electroMechanicalLaw` whose
  passive law is one of these; any other law is refused with an error. For
  `GuccioneElastic`, mu = k (cf + ct + cfs)/6 is a declared scalar stiffness
  measure, not the shear modulus of an isotropic material.

`impK` weights the momentum stabilisation and sets the pressure-stabilisation
diffusivity (`rAUf`, the reciprocal of the momentum diagonal), as well as the
preconditioner. Changing it therefore changes both stabilisation strengths
at fixed `scaleFactor`s: for a steady case with uniform properties, `twoMu`
is equivalent to `deviatoric` with the momentum scale factor multiplied by
3/2 and the pressure scale factors by 2/3. The keyword is refused without
`solvePressure`.

### Mixed system scaling

The displacement and pressure unknowns, and the momentum and pressure
equations, can differ by many orders of magnitude, which degrades the
matrix-free (MFFD) differencing step and the linear solves. The mixed PETSc
system can be scaled (opt-in; off by default, so existing cases are
unchanged):

```text
scaleMixedPETScFields yes;     // pressure unknown pHat = p/s_u
pressureUnknownScale  twoMu;   // s_u: twoMu (default), user or none
// pressureUnknownScaleValue 1e6;   // s_u with pressureUnknownScale user
// pressureScaleByTwoMu yes;        // default: as scaleMixedPETScFields
// pressureScaleFactor  1;          // s_r = pressureScaleFactor*(2 mu or 1)
```

With the PETSc unknowns (D, pHat), p = s_u pHat, and the pressure rows of the
residual multiplied by s_r, the system has the same roots, so the converged D
and p are unchanged; only the conditioning, and the norm used by the
nonlinear stopping test, change. The assembled matrix, which is an
approximate preconditioning matrix (the Newton operator is matrix-free with
`snes_mf_operator`), is transformed consistently: its A_Dp block is multiplied
by s_u, A_pD by s_r and A_pp by s_r s_u, in every cell, boundary and
processor coefficient. The transformation is exact; it does not make the
assembled blocks an exact Jacobian. The OpenFOAM `p` field is always the
physical pressure.

The scaling balances magnitudes; it is not a nondimensionalisation. The
extensive pressure rows have units of m³, so multiplied by 2μ they are in
N·m, a length factor (of the order of the cell size) away from the momentum
rows in N. Because the residual norm changes, SNES tolerances change meaning
when the scaling is switched on: the pressure rows are multiplied by 2μ, so
`snes_atol`, and `snes_rtol` when the initial pressure residual dominates,
may need retuning.

`twoMu` is the volume average of 2μ as the laws declare it (the `scalarShear`
tangent, as for `mixedImplicitStiffness twoMu`), independently of
`mixedImplicitStiffness`. A law that declares no 2μ is refused with a pointer
to an explicit scale: `pressureUnknownScale user` with
`pressureUnknownScaleValue`, and `pressureScaleByTwoMu no` with
`pressureScaleFactor`. The keywords and defaults are those of the original
implementation by Philip Cardiff on the `fix-pressure-solid` branch and its
port to the total-Lagrangian solver in #627 (Chanikya Valeti), with two
changes: the switch is spelled `scaleMixedPETScFields` (the earlier
`scaleMixedPetScFields` is still read, with a warning), and the undocumented
`pressureUnknownScale 2mu` alias is no longer accepted. There, `twoMu` was
1.5 `impK`, which equals 2μ only where `impK` is (4/3)μ.

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
- `solidSpringDashpot`
- `fixedDisplacementZeroShear`
- `symmetry`
- `slip`

On traction patches, the force is formed using the current deformed normals and
area measures. On spring-dashpot patches, the spring and dashpot force uses the
reference area and the traction and pressure use the current area. On
zero-shear, symmetry, and slip patches, the tangential force is projected away.

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
- stress is delegated to the `mechanicalConstitutiveLaw` framework;
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
