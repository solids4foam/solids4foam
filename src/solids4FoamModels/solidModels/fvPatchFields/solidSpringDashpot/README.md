# solidSpringDashpot: spring-dashpot (Robin) boundary condition

## Goal

Some solid boundaries rest on an elastic support: they are neither clamped nor
traction-free. `solidSpringDashpot` applies a traction proportional to the
displacement (spring) and to the velocity (dashpot), with separate normal and
tangential coefficients. Mathematically it is a Robin condition on the
displacement.

## Applications

The condition replaces surrounding material that is not meshed by a
distributed spring-dashpot support:

- **Heart:** the pericardium and the surrounding tissue hold the ventricles.
  The epicardium is usually given normal-only springs, which let it slide
  tangentially, often with stiffness varying from apex to base (Pfaller 2019;
  Strocchi 2020). The cut base, valve rims and vessels get springs in all
  directions instead of a fixed base, which would suppress base descent
  (Augustin 2021; Barnafi 2022).
- **Arteries:** the perivascular tissue supports the outer wall (Moireau
  2012). For cerebral arteries and aneurysms this is the surrounding brain
  tissue and cerebrospinal fluid, with a much stiffer support where the vessel
  touches bone (Shidhore 2023).
- **Elastic foundations:** any body resting on a compliant support.

## Formulation

$$
\mathbf{P}\mathbf{N} = -\mathbf{K}\,\mathbf{D} - \mathbf{C}\,\dot{\mathbf{D}}
+ \mathbf{t} - p\,\mathbf{n}
$$

$$
\mathbf{K} = k_N\,\mathbf{N}\otimes\mathbf{N} + k_T\,(\mathbf{I}-\mathbf{N}\otimes\mathbf{N}),\qquad
\mathbf{C} = c_N\,\mathbf{N}\otimes\mathbf{N} + c_T\,(\mathbf{I}-\mathbf{N}\otimes\mathbf{N})
$$

- `P` is the first Piola-Kirchhoff stress and `N` is the unit normal in the
  reference configuration. `K` [Pa/m] and `C` [Pa s/m] act per unit reference
  area, as in Regazzoni 2022 and Barnafi 2022. Pfaller 2019 uses the
  normal-only case, `t = N (k u·N + c u̇·N)`.
- The optional `traction` `t` and `pressure` `p` act per unit current area,
  with `n` the current normal. They are constant in time (no time series).
- `dD/dt` uses backward Euler, `(D - D.oldTime())/deltaT`, whatever the
  `ddtSchemes` and `d2dt2Schemes` entries are.
- Limits: `K → 0` gives a traction condition. `kNormal → ∞` with
  `kTangential = 0` gives `fixedDisplacementZeroShear` with zero displacement.

## Implementation

`solidSpringDashpot` derives from `solidDirectionMixed`:

- `valueFraction = K_eff (K_eff + impK·deltaCoeffs)^-1`, per normal and
  tangential direction, with `K_eff = (K + C/Δt)·dA/da`
- `refValue = 0`
- `refGrad = tractionBoundarySnGrad(t + (C/Δt)·D_old, p)`

This is a Robin condition with coefficients `K_eff` and `impK`, so the spring
and dashpot are implicit. `impK` is taken as the slope of
`tractionBoundarySnGrad`, so it matches the solid model's own implicit
stiffness, including `2μ` when `solvePressure` is used.

- **Matrix diagonal:** `directionMixedFvPatchField::snGradTransformDiag()`
  returns `sqrt(valueFraction_ii)`. For a soft spring (small valueFraction)
  this makes the matrix much stiffer than the spring, and the outer
  iterations converge slowly. This condition returns `valueFraction_ii`, the
  diagonal of the linearised boundary gradient. Only the coupling between
  components stays explicit when `kNormal ≠ kTangential`.
- **Imposed traction:** `linearGeometryTotalDisplacement` and
  `nonLinearGeometryTotalLagrangianTotalDisplacement` impose
  `-K·D - C·(D - D_old)/Δt + t - p n` directly on the boundary faces,
  through `springDashpotTraction()`, as they do for `solidTraction`. This is
  used by the implicit segregated, explicit and PETSc SNES algorithms.
- **Area ratio:** in the total Lagrangian model, `dA/da` is taken from `J` and
  `Finv`, because `tractionBoundarySnGrad` expects a traction per unit current
  area there.

Supported solid models: total-displacement models with linear geometry or a
total Lagrangian formulation (verified with `linearGeometryTotalDisplacement`
and `nonLinearGeometryTotalLagrangianTotalDisplacement`). Incremental models,
which solve for `DD`, and updated Lagrangian models, whose patch normals are
not reference normals, are refused with a fatal error.

The `nonOrthogonalCorrections` and `limitCoeff` entries of `solidDirectionMixed`
are read. `secondOrder yes` is refused, because the valueFraction above assumes
the first-order boundary value; with `secondOrder` the factor would be
`2·deltaCoeffs`.

## Usage

```text
outerWall
{
    type            solidSpringDashpot;
    kNormal         uniform 2e5;        // [Pa/m], required
    kTangential     uniform 0;          // [Pa/m], default 0
    cNormal         uniform 5e3;        // [Pa s/m], default 0
    cTangential     uniform 0;          // [Pa s/m], default 0
    traction        uniform (0 0 0);    // [Pa], default 0
    pressure        uniform 0;          // [Pa], default 0
    weightField     apicobasal;         // optional
    value           uniform (0 0 0);
}
```

- All coefficients can be `nonuniform` lists, one value per face.
- `weightField` names a `volScalarField`, read from the start time. Its patch
  values multiply all four coefficients once, at construction. The weighted
  coefficients are written out instead of the field name, so restarts and
  `decomposePar` do not apply the weight twice.
- **Static runs:** the dashpot uses `deltaT` as a physical time step. When
  `deltaT` is a pseudo time for load stepping, set `cNormal` and
  `cTangential` to 0; a warning is printed when `d2dt2Schemes` is
  `steadyState` and a dashpot coefficient is non-zero. Quasi-static runs in
  physical time, as in cardiac electromechanics, can keep the dashpot.

## Typical values from the literature

Unit conversions: 1 kPa/mm = 1e6 Pa/m; 1 dyne/cm³ = 10 Pa/m.

- **Pfaller 2019, epicardium:** normal springs and dashpots,
  k = 0.2 kPa/mm = 2e5 Pa/m (swept 0.1–5 kPa/mm),
  c = 5e-3 kPa s/mm = 5e3 Pa s/m.
- **Pfaller 2019, great vessels:** all directions, k = 2e3 kPa/mm,
  c = 1e-2 kPa s/mm.
- **Strocchi 2020, epicardium:** normal springs varying in space, up to
  50 kPa/mm at the apex (swept 1–50 kPa/mm), no damping.
- **Strocchi 2020, cut veins:** all directions, k = 10 kPa/mm.
- **Barnafi 2022, base:** K⊥ = 2e5, K∥ = 2e4 Pa/m; C⊥ = 1e4,
  C∥ = 2e3 Pa s/m.
- **Shidhore 2023, cerebral aneurysms:** all directions, no damping;
  k = 1e5, 1e6 and 1e7 Pa/m swept and 1e7 Pa/m used for tissue contact,
  1e10 Pa/m for bone contact (wall E = 1 MPa, ν = 0.49).
- **Moireau 2012, aorta:** the spring and dashpot coefficients are calibrated
  against measured wall motion.

## Verification

Tested with OpenFOAM v2412 against exact solutions (bar of length L = 10 mm,
E = 100 kPa, ν = 0.3, 20 cells along the bar):

- **Bar in uniaxial strain on a normal spring**, loaded at the far end,
  k = 1e2 … 1e8 Pa/m: spring-face and loaded-face displacements exact to
  ≤ 6e-9.
- **Load through the condition's own `traction` entry**, far end fixed:
  2e-10.
- **Simple shear with cyclic sides** (tangential spring): 1e-9; normal
  leakage 4e-21.
- **Mass-spring-dashpot** (dynamic, `d2dt2 Euler`, 20 % of critical damping):
  matches the discrete backward-Euler oscillator to 4e-5 (of order k L / M).
  The explicit algorithm (Δt = 1e-6 s) matches the implicit run
  (Δt = 1e-5 s) to 8e-4.
- **Total Lagrangian, free lateral faces, 20 % strain**, spring-face area
  change ≈ 13 %: force balance F = k·D·A₀ exact to 2e-7.
- **`weightField`, restart and a 2-processor run:** exact to 2e-10.
- **Bar with a body force** (quadratic solution): the loaded-end error falls
  by 4 with each mesh halving (2.3e-6, 5.7e-7, 1.4e-7, 3.5e-8 for 10–80
  cells), so the first-order boundary value gives second-order convergence.
- **Thick cylinder under internal pressure, normal springs on the curved
  outer surface** (plane strain, exact Lamé solution with the Robin outer
  condition): inner radial displacement errors 1.7e-3, 4.0e-4 and 9.8e-5 on
  10×20, 20×40 and 40×80 cells, with 51–64 outer iterations on all meshes
  and k = 1e5 … 1e8 Pa/m.
- **Refusals:** `secondOrder yes` and an updated Lagrangian (incremental)
  model stop with a fatal error; a dashpot in a `steadyState` run prints the
  warning.

Outer iterations to a step-norm tolerance of 1e-10, against the same
condition with the base-class diagonal and without the imposed traction:

| Case | Base class | This condition |
|---|---|---|
| Bar, k = 1e2 Pa/m | not converged in 10000 | 57 |
| Bar, k = 1e4 Pa/m | 4072 | 40 |
| Bar, k = 2e5 Pa/m | 982 | 64 |
| Bar, k = 1e8 Pa/m | 120 | 104 |
| Dynamic, 1000 steps (total) | 53155 | 5423 |
| Total Lagrangian, free faces | 2428 | 322 |
| Total Lagrangian, uniaxial strain | 452 | 86 |

The tutorial `tutorials/solids/linearElasticity/springSupportedBar` contains
the bar on a normal spring and the total Lagrangian force balance, with a
regression test against the exact solutions.

## References

- Pfaller M.R. et al. (2019). The importance of the pericardium for cardiac
  biomechanics: from physiology to computational modeling. Biomech Model
  Mechanobiol 18:503–529. doi:10.1007/s10237-018-1098-4
- Strocchi M. et al. (2020). Simulating ventricular systolic motion in a
  four-chamber heart model with spatially varying Robin boundary conditions to
  model the effect of the pericardium. J Biomech 101:109645.
  doi:10.1016/j.jbiomech.2020.109645
- Augustin C.M. et al. (2021). A computationally efficient physiologically
  comprehensive 3D–0D closed-loop model of the heart and circulation. Comput
  Methods Appl Mech Eng 386:114092. doi:10.1016/j.cma.2021.114092
- Regazzoni F. et al. (2022). A cardiac electromechanics model coupled with a
  lumped-parameter model for closed-loop blood circulation. Part I: model
  derivation. J Comput Phys 457:111083. doi:10.1016/j.jcp.2022.111083
- Barnafi N.A. et al. (2022). A comparative study of scalable multilevel
  preconditioners for cardiac mechanics. arXiv:2208.06191
- Moireau P. et al. (2012). External tissue support and fluid–structure
  simulation in blood flows. Biomech Model Mechanobiol 11:1–18.
  doi:10.1007/s10237-011-0289-z
- Shidhore T.C. et al. (2023). Comparative assessment of biomechanical
  parameters in subjects with multiple cerebral aneurysms using
  fluid–structure interaction simulations. J Biomech Eng 145:051003.
  doi:10.1115/1.4056317
