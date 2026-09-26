# solidRobin: spring-dashpot (Robin) boundary condition

## Goal

Some solid boundaries rest on an elastic support: they are neither clamped nor
traction-free. `solidRobin` applies a traction proportional to the
displacement (spring) and to the velocity (dashpot), with separate normal and
tangential coefficients. Before this, solids4foam had no spring, Robin or
elastic-foundation condition for the displacement.

## Motivation: cardiac mechanics

Heart meshes are cut off at the valves, and in the body the ventricles are
held by the pericardium and the surrounding tissue. The usual way to model
this is with Robin conditions:

- **Epicardium (pericardium):** normal-only springs, often with stiffness that
  varies from apex to base, which lets the epicardium slide tangentially
  (Pfaller 2019; Strocchi 2020).
- **Base / valve rims / cut vessels:** springs in all directions instead of a
  fixed base. A fixed base removes base descent, the dominant systolic motion
  (Augustin 2021; Barnafi 2022).
- **Damping:** dynamic solid solvers need dashpots so the springs do not ring
  (Pfaller 2019).

## Formulation

$$
\mathbf{P}\mathbf{N} = -\mathbf{K}\,\mathbf{D} - \mathbf{C}\,\dot{\mathbf{D}}
\;(+\; \mathbf{t},\ p \text{ as in } \texttt{solidTraction})
$$

$$
\mathbf{K} = k_N\,\mathbf{N}\otimes\mathbf{N} + k_T\,(\mathbf{I}-\mathbf{N}\otimes\mathbf{N}),\qquad
\mathbf{C} = c_N\,\mathbf{N}\otimes\mathbf{N} + c_T\,(\mathbf{I}-\mathbf{N}\otimes\mathbf{N})
$$

- `P` is the first Piola-Kirchhoff stress and `N` is the unit normal in the
  reference configuration. `K` [Pa/m] and `C` [Pa s/m] act per unit reference
  area, as in Regazzoni 2022 and Barnafi 2022. Pfaller 2019 uses the
  normal-only case, `t = N (k u·N + c u̇·N)`.
- The optional `traction` and `pressure` act per unit current area, exactly
  as in `solidTraction`.
- `dD/dt` uses backward Euler: `(D - D.oldTime())/deltaT`.
- Limits: `K → 0` gives `solidTraction`. `kNormal → ∞` with
  `kTangential = 0` gives `fixedDisplacementZeroShear` with zero displacement.

## Implementation

`solidRobin` derives from `solidDirectionMixed`. The spring and dashpot go
into the matrix implicitly, and the traction goes through the solid model's
`tractionBoundarySnGrad`:

- `valueFraction = K_eff (K_eff + impK·deltaCoeffs)^-1`, per normal and
  tangential direction, with `K_eff = (K + C/Δt)·dA/da`
- `refValue = 0`
- `refGrad = tractionBoundarySnGrad(t + (C/Δt)·D_old, p)`

At convergence this reproduces the traction above exactly. `impK` is taken as
the slope of `tractionBoundarySnGrad`, so it matches the solid model's own
implicit stiffness, including `2μ` when `solvePressure` is used.

Why not derive from `solidTraction`? With `t = -K·D` applied explicitly, the
spring force lags one outer iteration. For the global stretch/rigid mode of
the body that lag is amplified by roughly `k L / E`, which is well above 1
for stiff cardiac springs, so the explicit version would need heavy
relaxation. The implicit version converges for any `k` (verified up to
`k L / E ≈ 7400`).

Supported solid models: linear geometry and total Lagrangian
(e.g. `linearGeometryTotalDisplacement`,
`nonLinearGeometryTotalLagrangianTotalDisplacement`). In the total Lagrangian
case, the reference-to-current area ratio `dA/da` is taken from `J` and
`Finv`. Updated Lagrangian models are refused, because they move the mesh and
the patch normals are then not reference normals.

## Usage

```
EPI
{
    type            solidRobin;
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

## Typical values from the literature

| Source | Boundary | Form | Stiffness | Damping |
|---|---|---|---|---|
| Pfaller 2019 | epicardium | normal | k = 0.2 kPa/mm = 2e5 Pa/m (swept 0.1–5 kPa/mm) | c = 5e-3 kPa s/mm = 5e3 Pa s/m |
| Pfaller 2019 | great vessels | all directions | 2e3 kPa/mm | 1e-2 kPa s/mm |
| Strocchi 2020 | epicardium | normal, varies in space | up to 50 kPa/mm at the apex (swept 1–50) | — |
| Strocchi 2020 | cut veins | all directions | 10 kPa/mm | — |
| Barnafi 2022 | base | normal / tangential | K⊥ = 2e5, K∥ = 2e4 Pa/m | C⊥ = 1e4, C∥ = 2e3 Pa s/m |

Unit conversion: 1 kPa/mm = 1e6 Pa/m.

## Verification

Tested with OpenFOAM v2412, against analytical solutions:

| Test | Result |
|---|---|
| Bar in uniaxial strain on a normal spring, loaded at the far end, k = 1e4 … 1e8 Pa/m | spring-face and loaded-face displacements exact to ≤ 2e-8 |
| Stiff spring, k L / E ≈ 7400 | exact; 120 outer iterations vs 106 for `fixedDisplacement` |
| Load applied through the BC's own `traction` entry, far end fixed | 2e-10; 113–139 iterations for all k |
| Simple shear with cyclic sides (tangential spring) | 1e-9; normal leakage 4e-21 |
| Mass-spring-dashpot (dynamic, `d2dt2 Euler`) | matches the discrete oscillator to 1.6e-4 (≈ k L / M); period 44.40 ms vs 2π√(m/k) = 44.43 ms; damped (ζ = 0.2) period 45.50 ms vs 45.35 ms |
| Total Lagrangian, 25 % strain, spring face area changes ≈ 13 % | force balance F = k·D·A₀ exact to 7e-8 (14 % error when the area scaling is disabled) |
| `weightField`, restart, 2-processor run | identical to 1e-11 |

Note: a static body held only by soft springs and loaded with `solidTraction`
converges slowly in the segregated solver (about 1000 outer iterations at
k L / E ≈ 0.015). The cause is the nearly free rigid-body mode, not this
condition: with the rigid mode removed, the iteration count does not depend on
k. In dynamic runs, inertia regularizes that mode.

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
