# womersleyTube: investigation of the sub-nominal temporal order

Status: **complete.** Sections 1-4 are the code trace and the hypotheses as
written before any run (commit 8f8c8468e), unchanged. Sections 5-15 are the
results. Short answer: the sub-nominal order is caused by an O(dt)
mass-flux inconsistency on the artificial tube ends, where the exact
pressure and the exact normal velocity gradient are imposed (section 8); it
is not an FSI, interface, ALE, solid, start-up or tolerance effect. A
small, opt-in `pimpleFluid` option that makes the end flux consistent with
the boundary velocity flux (`fluxConsistentPatches`, section 6.10) removes
the O(dt) term and gives second order in every QoI of the production study
(m2, Robin-Neumann, exact end data unchanged): 1.87-2.08 from 100/200/400
steps per period. It leaves an O(dt h) boundary-flux term, so at fixed mesh
the order formally tends to one as dt -> 0, well beyond the steps used
(section 15). It also reduces most finest-mesh errors of the mesh study
(the wave speed by about 70 times; the profile is unchanged and the flow
phase, 1.5e-5 before, is 8.3e-5).

Scope: the `womersleyTube` verification case only. No change to the
production numerical methods; the manuscript is not edited here.

## 1. The anomaly as reported

PR #513 (OpenFOAM-v2412, macOS) reported, for the time-step study on mesh
factor 2 (64 axial cells) with Robin-Neumann coupling, 50/100/200 steps per
period, six periods, Fourier coefficient over periods 5-6:

| Quantity | Time order (successive differences) |
|---|---:|
| flow_amp | 1.64 |
| flow_phase | 1.84 |
| wallMid_amp | 1.71 |
| wallMid_phase | 1.74 |
| speed | 1.29 |
| attenuation | 1.51 |

plus a start-up transient of about 0.5% of the wall displacement in the
first period that *grows* as the time-step is reduced, and a 0.26%-of-
amplitude v2412/v2512 shift in the Robin-Neumann regression value.

## 2. What the code actually does (traced before forming hypotheses)

These points correct or sharpen premises of the task description.

1. **The case does not start from rest.** The fluid velocity and pressure,
   the solid displacement `D` and three old-time levels `D_0`, `D_0_0`,
   `D_0_0_0` (at -dt, -2dt, -3dt), and the solid point displacement `pointD`
   are all set to the exact solution with `#codeStream`. There is no ramp.
   "Exact initialisation" therefore already exists; the question is whether
   it is *discretely consistent*.
2. **But the initial state is not fully consistent.**
   - (a) The old-time solid fields set only their *internal* (cell) values.
     Their boundary values on `interface`, `outer` and `inlet|outlet` are the
     placeholder `value uniform (0 0 0)` and old-time fields are never
     re-evaluated. `solidModel::faceZoneAcceleration()` evaluates the
     `backward` `d2dt2` on the **boundary** values of `D`, `D_0`, `D_0_0`,
     `D_0_0_0` (and a lazily created fifth level), and this is the
     acceleration that the Robin condition (`elasticWallPressure`) imposes on
     the fluid, `dp/dn = -rho a_n` at convergence, and that enters the
     kinematically-consistent interface flux in `pimpleFluid`. For the first
     three to four steps the imposed acceleration is therefore of order
     `|D|/dt^2` instead of `omega^2 |D|`: about 250 times the true value at 100
     steps per period, scaling as `dt^-2`.
   - (b) The initial `pointD` is the exact displacement at the points, while
     from the first step on `pointD` is the solid model's interpolation of
     the cell values. The fluid interface moves by the increment
     `pointD - pointD.oldTime()`, so the first step carries the interpolation
     mismatch as a spurious wall motion, a spurious wall velocity of order
     `mismatch/dt`.
   - (c) The fluid has no `U_0`, `Uf_0` or old mesh: its first step is Euler
     in both OpenFOAM versions (a single first-order step, normally harmless).
3. **Interface velocity** (`newMovingWallVelocity` for IQN-ILS,
   `elasticWallVelocity` for Robin): tangential part = BDF2 of the face
   centres; normal part = `fvc::meshPhi(U)`, which for `backward` in
   OpenFOAM.com is `1.5 phi^{n+1/2} - 0.5 phi^{n-1/2}`, a second-order
   extrapolation of the swept-volume flux to `t^{n+1}`. Nominally second order
   after the first step.
4. **Robin pressure**: `p + c1 dp/dn = p_prev - c1 rho a_prev`; at coupling
   convergence `dp/dn = -rho a_n(solid boundary)`. The interface flux is
   replaced by `meshPhi - rAU a_n |Sf|` so the converged flux equals the mesh
   flux. `a_n` is the solid `backward` `d2dt2` of the boundary displacement:
   nominally second order once the history is genuine.
5. **Mesh motion**: `velocityLaplacian`, interface `pointMotionU = Delta d/dt`,
   so mesh positions follow the solid exactly (to the coupling tolerance);
   the fluid mesh sits at a constant offset `-d(0)` from the solid
   (5e-4 R, time-independent).
6. **OpenFOAM version difference already visible in the sources**: the
   installed v2412 has the solids4foam `optionalFixes/OpenFOAM-v2412/
   backwardDdtScheme.C` compiled in (source and `libfiniteVolume.so` dated
   2025-08-05). Its `deltaT0_(vf)` falls back to Euler when the field's old
   and old-old levels have the same time index. Stock v2512 falls back to
   Euler when `time().timeIndex() < 2`, for every field. The two differ only
   in which `fvc::ddt`/`meshPhi` calls are Euler during the first steps
   (and for lazily created old levels), i.e. exactly in the start-up.

## 3. Hypotheses (written before the study)

"Start-up" means anything that only acts in the first few steps; such an
error can still contaminate the analysis window if it excites a slowly
decaying free mode (the elastic wall is undamped apart from the fluid
viscosity and the weak BDF2 dissipation, 0.9996 per period at 100 steps).

| # | Hypothesis | Mechanism | QoIs affected | Supporting signature | Refuting signature | Cheapest test |
|---|---|---|---|---|---|---|
| H1 | BDF2 first-step (Euler) treatment | One first-order step: local O(dt^2) error once | All, through the transient | Transient amplitude falls as dt^2 | Transient grows as dt falls (already reported) | Per-period Fourier coefficients vs dt |
| H2a | Zero old-time boundary values of `D_0..D_0_0_0` | Spurious boundary acceleration O(D/dt^2) for 3-4 steps fed to the Robin pressure condition: an impulse whose velocity kick grows as 1/dt | All (Robin); not IQN-ILS | Transient grows as dt falls; Robin transient >> IQN-ILS transient; removed by exact old-time boundary values | Transient unchanged when boundary values are exact | Rerun with exact boundary values in the old-time files (an initial-condition fix only) |
| H2b | Initial `pointD` inconsistent with the solid interpolation | First-step mesh jump: spurious wall velocity mismatch/dt | All, both couplings | First-step wall velocity spike; transient in IQN-ILS too | No spike at step 1 | Inspect step-1 interface motion; IQN-ILS transient |
| H3 | Start-up ramp | None exists | – | – | – | (not applicable; H2 covers start-up) |
| H4 | Interface velocity first order in time | O(dt) wall velocity error | flow, profile, wall | Order -> 1 as dt -> 0 in a clean (transient-free) run | Order -> 2 in a clean run | Clean run; compare IQN-ILS and Robin orders |
| H5 | Robin pressure condition first order in time | `a_n` or `prevPressure` lagged a step | Pressure-related: speed, attenuation, wall | Robin order < IQN-ILS order at same dt | Robin and IQN-ILS converge together | Time study with IQN-ILS too |
| H6 | ALE / mesh-flux (GCL) time error | `backward` with extrapolated meshPhi is not GCL-exact | All, at relative size ~ d/R = 5e-4 | Error independent of the QoI scale, ~1e-5 | Effect below 1e-5 | Order-of-magnitude bound only |
| H7 | Predictor / interface extrapolation | Changes only the first iterate | Through tolerance (H10) | Orders change with the predictor | – | Folded into H10 |
| H8 | Extraction artefact | Fourier leakage of transient free modes; 10 fluid samples per period; log-pressure fit | Fluid-sample QoIs (speed, attenuation, profile) | Orders change with window / extraction | Same orders under alternative extraction | Re-analyse stored histories: other windows, 1-period vs 2-period, least-squares with transient term |
| H9 | Spatial-temporal interaction on mesh m2 | Non-separable terms (e.g. Rhie-Chow `ddtCorr`, `rAU ~ dt` in the pressure equation, wall-face extrapolation) give errors `C h^2 g(dt)` that do not cancel in successive differences | All; most where the spatial error is largest (wall phase, attenuation) | Orders improve from m1 -> m2 -> m4 | Orders identical on m1, m2, m4 | Time study on m1 (cheap) and m4 |
| H10 | Coupling / solver tolerance floor | Per-step residual errors accumulate over N = steps x periods steps; a floor grows relative to O(dt^2) as dt falls | All, most for small differences | Results move by a fraction of the successive differences when tolerances are tightened 100x | Movement << successive differences | One run at 200 steps with tolerances tightened 100x |
| H11 | Boundary data at the wrong time level | Coded BCs evaluated at t^n instead of t^{n+1} | All | Phase error O(dt) | Coded BCs use the current time | Code inspection (they use `time().value()` = t^{n+1}) |
| H12 | OpenFOAM detail making a term first order part of the time | e.g. `timeIndex() < 2` Euler fallback for every field (v2512), lazily created old levels, ddtCorr | Start-up; version differences | v2412/v2512 differences confined to the start-up | Differences persist in the periodic state | Version comparison of per-period coefficients |

Ranking before running (from the code trace only): H2a (Robin study) and
H2b first, then H9 and H10, then H8; H4/H5/H6 are not expected to be first
order from the code, and H11 is ruled out by inspection (the coded
conditions read `this->db().time().value()`, the new time).

## 4. Planned diagnostics (cheapest first)

1. Per-period analysis of the existing study's histories (and of a re-run on
   this machine): transient amplitude and decay per period, versus dt.
2. Window test: analyse periods 3-4, 5-6 and, with longer runs, later periods.
3. Initial-condition consistency test (H2a): exact old-time boundary values.
4. Mesh dependence (H9): time study on m1 (cheap), then m4.
5. Extend to 400 steps per period.
6. Tolerance test (H10) at one fine step.
7. IQN-ILS time study subset (H4/H5).
8. Alternative extraction on stored histories (H8).
9. Version difference (v2412 vs v2512), minimal.

---

# Results

All runs: solids4foam `development` at aee0c8e35 (no source change), built
privately; OpenFOAM-v2512 (stock) unless marked v2412 (which has the
solids4foam `optionalFixes` `backwardDdtScheme.C` compiled in); Linux
x86_64, PETSc 3.24; one core per run under Slurm. Driver:
`scripts/womersley_temporal_investigation.py` (reuses the set-up, readers
and Fourier extraction of `womersley_tube_verification.py`). Machine-readable
results are in `temporal/`:

- `temporal_runs.csv`: one row per run and two-period window (mesh factor,
  steps per period, dt, periods, variant, start-up, coupling, tolerances,
  OpenFOAM build, iterations, run time, every QoI error against the exact
  solution, including the alternative extractions `speed_alt`,
  `attenuation_alt` (complex-exponential fit of p) and `speed_ux`,
  `attenuation_ux` (wave number from the axial velocity));
- `temporal_orders.csv`: for each series and window, the values, successive
  differences and successive-difference orders, and the error-based orders;
- `runs/<run>.json`: per-period and per-window coefficients, the transient
  diagnostics and the residual spectrum of every run;
  `runs/robin_m2_n*_p12_history.csv.gz`: wall displacement and flow-rate
  histories of the main study series;
- `end_flux_inconsistency.csv`: the measured end-face flux inconsistency.

Run names are `<coupling>_m<mesh factor>_n<steps per period>_p<periods>[_variants]`.
Variants (all in the driver, none changes the production case):

| Variant | What changes |
|---|---|
| `exactOldBoundary` | boundary values of `D`, `D_0`, `D_0_0`, `D_0_0_0` set to the exact solution (initial data only) |
| `tight` | coupling tolerance 1e-8, Robin 1e-7/5e-5, SNES rtol/atol/stol 1e-10/1e-11/1e-8, fluid solvers 1e-14 |
| `pimple4`, `pimple10` | 4 or 10 PIMPLE outer correctors |
| `fluidOnly` | fluid region alone; the interface mesh velocity is the exact wall displacement increment / dt (the mesh then follows the exact wall from its undeformed start, as it follows the solid in the coupled case); `newMovingWallVelocity`, zero-gradient p |
| `transpiration` | (after `fluidOnly`) static mesh; the exact wall velocity imposed at the undeformed wall: no mesh motion, ALE flux or moving-wall condition |
| `ddtCoeff1`, `noDdtCorr` | constant Rhie-Chow `ddtCorr` coefficient (`backward 1`); `ddtCorr` off |
| `lsS4f` | fluid `gradSchemes leastSquaresS4f` |
| `rigid` | (self-convergence) no-slip wall |
| `endsZeroGrad` | (self-convergence) zero-gradient velocity at the tube ends |
| `endsFixedGradient` | the same exact end velocity gradient through `fixedGradient` (set each step by a coded function object) instead of `codedMixed` |
| `endsDirichletU` | (fluid-only) exact end velocity with `fixedFluxPressure` and a pressure reference cell |
| `of2412` | the OpenFOAM-v2412 build |

## 5. Reproduction of the original anomaly (A)

Robin-Neumann, m2, periods 5-6, v2512/Linux:

| Quantity | n50 | n100 | n200 | order (this work) | order (PR #513, v2412/macOS) |
|---|---:|---:|---:|---:|---:|
| flow_amp | -2.061e-3 | -1.701e-4 | 4.371e-4 | 1.64 | 1.64 |
| flow_phase | 1.243e-4 | 1.053e-3 | 1.311e-3 | 1.85 | 1.84 |
| wallMid_amp | 7.135e-3 | 8.427e-4 | -1.083e-3 | 1.71 | 1.71 |
| wallMid_phase | -7.872e-3 | -3.737e-3 | -2.494e-3 | 1.73 | 1.74 |
| speed | 2.969e-3 | 1.208e-3 | 4.889e-4 | 1.29 | 1.29 |
| attenuation | -1.988e-2 | -9.629e-3 | -6.040e-3 | 1.51 | 1.51 |

The stored values are reproduced to three or four digits and the orders to
0.01. The anomaly is reproducible and platform independent.

## 6. Diagnostics and results

### 6.1 Finer time steps (B) and analysis window (E)

Robin, m2, successive-difference orders from 50/100/200 and 100/200/400:

| Quantity | periods 5-6 | periods 11-12 | exact old-time boundary values (5-6) | IQN-ILS (5-6) |
|---|---|---|---|---|
| flow_amp | 1.64 / **1.12** | 1.64 / 1.12 | 1.64 / 1.12 | 1.94 / 0.41 |
| flow_phase | 1.85 / 1.79 | 1.85 / 1.83 | 1.85 / 1.83 | 1.53 / 4.54 |
| wallMid_amp | 1.71 / **1.31** | 1.71 / 1.32 | 1.71 / 1.32 | 1.78 / 0.66 |
| wallMid_phase | 1.73 / **1.14** | 1.73 / 1.14 | 1.73 / 1.14 | 1.61 / 1.43 |
| speed | 1.29 / **0.75** | 1.29 / 0.76 | 1.29 / 0.76 | 1.38 / 0.22 |
| attenuation | 1.51 / **0.94** | 1.52 / 0.93 | 1.51 / 0.93 | 1.46 / 0.92 |

- The orders **fall** as dt is refined, towards one. Pre-asymptotic
  higher-order terms would make them rise towards two; an emerging O(dt)
  term (or a floor) makes them fall.
- Moving the window from periods 5-6 to 11-12 changes no order by more than
  0.04 and no value by more than 2.3e-5 (wall amplitude at n50); at n100 the
  change is below 1e-6. The analysis window is not the cause (E).
- m1 (cheaper, 32 axial cells) shows the same behaviour: 1.51/1.22,
  2.55/-, 1.71/1.53, 1.49/1.16, 1.28/1.17, 1.20/0.97, and its successive
  differences are nearly those of m2 (wall amplitude 6.00e-3 and 1.83e-3 on
  m1, 6.29e-3 and 1.93e-3 on m2; speed 1.72e-3/7.09e-4 against
  1.76e-3/7.2e-4). The time error is almost mesh independent, so m1 was
  used for the cheap discriminating tests.
- A Robin m1 run at 800 steps per period failed at t = 34 s: the PETSc SNES
  diverged (`DIVERGED_DTOL`, floating-point exception) while the coupling
  residual was converging normally. Recorded as a solid-solver robustness
  failure at small dt on the coarsest mesh; not needed for the diagnosis.

### 6.2 Start-up and initial data (D)

The case already starts from the exact solution (section 2), so "start from
rest" is not the configuration studied. The initial data are, however, not
fully consistent: the old-time solid fields carried zero boundary values
(H2a). Largest residual of the wall displacement at L/2 after removing the
periods 11-12 fundamental, relative to the amplitude, in period 1:

| m2 | n50 | n100 | n200 | n400 |
|---|---:|---:|---:|---:|
| Robin, tutorial initial data | 1.8e-2 | 8.4e-3 | 1.5e-2 | **3.4e-2** |
| Robin, `exactOldBoundary` | 1.7e-2 | 4.2e-3 | 2.4e-3 | 2.6e-3 |
| IQN-ILS, tutorial initial data | 1.7e-2 | 4.3e-3 | 2.6e-3 | 2.5e-3 |

- **The reported "start-up transient that grows as dt is reduced" is
  explained**: it is the Robin coupling feeding the solid's *boundary*
  acceleration, the `backward` `d2dt2` of the boundary values of `D` and its
  old levels, to the fluid (`elasticWallPressure`, and the kinematically
  consistent interface flux in `pimpleFluid`). The old-time boundary values
  are the placeholder zero, so for the first three or four steps the imposed
  acceleration is of order |D|/dt^2. IQN-ILS does not use that acceleration
  and does not show the growth; exact boundary values remove it.
- It does not affect the analysis: `exactOldBoundary` changes periods 5-6
  by at most 2.5e-6 and periods 11-12 by about 1e-8, and leaves every order
  unchanged (table in 6.1). The free-mode content of the transient peaks at
  about 1.5 times the forcing frequency (the first standing pressure-wave
  mode of the finite tube, c/(2L) = 1.46 f), damped by about 0.4 per period.
- A ramped start from rest was therefore not run: the exact start is
  already available, a ramp would only lengthen the transient, and the
  window test shows that the transient is gone from the analysed periods.

### 6.3 Tolerances (G) and fluid inner iterations

| Change (periods 5-6) | flow_amp | wallMid_amp | wallMid_phase | speed | attenuation | coupling iterations/step |
|---|---:|---:|---:|---:|---:|---:|
| `tight` - base, m2 n200 | 1.3e-7 | -1.3e-6 | 1.5e-6 | -5.0e-7 | 5.2e-6 | 10.3 -> 15.5 |
| `tight` - base, m1 n400 | 2.6e-7 | 4.8e-7 | -1.5e-6 | 4.7e-8 | -1.2e-6 | 8.2 -> 16.9 |

The smallest successive difference at 200 -> 400 steps is 7.5e-5 (flow
phase); the tolerance effect is at least 15 times smaller, and the `tight`
m1 series has the same orders (1.51/1.22, 1.71/1.53, 1.49/1.16, 1.28/1.18,
1.20/0.98). `pimple4` reproduces the base to about 1e-8. Neither the
coupling and solver tolerances nor the PIMPLE splitting are the cause.

### 6.4 Coupling and interface conditions (F)

IQN-ILS (`newMovingWallVelocity`, zero-gradient pressure) and Robin-Neumann
(`elasticWallVelocity`, `elasticWallPressure`) show the same fall of the
order (6.1; m1 IQN-ILS: 1.71/1.09, 2.04/-, 1.63/1.31, 1.53/0.72, 1.29/0.91,
1.28/0.57) and differ by at most 2e-4 at n50. The behaviour is therefore
not specific to the Robin pressure condition (H5), and since the two
interface velocity conditions share their formulation (BDF2 of the face
centres plus the extrapolated mesh flux), the decisive test was to remove
the interface altogether.

### 6.5 Isolating the fluid (fluid-only bisection, m1, periods 5-6)

The fluid region alone, with the exact wall motion imposed, reproduces the
first-order behaviour. Each row removes or replaces one ingredient.
Orders of flow_amp from 50/100/200/400/800/1600 steps per period, and the
last difference:

| Fluid-only variant | flow_amp orders | d(800->1600) | profile orders |
|---|---|---:|---|
| moving mesh, exact wall motion (`fluidOnly`) | 1.60 1.33 1.13 **1.03** | 6.7e-5 | 0.53 0.55 0.83 0.74 |
| + 10 PIMPLE correctors (`pimple10`) | 1.60 1.33 1.13 1.02 | 6.7e-5 | 0.48 0.52 0.82 0.73 |
| static mesh, exact wall velocity (`transpiration`) | 1.70 1.47 1.26 **1.11** | 6.6e-5 | 1.73 1.44 1.43 1.06 |
| + constant `ddtCorr` coefficient (`ddtCoeff1`) | 1.73 1.49 1.27 1.14 | 6.2e-5 | 1.79 1.75 1.71 1.63 |
| + no `ddtCorr` (`noDdtCorr`) | 1.67 1.17 0.77 0.63 | 1.6e-4 | 2.13 1.60 - - |
| + `leastSquaresS4f` (`lsS4f`) | 1.70 1.47 1.26 1.12 | 6.6e-5 | 1.73 1.63 1.72 1.52 |
| + zero-gradient end velocity (`endsZeroGrad`, not exact) | **2.05 2.35** - - | 8.7e-5 (floor) | 2.00 1.64 1.95 1.37 |

- Mesh motion, ALE, the moving-wall condition, the solid and the coupling
  are not needed: the static-mesh problem with Dirichlet wall data has the
  same first-order term (6.6e-5 against 6.7e-5).
- Not the Rhie-Chow `ddtCorr` (H12): a constant coefficient leaves the
  orders unchanged, and switching it off makes them lower.
- Not the boundary gradient reconstruction: `leastSquaresS4f` (full
  boundary delta vectors) changes the QoIs by less than 3e-5 and not the
  orders. (The tube mesh is orthogonal at its boundaries, so
  non-orthogonal boundary corrections are not expected to matter here.)
- The **end velocity condition** is implicated: replacing the exact
  gradient by zero gradient restores second order of the flow rate and the
  profile until a floor of about 1e-4 is reached.
- The term is **mesh independent**, which rules out a Rhie-Chow-type
  O(dt h^2) momentum-interpolation error (H9): `transpiration`
  d(800->1600) of flow_amp is 6.6e-5 (m1), 8.0e-5 (m2) and 7.7e-5 (m4).

### 6.6 Mechanism, and direct measurement

At the tube ends the tutorial imposes the exact pressure (`codedFixedValue`)
and the exact normal velocity gradient (`codedMixed` with
`valueFraction 0`). In OpenFOAM, `mixedFvPatchField::assignable()` returns
`false` regardless of the value fraction, so `constrainHbyA` sets
`HbyA_b = U_b` on these patches, as for a fixed value. The pressure
equation then gives the end flux

    phi_b = U_b . S_f - rAU_b snGrad(p) |S_f|,

and since the pressure is fixed on the same patch, `snGrad(p)`, the axial
pressure gradient that drives the wave, does not vanish. The continuity
flux through the ends therefore differs from the flux of the boundary
velocity used by the momentum equation by `rAU dp/dn |S_f|`, with
`rAU ~ dt/1.5`: an O(dt) inconsistency, independent of the mesh and of the
coupling. Measured in the converged fields (final time, same phase in every
run; mean |phi_b - U_b . S_f| over the patch / largest |U_b . S_f|, leaving
out the end faces that touch the wall: on the moving mesh `phi` is relative
to the mesh motion, and those faces have a corner point that moves axially
with the solid):

| Run (outlet) | n50 | n100 | n200 | n400 | n800 | n1600 |
|---|---:|---:|---:|---:|---:|---:|
| fluid-only static, m1 |  | 3.1e-02 | 1.8e-02 | 1.0e-02 | 5.3e-03 | 2.7e-03 |
| fluid-only static, m2 |  | 1.9e-02 | 1.4e-02 | 8.8e-03 | 5.1e-03 | 2.7e-03 |
| fluid-only static, m4 |  | 7.3e-03 | 6.4e-03 | 5.1e-03 | 3.7e-03 | 2.3e-03 |
| **coupled Robin, m2 (the study)** | 2.2e-02 | 1.9e-02 | 1.4e-02 | 9.2e-03 |  |  |
| coupled IQN-ILS, m2 | 2.2e-02 | 1.9e-02 | 1.4e-02 | 9.3e-03 |  |  |
| coupled Robin, m2, `fluxConsistent2` (final fix, 6.10) | 2.3e-04 | 1.8e-04 | 1.3e-04 | 8.7e-05 |  |  |
| coupled Robin, m2, `fluxConsistent` (first fix, 6.10) | 1.6e-05 | 1.0e-05 | 1.0e-05 | 1.0e-05 |  |  |
| fluid-only static, m2, `fluxConsistent2` |  | 2.9e-03 |  | 2.4e-05 |  | 2.5e-05 |
| fluid-only, exact end velocity + `fixedFluxPressure` |  | 1.0e-10 |  | 1.0e-10 |  | 1.0e-10 |

(Inlet values are similar or smaller; all in `end_flux_inconsistency.csv`.)
Without the fix the inconsistency tends to halve with dt; at large dt on
fine meshes it falls more slowly because the viscous part of the momentum
diagonal then limits `rAU`. It is present, at about 1-2% of the end flux,
in the production study for both couplings.

### 6.7 Consistent alternatives: is second order recovered?

The exact end condition (pressure plus normal velocity gradient) has no
flux-consistent realisation in OpenFOAM's segregated pressure-velocity
algorithm:

1. `codedMixed` (the tutorial): end flux off by `rAU dp/dn`, O(dt) (6.6).
2. `fixedGradient` (`endsFixedGradient`, assignable, same continuum
   problem): `HbyA_b` is then extrapolated from the cell and ignores the
   imposed gradient, so the end flux misses `(U_b - U_P).S_f = O(h du/dn)`,
   and the pressure correction absorbs it with weight `1/rAU ~ 1/dt`.
   Fluid-only, m1 and m2: the pressure error roughly doubles at each
   halving of dt and the flow phase drifts by about 1e-3 per halving
   (worse than first order); unchanged with `noDdtCorr`. In the coupled
   Robin case (m1, 50/100/200/400): flow_amp 1.56/3.09, flow_phase
   0.46/-0.06, wallMid_amp 2.18/1.74, wallMid_phase 1.33/0.34, speed
   2.04/2.23, attenuation 1.95/2.42; on m2: flow_amp 1.78/1.07, flow_phase 0.66/-0.19,
   wallMid_amp 2.05/2.69, wallMid_phase 1.59/0.47, speed 1.91/1.83,
   attenuation 2.12/2.42. Changing *only* the end treatment therefore moves
   the coupled amplitude and wave-number QoIs to second order, while the
   phases now drift by the `1/dt` mechanism. Neither form is consistent.
3. Exact end velocity with `fixedFluxPressure` (`endsDirichletU`): the
   end flux is consistent to 1e-10, and the first-order term disappears.
   Fluid-only, flow_amp over 50 ... 1600 steps per period stays within
   2.1e-3 to 2.8e-3 (m1) and 6.3e-4 to 8.4e-4 (m2), and the change from 50
   to 100 steps falls from 2.9e-3 (tutorial ends, m2) to 1.1e-4; the
   gauge-free wave number from the axial velocity changes by at most
   3e-4 (m2) (against 2.6e-3 from 50 to 100 with the tutorial ends). But
   the pressure level is then set only by a reference cell (the pressure
   QoIs are meaningless, and in the coupled problem the wall traction
   would be wrong), and the windows wander by about 1e-4 (m2) to 5e-4
   (m1) from period to period, so no clean order can be measured below that
   level. It shows that the *time* error of the scheme itself is about 1e-4
   or less at 50 steps per period, but it is not a usable verification
   configuration.

None of these end conditions, available without a code change, is both
exact and flux consistent. Section 6.10 adds the missing piece in the fluid
solver.

### 6.8 Mesh factor 4 (C)

Robin, m4 (128 axial cells), tutorial end conditions, 100/200/400 steps per
period, 8 periods (about 2.5, 4 and 6 h on one core):

| Quantity | n100 | n200 | n400 | order (5-6) | order (7-8) | m2 order |
|---|---:|---:|---:|---:|---:|---:|
| flow_amp | -1.195e-3 | -7.167e-4 | -5.357e-4 | 1.40 | 1.41 | 1.12 |
| flow_phase | -2.518e-4 | 1.360e-5 | 9.431e-5 | 1.72 | 1.73 | 1.79 |
| wallMid_amp | 2.541e-3 | 8.611e-4 | 2.724e-4 | 1.51 | 1.52 | 1.31 |
| wallMid_phase | -2.833e-3 | -1.875e-3 | -1.531e-3 | 1.48 | 1.47 | 1.14 |
| speed | 1.401e-3 | 8.957e-4 | 6.287e-4 | 0.92 | 0.93 | 0.75 |
| attenuation | -6.884e-3 | -4.470e-3 | -3.479e-3 | 1.28 | 1.28 | 0.94 |

Refining the mesh does not restore second order: spatial contamination does
not explain the temporal order. The m4 differences are somewhat smaller
than on m2 (e.g. attenuation 9.9e-4 against 1.9e-3 from 200 to 400 steps),
consistent with the slower fall of the end-flux inconsistency on finer
meshes, where the viscous part of the momentum diagonal also limits `rAU`
(6.6).

### 6.9 Post-processing (H8)

- The wall displacement and the flow rate are sampled every step and their
  Fourier coefficients over whole periods are exact for the fundamental;
  their time stamps are the new time (a one-step lag would give a phase
  error of omega dt = 0.13 rad at n50, not observed).
- The fluid sets are sampled ten times per period at exact times; harmonics
  that could alias onto the fundamental are 9 and 11 times the frequency
  and negligible (the residual after removing the fundamental is 1.1e-3 of
  the amplitude in every run, independent of dt).
- Alternative extractions on the same runs: a complex-exponential fit of
  p (`speed_alt`, `attenuation_alt`) and the wave number of the axial
  velocity (`speed_ux`, `attenuation_ux`) show the same first-order
  behaviour as the log-pressure fit in the original set-up (e.g.
  `speed_ux`, m1 static: 1.69/1.54/1.36/1.29) and lose it with consistent
  ends. The low orders are not an artefact of the measurement; the
  reported definitions were not changed.

### 6.10 Fix: flux-consistent end patches (`fluxConsistentPatches`)

Added after the first version of this report (commit f0bab2730), on request,
as an opt-in option of `pimpleFluid` (OpenFOAM.com form,
`pimpleFluid.esi.C`):

```
PIMPLE
{
    ...
    fluxConsistentPatches (inlet outlet);
}
```

The listed patches must have a velocity condition that fixes the boundary
value (fixedValue, or mixed such as the tutorial's `codedMixed`, so that
`constrainHbyA` sets `HbyA_b = U_b`) and a fixed pressure; anything else
stops with a fatal error. The option is empty by default, so no other case
changes, and the end conditions themselves (exact pressure, exact normal
velocity gradient) are unchanged. A pure boundary condition cannot do this,
because `phiHbyA` is assembled in the solver and, on an assignable patch,
`HbyA_b` is the extrapolated `UEqn.H()` that a velocity condition cannot
set.

**Implementation (final).** On the listed patches `HbyA` is given the
pressure gradient of the boundary cell, as it has in the cell:

    phiHbyA_b = (U_b + rAtU_b grad(p)_P) . S_f,

so that the pressure equation gives

    phi_b = U_b . S_f + rAtU_b (grad(p)_P . S_f - snGrad(p) |S_f|),

the boundary velocity flux up to `O(rAtU h d2p/dn2)` (about 60 times
smaller than the original `rAtU dp/dn` on m2, and vanishing with the mesh),
while the fixed pressure stays an implicit Dirichlet condition of the
pressure equation.

**First implementation (superseded, recorded; never committed).** It set
`phiHbyA_b = U_b . S_f + rAtU_b snGrad(p)|S_f|` with `snGrad(p)` from the
previous iterate, which makes the converged end flux exactly `U_b . S_f`:

```
phiHbyA.boundaryFieldRef()[patchI] =
    (U().boundaryField()[patchI] & mesh().Sf().boundaryField()[patchI])
  + rAtU.boundaryField()[patchI]
   *p().boundaryField()[patchI].snGrad()
   *mesh().magSf().boundaryField()[patchI];
```

(in place of the final `(U_b + rAtU_b grad(p)_P) . S_f`; the results labelled
`fluxConsistent` in `temporal/` were produced with it).
It gave the same orders on m1 and m2 (variant `fluxConsistent`; Robin m2
2.03/1.86/1.96/2.14/1.94/2.10 from 50/100/200 and 2.00/1.90/1.96/2.09/
1.92/2.08 from 100/200/400; fluid-only m2 flow_amp 2.61 2.04 1.96 1.96),
but it failed twice: IQN-ILS m2 n50 diverged in the first time-step (PETSc
floating-point exception after the quasi-Newton iteration diverged), and
Robin m4 n200 stagnated at a coupling residual of about 2e-2 in the first
step (also with exact old-time boundary values). At convergence the fixed
end pressure no longer enters the end flux, so the pressure level near the
ends is held only through the momentum gradient, and the lagged correction
becomes unstable as the mesh is refined. The final implementation keeps the
pressure implicit and converges in both cases (m4 n200: 14 iterations in
the first step, as without the option; IQN-ILS m2 n50: at most 14).

Results with the final implementation (variant `fluxConsistent2`; private
build; the m2 and m4 coupled runs on 5 and 7 MPI ranks, decomposed into
axial slabs that keep processor boundaries away from the sampling
stations, see 6.11):

- Fluid-only static problem, orders 50/100/200/400/800/1600:

| m2 | tutorial ends | `fluxConsistent2` |
|---|---|---|
| flow_amp | 1.79 1.42 1.05 0.90 | **2.10 2.06 2.02 2.04** |
| flow_phase | 1.42 1.51 1.46 1.42 | **1.96 2.02 1.95 1.82** |
| speed | 0.21 -0.72 0.04 0.33 | **2.05 2.23 2.47 2.66** |
| attenuation | 0.62 0.60 0.38 0.96 | **2.75 2.76 2.71 2.92** |
| speed_ux | 1.79 1.39 1.05 0.92 | **2.24 2.19 2.13 2.18** |

  (m1: flow_amp 2.02 2.01 2.02 1.99.)

- **Coupled Robin-Neumann, m2 (the production study):**

| Quantity | tutorial 50/100/200 | tutorial 100/200/400 | **fix 50/100/200** | **fix 100/200/400** | fix, periods 11-12 |
|---|---:|---:|---:|---:|---:|
| flow_amp | 1.64 | 1.12 | 2.04 | **2.03** | 2.04 / 2.05 |
| flow_phase | 1.85 | 1.79 | 1.85 | **1.87** | 1.85 / 1.91 |
| wallMid_amp | 1.71 | 1.31 | 1.96 | **1.95** | 1.96 / 1.97 |
| wallMid_phase | 1.73 | 1.14 | 2.14 | **2.08** | 2.14 / 2.06 |
| speed | 1.29 | 0.75 | 1.94 | **1.93** | 1.94 / 1.97 |
| attenuation | 1.51 | 0.94 | 2.10 | **2.08** | 2.11 / 2.07 |

  Errors with the fix (periods 5-6; n50, n100, n200, n400): flow_amp
  -1.229e-3, 5.151e-4, 9.406e-4, 1.045e-3; flow_phase 1.430e-4, 1.074e-3,
  1.332e-3, 1.403e-3; wallMid_amp 5.283e-3, -6.653e-4, -2.196e-3,
  -2.592e-3; wallMid_phase -6.118e-3, -2.342e-3, -1.486e-3, -1.283e-3;
  speed 1.408e-3, -6.79e-5, -4.534e-4, -5.549e-4; attenuation -1.307e-2,
  -4.234e-3, -2.177e-3, -1.690e-3. The pressure and profile errors become
  independent of dt from n100 (changes 2e-6 to 2e-5). Coupling iterations
  9.2-11.0 per step (8.7-10.3 without the fix); period-to-period change
  <=7e-5.
- Coupled Robin, m1 (first implementation): 2.03/1.99, 1.85/1.85,
  1.96/1.95, 2.14/2.09, 1.94/1.90, 2.07/2.17.
- IQN-ILS against Robin-Neumann with the fix, m2 n100: at most 2.3e-4
  (wall amplitude), as without it (3.0e-4); 15.2 against 9.2 coupling
  iterations per step.

The time error at the study's 200 steps per period, estimated with order two
from the n200 -> n400 difference (e/3), is now 1.3e-4 (wall amplitude),
3.4e-5 (speed) and 1.6e-4 (attenuation), against 1.6e-3, 0.9e-3 and 3.8e-3
(order one) with the tutorial ends.

**Mesh study with the fix** (Robin, 200 steps per period, periods 5-6;
original values from PR #513 for comparison):

| QoI | m1 | m2 | m4 | mesh order | original m4 | original order |
|---|---:|---:|---:|---:|---:|---:|
| profile | 1.790e-2 | 4.960e-3 | 1.203e-3 | 2.04 | 1.18e-3 | 2.06 |
| flow_amp | 6.100e-3 | 9.406e-4 | -3.083e-4 | 2.05 | -7.15e-4 | 2.16 |
| flow_phase | 6.048e-3 | 1.332e-3 | 8.33e-5 | 1.92 | 1.5e-5 | 1.92 |
| wallMid_amp | -1.031e-2 | -2.196e-3 | -1.787e-4 | 2.01 | 8.57e-4 | 2.15 |
| wallMid_phase | -2.372e-3 | -1.486e-3 | -1.148e-3 | 1.40 | -1.88e-3 | 0.97 |
| speed | -2.333e-3 | -4.534e-4 | 1.24e-5 | 2.01 | 8.93e-4 | 2.48 |
| attenuation | -2.076e-3 | -2.177e-3 | -1.961e-3 | - | -4.47e-3 | 0.63 |

(Profile order from the errors of the two finest meshes, the others from
successive differences, as in `Allverify`.) The amplitude, flow and
wave-speed QoIs converge at second order, with smaller finest-mesh errors
(the wave speed to 1.2e-5). The wall phase improves from first order to
1.40. The attenuation error no longer converges at a low order: it is
about -2e-3 on all three meshes, a level that no longer depends on the
mesh or (6.10) the time-step, consistent with the linear-theory limit of
about 1e-3 (README) amplified about seven times in Im(k); it was about half
of the original finest-mesh attenuation error, the rest being the end-flux
time error.

### 6.11 Parallel sampling artefact

With scotch decomposition on 4 and 8 ranks, processor boundaries fell
exactly at x = L/4, L/2 and 3L/4: the axial pressure set then contains
duplicate points (merged by the driver; serial results are unchanged), and,
more seriously, the flow rate through the `midPlane` faceZone at L/2 is
wrong (flow phase error -1.58e-2 instead of 1.07e-3 at m2 n100), while the
wall, wave-speed and attenuation QoIs are unaffected. The offset does not
depend on dt, so it does not change the orders, but the values are
unusable. Parallel runs therefore use axial slabs (`simple`) with a number
of ranks (3, 5 or 7) that keeps every processor boundary at least 0.5 m
from the sampling stations; the driver refuses other counts. With 5 slabs
the flow-phase error at m2 n100 is 1.074e-3 (serial 1.079e-3).

## 7. Version difference

v2412 against v2512, m2, n100 (periods; relative errors):

| Period | Robin wallMid_amp | Robin flow_amp | IQN-ILS wallMid_amp | IQN-ILS attenuation |
|---|---:|---:|---:|---:|
| 1 | -1.7e-3 | -3.3e-4 | 1.3e-6 | - |
| 2 | 5.4e-4 | 7.1e-5 | -3.4e-5 | 9.3e-5 |
| 3 | -1.6e-4 | 3.4e-5 | -1.1e-5 | 9.9e-5 |
| 4 | 4.5e-5 | 1.1e-5 | -1.7e-5 | 1.2e-4 |
| 6 | 4.5e-6 | 5.3e-7 | -1.6e-5 | 1.2e-4 |
| 12 | 1.3e-6 | 1.5e-7 | -1.6e-5 | 1.1e-4 |

- **Reproduced exactly**: the regression value u_r(t = 25) is
  -2.229508061e-4 (v2512) and -2.217864674e-4 (v2412) for Robin, and
  -2.236660656e-4 / -2.237665584e-4 for IQN-ILS, to all ten digits of the
  PR #513 comment.
- **Robin-Neumann sensitivity is start-up only.** The two versions agree at
  the first step and separate from the second; the difference decays by
  3-4 times per period and is 4.5e-6 in period 6. The regression samples
  t = 25 s, half a period into the start-up, which is why it sees a 0.26%
  shift. The source is the only functional difference between the two
  installations in the time discretisation: the installed v2412 has the
  solids4foam `optionalFixes` `backwardDdtScheme.C`, whose Euler fall-back
  tests the old-time indices of each field, whereas stock v2512 falls back
  to Euler only when `timeIndex() < 2`. They differ in which `ddt`,
  `meshPhi` and boundary-acceleration terms are Euler during the first
  steps, which matters for Robin because it feeds the solid boundary
  acceleration, built from the (zero) old-time boundary values, to the
  fluid (6.2). (Identified by code inspection and by the start-up-only
  signature; not isolated by a rebuilt OpenFOAM.)
- **IQN-ILS sensitivity is persistent but small**: about 1.1e-4 in the
  attenuation and up to 6e-5 in the other QoIs, constant from period 4 on
  and the same at n100 and n200 on m1 (so not an O(dt) term). Its source
  was not identified; it is below every successive time difference used and
  below every finest-mesh error in the PR #513 table (the flow-phase
  offset, 9e-6, is below the finest-mesh flow-phase error, 1.5e-5).

## 8. Most likely explanation, ranked by evidence

1. **O(dt) mass-flux inconsistency on the tube ends (established).** The
   tutorial's end condition, exact pressure with an exact normal velocity
   gradient through `codedMixed`, is treated by OpenFOAM as fixing the
   velocity, so the continuity flux through the ends differs from the
   boundary velocity flux by `rAU dp/dn`, proportional to dt (6.6). It is
   present in every coupled and fluid-only configuration that keeps these
   end conditions, independent of the mesh, the coupling, ALE, `ddtCorr`,
   the gradient scheme, the tolerances and the start-up, and it disappears
   with a flux-consistent end condition (6.7), and in particular with the
   end flux made equal to the boundary velocity flux, which restores second
   order in the coupled study (6.10). The observed orders are the
   transition from the BDF2 O(dt^2) error to this O(dt) term: 1.3-1.8 from
   50/100/200 and 0.75-1.8 (mostly about 1) from 100/200/400.
2. **Start-up transient (established, but not the cause of the order).**
   Robin feeds the solid boundary acceleration built from zero old-time
   boundary values to the fluid; the period-1 transient grows as dt falls
   (3.4e-2 at n400), and is removed by exact boundary values. It is gone
   from the analysis window (<3e-6).

## 9. Hypotheses ruled out

| # | Hypothesis | Verdict and evidence |
|---|---|---|
| H1 | BDF2 first-step treatment | Ruled out: windows 5-6 and 11-12 agree; transient gone by period 5 |
| H2a | Zero old-time boundary values | Real start-up defect (explains the dt-growing transient), not the cause of the order (<3e-6 in the window, orders unchanged) |
| H2b | Initial `pointD` mismatch | Not the cause: IQN-ILS has no dt-growing transient; any first-step mesh jump is gone by period 5 |
| H3 | Ramp | Not applicable (exact start) |
| H4 | Interface velocity time accuracy | Ruled out: the static-mesh, Dirichlet-wall fluid problem has the same O(dt) term |
| H5 | Robin pressure condition | Ruled out: IQN-ILS and Robin behave the same; fluid-only (no Robin) too |
| H6 | ALE / GCL | Ruled out: static mesh has the same term (`backward` meshPhi is GCL-consistent by construction, `1.5 phi^{n+1/2} - 0.5 phi^{n-1/2}`) |
| H7 | Predictor | Ruled out with H10 (converged results unchanged by tolerances) |
| H8 | Extraction artefact | Ruled out: alternative extractions behave the same; windows irrelevant |
| H9 | Space-time interaction | Ruled out as cause: the O(dt) term is mesh independent (m1, m2, m4) |
| H10 | Tolerance floor | Ruled out: 100x tighter tolerances move QoIs by <5e-6 |
| H11 | Boundary data time level | Ruled out by inspection: coded conditions use the new time |
| H12 | OpenFOAM detail | Confirmed in a different form: `mixedFvPatchField::assignable() == false` with `constrainHbyA` (6.6); `ddtCorr` and the `timeIndex() < 2` Euler fall-back ruled out as causes of the order |

## 10. Is nominal second order recovered?

Without a code change, no (6.7). With `fluxConsistentPatches` on the tube
ends (6.10), **over the range tested, yes; formally, at fixed mesh, no**
(section 15). Every QoI of the production study (m2, Robin-Neumann, same
exact end data) converges at 1.87-2.08 from 100/200/400 steps per period and
1.85-2.14 from 50/100/200, in both analysis windows. The change does not
alter the discretisation of the equations; it removes a demonstrated O(dt)
inconsistency of the end flux. The boundary-flux term that remains is
O(dt h): the scheme is formally second order under joint refinement, and at
fixed mesh second order holds until the BDF2 error falls to the size of
that term (about 3200 steps per period on m1, 6400 on m2), beyond which the
order tends to one. The exact-flux form has no such term but is not robust
(15.6).

## 11. Statement for the paper

> With the exact travelling-wave data imposed on the artificial tube ends
> as a fixed pressure and a prescribed normal velocity gradient, the
> observed temporal orders fell from 1.3-1.8 (50/100/200 steps per period)
> to about one (100/200/400), on every mesh (m1, m2, m4). The cause was not
> the fluid-solid coupling: the behaviour was identical for IQN-ILS and
> Robin-Neumann, persisted in a fluid-only computation with the exact wall
> motion on a static mesh, and was independent of the mesh, the start-up,
> the analysis window and the tolerances. It was an O(Δt) inconsistency of
> the mass flux through the ends: the segregated pressure-velocity
> algorithm treats the mixed velocity condition as fixing the value, so the
> end flux differed from the boundary velocity flux by the momentum
> coefficient (∝ Δt) times the axial pressure gradient. With the end flux
> made consistent with the boundary velocity flux (an option of the fluid
> solver, leaving the boundary data unchanged), every quantity converges at
> second order in time over the steps used (1.87-2.08 from 100/200/400 steps
> per period), and the time error at 200 steps per period falls by an order
> of magnitude. The remaining boundary-flux inconsistency is O(Δt h), the
> same order as OpenFOAM's standard treatment of a fixed-pressure outlet:
> the scheme is second order under joint refinement, while at fixed mesh the
> order would tend to one only at much smaller time steps (beyond about
> 6400 steps per period on this mesh).

The exact reference is what made this boundary inconsistency visible: the
orders are reference-free, but the fluid-only and end-flux diagnostics, and
the confirmation that the fixed solution converges at second order towards
the exact one, are not.

With the fix, the mesh study (6.10) gives second order for the profile,
flow rate, wall amplitude and wave speed (2.01-2.05; the flow phase 1.92),
1.40 for the wall phase, and an attenuation error of about -2e-3 that no
longer depends on the mesh or the time-step (the linear-theory limit,
amplified in Im(k)). The finest-mesh errors fall to 1.2e-3 (profile),
3.1e-4 (flow amplitude), 1.8e-4 (wall amplitude), 1.15e-3 (wall phase),
1.2e-5 (wave speed) and 2.0e-3 (attenuation); the flow-phase error,
8.3e-5, is larger than the original 1.5e-5, which was smaller than the
time error it contained (probably a cancellation between the spatial error
and the end-flux time error). In the original set-up,
about half of the finest-mesh attenuation error and most of the wave-speed
error were the end-flux time error.

## 12. Bugs

- **No solids4foam source bug** was found in the coupling, the interface
  conditions, the solid, the time schemes or the mesh motion. The
  `fluxConsistentPatches` option fills a gap (no flux-consistent way to
  combine a fixed pressure with a prescribed velocity gradient) rather than
  fixing an error in existing code.
- **Verification set-up defects** (demonstrated, in the tutorial files, not
  in the library):
  - the end velocity condition (`codedMixed`, `valueFraction 0`, with a
    fixed pressure on the same patch) does not realise the intended
    gradient condition in the continuity flux: the end flux differs from
    the boundary velocity flux by an O(dt) term (6.6). This follows from
    OpenFOAM's treatment of mixed patches as non-assignable, which is a
    design choice of OpenFOAM rather than a defect, but it makes this
    combination time-inconsistent;
  - the boundary values of the old-time solid displacements are zero
    instead of the exact solution, contrary to the README's "exact old-time
    levels"; this produces the Robin start-up transient that grows as dt
    falls (6.2).
- Not a bug, recorded: the Robin m1 n800 run stopped with a PETSc SNES
  divergence (6.1).

## 13. Recommended final verification configuration

1. Add `fluxConsistentPatches (inlet outlet);` to the tutorial's
   `system/fluid/fvSolution` PIMPLE dictionary once the option is merged,
   and regenerate the regression references (both versions) and the stored
   results in the verification README; the numbers in 6.10 are those this
   configuration gives.
2. Fix the old-time boundary values in `0/solid/D*` (the `exactOldBoundary`
   change) at the same time. It removes the Robin start-up transient that
   grows as dt falls and changes the analysed periods by less than 3e-6;
   it also moves the regression value at t = 25 s.
3. Keep the six-period, last-two-periods analysis, the QoI definitions and
   the tolerances; none needed changing.
4. Extend the time study to 400 steps per period (50/100/200/400), and run
   any parallel cases with axial slabs that avoid the sampling stations
   (6.11).
5. Report the attenuation as converged to the linear-theory level
   (-2e-3), not as a first-order quantity.

## 14. Remaining open items

- The persistent IQN-ILS v2412/v2512 offset (1e-4 in the attenuation).
- The ~1e-4 floor of the self-convergence and Dirichlet-velocity fluid-only
  variants (not present with `fluxConsistent`).
- `fluxConsistentPatches` is implemented for OpenFOAM.com only (the tutorial
  runs only there).

## 15. Addendum: formal order of `fluxConsistentPatches`

Question: at fixed mesh, is the final treatment (`HbyA_b = U_b + rAtU_b
grad(p)_P`) formally second order in time, or does it keep an
O(dt h) boundary-flux error that eventually makes the time order one?

**Answer: it keeps an O(dt h) boundary-flux term.** At fixed mesh it is
formally first order in time as dt -> 0, with a coefficient proportional to
h; it is second order under joint refinement (dt proportional to h) and in
practice over the tested range, and it has the same boundary consistency as
OpenFOAM's standard fixed-pressure outlet. "Restores second order" must be
qualified accordingly (15.8).

### 15.1 Derivation

From `pimpleFluid::evolve()` (OpenFOAM.com form), with `consistent` off so
that `rAtU = rAU`:

- `UEqn.A()` is `D()/V`, the diagonal including the implicit boundary
  coefficients, with `extrapolatedCalculated` patches, so on the end patch
  `rAU_b = rAU_P = 1/A_P`.
- The pressure equation `laplacian(rAU, p) = div(phiHbyA)` gives on a
  fixed-value pressure patch (orthogonal end faces, no non-orthogonal
  correction) `pEqn.flux()_b = rAU_P |S_f| snGrad(p)`, with
  `snGrad(p) = (p_b - p_P)/d`, d the centre-to-face distance (h/2 for the
  axial cell size h).
- After the corrector, `U = HbyA - rAU grad(p)`, so at PISO convergence
  `HbyA_P = U_P + rAU_P g_P`, with `g_P = fvc::grad(p)_P` (the same
  `gradp()` the implementation uses).

With `phiHbyA_b = (U_b + rAU_P g_P) . S_f` (the implementation; `U_b` and
`g_P` from the previous corrector, equal to the current ones at
convergence):

    phi_b = phiHbyA_b - pEqn.flux()_b
          = U_b . S_f + rAU_P |S_f| (g_P . n - snGrad(p)).          (1)

Equivalently, since `HbyA_P - U_P = rAU_P g_P`,
`phiHbyA_b = (HbyA_P + U_b - U_P) . S_f`: HbyA is extrapolated to the face
with the increment that the velocity condition gives U. For a zero-gradient
velocity (`U_b = U_P`) this is exactly OpenFOAM's standard treatment of a
fixed-pressure outlet (an assignable velocity condition, `HbyA_b` the
extrapolated cell value), which therefore carries the same residual (1).
The original `codedMixed` treatment instead gives
`phi_b = U_b . S_f - rAU_P |S_f| snGrad(p)`.

### 15.2 Scaling of rAU with dt

`A_P = gamma/dt + a_nu` (+ the negligible convection), with gamma = 1
(Euler, first step) or 3/2 (BDF2) and `a_nu ~ nu sum |S_f| delta_f / V`
(about 1.05, 4.2 and 17 s^-1 in the end cells of m1, m2 and m4). Hence

    rAU_P = dt / (gamma + a_nu dt) = dt/gamma - a_nu dt^2/gamma^2 + ...

`A` is the same in every PIMPLE outer corrector (the problem is linear). So
`rAU = O(dt)` once `dt < gamma/a_nu` (about 35, 140 and 560 steps per period
on m1, m2 and m4); at larger steps it is limited by viscosity, which is why
the end-flux inconsistency fell slowly at 100-400 steps on m4.

### 15.3 Spatial order of `g_P . n - snGrad(p)`

Applying OpenFOAM's own operators (`leastSquares` gradient, the
fixed-value `snGrad`) to the *exact* pressure (values set at the cell
centres and on the end faces), on every end face:

| Mesh | h (axial) | (g_P.n - snGrad p) / (h d2p/dn2) | max misfit | max abs value |
|---|---:|---:|---:|---:|
| m1 | 0.469 | -0.244 | 0.5% | 2.34e-6 |
| m2 | 0.234 | -0.247 | 0.3% | 1.18e-6 |
| m4 | 0.117 | -0.249 | 0.1% | 5.94e-7 |

So `g_P . n - snGrad(p) = -(h/4) d2p/dn2 + O(h^2)`. The least-squares
gradient is accurate at the cell centre; the fixed-value `snGrad` is the
normal gradient at the midpoint between the cell centre and the face, h/4
away. The residual is a first-order *location* mismatch. (On interior faces
both are centred on the face and the corresponding Rhie-Chow difference is
O(h^2).) From (1),

    phi_b - U_b . S_f = -(dt/gamma) (h/4) |S_f| d2p/dn2 + O(dt h^2, dt^2 h).   (2)

The original treatment has `-(dt/gamma) |S_f| dp/dn`: the ratio is
`(h/4) |d2p/dn2| / |dp/dn| = (h/4)|k|` for the travelling wave
(|k| = 0.1455 m^-1), about 1/59, 1/117 and 1/234 on m1, m2 and m4.

### 15.4 Formal order of the time discretisation at fixed mesh

The boundary-flux error (2) is a source of mass on the end faces,
proportional to dt, which the solution responds to linearly. At fixed h the
solution therefore has an error

    e(dt) = B dt^2 + C h dt + ...,

with B the BDF2 coefficient and C independent of dt and h (proportional to
the pressure curvature on the ends):
**formally first order in time at fixed mesh**, with a coefficient that
vanishes as O(h). Under joint refinement (dt proportional to h) the term is
O(h^2), so the space-time scheme is formally second order. The order seen
at fixed mesh is two while `B dt > C h`, and tends to one for
`dt < C h / B`.

(By the same argument the interior momentum interpolation carries an
O(dt h^2) term, `rAU (interp(grad p) - snGrad p)`, one order higher in h
than the boundary term. It was not isolated here: with exact-flux ends the
fluid sub-problem converges at second order to 12800 steps per period on
m1 (15.5), so any such term is below about 1e-7 there.)

### 15.5 Numerical evidence

Fluid-only static problem (exact wall velocity, tutorial end data, m1),
periods 5-6, 200 ... 12800 steps per period; `scripts/womersley_boundary_order.py`.

**The term itself.** `D = Q(fluxConsistent2) - Q(fluxConsistent)` (complex
flow-rate coefficient) isolates (2): the two runs differ only in it.

| steps/period | 200 | 400 | 800 | 1600 | 3200 | 6400 | 12800 |
|---|---:|---:|---:|---:|---:|---:|---:|
| \|D\| | 1.25e-5 | 6.49e-6 | 3.37e-6 | 1.72e-6 | 8.68e-7 | 4.36e-7 | 2.19e-7 |

Halving ratios 1.93, 1.93, 1.96, 1.98, 1.99, 1.99 (constant `ddtCorr`
coefficient; 1.90-1.96 with the default): **exactly first order in dt**.
`|D/F|`, with F the first-order part of the original (`codedMixed`)
solution, tends to 0.027 on m1 and is still falling (0.069, 0.036, 0.022,
0.0175 from n200 to n1600) on m2, consistent with the scaling with h in (2)
(on m2 the viscous part of `rAU` matters to smaller dt, 15.2).

**The order it produces.** Successive-difference orders of the complex
flow-rate coefficient, with a constant `ddtCorr` coefficient (`backward 1`;
see the note below):

| triplet (steps/period) | 200/400/800 | 400/800/1600 | 800/1600/3200 | 1600/3200/6400 | 3200/6400/12800 |
|---|---:|---:|---:|---:|---:|
| exact end flux (`fluxConsistent`) | 2.00 | 2.01 | 2.02 | 2.04 | **2.08** |
| final form (`fluxConsistent2`) | 2.02 | 2.03 | 2.05 | 1.97 | **1.57** |

and for the flow phase 1.94, 1.98, 2.00, 2.03, 2.08 against 1.77, 1.70,
1.58, 1.42, **1.27**; the wave speed and attenuation of the final form
change sign in their differences beyond 1600 steps, where the O(dt h) term
crosses the BDF2 term. The BDF2 error of the exact-flux form at 3200 steps
is 9.1e-7, against `|D|` = 8.7e-7: on m1 the crossover is at about 3200
steps per period (by (2), about 6400 on m2 and 12800 on m4). The coupled
study (to 400 steps per period on m2) is well inside the second-order
range, which is why it shows 1.87-2.08.

**Spatial order of the residual** (15.3): `-(h/4) d2p/dn2` from OpenFOAM's
operators on the exact pressure, on m1, m2 and m4.

**Note: the default `ddtCorr` limits fixed-mesh convergence at about
1e-6.** With OpenFOAM's default Rhie-Chow `ddtCorr` coefficient,
`1 - min(|phi - U_f.S|/|phi|, 1)`, both implementations stop converging
below differences of about 1e-6 from 3200 steps per period (e.g. the
profile differences stay at -3.6e-6 per halving), which hides the turn to
first order. The floor is not the PIMPLE iteration (10 outer correctors
change nothing), nor the start-up transient (identical at every step); it
disappears with a constant coefficient. The coefficient switches as the
oscillating face fluxes pass through zero, so it is not smooth in time. This
is an interior effect of OpenFOAM's standard momentum interpolation, two
orders of magnitude below the differences used in the coupled study.

### 15.6 The first implementation (exact end flux)

`phiHbyA_b = U_b . S_f + rAU_P |S_f| snGrad(p^{k-1})`, with `snGrad` from
the previous corrector, gives at convergence `phi_b = U_b . S_f` exactly: no
dt-proportional term (it shows second order where the second form shows
the first-order term, 15.5). But with exact flux and an exact pressure on
the same patch, the pressure equation has two boundary conditions for one
second-order equation. The converged state is that of a *Neumann*
(prescribed-flux) pressure problem on the ends, in which the imposed
pressure value enters only through the momentum gradient `grad(p)` of the
end cells. The implicit Dirichlet coefficient is cancelled by the lagged
explicit term, so the iteration is a deferred correction whose contraction
depends on how strongly the pressure level is otherwise fixed (through a
Robin interface, weakly; with a zero-gradient interface, not at all). That
is why it stagnated on m4 and diverged with IQN-ILS at 50 steps. Making the
cancellation implicit (removing the patch from the pressure equation)
turns the ends into a genuine Neumann condition and leaves the pressure
level undetermined with IQN-ILS. There is therefore no well-posed treatment
that both imposes the flux exactly and keeps the fixed pressure as an
implicit Dirichlet condition of the pressure equation: the data are
over-specified for the projection, and one of the two can only be
satisfied to the accuracy of the discretisation.

The clean choice is to keep the pressure implicit and make the flux
consistent to the order of the pressure equation's own boundary gradient:

- the present form, `G = g_P . n`, gives O(dt h);
- a location-consistent form, `G = g_P . n + (d/2) (d2p/dn2)_P`, i.e. the
  cell gradient carried to where `snGrad` is centred, gives O(dt h^2), the
  same order as the interior faces. It needs a curvature estimate from the
  interior cells only (e.g. the normal derivative of `grad(p)` between the
  boundary cell and its interior neighbour), so that the fixed pressure
  stays implicit; using `snGrad` itself for the curvature reproduces the
  first implementation. Not implemented or tested here.

### 15.7 Where the treatment belongs

- OpenFOAM's convention is that a velocity condition on a fixed-pressure
  patch is *assignable* (`inletOutlet`, `pressureInletOutletVelocity`,
  `fixedNormalInletOutletVelocity` are mixed conditions that override
  `assignable()` to true), so `HbyA_b` is extrapolated, and the residual (1)
  is accepted. The tutorial's `codedMixed` (non-assignable) with a fixed
  pressure is outside that convention, which is what produced the O(dt)
  term.
- An assignable `fixedGradient` alone is not enough for a non-zero gradient:
  `HbyA_b` is then the cell value and misses `U_b - U_P`, and the pressure
  absorbs the mismatch with weight 1/rAU (the `endsFixedGradient` failure).
  The consistent rule is `HbyA_b = HbyA_P + (U_b - U_P)` on such patches,
  which the present implementation reproduces at convergence and which
  reduces to the standard treatment when `U_b = U_P`.
- Recommendation: a **specialised capability for over-specified
  (exact/manufactured) boundary data**, keyed to the velocity boundary
  condition rather than to a solver patch list: a dedicated assignable
  velocity condition that imposes a prescribed normal gradient, recognised
  by the solver, which applies `HbyA_b = HbyA_P + (U_b - U_P)` on it. This
  removes the possibility of pairing the option with an inconsistent
  condition (the patch list only checks the pairing at start-up). It should
  not be made a general change to all assignable patches: for inflow
  segments of `inletOutlet` it would change OpenFOAM's standard behaviour.
  Ordinary outlets (zero-gradient or `inletOutlet` velocity with a fixed
  pressure) already have the same O(dt h) consistency and need nothing.

### 15.8 Wording

Replace "restores second order" by: "removes the O(dt) end-flux
inconsistency; the remaining boundary term is O(dt h), so the time order is
two over the tested range (to 400 steps per period in the coupled study, to
about 1600-3200 steps per period on the fluid sub-problem) and the scheme is
formally second order under joint refinement, but at fixed mesh the order
tends to one as dt -> 0, with a coefficient that vanishes with h".

## 16. Final verification from the production path

The two changes are now in the tutorial (`fluxConsistentPatches (inlet
outlet)` in `system/fluid/fvSolution`; exact boundary values of `D`, `D_0`,
`D_0_0`, `D_0_0_0` through `womersleyPatchDisplacement` in
`system/womersleyCode`), and `Allverify` runs the time study at
50/100/200/400 steps per period and checks the order of the finest three
(threshold 1.7). The full `./Allverify --cores 7` (OpenFOAM-v2512, Linux,
serial cases, a clean build of the branch) passes every check;
`temporal/production/` holds its summary and CSV.

| Quantity | mesh order | time order 50/100/200 | 100/200/400 | finest-mesh error (m4, n200) |
|---|---:|---:|---:|---:|
| profile | 2.04 | - | - | 1.20e-3 |
| flow_amp | 2.05 | 2.04 | 2.05 | -3.09e-4 |
| flow_phase | 1.92 | 1.86 | 1.91 | 8.23e-5 |
| wallMid_amp | 2.01 | 1.96 | 1.97 | -1.76e-4 |
| wallMid_phase | 1.40 | 2.14 | 2.06 | -1.15e-3 |
| speed | 2.02 | 1.94 | 1.97 | 1.40e-5 |
| attenuation | - (-2.08e-3, -2.18e-3, -1.96e-3) | 2.10 | 2.07 | -1.96e-3 |

- IQN-ILS against Robin-Neumann (m2, n100): at most 2.16e-4 (wall
  amplitude); 14.6 and 9.3 coupling iterations per step.
- Periodicity (change from the fifth to the sixth period): at most 6.6e-5.
- First-period transient (largest residual over the amplitude, m2): 1.6e-2,
  2.7e-3, 2.2e-3, 2.7e-3 at 50, 100, 200, 400 steps per period, against 1.8e-2,
  8.4e-3, 1.5e-2, 3.4e-2 before the old-time boundary values were exact: it
  no longer grows as dt falls.
- These agree with the investigation runs (6.10): the mesh-study values to
  two or three digits, and the time orders to within 0.04. The differences
  come from the investigation's m2 runs having been decomposed on five
  ranks and run without exact old-time boundary values (which changes the
  analysed periods by less than 3e-6).
- Regression (`regressionTest.sh`, u_r at t = 25 s, 50 steps): new
  references -2.241790827e-4 m (IQN-ILS) and -2.242599027e-4 m (Robin) from
  v2512, which changed by -5.1e-7 and -1.3e-6 m with the two tutorial
  changes. v2412 now differs by 1.1e-7 (IQN-ILS) and 1.2e-8 m (Robin): the
  former 1.2e-6 m Robin difference came from the start-up and has gone, so
  the tolerance between versions is reduced from 2e-6 to 5e-7 m. Passes on
  both versions.
- The tutorial's default run (IQN-ILS, two periods) is now within 0.5%
  (profile), 0.05% and 0.0011 rad (flow rate), 0.08% (wall amplitude), 0.02%
  (wave speed) and 0.41% (attenuation) in the second period; the tutorial
  README is updated.
