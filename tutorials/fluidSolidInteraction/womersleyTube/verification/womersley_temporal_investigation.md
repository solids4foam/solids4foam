# womersleyTube: investigation of the sub-nominal temporal order

Status: **complete.** Sections 1-4 are the code trace and the hypotheses as
written before any run (commit 8f8c8468e), unchanged. Sections 5-14 are the
results. Short answer: the sub-nominal order is caused by an O(dt)
mass-flux inconsistency on the artificial tube ends, where the exact
pressure and the exact normal velocity gradient are imposed (section 8); it
is not an FSI, interface, ALE, solid, start-up or tolerance effect.

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
run; mean |phi_b - U_b . S_f| over the patch / largest |U_b . S_f|):

| Run | n50 | n100 | n200 | n400 | n800 | n1600 |
|---|---:|---:|---:|---:|---:|---:|
| fluid-only static, m1, outlet | | 3.3e-2 | 2.0e-2 | 1.1e-2 | 5.8e-3 | 3.0e-3 |
| fluid-only static, m2, outlet | | 2.0e-2 | 1.4e-2 | 9.1e-3 | 5.3e-3 | 2.8e-3 |
| fluid-only static, m4, outlet | | 7.3e-3 | 6.4e-3 | 5.2e-3 | 3.7e-3 | 2.4e-3 |
| **coupled Robin, m2 (the study), outlet** | 2.2e-2 | 1.8e-2 | 1.3e-2 | 9.1e-3 | | |
| coupled IQN-ILS, m2, outlet | 2.2e-2 | 1.8e-2 | 1.3e-2 | 9.2e-3 | | |
| fluid-only, exact end velocity + `fixedFluxPressure` | | 9.8e-11 | | 9.8e-11 | | 9.8e-11 |

(Inlet values are 3-5 times smaller; `end_flux_inconsistency.csv`.) The
inconsistency tends to halve with dt; at large dt on fine meshes it falls
more slowly because the viscous part of the momentum diagonal then limits
`rAU`. It is present, at about 1-2% of the end flux, in the production
study for both couplings.

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

Second order is therefore **not recovered under any configuration that
keeps the exact end data and is usable for the coupled problem**. With the
consistent but gauge-free alternative, the temporal error at 50 steps per
period drops by more than an order of magnitude, to the 1e-4 floor of that
set-up.

### 6.8 Mesh factor 4 (C)

Robin, m4, 100/200/400 steps per period, 8 periods: **running at the time of
this commit**. The runs cost far more than estimated (about 2.5 h, 4 h and
6 h on one core). They are not needed for the diagnosis: the O(dt) term is
already shown to be mesh independent on the coupled m1 and m2 series (6.1)
and on the fluid-only static problem on m1, m2 and m4 (6.5: 6.6e-5, 8.0e-5,
7.7e-5), and on m4 the end-flux inconsistency is still present (6.6). The
results will be added when the runs finish.

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
   with a flux-consistent end condition (6.7). The observed orders are the
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

No, not in any configuration that keeps the exact end data and remains
usable for the coupled problem (6.7). With a flux-consistent but gauge-free
end condition (fluid only), the first-order term disappears and the
remaining time dependence is at the 1e-4 level already at 50 steps per
period, which bounds the BDF2 time error of the fluid scheme but does not
give a measurable order.

## 11. Statement for the paper

> In womersleyTube the observed temporal orders at fixed mesh fall from
> 1.3-1.8 (50/100/200 steps per period) to about one (100/200/400). The
> cause is not the fluid-solid coupling: the behaviour is identical for
> IQN-ILS and Robin-Neumann, persists in a fluid-only computation with the
> exact wall motion and on a static mesh, and is independent of the mesh,
> the start-up, the analysis window and the solver and coupling
> tolerances. It is an O(Δt) inconsistency of the mass flux on the
> artificial tube ends, where the exact pressure and the exact normal
> velocity gradient are imposed: OpenFOAM treats the mixed velocity
> condition as fixing the value, so the end flux carries a term
> proportional to the momentum coefficient (∝ Δt) times the axial pressure
> gradient. With a flux-consistent end treatment the time dependence at 50
> steps per period falls by more than an order of magnitude. The time
> study of this case therefore verifies the boundary treatment of the
> truncated domain rather than the second-order accuracy of the coupled
> scheme, and we do not claim temporal second order from it.

The mesh orders are unaffected (successive differences at fixed dt cancel
the common time error), but the finest-mesh errors against the exact
solution at 200 steps per period contain this time error: extrapolating the
m2 sequence with order one, the time error at n200 is about twice the
200 -> 400 difference, e.g. 1.6e-3 in the wall amplitude, 0.9e-3 in the
speed and 3.8e-3 in the attenuation, comparable to or larger than the
finest-mesh errors in `tab:womersley`.

## 12. Bugs

- **No solids4foam source bug** was found in the coupling, the interface
  conditions, the solid, the time schemes or the mesh motion.
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

1. Keep the production algorithms and the six-period, last-two-periods
   analysis; keep the mesh study at 200 steps per period (its orders are
   valid) and the IQN-ILS/Robin comparison.
2. Report the time study as a diagnosed boundary-limited result (section
   11), with the fluid-only and end-flux evidence; do not tighten
   tolerances or change the window to improve it (neither helps).
3. Fix the old-time boundary values in `0/solid/D*` (the `exactOldBoundary`
   change). It removes the dt-growing start-up transient and changes the
   analysed periods by less than 3e-6, but it moves the regression value at
   t = 25 s, so it needs new regression references on both versions; not
   applied here.
4. State in the finest-mesh table that the errors at 200 steps per period
   contain a time error of 1e-3 order from the end treatment.
5. If a temporal-order verification of the coupled scheme is wanted, it
   needs a flux-consistent exact outflow condition (for example an
   assignable gradient condition whose `HbyA` boundary value carries the
   imposed gradient, or exact Dirichlet velocity with `fixedFluxPressure` and
   a time-dependent pressure datum), which is a code change outside this
   investigation.

## 14. Remaining open items

- The persistent IQN-ILS v2412/v2512 offset (1e-4 in the attenuation).
- The ~1e-4 floor of the self-convergence and consistent-end fluid-only
  variants.
