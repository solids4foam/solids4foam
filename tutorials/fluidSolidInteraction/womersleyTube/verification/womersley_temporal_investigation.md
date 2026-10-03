# womersleyTube: investigation of the sub-nominal temporal order

Status: **hypotheses written before any run** (section 2). Sections 3 onwards
are filled in as the diagnostics complete.

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
