# Robin coefficient study: results (2026-09-13)

OpenFOAM v2412, branch `feature/robin-auto-hs` (based on PR #450), serial
unless noted. All runs use `fixedRelaxation` with relaxation factor 1 unless
noted, and `robinFluxTolerance 1` (see README). "Mean/max" are FSI iterations
per time step. For `cerebralAneurysm` the iteration counts are capped by inner
solver floors (the solid solver skips solves once the load change falls below
its tolerance), so the iterations needed to bring both the displacement and
pressure residuals below 1e-3 are reported instead.

Variants:

- `pWaveSpeed`: the original default, `hs = c_p dt_eff`;
- `thicknessLimited`: `hs = l tanh(H_eff/l)`, ray-cast thickness, halved when
  wetted on both sides;
- best fixed: best `thicknessLimited` `hsScale` from a sweep (0.5-3);
- secant median/minimax: `hsModel secant; secantUpdate iteration;` with
  `secantFit median` or `minimax` (defaults otherwise).

## Summary table

Columns: `pWaveSpeed` (pWS), `thicknessLimited` (TL), best fixed (`hsScale`
in brackets), secant median (med) and secant minimax (mm).

| case | pWS | TL | best fixed | med | mm |
|---|---|---|---|---|---|
| 3dTube dt | 4.57/6 | 5.05/6 | 4.57 (1.03x) + | **4.08/6** | 4.12/6 |
| 3dTube 4dt | 52.4/57 | 6.78/10 | 6.78 (1x) | **5.90/11** | 6.43/10 |
| 3dTube 16dt | crash | 21.6/24 | 8.87/11 (2x) | **7.27/12** | 7.67/12 |
| beam modified | 52.2/194 | 50.0/96 | 27.4/52 (2x) | 28.8/38 | **24.0/34** |
| beam original | 9.80/12 | 108.6/143 | 61.0/79 (2x) | 9.45/22 | **9.05/22** |
| container | crash | 5.83/19 | 5.83 (1x) | **5.55/20** | 6.42/30 |
| aneurysm | crash | 7.5 | 5.8 (1.5x) | **5.2** | 8.9 |
| dam break | 11.5/22 | 14.1/26 | 11.3/25 (2x) | **10.9/21** | 11.4/21 |
| Turek-Hron | crash | 13.9/100 * | 13.9/100 * | crash | 13.9/100 * |

Cases (time steps): `3dTube` with the tutorial time step dt (60), 4 dt (40)
and 16 dt (15); `beamInCrossFlow` modified (20) and original (20);
`fillingElasticContainer` (2000); `cerebralAneurysm` (10), FSI iterations to
1e-3; `flexibleDamBreakRobin` (~400, adaptive dt); `HronTurekFsi3Robin` (300
coupled).

`+` best fixed choice: `pWaveSpeed` (1.03x).

`*` most steps reach `nOuterCorr = 100`: the iterations contract at about 0.95.

References: `fillingElasticContainer` as shipped (IQN-ILS, `constantHs 0.1`)
needs 8.75 / 30; with fixed relaxation 1 and `constantHs 0.1`, 8.70 / 27.
`cerebralAneurysm` as shipped (`constantHs 5e-4`) needs 5.7 iterations to 1e-3.

Accuracy: converged variants agree with each other to within 5e-6 relative
(`3dTube`), 4e-5 (`beamInCrossFlow` original), 4e-4 displacement / 8e-3 force
(`beamInCrossFlow` modified), 3e-3 (`fillingElasticContainer` apex, all within
the regression tolerance), 5e-4 (`cerebralAneurysm`) and 2% (dam break, whose
adaptive time steps differ between runs). The `beamInCrossFlow` regression test
passes for all four variants.

## Findings

1. **The theory holds quantitatively.** The measured solid impedance (from the
   interface pressure and acceleration changes between iterations) does not
   depend on the Robin coefficient used, and the predicted contraction factor
   `(S_f/S_s)|S_s - alpha|/(S_f + alpha)` matches the observed one (3dTube:
   0.025-0.042 predicted vs 0.031-0.043 observed; Turek-Hron: 0.91-0.96 vs
   0.956). The toy model predicts the 3dTube optimum and divergence limit
   (for 4 dt: optimum 1.34x the tanh estimate, divergence above 2.3x; observed:
   1x best, 2x diverging; for 16 dt: optimum 3.2x; observed 2-3x best).
2. **`pWaveSpeed` is not a safe default.** It crashes in 4 of 9
   configurations (thin or nearly incompressible walls, large time steps),
   because `c_p dt_eff` exceeds the wall thickness by 2-7x and alpha > 2 S_s
   diverges for strong added mass. It is only the best choice for the stiff,
   weakly added-mass beam, where the solid impedance is stiffness dominated.
3. **`thicknessLimited` is safe but not optimal.** It never crashed: the
   slab estimate `rho_s l tanh(H/l)` is a lower bound of the solid impedance
   (free back face, no stiffness). It is 1.5-11x slower than optimal when the
   stiffness term matters (large dt, soft or stiff beams, grid-scale modes
   with impedance about `rho_s l^2/dx`).
4. **The online secant estimate is the best automatic choice**, but only
   with per-iteration updates. Updating once per time step lags one step and
   oscillates with period two between the soft and stiff modes (a quiet mode
   for the current coefficient becomes the unstable one after the change).
   This was fixed by (a) skipping the first iteration pair of each step, which
   also contains the physical load change of the new step, (b) keeping impedance
   samples from recent steps, (c) limiting increases more than decreases, and
   (d) per-iteration updates.
5. **median vs minimax.** The median fit is fastest on average (best in 6 of 9)
   but crashed on Turek-Hron by tracking the stiff mode. The minimax fit
   (minimise the largest predicted rate over the recorded samples) never
   crashed and was within 20% of the best in 8 of 9, at worst 1.7x
   (`cerebralAneurysm`).
6. **The divergence safeguard only helps the static models** (it rescued
   `pWaveSpeed` at 16 dt) and caused oscillations with the secant model, so it
   is now ignored for `hsModel secant`.
7. **Spatially varying secant scale**: no gain (cerebralAneurysm 8.0 vs 5.2
   iterations to 1e-3); the spatial variation from the ray-cast wall thickness
   is enough.
8. **Limit of a scalar Robin coefficient.** In Turek-Hron the flap has a soft
   bending mode with `S_f/S_s ~ 9` (needs alpha ~ 10) and a stiff mode with
   `S_s ~ 520` (needs alpha ~ 520). The best single alpha gives a worst rate of
   0.89, which is what the minimax fit predicts and reaches. Adding IQN-ILS on
   top of the Robin condition did not help: it crashed or was unchanged on
   Turek-Hron and the modified beam, because IQN-ILS accelerates the interface
   displacement while the Robin pressure/acceleration state is not included in
   its unknowns.
9. **PR #450 criteria.** The leakage-flux residual levels off at a floor
   (0.019 in 3dTube, 0.026 in cerebralAneurysm, 0.036 in Turek-Hron) and the
   pressure residual can stall at inner-solver floors, so the default
   tolerances (equal to `outerCorrTolerance`) force `nOuterCorr` iterations in
   every step regardless of `hs`.
10. **Robustness fixes found on the way:** the solid density is now taken from
   the solid model (the `rho` field is not registered before
   `couplingStartTime`, which made the original lookup fail), rays skip their
   own warped origin face, and the Robin number report uses global reductions.

## Recommendation and next steps

- Use `hsModel secant; secantFit minimax; secantUpdate iteration;` (seeded by
  `thicknessLimited`) as the automatic choice, and `thicknessLimited` as the
  a-priori fallback. Consider making the secant model the default once the
  regression references are updated.
- Revisit the PR #450 flux residual: measure its change between iterations
  (like the pressure residual) rather than its absolute value.
- For problems like Turek-Hron, look beyond a scalar coefficient: a two-level
  (coarse-mode) Robin coefficient, or quasi-Newton acceleration that includes
  the interface pressure.
- Not yet verified: OpenFOAM-9 and foam-extend builds (the ray cast uses
  `triSurfaceSearch` on OpenCFD versions and a brute-force search elsewhere),
  larger parallel runs (the modified beam on 4 cores reproduces the serial
  iterations: 24.0 mean), and mesh-refinement studies.

## Part 2: interface leakage and Robin convergence criteria

### Why the converged Robin interface leaks

Tukovic et al. (2018) derive the interface volume flux from the Euler
momentum balance at the wall, Eqs. (28)-(30): `(1/a_P)_bi = dt`,
`(H/a_P)_bi = v_bi^(m-1)` and `Vdot = n.v^(m-1) S - dt (dp/dn) S`, with the
normal velocity BC taken from that flux, Eq. (31), and the Robin condition,
Eq. (20), using the transferred solid acceleration. At convergence Eq. (17)
gives `Vdot = S (v_n^(m-1) + dt a_s)`, which equals the mesh (GCL) flux only if
(a) the old fluid interface velocity equals the old wall velocity, and (b) the
swept-volume flux of the fluid mesh equals `S v_s,n` in the same discrete
kinematics as the solid face acceleration. The implementation followed the
paper (`Uf.oldTime()`, whose normal component is `Vdot/S`), but (b) does not
hold discretely: the fluid mesh is moved by vertex displacements while `a_s`
is a face-centre solid quantity mapped between meshes, faces change area and
orientation, and with the backward scheme `fvc::meshPhi` is a multi-level
combination. Because Eq. (29) carries the fluid's own old velocity forward,
each step's mismatch was carried into the next step (measured carry-over at
step n+1 = leak at step n), giving the leakage floors seen with PR #450
(0.019 in 3dTube, 0.04 in cerebralAneurysm, relative to the throughput).

Dirichlet-Neumann has no such error: the wall velocity sets the flux to the
mesh flux (3dTube: relative wall flux 4e-22, machine zero).

### Fixes

1. `robinFluxFromWallVelocity` (pimpleFluid, default on): build the Robin
   flux from the old wall velocity (set by `elasticWallVelocity` from the mesh
   motion) instead of the old fluid face velocity, removing the carry-over.
   interFluid already did this.
2. `robinKinematicConsistency` (pimpleFluid and interFluid, default on):
   build the explicit (old-velocity) part of the Robin interface flux from the
   mesh flux, `phiHbyA_b = meshPhi_b - rAU_b a_s |S|` (interFluid:
   `meshPhi_b + rAU_b |S| (phig_b - rho_b a_s)`), where `a_s` is the solid
   acceleration used by the Robin condition. The Robin condition itself, and
   with it the solid inertia and the added-mass coupling, is unchanged; the
   interface flux is `meshPhi + rAU (a_f - a_s) |S|`, which equals the mesh
   flux exactly once the iterations converge (`a_f = a_s`).

   A first version instead replaced the Robin right-hand-side acceleration by
   the acceleration implied by the fluid mesh motion. It gave the same leakage
   reduction but changed the Robin iteration: the mesh-derived acceleration is
   smoother (slower 3dTube convergence) and, with under-relaxation or IQN-ILS,
   follows the relaxed displacement, which made the fillingElasticContainer
   tutorial (IQN-ILS, relaxation 0.1) diverge. The flux-side form avoids both.

Leakage relative to the wall-motion flux at convergence:

| case | original | wall-velocity start | + kinematic consistency |
|---|---|---|---|
| 3dTube | 1e-3 to 0.25 | 1e-3 (0.25 first step) | ~1e-6 |
| beamInCrossFlow modified | ~1e-3 | ~1e-3 | 1e-6 to 1e-4 |
| cerebralAneurysm | ~5% | 5% -> 1.4% (decaying) | 2e-4 to 9e-4 |
| flexibleDamBreakRobin (interFluid) | 4.6e-4 (throughput) | - | 6.7e-9 |

With kinematic consistency the leakage is proportional to the Robin
iteration error, like Dirichlet-Neumann (whose leakage is exactly zero); it
reaches Dirichlet-Neumann levels only if the FSI iterations are converged to
machine precision. Iteration counts are unchanged except 3dTube (4.2 -> 6.3):
the mesh-derived acceleration is smoother (vertex interpolation), so
high-wavenumber interface modes get weaker Robin feedback. In 3dTube the
result moves towards the Dirichlet-Neumann solution (force difference 0.46% ->
0.024%).

The remaining 2% difference between Robin and Dirichlet-Neumann in the
modified beam is in the pressure force (viscous forces agree to 0.05%) and is
unaffected by the leakage: the Robin interface flux uses the inertial wall
momentum balance (`rAU_b = dt`, `phiHbyA_b = U_old.S`), dropping the
convective and viscous parts of `HbyA` at the wall face, which changes the
wall pressure gradient by a first-order (half-cell) term in this viscous
case. This is a separate discretisation question.

### Convergence criteria

PR #450 required the displacement residual, the relative interface pressure
change and the absolute leakage (normalised by the throughput) to satisfy
`outerCorrTolerance`. Problems found:

- the leakage floor (above) made every 3dTube step run to `nOuterCorr` (100),
  and most flexibleDamBreakRobin steps;
- the pressure criterion was applied before `couplingStartTime`, so every
  uncoupled Turek-Hron step ran to `nOuterCorr`;
- the pressure change is normalised by the interface pressure norm, which is
  dominated by the mean pressure in cerebralAneurysm (~13 kPa), so it is
  lenient there; with kinematic consistency the leakage is the sharper
  measure.

Final design (`robinConvergence residual`, default): displacement residual
`<= outerCorrTolerance`, pressure change `<= robinPressureTolerance` (default
`outerCorrTolerance`) and, if the Robin condition is kinematically consistent,
leakage normalised by the throughput plus the interface motion flux `<=
robinFluxTolerance` (default `10*outerCorrTolerance`); otherwise the leakage
is only reported. Before `couplingStartTime` only the displacement residual
is used. The residual file gains a convergence-state column.

An iteration-error criterion (`robinConvergence iterationError`) was also
implemented and compared: it estimates the remaining error of each residual
as `R rho/(1 - rho)` from the contraction rate over the iterations in which
the solid solution changed, requires it when the iterations contract slowly
(rate >= 0.9), stops early when it satisfies the tolerance, and accepts
stalled residuals within 10 times the tolerances. It helped where the
iterations contract quickly and cleanly (3dTube -17% iterations with errors at
the tolerance, stiff beam -7%) but inner-solver noise (interFluid pressure
solves, the solid solver in cerebralAneurysm skipping its solve once its own
tolerance is met, which produces residual plateaus) made the rate estimates
unreliable and added iterations (cerebralAneurysm 21 vs 17, container 7.2 vs
6.8) with no accuracy gain, so it is optional.

Mean FSI iterations per step (steps reaching `nOuterCorr` in brackets) and
final interface leakage, normalised by the throughput plus the interface
motion flux (the Robin interfaces are excluded from the throughput), with
the PR #450 criteria (absolute leakage normalised by the
throughput, tolerance `outerCorrTolerance`; run with the wall-velocity flux
start) and with the final defaults of this work (`hsModel secant`,
`robinKinematicConsistency`, `robinConvergence residual`):

| case | PR #450 criteria | final defaults | final leakage |
|---|---|---|---|
| 3dTube | 100 (60 of 60) | 5.87 | 8e-8 |
| beamInCrossFlow modified | 23.9 | 24.3 | 1e-6 |
| beamInCrossFlow original | 9.05 | 9.0 | 8e-8 |
| cerebralAneurysm | 30 (10 of 10) | 17.3 | 1.3e-4 |
| flexibleDamBreakRobin | 94.5 (384 of 410) | 11.8 | 5e-9 |
| container, relaxation 1 | 21.4 (1262 of 2000 unconverged) | 6.34 | 1.1e-4 |
| container, as shipped | - | 8.67 (0) | reported only |
| Turek-Hron, before coupling | 100 | 1 | - |
| Turek-Hron, coupled | 100 | 100 (scalar-coefficient limit) | 5e-6 |

Container: `fillingElasticContainer` with fixed relaxation 1, or as shipped
(IQN-ILS, relaxation 0.1). Turek-Hron: `HronTurekFsi3Robin`; the coupled steps
reach `nOuterCorr` because of the scalar-coefficient limit (see part 1).

In 3dTube the final defaults move the solution towards the
Dirichlet-Neumann one: maximum displacement difference 0.49% -> 0.12%, force
difference 0.46% -> 0.23%.

With under-relaxation or interface acceleration (Aitken, IQN-ILS) the fluid
mesh follows the relaxed displacement during the iterations; building the
interface flux from that mesh flux changed the iteration map and made the
container tutorial diverge, so the flux-side consistency is only applied when
the fluid interface is moved to the solid interface in every iteration
(`fluidMeshFollowsSolid()`: fixed relaxation with factor 1, the recommended
Robin setting). Otherwise the paper's flux is used and the leakage is only
reported.

Regression tests (OpenFOAM-v2412): 3dTube, beamInCrossFlow (all four
variants) and cerebralAneurysm pass with their existing references. The
fillingElasticContainer reference apex displacement was updated from -0.479
to -0.4822: the displacement-only criterion accepted most steps after about
one FSI iteration (mean 1.26), and the converged coupling (mean 8.7) gives
-0.4822 (-0.4817 to -0.4826 across coupling methods and Robin coefficient
models); the shift comes from the Robin convergence criteria of PR #450.

### Remaining limitations

- The flux-side consistency is only implemented for the OpenFOAM.com
  pimpleFluid and interFluid; the OpenFOAM.org and foam-extend variants keep
  the paper's formulation (leakage reported only), and they have not been
  compiled in this study.
- The start-up (first time step) mismatch and the first-order wall-pressure
  difference in the viscous beam (above) are not addressed.
- A scalar Robin coefficient cannot converge problems with widely separated
  interface modes quickly (Turek-Hron).
