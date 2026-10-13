# beamInCrossFlow: Richter (2012) 3-D steady FSI verification

This directory holds the opt-in verification campaign for the stationary
three-dimensional FSI benchmark of Richter (2012). It is not run by CI; the
tutorial regression test (`../regressionTest.sh`) is separate.

<!-- markdownlint-disable MD013 -->

## Campaign tools

| File | Purpose |
| --- | --- |
| `scripts/richter_case.py` | Builds one self-documenting case (`campaign_case.json`): mesh level `f`, solid thickness multiplier `m`, time step, coupling test, fluid/solid/mesh-motion solver tolerances, PIMPLE correctors, rigid-beam (fluid-only) mode, initialisation from a converged rigid-beam flow. |
| `campaign/run_case.sh` | Meshes, checks, decomposes and runs a case; continues from the latest time if the case is already decomposed. |
| `campaign/pack.slurm` | Runs several cases concurrently on MeluXina nodes (`srun --exact`). |
| `scripts/richter_analyse.py` | Quantities of interest, steadiness (changes over the last 1 and 2 s, exponential approach fit), coupling iterations, capped steps and the accepted `r1`/`r2` of every step. |
| `scripts/solid_discretisation.py` | Static solid-only diagnostic (uniform 100 Pa on the upstream face) with independent through-thickness and in-plane refinement. |
| `scripts/beam_in_cross_flow_verification.py`, `scripts/rigid_fluid_diagnostic.py` | Earlier drivers, kept for the graded-family definitions they export. |
| `reference/campaign/` | Compact evidence (CSV/JSON) of the campaign. |

All definitive runs use one platform and one build: MeluXina, OpenFOAM v2412
(foss-2024a, GCC 13.3), PETSc 3.22.0, and a frozen copy of the
`s4f-verif-v2412` build of solids4foam `aee0c8e350267` (-O2), stored in
`/project/home/p201287/philip/beam-richter-opus/user` (md5 sums identical to
the source build). The case files come from this branch; each case records
the commit that built it.

## Case-definition audit (Richter 2012, CMAME 223-224, section 7.2)

Re-read from the primary source (https://d-nb.info/1188567608/34, pp. 18-20).

| Item | Richter | Here | Note |
| --- | --- | --- | --- |
| Domain | `(0,1.5)x(0,0.4)x(-0.4,0.4)`, symmetric in the x-y plane, half computed | `z <= 0` half, `z = 0` symmetry | same |
| Solid | `(0.4,0.5)x(0,0.2)x(-0.2,0.2)`, StVK, `mu = 5e5`, `nu = 0.4` | StVK, `E = 1.4 MPa`, `nu = 0.4` | identical Lame constants |
| Clamp | `u = 0` on the base (y = 0), `n.u = 0` on the symmetry base | same | |
| Fluid | `rho = 1000`, `nu = 1e-3` | same | |
| Inlet | `v = 0.3 y(0.4-y)(0.4^2-z^2)/(0.2^2 0.4^2)`, "peak velocity vmax = 0.3" | `paraboloidInletVelocity`, `maxVelocity 0.3`: the same polynomial | see below |
| Outlet | do-nothing, `J sigma F^-T n = 0` | `p = 0`, `grad U . n = 0` | differs only by the normal viscous stress at an outlet 1 m downstream |
| Point A | `(0.45, 0.15, 0.15)` "close to the outer corner" | `(0.45, 0.15, -0.15)` | mirror point; a mesh vertex at every level; `u_z` changes sign |
| Drag | `J_drag = int_Gamma_i n.(J sigma_f F^-T) e1 ds` on the computed half | pressure + viscous traction on the deformed half interface (`forces`) | same functional; Richter evaluates it as a residual on the base |
| Problem | stationary monolithic ALE | transient-to-steady partitioned continuation | the steady state is shown to be independent of history and time step (below) |

Inlet speed. The printed polynomial and the sentence after it both give a
peak of 0.3 m/s. The next paragraph's "average inflow velocity of 0.2" and
`Re = 40` are inconsistent with it: the area mean of the printed bi-parabolic
profile is `(2/3)^2 0.3 = 0.133 m/s`; `0.2 = (2/3) 0.3` is the mean of a
single (2-D) parabola, as in the FSI-1 benchmark that this section follows.
We therefore read the "0.2" as a 2-D slip and use the printed profile
(peak 0.3), as decided for this study. Caveat recorded for the comparison: an
earlier solids4foam sequence with peak 0.2 m/s landed close to Richter's
Table 6 (`u_x(A)` L3 5.915e-5 vs 5.924e-5, `F_x` 1.331 vs 1.327), so Richter's
tabulated values may correspond to a peak of 0.2. Nothing here is tuned to
Table 6; the comparison is made only after our own sequence is resolved.

Steady interpretation. Richter solves the stationary equations; the solid
density does not enter. We march the transient equations (solid density
1000) to a steady state and accept values only when the steadiness and
history-independence checks below are met.

## Findings

### 1. The coupling-tolerance "floor" was a fluid-solver artefact

The October 10 production attempt (jobs 5323918-5323921, strict test
`max(r1, r2) <= 1e-6`, PIMPLE `nOuterCorrectors 1`) aborted at
`nOuterCorr 200` at every level with `r1` stuck at 1.6e-6 to 2.8e-5 while
`r2` was ~2e-11 (`reference/campaign/failed_r03_2026-10-10`). Per step,
IQN-ILS reduced `r1` to ~1e-4 in ~12 iterations and then entered a slow
linear tail (~0.6 % per iteration). With one PIMPLE corrector per FSI
iteration the fluid response depends on the iteration history (each FSI
iteration is one more PIMPLE iteration from the previous fluid iterate), so
the secant information is inconsistent: the same mechanism found for
Turek-Hron FSI3 (PR #472).

L1 diagnostic on MeluXina (`m = 4`, 0-5 s, strict 1e-6, `nOuterCorr 100`,
`reference/campaign/coupling_diagnostic_L1m4`):

| PIMPLE correctors | linear solvers | mean FSI it./step | capped steps | max accepted `r1` (last 1 s) | wall time 0-5 s (32 ranks) |
| ---: | --- | ---: | ---: | ---: | ---: |
| 1 | default (1e-6) | 24.2 | 6 | 5.5e-5 | 8876 s |
| 1 | tight (1e-9/1e-10) | 24.0 (to 3.2 s) | 1 | - | 9079 s (3.2 s only) |
| 3 | default | 8.4 | 0 | 9.6e-7 | 4977 s |
| 3 | tight | 8.1 | 0 | 9.7e-7 | 7409 s |

With three PIMPLE correctors the strict 1e-6 test is attained in 8-15
iterations at every step, including near the steady state, and the run is
1.8x faster. The `QoI` histories of all four variants agree to ~1e-6 in the
transient.

### 2. Loose fluid-solver tolerances leave history-dependent "steady states"

L0 (`m = 2`), workstation, strict 1e-6, PIMPLE 3. Started from rest (1 s
ramp) or from a converged rigid-beam flow (full inlet at t = 0):

| p/U tolerance | start | `F_x` (N) | `F_y` (N) | `u_x(A)` (m) |
| --- | --- | ---: | ---: | ---: |
| 1e-6 | rest | 2.158143 | 0.328095 | 8.139336e-5 |
| 1e-6 | rigid flow | 2.157357 | 0.326941 | 8.139184e-5 |
| 1e-6, motion 1e-12 | rest | 2.157361 | 0.326939 | 8.139200e-5 |
| 1e-6, motion 1e-12 | rigid flow | 2.158141 | 0.328108 | 8.139318e-5 |
| 1e-10, motion 1e-12 | rest (t = 20 s) | 2.157359 | 0.326939 | 8.139197e-5 |
| 1e-10, motion 1e-12 | rigid flow (t = 12 s) | 2.157357 | 0.326952 | 8.139180e-5 |

With the tutorial's tolerances (1e-6) a run freezes on one of two states that
differ by 3.6e-4 in `F_x`, 0.36 % in `F_y` (pressure only; the viscous force
is identical) and 1.9e-5 in `u_x(A)`; the interior ALE mesh differs by up to
13 % of the local displacement. Which state is reached is not tied to the
start. With p/U tolerance 1e-10 both starts agree to 1.4e-6 (`F_x`) and
2e-6 (`u_x`). Rigid-beam runs (no mesh motion) are unique to 2e-7 across
decompositions even with loose tolerances. The production settings therefore
use p/U 1e-10 and mesh motion 1e-12. The finding is local (workstation,
v2512) at L0: on MeluXina at L1 the loose-tolerance run (`P_L1_init_loose`)
agreed with the tight one to 1e-8, so the effect was not reproduced there; it
is reported as solver-tolerance-sensitive ALE trajectories leading to
different discrete equilibria (supported, not proven: a crossover test on a
common interior mesh was not run). The interior mesh differences were up to
1.6e-5 m, about 1e-3 of the local cell size. On MeluXina with the production
tolerances the rest and rigid-flow starts agree to 1.3e-6 at L0 and 1.6e-6
at L1.

### 3. Separating structural and fluid errors

Static solid-only diagnostic, standard solid, workstation (MeluXina rerun
running): `u_x(A)` (1e-4 m) under 100 Pa on the upstream face.

| cells through thickness `N_t` | 4 | 8 | 16 | 32 | 64 |
| --- | ---: | ---: | ---: | ---: | ---: |
| in-plane `8f`, `N_t = 4f` (original family) | 1.0522 | 1.2508 | 1.3566 | 1.3955 | - |
| in-plane `8f`, `N_t = 8f` (`m = 2`) | - | 1.2394 | 1.3510 | 1.3952 | 1.4087 |
| in-plane 8, `N_t` varied | 1.0522 | 1.2394 | 1.3404 | - | - |
| high-order solid, `N_t = 4f` | 1.2217 | 1.3699 | 1.4033 | - | - |

The displacement error is controlled by the through-thickness resolution;
in-plane refinement changes it by ~1 %. The standard and high-order sequences
extrapolate to 1.414-1.419e-4. The original family is 26/12/4.5/1.6 % too
stiff at L0-L3: this, not the fluid, made `u_x(A)` converge from below and
slowly in the October 0.2 m/s sequence. With `m = 2` the solid errors are
12.4/4.5/1.4/0.42 % (orders 1.34, 1.71). `m = 4` was rejected: it doubles the
maximum Courant number (thin fluid cells over the beam) and the rigid L3 run
became unstable at `Co ~ 3`.

Rigid-beam fluid (same graded fluid mesh, `m = 2`, dt 0.00625, t = 24 s):
`F_x` = 2.167855, 2.306750, 2.361026, 2.379550 N at L0-L3, changes
0.1389/0.0543/0.0185 N, observed orders 1.36 and 1.55. The steady rigid flow
is reached only after ~16-20 s (wake mode, tau ~ 2 s).

### 4. Time step

The steady state of the PIMPLE/Rhie-Chow discretisation depends weakly on
dt. Rigid `F_x`, t = 12 s: L1 2.307303/2.306835/2.306578 for dt
0.003125/0.00625/0.0125 (1-2e-4); L2 2.361104/2.361064 for 0.00625/0.0125
(1.7e-5). dt = 0.00625 s is held fixed at every level (max Courant 1.5 at
L3); the FSI dt check is `c1_L1_tol1e-6_dt0125`.

### 5. Coupling criterion (final) and its evidence

The relative measure `r1` normalises the residual by the largest residual of
the current step. As the steady state is approached the per-step change, and
with it that normaliser, goes to zero, so any fixed `r1` target eventually
falls below the solver noise floor. On L1 (from rest, cap 50) the strict test
needed 8-10 iterations per step until t ~ 10 s and then rose to the cap: at
1e-6, 29 capped steps after t = 15.5 s (`r1` floor 1.5e-6, `r2` ~ 1e-14); at
1e-7, 34 capped steps by t = 12 s; from the rigid flow, 67 capped steps with
`r1` up to 0.04 and `r2` <= 7e-9. At 1e-5 no step capped (max 18 iterations).

Production criterion: strict `max(r1, r2) <= 1e-5`, `nOuterCorr 30`,
`allowUnconvergedCoupling yes` as a counted safety net (capped steps are
listed with their `r1`, `r2` in `campaign_steps.csv`), PIMPLE 3 correctors,
p/U 1e-10, mesh motion 1e-12, SNES rtol/stol 1e-6.

Evidence that the coupling error is negligible (relative changes):

| test | `u_x(A)` | `F_x` |
| --- | ---: | ---: |
| L1 from rest, strict 1e-4 / 1e-5 / 1e-6 / 1e-7 at t = 16 s (12 s for 1e-7) | <= 1.0e-6 (1e-7: 2.4e-6, not yet steady) | <= 1.1e-6 |
| L1, continued from steady with strict 1e-8 (cap 100, every step capped, `r2` ~ 1e-14) | 3e-7 | 2e-6 |
| L1, continued with the legacy `min(r1, r2) <= 1e-6` for 6 s | +-1e-6 noise, no drift | +-3e-6 |
| L1, continued with PIMPLE 6 correctors | 6e-7 | 3e-6 |
| L1, SNES rtol/stol 1e-9 vs 1e-6 (from rest, t = 16 s) | < 1e-8 | < 1e-8 |
| L0, strict 1e-5 vs 1e-6 (rest t = 24 s, rigid-flow start t = 16 s) | <= 3e-7 | <= 4e-7 |
| L2, strict 1e-4 vs production (same start; t = 2 s / t = 8 s) | 1.9e-6 / 1.1e-7 | 1.6e-6 / 8e-8 |

The capped `1e-6`/`1e-7` runs are sensitivity evidence, not references; the
coupling error is below 2e-6 at every tested setting, two to three orders
below the L1-L2 mesh increments.

### 6. Time step

The steady state of this transient PIMPLE (Rhie-Chow) discretisation depends
on dt. It is reversible (L1: from the dt = 0.00625 steady state, dt = 0.0125
moves `u_x(A)` to 9.6168e-5 and the reverse branch returns to 9.6261e-5), so
it is a property of the discrete steady equations, not of the transient
path.

| level | dt = 0.0125 | dt = 0.003125 | (relative to dt = 0.00625) |
| --- | ---: | ---: | --- |
| rigid L1 `F_x` | -1.1e-4 | +2.0e-4 | t = 12 s |
| rigid L2 `F_x` | -1.7e-5 | +3.6e-5 | t = 12 s |
| rigid L3 `F_x` | - | +3.4e-6 | t = 24 s |
| FSI L1 `u_x(A)` | -9.6e-4 | +1.4e-3 | restarted from steady |
| FSI L1 `F_x` | -5.9e-4 | +8.7e-4 | restarted from steady |
| FSI L2 `u_x(A)` | -2.3e-4 | +3.9e-4 | restarted from steady |
| FSI L2 `F_x` | -1.4e-4 | +2.4e-4 | restarted from steady |

The dt effect falls by 6-10x per level for the rigid beam and by 3.6-4.2x from L1 to L2 for the FSI (so ~5e-5 expected at L3, ~2 % of the expected L2-L3 increment), consistent with a Rhie-Chow
term that vanishes with h at fixed dt, and grows as dt decreases (as the
`rAU`-weighted pressure smoothing weakens). dt = 0.00625 s is therefore held
fixed at all levels; the FSI dt effect is ~5x the rigid one because, on the
moving mesh, `ddtCorr` uses the face velocity `Uf`, which is consistent with
`phi` at steady state, so the `ddtCorr` contribution vanishes and only the
`rAU`-dependent smoothing remains (interpretation).

### 7. Steady state

L0-L2 histories are in `reference/campaign/production_histories`. From the
converged rigid-beam flow the FSI state is within 3.4e-4 (`u_x`) after 2 s,
1.6e-5 after 4 s and 2.5e-7 per 2 s after 12 s (L2). From rest the wake mode
(tau ~ 2 s) needs ~16-20 s.

## Production sequence (MeluXina, one build), `m = 2`, dt = 0.00625

| level | fluid cells | solid cells (`N_t`) | `u_x(A)` (m) | `u_y(A)` (m) | `u_z(A)` (m) | `F_x` (N) | `F_y` (N) | `F_z` (N) |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| L0 | 8,112 | 512 (8) | 8.13919e-5 | 1.95723e-6 | -2.03080e-6 | 2.15736 | 0.32694 | -0.25742 |
| L1 | 64,896 | 4,096 (16) | 9.62607e-5 | 1.94057e-6 | -2.15035e-6 | 2.30568 | 0.33773 | -0.27647 |
| L2 | 519,168 | 32,768 (32) | 1.024547e-4 | 1.92632e-6 | -2.20300e-6 | 2.36363 | 0.34428 | -0.28799 |
| L3 (mean over t = 4.0-4.5 s) | 4,153,344 | 262,144 (64) | 1.044776e-4 | 1.91580e-6 | -2.21539e-6 | 2.38305 | 0.34577 | -0.29121 |

`u_z` is reported for the computed `z <= 0` half (Richter's A is the mirror
point, `u_z` changes sign); `F` is the force on the half interface.

Relative increments (L0-L1 / L1-L2 / L2-L3) and observed orders (L0-L2,
L1-L3):

| quantity | increments | orders | conditional Richardson value (L1-L3 order) | correction to L3 | p = 2 instead |
| --- | --- | --- | ---: | ---: | ---: |
| `u_x(A)` | 15.4 / 6.0 / 1.9 % | 1.26, 1.61 | 1.0546e-4 m | +0.94 % | 1.0515e-4 m (+0.65 %) |
| `F_x` | 6.4 / 2.5 / 0.82 % | 1.36, 1.58 | 2.3929 N | +0.41 % | 2.3895 N (+0.27 %) |
| `F_y` | 3.2 / 1.9 / 0.43 % | 0.72, 2.14 | - | - | - |
| `F_z` | 6.9 / 4.0 / 1.1 % | 0.73, 1.84 | - | - | - |
| `u_z(A)` | 5.6 / 2.4 / 0.56 % | 1.18, 2.09 | - | - | - |
| `u_y(A)` | -0.86 / -0.74 / -0.55 % | 0.23, 0.44 | - | - | - |

Interpretation (wording agreed with the final Astra review): `u_x` and `F_x`
converge monotonically with effective orders increasing towards second order
(rigid-beam `F_x`: 1.36, 1.55; static solid: 1.34, 1.71). The finest three
meshes suggest an approach to the asymptotic regime, but order stability is
not established, so no GCI is claimed. The Richardson values are conditional
extrapolations; the `p_obs`/`p = 2` pair is an exponent-sensitivity range, not
a bound (a factor-of-safety 3 allowance would be 2.8 % and 1.2 %). `F_y` and
`u_z` converge monotonically with unstable orders; `u_y` converges slowly with
a low effective order. All non-spatial errors (coupling <= 2e-6, history <=
1.6e-6, steadiness, L3 oscillation ~1e-4 decaying) are far below the L2-L3
increments; the fixed-dt Rhie-Chow sensitivity is 2-4e-4 at L2 and falls ~4x
per level, so it is part of the mesh-dependent discretisation trend.

### L3 late-time oscillation

All levels show a ~8 Hz oscillation of the beam after the start, with the
same frequency at L2 (8.0 Hz) and L3 (7.9 Hz), consistent with a damped
fluid-loaded first bending mode (dry cantilever estimate ~15 Hz, lowered by
added mass); it is not vortex shedding (the rigid-beam flows are steady to
1e-7 and the frequency does not change with the mesh). It decays
geometrically, by a factor 2.1 per 0.25 s at L3 (computed damping ratio ~6 %,
physical plus numerical). L3 half-amplitude of `u_x(A)` per 0.25 s window
(linear trend removed): 4.6e-4 (2.50-2.75 s), 9.4e-5 (3.00-3.25), 2.0e-5
(3.50-3.75), 4.2e-6 (4.00-4.25), 2.1e-6 (4.25-4.50); `F_x` 3.5e-4 to 1.6e-6.
The plain restart of segment 2 at t = 3 s (same decomposition, `D_0_0`
written) did not re-excite it; the L2 restart with re-decomposition at
t = 4 s did (to 1.4e-4), after which it decayed again to the same mean. The
L3 values are means over 4.0-4.5 s. Successive 0.25 s window means still
change by +6e-7 (`u_x`), +7e-7 (`F_x`), +1.0e-5 (`F_y`) and -2e-5 (`F_z`)
relative, a slow wake drift as at L2; the estimated remainders (rate x
tau ~ 2 s) are ~2e-6 (`u_x`), ~3e-6 (`F_x`) and ~1e-4 (`F_y`, `F_z`), all
below 5 % of the L2-L3 increments.

## Comparison with Richter (2012), after the sequence above

Richter Table 6 (half interface, symmetric half domain, as here): drag
1.327 +- 0.01 N, `u_x(A)` 5.924e-5 +- 1e-7 m. For the implemented printed
inlet (peak 0.3 m/s; area mean 0.133 m/s, `Re = 0.133 x 0.2 / 1e-3 = 27`
with Richter's length scale) the sequence approaches 2.39 N and
1.05e-4 m (L3: 2.38305 N, 1.04478e-4 m), i.e. about 1.8x the tabulated values; the discrepancy greatly
exceeds the estimated numerical corrections (< 1 %). Earlier solids4foam
calculations with peak 0.2 m/s (coarser solid family, `N_t = 4f`) gave
1.3311 N and 5.915e-5 m at L3, close to Table 6. This suggests a possible
inlet-definition mismatch in the source (the stated mean 0.2 m/s and
`Re = 40` are inconsistent with the printed profile) but does not establish
which inlet produced Table 6. No setting was tuned to Table 6.

## Astra (GPT-6-Astra, xhigh) review, 2026-10-12

Prompt and full reply: `reference/campaign/astra/`. Main points and actions:

- Relative-reduction tests become unattainable at steady state; recommends a
  calibrated absolute floor, `R <= max(1e-5 R_step, R_abs)`. Not available
  in the frozen build; approximated by strict 1e-5 with a counted cap and an
  `r2` audit of capped steps, plus the sensitivity table above. Possible
  follow-up: implement the floor in solids4foam (separate PR).
- Capped 1e-6/1e-7 runs are sensitivity evidence, not references (adopted);
  check PIMPLE 3 vs 5-6 (done: PIMPLE 6, 6e-7) and SNES 1e-9 (done, < 1e-8).
- History mechanism: "supported, not proven"; do not claim reproduction on
  MeluXina (wording adopted; crossover test not run).
- Keep dt fixed, conditional on dt checks at L2 (running) and ideally L3; a
  rigid L4 at another dt cannot enter the extrapolation and does not bound
  the FSI error (rigid L4 not run).
- Keep L3 for credible coupled orders; rigid-flow initialisation at L3 is
  acceptable; stopping budget: transient and iterative errors each < 5 % of
  the finest mesh increment.
- `m = 2` family is systematic; diagnostic solid/fluid errors cannot simply
  be added; monotonicity alone is insufficient for GCI.
- Reported L1 rest/rigid-start difference corrected to 1.6-3.2e-6 relative.

## Astra final review (2026-10-13)

Prompt and reply: `reference/campaign/astra/astra_2_final_*`. Adopted: the
"suggested, not demonstrated" asymptotic wording; conditional Richardson
values without a headline GCI; fixed-dt wording ("the discretisation trend
includes mesh-dependent effects of the time-step-dependent Rhie-Chow
formulation"); window means with stated settling and amplitude; the
qualified Richter incompatibility statement. Open recommendation (decision
for Philip, not run): a coupled L3 dt/2 steady-sensitivity run (estimated
>= 60 node-hours) would close the remaining dt gap at L3.

## Cost (physical cores = node-hours x 128)

73.0 node-hours = 9,350 physical core-h in total since 2026-10-11, of which
L3 44.6 node-hours (segment 1, job 5333947: 36.0; segment 2, job 5336554:
8.6; the failed quota attempt 5333773: 0.1). Segment 2 ran ~1.7x faster per
simulated second than segment 1 (different nodes). Earlier totals in this
README and in reports were sacct CPUTimeRAW, which counts 256 logical CPUs
per node and so doubled the physical figure.

## STATUS (2026-10-13, campaign complete for L0-L3)

All production and check runs are finished; no beam jobs are queued. Open
items (decisions for Philip): a coupled L3 dt/2 steady-sensitivity run (Astra
recommendation, roughly 50-60 node-hours); an absolute-floor coupling
criterion in solids4foam (separate PR); clarification of the Richter inlet
with the author. Run directories: `/project/home/p201287/philip/beam-richter-opus/runs`
(L3 keeps its decomposed t = 4, 4.5 s state); archives in
`/project/home/p201136/philip-beam-archive/beam-richter-opus/`.

<!-- markdownlint-enable MD013 -->
