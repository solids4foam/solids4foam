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
use p/U 1e-10 and mesh motion 1e-12 (MeluXina confirmation at L1:
`P_L1_init_loose` vs `P_L1_init`, running).

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

## Production design (pending Astra review and L1 matrix)

- Mesh family: graded fluid family, every count x2 per level, `m = 2`.
- Coupling: strict `max(r1, r2) <= 1e-6`, `nOuterCorr 50`,
  `allowUnconvergedCoupling yes` only as a counted safety net, IQN-ILS as in
  the tutorial; PIMPLE `nOuterCorrectors 3`; p/U tolerance 1e-10, mesh
  motion 1e-12, SNES rtol/stol 1e-6 (1e-9 checked at L1).
- Steady state: L0-L2 from rest and from the converged rigid flow (history
  check); L3 from the converged rigid flow.

## STATUS (2026-10-12)

MeluXina jobs (account p201287, `/project/home/p201287/philip/beam-richter-opus`):

- 5329141 `beamO_c1`: L1 coupling-tolerance matrix from rest to t = 16 s
  (`runs/c1_L1_tol{1e-4,1e-5,1e-6,1e-7}`, `_loose`, `_dt0125`).
- 5329220 `beamO_p0`: `P_L0_rest`, `P_L0_init`, `P_L1_init`,
  `P_L1_rest_snes1e-9`, `P_L1_init_loose`, `P_L2_init` (L2 needs a
  continuation beyond the 14 h walltime).
- 5329151 `beamO_solid`: solid diagnostic on MeluXina (m = 1, m = 2,
  high order; m = 1 at 16x).

Next: Astra review of the coupling criterion and production plan; L2 from
rest; L3 from the rigid L3 flow (estimate 2700-4500 core-h, to be confirmed
with Philip if above 3000).

<!-- markdownlint-enable MD013 -->
