# bendingBeamInChannel verification

Opt-in verification of the forcing methods of `immersedBoundaryForce` for a
body that deforms while it moves through a stationary mesh:

- `cutLink`;
- `penalty`;
- legacy incremental, `method incremental`;
- `ghostCell` (documented for static bodies; its result here is exploratory).

The case is the beam of the tutorial, which bends with the prescribed
`customProfileBend` motion in a channel at $Re = 20$. The verification
asks three questions and keeps their answers apart. The verification is
not run by `tutorials/Alltest` or `tutorials/Alltest-regression`.

<!-- markdownlint-disable MD013 -->

<!-- markdownlint-disable MD013 -->

<!-- markdownlint-disable MD013 -->

<!-- markdownlint-disable MD013 -->

| | Question | Needs a reference? |
| --- | --- | --- |
| A | How does the net flux through the inlet and the outlet compare with the exact rate of change of the volume of the beam? | no: the analytic rate is exact |
| B | How far are $C_d$ and $C_l$ from the tracked body-fitted series? | yes, of uncertain provenance |
| F | What does a run cost? | no, it is observed |

<!-- markdownlint-enable MD013 -->

<!-- markdownlint-enable MD013 -->

<!-- markdownlint-enable MD013 -->

<!-- markdownlint-enable MD013 -->

Three further groups of results, all from the tracked data in
`reference/results/`, are kept separate from A, B and F in `results.md`:
C, the sensitivity to the finite-volume settings (time step, time scheme,
convection scheme, corrector counts, tolerances); D, the sensitivity to the
formulation (penalty coefficient, the velocity imposed inside the body); and
E, the limits of what the tutorial can show.

## What it shows

The bending changes the volume of the beam inside the channel. With the
prescribed motion the rate is exactly $dV/dt = c \omega \sin(2 \omega t)$, with
$c = (10/7) w L_z A_x A_y / H$ and $\omega = 2 \pi$ rad/s (the rate has twice
the bending frequency), peaking at 0.02831 m$^3$/s, 6.9 % of the
throughflow. For an incompressible flow it is the net flux out of the
channel. It is zero for $A_y = 0$ (Study B). Nothing is needed to know it: the
driver evaluates it in closed form and checks it against the area of the
deformed outline.

The prescribed bending exposes a difference in how the formulations represent
this flux:

- **cutLink** carries the flux. Its net exterior flux is 0.924, 0.950, 0.976
  and 0.988 of the exact peak on Levels 1 to 4; the residual halves with each
  level from Level 2 (4.96 %, 2.45 %, 1.22 % of the peak). It is a
  decay observed on four levels for one diagnostic, not a formal order of
  convergence. It does not depend on the time step, the time scheme or the
  number of ranks.
- **penalty**, legacy **incremental** and **ghostCell** have a net exterior flux
  of zero to the solver tolerance (at most $10^{-8}$ m$^3$/s) in every
  run, whatever the finite-volume settings. This is consistent with a
  continuity treatment that does not represent the flux of a body whose volume
  changes; the tutorial motion is one in which such a body would show it.

The $C_d$ differences between the methods and the tracked body-fitted series
are those of the table below. They are agreement with a numerical reference
whose provenance is not fully documented, not errors: a coarse body-fitted
solution changes by 2.9 % in $C_d$ and 6.2 % in $C_l$ between the time steps
0.01 and 0.0025 s (`reference/results/bodyfitted-R1-control.csv`), which is of
the size of the differences between formulations. The lift of the reference is
unverified and no lift ranking is made.

<!-- markdownlint-disable MD013 -->

<!-- markdownlint-disable MD013 -->

<!-- markdownlint-disable MD013 -->

<!-- markdownlint-disable MD013 -->

| Method (estimator) | $C_d$ nRMS, Level 1 / 2 / 3 | $C_d$ second harmonic, Level 1 / 2 / 3 (reference 5.89) |
| --- | --- | --- |
| cutLink, surface traction | 8.6 / 5.0 / 3.2 % | 5.33 / 5.72 / 5.82 |
| cutLink, momentum exchange | 7.3 / 4.7 / 3.4 % | 5.20 / 5.57 / 5.72 |
| penalty | 10.9 / 7.4 / 6.7 % | 1.81 / 1.84 / 1.91 |
| incremental | 14.0 / 7.6 / 6.8 % | 1.21 / 1.84 / 1.83 |
| ghostCell (exploratory) | 12.0 / 9.8 / 8.1 % | 2.08 / 2.15 / 2.03 |

<!-- markdownlint-enable MD013 -->

<!-- markdownlint-enable MD013 -->

<!-- markdownlint-enable MD013 -->

<!-- markdownlint-enable MD013 -->

The persistent deficit of the second harmonic of penalty, incremental and
ghostCell is
not changed by any of the finite-volume settings tested (a property of the
formulations as implemented, within those settings): it is 1.84 to 2.02 in every
penalty run, 11 of them with two
time steps, two time schemes, three convection schemes, two outer corrector
counts, a doubled pressure-corrector count, a tighter pressure tolerance and
three penalty coefficients over two decades.
That the missing volume flux causes it is not shown. `results.md` has the
evidence, the figures and the limitations.

<!-- markdownlint-disable MD013 -->

<!-- markdownlint-disable MD013 -->

<!-- markdownlint-disable MD013 -->

<!-- markdownlint-disable MD013 -->

| Figure (`images/`) | Content |
| --- | --- |
| `flux_vs_level.png` | net exterior flux against the exact rate, Levels 1-4, Studies A and B |
| `force_histories.png` | $C_d$ and $C_l$ of three methods against the body-fitted series; $C_d$ second harmonic of all four |
| `fv_sensitivity.png` | change of the $C_d$ nRMS and second harmonic with each finite-volume setting |
| `bodyfitted_control.png` | the body-fitted case at two time steps |
| `observed_cost.png` | observed time per step of matched cutLink and penalty runs |
| `mesh_refinement.png` | $C_d$ and $C_l$ nRMS, $C_d$ second harmonic and flux residual of the four methods, Levels 1-3 |
| `histories_L2.png`, `histories_L3.png` | $C_d(t)$ and $C_l(t)$, $1 < t < 3$ s, of the four methods against the body-fitted series |

<!-- markdownlint-enable MD013 -->

<!-- markdownlint-enable MD013 -->

<!-- markdownlint-enable MD013 -->

<!-- markdownlint-enable MD013 -->

## Running

The verification needs OpenFOAM v2412 or later (the immersed boundary library
is built only with OpenFOAM.com), solids4foam, and `python3`. Nothing else is
needed to run the exact tests and the dry run.

```bash
./Allverify unit                      # exact tests, no OpenFOAM needed
./Allverify dry-run --runs L2         # cases, ranks and estimated cost
./Allverify                           # smoke: Level 1, first 0.5 s
./Allverify --runs flux               # Levels 1-2, first 0.5 s: volume balance
./Allverify --runs L1 L2              # force coefficients
./Allverify check                     # PASS / FAIL / SKIPPED
```

| `--runs` | Cases | Estimated allocated rank-hours |
| --- | --- | --- |
| `smoke` (default) | cutLink and penalty, Level 1, 0 to 0.5 s | 0.01 |
| `flux` | cutLink and penalty, Levels 1-2, 0 to 0.5 s | 0.06 |
| `studyB` | cutLink, $A_y = 0$, Levels 1-3, 0 to 0.5 s | 0.4 |
| `L1` | the four methods, Level 1, 0 to 8 s | 0.3 |
| `L2` | cutLink, penalty, incremental, Level 2, 0 to 8 s | 2.0 |
| `L3` | cutLink, penalty, Level 3, 0 to 8 s | 7.7 |
| `L4flux` | cutLink, Level 4 (1.28 M cells), 0 to 0.5 s | 1.7 |

The estimates are scaled from the observed allocated rank-hours of single
runs on a shared node, with the number of ranks given by `--ranks` (default
by level: 1, 2, 4, 4). They are not a performance model. The whole
verification, `--runs L1 L2 L3 studyB L4flux`, is about 12 allocated rank-hours
and is never required.

Each case is a copy of the tutorial under `work/cases/` (not tracked) with
only the method settings, the mesh level, the end time, the outer correctors,
the write interval (0.125 s, a quarter of the period of the volume rate,
where the fluxes are read) and the decomposition changed. The tutorial is
never modified. `work/manifest.json` records every case with the SHA-256 of
its inputs and the load of the machine at the start. Parallel cases are run
with the solids4foam run functions and no CPU binding is requested
(`OMPI_MCA_hwloc_base_binding_policy=none`): by default Open MPI binds the
ranks of every `mpirun` to the first cores, so that two simultaneous runs
share cores. The driver refuses to start more ranks than the cores its
process may use.

## Method settings

Written explicitly into every case, those of the stationary and oscillating
cylinder verifications:

<!-- markdownlint-disable MD013 -->

<!-- markdownlint-disable MD013 -->

<!-- markdownlint-disable MD013 -->

<!-- markdownlint-disable MD013 -->

| Method | Settings | `nOuterCorrectors` |
| --- | --- | --- |
| cutLink | `forceEstimator surfaceTraction` (the momentum exchange is written to the same file) | 1 |
| penalty | `penaltyCoeff 1e3`, `weighting volumeFraction`, `surfaceRateCoeff 3`, `occupancy signedDistance` | 1 |
| incremental | `couplingCoeff 0.2`, `occupancy vertexFraction` | 5 |
| ghostCell | `penaltyCoeff 1e3`, `imageDistance 1.5`, `imageInterpolation cellPoint` | 1 |

<!-- markdownlint-enable MD013 -->

<!-- markdownlint-enable MD013 -->

<!-- markdownlint-enable MD013 -->

<!-- markdownlint-enable MD013 -->

The outer corrector count of incremental is that of the method and is not
equalised with the others. `nCorrectors` is 3 for all. The time step is the
adaptive one of the tutorial, with `maxCo 0.5`. This is the Courant number of
the controller, `fluidModel::CourantNo`, which is face based:
$\max_f (|\phi|/|S_f| \, \delta) \Delta t$ with $\delta$ the cell distance.
The Courant number printed in the log by the fluid model is the cell-sum
definition of OpenFOAM, $0.5 \max (\sum_f |\phi|/V) \Delta t$, which is about
twice as large (up to 1.0 in the tracked runs, 4.1 for the body-fitted case with
$\Delta t$ = 0.01 s). They are two definitions, not an overshoot of the
controller. The driver reports the second one as `max_Co_in_step`.

## Force reconstruction and normalisation

$C_d = 2 F_x / (\rho \bar U^2 H L_z)$ and $C_l = 2 F_y / (\rho \bar U^2 H L_z)$
with $\bar U$, $H$, $L_z$ read from the tutorial (the coefficient is 9.855 times
the kinematic force). The estimators are kept apart:

- cutLink writes the surface traction and the momentum exchange; both are the
  hydrodynamic force and are never averaged;
- penalty, incremental and ghostCell write a forcing and the inertia of the
  fluid inside the body; their sum is the hydrodynamic force, the forcing
  alone is not. The inertia is a first-order difference of the momentum of the
  fluid inside the body, whatever `ddtSchemes` says: a second-order time scheme
  in the fluid equations does not make this force second-order accurate.

The nRMS is the root mean square over time (trapezoid rule) of the difference
from the tracked series over $1 < t < 8$ s, normalised by the root mean square
of the series; the simulated force is interpolated linearly on the reference
times. The harmonics are a direct least-squares fit of a constant and the
first and second harmonics over the same window, without phase averaging or
alignment.

## Check semantics

`./Allverify check` prints one line per check and a summary:

- **PASS**: the check was made and its tolerance met;
- **FAIL**: the check was made and its tolerance not met (exit status 1);
- **SKIPPED**: the check could not be made: no completed case, a run shorter
  than the comparison window, no `phi` field, no tracked result for the case
  (never counted as a pass; `--strict` makes it a failure).

<!-- markdownlint-disable MD013 -->

<!-- markdownlint-disable MD013 -->

<!-- markdownlint-disable MD013 -->

<!-- markdownlint-disable MD013 -->

| Check | Tolerance |
| --- | --- |
| closed-form volume rate against the area of the deformed outline | $10^{-4}$ of the peak |
| cutLink flux ratio at the peaks against the tracked value | 0.01 |
| cutLink flux bias with $A_y = 0$ against the tracked value | 5 % |
| penalty, incremental, ghostCell net exterior flux | $10^{-6}$ of the analytic peak |
| $C_d$ nRMS against the tracked value | 1.0 point |
| $C_d$ second harmonic against the tracked value | 0.1 |

<!-- markdownlint-enable MD013 -->

<!-- markdownlint-enable MD013 -->

<!-- markdownlint-enable MD013 -->

<!-- markdownlint-enable MD013 -->

Both tolerances are the thresholds registered in advance for the sensitivity
study (a change of 1.0 point of $C_d$ nRMS, or of 0.1 in the second harmonic,
is material); they were not adjusted after looking at a result. The check is a
regression check against the tracked values, not a measure of accuracy. The
tracked values come from runs with the adaptive time step of the tutorial,
which is what the generated cases use. The archived fixed-time-step runs of the
campaign are not conforming reruns: for the penalty method they give a second
harmonic 0.06 to 0.18 above the tracked 1.844 (11 runs, for several time
steps, schemes and coefficients), so six of them would fail this tolerance.
The cause of this offset (the time-step control, or run-to-run variation of a
small amplitude) was not isolated, and no complete adaptive-step rerun of a
tracked penalty case exists. A failed A2 check of a penalty case is therefore
to be reported as a failure, not tuned away.

## Cost

`analyse` reports the allocated rank-hours (ranks times wall time), the
number of steps, the cost per cell-step and the ratio of CPU to wall time. They
are observed costs of single runs on a shared machine, not intrinsic costs of
the methods: the number of steps, the rank count, the load of the node and the
CPU binding all change them. In matched runs after the CPU-affinity
correction (below) the wall time per step of cutLink is 1.5 to 2.1 times that
of penalty for the same settings; the cutLink run also evaluates both of its
force estimators. The cost columns of `cases.csv` are `nan` where no wall time
of the tracked
run is available (its eight Study B rows). See `results.md`, Section 6, for what
is and is not
established.

## Provenance and data availability

- **Solver of the tracked results**: solids4foam `f33b244a2` (the
  `development` branch before this verification), OpenFOAM v2412, serial runs at
  Level 1 and mostly at Level 2, two to eight ranks at Levels 2 and 3, SLURM or
  direct MPI; the ghostCell Level 1 run is a repeat with the write interval
  0.125 s of the other runs (the first run wrote every 2 s).
- **Cases**: generated from the tutorial by this driver; Study A is the
  tutorial ($A_y$ = 0.8 m), Study B has $A_y$ = 0.
- **Body-fitted series** `reference/CdBodyFitted.dat`: see the tutorial README.
  Its original case dictionaries, mesh, logs, time-step and mesh sensitivity and
  force-integration set-up have not been recovered; a request to its author is
  outstanding. The comparison with it is a comparison with a numerical
  solution of uncertain accuracy.
- **Body-fitted control** (`reference/results/bodyfitted-R1-*.csv`): the
  archived body-fitted case with 24,720 cells (OpenFOAM v2512, `backward`, two
  outer correctors) run at $\Delta t$ = 0.0025 s and 0.01 s on one installation,
  with only `deltaT`, the write intervals and the decomposition changed.
  The solver is not part of this repository and the two runs are not
  reproducible by `Allverify`. It does not establish the uncertainty of the
  tracked series, which has about four times as many cells.
- **Tracked results** (`reference/results/`): `cases.csv` (one row per
  estimator and coefficient: nRMS, harmonics, noise, fluxes, cost),
  `history-*.csv` (display histories every 0.02 s), `fv-sensitivity.csv`,
  `cost-observed.csv` and `cost-matched-pairs.csv`, `bodyfitted-R1-*.csv`. They
  were computed by the functions of this driver from the complete saved runs;
  fields and logs of those runs are not tracked.

## Limitations

- one tutorial, one body motion, one Reynolds number, a laminar flow;
- the body-fitted series is the only force reference and its accuracy is
  unresolved, in particular that of the lift, which is not ranked;
- mesh and time step were not varied independently except at Level 2: no
  formal order of convergence is claimed for the forces or for the flux;
- the cost is observed, from single runs;
- ghostCell is exploratory (its forces are noisier and carry a decomposition
  scatter of about 1 % of their RMS), and incremental and ghostCell were run
  once per level;
- the surface-traction quadrature options of other pull requests are not
  evaluated here.
