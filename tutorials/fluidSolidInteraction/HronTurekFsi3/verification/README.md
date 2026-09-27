# HronTurekFsi3 verification studies

This directory contains opt-in verification studies for the `HronTurekFsi3`
tutorial. They compare the periodic response of the Turek-Hron FSI3
benchmark, computed with the partitioned Dirichlet-Neumann IQN-ILS and
Robin-Neumann couplings, with the published reference values of Turek and
Hron (2006). The studies are deliberately separate from `regressionTest.sh`:
the regression test checks that the tutorial remains numerically stable,
whereas these studies check convergence towards the benchmark. Nothing here is
run by `tutorials/Alltest` or `tutorials/Alltest-regression`.

Source a supported OpenFOAM environment, build solids4foam with PETSc, and run
from this directory:

```bash
cd tutorials/fluidSolidInteraction/HronTurekFsi3/verification
./Allverify                            # IQN-ILS mesh study, levels 1x and 2x
./Allverify --coupling robin           # the same sweep with Robin-Neumann coupling
./Allverify --levels 1,2,4             # add the 4x mesh (expensive)
./Allverify --study coupling           # Robin vs IQN-ILS, 1x and 2x meshes
./Allverify --quick                    # smoke run to t = 2.3 s, no checks
./Allverify --reuse                    # re-evaluate completed runs
```

The driver requires `python3`, `blockMesh` and `solids4Foam`; `gnuplot` is
optional and is used for the history plots. Each run is a complete copy of the
tutorial under `verification/work/`, so the tutorial itself and its regression
test are not modified. Results are written to `verification/postProcessing/`
as CSV files, `verification_summary.md`, and PNG history plots. Both
directories are ignored by Git and are retained to make a failed run
diagnosable.

## Benchmark quantities

The benchmark reports the displacement of the plate tip point A at
`(0.6, 0.2)` and the drag and lift on the cylinder and plate together, each as
`mean ± amplitude [frequency]`, where the mean is `(max + min)/2` and the
amplitude is `(max - min)/2` over the last full period of the plate motion.
The driver evaluates the same statistics over the closing `1 s` of each run
(about five `u_y` periods). Periods are delimited by upward crossings of the
window mean with a hysteresis band, so that the harmonics in the force signals
are not counted as periods, and the extrema of every quantity are taken over
the last full `u_y` period: `u_x` and the drag oscillate at twice the plate
frequency with alternating troughs, so their own, shorter period would miss
the deeper trough. The frequency of each quantity is the average over its own
full periods in the window.

The verification copies differ from the tutorial in four respects:

- the forces are integrated over both the `cylinder` and `plate` patches with
  `rhoInf 1000` and divided by the `0.015 m` mesh thickness, so that they are
  per unit depth like the published values;
- the plate uses `StVenantKirchhoffElastic`, the constitutive law specified by
  the benchmark, instead of the tutorial's `neoHookeanElastic`;
- the interface tolerance `outerCorrTolerance` is `1e-5` rather than the
  tutorial's `1e-6`, at which the IQN-ILS residual occasionally stalls just
  above the tolerance and a long run aborts (about once in five thousand
  steps);
- the run continues to `t = 7 s` (the tutorial stops at `6 s`) so that the
  closing window is well inside the periodic regime. The coupling is
  activated at `t = 2 s`, as in the tutorial.

The reference values are in `reference/HronTurekFsi3_verification_references.json`.
The primary values are the Featflow FSI3 results on level 4 with
`Δt = 0.00025 s`: `u_x = -2.88 ± 2.72 mm [10.93 Hz]`,
`u_y = 1.47 ± 34.99 mm [5.46 Hz]`, `F_D = 460.5 ± 27.74 N/m [10.93 Hz]` and
`F_L = 2.50 ± 153.91 N/m [5.46 Hz]`. This is the discretisation of the
published reference time history, `reference/TurekHron_fsi3_reference_history.csv`
(subsampled to `1 ms` from the Featflow `ref_fsi3.point` file): the driver's
extraction applied to that history reproduces the table to within the
tabulated digits, which checks the extraction itself. The frequently quoted
summary values of Turek and Hron (2006), `u_x = -2.69 ± 2.53 mm [10.9 Hz]`,
`u_y = 1.48 ± 34.38 mm [5.3 Hz]`, `F_D = 457.3 ± 22.66 N/m` and
`F_L = 2.22 ± 149.78 N/m`, and the solids4foam values of Tuković et al.
(2018) are reported alongside for information. The driver overlays the
published history on the closing window of each run, with the phases aligned
at the last `u_y` maximum.

## Mesh levels and time steps

Level 1 is the mesh shipped with the tutorial. Each further level doubles the
in-plane block divisions of both the fluid and the solid mesh, leaving the
single spanwise cell alone, and halves the time step so that the Courant number
is unchanged:

| Level | Refinement | Fluid cells | Solid cells | Δt (s) | Default cores |
|---:|---:|---:|---:|---:|---:|
| 1 | 1x | 5 336 | 630 | 0.001 | 1 |
| 2 | 2x | 21 344 | 2 520 | 0.0005 | 8 |
| 3 | 4x | 85 376 | 10 080 | 0.00025 | 16 |

Levels 1 and 2 form the default sweep; level 3 is reachable with
`--levels 1,2,4` and is expensive. Use `--cores N` to run every level on the
same number of MPI ranks (`--cores 1` runs in serial). Level 1 is serial by
default because the small meshes and the many short solid solves make a
four-rank run slower than a serial one. Fields are written only at the end
time; the point displacement and force histories are written every time step
regardless.

## Acceptance criteria

- Every run must complete to the end time.
- Every run must be periodic: the amplitude of each monitored quantity over
  the last two full `u_y` periods must agree to within `2%`.
- On the finest level of the sweep, each primary quantity must be within its
  tolerance of the Featflow level-4 value: `2%` for the mean drag, `3%` for
  the frequencies, `5%` for the `u_y` amplitude, `10%` for the `u_x` mean and
  amplitude and the lift amplitude, and `15%` for the drag amplitude. The
  `u_y` and lift means, which are close to zero, are reported as diagnostics
  only.
- The reference error of a primary quantity may not grow by more than one
  percentage point between the coarsest and the finest level of the sweep.

## Coupling study

The `coupling` study compares the IQN-ILS and Robin-Neumann transients after
the coupling start, sample by sample, on the 1x and 2x meshes up to
`t = 2.05 s`. Both variants restart from the state of a single uncoupled
Dirichlet-Neumann run to `t = 2 s`: the Robin variant's plate conditions,
`elasticWallPressure` and `elasticWallVelocity`, are active during the
uncoupled phase too and do not act as a zero-gradient wall while the plate is
at rest, so separate runs from `t = 0` enter the coupling from different flow
states (`2.8 N/m` apart in lift at `t = 2 s`), after which the transients
cannot be compared.

From the common start, both couplings converge every step, yet their histories
differ by more than the interface tolerance. The difference is not a time
discretisation error (it is unchanged when `Δt` is halved), nor does it come
from `elasticWallVelocity` (which gives results identical to
`newMovingWallVelocity` with a zero-gradient pressure) or from the
zero-gradient pressure itself (`fixedFluxPressure` gives identical results).
It comes from the converged Robin pressure condition, which prescribes the
interface pressure gradient from the interpolated solid acceleration, whereas
the Dirichlet-Neumann wall is kinematically exact in the discrete sense; the
two formulations are different spatial discretisations of the same interface
condition. The study therefore checks, for each of `u_x`, `u_y`, drag and
lift, with the maximum difference normalised by the maximum of the IQN-ILS
history:

- that the Robin vs IQN-ILS difference decreases from the 1x to the 2x mesh;
- that on the 2x mesh it is smaller than the change in the IQN-ILS history
  between the 1x and 2x meshes, i.e. smaller than the discretisation error of
  either solution.

It also reads the recorded residual history to ensure that every coupled
Robin time step terminates with the displacement, pressure-change and
leakage-flux residuals below their configured tolerances, and it reports the
number of FSI iterations of each coupling.

The Robin-Neumann coupling is expensive on this case, which is why the study
stops at `t = 2.05 s`. The plate is thin, wetted on both sides and as dense as
the fluid, and its response involves several bending modes with different
interface impedances, which a single Robin coefficient cannot match: on the
1x mesh the automatic `secant` coefficient settles at about twice its
`thicknessLimited` seed with a predicted contraction factor of `0.92` per
iteration, and every coupled step needs `150` to `180` unrelaxed fixed-point
iterations; on the 2x mesh about `52`. Larger coefficients (`hsModel
pWaveSpeed`, or a constant `hs` of `0.03 m` or more) diverge at the coupling
start, and smaller ones converge even more slowly. The `robin` mesh study is
therefore possible (`--coupling robin`) but is not expected to be run
routinely. The cost and the uncoupled-phase behaviour are the subject of
[issue #493](https://github.com/solids4foam/solids4foam/issues/493).

## Reference results

Recorded with OpenFOAM v2412 on an Apple M1 Ultra shared with other jobs.

### Coupling comparison

| Quantity | Robin vs IQN-ILS, 1x | Robin vs IQN-ILS, 2x | IQN-ILS 1x to 2x |
|---|---:|---:|---:|
| `u_x` | 5.08% | 1.92% | 8.51% |
| `u_y` | 1.72% | 1.22% | 26.5% |
| drag | 0.24% | 0.13% | 1.55% |
| lift | 3.53% | 2.76% | 37.0% |

| Mesh | Coupling | Mean FSI iterations per step | Maximum |
|---|---|---:|---:|
| 1x | IQN-ILS | 9.0 | 11 |
| 1x | Robin-Neumann | 151.5 | 180 |
| 2x | IQN-ILS | 9.0 | 13 |
| 2x | Robin-Neumann | 51.7 | 66 |

The worst converged Robin pressure-change and leakage-flux residuals were
`1.0e-5` and `4.7e-7`.

### Mesh study

To be recorded from the first complete sweep.

## FSI1 steady benchmark

`--benchmark fsi1` runs the steady FSI1 test of Turek and Hron instead of
FSI3. It has the same geometry and fluid, a mean inflow of `0.2 m/s`
(`Re = 20`) and a plate with `E = 1.4 MPa`, `ν = 0.4` and
`ρ = 1000 kg/m^3`, and the flow and plate settle to a steady state. The
steady point-A displacement, drag and lift are compared with the Featflow FSI1
table:

```bash
./Allverify --benchmark fsi1                   # IQN-ILS mesh study, 1x and 2x
./Allverify --benchmark fsi1 --levels 1,2,4    # add the 4x mesh (expensive)
./Allverify --benchmark fsi1 --study coupling  # Robin vs IQN-ILS, 1x mesh
./Allverify --benchmark fsi1 --quick           # smoke run, no checks
```

The default `--benchmark fsi3` behaviour is unchanged. The FSI1 settings are
in the `fsi1` entry of the reference JSON file. In addition to the FSI3 changes
listed above, the verification copies set:

- the inlet `maxValue` to `0.3 m/s`, ramped over `transitionPeriod 2 s`, in
  both `U.dirichletNeumann` and `U.robin`;
- `E = 1.4e6 Pa`;
- the coupling start at `t = 2 s`, the end of the ramp;
- fluid solver tolerances of `1e-9` instead of `1e-6`, which is needed because
  the FSI1 loads and displacements are 30 to 40 times smaller than those of
  FSI3. With the tutorial tolerances, the inner-solver residuals put a floor
  of about `1e-4` on the relative interface residual, and IQN-ILS stops
  converging on the 4x mesh. On the 1x and 2x meshes, which converge with
  either tolerance, the steady values agree to `2e-4`;
- `restart yes` and `writePrecision 12`, so that the coupling study can
  restart from the steady state.

The run is a pseudo-transient route to the steady state with
`Δt = 0.025 s` (`0.0125 s` on the 4x mesh) to `t = 30 s`. The time step is
25 times that of FSI3 and is limited by the IQN-ILS coupling rather than
by accuracy. These limits were found with the tutorial solver tolerances. At
`Δt = 0.05 s` and `0.1 s` on the 1x mesh, with coupling from `t = 0`, IQN-ILS
stalls near the interface tolerance at `t = 2.3 s` and `0.8 s`. On the 4x
mesh, `Δt = 0.025 s` (maximum Courant number 6) stalls in the first coupled
steps. Coupling from `t = 0` or `t = 0.5 s`, while the inflow is still small,
stalls on the 2x mesh, so the coupling starts at the end of the ramp. The
steady values depend slightly on the time step through the PIMPLE flux
interpolation: going from `Δt = 0.025 s` to `0.01 s` on the 1x mesh changes
`u_y` by `0.04%`, drag by `0.04%` and `u_x` by `0.26%`, all far below the
discretisation error.

### FSI1 reference and acceptance

The reference values are the finest level (7+0, about one million elements) of
the Featflow FSI1 table: `u_x(A) = 2.270493e-5 m`, `u_y(A) = 8.208773e-4 m`,
drag `14.29426 N/m` and lift `0.7637460 N/m`. All six tabulated levels are
recorded in the JSON file (`featflowLevels`). Between the coarsest level 2+0
and level 7+0, the reference changes by `0.73%` in `u_x`, `0.19%` in `u_y`,
`0.14%` in drag and `0.26%` in lift. Its own discretisation error is therefore
negligible against the tolerances below.

- Every level must reach a steady state: the relative spread of each quantity
  over the closing `20%` of the run must be below `1e-3`. A level that has not
  settled fails instead of being compared.
- `u_y(A)` and drag are the primary quantities. On the finest level of the
  sweep, drag must be within `0.5%` of the reference, and `u_y` within `6%`
  on the 2x mesh or `3%` on the 4x mesh. The tolerances follow from the
  recorded errors below. `u_y` converges at an observed order of about 1.6,
  so its error falls by a factor of about three per level. The 4x tolerance
  is an extrapolation, because the 4x level has not yet been run.
- `u_x(A)` (`23 μm`, the axial stretch of the plate) and the lift (`0.76 N/m`,
  set by the small asymmetry of the cylinder position) are reported with a
  `3%` indicative tolerance but do not fail the study. They are small
  resultants and are the most sensitive to the time step and to the remaining
  steady-state drift. `u_x`, for example, changes by `0.26%` between
  `Δt = 0.025 s` and `0.01 s`, and the lift settles last.
- The reference error of a primary quantity may not grow by more than `0.1`
  percentage points between the coarsest and the finest level.

### FSI1 coupling study

The opt-in coupling study compares the steady states of the two couplings on
the 1x mesh. It is informative only: the Robin-Neumann coupling on the
Turek-Hron cases is studied separately in
[issue #493](https://github.com/solids4foam/solids4foam/issues/493),
and the FSI1 acceptance rests on the IQN-ILS mesh study. Both couplings
restart from the steady IQN-ILS state of the mesh-study run at `t = 30 s`
and continue for `1 s`. IQN-ILS keeps its
Dirichlet conditions, and Robin-Neumann has the plate switched to
`elasticWallPressure` and `elasticWallVelocity`. The driver reports the
difference between the means of each quantity over the closing half of the
continuation. It fails only if the start state has not settled or a Robin step
does not meet its convergence criteria.

A restart replaces separate runs from rest for two reasons. First, the
Robin-Neumann fixed-point iterations are very expensive near rest: with
coupling from `t = 0`, the first two steps need 283 and 585 iterations. With
coupling from `t = 2 s`, the first coupled step needs 675 iterations. Second,
before the coupling starts, the Robin plate conditions do not act as a
zero-gradient wall, so the two couplings would enter the coupled phase from
different flow states.

The continuation uses `Δt = 0.001 s` (the FSI3 time step) for both couplings.
At `Δt = 0.025 s`, the Robin iterations contract by only about `0.99` per
iteration and do not reach the interface tolerance within 300 iterations, and
a constant `hs` of `0.1 m` diverges within a few iterations. The change of
time step moves the state by up to `0.2%` in drag and `1%` in `u_x`,
identically for both couplings, so the comparison is between two converged
couplings on the same transient.

### FSI1 recorded results

Recorded with OpenFOAM v2412 on an Apple M1 Ultra. The 1x mesh ran in serial
in `781 s`. The 2x mesh ran on two ranks in `3062 s`. IQN-ILS needed a mean of
`3.7` FSI iterations per coupled step on the 1x mesh (at most 8) and `4.1` on
the 2x mesh (at most 9). The largest steady spread over the closing `20%` was
`2.0e-4` on the 1x mesh and `3.1e-4` on the 2x mesh, both in the lift. The
errors are relative to Featflow level 7+0:

| Quantity | 1x | Error | 2x | Error | Observed order | Featflow 7+0 |
|---|---:|---:|---:|---:|---:|---:|
| `u_x(A)` (m) | 2.20672e-5 | 2.81% | 2.24843e-5 | 0.97% | 1.5 | 2.270493e-5 |
| `u_y(A)` (m) | 7.04258e-4 | 14.21% | 7.82820e-4 | 4.64% | 1.6 | 8.208773e-4 |
| drag (N/m) | 14.18027 | 0.80% | 14.26063 | 0.24% | 1.7 | 14.29426 |
| lift (N/m) | 0.803445 | 5.20% | 0.775620 | 1.55% | 1.7 | 0.7637460 |

The 1x mesh has 5 336 fluid and 630 solid cells, and the 2x mesh has 21 344
and 2 520. They are comparable in resolution to Featflow levels 2+0 and 3+0
(992 and 3 968 biquadratic elements), which are within `0.2%` of level 7+0
in `u_y` and drag. The second-order finite-volume discretisation needs
considerably finer meshes for the same plate deflection. All four errors fall
monotonically. The observed orders use the finest Featflow value as the exact
solution.

The 4x level (`Δt = 0.0125 s`) is reachable with `--levels 1,2,4`. On three
ranks it took about `46 s` per coupled step, which is more than 14 hours to
`t = 30 s`, so it was not run for this record.

In the steady coupling study, both couplings restarted from the
`t = 30 s` IQN-ILS state on the 1x mesh and ran to `t = 31 s` at
`Δt = 0.001 s`. The table gives the means over `t = 30.5` to `31 s`. The
differences and the step-to-step noise (the range over the same window) are
relative to the steady value:

| Quantity | IQN-ILS | Robin-Neumann | Difference | IQN-ILS noise | Robin noise |
|---|---:|---:|---:|---:|---:|
| `u_x(A)` (m) | 2.233184e-5 | 2.233136e-5 | 0.0022% | 0.54% | 0.074% |
| `u_y(A)` (m) | 7.045302e-4 | 7.045790e-4 | 0.0069% | 0.042% | 0.034% |
| drag (N/m) | 14.20647 | 14.20644 | 0.0002% | 0.088% | 0.006% |
| lift (N/m) | 0.8035243 | 0.8034405 | 0.0104% | 16.5% | 0.34% |

The Robin-Neumann coupling therefore converges to the same steady discrete
solution as IQN-ILS, to within `0.01%` in all four quantities. This is
consistent with the FSI3 coupling study, which attributes the transient
difference to the Robin pressure gradient built from the solid acceleration:
that term vanishes at a steady state. The first continued step differs by
`0.12%` in `u_x`, `0.22%` in drag and `2.8%` in lift, which is the size of
the IQN-ILS step-to-step noise.
The lift noise of IQN-ILS alternates from step to step at `Δt = 0.001 s` and
is the reason the comparison uses window means.

IQN-ILS needed a mean of `5.9` iterations per step (at most 8). Robin-Neumann
needed `17.6` (at most 56), with the automatic `secant` coefficient. The two
continuations took `1053 s` (IQN-ILS) and `1276 s` (Robin-Neumann) in serial.
Of the 1000 Robin steps, 976 met all three Robin criteria. The other 24 ended
on a stalled pressure residual, which the coupling accepts
(`robinConvergenceState 2`). Over all steps, the interface-displacement
residual was at most `4.2e-7` and the leakage-flux residual at most
`4.4e-7`, but in those 24 steps the pressure-change residual was up to
`9.8e-5`, against a `robinPressureTolerance` of `1e-5`.

## References

S. Turek and J. Hron, Proposal for numerical benchmarking of fluid-structure
interaction between an elastic object and laminar incompressible flow. In:
H.-J. Bungartz and M. Schäfer (eds), *Fluid-Structure Interaction*, Lecture
Notes in Computational Science and Engineering 53, Springer, 2006, 371-385.
Reference data: <https://wwwold.mathematik.tu-dortmund.de/~featflow/en/benchmarks/cfdbenchmarking/fsi_benchmark.html>.

Ž. Tuković, A. Karač, P. Cardiff, H. Jasak, A. Ivanković, OpenFOAM finite
volume solver for fluid-solid interaction, *Transactions of FAMENA*, 42(3),
1-31, 2018.

Ž. Tuković, M. Bukač, P. Cardiff, H. Jasak, A. Ivanković, Added mass
partitioned fluid-structure interaction solver based on a Robin boundary
condition for pressure. In: *OpenFOAM: Selected Papers of the 11th Workshop*,
Springer, 2019, 1-22.
