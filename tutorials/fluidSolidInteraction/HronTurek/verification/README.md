# HronTurek verification studies

This directory contains opt-in verification studies for the `HronTurek`
tutorial. By default they compare the periodic response of the Turek-Hron FSI3
benchmark, computed with the partitioned Dirichlet-Neumann IQN-ILS and
Robin-Neumann couplings, with the published reference values of Turek and
Hron (2006). The steady FSI1 and the periodic FSI2 benchmarks are described
in their own sections below. The studies are deliberately separate from `regressionTest.sh`:
the regression test checks that the tutorial remains numerically stable,
whereas these studies check convergence towards the benchmark. Nothing here is
run by `tutorials/Alltest` or `tutorials/Alltest-regression`.

Source a supported OpenFOAM environment, build solids4foam with PETSc, and run
from this directory:

```bash
cd tutorials/fluidSolidInteraction/HronTurek/verification
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

The reference values are in `reference/HronTurek_verification_references.json`.
The primary values are the Featflow FSI3 results on level 4 with
`Δt = 0.00025 s`: `u_x = -2.88 ± 2.72 mm [10.93 Hz]`,
`u_y = 1.47 ± 34.99 mm [5.46 Hz]`, `F_D = 460.5 ± 27.74 N/m [10.93 Hz]` and
`F_L = 2.50 ± 153.91 N/m [5.46 Hz]`. This is the discretisation of the
published reference time history, `reference/TurekHron_fsi3_reference_history.csv`
(subsampled to `1 ms` from the Featflow `ref_fsi3.point` file): the driver's
extraction applied to that history reproduces the table, which checks the
extraction itself. The displacement statistics, the drag mean and the lift
amplitude agree to within `0.7%`, and the frequencies to within `0.3%`
(`5.473` against `5.46 Hz` for `u_y`). The drag amplitude is `1.1%` low
(`27.43` against `27.74 N/m`), and the near-zero `u_y` and lift means differ
by `0.02 mm` and `0.04 N/m`. These differences are small against the
tolerances below. The frequently quoted
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
- Every run must be periodic: for `u_x`, `u_y` and the lift, the mean
  amplitude over the last two full `u_y` periods must agree with that over
  the two before to within `2%`. The drag amplitude scatters by about `2%`
  from cycle to cycle without a trend, so its change is reported only.
- On the finest level of the sweep, each primary quantity must be within its
  tolerance of the Featflow level-4 value: `2%` for the mean drag, `3%` for
  the frequencies, `5%` for the `u_y` amplitude, `10%` for the `u_x` mean and
  amplitude, and `15%` for the drag and lift amplitudes. The `u_y` and lift
  means, which are close to zero, are reported as diagnostics only. The lift
  amplitude converges slowest of all the quantities: its error falls from
  `49%` on the 1x mesh to `13%` on the 2x mesh, an observed order of about
  `1.9`, and the `15%` bound applies to the default 2x finest level.
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

IQN-ILS with `predictor yes`, run to `t = 7 s` and evaluated over the closing
`1 s`. Level 1 ran in serial and level 2 on eight ranks. Both have been
periodic since about `t = 4.2 s`.

| Quantity | 1x | 2x | Featflow level 4 | Error at 2x |
|---|---:|---:|---:|---:|
| `u_x` mean (mm) | -2.266 | -2.791 | -2.880 | 3.1% |
| `u_x` amplitude (mm) | 2.198 | 2.662 | 2.720 | 2.1% |
| `u_y` amplitude (mm) | 29.94 | 34.17 | 34.99 | 2.3% |
| `u_y` frequency (Hz) | 5.594 | 5.524 | 5.460 | 1.2% |
| drag mean (N/m) | 457.2 | 459.3 | 460.5 | 0.3% |
| drag amplitude (N/m) | 23.43 | 28.06 | 27.74 | 1.2% |
| lift amplitude (N/m) | 229.2 | 174.4 | 153.9 | 13.3% |

Every primary error decreases from the 1x to the 2x mesh. IQN-ILS needed
about 7 FSI iterations per coupled step on both levels; the runs took about
`1.1 h` (1x, serial) and `4.1 h` (2x, eight ranks).

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

## FSI2 periodic benchmark

`--benchmark fsi2` runs the periodic FSI2 test of Turek and Hron. It has the
same geometry and fluid as FSI3, a mean inflow of `1 m/s` (`Re = 100`) and a
heavy, soft plate with `ρ = 10000 kg/m^3`, `E = 1.4 MPa` and `ν = 0.4`. The
plate flaps with a large amplitude: `u_y(A)` is about `±80 mm` at `1.93 Hz`.
The FSI3 machinery is reused unchanged: the same periodic statistics,
periodicity test, history overlay and force integration.

```bash
./Allverify --benchmark fsi2                   # IQN-ILS mesh study, 1x and 2x
./Allverify --benchmark fsi2 --levels 1,2,4    # add the 4x mesh (expensive)
./Allverify --benchmark fsi2 --quick           # smoke run to t = 2.3 s, no checks
```

The FSI2 settings are in the `fsi2` entry of the reference JSON file. Runs,
CSV files and plots carry an `fsi2_` prefix. In addition to the FSI3 changes
listed above, the verification copies set:

- the inlet `maxValue` to `1.5 m/s`, ramped over the benchmark's
  `transitionPeriod 2 s`, in both `U.dirichletNeumann` and `U.robin`;
- `ρ = 10000 kg/m^3` and `E = 1.4e6 Pa`;
- the coupling start at `t = 2 s`, the end of the ramp;
- `outerCorrTolerance 1e-6`, the tutorial value, instead of the `1e-5` of the
  FSI3 copies. At `1e-5`, the partly converged interface leaves step-to-step
  noise in the forces: on the 1x mesh the lift noise was `14 N/m` rms, with
  peaks of `40 N/m`, and on the 2x mesh `30 N/m` rms, with peaks of
  `140 N/m`. Because the amplitude is taken from the extrema, this noise
  inflated the lift amplitude by `13%` on the 1x mesh and `30%` on the 2x
  mesh. At `1e-6` the noise falls by a factor of about six. IQN-ILS needs one
  more iteration per step, and the lift amplitude is within `3%` of that of
  the smoothed signal. At `1e-7` there was no further change. The stall that
  motivates `1e-5` for FSI3 is discussed below; it is not avoided by `1e-5`.

The time steps are those of FSI3, `Δt = 0.001 s` on the 1x mesh, `0.0005 s`
on the 2x mesh and `0.00025 s` on the 4x mesh. The maximum Courant number is
`0.40` on 1x and `0.41` on 2x, half that of FSI3 because the inflow is halved;
the moving plate adds little. A period of about `0.52 s` is resolved by
about 500 steps on the 1x mesh. The Featflow tables change by at most `1%`
between `Δt = 0.01 s` and `0.0005 s` on level 4, so the time-step error is
small against the mesh error.

### FSI2 end time and window

The run ends at `t = 10.5 s` and is evaluated over the closing `2.6 s`, from
`t = 7.9 s` on the 2x mesh. That is four full `u_y` periods, which is the
minimum for the two-against-two periodicity test. Both limits come from the
runs:

- The flow before the coupling starts is steady. After the coupling starts,
  the flapping grows from rest by about a factor of 1.6 every `0.5 s`, and
  the `u_y` amplitude saturates at `t ≈ 7.4 s` (2x) and `7.7 s` (1x).
- After saturation the amplitude drifts slowly downwards, by about `0.1%` per
  period on the 2x mesh and `0.2%` on the 1x mesh, while the fluid mesh
  distorts. The largest non-orthogonality rises steadily during the flapping,
  on the 2x mesh from `26°` at rest to `42°` at `t = 9 s`, `58°` at
  `10.5 s` and `61°` at `11 s`. On the 1x mesh it reaches `79°` at
  `13.3 s`. The `velocityLaplacian` motion solver is not reversible over a
  period, so the interior points drift. Eventually IQN-ILS stops converging
  within 30 iterations. With `outerCorrTolerance 1e-6`, this happened at
  `t = 11.36 s` on the 2x mesh and at `18.67 s` on the 1x mesh. With `1e-5`,
  it happened at `11.35 s` (2x) and `19.65 s` (1x), and on the 2x mesh the
  solid SNES then failed. The loosened tolerance therefore does not avoid it.
  The run is stopped at `10.5 s`, before the mesh degrades further.

The mesh drift is a limitation of the tutorial's mesh motion on this case. It
is not addressed here.

### FSI2 reference and acceptance

The reference values are the Featflow FSI2 results on level 4 with
`Δt = 0.0005 s`: `u_x = -14.85 ± 12.70 mm [3.86 Hz]`,
`u_y = 1.30 ± 81.6 mm [1.93 Hz]`, `F_D = 215.06 ± 77.65 N/m [3.86 Hz]` and
`F_L = 0.61 ± 237.8 N/m [1.93 Hz]`. All three published FSI2 tables (levels
2 to 4 at `Δt = 0.02`, `0.01` and `0.0005 s`) are recorded in the JSON file
(`featflowTables`).

This table is the discretisation of the published reference history,
`reference/TurekHron_fsi2_reference_history.csv`. It is subsampled to `1 ms`
from the Featflow `ref_fsi2.point` file. Applied to that history, the
driver's extraction reproduces all the amplitudes and the `u_x` and drag
means to within `0.3%`, and the frequencies to the tabulated digits. The
exceptions are the two near-zero means. The `u_y` mean comes out as `1.25`
rather than `1.30 mm`, and the lift mean as `0.27` rather than `0.61 N/m`.
Both differences are below `0.15%` of the corresponding amplitude. The lift
mean is `0.75 N/m` at the full `0.5 ms` sampling. The summary values of Turek
and Hron (2006), `u_x = -14.58 ± 12.44 mm [3.8 Hz]`,
`u_y = 1.23 ± 80.6 mm [2.0 Hz]`, `F_D = 208.83 ± 73.75 N/m` and
`F_L = 0.88 ± 234.2 N/m`, are reported alongside for information.

From level 3 to level 4 the Featflow reference changes by `2.5%` in the `u_x`
mean, `2.0%` in the `u_x` amplitude, `1.1%` in the `u_y` amplitude, `0.9%` in
the drag mean, `2.5%` in the drag amplitude and `1.3%` in the lift
amplitude. The frequencies do not change. The level-4 values are therefore
themselves uncertain by about `1%`, which sets the floor of the tolerances
below.

- Every run must complete to the end time.
- Every run must be periodic: for all four quantities, including the drag,
  the mean amplitude over the last two full `u_y` periods must agree with that
  over the two before to within `2%`. The largest change was `1.1%` (`u_x` on
  the 1x mesh).
- On the finest level of the sweep, each primary quantity must be within its
  tolerance of the Featflow level-4 value: `1.5%` for the drag mean; `3%` for
  the `u_y` amplitude and the `u_x` and `u_y` frequencies; `5%` for the `u_x`
  mean and amplitude and the drag amplitude; and `12%` for the lift
  amplitude. These are the recorded 2x errors, rounded up by about 1.5 times
  and by at least the reference's own spread. The `u_y` and lift means, and
  the drag and lift frequencies (which equal the `u_x` and `u_y`
  frequencies), are reported only.
- The reference error of a primary quantity may not grow by more than one
  percentage point between the coarsest and the finest level of the sweep.
  The lift amplitude is exempt (`errorGrowthChecked false`), as it does not
  converge monotonically: it is `8%` high on both meshes.

The lift amplitude converges slowest, as it does in FSI3. Its extrema-based
value is `256.7 N/m` on the 1x mesh and `256.3 N/m` on the 2x mesh, against
`237.8 N/m`. On the 2x mesh, the last period, smoothed over 21 steps, gives
`250.7 N/m` (`5.4%` high). The remaining `2.4%` is residual coupling noise.
The history overlay shows the lift peaks sharper than the reference and the
plateaus between them slightly lower. The frequencies converge next slowest:
`7%` high on the 1x mesh and `1.8%` high on the 2x mesh.

### FSI2 recorded results

Recorded with OpenFOAM v2412 on an Apple M1 Ultra shared with other jobs.
The 1x mesh ran in serial in `4173 s` (`1.2 h`). The 2x mesh ran on seven
ranks in `14346 s` (`4.0 h`), because one core was in use by the 1x run.
IQN-ILS needed a mean of `4.6` FSI iterations per coupled step on both meshes
(at most 10 on 1x and 12 on 2x). The errors are relative to Featflow level 4,
`Δt = 0.0005 s`:

| Quantity | 1x | Error | 2x | Error | Featflow level 4 |
|---|---:|---:|---:|---:|---:|
| `u_x` mean (mm) | -11.643 | 21.6% | -14.396 | 3.1% | -14.85 |
| `u_x` amplitude (mm) | 10.618 | 16.4% | 12.458 | 1.9% | 12.70 |
| `u_x` frequency (Hz) | 4.130 | 7.0% | 3.929 | 1.8% | 3.86 |
| `u_y` mean (mm) | 1.138 | 12.5% | 1.253 | 3.6% | 1.30 |
| `u_y` amplitude (mm) | 71.51 | 12.4% | 80.55 | 1.3% | 81.6 |
| `u_y` frequency (Hz) | 2.063 | 6.9% | 1.965 | 1.8% | 1.93 |
| drag mean (N/m) | 205.19 | 4.6% | 216.04 | 0.5% | 215.06 |
| drag amplitude (N/m) | 67.45 | 13.1% | 79.24 | 2.0% | 77.65 |
| lift mean (N/m) | 0.90 | - | -0.88 | - | 0.61 |
| lift amplitude (N/m) | 256.70 | 7.9% | 256.32 | 7.8% | 237.8 |

Every primary error except that of the lift amplitude falls by a factor of
four to ten from the 1x to the 2x mesh. The largest periodicity change on the
2x mesh was `0.9%` (drag). `reference/fsi2_iqnils_mesh_2x_history.png`
overlays the 2x history on the published one.

The 4x level (`Δt = 0.00025 s`, 16 ranks by default) is reachable with
`--levels 1,2,4`, but it has not been run. Its end time would first have to
be checked against the mesh drift described above: the finer mesh
distorted sooner.

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
