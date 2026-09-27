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
routinely.

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
