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
./Allverify --study coupling           # Robin vs IQN-ILS transient on the 1x mesh
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
amplitude is `(max - min)/2` over the last full period. The driver evaluates
the same statistics over the closing `1 s` of each run (about five `u_y` and
ten `u_x` periods), delimiting periods by upward crossings of the window mean
with a hysteresis band so that the harmonics in the force signals are not
counted as periods. The frequency is the average over the full periods in the
window.

The verification copies differ from the tutorial in three respects, each of
which brings the copy closer to the benchmark definition:

- the forces are integrated over both the `cylinder` and `plate` patches with
  `rhoInf 1000` and divided by the `0.015 m` mesh thickness, so that they are
  per unit depth like the published values;
- the plate uses `StVenantKirchhoffElastic`, the constitutive law specified by
  the benchmark, instead of the tutorial's `neoHookeanElastic`;
- the run continues to `t = 7 s` (the tutorial stops at `6 s`) so that the
  closing window is well inside the periodic regime. The coupling is
  activated at `t = 2 s`, as in the tutorial.

The reference values are in `reference/HronTurekFsi3_verification_references.json`.
The primary values are the Turek-Hron level-4 results: `u_x = -2.69 ± 2.53 mm
[10.9 Hz]`, `u_y = 1.48 ± 34.38 mm [5.3 Hz]`, `F_D = 457.3 ± 22.66 N/m` and
`F_L = 2.22 ± 149.78 N/m`. The solids4foam values of Tuković et al. (2018) are
recorded alongside for information. `reference/TurekHron_fsi3_reference_history.csv`
holds the published level-4 time history (subsampled to `1 ms` from the
Featflow `ref_fsi3.point` file), which the driver overlays on the closing
window of each run with the phases aligned at the last `u_y` maximum, and from
which it also evaluates the same statistics as a consistency check on the
extraction.

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
  the last two full periods must agree to within `2%`.
- On the finest level of the sweep, each primary quantity must be within its
  tolerance of the Turek-Hron value: `10%` for the `u_x` mean and amplitude,
  `5%` for the `u_y` amplitude, the lift amplitude and the frequencies, and
  `3%` for the mean drag. The `u_y` and lift means, which are close to zero,
  and the drag amplitude, which varies by more than `30%` between the
  benchmark's own mesh levels and time steps, are reported as diagnostics only.
- The reference error of a primary quantity may not grow by more than one
  percentage point between the coarsest and the finest level of the sweep.

The frequency tolerance accounts for the resolution of the published values:
the tabulated `5.3 Hz` and `10.9 Hz` come from an FFT of a short window,
whereas the published time history itself gives `5.47 Hz` and `10.95 Hz`.

## Coupling study

The `coupling` study runs the 1x mesh with IQN-ILS and with Robin-Neumann
coupling from the coupling start at `t = 2 s` to `t = 2.3 s` and compares the
two transients sample by sample. Both couplings converge the same interface
problem at every time step, so the histories must agree from the first coupled
step: the maximum difference in each of `u_x`, `u_y`, drag and lift, normalised
by the maximum of the IQN-ILS history over the window, must be within `1%`.
The driver also reads the recorded residual history to ensure that every
coupled Robin time step terminates with the displacement, pressure-change and
leakage-flux residuals below their configured tolerances, and it reports the
number of FSI iterations of each coupling.

### Recorded coupling result: FAIL

The first complete coupling study (OpenFOAM v2412, 1x mesh, four ranks each)
did not pass, and the failure is recorded here because it is a finding about
the Robin-Neumann implementation on this case rather than about the
benchmark. Both couplings converged every coupled step (IQN-ILS: 12.1
iterations per step on average, at most 14; Robin-Neumann: 160 on average, at
most 225, with the worst converged pressure-change and leakage-flux residuals
at `1e-5`), yet the transients differ from the first coupled step. At
`t = 2.001 s` the converged force on the plate is `+0.64 N` in the x direction
with IQN-ILS but `-9.38 N` with Robin-Neumann, and the Robin-Neumann plate
then oscillates with a `u_y` amplitude of about `18 mm` within `0.1 s`, whereas
the IQN-ILS oscillation grows gradually from below `1 mm`, as in the tutorial.
Over the window to `t = 2.3 s` the maximum differences are `18` times the
IQN-ILS maximum for `u_x`, `12` times for `u_y`, `40%` for the drag and
`4.3` times for the lift. Two converged partitioned couplings of the same
discrete problem should agree to the interface tolerance, so this points at
the Robin condition's treatment of the coupling start (the flow is developed
and the plate at rest when the coupling is switched on at `t = 2 s`) or of a
thin plate wetted on both sides; the `beamInCrossFlow` coupling study, where
the coupling is active from a ramped start, agrees to `1%`. The discrepancy is
left open here and the IQN-ILS results above are the verification of record.

The study also stops well short of the periodic regime because the
Robin-Neumann coupling is expensive on this case. The plate is thin, wetted on
both sides and as dense as the fluid, and its response involves several
bending modes with different interface impedances, which a single Robin
coefficient cannot match: the automatic `secant` coefficient settles at about
twice its `thicknessLimited` seed with a predicted contraction factor of `0.92`
per iteration, so every coupled step needs `150` to `225` unrelaxed
fixed-point iterations to reach the displacement tolerance, against about `12`
IQN-ILS iterations. Larger coefficients (`hsModel pWaveSpeed`, or a constant
`hs` of `0.03 m` or more) diverge at the impulsive coupling start, and smaller
ones converge even more slowly. The `robin` mesh study is therefore possible
(`--coupling robin`) but is not expected to be run routinely.

## Reference results

To be recorded from the first complete runs.

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
