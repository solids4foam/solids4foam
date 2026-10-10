# HronTurek verification studies

This directory contains opt-in verification studies for the `HronTurek`
tutorial. By default they compare the periodic response of the Turek-Hron FSI3
benchmark, computed with the partitioned Dirichlet-Neumann IQN-ILS and
Robin-Neumann couplings, with the Featflow reference values of the Turek and
Hron benchmark (see the provenance section below). The steady FSI1 and the periodic FSI2 benchmarks are described
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
./Allverify --levels 2 --delta-t 0.00025          # time-step diagnostic
./Allverify --levels 2 --outer-corr-tolerance 1e-6 --allow-unconverged-coupling
python3 scripts/fsi3_refinement_analysis.py \
    --diagnostic iqnils_mesh_2x_dt0.00025      # three-level analysis of 1x, 2x, 4x
```

`--delta-t` and `--outer-corr-tolerance` override the time step and the
interface tolerance of every level of an FSI3 or FSI2 mesh study; the runs
and results carry a `_dt<value>` or `_tol<value>` suffix.
`--allow-unconverged-coupling` lets a tighter-tolerance diagnostic accept a
step that stalls at `nOuterCorr` instead of aborting; the residual file
records it and the analysis counts such steps.

The driver requires `python3`, `blockMesh` and `solids4Foam`; `gnuplot` is
optional and is used for the history plots. Each run is a complete copy of the
tutorial under `verification/work/`, so the tutorial itself and its regression
test are not modified. Results are written to `verification/postProcessing/`
as CSV files, `verification_summary.md`, and PNG history plots. Each study
appends its section to `verification_summary.md`, so that consecutive studies,
such as FSI1 then FSI2, are kept together; delete it to start afresh. Both
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

Both the tutorial and verification copies integrate the forces over the
`cylinder` and `plate` patches with `rhoInf 1000`. The plotted and verification
forces are divided by the `0.015 m` mesh thickness to report force per unit
depth, like the published values.

Both also use `StVenantKirchhoffElastic`, the benchmark's constitutive law,
and `outerCorrTolerance 1e-5` for FSI3. At `1e-6`, the IQN-ILS residual
occasionally stalls just above the tolerance and a long run aborts (about
once in five thousand steps).

The FSI3 verification run continues to `t = 7 s` (the tutorial stops at
`6 s`) so that the closing window is well inside the periodic regime.
The coupling is activated at `t = 2 s`, as in the tutorial.

The benchmark is the FSI3 test of Turek and Hron as published: a channel of
`2.5 m × 0.41 m`, a cylinder of radius `0.05 m` centred at `(0.2, 0.2)`, a
plate of `0.35 m × 0.02 m` attached to it with its tip at point A,
`(0.6, 0.2)`, a parabolic inflow of mean `2 m/s`, a fluid of
`ρ = 1000 kg/m^3` and `ν = 0.001 m^2/s`, and a St. Venant-Kirchhoff plate of
`ρ = 1000 kg/m^3`, `E = 5.6 MPa` and `ν = 0.4` in plane strain. The inflow is
applied without the benchmark's `2 s` start-up ramp, and the coupling starts
at `t = 2 s`; neither affects the periodic state that is evaluated.

### FSI3 reference values and their provenance

The reference values are in `reference/HronTurek_verification_references.json`,
which also records all nine rows of the Featflow FSI3 tables
(`featflowTables`) and the 2006 tables (`turekHron2006Tables`). Three sets of
values are in circulation, and only the first is used for verification:

- **A. Featflow FSI3 table, level 4+0 (15 872 elements), `Δt = 0.00025 s`**
  (primary): `u_x = -2.88 ± 2.72 mm [10.93 Hz]`,
  `u_y = 1.47 ± 34.99 mm [5.46 Hz]`, `F_D = 460.5 ± 27.74 N/m [10.93 Hz]` and
  `F_L = 2.50 ± 153.91 N/m [5.46 Hz]`. Source: the "Results for FSI3 with
  timestep Δt=0.00025" table of the Featflow FSI tests page (retrieved
  2026-10-03). This is the discretisation of the published reference history
  `ref_fsi3.point` (`media/fsi/data/fsi3/0p00025`), from which
  `reference/TurekHron_fsi3_reference_history.csv` is subsampled to `1 ms`.
- **B. Turek and Hron (2006)** (historical only):
  `u_x = -2.69 ± 2.53 mm [10.9 Hz]`, `u_y = 1.48 ± 34.38 mm [5.3 Hz]`,
  `F_D = 457.3 ± 22.66 N/m` and `F_L = 2.22 ± 149.78 N/m`. These are level 4
  of the 2006 proceedings table at the smaller of its two time steps
  (`Δt = 0.0005 s`). The Featflow page keeps the 2006 tables only inside an
  HTML comment marked "old values"; they were superseded by table A. Its
  drag amplitude is `18%` below A, and its frequency is given to two digits.
- **C. Tuković et al. (2018)** (historical only): the earlier solids4foam
  result, `u_x = -2.72 ± 2.58 mm [11.07 Hz]`, `u_y = 1.67 ± 33.84 mm [5.53 Hz]`,
  `F_D = 459.18 ± 24.86 N/m` and `F_L = 1.59 ± 155.9 N/m`.

The driver's extraction, applied to the reference history at its full
`0.25 ms` sampling, reproduces table A: the displacement statistics, the drag
mean and the lift amplitude to within `0.6%`, the frequencies to within
`0.3%`, and the drag amplitude to `0.8%` (`27.51` against `27.74 N/m`;
`27.43 N/m` from the `1 ms` CSV). This checks the extraction itself.

The reference is not exact, and its uncertainty differs between the
quantities:

- From level 3 to level 4 at `Δt = 0.00025 s`, table A changes by `3.8%` in
  the `u_x` mean, `4.0%` in the `u_x` amplitude, `1.6%` in the `u_y`
  amplitude, `0.3%` in the drag mean, `4.5%` in the drag amplitude and `2.6%`
  in the lift amplitude. The frequencies do not change.
- Levels 2, 3 and 4 are not monotone in the `u_x` mean, the `u_x` and `u_y`
  amplitudes and the drag amplitude (`35.73`, `34.43`, `34.99 mm` in `u_y`).
  The lift amplitude rises by about `4 N/m` per level (`146.0`, `149.9`,
  `153.9 N/m`) with no sign of convergence.
- At level 4, the three tabulated time steps differ by up to `0.7%` in the
  lift amplitude and `1.0%` in the drag amplitude.
- The reference history itself is not strictly periodic: over its `7.9`
  periods (`t = 5` to `6.44 s`) the per-period lift amplitude falls from
  `157.7` to `153.5 N/m` (`2.7%`) and the drag amplitude from `27.8` to
  `27.5 N/m`.

An uncertainty of about `1%` is therefore defensible for the drag mean and
the frequencies only. For the `u_y` amplitude it is about `2%`, and for the
`u_x` mean and amplitude and the drag and lift amplitudes it is `3` to `5%`
at least.

The driver overlays the published history on the closing window of each run,
with the phases aligned at the last `u_y` maximum.

## Mesh levels and time steps

Level 1 is the mesh shipped with the tutorial. Each further level doubles the
in-plane block divisions of both the fluid and the solid mesh, leaving the
single spanwise cell alone, and halves the time step so that the Courant number
is unchanged. The block gradings are kept, so the meshes are a smooth family
with a refinement ratio of 2 rather than strictly nested; the plate has 6, 12
and 24 cells through its thickness. Because the time step is halved with the
cell size, the sequence is a combined space-time refinement path, and an
order observed along it is not a purely spatial order:

| Level | Refinement | Fluid cells | Solid cells | Δt (s) | Default cores |
|---:|---:|---:|---:|---:|---:|
| 1 | 1x | 5 336 | 630 | 0.001 | 1 |
| 2 | 2x | 21 344 | 2 520 | 0.0005 | 8 |
| 3 | 4x | 85 376 | 10 080 | 0.00025 | 16 |

Levels 1 and 2 form the default sweep; level 3 is reachable with
`--levels 1,2,4` and is expensive. Use `--cores N` to run every level on the
same number of MPI ranks (`--cores 1` runs in serial). The recorded 4x run
used 32 ranks, which gives each rank about as many cells as the 2x run on
eight; on one 128-core node it ran at `9.1 s` per step on 32 ranks, `11.1 s`
on 64 and `18.9 s` on 128, so more ranks do not help. Level 1 is serial by
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
  means, which are close to zero, are reported as diagnostics only. The
  tolerances were set from the 1x and 2x results and are not changed here.
  With `--levels 1,2,4` the study fails on the 4x drag amplitude (`18.1%`
  against `15%`); the lift amplitude passes at `15.0%`. See the three-level
  study below: both amplitudes are inflated at `Δt = 0.00025 s` by
  coupling noise, and the drag amplitude is still rising.
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

### Coupling comparison

Recorded with OpenFOAM v2412 on an Apple M1 Ultra shared with other jobs.

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

### Three-level mesh study (1x, 2x, 4x)

Recorded with OpenFOAM v2512 on a shared 192-core node (xenosim), solids4foam
`development` at `aee0c8e35`. IQN-ILS with `predictor yes` and
`outerCorrTolerance 1e-5`, run to `t = 7 s` and evaluated over the closing
`1 s`. The compact results are in `reference/fsi3_refinement/`:
`fsi3_iqnils_refinement_study.json` and the two CSV files (per level and per
quantity), the driver's `iqnils_mesh_sweep_1x2x4x.csv`, and the histories of
every run from `t = 4.5 s` at `1 ms`. `reference/iqnils_mesh_4x_history.png`
overlays the 4x history on the published one.

| Level | Fluid cells | Solid cells | Δt (s) | Ranks | Wall time | Mean (max) FSI iterations | Largest final residual |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 1x | 5 336 | 630 | 0.001 | 1 | 1.4 h | 6.96 (11) | 9.994e-6 |
| 2x | 21 344 | 2 520 | 0.0005 | 8 | 5.5 h | 7.01 (13) | 9.999e-6 |
| 4x | 85 376 | 10 080 | 0.00025 | 32 | 40.8 h | 6.69 (15) | 9.9996e-6 |

All three runs completed. Every one of the 5 000, 10 000 and 20 000 coupled
steps ended with the interface residual below `1e-5`; an unconverged step
aborts the run, so this holds by construction and is confirmed from the
residual files. The iteration count does not grow with refinement. In the
analysis window the largest two-against-two amplitude change was `1.9%`
(4x drag), and the 4x lift passed the `2%` periodicity test at `1.75%`.

The 1x and 2x results reproduce the earlier record (OpenFOAM v2412 on an
Apple M1 Ultra) to within `0.5%` on the 2x mesh, for example
`174.14` against `174.4 N/m` in the lift amplitude. On the 1x mesh they agree
to within `2.4%` (`22.87` against `23.43 N/m` in the drag amplitude). A 1x
run with OpenFOAM v2412 on MeluXina agrees with the v2512 1x run to within
`1.6%` in all primary quantities.

Benchmark statistics (extrema over the last full `u_y` period), and errors
against Featflow level 4 (signed, relative to the reference):

| Quantity | 1x | 2x | 4x | Featflow L4 | Error 1x | Error 2x | Error 4x |
|---|---:|---:|---:|---:|---:|---:|---:|
| `u_x` mean (mm) | -2.232 | -2.788 | -3.085 | -2.880 | +22.5% | +3.2% | -7.1% |
| `u_x` amplitude (mm) | 2.163 | 2.657 | 2.917 | 2.720 | -20.5% | -2.3% | +7.2% |
| `u_x` frequency (Hz) | 11.176 | 11.046 | 10.950 | 10.930 | +2.3% | +1.1% | +0.2% |
| `u_y` mean (mm) | 1.669 | 1.495 | 1.467 | 1.470 | - | - | - |
| `u_y` amplitude (mm) | 29.85 | 34.18 | 36.23 | 34.99 | -14.7% | -2.3% | +3.5% |
| `u_y` frequency (Hz) | 5.588 | 5.523 | 5.474 | 5.460 | +2.3% | +1.2% | +0.3% |
| drag mean (N/m) | 456.6 | 459.7 | 462.1 | 460.5 | -0.8% | -0.2% | +0.3% |
| drag amplitude (N/m) | 22.87 | 28.03 | 32.76 | 27.74 | -17.6% | +1.0% | +18.1% |
| lift mean (N/m) | 0.82 | 0.25 | 6.00 | 2.50 | - | - | - |
| lift amplitude (N/m) | 229.1 | 174.1 | 177.0 | 153.9 | +48.8% | +13.2% | +15.0% |
| lift frequency (Hz) | 5.588 | 5.523 | 5.493 | 5.460 | +2.4% | +1.2% | +0.6% |

The drag and lift frequencies equal the `u_x` and `u_y` frequencies; the 4x
lift frequency is biased by `0.3%` by the coupling noise described below,
which moves the zero crossings.

#### Coupling noise in the 4x forces

At `Δt = 0.00025 s` the `1e-5` interface tolerance leaves step-to-step noise
in the forces: an rms of `7.3 N/m` in the lift, with peaks of about
`35 N/m`, and `1.1 N/m` in the drag on the 4x mesh, against `1.7` and
`0.3 N/m` on the 2x mesh at `Δt = 0.0005 s`. The same noise appears on the 2x
mesh when the time step alone is halved (`7.2 N/m`), so it comes from the
time step, not the mesh. Because the benchmark amplitude is taken from the
extrema, the noise inflates the 4x drag and lift amplitudes. The
displacements are not affected. The analysis therefore also reports the
statistics of the histories smoothed by a centred `4 ms` moving average,
which attenuates the `5.5 Hz` fundamental by less than `0.1%`:

| Quantity | 1x | 2x | 4x | Featflow L4 | Error 4x | Ratio (4x-2x)/(2x-1x) | Order |
|---|---:|---:|---:|---:|---:|---:|---:|
| drag amplitude, smoothed (N/m) | 22.63 | 27.63 | 30.63 | 27.74 | +10.4% | 0.60 | 0.74 |
| lift amplitude, smoothed (N/m) | 228.78 | 172.27 | 165.56 | 153.91 | +7.6% | 0.12 | 3.1 |

#### Observed order along the refinement path

The table gives the successive differences, their ratio
`R = (f_4x - f_2x)/(f_2x - f_1x)`, and the observed order
`p = log2(1/R)` along the space-time refinement path. An order is reported
only where the differences have the same sign and decrease, the quantity is
not near zero, and the 2x-to-4x change exceeds the late-time variability. The
late-time variability is the range of the per-period values from `t = 4.5 s`
to the end of the run (13 periods). After the amplitudes saturate, the means,
amplitudes and frequency still wander slowly, by up to about `1%` over a
second or more, which the `1 s` window and the periodicity test do not see.
On the 4x mesh the late-time variability is `1.0%` in the `u_x` mean, `0.8%`
in the `u_x` amplitude, `0.6%` in the `u_y` amplitude, `0.1%` in the
frequency, `0.8%` in the drag mean and `1.0%` in the smoothed lift amplitude.

| Quantity | 2x - 1x | 4x - 2x | Ratio | Order | Status |
|---|---:|---:|---:|---:|---|
| `u_x` mean (mm) | -0.556 | -0.297 | 0.53 | 0.91 | monotone, not asymptotic |
| `u_x` amplitude (mm) | 0.494 | 0.260 | 0.53 | 0.93 | monotone, not asymptotic |
| `u_y` amplitude (mm) | 4.32 | 2.06 | 0.48 | 1.07 | monotone, not asymptotic |
| `u_y` frequency (Hz) | -0.065 | -0.049 | 0.75 | 0.42 | monotone, not asymptotic |
| drag mean (N/m) | 3.10 | 2.35 | 0.76 | - | 2x-to-4x change within the variability |
| drag amplitude (N/m) | 5.16 | 4.73 | 0.92 | 0.13 | monotone, not asymptotic; noise-inflated |
| drag amplitude, smoothed (N/m) | 5.00 | 3.00 | 0.60 | 0.74 | monotone, not asymptotic |
| lift amplitude (N/m) | -54.9 | 2.8 | -0.05 | - | sign change within the variability |
| lift amplitude, smoothed (N/m) | -56.5 | -6.7 | 0.12 | 3.1 | monotone, order above 2: not demonstrably asymptotic |
| `u_y` and lift means | | | | - | near zero: undefined |

The schemes are formally second order in space and time. No quantity shows
an observed order consistent with that (taken as `1.5 ≤ p ≤ 2.5`), so no
Richardson extrapolation is meaningful, and none is used.

#### Time-step and coupling contributions

The time step halves with the cell size, so the observed orders mix the two.
Four runs on the 2x mesh separate the contributions there: `Δt = 0.0005` and
`0.00025 s`, each with the interface tolerance `1e-5` and `1e-6`. The run at
`Δt = 0.00025 s` and `1e-6` stopped at `t = 6.73 s` (see below) and is
evaluated over `t = 5.73` to `6.73 s`; its 18 923 completed steps all
converged. The changes are relative, from the smoothed statistics for the
amplitudes and means; the changes in the means are in their magnitudes:

| Quantity | Halving Δt at 1e-5 | Halving Δt at 1e-6 | 1e-5 to 1e-6 at Δt = 0.0005 | 1e-5 to 1e-6 at Δt = 0.00025 | 2x to 4x |
|---|---:|---:|---:|---:|---:|
| `u_x` mean | +1.4% | -0.7% | +1.7% | -0.4% | +10.6% |
| `u_x` amplitude | +1.4% | -0.4% | +1.7% | -0.2% | +9.8% |
| `u_y` amplitude | +0.8% | -0.0% | +0.2% | -0.7% | +6.0% |
| `u_y` frequency | -0.1% | +0.4% | -0.4% | +0.0% | -0.9% |
| drag mean | +0.1% | -0.2% | +0.2% | -0.1% | +0.7% |
| drag amplitude | +3.3% | +2.5% | -1.4% | -2.2% | +10.9% |
| lift amplitude | +0.3% | +1.2% | -0.5% | +0.4% | -3.9% |

The changes from the time step and the tolerance are small against the
2x-to-4x change for the displacements. They are of the same size as the
change for the frequency, about a third of it for the smoothed lift
amplitude, and a quarter of it for the smoothed drag amplitude. The
frequency is therefore not resolved beyond about `0.4%`. On the
benchmark-defined extrema, halving the time step at `1e-5` raises the drag
amplitude by `11%` and the lift amplitude by `3.3%`; tightening the
tolerance to `1e-6` at `Δt = 0.00025 s` removes most of that
(`-9.8%` and `-3.2%`). The force noise is coupling error, not time
discretisation error.

The tighter tolerance cannot be used for the 4x level itself. At
`Δt = 0.00025 s` and `outerCorrTolerance 1e-6`, IQN-ILS twice stalled just
above the tolerance (`5` to `9 × 10^-6`, at `nOuterCorr 30`) and then
diverged within a step, after which the solid SNES failed: on the 2x mesh
at `t = 6.73 s`, and on the 4x mesh, restarted at `t = 4 s` from the `1e-5`
run, at `t = 4.73 s`. The fluid solver tolerances (`1e-6`, absolute) put a
floor of about `5e-6` on the relative interface residual at this time step,
as they do for FSI1 at `1e-4`. A clean `1e-6` 4x run would need tighter fluid
tolerances.

#### Interpretation

Self-convergence. Every primary quantity except the benchmark-defined lift
amplitude changes monotonically from 1x through 2x to 4x; for the drag mean
the 2x-to-4x change (`0.5%`) is within its late-time variability. The frequency converges at a decreasing rate. The displacements change
by `6` to `11%` from 2x to 4x with an observed order of about `1`, half the
formal order, so they are not yet in the asymptotic range. The smoothed lift
amplitude changes by only `3.9%` from 2x to 4x after `25%` from 1x to 2x,
but its apparent order of `3.1` exceeds the formal order and rests on three
points, so it does not demonstrate an asymptotic range either. The
smoothed drag amplitude is still rising (`+10.9%`).

Agreement with Featflow. The 2x displacements and drag amplitude were within
`1` to `3%` of Featflow level 4, but the 4x mesh moves past it: the `u_y`
amplitude to `+3.5%` and the `u_x` mean and amplitude to `7%` beyond, still
changing in the same direction. The 2x agreement was a crossing, not
convergence. The frequency converges onto Featflow (`+0.25%` at 4x, within
the `0.4%` time-step and coupling uncertainty), and the drag mean agrees to
`0.3%`, within its late-time variability of `0.8%`. The lift amplitude does
not close onto Featflow: smoothed, it is `11.9%` above on the 2x mesh and
`7.6%` above on the 4x mesh, and its small last change points to a limit of
about `164` to `166 N/m`, `6` to `8%` above `153.9 N/m`. The
benchmark-defined value at 4x is `15.0%` above, inflated by coupling noise.

Reference uncertainty. The remaining differences in the lift amplitude and
the displacements exceed the `1%` sometimes assigned to Featflow level 4 but
are comparable to the reference's own level-3-to-4 changes (`2.6%` in the
lift amplitude, `4%` in `u_x`), and the Featflow lift amplitude itself rises
by `4 N/m` per level without converging. Whether solids4foam converges to a
value different from the exact solution, or Featflow level 4 is itself
`5%` or more from it, cannot be settled from these data. It needs either a
finer solids4foam level or the CSM3 and CFD3 sub-benchmarks, which separate
the plate and the fluid discretisations.

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
- `outerCorrTolerance 1e-6`, also used by `Allrun fsi2`, instead of the
  `1e-5` of the
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

## CSM3 structural study

`csm3/` holds a structural-only study of the FSI3 plate (the Turek-Hron CSM3
benchmark and its steady CSM2 analogue, which has the FSI3 material), used to
isolate the spatial convergence of the solid from the fluid and the coupling.
See `csm3/README.md`; it is run separately from `Allverify`.

## FSI3 solid-only refinement at fixed fluid resolution

An error-decomposition study: the 2x fluid mesh, `dt = 0.0005 s`, the IQN-ILS
coupling and its `1e-5` tolerance are fixed, and only the solid mesh is
refined (2x = 210 x 12, 4x = 420 x 24, 8x = 840 x 48 cells). The fluid
(168 interface faces) is not refined, so this is not a spatial-order study and
no order is computed.

```bash
./Allverify --levels 2 --cores 8                              # solid 2x (baseline)
./Allverify --levels 2 --solid-refinement 4 --cores 14
./Allverify --levels 2 --solid-refinement 8 --cores 14
python3 scripts/fsi3_solid_refinement_analysis.py
```

`--solid-refinement` refines the solid block mesh by that factor
independently of `--levels` (runs carry a `_solid<N>x` suffix). The analysis
writes the QoIs, the successive changes, the interface face counts, the force
balance across the AMI (rms of fluid + solid interface force over the rms
fluid force, over the analysis window) and the coupling iterations; the
compact results are in `reference/fsi3_solid_refinement/`.

## CFD3 fluid-only study (rigid flag)

`scripts/hron_turek_cfd3.py` builds a single-region, static-mesh version of the
tutorial's fluid region (same `blockMeshDict` refined by the FSI3 factor,
same schemes and solver settings, same inlet and outlet), with the plate a
fixed no-slip wall: the Turek-Hron CFD3 benchmark (rigid flag, mean inflow
`2 m/s`, `Re = 200`). It is not part of `Allverify`.

```bash
python3 scripts/hron_turek_cfd3.py --level 2 --delta-t 0.0005 --end-time 25 --cores 8
python3 scripts/cfd3_analysis.py     # reads verification/work/cfd3_*
```

The impulsively started flow needs about 15 s (1x) to reach the periodic
state, so the runs go to 25 s and the closing 1 s is analysed (3-4 lift
periods). `cfd3_analysis.py` evaluates the matched path (1x/2x/4x with
`dt = 1e-3/5e-4/2.5e-4`), the fixed-`dt` path (`dt = 2.5e-4`) and a temporal
control on 2x, with the lift period delimiting the periods of all signals,
and writes `postProcessing/cfd3_study.json` and two CSV files; the compact
copies and 4 s force histories are in `reference/cfd3/`. The published
reference is Featflow level 4+0, `dt = 0.005`. Pressure difference is not a
CFD3 benchmark quantity and is not reported; the pressure and viscous parts
of drag and lift are.

## Prescribed-motion ALE study (fluid only)

`scripts/hron_turek_ale.py` restarts the developed rigid CFD3 state of a level
(`cfd3_<level>x_dt<dt>`, 25 s) and bends the flag with a prescribed
periodic motion, on the same fluid mesh with the FSI3 mesh motion
(`velocityLaplacian`, `newMovingWallVelocity` on the plate): the transverse
displacement is `A s(t) phi(x) sin(2 pi f (t - 25))` with `A = 35 mm`,
`f = 5.5 Hz`, `phi` the clamped-free beam mode 2 normalised to unit tip value,
`s` a smooth 1 s ramp; the point velocity is the finite difference of the
displacement, so the Euler-integrated mesh follows it exactly. The motion is
a `codedFixedValue` point-patch condition (compiled at run time); a `coded`
function object prints the pressure force on the flag projected on `phi`
(`GENF` log lines). Nothing in `src/` is changed.

```bash
python3 scripts/hron_turek_ale.py --level 2 --delta-t 0.0005 --duration 8 --cores 8
python3 scripts/ale_analysis.py     # reads verification/work/ale_*
```

`ale_analysis.py` takes the Fourier components at the forcing frequency over
the last five periods (in phase with the displacement, `a1`; in phase with the
velocity, `b1`), the extrema-based amplitude, the drag mean and second
harmonic, and the pressure/viscous split, with the change between the last two
five-period windows as the noise level; the matched (1x/2x/4x with
`dt = 1e-3/5e-4/2.5e-4`) and fixed-`dt` (`2.5e-4`) paths are written to
`postProcessing/ale_study.json` and two CSV files (compact copies and 2 s
histories in `reference/ale/`).

## Trajectory replay and energy balance (fluid only)

The ALE study above prescribes an assumed mode-2 motion. The replay
replaces it by the actual coupled motion:

1. `scripts/hron_turek_replay_source.py` restarts the converged FSI3 2x run
   (`iqnils_mesh_2x`, t = 7 s) for 1.2 s and records the fluid `plate`
   point positions at every step.
2. `scripts/hron_turek_replay.py` fits them over five periods (mean plus
   four harmonics in x and y per point, f = 5.5236 Hz) and replays them on
   the CFD3 fluid mesh of any level from its developed rigid state (25 s,
   1 s ramp, 8 s). The setup is the same as in the ALE study.
3. `scripts/replay_analysis.py` gives the force harmonics.

```bash
python3 scripts/hron_turek_replay_source.py --run iqnils_mesh_2x --cores 8 --duration 1.2
python3 scripts/hron_turek_replay.py --level 2 --delta-t 0.0005 --duration 8 --cores 8
```

`scripts/hron_turek_energy.py` restarts a finished replay from a written
time (1.2 s; 0.3 s blend and 1.6 s for the perturbed runs). It fixes the
replay boundary condition so that it keys the table on the `constant`
points, which a restart needs. It optionally scales the amplitude or the
frequency of the motion smoothly, with continuous phase. A `coded` function
object writes `energy.dat` at every step. The force is
`rho (p Sf + Sf . devReff)`, as in the `forces` function object; its plate
totals match `forcesPlate` to all printed digits. The file contains:

- the pressure and viscous power on the flag;
- the discrete work increments `sum_f f.(d^n - d^(n-1))` of each harmonic
  and direction;
- the generalised forces on the harmonic shapes.

`scripts/energy_analysis.py` evaluates, over the last five periods:

- the net work per cycle W, which is positive into the flag;
- its split by harmonic;
- the first-harmonic generalised force on the replayed y shape, normalised
  to the tip amplitude: `Q_in` is in phase with the displacement
  (added stiffness or mass), and `Q_quad` is in phase with the velocity,
  with `W_1y = pi A Q_quad`;
- a restart check against the `POWR` trace of the original replay.

`scripts/energy_summary.py` combines this with the coupled FSI3 levels.

All runs use the coupled 2x motion (tip amplitude 34.30 mm, 5.5236 Hz),
v2512, GCC 11.4 `-O3`, xenosim. The work per cycle is in J/m, the forces
in N/m. The work increments are weighted by the shape of each face, and
their sum (`W_dd`) differs from the power-based W by `O(omega dt)`: the
fluid wall velocity is a BDF2 mesh-flux velocity. That difference is
0.12 J/m at dt = 1e-3 and 0.03 J/m at dt = 2.5e-4.

| run | W | W_dd | W_p | W_v | Q_in | Q_quad | gross \|P\| work |
|---|---:|---:|---:|---:|---:|---:|---:|
| 1x, dt 1e-3 | -1.987 | -1.870 | -2.073 | 0.086 | 50.34 | -16.97 | 8.38 |
| 1x, dt 2.5e-4 | -1.745 | -1.717 | -1.874 | 0.129 | 47.75 | -14.70 | 8.38 |
| 2x, dt 5e-4 (source mesh) | 0.082 | 0.161 | -0.123 | 0.205 | 72.17 | 2.00 | 7.40 |
| 2x, dt 2.5e-4 | 0.141 | 0.181 | -0.071 | 0.212 | 72.67 | 2.55 | 7.34 |
| 4x, dt 2.5e-4 | 0.125 | 0.170 | -0.103 | 0.229 | 81.02 | 2.34 | 6.89 |
| 2x, amplitude x0.95 | 0.018 | 0.102 | | | 81.03 | 1.14 | 5.98 |
| 2x, amplitude x1.05 | -0.031 | 0.041 | | | 60.83 | 1.19 | 9.02 |
| 2x, frequency x0.99 | 1.326 | 1.402 | | | 69.66 | 13.39 | 7.24 |

### Validation

- The replayed coupled motion on its own mesh gives W = 0.08 to 0.18 J/m,
  depending on dt and on the work definition. That is 1 to 2% of the gross
  exchange, but 10 to 20% of E1 = 0.81 J/m, so zero is resolved only to
  about 0.1 J/m.
- The plate lift fundamental (179 N/m) matches the coupled lift amplitude
  (174 N/m).
- Each restart reproduces the original pressure-power trace: the largest
  difference is 0.38 W/m at a rms of 41 to 51 W/m, and the work difference
  is at most 0.0015 J/m per cycle.
- W changes by less than 0.003 J/m between windows shifted by half a period.

### Findings

- **1x mesh.** At fixed motion the 1x mesh extracts 1.9 to 2.1 J/m per
  cycle more than 2x. Nearly all of it is pressure work on the first
  transverse harmonic. The sign matches the smaller 1x coupled amplitude.
  Read literally as numerical dissipation, this is an overstatement: the
  1x fluid response has a different phase.
- **Energy route is indeterminate for 2x to 4x.** The 4x mesh changes W by
  -0.016 J/m (fixed dt) or +0.044 J/m (matched dt). This is within the
  resolution of W, about ±0.05 J/m from dt and from the definition of the
  discrete work.
  - At fixed shape and frequency, W is almost flat in amplitude: ±5%
    amplitude gives +0.02 and -0.03 J/m, against +0.08 J/m at the baseline.
    The secant is about -0.014 J/m per mm, but curvature or noise is not
    resolved.
  - With that slope, 0.03 J/m would already account for the +2 mm 2x to
    4x change. The energy argument therefore can neither confirm nor
    exclude the mechanism there.
  - For 1x it fails quantitatively: it predicts a change of about 100 mm.
  - The implied dA/dW is incompatible between levels: 2.1 to 2.3 mm per
    J/m from 1x, against +47 or -129 mm per J/m from 4x.
- **W depends strongly on frequency** (-22.5 J/m per Hz). Extrapolated to
  the observed 0.9% frequency decrease from 2x to 4x, the frequency shift
  adds about 1 J/m, which the amplitude cannot remove. Shape and harmonic
  changes must close the balance, and the replay of a single motion
  cannot show how.
- **In-phase force (hypothesis).** What changes measurably with the mesh at
  fixed motion is the in-phase generalised force: Q_in = 50.3, 72.2 and
  81.0 N/m (+12% from 2x to 4x). Q_in also falls steeply with amplitude
  (-5.9 N/m per mm).
  - A one-mode in-phase balance at fixed frequency,
    `dA = -dQ_in/(dQ_in/dA)`, optionally with a secant structural term
    `Q_in/A`, gives -2.7 to -4.2 mm for 1x (observed -4.32) and +1.0 to
    +1.5 mm for 4x (observed +2.05).
  - This is a consistent hypothesis, not an identified mechanism. At one
    frequency, Q_in cannot separate added mass, stiffness and reactive
    vortical load. The structural tangent is not measured. The energy
    residual is not closed.
- **u_x mean.** The time-mean tip u_x is kinematic: the inextensible
  estimate `-1/4 int |W'|^2 dx` gives -2.66 mm against -2.70 mm fitted.
  Across the coupled levels the time-mean u_x scales as A^1.88 to 1.93,
  and A^2 is within 1.6% (1x) and 0.4% (4x). The benchmark midrange
  (max+min)/2 scales as A^1.64 to 1.70, because the 2f waveform changes
  too.

`reference/energy_balance/` holds:

- `energy_runs.json`, the per-run analysis;
- `energy_balance.csv` and `energy_balance.json`, the summary and
  predictions;
- `defect_exposure.json`, the PR #546 exposure check.

### PR #546 exposure

The study binary was built at aee0c8e35 with GCC 11.4 `-O3`; it has the
same `src/` as 27c9bfdf7.

- **Aliasing defect.** The alias reproducer fails at `-O3` on xenosim.
  For the original vector & symmTensor product the error is in the z
  component only, and it vanishes for 2-D data
  (`platform/Test-tmpDotAlias2D.C`).
- **Barycentric weights.** The degenerate-weight branch cannot trigger
  here: the smallest interface fan triangle is 9.8e3 times above the
  threshold, at solid 8x (`scripts/plate_interface_faces.py`).
- **A/B against the fixed tree.** 27c9bfdf7 + 28751407c + c8efdd390 was
  built with the same toolchain. Against the study binary, all of these
  are identical to every written digit:
  - an FSI3 2x coupled restart (600 steps): forces, point A, every plate
    point and the residual file;
  - a CFD3 2x restart (500 steps), run with the 28751407c build only (the
    rigid flag uses `fixedValue`, not `newMovingWallVelocity`);
  - a replay 2x restart (600 steps), whose run exercises
    `newMovingWallVelocity`.

  The comparison script is `scripts/ab_compare.py`.

## Fluid-side variants of the trajectory replay (Q_in convergence)

At fixed coupled 2x motion the in-phase generalised force on the flag,
`Q_in`, is 50.3, 72.2 and 81.0 N/m at 1x, 2x and 4x (dt 1e-3, 5e-4, 2.5e-4;
observed order 1.3). This section tests which fluid treatment limits that
convergence. Each variant changes one treatment of the replay restart cases
of the section above, the rest of the setup is unchanged. The analysis is the
same (`energy_analysis.py`, last five periods, windows ending 33.2 s at 1x
and 2x and 30.2 s at 4x).

Tools: `scripts/hron_turek_variants.py` (builds the case with
`hron_turek_energy.py`, then applies the variants), `scripts/run_replay_variant.sh`,
`scripts/variants_summary.py`, and `platform/htVariants/` (a small library with
the GCL mesh class). Results: `reference/replay_variants/`
(`variants.csv` and `.json`: all runs; `convergence.csv`: per variant).

Setup: parent commit c821aeb8b, study binary built from aee0c8e35 (PR #546
has no effect), OpenFOAM v2512 (OpenCFD), GCC 11.4.0 `-O3`, Ubuntu 22.04,
xenosim, run on 1 (1x), 8 (2x) and 24 (4x) MPI ranks. The baseline restart at
2x reproduces `Q_in = 72.172` of the earlier study to all digits.

### Variants, as checked in the source

The FSI3 fluid has `p` zeroGradient and `U` `newMovingWallVelocity` on the
flag, `default leastSquares`, `ddt backward`, `velocityLaplacian` with
`diffusivity quadratic inverseDistance`. There is no `fixedValueCorrected`
in v2512 or solids4foam.

| name | change (dictionary) | what it does (source) |
|---|---|---|
| pwall | `0/p` plate: `type movingWallPressure;` | solids4foam: `gradient = -n.a_wall`, `a_wall` the BDF2 wall acceleration from `newMovingWallVelocity`; `pimpleFluid` adds `rAU*snGrad(p)*|Sf|` to `phiHbyA`, so the wall flux stays the mesh flux |
| uwall | `0/U` plate: `type movingWallVelocity;` | v2512: Euler face-centre velocity, normal part from meshPhi; drops the non-orthogonal `snGrad` of `newMovingWallVelocity` |
| gradS4f | `fvSchemes`: `grad(U) leastSquaresS4f; grad(p) leastSquaresS4f;` | solids4foam: the least-squares fit uses boundary values only on patches that fix the value, so a zeroGradient `p` wall is extrapolated; `default` cannot be changed because `cellMotionU` has no face-value list. It also changes the viscous force evaluation |
| gcl | `dynamicMeshDict`: `dynamicFvMesh dynamicMotionSolverBDF2FvMesh;`, `controlDict`: `libs ("libhtVariants.so");` | mesh flux `1.5 phiE^n - 0.5 phiE^(n-1)`, so that `div(meshPhi)` equals the backward volume derivative `(3V - 4V0 + V00)/(2 dt)`. The standard flux is the Euler swept-volume flux, a GCL violation of `O(dt)` for backward. One Euler step at a restart |
| mesh_inv | `dynamicMeshDict`: `diffusivity inverseDistance 2(plate cylinder);` | linear instead of quadratic inverse distance |

`pwall` needs `newMovingWallVelocity`, so it cannot be combined with
`uwall`. The mesh variant is only a screening: the velocity-Laplacian mesh is
path dependent, so a restart changes the interior mesh slowly, and it was not
run at 4x.

### Results

`Q_in` in N/m (1x and 2x are full runs from t = 25 s or restarts from 32 s;
4x are restarts from 29 s):

| variant | 1x | 2x | 4x | 1x to 2x | 2x to 4x | order |
|---|---:|---:|---:|---:|---:|---:|
| baseline | 50.3 | 72.2 | 81.0 | +21.9 | +8.8 | 1.31 |
| pwall | 30.3 | 62.7 | 77.3 | +32.4 | +14.6 | 1.15 |
| uwall | 47.1 | 66.8 | | +19.6 | | |
| gradS4f | 11.0 | 50.8 | 69.3 | +39.7 | +18.5 | 1.10 |
| gcl | 57.1 | 76.1 | 82.7 | +19.0 | +6.6 | 1.51 |
| mesh_inv | 48.8 | 69.3 | | +20.6 | | |
| pwall + gradS4f + gcl | 39.4 | 65.0 | 76.2 | +25.6 | +11.2 | 1.19 |

Other 1x/2x combinations (restart from 32 s at 2x): pwall + gradS4f 32.8 /
61.5, pwall + gcl 36.9 / 66.4, gradS4f + gcl 17.6 / 54.3.
`Q_quad` (N/m), `W` (J/m) and the force amplitudes are in `variants.csv`;
for the 4x runs they are, baseline / pwall / gradS4f / gcl / all three:

| 4x | Q_quad | W | lift amplitude (N/m) | mean drag (N/m) |
|---|---:|---:|---:|---:|
| baseline | 2.34 | 0.125 | 163.3 | -18.6 |
| pwall | 2.48 | 0.140 | 161.8 | -17.3 |
| gradS4f | 1.90 | 0.084 | 167.5 | -14.4 |
| gcl | 1.18 | -0.033 | 161.1 | -19.2 |
| pwall + gradS4f + gcl | 0.77 | -0.066 | 163.6 | -15.6 |

The order is `log2` of the ratio of the two differences, for sequences that
are monotone in every row with a value (spatial and temporal refinement are
refined together, so it is a combined order).

### Validation

- A restart from 32 s with a variant switched on reproduces the same
  variant run from the developed rigid state at t = 25 s (2x): `Q_in` 62.90
  against 62.74 (pwall), 75.66 against 76.07 (gcl), 50.88 against 50.75
  (gradS4f). The switching transient is at most 0.4 N/m, so the 4x restarts
  are valid. This does not hold for `mesh_inv` (72.14 against 69.35), which
  is why the mesh variant is only a screening.
- The time step alone is minor (earlier study): at 1x, dt 1e-3 to 2.5e-4
  changes `Q_in` by -2.6 N/m, at 2x, dt 5e-4 to 2.5e-4 by +0.5 N/m.
- The `gcl` run selects `dynamicMotionSolverBDF2FvMesh` (log) and gives a
  first step identical to the baseline, as it should (no old flux).

### Findings

- No variant removes the slow convergence of `Q_in`. Every one that was run at
  all three levels is monotone with an observed order of 1.1 to 1.5.
- The `gcl` variant is the only one that helps, modestly: the 2x to 4x change
  falls from +8.8 to +6.6 N/m (-25%) and the order rises to 1.5. It moves the
  in-phase force up and `Q_quad` and `W` down at every level.
- `pwall` and `gradS4f` change the 1x value strongly (to 30 and 11 N/m) and
  give a larger 2x to 4x change than the baseline. They remove the
  zero-normal-gradient wall pressure error of the standard treatment, but
  `Q_in` then approaches the common value from further below.
- Across variants the spread of `Q_in` is 46, 25 and 14 N/m at 1x, 2x and
  4x: the treatments differ by first-order terms and the spread itself
  converges at about first order. The indicative limits from the three-level
  fits (not a Richardson extrapolation to be quoted) are 85 to 89 N/m for
  all five variants with 4x values, a tight cluster.
- `uwall` and `mesh_inv` move `Q_in` by at most 8% and 4% at 1x and 2x, and
  change the 1x to 2x difference by -10% and -6%.

Conclusion: no single fluid-side treatment limits the convergence of `Q_in`.
The first-order behaviour is shared by all treatments tested (wall pressure,
wall velocity, gradient, mesh motion, GCL), so it is a property of the
discretisation as a whole (a hypothesis: the singular flag tip and the
vortical in-phase load) rather than of one of these options. The cheapest
option that is better than the baseline is `gcl`, and its gain is 25%. The
untested items are the `ddtCorr`/PISO splitting, the convection scheme and
the tip geometry.

**Correction (coordinator review, 2026-10-10): the `gcl` variant is not a
GCL fix, and its result should not be used as evidence.** OpenFOAM v2512
`backwardDdtScheme::meshPhi` already returns
`(1 + c) mesh.phi() - c mesh.phi().oldTime()` with `c = dt/(dt + dt0)`, which is
`1.5 phi^n - 0.5 phi^(n-1)` at constant dt. solids4foam takes the mesh flux
through `fvc::meshPhi` (`pimpleFluid.esi.C` via `fvc::makeRelative`, and
`newMovingWallVelocity`). The baseline therefore already satisfies the
space-conservation law for backward ddt; the standard flux is not an `O(dt)`
GCL violation. `dynamicMotionSolverBDF2FvMesh` overwrites `mesh.phi()` itself
with `1.5 phiE^n - 0.5 phiE^(n-1)`, so `fvc::meshPhi` applies the BDF2
combination a second time. That variant is thus GCL-inconsistent, and its
modest gain (2x to 4x change -25%) is unexplained and possibly an artefact.
The conclusion that all valid treatments share an observed order of about
1.1-1.3 is unchanged. The statements above that `gcl` "helps" and is "the
cheapest option that is better than the baseline" are withdrawn.

Cost: about 132 core-hours of wall time times ranks on xenosim, which was
heavily loaded (1x about 2, 2x about 47 of which 12 are restarts, 4x about 78,
about 4 lost to a run killed by mistake). No 8x run was made.

## Where the mesh dependence of Q_in sits on the flag (tip localisation)

Question: is the slow convergence of `Q_in` (50.3, 72.2, 81.0 N/m at 1x, 2x,
4x, order 1.3) carried by the geometric singularities, the free-end corners
and the flag-cylinder corners? The baseline replays were rerun unchanged
(same restart times, ranks, binary and window as the section above: 1x
restart 25 s to 33.2 s on 1 rank, 2x restart 32 s to 33.2 s on 8, 4x restart
29 s to 30.2 s on 24) with one added coded function object (`plateTraction`,
`scripts/hron_turek_tip.py`) that writes the face traction (pressure and
viscous, N/m, same definition as `plateEnergy`) of every flag face at every
step of the last 5.2 periods. `scripts/tip_localisation.py extract` reduces this to the
first-harmonic phasor of each face over the five-period window of
`energy_analysis.py` (`reference/tip_localisation/faces_{1x,2x,4x}.csv`, 84, 168
and 336 faces), and `analyse` sums the face contribution
`Re[G1y_f (a_f + i b_f)]/|D_tip|` over geometric regions. The restarts
reproduce `Q_in` of the baselines to all digits, and the sum over the regions
reproduces `Q_in`, the lift amplitude and the mean drag of
`energy_analysis.py` (50.3119, 72.1723, 81.0181 N/m; 176.50, 178.99,
163.26 N/m; -13.406, -16.195, -18.643 N/m).

Regions are defined by x on the undeformed flag (x from 0.24899 to 0.6 m,
thickness t = 0.02 m); a face is split between regions in proportion to its
x-extent, so the partition is the same on the three meshes (the 0.1 t bands,
0.002 m, are below the face size at 1x/2x and about one face at 4x, so they
only carry a fraction of one face). `tip_corner_1t` is the end face plus the
last 1 t of both long faces; `cyl_corner_1t` the first 1 t at the
cylinder; `smooth_mid` the rest, split into five equal segments `mid_1` to
`mid_5` (x from 0.269 to 0.58 m).

| Region | Q_in 1x | 2x | 4x | 1x to 2x | 2x to 4x | share of 2x to 4x | order |
|---|---|---|---|---|---|---|---|
| tip end face | 4.27 | 4.61 | 4.63 | +0.34 | +0.02 | 0.3% | 3.8 |
| tip corner, last 1 t (incl. end face) | 2.33 | 6.89 | 8.84 | +4.56 | +1.95 | 22% | 1.2 |
| cylinder corner, first 1 t | 0.015 | 0.013 | 0.007 | -0.002 | -0.006 | -0.1% | non-convergent |
| mid_1 (0.27-0.33) | 1.74 | 1.73 | 1.58 | -0.02 | -0.15 | -1.7% | no |
| mid_2 | 11.16 | 11.49 | 11.52 | +0.33 | +0.03 | 0.4% | 3.3 |
| mid_3 | 18.03 | 19.77 | 20.69 | +1.74 | +0.92 | 10% | 0.9 |
| mid_4 (0.46-0.52) | 5.78 | 9.33 | 11.83 | +3.55 | +2.50 | 28% | 0.5 |
| mid_5 (0.52-0.58) | 11.24 | 22.94 | 26.54 | +11.70 | +3.60 | 41% | 1.7 |
| all smooth faces (mid_1..5) | 47.96 | 65.27 | 72.17 | +17.31 | +6.90 | 78% | 1.3 |
| total | 50.31 | 72.17 | 81.02 | +21.86 | +8.85 | 100% | 1.3 |

(N/m; share is of the 2x to 4x change; order only where the 1x-2x and 2x-4x changes have the same sign
and the ratio exceeds 1. `reference/tip_localisation/tip_localisation_regions.csv` and `.json`
also hold the per-length values, the lift amplitude and the mean drag of
every region, and the nested sub-bands.) The mean drag changes by -1.2 N/m
(tip end face), -2.1 N/m (smooth faces) and -0.35 N/m (tip corner, last 1 t) from 2x to 4x, of a total of
-2.45 N/m.

Reading. The singularity hypothesis is not supported. The cylinder corners
contribute nothing (the traction there is small and the flag is nearly at
rest). The free-end corner carries 22% of the 2x-to-4x change, about its
share of the 1x-to-2x change (21%), and its order (1.2) equals that of the
total, so it does not converge more slowly than the rest; the end face itself
is converged (order 3.8). The smooth long faces carry 78% of the change and converge at
the same order as the total (1.3). The change comes from the aft third of the flag
(`mid_4`, `mid_5`, x of 0.46 to 0.58 m, 69% of the 2x-to-4x change): the
traction profile (`traction_distribution.png`) has a second pressure lobe near
x of 0.56 m, whose amplitude and in-phase part keep growing with refinement
(Q_in density at the lobe peak, upper face, 140, 248, 283 N/m per m), with a minimum near
x of 0.50 m that is poorly resolved at 1x. This is a smooth flow feature of the aft flag, not a corner singularity. Per unit length the tip corner
is the most sensitive region together with `mid_5`, so a tip contribution of
the expected kind exists, but it is a minority of the total. Within the last
0.5 t the cumulative Q_in (`tip_cumulative.png`) is 2.5, 4.8, 6.0 N/m, so the
refinement change is small in absolute terms there.
Caveats: the geometric bands below one cell width (0.1 t) are not resolved by any of the
meshes, so nothing is concluded about the corner at that scale; three
levels with orders of 0.5 to 1.7 in the segments (not asymptotic) mean orders
are indicative.

Run: `hron_turek_tip.py` (case set-up), `run_replay_variant.sh <case>` (as
above), `tip_localisation.py extract|analyse`; run directories on xenosim
`~/ht_tip/work/t_{1x,2x,4x}`, about 19 core-hours in total. Plots:
`reference/tip_localisation/traction_distribution.png` (amplitude, phase and
Q_in density of the first-harmonic traction along the upper and lower faces) and
`tip_cumulative.png` (cumulative Q_in from the tip).
