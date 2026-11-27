# Verification of the oscillating-wall Stokes layer

This directory verifies `oscillatingWallStokesLayer` against the **exact**
solution of the Stokes second problem, for the four immersed-boundary methods of
`immersedBoundaryForce`: `cutLink` (with the momentum exchange and the surface
traction as two estimators of the same run, and the `forcing` estimator as a
third), `penalty`, the legacy `incremental` method and (exploratory)
`ghostCell`. Nothing in the tutorial or in `src/` is changed: every case is a
copy of the tutorial with the method settings written into `constant/fvOptions`.

Three kinds of statement are kept apart in `results.md`:

- **exact**: comparisons with the closed-form solutions of
  [`results.md`](results.md#exact-solutions) (the periodic solution, the
  periodic solution with the lid of the mesh, and the exact solution of the
  problem that is simulated, which starts from rest);
- **observed**: the cost, the sensitivity to the time step, the behaviour of the
  force estimators, taken from single runs;
- **exploratory**: everything about `ghostCell`, which is documented for static
  bodies (here the slab moves in its own plane, so its position is invariant;
  see the applicability audit in `results.md`).

No result of a body-fitted solver is used here.

## Quick start

```bash
cd tutorials/fluids/immersedBoundary/oscillatingWallStokesLayer/verification
./Allverify unit            # exact tests, no OpenFOAM needed (seconds)
./Allverify dry-run --runs matrix    # the cases, ranks and observed cost
./Allverify                 # smoke: cutLink, Level 1 (about 10 s); runs only
./Allverify analyse         # the result rows of the completed cases
./Allverify check           # PASS / FAIL / SKIPPED (analyses, compares)
./Allverify tables          # the tables of results.md from the tracked data
```

`./Allverify` needs OpenFOAM.com (v2412 was used for the tracked results) with
`solids4Foam` and `libimmersedBoundary` in the path. Without them the run is
reported as `SKIPPED`, never as `PASS`.

<!-- markdownlint-disable MD013 -->

| Command | What it does |
| ------- | ------------ |
| `unit` | the unit tests of the exact solutions, the force and shear conventions, the error analysis and the driver (no solver) |
| `dry-run --runs R` | lists the cases of the named runs with their ranks; starts nothing |
| `generate --runs R` | writes the cases to `work/cases/` and records the SHA-256 of their inputs in `work/manifest.json` |
| `run --runs R` | generate and run (completed cases are kept unless `--force`); the default run is `smoke`; the subcommand `run` may be left out (`./Allverify --runs L1 L2`); it does not analyse or check, follow with `analyse` and `check` |
| `status` | state of the cases in the manifest |
| `analyse` | the result rows of the completed cases, `work/results.csv` |
| `check` | the exact-solution checks, the consistency of the tracked results, and the completed cases against the tracked rows |
| `tables` | the tables of `results.md` from `reference/results/cases.csv` |

<!-- markdownlint-enable MD013 -->

Named runs: `smoke` (Level 1), `L1` to `L4` (all variants of a level), `matrix`
(Levels 1-4, the tracked matrix) and `dt` (the time-step control, 36 cases).
`--methods`, `--levels` and `--ranks` select within them. The control runs of
`reference/results/controls/` have named runs too (see Controls below). The
matrix is about 12 serial minutes on one core of a recent node (the observed
wall times are in `reference/results/cases.csv`), and the time-step run about an
hour; both are opt-in and neither is part of the CI of the repository.

## What is computed

For each case, over the fourth period (3 < t <= 4 s):

- the **wall shear stress**, the force on the slab in x divided by the wall area
  (0.05 m x 0.1 m) for rho = 1, from the estimators below; amplitude, phase, RMS
  of the time history and the fluctuation about its mean and first three
  harmonics;
- the **velocity** u_x in the fluid cells (centre above the wall), as
  volume-weighted L1, L2 and Linf norms relative to U0, over the eight snapshots
  of the period, and separately in the cells within one Stokes thickness of the
  wall (`near`) and in the cells of the first row above the wall (`iface`).

Each is compared with two exact solutions: the periodic solution of the
tutorial, and the exact solution of the simulated initial-value problem (the
columns with the suffix `_ivp`; the columns with the suffix `_al` are the
diagnostic with the inertia averaged over two steps, `results.md`, Alignment),
because the start-up transient has not decayed in the fourth period
(`results.md`, Exact solutions).

The estimator of the shear is `momentumExchange` and `surfaceTraction` (two
columns of the same cutLink run), `forcing` (a separate run; the tutorial's
`forceEstimator forcing`) for cutLink, and the forcing plus the rate of change
of the momentum of the fluid inside the slab (columns 2 and 8 of the force file,
as in the tutorial's plot) for `penalty`, `incremental` and `ghostCell`.

## Settings

Those of the stationary, oscillating and bending verifications, written
explicitly into every case:

<!-- markdownlint-disable MD013 -->

| Method | Settings |
| ------ | -------- |
| `cutLink` | tutorial (`writeSurfaceTraction yes`); the `forcing` run adds `forceEstimator forcing` |
| `penalty` | `penaltyCoeff 1e3`, `weighting volumeFraction`, `surfaceRateCoeff 3`, `occupancy signedDistance` |
| `incremental` | `couplingCoeff 0.2`, `occupancy vertexFraction`, `nOuterCorrectors 5` |
| `ghostCell` | `penaltyCoeff 1e3`, `imageDistance 1.5`, `imageInterpolation cellPoint` |

<!-- markdownlint-enable MD013 -->

The time step is the 0.001 s of the tutorial (`backward`), also for the controls
except where stated: the `dt` run changes only `deltaT`.

## Tracked data

- `reference/results/cases.csv`: one row per method, estimator, level and time
  step, with the quantities above, the volume of the body in the mesh, the wall
  time and the allocated rank-hours of the run that produced it;
- `reference/results/shear-L<k>.csv`: the shear of every estimator at every time
  step of the fourth period, with the two exact solutions;
- `reference/results/profiles-L<k>.csv`: the row-averaged velocity at t/T = 0,
  0.25, 0.5, 0.75 of the fourth period, with the two exact solutions;
- `reference/results/scaling.csv`: wall time of the Level 4 `penalty` case on 1,
  2, 4 and 8 ranks over the first second;
- `images/`: previews written by `scripts/plot_verification.py`; the manuscript
  figures (PDF) are written by the same script without `--png-only`.
- `reference/results/controls/`: the control runs (time steps down to 3.1e-5 s,
  the wall placed on a mesh face, the tables derived from them by
  `scripts/stokes_controls.py`); their discussion is in
  [`controls.md`](controls.md), regenerated by `scripts/make_controls_doc.py`,
  and the figures by `scripts/plot_controls.py`. The tests of the controls run
  with `python3 -m unittest discover -s scripts` (`./Allverify unit` runs the
  driver tests only).

`analyse` writes only `work/results.csv`. `./Allverify tables` regenerates the
tables of `results.md`, and `python3 scripts/plot_verification.py --out images
--dpi 110 --png-only` regenerates the previews, from the tracked data only.

### Where every tracked file comes from

<!-- markdownlint-disable MD013 -->

| Tracked file | Generator | Runs it needs |
| ------------ | --------- | ------------- |
| `cases.csv`, `shear-L*.csv`, `profiles-L*.csv` | `scripts/export_results.py cases` | `matrix`, `dt` |
| `controls/cases-timestep.csv` | `scripts/export_results.py cases` | `dt-L1`, `dt-fine`, `dt-6e-5`, `dt-3e-5` |
| `controls/shear-yw0125-L*.csv` | `scripts/export_results.py wall` | `wall-face` |
| `controls/velocity-regions.csv` | `scripts/velocity_regions.py` | `matrix`, `dt-fine` |
| `controls/temporal-convergence.csv`, `offset-audit.csv`, `cell-average.csv` | `scripts/stokes_controls.py all` | none (tracked data only) |
| `scaling.csv` | `scripts/export_results.py scaling` | Level 4 `penalty`, `--end-time 1`, on 1, 2, 4, 8 ranks |
| `controls.md`, `images/` | `make_controls_doc.py`, `plot_controls.py`, `plot_verification.py` | none (tracked data only) |

<!-- markdownlint-enable MD013 -->

The export scripts read the completed cases of one or more work directories
(`--work-dir`, so that independent jobs do not share a manifest) and call the
analysis functions of the driver; they write to `work/export` unless `--out` is
given. What does not regenerate exactly: the column `rank_h` of the tracked rows
is the allocated rank-hours of the scheduler job, the exporter writes ranks x
the wall time of `Allrun.verify`, and wall times depend on the node;
`scaling.csv` is a timing and is not expected to repeat to the second. The
eight-period control (`T8`) has no tracked rows beyond the statement in
`results.md`: `check` compares its fourth period with the tracked row and its
eighth period with its own fourth.

## Controls

Every control of `controls.md` is a named run or an option of the driver, with
the settings of the tutorial otherwise:

<!-- markdownlint-disable MD013 -->

| Control | Command | Notes |
| ------- | ------- | ----- |
| time step 2.5e-3 to 2.5e-4 s, Level 1 | `./Allverify run --runs dt-L1` | 15 cases |
| time step 1.25e-4 s, Levels 1-4 | `./Allverify run --runs dt-fine` | 20 cases |
| time step 6.25e-5 s | `./Allverify run --runs dt-6e-5` | `penalty`, `incremental`, `ghostCell` on Level 3; all variants on Level 4 |
| time step 3.125e-5 s | `./Allverify run --runs dt-3e-5` | `incremental`, Levels 3-4 (about 0.5 and 1.8 rank-h) |
| the wall on a mesh face | `./Allverify run --runs wall-face` | `cutLink` momentum exchange and forcing, Levels 2-4; `--wall-y 0.125` moves only the top vertices of `slab.stl` to y = 0.125 m (the wall of the tutorial is 0.1234 m); the exact solutions use the moved wall |
| eight periods | `./Allverify run --runs T8` | `cutLink` momentum exchange and `penalty`, Level 3; `--end-time 8` is an option of any run (the cases are tagged `-T8`) |
| all of the above | `./Allverify run --runs controls` | about 7 rank-h in total (`dry-run --runs controls` lists the observed cost of each case) |

<!-- markdownlint-enable MD013 -->

Cases with a changed step, end time, wall position or rank count carry it in
their name (`penalty-L4-dt0.00025`, `cutLink-ME-L3-yw0.125`, `penalty-L3-T8`,
`penalty-L4-T1-np2`). `check` compares a time-step case with the tracked row of
the same method, estimator, level and step; a wall case with the tracked shear
histories of `controls/shear-yw0125-L*.csv` (tolerance `CHECK_SERIES`, 2e-7 Pa,
twice the rounding of the files); an end-time case as described above.

## Pre-registered criteria

Fixed before the matrix was analysed and not changed afterwards:

- the comparison window is the fourth period, 3 < t <= 4 s, for every method;
- the noise of the shear is the RMS of the residual of the least-squares fit of
  the mean and the first three harmonics, relative to the RMS of the exact shear
  (nothing is filtered or shifted in any reported result);
- `check` tolerances for a re-run against the tracked row of the same case:
  `CHECK_SHEAR_AMP` = 0.005 (absolute difference of the amplitude error) and
  `CHECK_VELOCITY` = 2e-4 (absolute difference of the velocity L2 norm); the
  closed forms are checked to 1e-9 of U0 (shear, velocity) and 1e-8 Pa
  (amplitude).

## Limitations

See `results.md`, Limitations. In short: single runs on a shared node (the costs
are observations, not performance models), a periodic-target and an
initial-value-target comparison that differ by 0.43 % of the shear (0.012 of U0
in velocity), no formal convergence order is claimed from the Level sequences
(the time step is not refined with the mesh), `ghostCell` is exploratory, and
the shear of `penalty`, `incremental` and `ghostCell` is a small difference of
two large terms.
