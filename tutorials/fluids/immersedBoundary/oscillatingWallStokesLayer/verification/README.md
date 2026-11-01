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
./Allverify                 # smoke: cutLink, Level 1 (about 10 s), then check
./Allverify check           # PASS / FAIL / SKIPPED
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
| `run --runs R` | generate and run (completed cases are kept unless `--force`); the default run is `smoke`; follow with `analyse` and `check` |
| `status` | state of the cases in the manifest |
| `analyse` | the result rows of the completed cases, `work/results.csv` |
| `check` | the exact-solution checks, the consistency of the tracked results, and the completed cases against the tracked rows |
| `tables` | the tables of `results.md` from `reference/results/cases.csv` |

<!-- markdownlint-enable MD013 -->

Named runs: `smoke` (Level 1), `L1` to `L4` (all variants of a level), `matrix`
(Levels 1-4, the tracked matrix) and `dt` (the time-step control, 36 cases).
`--methods`, `--levels` and `--ranks` select within them. The matrix is about 12
serial minutes on one core of a recent node (the observed wall times are in
`reference/results/cases.csv`), and the time-step run about an hour; both are
opt-in and neither is part of the CI of the repository.

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
  the wall placed on a mesh face) and the tables derived from them by
  `scripts/stokes_controls.py`; their discussion is in
  [`controls.md`](controls.md), regenerated by `scripts/make_controls_doc.py`,
  and the figures by `scripts/plot_controls.py`. The tests of the controls run with
  `python3 -m unittest discover -s scripts` (`./Allverify unit` runs the driver
  tests only).

The shear, profile and scaling files are written by research export scripts that
call the analysis functions of the driver on the complete runs; `analyse` writes
only `work/results.csv`. `./Allverify tables` regenerates the tables of
`results.md`, and `python3 scripts/plot_verification.py --out images --dpi 110
--png-only` regenerates the previews, from the tracked data only.

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
