# Verification of the immersed Taylor-Couette flow

This directory verifies `immersedTaylorCouette` against the **exact** solution
of the circular Couette flow: the azimuthal velocity in the gap and the torque
on the rotating rotor, for the four immersed-boundary methods of
`immersedBoundaryForce`: `cutLink` (with the momentum exchange, the surface
traction and the `forcing` estimate as three estimators, each from its own run),
`penalty`, the legacy `incremental` method and (exploratory) `ghostCell`.
Nothing in the tutorial or in `src/` is changed: every case is a copy of the
tutorial with the method settings written into `constant/fvOptions`.

Three kinds of statement are kept apart in `results.md`:

- **exact**: comparisons with the closed-form solution (`results.md`, Exact
  solution), whose torque is derived and tested independently of the README;
- **observed**: the cost, the sensitivity to the time step, the behaviour of the
  torque estimators, taken from single runs;
- **exploratory**: everything about `ghostCell`, which is documented for static
  bodies (here the rotor turns about its own axis, so its position is invariant;
  see the applicability audit in `results.md`).

No result of another solver is used here.

## Quick start

```bash
cd tutorials/fluids/immersedBoundary/immersedTaylorCouette/verification
./Allverify unit            # exact tests, no OpenFOAM needed (seconds)
./Allverify dry-run --runs matrix    # the cases, ranks and observed cost
./Allverify                 # smoke: cutLink, Level 1 (about 20 s), then check
./Allverify check           # PASS / FAIL / SKIPPED
./Allverify tables          # the tables of results.md from the tracked data
```

`./Allverify` needs OpenFOAM.com (v2412 was used for the tracked results) with
`solids4Foam` and `libimmersedBoundary` in the path. Without them the run is
reported as `SKIPPED`, never as `PASS`.

<!-- markdownlint-disable MD013 -->

| Command | What it does |
| ------- | ------------ |
| `unit` | the unit tests of the exact solution, the torque convention, the error analysis and the driver (no solver) |
| `dry-run --runs R` | lists the cases of the named runs with their ranks; starts nothing |
| `generate --runs R` | writes the cases to `work/cases/` and records the SHA-256 of their inputs in `work/manifest.json` |
| `run --runs R` | generate and run (completed cases are kept unless `--force`); the default run is `smoke`; follow with `analyse` and `check` |
| `status` | state of the cases in the manifest |
| `analyse` | the result rows of the completed cases, `work/results.csv` |
| `check` | the exact-solution checks, the consistency of the tracked results, and the completed cases against the tracked rows |
| `tables` | the tables of `results.md` from `reference/results/cases.csv` |

<!-- markdownlint-enable MD013 -->

Named runs: `smoke` (Level 1), `L1` to `L4` (the variants of a level; Level 4
runs on four ranks), `matrix` (Levels 1-3, the tracked matrix) and `dt` (the
time-step control, 40 cases). `--methods`, `--levels` and `--ranks` select
within them. Levels 1-3 of the matrix are about 1.1 serial hours in all (5 s to
15 minutes per case); a Level 4 case is 0.6 to 2.5 rank-hours on four ranks; the
time-step run is about 9 serial hours. All are opt-in and none is part of the CI
of the repository.

Parallel runs are decomposed in slabs along y. The pressure reference point of
the tutorial, `(0 0.075 0.05)`, lies on the plane x = 0, so a decomposition
along x (the default of many cases) puts it on a processor boundary and
`setRefCell` stops with "found on 2 domains".

## What is computed

For each case:

- the **rotor torque** about the axis, the mean over the window 8 <= t <= 10 s
  (time-weighted: the step is adaptive), the standard deviation in the window,
  and the difference from the exact torque, (T - T_exact)/T_exact, so that a
  positive value is a torque that is too large in magnitude;
- the **steady state**: the relative change of the torque over the window and
  the relative L2 change of the velocity between t = 5 and 10 s, both below 1e-3
  (predefined). A case that is not steady is reported with its mean and variance
  and never as converged;
- the **velocity** at t = 10 s in the fluid cells (those with a centre in R1 < r
  < R2): volume-weighted L1, L2 and Linf norms of the azimuthal difference from
  the exact velocity at the cell centres, and of the radial velocity (exact
  value zero), relative to omega R1; also for the cells at least two cell widths
  from both walls (`far`) and for the others (`iface`);
- an independent torque from the velocity field alone (the viscous torque
  transmitted across circles inside the gap), which in a steady flow equals the
  torque on the rotor and uses no force output of the solver;
- the torque on the stator (reported, not used: see `results.md`).

## Settings

Those of the stationary, oscillating and bending verifications, written
explicitly into every case:

<!-- markdownlint-disable MD013 -->

| Method | Settings |
| ------ | -------- |
| `cutLink` | tutorial; `forceEstimator surfaceTraction` or `forcing` for the other two estimators |
| `penalty` | `penaltyCoeff 1e3`, `weighting volumeFraction`, `surfaceRateCoeff 3`, `occupancy signedDistance` |
| `incremental` | `couplingCoeff 0.2`, `occupancy vertexFraction`, `nOuterCorrectors 5` |
| `ghostCell` | `penaltyCoeff 1e3`, `imageDistance 1.5`, `imageInterpolation cellPoint` |

<!-- markdownlint-enable MD013 -->

The time step of the tutorial is limited by a Courant number of 0.5; because the
cells halve with each level, so does the step. The `dt` run adds `maxDeltaT`
(0.01, 0.005 and 0.0025 s) to separate the time step from the mesh.

## Tracked data

- `reference/results/cases.csv`: one row per method, estimator, level and time
  step cap, with the quantities above, the wall time and the allocated
  rank-hours of the run that produced it;
- `reference/results/history-<variant>-L<k>.csv`: the torque on the rotor at
  every time step of the tutorial runs;
- `reference/results/profiles-L<k>.csv`: azimuthal averages of the velocity in
  radial bins of half a cell width, with the exact velocity;
- `reference/results/scaling.csv`: wall time of the Level 4 `penalty` case on 1,
  2, 4 and 8 ranks over the first 1.5 s;
- `images/`: previews written by `scripts/plot_verification.py`; the manuscript
  figures (PDF) are written by the same script without `--png-only`.
- `reference/results/controls/`: the control runs (step caps, 20 s runs, the
  surface-traction audit) and the tables derived from them by
  `scripts/taylor_couette_controls.py`; their discussion is in
  [`controls.md`](controls.md), regenerated by `scripts/make_controls_doc.py`,
  and the figures by `scripts/plot_controls.py`. The tests of the controls run with
  `python3 -m unittest discover -s scripts` (`./Allverify unit` runs the driver
  tests only).

The histories, profiles and scaling files are written by research export scripts
that call the analysis functions of the driver on the complete runs; `analyse`
writes only `work/results.csv`. `./Allverify tables` regenerates the tables of
`results.md`, and `python3 scripts/plot_verification.py --out images --dpi 110
--png-only` regenerates the previews, from the tracked data only.

## Pre-registered criteria

Fixed before the matrix was analysed (the steady criterion compares the end
values of the window, which is noisy for oscillating cases; it was kept as
defined and the flag is not compared by `check` when it is near the threshold):

- the comparison window is 8 <= t <= 10 s for every method; the steady criteria
  are 1e-3 for the drift of the torque and for the change of the velocity;
- `check` tolerances for a re-run against the tracked row of the same case:
  0.005 (absolute difference of the torque difference) and 1e-3 (absolute
  difference of the velocity L2 norm); the closed form of the torque is checked
  against the stress integral to 1e-9 relative, and the torque of the velocity
  field against the force output to 0.05 (Levels 2-4, all but the legacy
  `incremental` method, whose difference is reported).

## Limitations

See `results.md`, Limitations. In short: single runs on a shared node (the costs
are observations, not performance models); no formal convergence order is
claimed from the Level sequences, whose time step also changes; `ghostCell` is
exploratory; the torque of `penalty`, `incremental` and `ghostCell` is the
forcing torque, which equals the torque on the rotor only in the steady state
used here.
