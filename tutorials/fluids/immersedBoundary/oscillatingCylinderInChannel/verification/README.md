# oscillatingCylinderInChannel verification

Opt-in verification of the four forcing methods of `immersedBoundaryForce`
for a body that moves and accelerates through a stationary mesh:

- `cutLink`;
- `ghostCell` (documented for static bodies; its moving-body results are
  informative, not a general failure of the method);
- `penalty`;
- legacy incremental, `method incremental`.

The case is the cylinder (D = 0.1 m) oscillating in a closed channel of
2.2 x 0.41 m filled with fluid initially at rest, with its centre at
`x_c(t) = 1.1 + 0.25 sin(2 pi t/4)` m, `y_c = 0.2` m. The results are
compared with two moving body-fitted mesh solutions and, as a secondary
comparison, with the result of Wan and Turek (2007) digitised from a
figure of Erzincanli and Sahin (2013). The verification is not run by
`tutorials/Alltest` or `tutorials/Alltest-regression`.

## What it shows

Accurate steady drag does not by itself establish accurate unsteady forces.
For a moving body, the force estimator, the surface quadrature, the mesh
and time resolution, the time scheme and the crossing of cells by the body
all affect the drag and, more strongly, the lift.

The largest effect is a limitation of the greedy quadrature, the only
quadrature of the tested revision. The cutLink surface traction with it
gives an accurate drag but a lift that differs from body-fitted reference B
by 140 %, 149 % and 140 % on D/h = 10, 20 and 40, because the thinning of
the surface points in STL-facet order breaks the cancellation of the
transverse pressure force. On D/h = 80 every candidate point is kept and
the lift agrees to 4.6 %; this does not establish that the greedy
quadrature converges in general. This defect is a property of that
quadrature, not of the cutLink method as such; solids4foam pull request 612
changes the default quadrature and is not evaluated here.

Three force estimates are kept apart: the cutLink greedy surface traction;
the cutLink momentum exchange of the same runs; and, for penalty, the
forcing plus the inertia of the fluid inside the body. The differences from
reference B are read against the difference between the two references, A
and B (0.77 % in Cd, 3.0 % in Cl), and smaller differences are not ranked.
By that rule cutLink and penalty are not distinguished in drag at
D/h = 40 and 80 (1.10 and 0.65 % against 1.41 and 1.08 %), and penalty and
the cutLink momentum exchange are not distinguished in lift (3.7 and 3.5 %
against 5.3 and 4.8 %). Both are clearly closer to B than ghostCell and the
legacy incremental method in drag, and penalty and the momentum exchange in
lift. Under the tested configuration, penalty combines a drag and a lift
close to the reference spread with the lowest cost of the four methods at
the common D/h = 40 (1.493 allocated CPU hours); this is a result for this
case and these settings, not a general ranking. ghostCell is documented for
static bodies: its moving-body results are informative and not a general
failure of the method. `results.md` has the evidence, the figures and the
limitations.

| Figure (`images/`) | Content |
| --- | --- |
| `force_histories.png` | four methods vs B; penalty vs B and Wan-Turek |
| `error_vs_resolution.png` | raw Cd and Cl errors against D/h |
| `cutlink_estimators.png` | cutLink lift: greedy traction, momentum exchange |
| `error_vs_cost.png` | further diagnostic |
| `cell_crossing.png` | diagnostic; not regenerable from tracked files |

The files:

| File | Content |
| --- | --- |
| `Allverify` | runs the verification (wrapper of the driver) |
| `Allclean` | removes `work/` and `postProcessing/` |
| `results.md` | settings, results, limitations |
| `reference/` | reference histories and the compact results |
| `scripts/oscillating_cylinder_verification.py` | the driver |
| `scripts/plot_verification.py` | figures from the compact results |
| `scripts/test_oscillating_cylinder_verification.py` | unit tests |
| `images/` | previews of the figures (PNG, 110 dpi) |

The cases are generated in `work/` and the results are written to
`postProcessing/`. Neither is tracked, and `./Allclean all` in the tutorial
removes them.

## Running

The driver needs Python 3 with numpy. scipy is optional: it is needed only
for the filtered diagnostics (the 5 Hz low-pass errors and the force above
5 Hz), which are nan without it, and the primary raw-force comparison does
not use it. matplotlib is needed for the figures (no LaTeX installation).
Source OpenFOAM and solids4foam, then:

```bash
./Allverify                      # cutLink, Level 1: run, analyse, check
./Allverify --runs smoke         # the first 0.5 s of it (a few seconds)
./Allverify dry-run --runs R1A R1B R1C R1D
./Allverify --runs R1A R1B R1C R1D   # the four methods on Level 1
./Allverify status               # the cases of the manifest
./Allverify --help
```

Each run of the campaign is a copy of the tutorial in
`work/cases/<method>-L<level>-dt<deltaT>`; `results.md` lists the runs
(`R1A` ... `R5D`, `N5`-`N20`, `O4`). Each copy has:

- the same mesh, geometry, motion, fluid and boundary conditions;
- a fixed time step, halved with the cell size (`--ranks` to change the
  decomposition);
- the method settings below, written explicitly. They are the settings of
  the static-cylinder verification (solids4foam pull request #595), not
  all the defaults of `immersedBoundaryForce`.

Serial cases run with the tutorial's `Allrun`, parallel cases with
`solids4Foam::runParallel`. A manifest, `work/manifest.json`, records each
case and its job, so that a campaign can be resumed with `status`.

The method settings are:

- `cutLink`: `penaltyCoeff 1e3`, `pinnedRateCoeff 100`, `linkCorrection yes`,
  `apertureCoupling yes`, `forceEstimator surfaceTraction`;
- `ghostCell`: `penaltyCoeff 1e3`, `imageDistance 1.5`,
  `imageInterpolation cellPoint`;
- `penalty`: `penaltyCoeff 1e3`, `weighting volumeFraction`,
  `surfaceRateCoeff 3`, `occupancy signedDistance`;
- incremental: `couplingCoeff 0.2`, `occupancy vertexFraction`,
  `nOuterCorrectors 5`. `couplingCoeff 0.2` is the value chosen for
  consistency with #595; it is not the source default, which is 0.8.

The other methods use `nOuterCorrectors 1`. All use `nCorrectors 3` and the
tutorial's schemes (Euler, linearUpwind).

The other commands:

```bash
./Allverify generate --runs ...  # generate the cases only
./Allverify analyse              # postProcessing/cases.csv and histories
./Allverify check                # compare with reference/results/cases.csv
./Allverify export               # compact results to reference/results/
./Allverify plot                 # figures from reference/results/
./Allverify submit --runs ...    # optional: one SLURM job per case
```

- `analyse` checks each case for its numerical integrity and reports PASS
  or FAIL: normal end, complete force history, finite forces, Courant number
  at most 1.
- `check` (also run by `./Allverify`) compares each analysed case with the
  tracked result of the same case. It checks the normalised RMS differences
  from references A and B and the fluctuation level, within 5 % relative or
  0.002 absolute, and reports PASS or FAIL. The fluctuation level
  (`Cd_hf_rms_above_5Hz`) needs scipy: without it that comparison is
  reported as SKIP and counted as skipped, and the raw-force comparisons are
  made. A nonfinite value of any other metric, or of that one with scipy
  installed, is a FAIL. It also prints the OpenFOAM version and solids4foam
  build recorded for the analysed cases next to the tracked baseline.
  These tolerances detect a change of a method, its settings or the force
  estimator, beyond the differences between compilers, MPI libraries and
  decompositions (about 1e-5). They are reproduction checks, not accuracy
  targets.
- `submit` is optional: `./Allverify` runs the cases locally, or in a batch
  job of your own. Its options are `--env-script` (default
  `$S4F_ENV_SCRIPT`), `--partition`, `--max-cpus`, `--time`, `--mem` and
  `--sbatch-args`.

The tests need no OpenFOAM, SLURM (`sbatch` and `squeue` need not be
installed) or scipy; the tests that need scipy are skipped without it. They
check the driver with synthetic data. They also check that the tracked
references, `reference/results/` and the tables of `results.md` agree, and
that the figures are regenerated from `reference/results/`; this last test
is skipped without matplotlib and scipy:

```bash
python3 -m unittest discover -s scripts -p 'test_*.py'
```

The previews in `images/` are written by
`python3 scripts/plot_verification.py --out images --dpi 110 --png-only`;
without `--png-only`, PDF and SVG versions are also written. The force
histories in the figures are drawn from the compact histories of
`reference/results/` (every 0.02 s, for display only); the errors are
computed from the complete histories. `cell_crossing.png` was drawn from
full-resolution extracts of the finest runs, which are not tracked (see
Repository size); `plot` skips it without them.

## Quadrature of the cutLink surface traction

The tracked results use the quadrature of the revision that produced them:
the greedy thinning of the candidate surface points, the only quadrature of
that revision. The cutLink cases always contain `tractionQuadrature greedy;`.
The driver does not read the solver source to decide this. The entry is
explicit, so a solver whose default has changed (solids4foam pull request
612 makes `clustered` the default) still runs the historical quadrature,
and the cases do not depend on the solver default.

Limitations by solver version:

- a solver without the `tractionQuadrature` keyword (the tested revision,
  and `development` at the base of this pull request) ignores the entry and
  uses its only quadrature, greedy;
- a solver with the keyword runs greedy as requested, and stops with an error
  on a name it does not know, so it never silently falls back to another
  quadrature for the historical cases.

`--traction-quadrature none|clustered` writes another quadrature, for a
solver that has the keyword. Those cases get a tagged name (for example
`cutLink-L1-dt0p0025-clustered`) and are not checked against the tracked
results. On a solver without the keyword they would run greedy under a
tagged name, so use them only with a solver that has it. The other methods
do not evaluate the surface traction and get no quadrature entry.

## Force definitions

The force file of `immersedBoundaryForce`,
`postProcessing/immersedBoundary/<time>/cylinder.dat`, is parsed by its
header:

- `cutLink`: columns 2-4 are the selected estimator, here
  `surfaceTraction`, the pressure and viscous traction integrated over the
  surface of the cylinder, which is the whole hydrodynamic force; columns
  8-10 are zero; columns 11-13 are the alternative momentum-exchange
  estimate, which accounts for the change of the momentum of the penalised
  cells. The two estimates are reported separately, never added.
- `ghostCell`, `penalty`, legacy incremental: columns 2-4 are the force of
  the forcing and columns 8-10 the inertia of the fluid inside the cylinder;
  their sum is the hydrodynamic force.

`C = 2F/(rho U_ref^2 D L_z)`, with `U_ref = 2 pi A/T`, `rho`, `D` and `L_z`
read from the case.

## Tracked results

`reference/results/` was exported by `./Allverify export` after
`./Allverify analyse` of the 24 runs of `results.md`. The Level 4 runs and
the ghostCell, Level 3 incremental and outer-corrector runs were added to
the earlier runs of the same revision; the earlier results are unchanged.

### Revision and provenance

- **Historical baseline:** all 24 runs used the solids4foam build
  `f33b244a2` with OpenFOAM v2412. This is stated from the build record of
  those runs; the runs themselves wrote no per-case provenance, and none is
  added to the tracked files afterwards.
- **Revision anchor:** the immersed-boundary, fluid-model and tutorial
  sources of `f33b244a2` are identical to those of `development` at
  `8cf077860`, the base of this pull request (checked with `git diff`).
  `8cf077860` is the only later revision named in these documents. Later
  commits of `development` were not used for any result.
- **Newly generated cases** record, in `case.json` (`provenance`) and in
  the columns of the analysed `postProcessing/cases.csv`:
  `openfoam_version` (`WM_PROJECT_VERSION`), `openfoam_variant`,
  `solver_build` (the `SHA=` line of the `BUILD_INFO.txt` in the root of the
  build of the `solids4Foam` executable found on `PATH`, if there is one),
  `solver_binary_sha256` (the executable itself), and, separately,
  `tutorial_revision`, the revision of the tutorial tree that supplied the
  case (`git_sha`, as before). An entry that cannot be determined is empty;
  nothing is guessed. `case.json` also records the tracked baseline.
- **OpenFOAM v2512 reproduction of R1A:** `./Allverify` (cutLink Level 1,
  greedy quadrature, one rank, then `check`) was run once on OpenFOAM
  v2512 with a clean single-core build of the baseline source `f33b244a2`
  (gcc 11.4.0, Open MPI 4.1.2, no PETSc; the source compiled without
  change). The run ended normally, passed the integrity checks, and the
  five compared metrics (`Cd_nrms_B`, `Cl_nrms_B`, `Cd_nrms_A`,
  `Cl_nrms_A`, `Cd_hf_rms_above_5Hz`) equal the tracked v2412 values to
  all stored digits (5 of 5 PASS, tolerances unchanged); 102 s wall time.
  This covers one case on one platform: Levels 2-4, the other methods and
  parallel runs were not repeated on v2512. The `mesh_hash` of the v2512
  case differs from the tracked one; it hashes the mesh files including
  their headers, is not used by `check`, and was not investigated further.

### Repository size

`reference/results/` holds `cases.csv` (rows of all 24 runs) and the 0.02 s
histories of the 8 s runs that the figures draw, 19 files; `check` reads
only `cases.csv` and references A and B. The histories of the serial repeat
O4 (identical metrics to R2A), of the three outer-corrector runs and of the
smoke run, and the four full-resolution extracts of `cell_crossing`, are not
tracked (0.40 MB of 1.00 MB removed); the rows of those runs remain in
`cases.csv`. `./Allverify export --full-resolution` writes the extracts
again, after rerunning the finest cases of each method. `cell_crossing.png`
is therefore a preview of the original runs: it is not regenerated from
the tracked files, and `plot` writes the other four figures only.

## References

- Wan, D., Turek, S. (2007). Fictitious boundary and moving mesh methods for
  the numerical simulation of rigid particulate flows. Journal of
  Computational Physics, 222(1), 28-56.
  <https://doi.org/10.1016/j.jcp.2006.06.002>
- Erzincanli, B., Sahin, M. (2013). An arbitrary Lagrangian-Eulerian
  formulation for solving moving boundary problems with large
  displacements and rotations. Journal of Computational Physics, 255,
  660-679. <https://doi.org/10.1016/j.jcp.2013.08.038>
