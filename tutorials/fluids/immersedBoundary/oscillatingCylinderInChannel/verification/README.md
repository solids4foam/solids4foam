# oscillatingCylinderInChannel verification

Opt-in verification of the four forcing methods of `immersedBoundaryForce`
for a body that moves and accelerates through a stationary mesh:

- `cutLink`;
- `ghostCell` (documented for static bodies; included to show its
  moving-body behaviour);
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
all affect the drag and, more strongly, the lift. In particular, the
cutLink surface traction with the greedy quadrature (the only quadrature
of the tested revision) gives an accurate drag but a lift that differs
from body-fitted reference B by 140 %, 149 % and 140 % on D/h = 10, 20 and
40, because the thinning of the surface points in STL-facet order breaks
the cancellation of the transverse pressure force. On D/h = 80 every
candidate point is kept and the lift agrees to 4.6 %; this does not
establish that the greedy quadrature converges in general. At the common
D/h = 40, `penalty` gives both an accurate drag (1.41 %) and an accurate
lift (3.7 %) at the lowest cost of the four methods; this is a result for
this case and these settings, not a general ranking. `results.md` has the
evidence, the figures and the limitations.

| Figure (`images/`) | Content |
| --- | --- |
| `force_histories.png` | four methods vs B; penalty vs B and Wan-Turek |
| `error_vs_resolution.png` | raw Cd and Cl errors against D/h |
| `cutlink_estimators.png` | cutLink lift: greedy traction, momentum exchange |
| `error_vs_cost.png`, `cell_crossing.png` | further diagnostics |

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

The driver needs Python 3 with numpy. scipy is needed for the filtered
diagnostics, which are reported as nan without it, and matplotlib for the
figures (no LaTeX installation is needed). Source OpenFOAM and solids4foam, then:

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
- the method settings below, written explicitly so that a change of a
  solver default does not change the cases.

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
  `nOuterCorrectors 5`.

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
  0.002 absolute, and reports PASS or FAIL. These tolerances detect a change
  of a method, its settings or the force estimator, beyond the differences
  between compilers, MPI libraries and decompositions (about 1e-5). They
  are reproduction checks, not accuracy targets.
- `submit` is optional: `./Allverify` runs the cases locally, or in a batch
  job of your own. Its options are `--env-script` (default
  `$S4F_ENV_SCRIPT`), `--partition`, `--max-cpus`, `--time`, `--mem` and
  `--sbatch-args`.

The tests need no OpenFOAM and never call SLURM. They check the driver
with synthetic data. They also check that the tracked references,
`reference/results/` and the tables of `results.md` agree, and that the
figures are regenerated from `reference/results/`; this last test is
skipped without matplotlib and scipy:

```bash
python3 -m unittest discover -s scripts -p 'test_*.py'
```

The previews in `images/` are written by
`python3 scripts/plot_verification.py --out images --dpi 110 --png-only`;
without `--png-only`, PDF and SVG versions are also written. The force
histories in the figures are drawn from the compact histories of
`reference/results/` (every 0.02 s, for display only); the errors are
computed from the complete histories.

## Quadrature of the cutLink surface traction

The tracked results use the quadrature of the revision that produced them:
the greedy thinning of the candidate surface points, the only quadrature of
that revision. When the solver source of this repository has a
`tractionQuadrature` keyword, the driver writes `tractionQuadrature greedy;`
in the cutLink cases. That way the cases do not depend on the solver
default. If the keyword exists but `greedy` is no longer one of its
choices, the driver stops rather than compare another quadrature with the
tracked results.

`--traction-quadrature none|clustered` runs another quadrature where the
solver offers it. Those cases get a tagged name (for example
`cutLink-L1-dt0p0025-clustered`) and are not checked against the tracked
results. The other methods do not evaluate the surface traction and get no
quadrature entry.

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
`./Allverify analyse` of the 24 runs of `results.md`, all run with
solids4foam `f33b244a2`. The Level 4 runs and the ghostCell, Level 3
incremental and outer-corrector runs were added to the earlier runs of the
same revision; the earlier results are unchanged.

## References

- Wan, D., Turek, S. (2007). Fictitious boundary and moving mesh methods for
  the numerical simulation of rigid particulate flows. Journal of
  Computational Physics, 222(1), 28-56.
  <https://doi.org/10.1016/j.jcp.2006.06.002>
- Erzincanli, B., Sahin, M. (2013). An arbitrary Lagrangian-Eulerian
  formulation for solving moving boundary problems with large
  displacements and rotations. Journal of Computational Physics, 255,
  660-679. <https://doi.org/10.1016/j.jcp.2013.08.038>
