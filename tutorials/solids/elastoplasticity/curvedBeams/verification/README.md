# curvedBeams verification study

This opt-in study turns the curved-beams contact benchmark from
`solid-benchmarks/elastoPlasticity/curvedBeams` into a mesh-convergence study
inside the tutorial itself. It refines the tutorial `blockMesh`, extracts the
total reaction force history on the `fixed` patch of the lower beam from the
`solidForces` function object that the tutorial already runs, and compares the
history with the published curves of Neto et al. (2016).

The `solid-benchmarks` case varies the Coulomb friction coefficient over
`0.0`, `0.3` and `0.6` on a single mesh. Those three friction coefficients are
retained here as the study variants, and the mesh sweep is added on top: each
variant is run through the mesh family below and compared with the published
curve for its own friction coefficient.

It is deliberately separate from `regressionTest.sh`: the regression test
checks that the tutorial remains numerically stable, whereas this study checks
convergence towards a published benchmark. Nothing here is run by
`tutorials/Alltest` or `tutorials/Alltest-regression`.

## Running

The curvedBeams tutorial currently runs only in foam-extend, so the study does
too. Source a foam-extend environment with solids4foam built, and run:

```bash
source ~/bin/load-openfoam fe41
cd tutorials/solids/elastoplasticity/curvedBeams/verification
./Allverify
```

Useful options:

```bash
./Allverify --quick                    # two coarsest levels, smoke test only
./Allverify --levels 1,2,3             # the default sweep, stated explicitly
./Allverify --variants mu0.0,mu0.6     # frictionless and high-friction cases
./Allverify --reuse                    # resume a sweep without re-running cases
./Allverify --keep-going               # continue after an individual case fails
```

Each level is a complete copy of the tutorial under the ignored
`verification/work/` directory, so the tutorial itself and its regression test
are never modified. Results are written to the ignored
`verification/postProcessing/` directory as `mesh_convergence.csv`,
`verification_summary.md`, the extracted force histories under `profiles/`,
and, when `gnuplot` is available, `curvedBeams_reaction.pdf`.

No extra post-processing utility is needed: the reaction force is already
written by the `solidForces` function objects in the tutorial `controlDict`,
and the driver simply reads `postProcessing/0/solidForcesfixed.dat`.

## Mesh levels

Each beam is meshed with `radial x circumferential` divisions set by two m4
macros in `system/blockMeshDict.m4`, which the driver rewrites. The
circumferential count, which sets the number of contact faces, doubles between
levels:

| Level | Divisions per beam | Cells |
|---:|:---:|---:|
| 1 | `3 x 25` | 150 |
| 2 | `5 x 50` | 500 |
| 3 | `10 x 100` | 2 000 |
| 4 | `20 x 200` | 8 000 |

Level 2 is the mesh shipped with the tutorial; the radial count is tabulated
rather than computed because five cells cannot be halved exactly. Levels 1 to
3 form the default sweep, which takes under four minutes in serial on an Apple
M1 Ultra.

Level 4 is defined but **not** part of any recommended sweep: with the penalty
contact settings of the tutorial it diverges in the first displacement
increment, and the driver reports that failure rather than silently accepting
it. The divergence follows the circumferential count alone — `20 x 100` runs,
`10 x 200` does not — and it is insensitive to `penaltyScale` and to the
contact `relaxationFactor`, so it reflects the contact discretisation rather
than the choice of under-relaxation.

## Acceptance criteria

The reference curves in `reference/neto_mu*_reaction_*.txt` are digitised from
the published figures and therefore carry a digitisation uncertainty of roughly
one percent of the peak force. The difference from them stops falling once the
discretisation error drops below that floor, so it cannot serve as the
convergence measure on its own. The study therefore separates the two
questions:

- **Convergence** is measured by `force_rms_change_n`, the RMS change in the
  reaction force history, sampled at 64 equally spaced displacement stations
  in both components, between successive mesh levels. A full sweep requires
  this change to decrease monotonically with a positive net order.
- **Accuracy** is measured against the published curves on the finest mesh:
  `reference_rms_relative_error` and `reference_peak_relative_error`, each
  normalised by the peak reference force because the reaction passes through
  zero at the start and the end of the sliding, must both be within 5%.

A `--quick` run only exercises the two coarsest meshes, which are not expected
to meet the accuracy tolerances, so it checks only that the sweep produces
usable numbers.

The tutorial `Allrun` returns zero even when the solver diverges, so the driver
also inspects `log.solids4Foam` for a fatal error, a stack trace, or a run that
never reached its end time.

## Reference results

Recorded with foam-extend-4.1 and the default `mu0.3` variant on an Apple M1
Ultra, in serial:

| Level | Cells | Peak F<sub>x</sub> (N) | Peak F<sub>y</sub> (N) | RMS error vs. reference | Peak error vs. reference | Force RMS change (N) | Solver time (s) |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | 150 | 17.340 | 40.532 | 0.0437 | 0.0392 | – | 26 |
| 2 | 500 | 17.218 | 39.317 | 0.0196 | 0.0080 | 0.706 | 52 |
| 3 | 2 000 | 16.971 | 38.872 | 0.0225 | 0.0179 | 0.460 | 141 |

The reference peaks are 17.280 N in x and 39.005 N in y. The change in the
force history decreases monotonically and gives an observed net order of 0.62,
which is the order expected of a penalty contact solution with a sliding
contact front rather than of a smooth field. On the finest mesh the computed
history agrees with the published curves to 2.3% RMS and 1.8% at the peak; that
difference no longer falls between levels 2 and 3, which is consistent with the
digitisation uncertainty of the reference curves.

The frictionless `mu0.0` variant passes the same criteria over the same three
levels: RMS changes of 0.894 N and 0.596 N, a net order of 0.59, and a finest
mesh agreeing with the published curves to 3.2% RMS and 4.0% at the peak.

## References

D. Neto, M. Oliveira, L. Menezes, J. Alves, A contact smoothing method for
arbitrary surface meshes using Nagata patches, *Computer Methods in Applied
Mechanics and Engineering*, 299, 283–315, 2016.

I. Batistić, P. Cardiff, Ž. Tuković, A finite volume penalty based
segment-to-segment method for frictional contact problems, *Applied
Mathematical Modelling*, 101, 673–693, 2022.
