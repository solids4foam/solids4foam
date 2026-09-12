# narrowTmember verification study

This opt-in study migrates the mesh study from
`solid-benchmarks/linearElasticity/narrowTmember` into the tutorial itself. It
refines the tutorial `blockMesh` through the mesh family of Demirdžić et al.
(1997), samples the equivalent (von Mises) stress along the arc `r = 1.5R` in
the `z = 0` plane, where `R = 5 mm` is the fillet radius, and compares the
result with the published curve.

It is deliberately separate from `regressionTest.sh`: the regression test
checks that the tutorial remains numerically stable, whereas this study checks
convergence towards a published benchmark. Nothing here is run by
`tutorials/Alltest` or `tutorials/Alltest-regression`.

## Running

```bash
source ~/bin/load-openfoam v2512
cd tutorials/solids/linearElasticity/narrowTmember/verification
./Allverify
```

Useful options:

```bash
./Allverify --quick                 # two coarsest levels, smoke test only
./Allverify --variants segregated   # implicit segregated solution algorithm
./Allverify --levels 2,3,4          # a subset of the family
./Allverify --reuse                 # resume a sweep without re-running cases
```

Each level is a complete copy of the tutorial under the ignored
`verification/work/` directory, so the tutorial itself and its regression test
are never modified. Results are written to the ignored
`verification/postProcessing/` directory as `mesh_convergence.csv`,
`verification_summary.md`, the sampled profiles under `profiles/`, and, when
`gnuplot` is available, `narrowTmember_sigmaEq.pdf`.

## Mesh levels

The tutorial ships level 3 of the family. Every level scales its block
divisions by a power of two and retains the tutorial grading:

| Level | Cells |
|---:|---:|
| 1 | 624 |
| 2 | 4 992 |
| 3 | 39 936 |
| 4 | 319 488 |

Levels 1 to 4 form the default sweep and take about two minutes in total on an
Apple M1 Ultra. Level 5 is reachable with `--levels`, but at roughly 2.6
million cells it exceeded the memory of the 64 GB machine used here, so choose
it only on a larger machine.

## Acceptance criteria

As in the `ellipticPlate` study, the reference curve is read from a published
figure and carries a digitisation uncertainty of roughly one percent. The
difference from it stops falling once the discretisation error drops below that
floor, so the two questions are separated:

- **Convergence** is measured by `profile_rms_change_pa`, the RMS change in the
  sampled stress profile between successive mesh levels. A full sweep requires
  this change to decrease monotonically with a positive net order.
- **Accuracy** is measured against the published curve on the finest mesh:
  `reference_rms_relative_error` and `reference_peak_relative_error` must both
  be within 5%.

A `--quick` run only exercises the two coarsest meshes, which are not expected
to meet the accuracy tolerances, so it checks only that the sweep produces
usable numbers.

## Reference results

Recorded with OpenFOAM v2512 and the `petscSnes` variant on an Apple M1 Ultra:

| Level | Cells | Peak σ<sub>eq</sub> (MPa) | RMS error vs. reference | Profile RMS change (Pa) | Solver wall clock (s) |
|---:|---:|---:|---:|---:|---:|
| 1 | 624 | 6.866 | 0.1466 | – | <1 |
| 2 | 4 992 | 6.397 | 0.0456 | 2.35e+05 | 1 |
| 3 | 39 936 | 6.338 | 0.0206 | 7.66e+04 | 6 |
| 4 | 319 488 | 6.313 | 0.0170 | 2.05e+04 | 93 |

The profile change falls by a factor of 3.1 and then 3.7 for each halving of
the cell size, giving a net order of 1.76, and the finest mesh agrees with the
published curve to 1.7% RMS and 0.68% at the peak stress, which is within the
digitisation uncertainty of the reference.

## Sampled quantity

The arc is the one used in the published figure: radius `1.5R = 7.5 mm` about
the fillet centre, in the `z = 0` symmetry plane, sampled at 30 points from
`theta = -90°` to `theta = 180°`. The published angle is measured from the
point `(0, -1.5R)`, which is `theta + 90°`; the driver applies that offset when
comparing.

## References

I. Demirdžić, S. Muzaferija, M. Perić, Benchmark solutions of some structural
analysis problems using finite-volume method and multigrid acceleration,
*International Journal for Numerical Methods in Engineering*, 40(10),
1893–1908, 1997.
