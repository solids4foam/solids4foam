# ellipticPlate verification study

This opt-in study migrates the mesh-convergence study from
`solid-benchmarks/linearElasticity/ellipticPlate` into the tutorial itself. It
refines the tutorial `blockMesh` through the mesh family of Demirdžić et al.
(1997), samples the equivalent (von Mises) stress along the line
`r = 2.1 m`, `z = 0.3 m` at mid-thickness, and compares the result with the
published curve.

It is deliberately separate from `regressionTest.sh`: the regression test
checks that the tutorial remains numerically stable, whereas this study checks
convergence towards a published benchmark. Nothing here is run by
`tutorials/Alltest` or `tutorials/Alltest-regression`.

## Running

Source a supported OpenFOAM environment, build solids4foam with PETSc, and run:

```bash
source ~/bin/load-openfoam v2512
cd tutorials/solids/linearElasticity/ellipticPlate/verification
./Allverify
```

Useful options:

```bash
./Allverify --quick                    # two coarsest levels, smoke test only
./Allverify --levels 1,2,3,4,5         # add the 294,912-cell level
./Allverify --variants segregated      # implicit segregated solution algorithm
./Allverify --reuse                    # resume a sweep without re-running cases
```

Each level is a complete copy of the tutorial under the ignored
`verification/work/` directory, so the tutorial itself and its regression test
are never modified. Results are written to the ignored
`verification/postProcessing/` directory as `mesh_convergence.csv`,
`verification_summary.md`, the sampled profiles under `profiles/`, and, when
`gnuplot` is available, `ellipticPlate_sigmaEq.pdf`.

## Mesh levels

Every level scales the tutorial block divisions by a factor of two and retains
the tutorial `edgeGrading`, reproducing the mesh family used in the reference:

| Level | Block divisions | Cells |
|---:|:---:|---:|
| 1 | `6 x 4 x 3` | 72 |
| 2 | `12 x 8 x 6` | 576 |
| 3 | `24 x 16 x 12` | 4 608 |
| 4 | `48 x 32 x 24` | 36 864 |
| 5 | `96 x 64 x 48` | 294 912 |

Level 3 is the mesh shipped with the tutorial. Levels 1 to 4 form the default
sweep; level 5 is available explicitly. The full 1–5 `petscSnes` sweep takes
roughly three minutes in serial on an Apple M1 Ultra, with level 5 accounting
for most of that.

## Acceptance criteria

The reference curve in `reference/demirdzic_sigmaEq.txt` is digitised from
Figure 7 of Demirdžić et al. (1997) and therefore carries a digitisation
uncertainty of roughly one percent. The difference from it stops falling once
the discretisation error drops below that floor, so it cannot serve as the
convergence measure on its own. The study therefore separates the two
questions:

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

| Level | Cells | Peak σ<sub>eq</sub> (MPa) | RMS error vs. reference | Profile RMS change (Pa) |
|---:|---:|---:|---:|---:|
| 1 | 72 | 3.920 | 0.1117 | – |
| 2 | 576 | 3.438 | 0.0534 | 3.23e+05 |
| 3 | 4 608 | 3.573 | 0.0161 | 1.60e+05 |
| 4 | 36 864 | 3.563 | 0.0231 | 4.70e+04 |
| 5 | 294 912 | 3.563 | 0.0251 | 1.19e+04 |

The peak stress is converged to four significant figures between levels 4 and
5, the profile change between the two finest levels gives an observed order of
1.98, and the finest mesh agrees with the published curve to 2.5% RMS and 1.0%
at the peak. The residual difference is consistent with the digitisation
uncertainty of the reference curve.

## References

I. Demirdžić, S. Muzaferija, M. Perić, Benchmark solutions of some structural
analysis problems using finite-volume method and multigrid acceleration,
*International Journal for Numerical Methods in Engineering*, 40(10),
1893–1908, 1997.
