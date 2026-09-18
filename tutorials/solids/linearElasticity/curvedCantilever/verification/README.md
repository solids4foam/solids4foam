# curvedCantilever verification study

This opt-in study adds the mesh-convergence sweep that the tutorial README
already describes in words but does not ship. The tutorial builds the
Timoshenko curved-beam analytical solution as the
`curvedCantileverAnalyticalSolution` function object and samples both the
computed and the analytical stress along the radial line at `theta = 45°`. The
study refines the mesh, measures the error against that analytical solution on
the sampled line, and checks the observed order of accuracy.

It is deliberately separate from `regressionTest.sh`: the regression test
checks that the tutorial remains numerically stable, whereas this study
measures convergence. Nothing here is run by `tutorials/Alltest` or
`tutorials/Alltest-regression`.

## Running

```bash
source ~/bin/load-openfoam v2512
cd tutorials/solids/linearElasticity/curvedCantilever/verification
./Allverify
```

Useful options:

```bash
./Allverify --quick          # two coarsest levels, smoke test only
./Allverify --levels 2,3,4,5 # finer family
./Allverify --reuse          # resume a sweep without re-running cases
```

Each level is a complete copy of the tutorial under the ignored
`verification/work/` directory, so the tutorial itself and its regression test
are never modified. Results are written to the ignored
`verification/postProcessing/` directory as `mesh_convergence.csv`,
`verification_summary.md`, the error norms under `profiles/`, and, when
`gnuplot` is available, `curvedCantilever_errorNorms.pdf`.

## Mesh levels

The case is two dimensional, so each level scales the circumferential and
radial divisions of the single block by a power of two and leaves the one
spanwise cell alone. The tutorial ships level 2:

| Level | Circumferential × radial | Cells |
|---:|:---:|---:|
| 1 | 50 × 5 | 250 |
| 2 | 100 × 10 | 1 000 |
| 3 | 200 × 20 | 4 000 |
| 4 | 400 × 40 | 16 000 |

The whole default sweep takes about 35 seconds in serial on an Apple M1 Ultra,
including rebuilding the tutorial's analytical-solution library once per level.

## Acceptance criteria

Unlike the `ellipticPlate` and `narrowTmember` studies, the reference here is
an exact analytical solution rather than a digitised figure, so the error
against it *is* the convergence measure and no separate self-convergence
metric is needed:

- The relative L2 error must decrease monotonically with refinement.
- Its net order must exceed 0.5.
- The finest mesh must be within 5%.

The error combines the `sigma_xx`, `sigma_xy`, and `sigma_yy` components
sampled along the line, normalised together by the L2 norm of the analytical
stress over the same points, so one figure describes the whole profile rather
than a single component that happens to be small.

A `--quick` run only exercises the two coarsest meshes, which are not expected
to meet the accuracy tolerance, so it checks only that the sweep produces
usable numbers.

## Reference results

Recorded with OpenFOAM v2512 on an Apple M1 Ultra:

| Level | Cells | Relative L2 error | Relative L∞ error |
|---:|---:|---:|---:|
| 1 | 250 | 0.2522 | 0.2849 |
| 2 | 1 000 | 0.07906 | 0.1115 |
| 3 | 4 000 | 0.02279 | 0.04366 |
| 4 | 16 000 | 0.006626 | 0.01834 |

The L2 error falls by a factor of 3.2, then 3.5, then 3.4 for each halving of
the cell size, giving a net order of 1.75, close to second order. The L∞ error
converges more slowly, at a net order of 1.32, which is expected since it is
dominated by the highest-gradient point on the line.

## References

S. Timoshenko, J. N. Goodier, *Theory of Elasticity*, 3rd edition,
McGraw-Hill, 1970, pure-bending solution for a curved bar.
