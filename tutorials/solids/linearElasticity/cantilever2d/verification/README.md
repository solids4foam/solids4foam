# cantilever2d verification study

This opt-in study adds the order-of-accuracy sweep that the tutorial README
describes in words but does not ship. It is a port of the cantilever
order-of-accuracy study in the `solid-benchmarks` repository
(`linearElasticity/cantilever`) onto the block-structured tutorial mesh, so
that the same measurement can be made with any of the tutorial's solution
approaches rather than only the vertex-centred one.

The tutorial already builds the Timoshenko slender-cantilever analytical
solution as the `cantileverAnalyticalSolution` function object, which writes
the mean L1, mean L2 and LInf norms of the displacement and stress error over
the whole domain to the solver log at the end of the run. This study refines
the mesh, reads those norms back, and checks the observed order of accuracy.

It is deliberately separate from `regressionTest.sh`: the regression test
checks that the tutorial remains numerically stable, whereas this study
measures convergence. Nothing here is run by `tutorials/Alltest` or
`tutorials/Alltest-regression`.

## Running

```bash
source ~/bin/load-openfoam v2512
cd tutorials/solids/linearElasticity/cantilever2d/verification
./Allverify
```

Useful options:

```bash
./Allverify --quick                    # two coarsest levels, smoke test only
./Allverify --levels 2,3,4,5           # finer family
./Allverify --variants vertexCentred   # a different solution approach
./Allverify --reuse                    # resume a sweep without re-running cases
./Allverify --keep-going               # continue after an individual case fails
```

The variants are the tutorial's own solution approaches: `petscSnes` (the
default), `vertexCentred`, `highOrder` and `segregated`. The default is the
Jacobian-free Newton-Krylov approach, because the tutorial README records that
the segregated approach needs about 123 s on 5 120 cells where the coupled
approaches need about 1 s; `segregated` is therefore left opt-in. The
`unsCoupled` approach is not offered, since it only runs with foam-extend.

Each level is a complete copy of the tutorial under the ignored
`verification/work/` directory, so the tutorial itself and its regression test
are never modified. The only edits made to a copy are the block divisions in
`system/blockMeshDict` and the relative `SOLIDS4FOAM_ROOT` in
`src/Make/options`, which is rewritten as an absolute path because the copy
sits deeper in the tree than the tutorial. Results are written to the ignored
`verification/postProcessing/` directory as `mesh_convergence.csv`,
`verification_summary.md`, the error norms under `profiles/`, and, when
`gnuplot` is available, `cantilever2d_errorNorms.pdf`.

## Mesh levels

The case is two dimensional, so each level scales the axial and
through-thickness divisions of the single block by a power of two and leaves
the one spanwise cell alone. The cells are square at every level. The tutorial
ships level 4:

| Level | Axial × through-thickness | Cells | Spacing (m) |
|---:|:---:|---:|---:|
| 1 | 40 × 2 | 80 | 0.05 |
| 2 | 80 × 4 | 320 | 0.025 |
| 3 | 160 × 8 | 1 280 | 0.0125 |
| 4 | 320 × 16 | 5 120 | 0.00625 |
| 5 | 640 × 32 | 20 480 | 0.003125 |

Levels 1 to 4 are the default sweep; level 5 is available through `--levels`.
The whole default sweep takes about 100 seconds in serial on an Apple M1 Ultra,
almost all of it rebuilding the tutorial's analytical-solution library once per
level; the solver itself needs at most a second per level.

## Acceptance criteria

As in the `curvedCantilever` study, and unlike `ellipticPlate` and
`narrowTmember`, the reference is an exact analytical solution rather than a
digitised figure, so the error against it *is* the convergence measure and no
separate self-convergence metric is needed:

- The relative displacement L2 error and the relative stress L2 error must
  both decrease monotonically with refinement.
- Both net orders must exceed 1.4.
- On the finest mesh both must be within 1%.

The absolute norms written by the function object are normalised by the peak
analytical response, so that the reported figures are relative: the
displacement by the analytical tip deflection $$u_y = 0.0160275$$ m, and the
stress by the peak analytical bending stress $$\sigma_{xx} = 120$$ MPa. The
three reported stress components, `xx`, `xy` and `yy`, are combined into a
single root-sum-square figure before normalising, so one number describes the
whole stress field rather than one component that happens to be small.

A `--quick` run only exercises the two coarsest meshes, which are not expected
to meet the accuracy tolerance, so it checks only that the sweep produces
usable finite positive numbers.

One exception applies to the `highOrder` variant: below a roundoff floor of
$$10^{-8}$$ the monotone and order checks are skipped, because that
discretisation reproduces the analytical solution exactly and its remaining
error is arithmetic roundoff, which neither decreases nor has an order.

## Reference results

Recorded with OpenFOAM v2512 on an Apple M1 Ultra.

### Default `petscSnes` variant

| Level | Cells | Displacement L2 | Displacement L∞ | Stress L2 | Stress L∞ |
|---:|---:|---:|---:|---:|---:|
| 1 | 80 | 0.1703 | 0.3441 | 0.1031 | 0.1465 |
| 2 | 320 | 0.07135 | 0.1455 | 0.04720 | 0.1056 |
| 3 | 1 280 | 0.02228 | 0.04566 | 0.01527 | 0.03845 |
| 4 | 5 120 | 0.006043 | 0.01242 | 0.004308 | 0.01193 |

The net orders over the family are 1.61 for the displacement L2 error and 1.53
for the stress L2 error. Both are still approaching their asymptotic values
over this family: the displacement L2 error falls by a factor of 2.39, then
3.20, then 3.69 for each halving of the cell size, an order of 1.25, then 1.68,
then 1.88, so the last refinement is close to second order. The coarsest mesh
resolves the thickness with only two cells, which is what holds the net order
down: `./Allverify --levels 3,4,5` skips it and gives a displacement L2 error
of 0.02228, 0.006043 and 0.001559 for a net order of 1.92. The L∞ errors
converge more slowly, as expected, since they are dominated by the single worst
cell.

### Other variants

| Variant | Displacement L2 net order | Stress L2 net order | Finest displacement L2 |
|:---|---:|---:|---:|
| `petscSnes` | 1.61 | 1.53 | 0.006043 |
| `vertexCentred` | 2.07 | 2.00 | 0.0008479 |
| `highOrder` | exact | exact | 1.6e-10 |

The vertex-centred approach is cleanly second order on this family and is about
seven times more accurate than the cell-centred approach on the same mesh. The
high-order approach reproduces the analytical solution to machine precision
from level 2 onwards: the Timoshenko displacement field is a cubic polynomial,
which the high-order reconstruction represents exactly. Level 1 is the
exception, at an error of 0.19, because two cells through the thickness do not
support the high-order stencil. Running the `highOrder` and `vertexCentred`
sweeps together takes about 3.5 minutes.

## References

[1] C.E. Augarde, A.J. Deeks, The use of Timoshenko's exact solution for a
cantilever beam in adaptive analysis. *Finite Elements in Analysis and Design*,
44, 2008, 595–601, 10.1016/j.finel.2008.01.010.

[2] S. Timoshenko, J. N. Goodier, *Theory of Elasticity*, 3rd edition,
McGraw-Hill, 1970.
