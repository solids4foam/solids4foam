# annularPlate verification study

This opt-in study adds the mesh-convergence sweep that the tutorial README
describes, for the three configurations of the annular plate. The reference
is the closed-form Kirchhoff solution of each configuration, which the script
evaluates itself: the axisymmetric solution for the uniformly loaded plate
clamped at one edge and free at the other, and the $$\cos\varphi$$ mode for the
cosine load. The study refines the mesh, measures the error of the free-edge
deflection, the relative $$L_2$$ errors of $$w$$ and $$M$$ over the plate and,
for the cosine load, the error of the normal gradient of $$M$$ on the free
edge, and checks the observed order of accuracy.

It is deliberately separate from `regressionTest.sh`: the regression test
checks that the default configuration remains numerically stable, whereas
this study measures convergence. Nothing here is run by `tutorials/Alltest` or
`tutorials/Alltest-regression`.

## Running

Source an OpenFOAM.com environment (v2412 or newer) with solids4foam built,
and run:

```bash
cd tutorials/solids/beamsPlatesShells/annularPlate/verification
./Allverify
```

Useful options:

```bash
./Allverify --quick                            # two coarsest levels only
./Allverify --variants clampedInnerCosineLoad  # one configuration
./Allverify --levels 1,2,3                     # skip the finest level
./Allverify --reuse                            # resume without re-running
```

Each level of each configuration is a complete copy of the tutorial under the
ignored `verification/work/` directory, run with `./Allrun <configuration>`,
so the tutorial itself and its regression test are never modified. Results are
written to the ignored `verification/postProcessing/` directory as
`mesh_convergence.csv`, `verification_summary.md`, the error norms under
`profiles/`, and, when `gnuplot` is available, `annularPlate_convergence.pdf`.

## Mesh levels

Each level scales the radial and circumferential divisions of the four
quadrant blocks by a power of two and leaves the one through-thickness cell
alone. The tutorial ships level 2. The loop gain of the explicitly coupled
`clampedMoment` condition grows with refinement, so the study sets its
`relaxationFactor` per level from the reference file:

| Level | Radial x circumferential cells | `clampedMoment` relaxation |
|---:|:---:|---:|
| 1 | 8 x 76 | 0.01 |
| 2 | 16 x 152 | 0.01 |
| 3 | 32 x 304 | 0.005 |
| 4 | 64 x 608 | 0.002 |

The whole default sweep takes about 25 minutes in serial on an Apple M1 Ultra:
the three level-4 cases take 8 to 9 minutes each and the rest under a minute.

## Metrics

- Free-edge deflection: the mean of the boundary values of `w` on the free
  edge for the axisymmetric configurations, and the amplitude of the
  $$\cos\varphi$$ mode, $$2 \langle w \cos\varphi \rangle$$, for the cosine
  load, each compared with the exact value at the free edge.
- Relative $$L_2$$ errors of `w` and `M` over all cells against the exact
  solution at the cell centres.
- For the cosine load, the normal gradient of `M` on the free edge, taken as
  the difference between the edge value and the adjacent cell value divided by
  their distance and projected on the $$\cos\varphi$$ mode, against the exact
  $$\partial M / \partial r$$. This is the quantity the curvature term of the
  `freeEdgeMoment` condition determines.

The patch edge values are read in the order of the polyMesh patch faces, which
is the order the finite area mesh builds its patch edges in; the script checks
that each edge value is close to the value of its adjacent cell before using
them.

## Acceptance criteria

The reference is exact, so the error against it is the convergence measure.
For every configuration:

- The relative error of the free-edge deflection must decrease monotonically
  with refinement.
- Its net order must be at least 0.7 (the scheme is first order because the
  clamped and free edge conditions are imposed from one-sided estimates, and
  the coarsest annulus mesh is still in the pre-asymptotic range).
- The finest mesh must be within 5%.

A `--quick` run only exercises the two coarsest meshes, which are not expected
to meet the accuracy tolerance, so it checks only that the sweep produces
usable numbers.

`verification/Allclean` removes the generated `work/` and `postProcessing/`
directories, as does `./Allclean all` in the case directory; the plain case
`./Allclean` keeps them.

## Reference results

Recorded with OpenFOAM v2606 on an Apple M1 Ultra:

Free-edge deflection $$w_f$$ (exact value given above each table), its
relative error, and the relative $$L_2$$ errors of $$w$$ and $$M$$.

### clampedInner

Exact free-edge deflection 7.50303e-06 m.

| Level | Cells | Corr. | w_f (m) | Error | L2 w | L2 M |
|---:|---:|---:|---:|---:|---:|---:|
| 1 | 608 | 3153 | 8.53410e-06 | +13.74% | 16.44% | 0.58% |
| 2 | 2432 | 3892 | 8.10394e-06 | +8.01% | 9.39% | 0.66% |
| 3 | 9728 | 6578 | 7.82577e-06 | +4.30% | 5.01% | 0.44% |
| 4 | 38912 | 12432 | 7.67009e-06 | +2.23% | 2.59% | 0.24% |

### clampedInnerCosineLoad

Exact free-edge deflection amplitude 7.21429e-06 m; the last column is the
error of the free-edge normal gradient of $$M$$.

| Level | Cells | Corr. | w_f (m) | Error | L2 w | L2 M | dM/dn |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | 608 | 3137 | 8.25431e-06 | +14.42% | 16.88% | 0.60% | +7.18% |
| 2 | 2432 | 3700 | 7.82408e-06 | +8.45% | 9.73% | 0.46% | +4.37% |
| 3 | 9728 | 6771 | 7.54286e-06 | +4.56% | 5.21% | 0.31% | +2.39% |
| 4 | 38912 | 12822 | 7.38465e-06 | +2.36% | 2.70% | 0.18% | +1.25% |

### clampedOuter

Exact free-edge deflection 4.60288e-06 m.

| Level | Cells | Corr. | w_f (m) | Error | L2 w | L2 M |
|---:|---:|---:|---:|---:|---:|---:|
| 1 | 608 | 1923 | 5.31784e-06 | +15.53% | 16.85% | 5.23% |
| 2 | 2432 | 2514 | 4.99829e-06 | +8.59% | 9.50% | 2.75% |
| 3 | 9728 | 4703 | 4.80987e-06 | +4.50% | 5.02% | 1.41% |
| 4 | 38912 | 8854 | 4.70866e-06 | +2.30% | 2.57% | 0.72% |

The free-edge deflection error halves with every refinement in all three
configurations: the successive orders are 0.78, 0.90, 0.95 (`clampedInner`),
0.77, 0.89, 0.95 (`clampedInnerCosineLoad`) and 0.85, 0.93, 0.97
(`clampedOuter`), approaching first order, with net orders of 0.88, 0.87 and
0.92. The $$L_2$$ error of $$w$$ converges at the same rate (net orders 0.89,
0.88, 0.90). The $$L_2$$ error of $$M$$ is below 1% on every mesh for the
inner-clamped configurations, since $$M$$ is nearly uniform there, and
converges at first order (net 0.96) for `clampedOuter`. For the cosine load the
free-edge normal gradient of $$M$$, which the curvature term of
`freeEdgeMoment` sets, converges with a net order of 0.84. The number of
correctors grows with refinement because the `clampedMoment` relaxation factor
has to be reduced.

## References

S. Timoshenko, S. Woinowsky-Krieger, *Theory of Plates and Shells*, 2nd
edition, McGraw-Hill, 1959, circular and annular plates.
