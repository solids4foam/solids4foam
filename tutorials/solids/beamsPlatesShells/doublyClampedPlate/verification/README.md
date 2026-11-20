# doublyClampedPlate verification study

This opt-in study adds the mesh-convergence sweep that the tutorial README
describes. The plate is clamped on two opposite edges and free on the other
two, and there is no closed-form Kirchhoff solution for it, so the reference
is a converged Kirchhoff plate finite element solution of the same problem
(12-DOF ACM rectangular element, Richardson-extrapolated from 64 x 64 and
128 x 128 elements), recorded in
`reference/doublyClampedPlate_verification_references.json`. The study refines
the mesh, reads the deflection at the midpoint of a free edge and at the plate
centre from the written finite area field `w`, measures the errors against
that reference and checks the observed order of accuracy. It also reports the
order and the extrapolated limit obtained from the computed values alone, by
Richardson extrapolation of the three finest meshes.

It is deliberately separate from `regressionTest.sh`: the regression test
checks that the tutorial remains numerically stable, whereas this study
measures convergence. Nothing here is run by `tutorials/Alltest` or
`tutorials/Alltest-regression`.

## Running

Source an OpenFOAM.com environment (v2412 or newer) with solids4foam built,
and run:

```bash
cd tutorials/solids/beamsPlatesShells/doublyClampedPlate/verification
./Allverify
```

Useful options:

```bash
./Allverify --quick          # two coarsest levels, smoke test only
./Allverify --levels 1,2,3   # skip the finest level
./Allverify --reuse          # resume a sweep without re-running cases
./Allverify --coupled        # block-coupled solver, results under work/coupled
```

Each level is a complete copy of the tutorial under the ignored
`verification/work/` directory, so the tutorial itself and its regression test
are never modified. Results are written to the ignored
`verification/postProcessing/` directory as `mesh_convergence.csv`,
`verification_summary.md`, the error norms under `profiles/`, and, when
`gnuplot` is available, `doublyClampedPlate_convergence.pdf`.

## Mesh levels

Each level scales the two in-plane divisions of the single block by a power of
two and leaves the one through-thickness cell alone. The tutorial ships
level 2. The loop gain of the explicitly coupled `clampedMoment` condition
grows with refinement, so the study sets its `relaxationFactor` per level on
both clamped edges from the reference file:

| Level | Cells | `clampedMoment` relaxation |
|---:|:---:|---:|
| 1 | 25 x 25 | 0.01 |
| 2 | 50 x 50 | 0.01 |
| 3 | 100 x 100 | 0.002 |
| 4 | 200 x 200 | 0.001 |

The whole default sweep takes about an hour in serial on an Apple M1 Ultra:
levels 1 to 3 take about four minutes together and level 4 about 53 minutes,
since the number of correctors grows as the `clampedMoment` relaxation factor
is reduced. With `--coupled` the whole sweep takes a few seconds.

## Metrics

- `free_edge_midpoint`: the deflection at the midpoint of the free `top`
  edge, taken from the boundary values of `w` on the edges nearest the middle
  of the patch (the two central edges for an even number of cells). Reference
  0.0015882 m.
- `centre`: the deflection at the plate centre, the mean of the cells nearest
  to it (the four central cells for an even number of cells). Reference
  0.0013976 m.

The patch edge values are read in the order of the polyMesh patch faces, which
is the order the finite area mesh builds its patch edges in; the script checks
that each edge value is close to the value of its adjacent cell before using
them.

## Acceptance criteria

The reference is a converged solution of the same equations, so the error
against it is the convergence measure:

- The relative error of `free_edge_midpoint` must decrease monotonically with
  refinement.
- Its net order must be at least 0.8 (the free edge conditions are imposed
  from one-sided estimates, so the observed order falls towards one on fine
  meshes).
- The finest mesh must be within 2%.

A `--quick` run only exercises the two coarsest meshes, which are not expected
to meet the accuracy tolerance, so it checks only that the sweep produces
usable numbers.

`verification/Allclean` removes the generated `work/` and `postProcessing/`
directories, as does `./Allclean all` in the case directory; the plain case
`./Allclean` keeps them.

## Reference results

Recorded with OpenFOAM v2606 on an Apple M1 Ultra:

Deflection at the free-edge midpoint $$w_e$$ (reference 0.0015882 m) and at
the plate centre $$w_c$$ (reference 0.0013976 m), with their relative errors.

| Level | Cells | Corr. | w_e (m) | Error | w_c (m) | Error |
|---:|---:|---:|---:|---:|---:|---:|
| 1 | 625 | 7316 | 0.00163968 | +3.24% | 0.00141967 | +1.58% |
| 2 | 2500 | 7000 | 0.00160509 | +1.06% | 0.00140174 | +0.30% |
| 3 | 10000 | 28093 | 0.00159506 | +0.43% | 0.00139849 | +0.06% |
| 4 | 40000 | 49729 | 0.00159123 | +0.19% | 0.00139776 | +0.01% |

The free-edge midpoint error falls by a factor of two to three with every
refinement (successive orders 1.61, 1.30, 1.18; net order 1.36), falling
towards the first order of the one-sided free edge estimates, while the
centre deflection, which the second-order clamped edge treatment governs,
converges faster. Richardson extrapolation of the three finest meshes gives
an observed order of 1.39 and an extrapolated free-edge deflection of
0.00158885 m, 0.04% above the finite element reference, which confirms that
the scheme converges to the reference solution. The same values are obtained
with `--coupled`, in two outer iterations per level.

## References

A. Adini, R. W. Clough, Analysis of plate bending by the finite element
method, NSF Report G7337, 1960 (the ACM rectangular element of the reference
solution).
