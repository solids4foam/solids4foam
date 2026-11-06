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

The whole default sweep takes about 20 minutes in serial on an Apple M1 Ultra:
levels 1 to 3 take under two minutes together and level 4 about 17 minutes,
since the number of correctors grows as the `clampedMoment` relaxation factor
is reduced.

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
- Its net order must be at least 0.8 (the scheme is first order because the
  clamped and free edge conditions are imposed from one-sided estimates).
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
| 1 | 625 | 4380 | 0.00178226 | +12.22% | 0.00154300 | +10.40% |
| 2 | 2500 | 4131 | 0.00168758 | +6.26% | 0.00147415 | +5.48% |
| 3 | 10000 | 15758 | 0.00163917 | +3.21% | 0.00143750 | +2.86% |
| 4 | 40000 | 28117 | 0.00161403 | +1.63% | 0.00141799 | +1.46% |

The free-edge midpoint error halves with every refinement (successive orders
0.97, 0.96, 0.98; net order 0.97), and the centre deflection converges at the
same rate. Richardson extrapolation of the three finest meshes gives an
observed order of 0.945 and an extrapolated free-edge deflection of
0.00158685 m, 0.09% below the finite element reference, which confirms that
the scheme converges to the reference solution.

## References

A. Adini, R. W. Clough, Analysis of plate bending by the finite element
method, NSF Report G7337, 1960 (the ACM rectangular element of the reference
solution).
