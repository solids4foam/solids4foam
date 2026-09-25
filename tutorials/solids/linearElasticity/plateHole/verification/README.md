# plateHole verification study

This opt-in study migrates the mesh studies from
`solid-benchmarks/linearElasticity/plateHole` into the tutorial itself. It
supports both the structured `blockMesh` family and the unstructured
triangular-prism Gmsh family from the benchmark. The driver reads the error
norms that the `plateHoleAnalyticalSolution` function object already prints
against the analytical plate-with-hole solution, so no extra sampling is
needed.

The reference here is an exact analytical solution rather than a digitised
curve, so the error itself is meaningful and the study checks the observed
order of accuracy directly.

It is deliberately separate from `regressionTest.sh`: the regression test
checks that the tutorial remains numerically stable, whereas this study checks
convergence towards the analytical solution. Nothing here is run by
`tutorials/Alltest` or `tutorials/Alltest-regression`.

## Running

Source an OpenFOAM environment with solids4foam built (with PETSc for the
`petscSnes` variants), and run:

```bash
cd tutorials/solids/linearElasticity/plateHole/verification
./Allverify
```

Useful options:

```bash
./Allverify --quick                 # two coarsest levels, smoke test only
./Allverify --mesh triangular       # unstructured triangular-prism family
./Allverify --variants petscSnes    # PETSc SNES solution algorithm
./Allverify --variants highOrder    # high-order least-squares discretisation
./Allverify --variants vertexCentred # vertex-centred formulation
./Allverify --mesh triangular --variants vertexCentred
./Allverify --levels 2,3,4          # a subset of the family
./Allverify --reuse                 # resume a sweep without re-running cases
./Allverify --keep-going            # continue after an individual case fails
```

The variants are the solution approaches the tutorial `Allrun` accepts on an
OpenFOAM.com installation: `segregated` (the default), `petscSnes` and
`highOrder`, plus the legacy benchmark's `vertexCentred` formulation. The
`petscSnes`, `highOrder`, and `vertexCentred` variants require a solids4foam
build with PETSc. The two cell-centred PETSc approaches are skipped silently
by the tutorial when PETSc is unavailable, which the driver reports as a
missing solver log. The foam-extend-only `pressureDisplacement*` approaches
are out of scope.

Each level is a complete copy of the tutorial under the ignored
`verification/work/` directory, so the tutorial itself and its regression test
are never modified. For the structured family the driver changes the block
divisions in `system/blockMeshDict`. For the triangular family it installs the
benchmark Gmsh geometry, converts the mesh with `gmshToFoam`, and restores the
OpenFOAM boundary types with `changeDictionary`. For the segregated variant,
the solidModel convergence tolerance `rTol` is tightened from its default of
`1e-6` to `1e-10` so that the iterative error stays well below the
discretisation error on the finer meshes. With the shipped tolerance the
structured displacement error stalls near `2e-9` m from level 4 onwards and
the convergence curve flattens.

Results are written to the ignored `verification/postProcessing/` directory.
The structured sweep produces `mesh_convergence.csv` and
`verification_summary.md`; triangular output names have a `triangular_`
prefix. Per-variant convergence data are stored under `profiles/`, and a PDF
plot is generated when `gnuplot` is available.

## Mesh levels

The tutorial ships level 2 of the family. Every level scales the two in-plane
block divisions by a power of two; the third division spans the single cell
between the `empty` `frontAndBack` patches and is left alone, as is the
grading. The five levels reproduce `blockMeshDict.1` to `blockMeshDict.5` of
the benchmark exactly.

| Level | Block divisions of the first block | Cells |
|---:|:---|---:|
| 1 | (5 5 1) | 250 |
| 2 | (10 10 1) | 1 000 |
| 3 | (20 20 1) | 4 000 |
| 4 | (40 40 1) | 16 000 |
| 5 | (80 80 1) | 64 000 |

Levels 1 to 5 form the default sweep and take about four and a half minutes in
total, serial, on an Apple M1 Ultra.

The unstructured family uses the benchmark's nominal Gmsh spacings. Its first
four levels form the default triangular sweep; level 5 is retained as an
explicit extended option because its roughly quarter-million-cell segregated
solve is substantially more expensive.

| Level | Nominal spacing (m) | Cells |
|---:|---:|---:|
| 1 | 0.1 | 978 |
| 2 | 0.05 | 4 026 |
| 3 | 0.025 | 16 296 |
| 4 | 0.0125 | 64 257 |
| 5 | 0.00625 | approximately 256 000 |

## Acceptance criteria

The four metrics are the mean L2 and LInf norms of the cell displacement
difference and of the XX component of the cell stress difference. For a full
sweep each metric must

- decrease monotonically with every refinement, and
- give a net order of accuracy, measured between the coarsest and finest
  meshes against the effective cell spacing, above the minimum recorded in
  `reference/plateHole_verification_references.json`.

For the structured family the minimum orders are 1.5 for both displacement
norms, 1.0 for the stress L2 norm, and 0.7 for the stress LInf norm. The
triangular minima are respectively 1.5, 1.5, 0.9, and 0.55, reflecting the
lower stress convergence order on the unstructured cells next to the hole.

The stress minima are lower than the displacement ones because the largest
stress error sits in the cells at the hole, where the gradient is steepest, so
the LInf stress error converges at roughly first order.

A `--quick` run only exercises the two coarsest meshes, which are still well
outside the asymptotic range, so it checks only that the sweep produces usable
finite positive numbers.

## Reference results

Recorded with OpenFOAM v2512 and the default `segregated` variant on an Apple
M1 Ultra:

| Level | Cells | Spacing | D L2 | D LInf | Stress L2 | Stress LInf | Time |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | 250 | 0.1233 | 2.966e-08 | 1.072e-07 | 24 036 | 154 295 | 1 |
| 2 | 1 000 | 0.0617 | 1.349e-08 | 4.216e-08 | 10 420 | 105 030 | <1 |
| 3 | 4 000 | 0.0308 | 4.323e-09 | 1.283e-08 | 4 027 | 58 268 | 3 |
| 4 | 16 000 | 0.0154 | 1.206e-09 | 3.936e-09 | 1 475 | 30 174 | 21 |
| 5 | 64 000 | 0.0077 | 3.199e-10 | 1.095e-09 | 529 | 15 263 | 158 |

The net orders are 1.63 and 1.65 for the displacement L2 and LInf norms and
1.38 and 0.83 for the stress norms. The level-to-level orders of the
displacement L2 norm are 1.14, 1.64, 1.84 and 1.92, that is, the expected
second order once the coarsest meshes are left behind; the stress L2 norm
reaches 1.48 and the stress LInf norm 0.98 between the two finest meshes.

The `petscSnes` variant solves the same discretisation to a tighter tolerance
and needs no `rTol` adjustment, so it is slightly cleaner and faster: net
orders of 1.83, 1.81, 1.39 and 0.84 in the same order, with level-to-level
displacement L2 orders of 1.54, 1.86, 1.95 and 1.98, in one minute and
47 seconds for the whole sweep. The `highOrder` variant reaches a net order of
2.68 for the displacement L2 norm on the two coarsest meshes alone.

The added legacy variants were recorded with OpenFOAM v2512 on the same
machine. The table lists net orders over each default mesh family; all values
pass the acceptance criteria above.

| Mesh family | Variant | D L2 | D LInf | Stress L2 | Stress LInf |
|:---|:---|---:|---:|---:|---:|
| structured | `vertexCentred` | 1.98 | 1.98 | 1.96 | 1.76 |
| triangular | `segregated` | 1.90 | 1.66 | 1.02 | 0.68 |
| triangular | `petscSnes` | 1.92 | 1.73 | 1.12 | 0.76 |
| triangular | `highOrder` | 3.39 | 3.28 | 2.76 | 2.03 |
| triangular | `vertexCentred` | 1.93 | 1.82 | 1.00 | 0.61 |

## References

[1] S. Timoshenko, J. N. Goodier, *Theory of Elasticity*, McGraw-Hill, 1951.

[2] I. Demirdžić, S. Muzaferija, M. Perić, Benchmark solutions of some
structural analysis problems using finite-volume method and multigrid
acceleration, *International Journal for Numerical Methods in Engineering*,
40(10), 1893–1908, 1997.
