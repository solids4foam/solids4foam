# idealisedVentricle verification study

This opt-in study migrates the mesh-convergence study from
`solid-benchmarks/hyperElasticity/idealisedVentricle` into the tutorial itself.
It refines the tutorial `blockMesh` and its rotational `extrudeMesh` together,
samples the deformed mid-wall line of Problem 2 of Land et al. (2015), and
checks that the solution converges under uniform refinement.

It is deliberately separate from `regressionTest.sh`: the regression test
checks that the tutorial remains numerically stable, whereas this study checks
mesh convergence. Nothing here is run by `tutorials/Alltest` or
`tutorials/Alltest-regression`.

## Running

Source an OpenFOAM.com environment, build solids4foam with PETSc, and run:

```bash
source ~/bin/load-openfoam v2512
cd tutorials/solids/hyperelasticity/idealisedVentricle/verification
./Allverify
```

Useful options:

```bash
./Allverify --quick             # levels 1 and 2 only, smoke test
./Allverify --levels 1,2,3,4    # add the 829,440-cell level
./Allverify --cores 16          # fixed rank count for every level
./Allverify --levels 1,2,3,4 --cores 1,8,32,64   # one rank count per level
./Allverify --reuse             # resume a sweep without re-running cases
```

With `--cores auto` (the default) level 1 runs in serial and levels 2, 3, and 4
run on 8, 16, and 16 ranks respectively, through the tutorial's own
`./Allrun petsc parallel` path. Use `--cores N` when a scheduler allocation
requires a fixed rank count, or pass one value per requested level when the
allocation should grow with the mesh.

Each level is a complete copy of the tutorial under the ignored
`verification/work/` directory, so the tutorial itself and its regression test
are never modified. Results are written to the ignored
`verification/postProcessing/` directory as `mesh_convergence.csv`,
`verification_summary.md`, the sampled mid-wall lines under `profiles/`, and,
when `gnuplot` is available, `idealisedVentricle_midLine.pdf`.

## Mesh levels

Each level doubles the block divisions and the number of rotational extrusion
layers together, so the cell count grows by a factor of eight:

| Level | Block divisions | Extrusion layers | Cells |
|---:|:---:|---:|---:|
| 1 | `15 x 3 x 1` | 36 | 1 620 |
| 2 | `30 x 6 x 1` | 72 | 12 960 |
| 3 | `60 x 12 x 1` | 144 | 103 680 |
| 4 | `120 x 24 x 1` | 288 | 829 440 |

Level 1 is the mesh shipped with the tutorial. Levels 1 to 3 form the default
sweep; level 4 is available explicitly and needs substantially more memory and
wall-clock time. Choose levels appropriate for the available machine.

## Acceptance criteria

Land et al. publish Problem 2 as a cross-code comparison rather than as a
closed-form solution, so there is no analytical value to converge to. The
acceptance criterion is therefore self-convergence:

- `mid_line_rms_change_m`, the RMS distance between the deformed mid-wall lines
  of successive mesh levels, must decrease monotonically with a positive net
  order.
- Every level must produce a finite, negative apex position, that is, an
  inflated ventricle.

`diagnostic_apex_z_m` records the apex position obtained from an independent
solids4foam solution of the same configuration using the block-coupled mixed
pressure-displacement solid model. It is reported in the CSV and the summary as
a cross-formulation check, and is deliberately **not** a pass criterion, since
it is itself a numerical solution rather than a published reference.

A `--quick` run only exercises the two coarsest meshes and checks only that
the sweep produces usable numbers.

## Reference results

Recorded with OpenFOAM v2512 on a 20-core Apple M1 Ultra using the default
`--cores auto` allocation. The whole default sweep took 17 minutes.

| Level | Cells | Ranks | Apex position (mm) | Mid-wall RMS change (m) | Solver wall clock (s) |
|---:|---:|---:|---:|---:|---:|
| 1 | 1 620 | 1 | -26.487 | – | 50 |
| 2 | 12 960 | 8 | -27.268 | 5.95e-04 | 111 |
| 3 | 103 680 | 16 | -27.324 | 1.15e-04 | 817 |

The change in the mid-wall line falls by a factor of 5.2 for a halving of the
cell size, giving a net order of 2.38, and the apex position is converged to
three significant figures between the two finest meshes. The finest apex
position differs from the independent block-coupled mixed pressure-displacement
solution by 1.0%.

## Sampled quantity

The mid-wall line is the curve used in the Land et al. Problem 2 figures. It is
the ellipse midway between the endocardial and epicardial surfaces,

```
x = r_s sin(u),    y = 0,    z = r_l cos(u)
```

with `r_s = (7 + 10)/2 mm`, `r_l = (17 + 20)/2 mm`, and `u` running from `-pi`
at the apex to `-acos(1/4)` at the truncated base, sampled at 100 points. The
driver samples the displacement field `D` at those points with `cellPoint`
interpolation and adds it to the undeformed coordinates.

## References

S. Land, V. Gurev, S. Arens, C. M. Augustin, L. Baron, R. Blake, C. Bradley,
S. Castro, A. Crozier, M. Favino, T. E. Fastl, T. Fritz, H. Gao, A. Gizzi,
B. E. Griffith, D. E. Hurtado, R. Krause, X. Luo, M. P. Nash, S. Pezzuto,
G. Plank, S. Rossi, D. Ruprecht, G. Seemann, N. P. Smith, J. Sundnes,
J. J. Rice, N. Trayanova, D. Wang, Z. J. Wang, S. A. Niederer, Verification of
cardiac mechanics software: benchmark problems and solutions for testing active
and passive material behaviour, *Proceedings of the Royal Society A*, 471(2184),
20150641, 2015. https://doi.org/10.1098/rspa.2015.0641
