# manufacturedSolutionPlate verification study

This opt-in study adds the mesh-convergence sweep that the tutorial README
describes. The reference is the manufactured solution defined in
`writeManufacturedSolution.py` in the case directory, which the driver imports,
so the plate data and the polynomial are kept in one place. The study refines
the mesh for each mesh type, reads the written finite area fields `w`, `M` and
`theta`, measures their errors against the exact solution and checks the
observed order of accuracy where the discretisation is expected to converge.

It is deliberately separate from `regressionTest.sh`: the regression test
checks that the tutorial remains numerically stable, whereas this study
measures convergence. Nothing here is run by `tutorials/Alltest` or
`tutorials/Alltest-regression`.

## Running

Source an OpenFOAM.com environment (v2412 or newer) with solids4foam built,
and run:

```bash
cd tutorials/solids/beamsPlatesShells/manufacturedSolutionPlate/verification
./Allverify
```

Useful options:

```bash
./Allverify --quick                 # two coarsest levels, smoke test only
./Allverify --variants hex,distHex  # skip the triangles (no gmsh needed)
./Allverify --levels 1,2,3          # skip the finest level
./Allverify --reuse                 # resume a sweep without re-running cases
./Allverify --coupled               # block-coupled solver, results under work/coupled
```

Each level is a complete copy of the tutorial under the ignored
`verification/work/` directory, so the tutorial itself and its regression test
are never modified. Results are written to the ignored
`verification/postProcessing/` directory as `mesh_convergence.csv`,
`verification_summary.md`, the error norms under `profiles/`, and, when
`gnuplot` is available, `manufacturedSolutionPlate_convergence.pdf`.

The `tri` variant needs `gmsh`. The segregated sweep takes about an hour for
the quadrilateral variants and longer for the triangles, most of it on the two
finest levels; the `--coupled` sweep takes a few seconds per variant.

## Mesh levels and variants

Level $$n$$ uses $$25 \cdot 2^{n-1}$$ cells per side for the quadrilateral
meshes and the same target spacing for the triangles; the tutorial ships
level 2. The three variants are the mesh types of the tutorial `Allrun`:

| Variant | Mesh | Orders checked |
|:---|:---|:---:|
| `hex` | uniform quadrilaterals from `blockMesh` | yes |
| `distHex` | interior points moved randomly by up to 15% of the spacing | no |
| `tri` | unstructured triangles from Gmsh | no |

The loop gain of the explicitly coupled `clampedMoment` condition grows with
refinement, so the study sets its `relaxationFactor` per level from the
reference file (0.01, 0.01, 0.002 and 0.001 for levels 1 to 4); the coupled
solver does not use it.

## Metrics

All errors are relative. The $$L_2$$ norms are area weighted over the cells
and normalised by the norm of the exact field; the $$L_\infty$$ norms are
normalised by the largest exact value:

- `w_l2_error`, `w_linf_error`: deflection at the cell centres;
- `M_l2_error`, `M_linf_error`: moment sum at the cell centres;
- `theta_l2_error`: rotation $$\theta = -\nabla w$$, the least-squares
  gradient of the computed deflection;
- `free_edge_w_linf_error`, `free_edge_M_linf_error`: boundary values on the
  free edge;
- `clamped_edge_M_linf_error`: boundary values of `M` on the clamped edge.

The finite area boundary values are compared with the exact solution at the
foot of the normal from the adjacent cell centre, which is where OpenFOAM.com
places them (`faPatch::delta()` is the normal projection of the vector from
the cell centre to the edge centre). The patch edge values are read in the
order of the polyMesh patch faces, which is the order the finite area mesh
builds its patch edges in; the script checks that each edge value is close to
the value of its adjacent cell before using them.

## Acceptance criteria

For the variants with `check_orders` in the reference file (`hex`):

- The $$L_2$$ error of `w` must decrease monotonically with refinement and be
  at most 0.5% on the finest level.
- Its net order must be at least 0.8 and that of the $$L_2$$ error of `M` at
  least 0.7. The free edge conditions use one-sided estimates, so first order
  is the expected behaviour.

The other variants must only complete with finite errors; their results are
reported so that the behaviour on skewed meshes is visible.

Quick runs (`--quick`) check completion and finite errors only.

## Reference results

Recorded with OpenFOAM.com v2606 on an Apple M1 Ultra with `--coupled`, which
converges in one solve (the second iteration confirms it); the segregated
solver gives the same values to better than $$10^{-5}$$ relative (checked on
the two coarsest levels of every variant). The $$L_\infty$$ norms are in the
summary file; they are about twice the $$L_2$$ norms on the uniform mesh.

### hex

| Level | Cells | L2 w | L2 M | L2 theta | Free w | Free M | Clamped M |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | 625 | 1.22% | 1.07% | 1.68% | 2.00% | 4.54% | 0.44% |
| 2 | 2500 | 0.57% | 0.59% | 0.86% | 0.98% | 2.23% | 0.13% |
| 3 | 10000 | 0.28% | 0.31% | 0.44% | 0.49% | 1.11% | 0.04% |
| 4 | 40000 | 0.14% | 0.16% | 0.22% | 0.24% | 0.55% | 0.01% |

Successive orders of the $$L_2$$ error of `w`: 1.09, 1.05, 1.03; net orders
1.06 (`w`), 0.91 (`M`), 0.98 (`theta`). The clamped edge moment converges at
second order, the free edge values at first order. The segregated solver
needs 2873 and 4991 correctors on levels 1 and 2.

### distHex

| Level | Cells | L2 w | L2 M | L2 theta | Free w | Free M | Clamped M |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | 625 | 1.66% | 1.30% | 1.94% | 3.13% | 6.97% | 2.81% |
| 2 | 2500 | 2.87% | 1.41% | 2.18% | 3.43% | 8.72% | 3.90% |
| 3 | 10000 | 0.47% | 0.59% | 0.69% | 0.90% | 8.11% | 4.16% |
| 4 | 40000 | 1.05% | 0.51% | 0.77% | 1.13% | 9.02% | 4.00% |

Each level has a different random perturbation, and the errors vary with it
instead of decreasing: the net orders are 0.2 (`w`) and 0.5 (`M`).

### tri

| Level | Cells | L2 w | L2 M | L2 theta | Free w | Free M | Clamped M |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | 1472 | 42.7% | 51.7% | 69.2% | 83.2% | 116.2% | 16.7% |
| 2 | 5828 | 44.4% | 52.5% | 70.1% | 89.1% | 122.4% | 18.8% |
| 3 | 23266 | 42.4% | 50.7% | 67.3% | 83.2% | 135.8% | 8.9% |
| 4 | 92578 | 38.7% | 48.0% | 63.6% | 70.8% | 134.1% | 16.2% |

The errors do not decrease with refinement.

## What limits the accuracy on skewed meshes

The error on the perturbed and triangular meshes was traced by replacing each
edge condition in turn by the exact boundary values (evaluated at the foot of
the normal) and by distorting only parts of the mesh:

- With exact values on all edges, the interior discretisation on a smoothly
  skewed mesh gives the same 0.17% error as the uniform mesh.
- The clamped condition, the symmetry condition and the `freeEdgeMoment`
  condition keep their accuracy on skewed boundary cells.
- The `freeEdgeDisplacement` condition does not: distorting only the row of
  cells along the free edge reproduces the whole error, and the error grows in
  proportion to the tangential offset between the cell centres and their edge
  centres while it does not change with the mesh spacing.

The error is a smooth deflection mode, uniform along the free edge and shaped
like a cantilever under a tip shear force. The displacement condition
extrapolates the edge value of `w` from the cell gradient, which imposes the
vanishing bending moment only through a term of the order of the mesh spacing,
so the discrete system is weakly constrained against this mode. On a uniform
mesh the truncation errors of the condition are symmetric and do not excite it;
on a skewed boundary they do. The condition is consistent, but the accuracy on
general meshes requires a formulation that imposes the vanishing bending moment
through the edge value of `M`, which is left for further work.

The study also exposed, and fixed, a scaling error of the explicit
non-orthogonal correction in the block-coupled solver, which affected skewed
meshes only; the correction is now assembled implicitly, so the coupled solver
reproduces the converged segregated solution in one solve on any mesh.
