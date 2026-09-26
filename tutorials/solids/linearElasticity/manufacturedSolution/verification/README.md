# Manufactured-solution verification study

This opt-in study migrates the active convergence variants from
`solid-benchmarks/linearElasticity/manufacturedSolution`. The default sweep
combines the segregated and PETSc SNES solution procedures with regular
hexahedral, dual-polyhedral, and distorted-hexahedral meshes. Their tetrahedral
variants remain available explicitly. Both high-order approaches,
`highOrder-movingLeastSquares` and `highOrder-kExactLeastSquares`, are included
on regular hexahedral and tetrahedral meshes with polynomial degrees p=1, p=2,
and p=3. The JSON `p` entry sets `polynomialOrder` in each copied case. The per-variant
JSON `stencil_extra_cells` input sets the extra cells for face and cell stencils.
It defaults to 45 for p=1, 55 for p=2, and 65 for p=3 in the supplied JSON, for
both reconstruction methods. These correspond to nominal displacement orders 2,
2 (instead of 3), and 4, respectively; stress involves a derivative,
with nominal orders 1, 2, and 3. High-order variants require a PETSc-enabled build.
Variant names include the polynomial degree, for example
`hex-highOrder-movingLeastSquares-p3`. The `tet-structural` mesh label identifies
the structured tetrahedral mesh from `gmsh/tet-structured.geo`; the driver maps
it to the `Allrun` mesh argument `tet`.

Source an OpenFOAM.com, OpenFOAM.org, or foam-extend environment, ensure the
tutorial library can be built, and run:

```bash
cd tutorials/solids/linearElasticity/manufacturedSolution/verification
./Allverify
```

The default levels use 5, 10, 20, and 40 cells per coordinate direction for
the hexahedral meshes and equivalent target spacings for Gmsh. Use `--quick`
for the first two levels, `--variants` for a comma-separated subset, `--levels`
for custom cell counts, and `--reuse` when resuming a sweep. For example:

```bash
./Allverify --quick
./Allverify --variants hex-segregated,poly-petscSnes
./Allverify --variants tet-structural-segregated,tet-structural-petscSnes --quick
./Allverify --variants hex-highOrder-movingLeastSquares-p3,tet-structural-highOrder-kExactLeastSquares-p2
```

Results are written under the ignored `verification/postProcessing/`
directory. Each variant must produce finite positive errors and, in a full
sweep, lower finest-mesh errors for the displacement and stress L2 and L-infinity
norms. Each variant has JSON `minimum_net_order` inputs for `displacement` and
`stress`, applied to both L2 and L-infinity. Full sweeps require the measured net
orders to meet or exceed these values. The net order uses the coarsest and
finest errors and effective mesh spacings; intermediate errors need not decrease
monotonically. The summary reports each measured order, threshold, and result.

Initial thresholds are 1.5 for displacement and 0.5 for stress for standard
methods. High-order thresholds are p+0.5 and p-0.5, respectively, allowing a 0.5
margin below the nominal orders. These are editable expectations, not thresholds
calibrated against a completed full sweep.
Quick runs check finite positive errors and successful solver completion,
without enforcing the order thresholds. The CSV records p for high-order runs.
`--reuse` requires matching degree and stencil metadata; older runs without
it are rerun.
The kExact displacement errors compare cell averages; MLS errors compare
point values at cell centres. Stress errors use cell-centre values in both.
The study is not run by
`tutorials/Alltest` or `tutorials/Alltest-regression`.
