# Manufactured-solution verification study

This opt-in study migrates the active convergence variants from
`solid-benchmarks/linearElasticity/manufacturedSolution`. The default sweep
combines the segregated and PETSc SNES solution procedures with regular
hexahedral, dual-polyhedral, and distorted-hexahedral meshes. Tetrahedral
variants, which were disabled in the legacy driver, remain available
explicitly. Both cubic high-order approaches, `highOrder-movingLeastSquares`
and `highOrder-kExactLeastSquares`, are included on the regular hex mesh in the
default sweep. They are also available explicitly for `tet` and `poly` meshes
and require a PETSc-enabled build.

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
./Allverify --variants tet-segregated,tet-petscSnes --quick
./Allverify --variants hex-highOrder-movingLeastSquares,hex-highOrder-kExactLeastSquares
```

Results are written under the ignored `verification/postProcessing/`
directory. Each variant must produce finite positive errors and, in a full
sweep, lower finest-mesh errors with positive net convergence order for the
displacement and stress L2 and L-infinity norms. Quick runs check finite positive
errors and successful solver completion without asserting convergence order.
The kExact displacement errors compare cell averages; MLS errors compare
point values at cell centres. Stress errors use cell-centre values in both.
The study is not run by
`tutorials/Alltest` or `tutorials/Alltest-regression`.
