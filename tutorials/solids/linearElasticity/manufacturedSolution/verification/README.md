# Manufactured-solution verification study

This opt-in study migrates the active convergence variants from
`solid-benchmarks/linearElasticity/manufacturedSolution`. The default sweep
combines the segregated and PETSc SNES solution procedures with regular
hexahedral, dual-polyhedral, and distorted-hexahedral meshes. Tetrahedral
variants, which were disabled in the legacy driver, remain available
explicitly.

Source an OpenFOAM.com environment, ensure the tutorial library can be built,
and run:

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
```

Results are written under the ignored `verification/postProcessing/`
directory. Each variant must produce finite positive errors and, in a full
sweep, lower finest-mesh errors with positive net convergence order for the
displacement and stress L2 and L-infinity norms. The study is not run by
`tutorials/Alltest` or `tutorials/Alltest-regression`.
