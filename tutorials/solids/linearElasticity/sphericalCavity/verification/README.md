# Spherical-cavity verification study

This opt-in study migrates the polyhedral mesh-convergence sweep previously
kept in `solid-benchmarks/linearElasticity/sphericalCavity`. Each level is a
complete copy of the parent tutorial under `verification/work/`; the tutorial
and its normal regression tests are not modified.

Source a supported OpenFOAM.com environment with PETSc, ensure Gmsh is
available, and run:

```bash
cd tutorials/solids/linearElasticity/sphericalCavity/verification
./Allverify
```

The driver uses the tutorial's `petscSnes poly` configuration and Gmsh minimum
spacings of `0.02`, `0.01`, `0.005`, and `0.0025` m. These are the first four
refinements from the legacy study; the much larger next level remains available
with `--levels 0.02,0.01,0.005,0.0025,0.00125`. Use `--quick` for the first two
levels, `--levels` for a custom refinement sequence, and `--reuse` to retain
completed cases when resuming a sweep.

Results are written to the ignored `verification/postProcessing/` directory as
CSV plus `verification_summary.md`. The CSV records cell count, effective
spacing, and the mean-L2 and L-infinity errors for displacement and
`sigma_zz`, as calculated by the tutorial's analytical-solution function
object.

The full study passes when all four error norms decrease from the coarsest to
the finest mesh and have positive net convergence order. This checks
convergence directly toward the analytical solution without fixing
machine-dependent regression values. Quick mode checks only that both cases
finish and produce finite positive errors. The study is not run by
`tutorials/Alltest` or `tutorials/Alltest-regression`.
