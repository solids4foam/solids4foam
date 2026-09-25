# Contact patch test verification

This opt-in study migrates the solver-formulation cases previously kept in
`solid-benchmarks/linearElasticity/contactPatchTest`. It runs the
non-conformal contact patch test with `linearGeometryTotalDisplacement`, then
checks the transmitted stress against the analytical solution. The original
study also ran `unsLinearGeometry`, which has been removed together with the
legacy mechanical model. The parent tutorial and its normal regression test
are not modified.

With OpenFOAM and solids4foam loaded, run:

```bash
cd tutorials/solids/linearElasticity/contactPatchTest/verification
./Allverify
```

Pass `--reuse` to retain a completed case when resuming or regenerating the
summary.

Complete case copies and logs are retained in the ignored `work/` directory.
The ignored `postProcessing/` directory receives `results.csv` and a Markdown
summary.

## Acceptance criteria

The run must converge and write all three quantities
used by the existing tutorial regression test. The average relative error in
the transmitted vertical stress must be no greater than 5%, and the maximum
von Mises stress must be within 3% of the analytical magnitude of 10 kPa. The
tolerances are defined in
`reference/contact_patch_test_verification_references.json`.

This is a formulation verification rather than a mesh-convergence study: it
keeps the tutorial's non-matching 5-by-5 and 8-by-8 interface meshes fixed.
