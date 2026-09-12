# Elastoplastic Cook's membrane verification

This opt-in study migrates the structured quadrilateral mesh sweep previously
kept in `solid-benchmarks/elastoPlasticity/cooksMembrane`. It copies the parent
tutorial into `verification/work/`, runs the updated-Lagrangian elastoplastic
formulation on successively refined meshes, and records the top-right vertical
displacement. The tutorial and normal regression tests are not modified.

With OpenFOAM and solids4foam loaded, run:

```bash
cd tutorials/solids/elastoplasticity/cooksMembrane/verification
./Allverify
```

The default cell counts per side are `3, 6, 12, 24, 48`, matching the first
five levels of the legacy study. `./Allverify --quick` runs the first three as
a smoke test, and `--levels 12,24,48` selects an explicit subset. The legacy
96 x 96 level remains available with `--levels 3,6,12,24,48,96`; expect it to
be substantially more expensive once widespread yielding occurs.
Results are written to the ignored `postProcessing/` directory as CSV and a
Markdown summary; complete case copies remain under the ignored `work/`
directory for diagnosis. Pass `--reuse` to retain completed levels when
resuming or post-processing a sweep.

The full study passes when the finest displacement is within 8% of the 7.0 mm
reference, which represents the converged 6.97--7.09 mm range from Simo and
Armero, Areias, and César de Sá et al.; its reference error must also decrease
over the sweep with positive net order. The default full study took about 8
minutes on an Apple M1 Max reference machine with OpenFOAM-v2512. The quick
mode checks only that all requested cases finish and produce finite
displacement values.
