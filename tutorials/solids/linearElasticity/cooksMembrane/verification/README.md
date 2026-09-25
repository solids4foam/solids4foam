# Linear-elastic Cook's membrane verification

This opt-in study migrates the structured quadrilateral mesh sweep previously
kept in `solid-benchmarks/linearElasticity/cooksMembrane`. It copies the parent
tutorial into `verification/work/`, runs the segregated formulation on
successively refined meshes, and records the top-right vertical displacement.
The tutorial and normal regression tests are not modified.

With OpenFOAM and solids4foam loaded, run:

```bash
cd tutorials/solids/linearElasticity/cooksMembrane/verification
./Allverify
```

The default cell counts per side are `3, 6, 12, 24, 48, 96`, matching the
first six levels of the legacy study. `./Allverify --quick` runs the first
three as a smoke test, and `--levels 12,24,48` selects an explicit subset.
Results are written to the ignored `postProcessing/` directory as CSV and a
Markdown summary; complete case copies remain under the ignored `work/`
directory for diagnosis. Pass `--reuse` to retain completed levels when
resuming or post-processing a sweep.

The full study passes when the finest displacement is within 1% of 32.24 mm,
the midpoint of the 32.20--32.28 mm very-fine finite-element range reported in
the tutorial README, and its reference error decreases over the sweep with
positive net order. The quick mode checks only that all requested cases finish
and produce finite displacement values.
