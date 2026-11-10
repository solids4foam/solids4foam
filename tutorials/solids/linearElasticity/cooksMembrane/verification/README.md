# Linear-elastic Cook's membrane verification

This opt-in study is a structured quadrilateral mesh sweep. It copies the
parent tutorial into `verification/work/`, runs the segregated formulation on
successively refined meshes, and records the top-right vertical displacement.

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

The per-mesh tip displacements of the two published sources of the legacy
study are stored in `reference/` as CSV files: Zienkiewicz and Taylor (Q4
elements) and Abaqus (CPE4H). The legacy study applied a traction 1000 times
smaller than the tutorial, so its displacements have been scaled by 1000 to
match the tutorial load. The summary lists the finest value of each series
and its difference from the finest solids4foam result, and, when `gnuplot` is
available, `cooksMembrane_mesh_convergence.pdf` plots the sweep against both
series as a function of the average cell spacing. These series are reported
for comparison only and do not enter the pass criterion.

`verification/Allclean` removes the generated `work/` and `postProcessing/`
directories, as does `./Allclean all` in the case directory; the plain case
`./Allclean` keeps them.
