# Elastoplastic Cook's membrane verification

This opt-in study is a structured quadrilateral mesh sweep. It copies the
parent tutorial into `verification/work/`, runs the updated-Lagrangian
elastoplastic formulation on successively refined meshes, and records the
top-right vertical displacement.

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

The default sweep uses structured quadrilateral meshes. The tutorial can also
be run on an unstructured triangular-prism mesh created by Gmsh from
`cooksMembrane.geo` (`./Allrun triangular`), and `./Allverify --mesh
triangular` sweeps that family instead. Its levels are the nominal Gmsh
spacings of the legacy study, 25 mm halved at each level, and the default
levels `2, 3, 4, 5, 6` give 31, 119, 413, 1 591 and 6 323 cells with Gmsh
4.15. The triangular family requires `gmsh`, `gmshToFoam` and
`changeDictionary`, uses the same pass criterion as the structured family, and
writes its outputs with a `triangular_` prefix.

The per-mesh tip displacements of the four published sources of the legacy
study are stored in `reference/` as CSV files: Simo and Armero (Q1/P0 and
Q1/E4 elements), César de Sá et al. (Q4, Q1E4 and Qi6), Areias (Simplas,
δ and α), and Abaqus (CPE4H). The summary lists the finest value of each
series and its difference from the finest solids4foam result, and, when
`gnuplot` is available, `cooksMembrane_mesh_convergence.pdf` plots the sweep
against all series as a function of the average cell spacing. These series
are reported for comparison only and do not enter the pass criterion.

`verification/Allclean` removes the generated `work/` and `postProcessing/`
directories, as does `./Allclean all` in the case directory; the plain case
`./Allclean` keeps them.
