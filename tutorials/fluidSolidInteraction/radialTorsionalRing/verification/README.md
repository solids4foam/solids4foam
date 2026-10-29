# radialTorsionalRing verification

This opt-in study runs the tutorial on up to three meshes and three drive
configurations and analyses the results against the semi-analytical reference. It
is deliberately separate from `regressionTest.sh`, and nothing here is run by
`tutorials/Alltest` or `tutorials/Alltest-regression`. Scope and claim boundaries
are in the tutorial `README.md`; they apply here unchanged.

## Contents

- `Allverify`: runs the selected cases in `work/<config>_<mesh>` and then the
  convergence analysis. `Allclean` removes `work`.
- `scripts/rt_reference.py`: the finite-amplitude reduced model (neo-Hookean
  solid on Gauss-Lobatto-Legendre (GLL) nodes, Chebyshev-collocated fluid, SciPy Radau).
- `scripts/make_reference.py`: writes `reference.json`, the reference at every
  solver time step (Python 3.12.3, numpy 1.26.4, scipy 1.11.4; about 2 s).
- `scripts/makeCase.py`: writes the whole case for a mesh, a configuration, a
  rank count, `nNonOrthogonalCorrectors` and an absolute coupling tolerance.
  The default tolerance is per mesh (1e-12 for M1 and M2, 1e-10 for M3).
- `thresholds.json`: every threshold with its provenance, in four sections (`layerA`, `layerB_deterministic`,
  `numerical_noise`, `family_consistency`). The noise allowances are never used to widen a band.
- `scripts/verify.py WORK [--log LOG_OR_CASE ...]`: the driver (`--log` takes a case directory or a log file, absolute or relative; the run label is sanitised so that every output file stays under `--out`; a missing log is INCOMPLETE, exit code 1, no traceback). Layer A (numerical acceptance, `scripts/acceptance.py`)
  and layer B (scientific reference-error expectations, after the family-consistency validation of
  `scripts/family.py`), explicit statuses, JSON and CSV output, the NOT-VERIFIED list in every run,
  exit code 1 on layer-A FAIL or INCOMPLETE.
- `test_verify.py`: tests without CFD (synthetic logs, family-consistency logic, `Allrun` stored/override paths with
  stubbed OpenFOAM helpers, corrected criteria, hygiene; known good and bad logs when `RT_KNOWN_RUNS` is set).
  Run with `python3 -B -I` so that no bytecode is written into the tree.
- `scripts/status.py`: the single ordered definition of the status severities (an unrecognised status is UNKNOWN, non-PASS).
- `scripts/family.py`: the family-consistency validation, semantic (parsed dictionaries; comments, whitespace, number formats and entry order do not matter). It checks physics, drive tables, the mesh of both regions (cell counts, refinement ratio, all vertex and arc coordinates, block topology, grading, patches), the whole `controlDict` before and after the `functions` block (except `endTime`, `writeInterval`) with the function-object parameters and code versions, PIMPLE and solver settings, the coupling algorithm, transfer method and predictors, the solid model and law, the dynamic mesh, the boundary conditions including the coded conditions, and a semantic fingerprint of every small dictionary; it does not check the solver binary, the generated `polyMesh` points or the logs. Documented-equivalent differences give CONDITIONAL, anything else INVALID-INPUT, unreadable input INCOMPLETE.
- `scripts/analyse_case.py CASE`: the quantities of interest, errors, balances and
  coupling statistics for one case (`analysis.json`, `series.csv`).
- `scripts/convergence.py WORK`: errors (max-norm and L2-in-time), observed and self-convergence
  orders and mixed terms for the cases in `WORK` (`convergence.json`, `convergence.csv`); used by `verify.py`.

## Consistency check of the stored dictionaries

`./Allclean` regenerates the stored dictionaries with `makeCase.py --mesh M1
--config combined`, so a clean tutorial is always the M1 combined 4-rank set. To
check that a working tree matches the generator, run
`python3 scripts/makeCase.py $TMPDIR/chk --mesh M1 --config combined --no-reference`
and `diff -r` the `0`, `constant` and `system` directories against the tutorial.

## What is not included

The field-level audits (continuity location, contraction factor, mode growth,
mesh quality, swept-volume GCL) and the temporal studies were run by project
scripts on the written fields and are not part of this tutorial. The ALE/GCL
consecutive-save checks are incomplete.
