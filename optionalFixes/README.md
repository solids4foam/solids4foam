---
sort: 1
---

# Optional fixes for the OpenFOAM installation

As described below, there are different optional fixes for OpenFOAM.org and
foam-extend; none are required for OpenFOAM.com. The replacement files are in
the `OpenFOAM-9` and `foam-extend-4.1` directories. Run `./Allcheck` (with
`export S4F_NO_FILE_FIXES=0`) to list the files that differ from your
installation, or `./Allreplace` to copy them into `$FOAM_SRC`, then rebuild the
affected libraries.

The replacement files are copies of OpenFOAM and foam-extend files and keep
their original copyright and licence headers. A "Modified by solids4foam" note
in each header describes the change.

Earlier versions also replaced the `EulerDdtScheme.C` and `backwardDdtScheme.C`
ddt schemes. These are no longer needed: the foam-extend `pimpleFluid` computes
its `ddtCorr` flux correction itself. If you replaced these files, the
foam-extend files should be restored, as the replacement `ddtCorr` does not
suit the foam-extend solvers, such as `pimpleDyMFoam`.

## OpenFOAM.org

- `GeometricField.C`: this fix is required for consistent time discretisation.
  Some FSI cases may crash without this fix.

## foam-extend

- `pointBoundaryMesh.C`: without this fix, cases involving topological mesh
  changes will have a segmentation fault. For example, when using
  `crackerFvMesh`.

- `GeometricField.C`: this fix is required for consistent time discretisation.
  Some FSI cases may crash without this fix.
