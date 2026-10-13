---
sort: 1
---

# Optional fixes for the OpenFOAM installation

As described below, there are different optional fixes depending on whether you
are using OpenFOAM.org, OpenFOAM.com, or foam-extend. The replacement files are
in the `OpenFOAM-9`, `OpenFOAM.com` and `foam-extend-4.1` directories. Run
`./Allcheck` (with `export S4F_NO_FILE_FIXES=0`) to list the files that differ
from your installation, or `./Allreplace` to copy them into `$FOAM_SRC`, then
rebuild the affected libraries.

The replacement files are copies of OpenFOAM and foam-extend files and keep
their original copyright and licence headers. A "Modified by solids4foam" note
in each header describes the change.

## OpenFOAM.org

- `GeometricField.C`: this fix is required for consistent time discretisation.
  Some FSI cases may crash without this fix.

- `backwardDdtScheme.C`: this corrects the scheme for a moving mesh, for
  example, for the fluid domain in a fluid-solid interaction simulation.

## OpenFOAM.com

One copy of each file in `OpenFOAM.com` is used for all supported OpenFOAM.com
versions (v2312 to v2606).

- `backwardDdtScheme.C`: this corrects the scheme for a moving mesh, for
  example, for the fluid domain in a fluid-solid interaction simulation. The
  scheme falls back to Euler while a field has no old-old time level, instead
  of only on the first time step of the run. Apart from this change, the file
  is the same as the upstream file, which is identical in v2312 to v2606.

## foam-extend

- `pointBoundaryMesh.C`: without this fix, cases involving topological mesh
  changes will have a segmentation fault. For example, when using
  `crackerFvMesh`.

- `EulerDdtScheme.C`: This corrects the scheme for a moving mesh, for example,
  for the fluid domain in a fluid-solid interaction simulation.

- `backwardDdtScheme.C`: this corrects the scheme for a moving mesh, for
  example, for the fluid domain in a fluid-solid interaction simulation.

- `GeometricField.C`: this fix is required for consistent time discretisation.
  Some FSI cases may crash without this fix.
