---
sort: 1
---

# Optional fixes for the OpenFOAM installation

As described below, there are different optional fixes depending on whether you are using OpenFOAM.org, OpenFOAM.com, or foam-extend.

## OpenFOAM.com

  * `GeometricField.C`: this fix is required for consistent time discretisation. Some FSI cases may crash without this fix.

  * `backwardDdtScheme.C`: this corrects the scheme for a moving mesh, for example, for the fluid domain in a fluid-solid interaction simulation.


## OpenFOAM.com

  * `backwardDdtScheme.C`: this corrects the scheme for a moving mesh, for example, for the fluid domain in a fluid-solid interaction simulation.


## foam-extend

  * `EulerDdtScheme.C`: This corrects the scheme for a moving mesh, for example, for the fluid domain in a fluid-solid interaction simulation.

  * `GeometricField.C`: this fix is required for consistent time discretisation. Some FSI cases may crash without this fix.

  * `meshObjectBase.H`: without this fix, all runs will end in a segmentation. The solids4foam solver will work correctly; however, you may like to fix this if you plan to catch the return valve from the solver.

  * `pointBoundaryMesh.C`: without this fix, cases involving topological mesh changes will have a segmentation fault. For example, when using `crackerFvMesh`.

