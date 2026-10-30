# beamSolid

`beamSolid` runs the geometrically exact beam formulation of
[beamFoam](https://github.com/solids4foam/beamFoam) as a solids4foam solid
model. The beam is solved by beamFoam's `coupledTotalLagNewtonRaphsonBeam`
model; `beamSolid` connects it to the `solidModel` interface, so a beam case
runs with the `solids4Foam` application.

`beamSolid` and the beam models are built for OpenFOAM.com when Eigen is
available, that is when `S4F_NO_USE_EIGEN` is not set.

## Libraries

- `libS4FBeamModels` (`src/beamModels`): the beam models, beam boundary
  conditions, cross-sections and beam function objects, copied unchanged from
  beamFoam.
- `libS4FBeamSolid` (`src/beamSolid`): the `beamSolid` solid model.

The beam mesh utilities `createBeamMesh` and `setInitialPositionBeam` are also
copied unchanged from beamFoam, in `applications/utilities`.

## Case Setup

A beamFoam case runs with `solids4Foam` after three additions:

- `system/controlDict`: `application solids4Foam;` and
  `libs ("libS4FBeamSolid.so");`
- `constant/physicsProperties`: `type solid;`
- `constant/solidProperties`:

  ```c++
  solidModel beamSolid;

  beamSolidCoeffs
  {
  }
  ```

The solid model base class also needs `constant/g`; use zero for no gravity.

Everything else stays as in beamFoam. The beam model, its solver settings, the
cross-sections and the material are read from `constant/beamProperties`, the
mesh is created with `createBeamMesh` (and positioned with
`setInitialPositionBeam` if needed), and the unknowns are the beam displacement
`W` and rotation `Theta`. The beam boundary conditions, for example
`forceBeamDisplacementNR` and `momentBeamRotationNR`, apply the loads.
Continuum traction boundary conditions are not supported.

See the `beamCantilever`, `beamDynamicCantilever` and `beamFollowerForce`
tutorials in `tutorials/solids/beamsPlatesShells`.

## Interface Fields

The beam model works on the mesh that the solid model reads, so the two share
one mesh. After each time step, `beamSolid` sets the `solidModel` displacement
fields from the beam:

- `D` from `W`, and `DD` from the change in `D`;
- `pointD` from `pointW`, which includes the cross-section rotation, and
  `pointDD` from the change in `pointD`.

These fields, and `sigma`, are not written: the beam writes its own fields
(`W`, `Theta`, `U`, the section forces and moments, and so on).

## Synchronising with beamFoam

The beam sources in `src/beamModels`, `applications/utilities/createBeamMesh`
and `applications/utilities/setInitialPositionBeam` are copies of beamFoam
files. The beamFoam commit they come from is recorded in
`src/beamModels/BEAMFOAM_VERSION`. Only the build files (`Make/` and
`Allwmake`) belong to solids4foam. Change the beam sources in beamFoam, then
copy them across:

```bash
src/beamModels/syncBeamFoam.sh <beamFoam repository> <commit>
```

The script reports changes to beamFoam's list of sources, which `Make/files`
must then follow.

## Limitations

- Only OpenFOAM.com has been tested (v2412), in serial, with one beam region.
- Fluid-beam coupling is not available yet.
- The optional beam momentum contributions (for example Morison drag and ground
  contact) are compiled but have not been tested through `beamSolid`.
- Continuing a run from a later time is stopped by the `solidModel` restart
  check unless `restart no;` is set in `beamSolidCoeffs`. The beam writes the
  state that it needs to continue, and a run continued with `restart no;` gives
  the same answer as an uninterrupted run, to the written precision.
