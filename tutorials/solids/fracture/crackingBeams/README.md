---
sort: 1
---

# Cracking beams: `crackingBeams`

## Tutorial Aims

- Demonstrate crack propagation using a cohesive zone model, where the crack
  grows along internal mesh faces.
- Show the use of the `crackerFvMesh` dynamic mesh, which changes the mesh
  topology as the crack advances.

```note
This case only runs with foam-extend. With other OpenFOAM versions, the
`Allrun` script exits without running the case.
```

## Case Overview

Two beams, each `60 mm` long and `10 mm` deep, are joined along their common
mid-plane, apart from an initial crack over the first `10 mm` at the left end.
This resembles a double cantilever beam test. The case is two-dimensional
(plane stress), with a thickness of `1 mm`.

The left end of the upper beam (`topLoading`) is displaced upwards and the
left end of the lower beam (`bottomLoading`) is displaced downwards, each at a
rate of `1e-6 m/s`, pulling the beams apart. The remaining outer surfaces are
traction-free. The beams are linear elastic with a Young's modulus of
`200 GPa` and a Poisson's ratio of `0.3`. The case is run for `20` time steps
of `1 s`.

The crack is handled by two components:

- The `crackerFvMesh` dynamic mesh, selected in `constant/dynamicMeshDict`,
  converts internal faces into pairs of boundary faces on the initially empty
  `crack` patch when the `cohesiveZoneInitiation` law is satisfied. The
  `crackPathLimiter` restricts the crack to the mid-plane between the beams.
- The `solidCohesive` boundary condition on the `crack` patch in `0/D` applies
  the cohesive tractions between the two crack flanks. The `variableMixedMode`
  cohesive zone model is used, with a maximum normal and shear traction of
  `10 MPa` and mode-I and mode-II fracture energies of `50 J/m2`.

## Expected Results

The reaction force on the loaded patches initially rises elastically, reaches
a peak as the cohesive zone at the crack tip first fails, and then falls as the
crack runs along the mid-plane and the beams become more compliant. With
foam-extend-4.1, the peak vertical force on `topLoading` is approximately
`81.6 N` at `t = 6 s`, and `22` internal faces have been broken by
`t = 10 s`.

The force history is written to the `postProcessing` directory, and the crack
can be seen in ParaView by viewing the `crack` patch.

## Running the Case

The tutorial case is located at
`solids4foam/tutorials/solids/fracture/crackingBeams`. Run the included script
from this directory, with foam-extend loaded:

```bash
./Allrun
```

The script creates the mesh with `blockMesh` and runs the `solids4Foam` solver.
The case can be cleaned with `./Allclean`.
