---
sort: 2
---

# Cracking plate with a hole: `crackingPlateHole`

## Tutorial Aims

- Demonstrate crack propagation along a symmetry plane using a simple cohesive
  zone approach, which does not change the mesh topology.
- Show the use of the `simpleCrackerFvMesh` dynamic mesh with the
  `simpleCohesiveZone` boundary condition.

```note
This case only runs with foam-extend. With other OpenFOAM versions, the
`Allrun` script exits without running the case.
```

## Case Overview

The mesh is the same as in the linear elastic
[`plateHole`](../../linearElasticity/plateHole/README.md) tutorial: a `2 m` by
`2 m` region with a quarter of a circular hole of radius `0.5 m` at the
origin. The case is two-dimensional (plane stress). The plate is linear elastic
with a Young's modulus of `200 GPa` and a Poisson's ratio of `0.3`.

The `left` edge ($$x = 0$$) is displaced in the positive $$y$$ direction,
reaching `4e-5 m` at the end time of `20 s`, while the `right` edge
($$x = 2$$ m) is fixed. The `up` edge and the `hole` surface are
traction-free. The case is run for `20` time steps of `1 s`.

The `down` edge ($$y = 0$$) is a symmetry plane along which the crack grows.
Unlike the [`crackingBeams`](../crackingBeams/README.md) tutorial, no faces are
added to or removed from the mesh:

- The `simpleCohesiveZone` boundary condition on the `down` patch in `0/D`
  initially acts as a symmetry plane. Once the normal traction on a face
  reaches `sigmaMax`, the face is released and the cohesive traction from the
  chosen law is applied instead.
- The `Dugdale` cohesive law is used, with a maximum traction of `1 MPa` and a
  mode-I fracture energy of `10 J/m2`.
- The `simpleCrackerFvMesh` dynamic mesh, selected in
  `constant/dynamicMeshDict`, tracks and reports the released faces.

## Expected Results

Faces on the `down` patch are released one after another as the loading
increases, starting from the hole, where the stress is concentrated. The
reaction force on the `left` edge rises with the applied displacement, with
drops when the crack grows quickly enough to shed load. With foam-extend-4.1,
`30` faces are released over the run, and the reaction force first falls at
`t = 8 s`.

The force history is written to the `postProcessing` directory.

## Running the Case

The tutorial case is located at
`solids4foam/tutorials/solids/fracture/crackingPlateHole`. Run the included
script from this directory, with foam-extend loaded:

```bash
./Allrun
```

The script creates the mesh with `blockMesh` and runs the `solids4Foam` solver.
The case can be cleaned with `./Allclean`.
