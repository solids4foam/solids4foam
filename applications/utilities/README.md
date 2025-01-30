---
sort: 4
---

# Utilities

---

Prepared by Ivan Batistić with additions from Philip Cardiff

---

## Section Aims

- This document describes the `solids4foam` utilities located in `applications/utilities`;
- Utilities are executables that provide various pre and post-processing functionalities;
- These `solids4foam` utilities provide functionalities not available in the standard `OpenFOAM` utilities ([OpenFOAM-wiki](https://openfoamwiki.net/index.php/Main_OFUtilities)).

---

## `abaqusMeshToFoam`

- __Utility purpose__
  Mesh converter: converts [Abaqus](https://www.3ds.com/products-services/simulia/products/abaqus/) mesh (in `*.inp` format) into the FOAM mesh format.
  Currently, this utility only supports 3-D hexahedral cells/elements.
  Details regarding the FOAM mesh format can be found, for example, [here](https://www.openfoam.com/documentation/user-guide/4-mesh-generation-and-conversion/4.1-mesh-description#:~:text=By%20default%20OpenFOAM%20defines%20a,any%20restriction%20on%20its%20alignment.).
  Note that each distribution of `OpenFOAM` comes with a set of mesh converters (but not for Abaqus), see the available one [here](https://www.openfoam.com/documentation/user-guide/4-mesh-generation-and-conversion/4.5-mesh-conversion).

- __Arguments__

  - `<mesh.inp>` name of the Abaqus mesh file.

- __Options/parameters__

  None

- __Example of usage__
  ```bash
  $ abaqusMeshToFoam mesh.inp
  ```

  ```note
  - Only the following Abaqus element types are supported: C3D8 and C3D8R.
  - Only the first PART is used and the rest are ignored.
  - Node sets, element sets and surfaces are not converted.
  ```

---

## `addTinyPatch`

- __Utility purpose__
  For a chosen patch, find the closest face to the specified location and separate it into a new patch.
  The utility can be used, for example, to create patches for specifying point loads.

- __Arguments__

  - `<currentPatchName>` chosen patch name;
  - `<newTinyPatchName>` name of the one-face patch to be created;
  - `"(x y z)"` location vector.

- __Options/parameters__

  None

- __Example of usage__
  ```bash
  $ addTinyPatch Top TopNew "(30 30 0)"
  ```

<div style="text-align: center;">
<img src="./images/addTinyPatch.png" alt="Image" width="800">
    <figcaption>
     <strong>addTinyPatch: top patch before and after adding one-face patch</strong>
    </figcaption>
</div>

```note
When using `addTinyPatch` the original mesh is overwritten!
```

---

## `foamMeshToAbaqus`

- __Utility purpose__
  Mesh converter: converts FOAM mesh into [Abaqus](https://www.3ds.com/products-services/simulia/products/abaqus/) mesh (`*.inp` format).
  Currently, this utility only supports 3-D hexahedral cells/elements.

- __Arguments__

  None

- __Options/parameters__

  None

- __Example of usage__

  ```bash
  $ foamMeshToAbaqus
  ```

  Converted mesh is written to the`abaqusMesh.inp` file.
  Creates a node set and and element set and a surface for each boundary  patch.
  Also creates a element set for each material in the materials file (if it is exists).

  ```note
  - Only works for hexahedral cells as yet.
  - Created for Abaqus-6.9-2, but should work for new versions too.
  ```

---

## `perturbMeshPoints`

- __Utility purpose__
  Add a random perturbation to each interior mesh point. Boundary points are not perturbed, except for in-plane motion on `empty` and `wedge` patches.
  The utility can be used to create a distorted mesh to test the behavior (accuracy, order of accuracy, stability, etc.) of a discretisation procedure.

- __Arguments__
  None

- __Options/parameters__

  None

- __Dictionary__
  Inputs are defined in dictionary named `perturbMeshPointsDict` and located in `system` directory:

  ```c++
  seed        1;

  scaleFactor (5e-3 5e-3 5e-3);

  Gaussian    no;
  ```

  - `seed` is a scalar value used for the random number generator;
  - `scaleFactor` is a scaling vector which scales the point perturbation in each direction separately;
  - `Gaussian` enforces Gaussian distribution when perturbing points, otherwise uniform distribution is expected. Only used in combination with [OpenFOAM.com](https://www.openfoam.com/).

- __Example of usage__

  ```bash
  $ perturbMeshPoints
  ```

  The figure below shows mesh before and after using `perturbMeshPoints` utility.

<div style="text-align: center;">
  <img src="./images/perturbMeshPoints.png" alt="Image" width="800">
    <figcaption>
     <strong>perturbMeshPoints: mesh before and after point perturbation</strong>
    </figcaption>
</div>

```note
Perturbed mesh (`polyMesh`) is stored in the `0` directory and needs to be moved to `constant` before running the simulation!
```

```tip
For 2-D simulations, there is no need to perturb points in the `empty` direction. For an empty direction, zero scaling should be used, e.g.:
`scaleFactor (5e-3 5e-3 0);`
```

---

## `projectPatchToSphere`

- __Utility purpose__
  Project points of the specified patch onto the surface of a sphere.
  It does not alter internal points.

- __Arguments__

  - `<patchName>` patch to be projected;
  - `(x y z)` origin vector of the sphere.
  - `radius` radius of the sphere.

- __Options/parameters__

  None

- __Example of usage__

  ```bash
  $ projectPatchToSphere leftPatch (0 0 0) 0.5
  ```

```note
The mesh is overwritten!
```

---

## `splitPatch`

- __Utility purpose__
  Splits a patch into two patches by putting faces in the given bounding box in a new patch.

- __Arguments__
  None

- __Options/parameters__

  `-overwrite`  overwrite the original mesh when storing mesh after patch splitting.

- __Dictionary__
  ```c++
  patchToSplitName    Top;

  newPatchName        TopNew;

  boundBoxes
  (
      (0 0 0.01) (1 1 1)
  );
  ```

  - `patchToSplitName` is a name of the patch to be split;
  - `newPatchName`is the name of the splitted patch part;
  - `boundBoxes` is list of bounding boxes; each defined with two vectors:  `(xmin ymin zmin)` and `(xmax ymax zmax)`. Boundary patch faces inside this bounding box will be put in a new patch.
    _Note:_ The face's centre point is tested to see if it is inside the bounding boxes!

- __Example of usage__

  ```bash
  $ splitPatch
  ```

  <div style="text-align: center;">
    <img src="./images/splitPatch.png" alt="Image" width="800">
      <figcaption>
       <strong>splitPatch: top patch before and after split</strong>
      </figcaption>
  </div>

  ```note
  The mesh with the new patch is stored in a new time-step directory (`1` for example) and should be moved to the `constant` before running the simulation. Alternatively, the `-overwrite` option can be used to overwrite the original mesh.
  ```

---
