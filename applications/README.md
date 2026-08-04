---
sort: 8
---

# Applications

Everything that is compiled into an executable or is run from the shell lives
in `applications`. There is one solver, a small set of mesh and case
utilities, and a library of Bash functions used by the tutorials.

---

## Contents

| Item | What it is |
| ---- | ---------- |
| [`solids4Foam`](https://www.solids4foam.com/documentation/applications/solids4Foam.html) | The solver. Solid, fluid and fluid-solid interaction analyses, with the physics selected at run time. |
| [Utilities](https://www.solids4foam.com/documentation/applications/utilities.html) | Mesh and case pre- and post-processing: `abaqusMeshToFoam`, `addTinyPatch`, `foamMeshToAbaqus`, `perturbMeshPoints`, `projectPatchToSphere` and `splitPatch`. |
| [`decomposeParMonolithic`](https://www.solids4foam.com/documentation/applications/decomposeParMonolithic.html) | Decomposes the regions of a monolithic multi-region case as one graph, so that the partitioning is consistent across regions. |
| [`solids4FoamScripts.sh`](https://www.solids4foam.com/documentation/applications/scripts.html) | Bash functions used by the tutorial `Allrun` scripts, principally for making a case run on any OpenFOAM fork. |

---

## Building

`applications/Allwmake` builds the solver, the utilities, the tests and the
scripts; it is called in turn by the top-level `Allwmake`.
`applications/Allwclean` removes the resulting executables and dependencies.

```bash
cd applications
./Allwmake
```

---

## Source

[`applications`](https://github.com/solids4foam/solids4foam/tree/development/applications)
