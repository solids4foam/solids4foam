---
sort: 10
---

# coupledUnsLinGeomLinearElasticSolid

This page documents the block-coupled linear-elastic solid model. The runtime
type is:

```text
coupledUnsLinearGeometryLinearElastic
```

The model assumes small strains and small rotations and solves the momentum
equation for the total displacement `D` with a **block-coupled** methodology:
all three displacement components are solved simultaneously in one linear
system, rather than segregated component by component. Because linear
elasticity with linear boundary conditions gives a linear problem, no outer
iteration is required at all — one linear solve per time step is the answer.

The method is described in:

P. Cardiff, Z. Tukovic, H. Jasak and A. Ivankovic (2016), _A block-coupled
finite volume methodology for linear elasticity and unstructured meshes_,
Computers and Structures, 175, 100-122,
[10.1016/j.compstruc.2016.07.004](https://doi.org/10.1016/j.compstruc.2016.07.004).

```warning
This model is **only available with foam-extend**. On OpenFOAM.com and
OpenFOAM.org builds the class compiles, but `evolve()` calls `notImplemented`
and the run aborts. It relies on foam-extend's `BlockLduMatrix` machinery and
on `src/blockCoupledSolids4FoamTools`.
```

---

## User Guide

### What this model solves

`coupledUnsLinGeomLinearElasticSolid`:

- assembles the linear-geometry momentum equation directly as a block matrix,
  with the Laplacian, Laplacian-transpose and Laplacian-trace terms of
  `div(sigma)` all treated implicitly;
- solves for total displacement `D` in a single block solve per time step;
- uses the `uns` face-Gauss gradient, so the discretisation matches
  [unsLinGeomSolid](https://www.solids4foam.com/documentation/solid-models/unsLinGeomSolid.html);
- reads `mu` and `lambda` once at construction from the mechanical law.

The advantage over the segregated models is convergence: the inter-component
coupling that segregated solvers defer to the outer iterations is here inside
the matrix, which matters most on high-aspect-ratio and bending-dominated
problems. The cost is memory and a larger, denser linear system.

### Restrictions

The constructor enforces two hard restrictions:

- exactly **one** material — multi-material cases abort;
- the mechanical law must be `linearElastic` — any other law aborts.

The linear definition of stress is hard-coded into the block assembly, so
these are not incidental limits. Equation under-relaxation is also disabled:
if `fvSolution` relaxes `DEqn`, the run aborts.

### Supported solution algorithms

There is a single block-coupled path and no outer loop. This model does not
read `solutionAlgorithm`, `nCorrectors` or the convergence tolerances.

### Model options

`coupledUnsLinGeomLinearElasticSolid` reads no entries of its own, and the
inherited iteration controls do not apply. The coefficients sub-dictionary
must still be present, but it may be empty:

```text
solidModel       coupledUnsLinearGeometryLinearElastic;

coupledUnsLinearGeometryLinearElasticCoeffs
{}
```

### Linear solver setup

The block system is solved with a `blockD` solver selected from `fvSolution`,
_not_ with an ordinary `D` solver:

```text
solvers
{
    blockD
    {
        solver              EigenSparseLU;
    }
}
```

The block solvers available come from
`src/blockCoupledSolids4FoamTools/BlockLduSolvers`.

### Boundary conditions

Standard solids4foam patch types do not work here: the boundary conditions are
inserted directly into the block matrix, so they need block-aware
implementations. Use the `block*` patch types from
`src/blockCoupledSolids4FoamTools/blockFvPatchVectorFields`:

| Patch type | Purpose |
| --- | --- |
| `blockFixedDisplacement` | Prescribed displacement |
| `blockFixedDisplacementZeroShear` | Normal displacement fixed, zero shear |
| `blockFixedGradient` | Prescribed displacement gradient |
| `blockSolidTraction` | Prescribed traction and pressure |
| `blockSolidVelocity` | Prescribed velocity |
| `blockGlobalSolid` | Global patch, used for coupling |

A symmetry plane is handled by `blockFixedDisplacementZeroShear`; there is no
block-aware `symmetry` type.

### Required input files

- `D` in the time directory, with `block*` boundary conditions;
- `constant/solidProperties`;
- `constant/mechanicalProperties`, selecting a single `linearElastic` law;
- `constant/g`;
- a `blockD` entry in `fvSolution`.

### Field glossary

- `D`, `DD`: total and incremental displacement.
- `pointD`, `pointDD`: displacement at mesh points.
- `grad(D)`, `grad(DD)`: displacement gradients.
- `sigma`: symmetric stress tensor, evaluated once per time step from the
  mechanical law.
- `U`: velocity field, computed as `ddt(D)`.
- `solutionVec`: the block solution vector, written with `AUTO_WRITE`. It has
  one entry per _variable_ of the extended mesh, which includes boundary
  unknowns as well as cells, so it is longer than the cell count.

---

## Developer Notes

### Class role

`coupledUnsLinGeomLinearElasticSolid` inherits directly from `solidModel`, and
holds a `solidPolyMesh` (`extendedMesh_`) alongside the ordinary `fvMesh`. The
extended mesh carries the larger addressing needed by the block system, in
which boundary faces as well as cells are unknowns.

- `D` is the primary solution variable;
- `nonLinGeom()` returns `LINEAR_GEOMETRY`;
- `muf_` and `lambdaf_` are surface fields, set once at construction, which is
  why only `linearElastic` is permitted.

Everything inside `evolve()` is guarded by `#ifdef FOAMEXTEND`.

### Construction

The constructor writes a banner citing the reference paper, calls
`DisRequired()`, then validates the mechanical model: one law, and that law
`linearElastic`. It casts the law and copies `mu()` and `lambda()` into
`muf_` and `lambdaf_`. On foam-extend it also builds `extendedMesh_` and the
`solutionVec_` field, sized from `extendedMesh_.nVariables()`.

### `evolve()`

There is no loop. Each time step:

1. clear the extended mesh's cached global coefficients;
2. build an empty `BlockLduMatrix<vector>` and zero its diagonal and
   off-diagonals;
3. assemble three separate block matrices through `BlockFvm::laplacian`,
   `BlockFvm::laplacianTranspose` and `BlockFvm::laplacianTrace`, and add
   their diagonal, off-diagonal and processor-coupling contributions together
   by hand — a `BlockFvMatrix` class would tidy this up;
4. insert the boundary-condition equations with
   `extendedMesh_.insertBoundaryConditions()`;
5. add the transient and gravity terms via `extendedMesh_.addFvMatrix()` with
   `rho*d2dt2(D) - rho*g`;
6. abort if `DEqn` under-relaxation is requested;
7. solve with the `blockD` solver and copy the result back into `D` through
   `extendedMesh_.copySolutionVector()`;
8. under-relax `D`, interpolate to `pointD`, and evaluate `grad(D)` — twice,
   because fixed-displacement boundaries need the patch internal field to set
   `snGrad` on the gradient's boundary;
9. update `DD`, `grad(DD)`, `sigma`, `pointDD` and `U`.

The non-orthogonal correction of the Laplacian is treated implicitly, which is
one of the advantages the block approach has over the segregated models.

### Extension points

- the three `BlockFvm` terms if a different stress form is wanted; note that
  removing the `linearElastic` restriction means these must become functions
  of the mechanical law rather than of fixed `mu` and `lambda` fields;
- `solidPolyMesh::insertBoundaryConditions()` for a new block patch type.

There is commented-out code in `evolve()` that would enforce zero normal
displacement on symmetry patches at the point level; the note there suggests a
`blockSymmetry` boundary condition as the proper fix.

---

## Tutorials

No tutorial selects `coupledUnsLinearGeometryLinearElastic` by default, since
it needs foam-extend. Three cases ship a ready-made `unsCoupled` variant —
`constant/solidProperties.unsCoupled`, `0/D.unsCoupled` and
`system/fvSolution.unsCoupled`:

- `solids/linearElasticity/cantilever2d`
- `solids/linearElasticity/narrowTmember`
- `solids/linearElasticity/ellipticPlate`
