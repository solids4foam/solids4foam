---
sort: 15
---

# kirchhoffPlateSolid

This page documents the Kirchhoff thin-plate solid model. The runtime type is:

```text
kirchhoffPlate
```

Unlike every other solid model in solids4foam, this one does not discretise a
three-dimensional continuum. It solves the Kirchhoff plate equations on a
**finite area** mesh — a two-dimensional surface mesh — for the transverse
deflection `w` and the moment sum `M`. The formulation is the rotation-free
one of Torlak (2006), which splits the fourth-order plate equation into two
coupled second-order equations so that standard Laplacian discretisations can
be used.

```warning
This model requires a **finite area** mesh and OpenFOAM.com **v2412 or newer**.
The shipped tutorial explicitly excludes OpenFOAM.org and foam-extend. The
class itself is compiled for all forks, but the finite-area machinery it needs
is not available everywhere.
```

---

## User Guide

### What this model solves

The plate equations are split as:

```text
M equation:  rho*h*d2dt2(w) = laplacian(M) + p
w equation:  laplacian(bendingStiffness, w) + M = 0
```

where `w` is the transverse (out-of-plane) deflection, `M` is the moment sum,
`p` is the net transverse pressure, `h` is the plate thickness, `rho` is the
density, and the bending stiffness is

```text
bendingStiffness = E*h^3/(12*(1 - nu^2))
```

The two equations are solved segregated, one after the other, inside an outer
loop. The rotation field `theta = -grad(w)` is updated each iteration and used
by the clamped boundary condition.

### Restrictions

The constructor enforces:

- exactly **one** material — multi-material cases abort;
- the mechanical law must be `linearElastic`, from which `rho`, `E` and `nu`
  are read once at construction.

Beyond that:

- `tractionBoundarySnGrad()` is `notImplemented`, so ordinary solid traction
  boundary conditions cannot be used;
- the stress field `sigma` is never calculated, and its write option is
  switched to `NO_WRITE` at the end of `evolve()`;
- the volume mesh must be exactly one cell thick, with two opposite patches of
  equal face count. The constructor locates the finite-area patch and its
  shadow patch automatically, and aborts if the mesh does not satisfy this.

### Supported solution algorithms

A single segregated implicit loop over the two scalar equations. There is a
commented-out sketch in the source for solving them block-coupled, which is
not implemented. This model does not read `solutionAlgorithm`.

### Model options

| Entry | Dimensions | Description |
| --- | --- | --- |
| `plateThickness` | `[0 1 0 0 0 0 0]` | Plate thickness `h`; required |

`plateThickness` is read with `lookup()`, so the model fails at construction if
it is missing.

The relevant inherited `solidModel` entries are:

| Entry | Default | Relevance |
| --- | --- | --- |
| `nCorrectors` | `10000` | Maximum number of outer correctors |
| `solutionTolerance` | `1e-06` | Primary convergence tolerance, `w` and `M` |
| `alternativeTolerance` | `1e-07` | Secondary convergence tolerance |
| `infoFrequency` | `100` | Frequency for solver progress output |

Tight tolerances are usual here; the shipped tutorial uses `1e-10` for both.

### Required input files

Because the primary fields live on the finite-area mesh, the case layout is
different from the other solid models:

- `0/finite-area/w`, `0/finite-area/M` and `0/finite-area/p`, all `MUST_READ`;
- `0/finite-area/D`, normally `calculated`, since it is mapped from `w`;
- `system/finite-area/faMeshDefinition`, `faSchemes` and `faSolution`;
- `constant/solidProperties`;
- `constant/mechanicalProperties`, selecting a single `linearElastic` law;
- `constant/g` and `constant/dynamicMeshDict`.

The finite-area mesh is generated with `makeFaMesh` before running
`solids4Foam`.

### Recommended dictionary setup

Minimal example for `constant/solidProperties`:

```text
solidModel     kirchhoffPlate;

kirchhoffPlateCoeffs
{
    plateThickness       plateThickness [0 1 0 0 0 0 0] 0.1;

    nCorrectors          1000;
    solutionTolerance    1e-10;
    alternativeTolerance 1e-10;
}
```

`faSchemes` and `faSolution` provide the schemes and linear solvers for `w`
and `M`; the ordinary `fvSchemes` and `fvSolution` are still read but the plate
equations do not use them.

### Boundary conditions

The boundary conditions are finite-area patch fields on `w` and `M`, not the
usual `fvPatchField` types. Two are provided by this model:

| Patch type | Field | Purpose |
| --- | --- | --- |
| `clampedMoment` | `M` | Clamped edge; sets the moment from the rotation |
| `freeEdgeDisplacement` | `w` | Free edge |

A simply supported edge is `fixedValue` on `M` with `fixedValue` on `w`; a
symmetry edge is `zeroGradient` on both. The transverse pressure `p` is
normally `zeroGradient` on all edges, with its internal field setting the
load.

### Field glossary

Finite-area fields, on the surface mesh:

- `w`: transverse deflection; primary unknown.
- `M`: moment sum; primary unknown.
- `p`: net transverse pressure load.
- `theta`: rotation, `-grad(w)`.

Volume fields, mapped from the above onto the single-layer volume mesh so that
ordinary post-processing tools can read them:

- `wVf`, `MVf`, `pVf`, `thetaVf`: the mapped counterparts.
- `D`: displacement, mapped from `w*faceAreaNormals`, i.e. the deflection
  applied along the plate normal.
- `pointD`, `pointDD`, `DD`, `U`: derived from `D` as usual.

`sigma` exists but is never calculated and is not written.

---

## Developer Notes

### Class role

`kirchhoffPlateSolid` inherits from `solidModel` and holds a `faMesh`
(`aMesh_`) built from the volume mesh. It satisfies the `solidModel` interface
by mapping its finite-area solution onto the volume fields the rest of
solids4foam expects.

- `nonLinGeom()` returns `LINEAR_GEOMETRY`;
- the primary unknowns are `w_` and `M_`, not `D`;
- `D` is an output, reconstructed from `w`.

### Construction

The constructor builds `aMesh_`, reads `w`, `M` and `p` on the area mesh
(`MUST_READ`) and creates their volume counterparts, creates `theta` and
`gradTheta`, and reads `plateThickness`. It then validates the mechanical
model — one law, and that law `linearElastic` — extracts `rho`, `E` and `nu`,
and forms `bendingStiffness_`.

`calcAreaPatches()` identifies the finite-area patch and, from the first cell,
the opposite "shadow" patch by finding the face whose normal is most
anti-parallel. It aborts if the two patches do not have equal face counts,
which is the check that the mesh is one cell thick.

### `evolve()`

An outer mesh-update loop wrapping an inner equation loop. Each inner
iteration:

1. assembles and solves the `M` equation. `fac::d2dt2` is not available, so
   the second time derivative is formed manually as
   `(fac::ddt(w) - fac::ddt(w.oldTime()))/deltaT`, and all terms are moved to
   the left-hand side because `==` is not supported for `faScalarMatrix` here;
2. relaxes and solves the `w` equation
   `fam::laplacian(bendingStiffness, w) + M`;
3. updates `theta = -fac::grad(w)` and `gradTheta`, the latter being used for
   the non-orthogonal correction in the clamped boundary condition.

After convergence, `M`, `w`, `theta` and `p` are mapped to their volume
counterparts with `mapAreaFieldToSingleLayerVolumeField()`, `D` is built as
`w*aMesh_.faceAreaNormals()` and mapped, and `pointD`, `DD`, `pointDD` and `U`
follow. Finally `sigma().writeOpt()` is set to `NO_WRITE`.

### Extension points

- the commented-out block-coupled branch in `evolve()`, which would solve `w`
  and `M` simultaneously as a `BlockLduSystem<vector2, vector2>`; this is
  sketched but not implemented;
- `setTraction()` and the finite-area patch fields for new edge conditions;
- `tractionBoundarySnGrad()`, currently `notImplemented`, if the model is ever
  to participate in fluid-solid interaction.

Extending this class to Mindlin-Reissner plates would mean carrying the
rotations as independent unknowns rather than deriving them from `grad(w)`.

---

## Tutorials

Cases that select `kirchhoffPlate`:

- `solids/beamsPlatesShells/squarePlate` — a clamped square plate under
  uniform pressure, compared against an analytical solution built by the
  case's own `squarePlateAnalyticalSolution` utility.
