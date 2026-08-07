---
sort: 3
---

# Material models

A material model describes how a material responds, as opposed to a solid
model, which describes how the governing equations are discretised and
solved. solids4foam splits this into two independent models:

- the **mechanical model** (`mechanicalModel`), which turns deformation into
  stress and provides the density and the implicit stiffness used by the
  momentum equation. It is read from `constant/mechanicalProperties` and every
  solid model creates one.
- the **thermal model** (`thermalModel`), which provides the specific heat
  capacity and thermal conductivity used by the heat equation. It is read from
  `constant/thermalProperties` and is only created by the solid models that
  solve for temperature.

The two are deliberately separate. A thermal analysis therefore needs both
dictionaries: the density that multiplies the specific heat capacity is taken
from the mechanical model, not the thermal one.

---

## User Guide

### The `constant/mechanicalProperties` dictionary

`mechanicalProperties` is read `MUST_READ`, so it must exist for every solid
case. It has these top-level entries:

| Entry | Required | Description |
| --- | --- | --- |
| `planeStress` | yes | `yes`/`no`; plane stress or plane strain/3-D |
| `mechanical` | yes | List of named material sub-dictionaries |
| `writeSubMeshes` | no | Default `no`; write the material sub-meshes |

`planeStress` is read with `lookup()` and so has no default; omitting it is a
fatal error even for a genuinely three-dimensional case, where `no` is the
correct answer. It is read once by `mechanicalModel` and is then available to
every law through the protected `planeStress()` accessor, which searches the
current region first and then the `region0` or `solid` sub-registry — this is
what makes it work in fluid-solid interaction cases where the solid is a named
region.

`writeSubMeshes` is only meaningful for multi-material cases; it is read with
`lookupOrDefault<Switch>` and, when enabled, writes each material sub-mesh to
disk so that the partitioning can be inspected.

### Selecting a mechanical law

The `mechanical` entry is a **list**, not a dictionary. Each element is a
named sub-dictionary, and the `type` keyword inside it selects the mechanical
law at runtime:

```text
planeStress     no;

mechanical
(
    steel
    {
        type            linearElastic;
        rho             rho [1 -3 0 0 0 0 0] 7854;
        E               E [1 -1 -2 0 0 0 0] 200e9;
        nu              nu [0 0 0 0 0 0 0] 0.3;
    }
);
```

For a single-material case the name (`steel` above) is arbitrary and is used
only in solver output. It becomes significant as soon as there is more than
one entry — see below.

There are **two separate run-time selection tables**, and a law is registered
in exactly one of them:

- `linGeomMechLaw`, used when the solid model reports
  `LINEAR_GEOMETRY`;
- `nonLinGeomMechLaw`, used when the solid model reports
  `TOTAL_LAGRANGIAN` or `UPDATED_LAGRANGIAN`.

A law can therefore not be paired with an incompatible solid model. If a
nonlinear-geometry law is requested from a linear-geometry solid model, the
selector detects this and reports that the law "can only be used with a
nonlinear geometry solid model" rather than simply saying the type is unknown.
The two catalogues are documented at:

- [Linear geometry mechanical laws](https://www.solids4foam.com/documentation/material-models/linear-geometry-laws.html)
- [Nonlinear geometry mechanical laws](https://www.solids4foam.com/documentation/material-models/nonlinear-geometry-laws.html)

### Entries common to every law

Four entries are handled by the `mechanicalLaw` base class and so are
available inside any law sub-dictionary, whichever `type` is selected:

| Entry | Default | Description |
| --- | --- | --- |
| `type` | none | Required; the law's run-time type name |
| `rho` | none | Density; required by laws using the base `rho()` |
| `solvePressureEqn` | `no` | Solve a smoothing Laplacian for the pressure |
| `pressureSmoothingScaleFactor` | `100` | Scales that smoothing |

`rho` is read as a `dimensionedScalar` by `mechanicalLaw::rhoScalar()`. A law
that overrides `rho()` need not provide it, but almost all do use the base
implementation. Note that `rho()` uses `READ_IF_PRESENT`, so a `rho` field
placed in the time directory will be read in preference to the uniform value.

`solvePressureEqn` and `pressureSmoothingScaleFactor` are read with
`lookupOrAddDefault`, which means the defaults are written back into the
in-memory dictionary and appear in the `mechanicalProperties.withDefaultValues`
file the solver writes at the end of a run. That file is a useful way to see
exactly what was used.

A fifth, optional entry is `regionName`. The base class needs to know which
mesh region is the "base" solid region; it uses `regionName` if present, else
`solid` if that region exists, else `region0`. Only unusual multi-region
set-ups need to set it explicitly.

### Multi-material cases

Give the `mechanical` list more than one entry and solids4foam switches to its
multi-material path:

```text
planeStress     yes;

mechanical
(
    outer
    {
        type            linearElastic;
        rho             rho [1 -3 0 0 0 0 0] 1000;
        E               E [1 -1 -2 0 0 0 0] 200e+9;
        nu              nu [0 0 0 0 0 0 0] 0.3;
    }
    inner
    {
        type            linearElastic;
        rho             rho [1 -3 0 0 0 0 0] 1000;
        E               E [1 -1 -2 0 0 0 0] 20e+9;
        nu              nu [0 0 0 0 0 0 0] 0.35;
    }
);
```

The rules that follow from the implementation are:

1. **The name of each entry must be the name of a `cellZone`.** The keyword of
   each `mechanical` entry (`outer`, `inner` above) is taken verbatim as a
   cellZone name. If no cellZone of that name exists, construction aborts with
   "cellZone not found for material \<name\>".
2. **Every cell must be in exactly one of those cellZones.** Before the
   sub-meshes are built, `checkCellZones()` counts how many of the listed
   zones each cell belongs to and aborts if any cell is in none of them
   ("There are cells that are not in a material cellZone!") or in more than
   one. Cells in other, unlisted cellZones are not exempt: the check is over
   the whole mesh.
3. **Each material gets its own sub-mesh.** `solidSubMeshes` builds one
   `fvMeshSubset` per material from the corresponding cellZone, and the
   mechanical law for that material is constructed **on the sub-mesh**, not on
   the base mesh. Every field the law creates therefore lives on the sub-mesh,
   and the results are mapped back to the base mesh by `mechanicalModel`.
4. **Interfaces are internal faces, not patches.** A face shared by two
   materials stays an internal face of the base mesh; in each sub-mesh it
   appears on the auto-generated `oldInternalFaces` patch. There is no
   coupling dictionary and no patch to set up.

Cell zones are normally created after meshing with `setSet` followed by
`setsToZones` (or with `topoSet`). The `layeredPipe` tutorial does exactly
this, with a `batch.setSet` file that makes an `outer` and an `inner` cellSet
and subtracts one from the other so the two do not overlap.

```note
Because the sub-mesh partition is by cellZone, a material region does not need
to be contiguous, and the base mesh does not need to be conformal in any
special way. It does, however, have to be a single mesh: separately meshed
parts must be merged first, as in the `punch` tutorial.
```

### What happens at a bi-material interface

Stress is discontinuous across an interface between dissimilar materials,
which is exactly what a smooth cell-to-face interpolation gets wrong. Rather
than interpolate the base displacement onto each sub-mesh naively,
solids4foam computes an interface displacement that drives the two sides
towards traction equilibrium.

For each interface face, the stiffness-weighted update is

```text
Di = Di.prevIter
   + wab*(tractionB - tractionA)
   + wa*(Da - Da.prevIter) + wb*(Db - Db.prevIter)
```

where `A` and `B` are the two sides, `da` and `db` are the normal distances
from the face to the two cell centres, `Ka` and `Kb` are the implicit
stiffnesses there, and

```text
wab = da*db/(db*Ka + da*Kb)
wa  = db*Ka/(db*Ka + da*Kb)
wb  = 1 - wa
```

The consequences for the user are:

- The correction is **iterative**: it uses the previous outer-iteration values,
  so traction continuity is achieved at convergence, not in a single
  iteration. Expect multi-material cases to need more outer correctors than
  the equivalent single-material case, particularly when the stiffness ratio
  is large.
- It uses the implicit stiffness `impK`, so the stiffness contrast between
  materials is what sets the weights. This is the mechanism that keeps the
  solution oscillation-free across a large jump in Young's modulus.
- For nonlinear-geometry laws the face normal is rotated with Nanson's
  formula, using the relative deformation gradient for updated Lagrangian
  solid models and the total one for total Lagrangian solid models.

If no face is shared by two materials — for example, materials separated by an
internal boundary — the correction is skipped entirely and each sub-mesh is
filled by plain interpolation.

### Multi-material limitations

Some capabilities are single-material only, and fail with a clear
`notImplemented` or fatal error rather than silently giving a wrong answer:

| Capability | Status with more than one material |
| --- | --- |
| `correct(pointSymmTensorField&, ...)` | Not implemented |
| Face-quadrature `correct(...)`/`grad(...)` | Not implemented |
| `crackerFvMesh` | Fatal error at construction (foam-extend) |

The vertex-centred solid models are not affected by the first two rows in the
same way: they use `dualMechanicalModel`, which handles multiple materials by
constructing every law on the whole dual mesh and masking each law by the
cells belonging to its material.

### Thermal properties

Solid models that solve a temperature equation additionally read
`constant/thermalProperties`, which selects a `thermalLaw` giving the specific
heat capacity and thermal conductivity. See the thermal model page in this
section.

---

## Developer Notes

### Class layout

`mechanicalModel` derives from both `IOdictionary` (it *is* the
`mechanicalProperties` dictionary) and `PtrList<mechanicalLaw>` (it *is* the
list of laws). The number of laws is what selects between the single- and
multi-material code paths throughout the class; there is no separate flag.

`solidSubMeshes` is created lazily, on the first call to `solSubMeshes()`, and
deliberately aborts if constructed for a single material. Note the ordering
trap handled in `clearOut()`: the laws must be cleared before the sub-meshes,
because the laws hold `GeometricField`s registered on the sub-meshes.

### Accumulation pattern

`rho()`, `impK()`, `bulkModulus()` and `shearModulus()` all follow the same
shape. With one law they forward straight to it; with several they build a
`PtrList` of the per-sub-mesh fields and call
`solidSubMeshes::mapSubMeshVolFields` to assemble the base-mesh field. The
same pattern with `mapSubMeshSurfaceFields` is used for `sigma`, with one
difference: the surface field is zeroed first, because interface faces receive
a contribution from both materials, whereas volume fields store nothing on the
interface.

`impKf()` is not accumulated from per-law surface fields; it is a plain linear
interpolation of the assembled `impK`, with the comment that this gives the
best convergence in practice.

### Gradients and vol-to-point interpolation

In the multi-material path, `grad()` does not take the gradient of the base
field. It interpolates the base `D` to each sub-mesh (applying the interface
correction), takes the gradient on each sub-mesh, calls `correctBoundarySnGrad`
to restore non-orthogonal correction on boundaries that the sub-setting turned
into `calculated` patches, maps back, and finally re-applies
`gaussGrad::correctBoundaryConditions`. The surface-field variant additionally
overwrites the normal component with `n*fvc::snGrad(D)`, which the code notes
is necessary for convergence in many cases.

`interpolate()` follows the same shape using one
`enhancedVolPointInterpolation` (`newLeastSquaresVolPointInterpolation` on
foam-extend) per sub-mesh, then maps the point fields back.
`correctSymmPlanes()` then zeroes the normal component of `pointD` on symmetry
and empty patches by testing the average patch normal against the coordinate
axes, and applies the solid model's 2-D corrector.

### Adding a new mechanical law

See the `mechanicalLaw` base-class page in this section for the interface a new
law must implement and the concrete steps to register it. The base class lives
at
[`mechanicalLaw.H`](https://github.com/solids4foam/solids4foam/blob/development/src/solids4FoamModels/materialModels/mechanicalModel/mechanicalLaws/mechanicalLaw/mechanicalLaw.H).

### Dictionary write-back

`writeDict()` is called from `solidModel::end()`. Because each law's
dictionary is a copy rather than a reference into `mechanicalProperties`, the
law dictionaries have to be reassembled into a `mechanical` entry before being
written, which is why the function exists at all. The result is
`constant/mechanicalProperties.withDefaultValues`, showing every entry
including the ones filled in by `lookupOrAddDefault`.

---

## Tutorials

Cases with more than one entry in the `mechanical` list:

- `solids/multiMaterial/layeredPipe`: a bi-material thick-walled cylinder;
  the reference case for the multi-material machinery, with cellZones built by
  `setSet`/`setsToZones` and an analytical solution to compare against.
- `solids/linearElasticity/punch`: two separately meshed parts merged into one
  mesh, then split into the `punch_top` and `punch_bottom` cellZones.
- `solids/thermoelasticity/hotCylinder/hotCylinderPredefinedTFieldMultipleMaterials`:
  `steel` and `aluminium` materials, each a `thermoMechanicalLaw` wrapping a
  `linearElastic` law.
