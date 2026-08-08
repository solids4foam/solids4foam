---
sort: 27
---

# electroMechanicalLaw

This law adds fibre-directed active stress to the stress calculated by a
run-time selectable passive nonlinear mechanical law. The runtime type is:

```text
electroMechanicalLaw
```

---

## User Guide

### What it computes

The cell-centred overload first asks the passive law to calculate the Cauchy
stress, then adds an active contribution:

```text
A0           = f0 (x) f0
J            = det(F)
sigma        = sigma_passive + symm(F & (Ta_current*A0) & F^T)/J
```

When a cell-centred `Ta` field exists in the mesh object registry at the first
call to `correct`, `Ta_current` is that field. Otherwise the dictionary value
is ramped according to:

```text
Ta_current = (t/rampTime)*activeTension    for t < rampTime
Ta_current = activeTension                 otherwise
```

The field-based choice is made once and retained. For the face-centred
overload, `Ta` is interpolated from cells to faces. The implemented face
relation differs from the cell relation by multiplying by `J`:

```text
A0f          = f0f (x) f0f
Jf           = det(Ff)
sigmaf       = sigmaf_passive
             + Jf*symm(Ff & (Taf_current*A0f) & Ff^T)
```

The `volSymmTensorField` and `surfaceSymmTensorField` overloads are
implemented. The inherited `pointSymmTensorField` overload aborts with
`notImplemented`. On OpenFOAM.com and OpenFOAM.org, the inherited
`CompactListList` overload also aborts with `notImplemented`; that overload is
not compiled for foam-extend.

### Model options

| Entry | Required | Description |
| --- | --- | --- |
| `rho` | yes | Density, `[1 -3 0 0 0 0 0]` |
| `activeTension` | yes | Constant active tension, `[1 -1 -2 0 0 0 0]` |
| `rampTime` | yes | Ramp duration in simulation-time units |
| `passiveMechanicalLaw` | yes | Passive nonlinear-law sub-dictionary |
| `passiveMechanicalLaw.type` | yes | Runtime type of the passive law |
| `solvePressureEqn` | no | Base-class switch; default `false` |
| `pressureSmoothingScaleFactor` | no | Base-class scalar; default `100` |
| `regionName` | no | Base mesh region used to find the solid model |

`activeTension` and `rampTime` remain required when a `Ta` field will be used,
because the constructor reads both before checking the registry. A negative
`rampTime` is rejected. A zero value applies the full constant tension without
division by zero because the ramp branch is then not entered for non-negative
simulation time.

The top-level base-class pressure options are constructed, but this wrapper
does not call `updateSigmaHyd()`. Options needed by the selected passive law
belong inside `passiveMechanicalLaw`. If `regionName` is absent, the base class
selects `solid`, then `region0`; construction fails if neither region exists.

The law also requires `volVectorField` `f0` and `surfaceVectorField` `f0f` in
the initial time directory. They are read independently and are not normalised
by this class.

### Recommended dictionary setup

This setup follows the `idealisedVentricle` tutorial and uses a Guccione law
for passive heart tissue:

```text
planeStress     no;

mechanical
(
    heartTissue
    {
        type            electroMechanicalLaw;
        rho             rho [1 -3 0 0 0 0 0] 3000;
        activeTension   activeTension [1 -1 -2 0 0 0 0] 0;
        rampTime        0;

        passiveMechanicalLaw
        {
            type            GuccioneElastic;
            k               k [1 -1 -2 0 0 0 0] 10e3;
            cf              1;
            ct              1;
            cfs             1;
            uniformFibreField yes;
            f0              (0 0 1);
            bulkModulus     bulkModulus [1 -1 -2 0 0 0 0] 1e6;
        }

        // Optional base-class entries
        // regionName                    solid;
        // solvePressureEqn              no;
        // pressureSmoothingScaleFactor  100;
    }
);
```

The outer law still requires separate `f0` and `f0f` field files even when the
passive law uses its own uniform fibre setting. Set `activeTension` to the
fallback value required when no coupling model registers `Ta`.

### Field glossary

- `f0`: reference fibre direction at cell centres, read with `MUST_READ` and
  written automatically.
- `f0f`: reference fibre direction at faces, read separately with `MUST_READ`
  and written automatically.
- `f0f0`, `f0f0f`: cell- and face-centred structural tensors formed by
  `sqr(f0)` and `sqr(f0f)`.
- `Ta`: optional cell-centred active-tension field supplied through the object
  registry; this law does not create it.
- `F_`, `Ff_`: deformation gradients accessed through the base class for the
  cell- and face-centred active-stress transformations.
- Fields belonging to the passive law are maintained by that nested law.

---

## Developer Notes

### Class role

`electroMechanicalLaw` derives directly from `mechanicalLaw` and is registered
in the `nonLinGeomMechLaw` selection table. It owns an `autoPtr<mechanicalLaw>`
for the passive law, the `f0`, `f0f`, `f0f0` and `f0f0f` fields, the constant
`Ta_`, the scalar `rampTime_`, and two mutable flags that cache the registry
check for field-based active tension.

There are no fork-specific guards in this class. Its source is listed in both
`Make/files.openfoam` and `Make/files.foamextend`. The base-class
`CompactListList` interface is guarded by `#ifndef FOAMEXTEND`.

### Construction

Construction proceeds in this order:

1. The `mechanicalLaw` base reads `solvePressureEqn` with default `false` and
   `pressureSmoothingScaleFactor` with default `100`, then resolves
   `regionName`.
2. The `passiveMechanicalLaw` sub-dictionary and its required `type` entry are
   used to construct a nonlinear mechanical law through the runtime table.
3. `f0` is read and `f0f0` is formed as `sqr(f0)`.
4. `f0f` is read independently and `f0f0f` is formed as `sqr(f0f)`.
5. The required `activeTension` dimensioned scalar and `rampTime` scalar are
   read, and the registry-check flags are initialised to `false`.
6. A negative `rampTime` triggers `FatalError` with the message that it should
   be greater than or equal to zero.

The base class also emits a fatal error if it cannot identify `solid` or
`region0` and no `regionName` is supplied. It prints an informational message
when `solvePressureEqn` is enabled. On the first stress correction, this law
prints whether it selected field-based or constant active tension.

### Key methods

- `impK()` delegates to the passive law. It is the diffusivity used for the
  solid model's implicit Laplacian term and affects the outer-iteration
  convergence rate rather than the converged answer.
- `correct(volSymmTensorField&)` obtains passive stress, accesses `F`, and adds
  either registry-field or ramped constant active stress divided by `det(F)`.
- `correct(surfaceSymmTensorField&)` follows the same sequence with `Ff`,
  interpolates a registry `Ta` to faces, and multiplies the active term by
  `det(Ff)` as implemented.
- `materialTangentField()` delegates the complete face list to the passive
  law. The singular inherited `materialTangent()` is not overridden and aborts
  with `notImplemented`.
- `bulkModulus()` and `shearModulus()` delegate to the passive law.

### Extension points

A related active-passive wrapper can be made by copying this class, retaining
the nonlinear runtime construction of the passive law, and replacing the
active tensor assembled in the two `correct` overloads. New implementations
should decide explicitly which registry fields are detected, whether that
choice may change after the first correction, and which cell, face, point and
quadrature interfaces they support. A consistent Newton implementation also
needs the active contribution in its material tangent rather than only the
passive-law tangent.

The source is at
[electroMechanicalLaw.C](https://github.com/solids4foam/solids4foam/blob/master/src/solids4FoamModels/materialModels/mechanicalModel/mechanicalLaws/nonLinearGeometryLaws/electroMechanicalLaw/electroMechanicalLaw.C).

---

## Tutorials

- `solids/hyperelasticity/idealisedVentricle/caseOptions/petsc`
