# Mechanical Constitutive Laws in solids4foam

This directory contains the **next-generation mechanical constitutive
modelling framework** for solids4foam.

The design separates:

- **material behaviour** (constitutive laws),
- **constitutive state** (history variables),
- **kinematics** (deformation measures),
- **topology** (where stresses are evaluated),
- **solver logic** (FV, higher-order, mixed formulations).

The goal is to provide a **general, extensible, and high-performance**
infrastructure that:

- avoids duplicated code,
- avoids copying deformation and stress data,
- supports multiple materials,
- supports both small- and finite-strain formulations,
- supports both low-order (cell-centred) and higher-order
  (multi-integration-point) discretisations,
- supports segregated and mixed (e.g. displacement–pressure) solvers.

---

## High-level overview

```text
mechanicalConstitutiveLawManager
   |
   +-- mechanicalConstitutiveLaw (one per material)
   |
   +-- mechanicalConstitutiveLawState (history, one per material)
   |
   +-- integrationPointTopology
          |
          +-- cellCentredIntegrationPointTopology
          +-- compactCellIntegrationPointTopology
```

At runtime, a client of the framework:

1. constructs a `mechanicalConstitutiveLawManager` with the mesh and the
   `mechanicalProperties` dictionary;
2. selects an `integrationPointTopology` appropriate to its storage;
3. owns the kinematic, stress, and any requested tangent storage; and
4. calls a matching manager stress-update function. The manager then:
   - constructs kinematic views,
   - passes them to the constitutive laws,
   - updates stresses and optional tangents in-place.

## Solid-model integration boundary

The framework is currently **not connected to any `solidModel`**. This is
intentional: the manager and constitutive laws must remain independently
buildable while the solid-model integration is redesigned against the current
solver interfaces.

A future solid model should use the manager only as a constitutive-update
service. It is responsible for:

- constructing and owning the manager for its mesh and `mechanicalProperties`;
- selecting the integration-point topology that matches its field storage;
- maintaining kinematics, stress, and solver state at the chosen locations;
- calling the appropriate small- or finite-strain update at each solver stage;
- applying boundary conditions, interpolation, assembly, and convergence
  control itself.

The manager is responsible for material selection, constitutive state, and
writing the requested stress response to caller-owned storage. It must not
depend on `solidModel`, create solver fields, select a solution algorithm, or
drive a solver lifecycle. This one-way dependency keeps the material framework
usable by cell-, face-, point-, and higher-order discretisations.

---

## Core classes

### `mechanicalConstitutiveLaw`

Abstract base class representing a **mechanical constitutive law**.

Responsibilities:

- compute stress from supplied kinematics,
- optionally compute approximate scalar tangents,
- provide material properties (e.g. density, bulk modulus).

Key design points:

- Constitutive laws **do not store history**.
- Constitutive laws are **stateless functions** of:
  - kinematics,
  - state,
  - output response.
- Laws operate on generic list views (`UIndirectList`), not fields.

Derived examples:

- `linearElasticMechanicalConstitutiveLaw`
- `neoHookeanElasticMechanicalConstitutiveLaw`

---

### `mechanicalConstitutiveLawState`

Stores **constitutive history variables** for a single material.

Examples:

- plastic strain,
- internal variables,
- hardening parameters,
- deformation gradients (for history-dependent laws).

Key points:

- Owned and managed by `mechanicalConstitutiveLawManager`.
- One state object per material (and per boundary patch if required).
- Sized according to the number of integration points belonging
  to the material.

---

### `mechanicalConstitutiveLawResponse`

Lightweight wrapper providing **write access** to:

- stress values,
- optional scalar tangents.

Key points:

- Does **not allocate memory**.
- Holds views (`UIndirectList`) into solver-owned storage.
- Tangent computation is enabled only when explicitly requested.

Supported tangent modes:

- `none`
- `scalar`
- `scalarDeviatoric`
- `fourthOrder`
- `fourthOrderFiniteDifference`

Fourth-order tangents are only available on integration-point topologies that
report `supportsFourthOrderTangent()`, i.e. those whose integration points are
the locations at which a Jacobian operator evaluates fluxes.

---

### `mechanicalConstitutiveLawKinematics`

Encapsulates deformation measures passed to constitutive laws.

Two main variants:

- `smallStrainMechanicalConstitutiveLawKinematics`
- `finiteStrainMechanicalConstitutiveLawKinematics`

Key points:

- Kinematics are **read-only** views.
- Laws may operate element-by-element (faster) or field-wise.
- Avoids repeated reconstruction of strain measures.

---

### `mechanicalConstitutiveLawManager`

Central orchestrator that:

- owns constitutive laws,
- owns constitutive states,
- manages material-to-integration-point mappings,
- assembles material property fields,
- dispatches stress evaluations.

Responsibilities:

- Read `mechanicalProperties`.
- Construct laws via run-time selection.
- Maintain:
  - cell-based material definitions,
  - integration-point mappings.
- Provide stress update functions for:
  - small strain / finite strain,
  - low-order (`volField`) and
  - higher-order (`CompactListList`) storage,
  - flat `UList` storage at the integration points of any topology.

The flat-list form is the primitive: the `volField` and `CompactListList`
overloads are implemented in terms of it, and it is the only form that works on
storage belonging to a mesh other than the manager's own, such as a dual mesh.
There are also tangent-only updates, `updateTangentSmallStrain` and
`updateTangentFiniteStrain`, for solid models that take their Jacobian from the
manager while their residual stress still comes from elsewhere.

The manager is **discretisation-agnostic**.

---

### `integrationPointTopology`

Abstract base class defining how **cells map to integration points**.

Responsibilities:

- describe how many integration points exist,
- map a cell to its associated integration point IDs,
- state whether integration points are shared between cells, and whether they
  can carry a fourth-order (`mat66`) material tangent.

Current implementations:

- `cellCentredIntegrationPointTopology`
- `compactCellIntegrationPointTopology`
- `faceCentredIntegrationPointTopology`
- `pointCentredIntegrationPointTopology`
- `dualFaceIntegrationPointTopology`

A topology that can be built from the mesh alone is obtained by type name from
`mechanicalConstitutiveLawManager::topologyFor()`. One that needs more, such as
`dualFaceIntegrationPointTopology` with its dual-face-to-cell map, is
constructed by the caller and handed over with `registerTopology()`, which
takes ownership.

---

## Multiple materials

Materials are always **defined per cell**, via cellZones.

- Each constitutive law corresponds to one material.
- The manager builds:
  - a list of cells per material,
  - a list of integration points per material.
- Constitutive laws are applied only to their assigned integration points.

---

## Small vs finite strain

Two independent stress update pathways exist:

- **Small strain**: engineering strain, displacement gradient based
- **Finite strain**: true (Cauchy) stress, deformation gradient based

Both pathways support:

- multiple materials,
- history-dependent laws,
- optional scalar tangents,
- low- and high-order discretisations.

---

## Tangents and mixed formulations

The framework can return optional approximate scalar or fourth-order tangents
when explicitly requested by a client. The manager does not retain a universal
"current tangent" field and does not define an implicit or stabilisation
coefficient for a solid model. The caller owns the requested tangent storage,
chooses when it is refreshed, and decides how it enters its linearisation.

This question is resolved in `DESIGN-tangents.md`. In summary: the solid model
owns its implicit stiffness fields (`impK`, `impKf`, `rImpK`) and obtains their
values from this framework via `tangentRequest::scalar`; the fidelity of the
material tangent used in the Jacobian is a separate, per-solid-model dictionary
choice. The stabilisation model supplies the discrete operator, never the
coefficient, so the two are never additively combined.

Supported current requests include:

- requested explicitly by the solver,
- computed only when needed,
- scalar tangents suitable for segregated or mixed formulations,
- fourth-order tangents where the selected update interface supports them.

---

## Higher-order discretisations

Higher-order solid models store stress and kinematics as:

```c++
CompactListList<Type>
```

The constitutive framework operates on the **packed storage only**;
topology remains solver-owned.

---

## Testing

`Test-mechanicalConstitutiveLaw` (in `applications/test/`) exercises the
manager on a case mesh and its `constant/mechanicalProperties`. It checks the
closed-form linear elastic stress and tangent, agreement between the flat-list,
`CompactListList` and `GeometricField` paths, the tangent-only update, the
dual-face topology addressing and its fourth-order tangent, and the misuse
guards.

It requires every material to be `linearElastic`, and is run by the
`layeredPipe` tutorial's `regressionTest.sh`, that being the only tutorial with
more than one material. No solid model uses this framework yet, so that is
currently its only runtime coverage.

---

## Design philosophy

This framework is designed to:

- minimise code duplication,
- maximise extensibility,
- maintain high performance,
- support future research directions.

---

## Status and future extensions

Planned extensions include:

- additional tangent types,
- a dual-face integration point topology for the vertex-centred solid models,
- a law-defined hydrostatic response for mixed displacement-pressure solvers.

See `DESIGN-tangents.md` for the design that resolves how solid models obtain
material stiffness for their residual and their Jacobian, and for the staged
plan by which solid models adopt this framework.
