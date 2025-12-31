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
solidModel
   |
   v
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

At runtime:

1. The **solid model** constructs a `mechanicalConstitutiveLawManager`.
2. The manager reads and constructs one or more
   `mechanicalConstitutiveLaw` objects.
3. The manager builds a mapping from **cells → integration points**
   using an `integrationPointTopology`.
4. During each stress update, the manager:
   - constructs kinematic views,
   - passes them to the constitutive laws,
   - updates stresses and optional tangents in-place.

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
  - higher-order (`CompactListList`) storage.

The manager is **discretisation-agnostic**.

---

### `integrationPointTopology`

Abstract base class defining how **cells map to integration points**.

Responsibilities:

- describe how many integration points exist,
- map a cell to its associated integration point IDs.

Current implementations:

- `cellCentredIntegrationPointTopology`
- `compactCellIntegrationPointTopology`

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

The framework supports **optional approximate scalar tangents**:

- requested explicitly by the solver,
- computed only when needed,
- suitable for segregated or mixed formulations.

---

## Higher-order discretisations

Higher-order solid models store stress and kinematics as:

```c++
CompactListList<Type>
```

The constitutive framework operates on the **packed storage only**;
topology remains solver-owned.

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
- face- and vertex-based integration point topologies.
