# Mechanical Constitutive Laws in solids4foam

This directory contains the **mechanical constitutive modelling framework**
of solids4foam: every solid model takes its stress, implicit stiffness and
density from it.

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
          +-- faceCentredIntegrationPointTopology
          +-- compactFaceIntegrationPointTopology
          +-- pointCentredIntegrationPointTopology
          +-- dualFaceIntegrationPointTopology
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

The manager and constitutive laws are independent of `solidModel`: a solid
model uses the manager only as a constitutive-update service. `solidModel`
constructs one manager, on first use, from `constant/mechanicalProperties`,
and the solid model is responsible for:

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

A law is selected by its runtime type, `linearElastic` say, and its class
carries the `MechanicalConstitutiveLaw` suffix. The suffix is not only a
leftover from when the legacy laws of the same names existed: foam-extend
4.1's `solidModels` library has a `Foam::linearElastic` with a
`linearElastic.H`, so the short name would clash in both the namespace and the
include path, and the suffix keeps every law named the same way.

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
`updateTangentFiniteStrain`, plus `updateScalarTangent` for the cell-centred
case, for solid models that take their Jacobian from the manager while their
residual stress still comes from elsewhere.

A tangent query never disturbs constitutive state. Each law is evaluated
against a **shadow** of its state, which aliases the parent's old-time fields -
read-only through it - and owns its own current-time fields. This works because
a law is a pure function of the kinematics and the old-time state: it reads
only `*0` fields, through the `get*` accessors, and writes only current-time
fields. History is never copied, only aliased, so the cost is one set of
current-time fields.

The same mechanism makes the `fourthOrderFiniteDifference` tangent safe: it
evaluates the law once per Voigt strain component at perturbed kinematics, and
those perturbed results must not be left where `endTimeStep()` can read them or
`storeOldTime()` can commit them.

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
- `compactFaceIntegrationPointTopology`
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

## Scalar inputs from another case directory

A law that reads a scalar field solved by some other equation - the
temperature of `thermoMechanicalLaw`, say - declares it with
`requiredScalarInputs()`, and the manager gathers it on every evaluation. It
looks, in order, for:

1. a field of that name registered by another model, such as the temperature
   a thermal solid model solves for;
2. a case directory named for that input in the law's dictionary;
3. a field of that name in this case's current, or initial, time directory.

The second is for a field computed beforehand, e.g. a temperature history
from a separate thermal run, and works for any input of any law:

```text
steel
{
    type            thermoMechanicalLaw;
    ...
    inputCaseDirectories
    {
        T           "hotCylinderTemperatureField";
    }
}
```

`<input>caseDirectory` means the same, so a `TcaseDirectory` entry works
unchanged. The manager (via
`mechanicalConstitutiveLawCaseInputs`) builds a `Time` and mesh for the
directory, relative to the case, checks that the mesh matches, and reads the
field once per time step at the current time. A step with no such field keeps
the last one read; none at all is fatal. The copy is registered under the
input's name and written. In parallel the directory is looked for inside each
`processorN`.

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

The division of responsibility is: the solid model owns its implicit stiffness
fields (`impK`, `impKf`, `rImpK`) and obtains their values from this framework
via `tangentRequest::scalar`; the fidelity of the material tangent used in the
Jacobian is a separate, per-solid-model dictionary choice. The stabilisation
model supplies the discrete operator, never the coefficient, so the two are
never additively combined.

The available requests are the `tangentRequest` enumeration:

- `none`: stress only,
- `scalar` and `scalarDeviatoric`: scalar tangents suitable for segregated or
  mixed formulations,
- `fourthOrder` and `fourthOrderFiniteDifference`: full 6x6 Voigt tangents,
  where the selected update interface and topology support them.

A tangent is computed only when it is asked for, and the caller owns the
storage it is written into.

---

## Prestress

A geometry taken from an image is already loaded: an artery imaged at
diastole carries the stress of the diastolic pressure. The `prestress` switch
in `mechanicalProperties` lets every finite-strain law see that stress without
the mesh being moved to a stress-free configuration first:

```c++
prestress       yes;

mechanical
(
    ...
);
```

The manager then reads a prestress deformation gradient `Fpre` (a
dimensionless `volTensorField`)
from the start time. `Fpre` maps the stress-free configuration onto the mesh,
and the laws are given the solver's kinematics composed with it:

- `F = Fstep & Fpre`, with the old-time `F` composed the same way,
- `Finv = inv(Fpre) & Finv_step`,
- `J = J_step*det(Fpre)`.

`Fpre` is on the right because it acts first. The solver keeps its own
mesh-relative `F` and `J` for the momentum equation; the laws return the Cauchy
stress, which belongs to the current configuration whichever reference it was
computed from, so the two agree without the solver knowing about `Fpre`.

That holds for a law whose stress is a function of the total `F`, and a law
says so through `admitsPrestress()`: the hyperelastic laws (`neoHookeanElastic`,
`MooneyRivlinElastic`, `OgdenElastic`, `StVenantKirchhoffElastic`,
`GuccioneElastic`, `HolzapfelGasserOgdenElastic`) do, and `electroMechanicalLaw`
does when its passive law does. The answer is no by default, and the manager
refuses to start with prestress on and a law that answers no, naming it. A law
that integrates an increment `F & Finv0` onto a history variable, such as
`neoHookeanElasticMisesPlastic`, would otherwise cancel `Fpre` out of the
increment and start its history from the stress-free value, so it would run as
if there were no prestress. `OgdenElastic` also takes a prescribed initial
stress `sigma0`, and with prestress on the two add: `sigma0` is pushed forward
with the composed `F`, so it is read as a stress in the stress-free
configuration, on top of the stress `Fpre` gives. A prestressed case normally
leaves it at zero.

`Fpre` is held at the cells. A face evaluation uses it interpolated to the
faces, so only cell values (and the patch values that come with a
`volTensorField`) need to be written. It is written back at every write time,
at full precision whatever the `writePrecision`, as the constitutive state is,
so a continued run reads the `Fpre` it started with rather than a rounding of
it.

What prestress does not cover, and refuses rather than evaluates without
`Fpre`:

- laws that do not admit prestress, as above,
- `coupledPressureDisplacementSolid`, whose pressure equation is an
  incremental relation on the mesh and does not see `Fpre`,
- a topology change, through which `Fpre` is not carried (no solid model
  that changes topology is a finite-strain one),
- an `Fpre` that is not `calculated` on every patch but the constraint ones:
  a patch type that sets its own value would replace the composed patch
  values `updatePrestress` writes,
- small-strain evaluations, which have no deformation gradient to compose,
- flat-list evaluations on any topology but the cell-centred one, including
  the `CompactListList` (high-order quadrature) path,
- fourth-order tangents, which the laws give with respect to the composed
  deformation gradient rather than the solver's.

The mixed formulation of `nonLinearGeometryTotalLagrangianTotalDisplacement`
(`solvePressure`) is written on the `J` the laws see, `J_step*det(Fpre)`, in
its pressure residual and in its check of the law's volumetric response, since
the pressure stands in for that response; its approximate Jacobian is
unchanged. Scalar tangents are given at the composed kinematics; they are preconditioner
estimates and need nothing more. Laws with fibre directions take those
directions in the stress-free configuration, so directions defined on the
imaged geometry should be mapped back through `inv(Fpre)` first.

`Fpre` is found with the modified updated Lagrangian formulation (Gee, Förster
and Wall, 2010): the load is applied to the imaged mesh, `Fpre` becomes
`Fstep & Fpre` and the displacement is reset to zero, and this is repeated
until the imaged mesh is in equilibrium with the load at `Fstep = I`. The
`updatePrestress` utility performs one such update between runs; each run
should be a static solve, since `Fstep` from a transient one is a snapshot of a
vibration rather than the equilibrium. The `F` a solid model writes is its own,
mesh-relative deformation gradient, so copying it into `Fpre` discards the
prestress of every earlier pass and the passes do not converge.

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

The closed-form checks need every material to be `linearElastic`; the rest
apply to any law. It is run by ten tutorials' `regressionTest.sh`, each on the
laws and mesh of its own case: `layeredPipe` (more than one material),
`perforatedPlate` and `cylinderExpansion` (history-dependent plasticity, where
a tangent query has history to disturb), `viscoTube` (viscoelasticity),
`blockPunch`, `longWall`, `rotatingCylinder`, `ratCarotid`,
`idealisedVentricle` and `LandEtAl2015/problem3` (hyperelasticity).

`Test-prestress` checks the `prestress` switch on a case mesh with
finite-strain laws: that the prestressed volField, surfaceField, flat-list and
scalar tangent updates match a manager without prestress given the composed
kinematics, that the order of composition is visible to that check, and that
the unsupported evaluations are refused.

Every solid model takes its stress and tangent from this framework, so every
tutorial's regression test is runtime coverage of it as well.

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
- a law-defined hydrostatic response for mixed displacement-pressure solvers.
