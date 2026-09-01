# Design: persistent and prescribed constitutive state

Status: proposal, revised after review. Nothing here is implemented yet.

The first draft was reviewed and came back "do not implement as written". The
substantive objections are incorporated below and summarised in §10.

Applies to the `mechanicalConstitutiveLaw` framework. Resolves OQ-2 and the
wider requirement behind it.

---

## 1. This is three requirements, not one

They are usually discussed together, and they pull in different directions.

<!-- markdownlint-disable MD013 -->

| # | Requirement | Direction | Lifetime | Loss tolerance |
|---|---|---|---|---|
| R1 | Restart: history survives stop/start | read **and** write | every write interval | **none** — must round-trip exactly |
| R2 | Initial stress / initial strain, for any law | read only | once, at construction | none, but it is a per-cell specification by nature |
| R3 | Material orientation: fibre directions and similar | read only | once, at construction | as R2 |

<!-- markdownlint-enable MD013 -->

R2 and R3 are the same mechanism: read a field the user prepared, never write
it back. R1 is different: it is a machine-to-machine round trip.

Conflating them is what makes the problem look hard. Separated, R2 and R3 are
easy and R1 is easy for most topologies and awkward for one.

---

## 2. What the legacy laws already do

This is not a new problem, and the existing answer is good:

```c++
// linearElasticMisesPlastic.C:341-353
epsilonP_
(
    IOobject
    (
        "epsilonP",
        mesh.time().timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedSymmTensor("zero", dimless, symmTensor::zero)
),
```

A `volField` on the primary mesh with `READ_IF_PRESENT` and `AUTO_WRITE`. That
one declaration satisfies R1 and R2 at once, is inspectable, decomposes and
reconstructs with the standard tools, and is directly writable by `setFields`.
`HolzapfelGasserOgdenElastic` reads its fibre fields the same way, as
`MUST_READ` volFields, which is R3.

**The new framework should not invent a different answer where the old one
works.** It only needs a different answer where the old one cannot reach.

---

## 3. Why the new framework cannot simply do the same

State in the new framework lives at the integration points of a *topology*, and
only one topology is one-to-one with cells:

<!-- markdownlint-disable MD013 -->

| Topology | Integration points per cell | Cell-shaped? |
|---|---|---|
| `cellCentred` | 1 | yes |
| `compactCell` (quadrature) | fixed n | no, but regular |
| `faceCentred` | shared between cells | no, but `surfaceField` fits |
| `pointCentred` | shared between cells | no, but `pointField` fits |
| `dualFace` | variable, many | no |

<!-- markdownlint-enable MD013 -->

Two qualifications, both from review.

`faceCentred` and `pointCentred` are **not** short of a container: OpenFOAM has
mesh-aware surface and point fields with proper parallel handling. What they are
short of is *semantics* - an integration point shared between two cells can be
reached from two materials, so a single value per face is not obviously the
right model. That is a different problem from having nowhere to put the number.

`setFields` sets the internal and boundary values of `volField`s, and finite-area
fields. It does not set `surfaceField`s, `pointField`s or arbitrary `IOField`s.
So "a per-cell `volField` is what `setFields` can write" is right; "and nothing
else" was too strong.

---

## 4. Decision

**The user-facing form is always a `volField` on the primary mesh. The manager
maps between that and integration points.**

Materials are already defined per primary cell, `setFields` operates on
volFields, and a per-cell initial condition is what a user actually wants to
specify. The manager broadcasts a cell's value to every integration point of
that cell.

**This is exact only where an integration point belongs to exactly one cell**,
i.e. where `requiresUniqueIntegrationPointsPerMaterial()` is false -
`cellCentred`, `compactCell`, `dualFace`. For `faceCentred` and `pointCentred`
an integration point is reachable from two cells that may carry different
prescribed values, and the broadcast is then ambiguous. The first draft called
the cell-shaped form "not a compromise"; that is false for shared topologies.
Those need a stated mapping policy - reject a discontinuity across the point,
take the owner cell, or interpolate - and the policy has to be part of the
declaration rather than an implicit choice.

For R1 the direction reverses and the mapping has to be lossless, so it depends
on the topology:

<!-- markdownlint-disable MD013 -->

| Topology | Restart representation | Standard tools? |
|---|---|---|
| `cellCentred` | one `volField`, `READ_IF_PRESENT` / `AUTO_WRITE`. Identical to legacy | yes, fully |
| `compactCell`, fixed n per cell | n `volField`s, suffixed `_0` … `_n-1` | yes, fully |
| `dualFace`, `faceCentred`, `pointCentred` | `IOField` under `<time>/mechanicalState/<topology>/<law>/<name>` | **no** — see §7 |

<!-- markdownlint-enable MD013 -->

---

## 5. Laws declare their state; the manager owns the IO

Today each legacy law constructs its own `IOobject`s. That is exactly the
coupling this framework exists to remove: a law is supposed to be a pure
function with no knowledge of meshes, files or registries.

So laws *declare* what they need and the manager does the rest:

```c++
        //- Declare the state fields this law uses.
        //  Called once, before any evaluation. The manager allocates the
        //  fields, applies the defaults, reads any prescribed values, and
        //  arranges for the persistent ones to be written
        virtual void declareState(mechanicalConstitutiveLawStateSpec& spec) const
        {}
```

with

```c++
enum class stateRole
{
    //- History. Written every write interval, read back on restart (R1)
    persistent,

    //- Prescribed once from a field the user prepared, never written back
    //  (R2, R3). Absent means the default is used
    prescribed
};
```

and a law declaring, for example:

```c++
void linearElasticMisesPlasticMechanicalConstitutiveLaw::declareState
(
    mechanicalConstitutiveLawStateSpec& spec
) const
{
    spec.add<symmTensor>("epsilonP", stateRole::persistent, symmTensor::zero);
    spec.add<scalar>("epsilonPEq", stateRole::persistent, 0.0);
    spec.add<scalar>("sigmaY", stateRole::persistent, yieldStress(0.0));
}
```

The manager then, per law and topology: allocates the `Field<Type>` at the
right size, looks for a `volField` of that name on the primary mesh, broadcasts
it to the law's integration points if found, and registers the persistent ones
for writing.

This keeps every law free of IO, and it makes the set of state variables
introspectable, which is worth having on its own.

---

## 5a. Two things the first draft omitted entirely

**Boundary integration points.** A `boundaryAware()` topology carries a
*separate* `mechanicalConstitutiveLawState` per law per patch
(`mechanicalConstitutiveLawManager.C:275-306`), and commits it alongside the
internal one (`:381-390`). `cellCentred` is boundary-aware, so this is the
common case, not an exotic one. A single `volField` covers it only if the
mapping explicitly carries the internal field to the internal state and each
patch's boundary values to that patch's state, in both directions. The
declaration therefore has to name the boundary state as part of the same
variable, not leave it implicit.

**Mesh topology change.** `crackerFvMesh` calls `changeMesh()`
(`crackerFvMesh.C:561-579`). Constitutive state is a raw `Field` whose
`setSize()` zero-fills any growth (`mechanicalConstitutiveLawState.C:215-253`),
and the manager caches topologies and their state for its whole lifetime
(`mechanicalConstitutiveLawManager.H:142-191`). So a topology change today
would silently zero the history of every new cell and misalign the rest.

This is a blocker rather than a detail, and it affects *continuing* runs, not
just restarts. The framework needs a `mapPolyMesh` path that rebuilds the
topology and its addressing and maps **both current and old-time** state,
including onto newly created crack faces and points. Until that exists, a
history-dependent law on a changing mesh must be refused.

---

## 6. Initial stress and initial strain are manager-level, not per-law

The requirement is that they work for **any** law, including `linearElastic`.
Adding a `sigma0` term to every law would be repetitive, easy to get wrong, and
would have to be repeated for every law written in future.

Instead the manager applies them around the law:

```text
initial strain:   the law is evaluated at (gradD - gradD_initial)
initial stress:   sigma += sigma_initial, after evaluation
```

Both are `prescribed` fields with the reserved names `initialStrain` and
`initialStress`, allocated only if the corresponding `volField` is present.
When absent there is no field, no branch in the inner loop, and no cost.

Note the name: `gradD0` already means the *previous-time* displacement gradient
(`smallStrainMechanicalConstitutiveLawKinematics.H:86-92`) and must not be
reused for a prescribed initial value.

**Limitations, stated rather than hidden.** Both were understated in the first
draft.

* The additive treatment of initial strain is a **small-strain, total-form**
  construction. It is not automatically valid for a law written in incremental
  or rate form, and finite strain needs a multiplicative split
  `F = F_e * F_initial`. The reference convention has to be declared per
  kinematic formulation rather than assumed, and anything outside the case it
  is defined for must be refused rather than silently applied.
* `initialStress` needs a stated stress measure and frame. Adding a prescribed
  field to whatever a law returns is only meaningful once "Cauchy, current
  configuration" - or whatever the convention is - is written down.

---

## 7. The awkward case, and what to do about it

For `dualFace` and the other shared-point topologies there is no cell-shaped
lossless representation, because the number of integration points per cell
varies. Three options:

1. **A topology-aware persistent field per topology per law.** Not a bare
   `IOField`: that is an ordered `Field` plus IO and nothing more
   (`IOField.H:53-57`), with no mesh addressing and no identity for an
   integration point, so it is only safe on an unchanged mesh with an identical
   processor layout and ordering. A real solution carries global or topological
   identity for each integration point and supports decompose, reconstruct and
   redistribute - which means a purpose-built `regIOobject`, not a file in a
   subdirectory.
2. **Reduce to one value per cell.** Standard tools work, everything is
   inspectable, and restart is *approximate*: the loss is exactly the sub-cell
   variation of history, i.e. a restarted run behaves as though history had
   been homogenised within each cell.
3. **Refuse.** Error if a history-dependent law is used on a non-cell-shaped
   topology, until 1 is implemented.

**Recommendation: 3 now, 1 next, never 2.** With one widening from review: the
guard in 3 must cover *every* combination that lacks an exact defined mapping,
not just the multi-point topologies. That includes boundary state, shared
integration points with conflicting cell values, and any mesh that can change
topology. Option 2 produces a run that
restarts to a subtly different state with no indication that anything was lost,
which is the worst property a restart can have. Option 3 costs one guard and is
honest. Option 1 is the real answer and can be built when someone needs a
parallel restart of a vertex-centred elastoplastic case.

Note that the vertex-centred solvers currently in the tree only run
`linearElastic`, which has no history at all, so option 3 blocks nothing that
works today.

---

## 8. Staging

<!-- markdownlint-disable MD013 -->

| Step | Content | Unblocks |
|---|---|---|
| 1 | `declareState` and the spec object; manager allocates and applies defaults. No IO yet | the rest |
| 2 | `prescribed` fields read from `volField`s and broadcast to integration points | R2, R3, fibre directions, `setFields` |
| 3 | `initialStress` and `initialStrain` applied by the manager, small strain only | the specific request |
| 4 | `persistent` fields written and read as `volField`s for cell-shaped topologies | R1 for cell-centred solid models |
| 5 | Guard refusing history-dependent laws on non-cell-shaped topologies | honest failure instead of silent loss |
| 6 | Topology-aware persistent state with integration-point identity, plus decompose/reconstruct/redistribute support | R1 everywhere |
| 7 | `mapPolyMesh` handling so history survives a topology change | `crackerFvMesh` with a history-dependent law |

<!-- markdownlint-enable MD013 -->

Steps 1 to 3 are worth doing on their own: they deliver initial stress, initial
strain and fibre directions, none of which exist in the new framework today,
and none of which depend on the restart question being settled.

---

## 9. Open questions

* **Q1. Answered, and it could not be deferred.** A `prescribed` field that is
  never written would simply vanish on restart, because construction searches
  the restart time directory and not `0/`. So prescribed inputs are read from
  `constant/`, which is where time-independent material data belongs and where
  the fibre-angle inputs of the legacy anisotropic law would also sit. That
  keeps them present on restart without writing them every interval.
* **Q2.** Should `persistent` fields be written every write interval, or only
  on demand? Every interval matches the legacy `AUTO_WRITE` behaviour and makes
  restart automatic, at the cost of larger output for laws with several state
  variables. Whichever is chosen, only converged end-of-time-step state may be
  written: an intermediate nonlinear-iteration value must never be persisted as
  restart history. On load, current and old-time state must both be initialised
  coherently before the first evaluation, since a law reads old-time state and
  the manager only rolls current into old when the time index advances
  (`mechanicalConstitutiveLawManager.C:354-395`).
* **Q4.** What is the reference configuration and stress measure for
  `initialStress`, and which kinematic formulations may `initialStrain` be
  offered for? §6 needs these written down before it is implemented.

* **Q3.** Where a law's state has a natural physical name that collides across
  laws in a multi-material case (`epsilonP` for two different plastic laws),
  does the on-disk name need a per-material qualifier? The legacy code sidesteps
  this with sub-meshes, which this framework removes.

---

## 10. What the review changed

The first draft was reviewed by an independent agent, which recommended against
implementing it as written. The substantive corrections, all incorporated above:

1. `setFields` writes volFields *and* finite-area fields; "and nothing else"
   was too strong (§3).
2. `faceCentred` and `pointCentred` are not short of a container - mesh-aware
   surface and point fields exist. Their difficulty is shared integration
   points, which is a semantic problem, not a storage one (§3).
3. The per-cell broadcast is exact only where an integration point belongs to
   one cell. Calling it "not a compromise" was wrong for shared topologies,
   which need a declared mapping policy (§4).
4. **Boundary state was omitted entirely.** `boundaryAware()` topologies keep a
   separate state per law per patch, and `cellCentred` is boundary-aware, so
   this is the common case (§5a).
5. **Mesh topology change was omitted entirely.** `crackerFvMesh` changes the
   mesh; state is a raw `Field` that zero-fills on growth, so history would be
   silently corrupted on a continuing run, not merely on restart (§5a).
6. A bare `IOField` is not a restart format: no addressing, no integration-point
   identity, unsafe across decomposition. The lossless option needs a
   purpose-built `regIOobject` (§7).
7. Prescribed fields would vanish on restart, since construction reads the
   restart time and not `0/`. They move to `constant/`, which answers Q1 (§9).
8. §6 reused the name `gradD0`, which already means the previous-time gradient,
   and over-claimed generality: the additive treatment is small-strain
   total-form only, and initial stress needs a stated measure and frame (§6).

The review endorsed keeping `declareState` rather than letting laws own their
own IO, on the condition that the descriptor carries dimensions, persistence,
scope, mapping policy and validation. It also suggested keeping immutable
orientation inputs conceptually distinct from evolving state even where they
share a loading service, which is worth doing.
