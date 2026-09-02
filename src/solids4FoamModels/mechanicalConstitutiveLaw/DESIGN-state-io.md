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

## 3. The general shape of an integration-point topology

State lives at the integration points of a *topology*, and the first draft
treated those as either "cell-shaped" or "awkward". That is too coarse. Every
topology in use or planned is of one form:

> **n integration points attached to each entity of some mesh**, where the
> entity is a cell, a face or a point.

Writing it that way makes the IO question answerable uniformly, because
OpenFOAM already has a mesh-aware, parallel-aware, decomposable container for
each entity kind.

<!-- markdownlint-disable MD013 -->

| Topology | Mesh | Entity | n per entity | Natural container |
|---|---|---|---|---|
| `cellCentred` | primary | cell | 1 | `volField` |
| `faceCentred` | primary | face | 1 | `surfaceField` |
| `pointCentred` | primary | point | 1 | `pointField` |
| Gauss points in cells (`compactCell`) | primary | cell | fixed n | n `volField`s |
| Gauss points on faces | primary | face | fixed n | n `surfaceField`s |
| `dualFace` | **dual** | face | variable | none - see §7 |

<!-- markdownlint-enable MD013 -->

So the classification that matters for IO is not "is it cells" but three
independent properties:

1. **Which mesh** the entities belong to. The primary mesh is decomposed,
   reconstructed and written by the standard tools; a derived mesh such as the
   dual is not.
2. **Which entity kind**, which selects the container.
3. **Whether n is fixed**. A fixed n becomes n containers with an index suffix;
   a variable n has no entity-shaped representation at all.

**A topology is standard-IO-capable when its mesh is the primary mesh and n is
fixed.** Everything in the table except `dualFace` satisfies that, including
both Gauss-point cases. This is worth stating as an interface:

```c++
        //- The mesh entity that integration points are attached to
        enum class entityKind { cell, face, point };

        //- Entity kind for this topology
        virtual entityKind entity() const = 0;

        //- Integration points per entity, or -1 where it varies
        virtual label pointsPerEntity() const = 0;

        //- True if the entities belong to the mesh the manager was built on,
        //  rather than to a derived mesh such as a dual
        virtual bool onPrimaryMesh() const = 0;
```

The manager then chooses the container generically, and a new topology gets IO
for free by answering three questions rather than by having IO written for it.

Two qualifications carried over from review.

`setFields` sets the internal and boundary values of `volField`s, and
finite-area fields. It does not set `surfaceField`s, `pointField`s or arbitrary
`IOField`s. So the `setFields` route exists for cell-entity topologies and not
for the others - which matters for prescribing inputs (§4), not for restart.

A shared integration point - a face or point reachable from two cells - can be
reached from two *materials*. That is a semantic problem, not a storage one,
and it is what §4 has to deal with.

---

## 4. Decision

**Prescribed inputs are always specified as a `volField` on the primary mesh.
Persistent state is stored in the container its topology's entity implies.**

The two directions genuinely differ, which is why they get different answers.

**Prescribed (R2, R3).** A user specifies an initial stress, an initial strain
or a fibre direction per *region of material*, and materials are already
defined per primary cell. A cell-shaped input is therefore the natural
specification and not a loss of fidelity, and it is the one form `setFields`
can write. The manager broadcasts a cell's value to every integration point of
that cell.

This is exact where an integration point belongs to exactly one cell, i.e.
where `requiresUniqueIntegrationPointsPerMaterial()` is false - `cellCentred`,
`compactCell`, `dualFace`. For `faceCentred`, `pointCentred` and face Gauss
points an integration point is reachable from two cells that may carry
different values, and the broadcast is then ambiguous. The first draft called
the cell-shaped form "not a compromise"; that is false for shared topologies.
Those need a mapping policy stated in the declaration - reject a discontinuity
across the point, take the owner cell, or interpolate - rather than an implicit
choice.

**Persistent (R1).** The direction reverses and the mapping must be lossless,
so the container follows the entity:

<!-- markdownlint-disable MD013 -->

| Topology | Restart representation | Standard tools? |
|---|---|---|
| `cellCentred` | one `volField`, `READ_IF_PRESENT` / `AUTO_WRITE`. Identical to legacy | yes |
| `faceCentred` | one `surfaceField` | yes |
| `pointCentred` | one `pointField` | yes |
| Gauss points in cells | n `volField`s, suffixed `_0` … `_n-1` | yes |
| Gauss points on faces | n `surfaceField`s, suffixed likewise | yes |
| `dualFace` | none available | **no** - see §7 |

<!-- markdownlint-enable MD013 -->

The suffix convention is deliberate: n separate fields of a standard type
decompose, reconstruct and plot with no new tooling, where one field of n
values per entity would need all three written from scratch.

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

§4 leaves exactly one topology without a representation: `dualFace`, and any
future topology on a derived mesh or with a variable count. `faceCentred`,
`pointCentred` and both Gauss-point cases are now covered by standard
containers, so the awkward set is smaller than the first draft assumed.

For that remaining case there are three options:

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
not just topologies on a derived mesh. That includes boundary state, shared
integration points with conflicting prescribed values, and any mesh that can
change topology.

Note how little this now blocks. The cell-centred solid models and their
high-order Gauss-point variants - the ones in common use, and the priority for
migration - are all standard-IO-capable by §4, so restart works for them
through the same mechanism the legacy laws already use. Only the vertex-centred
path is held back, and it currently runs `linearElastic`, which has no history
to lose. Option 2 produces a run that
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
| 2 | `prescribed` fields read from `volField`s in `constant/` and broadcast to integration points, with a declared mapping policy for shared points | R2, R3, fibre directions, `setFields` |
| 3 | `initialStress` and `initialStrain` applied by the manager, small strain only | the specific request |
| 4 | `persistent` fields written and read in the container the topology's entity implies: `volField`, `surfaceField`, `pointField`, or n of them for a fixed Gauss-point count | R1 for the cell-centred solid models and their high-order variants |
| 5 | Guard refusing history-dependent laws wherever the mapping is not exact: derived meshes, variable counts, unresolved shared points, changing topology | honest failure instead of silent loss |
| 6 | Topology-aware persistent state with integration-point identity, plus decompose/reconstruct/redistribute support, for topologies on a derived mesh or with a variable count | R1 everywhere |
| 7 | `mapPolyMesh` handling so history survives a topology change | `crackerFvMesh` with a history-dependent law |

<!-- markdownlint-enable MD013 -->

Steps 1 to 3 are worth doing on their own: they deliver initial stress, initial
strain and fibre directions, none of which exist in the new framework today,
and none of which depend on the restart question being settled.

---

## 8a. Two categories this design does not yet cover

Migrating the solid models one at a time has now surfaced two shapes of law
that the sections above cannot express. Both were found the same way: by
asking which law each remaining tutorial actually uses.

### 8a.1 Live coupling fields

`thermoMechanicalLaw` subtracts `3 K alpha (T - T0) I` from the stress its
sub-law computed; `poroMechanicalLaw` subtracts `b (p + p0) I`. Neither `T` nor
`p` is kinematics, and neither is prescribed state in the sense of §4: they are
*solved every time step*, by another equation in the same solid model, and must
be current at the moment the law is evaluated.

So there is a third category alongside persistent (R1) and prescribed (R2, R3):

**R4. A field the solid model owns and re-supplies each update, which a law
reads but never writes.**

Temperature and pore pressure are the two present cases;
`electroMechanicalLaw` needs an activation field, and `viscousHookeanElastic`
needs the time increment it already gets.

**Proposal.** Do not change `evaluate`. A law declares a coupling field the
same way it declares any other state field, and the solid model hands the
manager the current values before each update:

```c++
    manager.setCouplingField("T", T);      // volScalarField on the primary mesh
```

The manager scatters cell values to the integration points of every law that
declared `T`, using the same broadcast and the same shared-topology mapping
policy §4 already defines for prescribed inputs. The law then reads `T` from
its state exactly like anything else, and remains a pure function of
(kinematics, state).

**Alternative rejected.** Passing a dictionary of auxiliary fields into
`evaluate`. That changes the signature every law implements in order to serve a
minority of them, and it puts field-shaped objects back into the one interface
this framework exists to keep free of them.

**Open point.** A coupling field must be *stale-proof*: if a solid model
forgets to set it, the law silently uses whatever was there last, which is the
kind of error that produces plausible wrong answers. The manager should track
whether each declared coupling field was set in the current update and fail
loudly if not, rather than defaulting.

### 8a.2 Composite laws

`thermoMechanicalLaw` and `poroMechanicalLaw` are decorators: they construct a
sub-law, delegate the stress to it, and then correct the result. So are
`electroMechanicalLaw` and, in effect, the plastic laws' elastic predictors.
Nothing in this framework yet says how a law owns another law.

**Proposal.** A composite constructs its sub-law with
`mechanicalConstitutiveLaw::New` on a sub-dictionary and owns it outright. The
manager continues to see exactly one law per material, with one state. The
composite forwards `declareState` to its sub-law, prefixing the names it
declares so two laws in a stack cannot collide, and forwards `evaluate`, then
modifies `response.stress()` in place, which is already a mutable view.

Two consequences worth stating:

* The finite-difference fourth-order tangent of the base class then measures
  the *composite's* response, not the sub-law's, which is what a solver wants.
* The scalar tangent must come from the composite, not be silently inherited:
  a thermal correction that is a pure pressure shift does not change the
  tangent, but a poro correction with a stress-dependent sub-law can.

**Open point.** Whether the state should be namespaced by prefixing names, as
proposed, or by giving each law in a stack its own state object. Prefixing is
less machinery; separate objects are harder to get wrong. This is the same
question the shadow state answered for tangent queries, and it should probably
be answered the same way.

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
   points, which is a semantic problem, not a storage one (§3). Following this
   through is what produced the general classification now in §3: the question
   is not whether a topology is cell-shaped but which mesh and entity its
   points attach to, and whether the count per entity is fixed. Gauss points in
   cells and Gauss points on faces both come out standard-IO-capable, and only
   a derived mesh or a variable count does not.
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
