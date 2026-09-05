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

Migrating the solid models one at a time surfaced two shapes of law that the
sections above cannot express. Both were found the same way: by asking which
law each remaining tutorial actually uses.

The first draft of this section proposed carrying live inputs in the state
object. **That was wrong**, and review killed it on a concrete point rather
than a stylistic one. It is rewritten below, with the refuted version kept in
§8a.4 because the reason it fails is the most useful thing here.

### 8a.1 Live coupling inputs

`thermoMechanicalLaw` subtracts `3 K alpha (T - T0) I` from the stress its
sub-law computed; `poroMechanicalLaw` subtracts `b (p + p0) I`. Neither `T` nor
`p` is kinematics, and neither is prescribed state in the sense of §4: they are
*solved every time step*, by another equation, and must be current at the
moment the law is evaluated.

So there is a third category alongside persistent (R1) and prescribed (R2, R3):

**R4. A value the solid model owns and re-supplies each evaluation, which a law
reads and never writes.**

`electroMechanicalLaw` needs an activation field, and `viscousHookeanElastic`
needs the time increment, which it already receives.

**Decision: an inputs object, passed to `evaluate`. Implemented.**

```c++
    virtual void evaluate
    (
        const smallStrainMechanicalConstitutiveLawKinematics& kin,
        const mechanicalConstitutiveLawInputs& inputs,
        mechanicalConstitutiveLawState& state,
        mechanicalConstitutiveLawResponse& response
    ) const;
```

`mechanicalConstitutiveLawInputs` is immutable and per-call. It holds
integration-point views, not `GeometricField`s: the manager scatters a
`volScalarField` once and hands the law a `UIndirectList`-shaped view, so the
law still never sees a mesh, a registry or a field type. It also carries the
global scalars, `dt` among them - see §8a.5, which moves `dt` out of the
kinematics object where it does not belong.

This costs a one-time signature change on every law. It buys the property that
matters: **every evaluation, including every shadow and every
finite-difference perturbation, receives the live values by construction.**
There is no setter, no ordering protocol, no epoch counter, and no way to be
stale.

**What is implemented, and what is not.** The object exists, carries `dt`, and
is threaded through every `evaluate`, every manager call path and every
finite-difference perturbation. `dt` has moved out of both kinematics classes,
which are now purely geometric. Laws ask for a live field by name, with
`findScalar` for an optional input and `getScalar` for a required one, the
latter failing rather than reading zero.

What is *not* implemented is the manager side: nothing yet scatters a
`volField` into a per-law input, because no law reads one yet and the mapping
policy for shared and boundary points, §8a.3, is still open. The first
consumer, `thermoMechanicalLaw` or `poroMechanicalLaw`, settles that. Until
then the `dt` path through the manager is compiled but unexercised, since no
law reads it either - `viscousHookeanElastic` will be the first.

### 8a.2 Composite laws

`thermoMechanicalLaw` and `poroMechanicalLaw` are decorators: they construct a
sub-law, delegate, and correct the result. So is `electroMechanicalLaw`.
Nothing in this framework says how a law owns another law.

A composite constructs its sub-law with `mechanicalConstitutiveLaw::New` on a
sub-dictionary and owns it. The manager continues to see one law per material.
Modifying `response.stress()` after delegating works as intended: it is a
non-const view over the caller's storage, so there is no copy to lose.

Three things the first draft got wrong or missed.

**Prefixed state names do not work.** Declaring `sub.epsilonP` on the parent's
behalf does not make the sub-law's own `state.getSymmTensorField("epsilonP")`
lookups qualified. A declaration-only prefix allocates fields nobody reads.
What is needed is either a scoped facade that prefixes every *lookup*, or a
child state object per layer.

**Child states need recursive lifecycle support.** They are the safer model,
and they are now implemented. `mechanicalConstitutiveLawState::child(name)`
returns a state owned by this one, created on first use at the parent's size.
A composite hands its sub-law that child, so the sub-law looks its history up
by its own unqualified name and nothing collides. Prefixing was the alternative
and does not work, for the reason above.

Two lifecycle operations recurse. `setSize` resizes the children, and
`storeOldTime` rolls them over - without which a sub-law would read this step's
values as though they were last step's. Restart IO does not recurse yet,
because it does not exist at all; that is §5.

Shadowing recurses too, and this is the part worth stating plainly: `child()`
on a shadow returns *a shadow of the parent's corresponding child*, created
lazily. Anything else and a tangent query evaluated into a shadow would reach
shadowed history at the top level and the sub-law's real history underneath -
writing to it, and committing perturbed values. `Test-mechanicalConstitutiveLaw`
checks exactly that: a shadow's child is not the parent's child, is itself a
shadow, reads the real child's history, and writing through it leaves the real
child untouched.

One deliberate asymmetry: the rollover walks the old-time table, so a field
that has a current entry and no old-time one is not history and is not rolled
over. That is how a law says a field is scratch rather than state.

**The sub-law's material properties are not reachable.** `kappa()` does exist
on the interface, so the first draft was wrong to say no properties are
exposed, but it returns a single `dimensionedScalar`. The legacy thermal
correction uses a sub-law bulk modulus that may be field-valued, and
interpolates it to faces. A general decorator therefore cannot express a
state-, point- or material-dependent thermal stiffness through today's
interface. Either material properties become a per-integration-point response,
or the thermal correction is admitted to be law-specific rather than a general
decorator.

### 8a.2c The first composite: thermoMechanicalLaw

`thermoMechanicalLaw` is migrated and reproduces the legacy law bit for bit on
`hotSphere`. It is worth recording what the composite actually needed, because
only one of the three costs listed above turned out to bite.

**Child states, as designed.** The sub-law is handed `state.child("mechanicalLaw")`
and goes on looking its history up by its own names. `applyStateSpec` recurses
through `childStateNames()` and `childLaw()`, so a sub-law's declared defaults
and prescribed fields are applied to its child exactly as they would be if it
were the top-level law. Without that recursion a sub-law's `sigma0` would be
silently absent.

**The bulk modulus did not bite, but only here.** The thermal term uses the
sub-law's bulk modulus as a single value. Every tutorial that uses this law
has `linearElastic` underneath, whose bulk modulus is uniform, so `kappa()` is
exact. A sub-law whose bulk modulus varies point to point still cannot be
expressed through this interface; that limitation is unchanged and is noted at
the point of use.

**The temperature is an input, not state, and that distinction earned its
keep.** `T` is declared through `requiredScalarInputs()` and gathered by the
manager on every evaluation, per law, as a view of that law's own integration
points. It is never stored and never rolled over.

Two things about the gather were not obvious.

The registry comes first and the file second - the opposite order from a
prescribed field. A coupling input is solved for, so a live field must never be
shadowed by a file. But the file is still needed, because of *when* the first
evaluation happens: a solid model builds its implicit stiffness while it is
being constructed, and a derived model that solves for the coupling field has
not reached its own members yet, since the base class constructor runs first.
At that moment `T` exists only as the initial condition on disk, which is the
right value to evaluate against anyway.

Boundary faces gather from the field's *patch* values, into storage of their
own. Sharing storage with the internal gather would have left whichever ran
last in place for both.

**Paths that do not gather say so.** Five evaluation paths do not gather
coupling inputs yet. Rather than hand a law an empty inputs object and let it
read a zero it cannot tell from a real value, they build their inputs through
`inputsWithoutCoupling`, which is fatal if any law asks for one. That guard
found a real gap during this work: the volField path gathered for its internal
points but not for its boundary faces.

### 8a.2d Wiring thermalLinGeomSolid, and what the two arms agree to

`thermalLinGeomSolid` now offers the framework, which is what lets the thermal
tutorials carry the comparison themselves rather than it being done by hand on
a substituted solid model.

The implicit stiffness is the one piece worth noting. `solidModel` grew a
shared `frameworkImpK`, because the field is not free to differ: it is
registered under `impK` and looked up by that name by the contact and cohesive
zone models, so its name, dimensions and boundary types have to match the
legacy one exactly. Both solid models that offer the framework now build it
there rather than each keeping a copy.

The face value differs in how it is reached. The framework has no separate face
tangent, so `impKf` is the interpolate of the cell field, where the legacy
model forms a face value directly. For a law whose stiffness does not vary
within a material these agree; what they do not do is steer the iteration
identically.

That shows up in what the two arms agree to, and it is worth being precise
because two of the three thermal tutorials behave differently:

* `slabCooling` runs a single time step and the two arms agree exactly.
* `hotSphere` runs five, and they agree to about 1.8e-7 of the displacement -
  the solution tolerance, not the model. Tightening `solutionTolerance` from
  1e-6 to 1e-10 takes the difference from 8.6e-8 to 1e-8 of the displacement,
  which is what a convergence-path difference looks like and not what a
  constitutive difference looks like.

Both tutorials now carry the comparison. `hotSphere` is the one that matters
for state, being the only one that runs enough steps to roll it over, and it
writes 14 significant figures while doing so: at the default six the two arms
differ by 3e-6 simply because that is the last digit written, and the check
would be measuring the file format.

### 8a.2e The second composite: poroMechanicalLaw

`poroMechanicalLaw` is migrated and reproduces the legacy law to round-off on
`rodAndSeabed`. `poroLinGeomSolid` is wired the same way as the thermal one.

The note in §8a.2 said this law "is not sub-law stress minus pressure", and
that is right, but it is worth being exact about when it matters. Legacy seeds
an effective stress once, `sigmaEff = sigma + b(p + p0)I`, and gives the
sub-law that field to work in. Whether that differs from subtracting the
pressure from whatever the sub-law writes depends on one thing: whether the
sub-law reads its incoming stress. No framework law does - every one of them
computes its stress from the kinematics and its own history - so for the laws
that exist today the two are numerically identical and the seeding is
unobservable.

It stops being identical for a sub-law whose strength depends on the stress
state, which is the Mohr-Coulomb case the legacy comment calls out. So the
seeding is reproduced rather than simplified away: `sigmaEff` is declared as
persistent state, seeded per point on first evaluation, and the sub-law is
given it to work in. "Identical until someone migrates Mohr-Coulomb, then
silently wrong" is not a trade worth taking.

The sub-law is handed the effective stress through the caller's storage rather
than through a list of its own - copy in, evaluate, copy out - because the
response wraps an indirect list and building an identity-indexed view over the
state field would cost an allocation per evaluation to save two passes.

The shipped `rodAndSeabed` now runs both ways and agrees bit for bit, which
took one more thing: `anisotropicBiotElastic`, its effective-stress law,
selects its two-dimensional form from `mesh.solutionD()`, and a framework law
is constructed from a dictionary with no mesh. The manager now injects
`solutionD` alongside `planeStress`, for the same reason and by the same route.

That case is also the one that shows why the effective stress is carried rather
than recomputed. `anisotropicBiotElastic` does not write the zz, yz and xz
components of the stress in the branch this case takes, so they are whatever
was in the storage it was handed. Give it the caller's total stress and those
three components are wrong; give it the effective stress and they are right.
The composite was written that way on the strength of an argument about
Mohr-Coulomb, and it turned out to be load-bearing for the very first sub-law
tried.

One thing this does not do.

**`sigmaEff` is no longer written for post-processing.** The legacy law
registers it `AUTO_WRITE`, so a case could look at the effective stress
directly. In the framework it is per-integration-point state and is not a
registered field, so nothing is written. That is a real capability difference
rather than an oversight, and it is what §5's state IO would restore.

### 8a.2b Corrections around a pure law: what survives, and what does not

This section proposed replacing decorator laws entirely: a material would
declare an input correction (an eigenstrain subtracted from the strain) and an
output correction (an additive stress), the manager would apply them, and
§8a.2's child-state problem would evaporate.

**Review rejected it, and was right.** The taxonomy is worth keeping; the claim
that it *replaces* composition is false. What follows is what stands, what
fell, and why - the failures being the more useful half.

#### What stands

* **The thermal case really is an eigenstrain.** For three-dimensional
  isotropic linear elasticity, `C : eps0 = 2 mu dev(eps0) + K tr(eps0) I`, and
  `eps0 = alpha (T - T0) I` is spherical, so this is `3 K alpha (T - T0) I` -
  identically the legacy stress correction. The framework and legacy laws build
  `K = lambda + (2/3) mu` from the same plane-stress and plane-strain formulae,
  so the equivalence holds under plane stress too.
* **Applying it on the input is more correct than subtracting a stress.** The
  legacy wrapper accepts *any* linear-geometry sub-law and subtracts
  `3 K alpha dT I` regardless. For a nonlinear or history-dependent sub-law
  that is not the same thing, and the eigenstrain form is the defensible one.
* **The two directions are genuinely different.** An eigenstrain changes what
  the material sees; an additive stress does not. That distinction is real and
  the nested dictionary cannot express it.

#### What fell

**1. "Poro needs nothing more" is false.** `checkSigmaEffReady` is not lazy
allocation with a throwaway value: it seeds the effective stress as
`sigma_total_incoming + b (p + p0) I`, and that seed is the *initial effective
stress*.

`stripFooting` nests `linearElasticMohrCoulombPlastic` inside
`poroMechanicalLaw`. That law is **pressure-sensitive**, yields on principal
effective stresses, starts its trial stress from an accumulated increment plus
`sigma0`, and warns explicitly that starting from zero stress is problematic.
Replacing the wrapper with "law output minus pore pressure" starts the law from
its own default state instead of the incoming total stress converted to
effective stress, which changes the first return mapping and the restart path.

The requirement does not vanish; it **moves** to the prescribed initial stress
of §4 and §6, which is designed and not implemented. Poro cannot be migrated
before that is.

**2. "Neither correction perturbs the Jacobian" is false in general.** The
correct statement is narrower. An additive stress leaves the constitutive
derivative alone only if it is deformation-independent. An eigenstrain adds no
term, since `d eps0 / d eps = 0`, but the law's tangent is then evaluated at
the *mechanical* strain, so it can change branch. A spherical thermal
eigenstrain does not disturb a J2 yield test, which is deviatoric - that
special case survives - but the Mohr-Coulomb law in this repository is
pressure-sensitive, and there a spherical eigenstrain can turn yielding on and
off.

**3. `electroMechanicalLaw` does not fit the additive-stress correction.** Its
active stress is `F (Ta f0 (x) f0) F^T / J`: a push-forward of a reference
structural tensor, so it depends on the deformation gradient, carries its own
tangent, and needs reference fibre data. `idealisedVentricle` nests it over
`GuccioneElastic`. That is real nesting, in a real case, that one
strain-independent eigenstrain plus one strain-independent additive stress
cannot represent. Active tension is a *material contribution*, not a coupling
correction.

**4. The legacy shim cannot be generic.** The three wrappers do not even agree
on the key holding the inner law: thermo uses `mechanicalLaw`, poro uses
`effectiveStressMechanicalLaw`, electro uses `passiveMechanicalLaw`. Any
translation is wrapper-specific, and for thermo over a nonlinear sub-law it
would silently change the model. Conversion must be opt-in and only where
equivalence is proven, not a blanket promise that old cases keep running.

#### Where this leaves the design

Both mechanisms are needed, and the boundary between them is the point:

* **Manager-applied corrections** for effects that are affine, stateless, and
  small-strain - the thermal eigenstrain, and swelling. These need no sub-law
  and no child state.
* **Composite laws** for the rest, because two real cases need them: poro,
  which needs an initial effective stress and a pressure-sensitive sub-law, and
  electro, whose active stress is finite-strain and kinematics-dependent. So
  §8a.2's child-state question is **not** avoided and must be answered.

A correction is therefore not "something the manager does"; it is a declared
component with a contract: which strain measures it supports, which inputs it
requires and their dimensions, whether it carries state, and whether it
contributes to the tangent. Anything that needs more than that contract is a
law, not a correction.

The dictionary split of `eigenstrain` and `additiveStress` survives for the
first category. The `inputs` object will also need tensor-valued views before
electro can be expressed at all; today it offers scalars and vectors only.

### 8a.3 What is still unresolved

* **Coupling on a shared topology is not an input policy, it is physics.** §4
  lets a prescribed value on a `faceCentred` or `pointCentred` point be taken
  from the owner cell when two materials disagree. For a temperature or a pore
  pressure that choice changes the constitutive response on one side of a
  material interface. It needs its own answer, not a reused one.
* **Boundary states are not covered.** The manager keeps per-patch states, and
  a coupling input must reach them with the field's *patch* values. Broadcasting
  internal cell values to a boundary is wrong for exactly the thermal and
  pressure cases this is for. §5a already flags boundary state as a blocker.
* **Not every source is a field on the primary mesh.** The legacy poro law
  looks its pressure up in a configured registry or region, and the legacy
  thermal law can read and map a field from another mesh entirely. Assuming
  "a `volField` on the solid mesh" silently removes those capabilities unless
  the solid model owns the inter-region mapping and hands over the result.
* **`poroMechanicalLaw` has hidden state.** It is not "sub-law stress minus
  pressure": it owns a separate effective-stress field, seeds it once from the
  incoming *total* stress as `sigmaEff = sigma + b (p + p0) I`, and passes the
  effective stress to the sub-law, precisely so that a stress-dependent
  strength sees effective and not total stress. A coupling input alone does not
  reproduce that initialisation contract, and migrating it naively would change
  the behaviour of any stress-dependent sub-law.
* **Restart and iteration ordering.** Knowing an input was supplied this update
  does not establish that it was supplied *after* the equation that produces it
  was solved. The inputs object makes staleness impossible within an
  evaluation; it does not by itself order the coupling loop.

### 8a.8 Prescribed state: what is implemented

A law declares the state it needs by overriding `declareState`, which is handed
a `mechanicalConstitutiveLawStateSpec`. Each declaration gives a name, a type,
a default, and a role. Two roles exist: `persistent` for history the law writes
and reads back, and `prescribed` for a value supplied from outside that the law
only reads. Only `prescribed` is implemented; `persistent` is declared and
defaulted but its restart IO is still the open item of §5.

The manager applies a declaration in `applyStateSpec`. Every declared field is
allocated and set to its default at both current and old time. A `prescribed`
field is then looked for as a `volField` of the same name, first in the
registry and then as a file in the current time directory, and broadcast to the
law's integration points. Absence is not an error: the default stands, which is
what lets a law declare a prescribed field without changing a single existing
case.

Two consequences are worth stating because neither is obvious.

**The default is the dictionary value.** `linearElastic` reads a uniform
`sigma0` from its own dictionary and declares *that* as the default. The
supplied-field route and the dictionary route are therefore the same mechanism
seen at two points: the dictionary sets the value everywhere, and a field, if
there is one, overrides it point by point. This is what makes the legacy
dictionary work unchanged while a spatially varying initial stress becomes
available to any law that declares one.

**Boundary states need the same treatment, and by a different route.** The
manager keeps a separate state per law per patch, evaluated exactly as the
internal points are. Those states were initially left out of the state spec,
and the result was not an error but silence: `mechanicalConstitutiveLawState`
creates a field on first access, so a law asking for `sigma0` on a patch got a
zero field and the initial stress was quietly dropped on every boundary. The
symptom was small - a fraction of a percent in the displacement - and it was
only visible against a legacy run over a real load path. A boundary state is
now filled by `applyStateSpecPatch`, which takes the field's *patch* values
rather than broadcasting cell values, per §8a.3.

That silent creation on access is worth revisiting. It is convenient for a law
initialising its own history, and it turns a mistyped state name into a field
of zeros rather than an error. The state spec is the beginning of the answer,
because a declared field no longer relies on it.

Four further points, each of which was a way to get a different answer from
legacy:

* **The dictionary wins over the field, because that is what legacy does.**
  `mechanicalLaw::makeSigma0` reads a `sigma0` field `READ_IF_PRESENT`, and the
  `linearElastic` constructor then assigns any dictionary `sigma0` over the
  whole of it. A case carrying both therefore runs on the dictionary value.
  The framework says this with a third role, `fixed`: the law declares `sigma0`
  as `fixed` when its dictionary supplied one and as `prescribed` when it did
  not, so the field is consulted only in the second case.
* **A face integration point takes the interpolated value, not a cell's.**
  Legacy forms its face initial stress as `linearInterpolate(sigma0)`. A
  broadcast would instead hand a face the value of whichever adjacent cell was
  visited last, which differs wherever the field is not uniform and depends on
  cell ordering and on the decomposition. A topology answers
  `integrationPointsAreFaces()` when its integration point index is the face
  index, and is then filled by interpolation. A topology holding several points
  per face - the compact face one - cannot answer that, and still broadcasts;
  §4's mapping question stands for it. Note that no solid model currently
  evaluates through the face-centred topology - `linGeomTotalDispSolid` uses
  the cell-centred and quadrature paths - so this rule is written to match
  legacy but is not yet exercised end to end by a tutorial.
* **A prescribed field is read at old time.** A shadow state aliases its
  parent's old-time fields and owns current-time fields that begin empty, on
  the premise that a law reads history and writes current values. A prescribed
  input breaks that premise: it is read, not written. Holding it at both times
  and reading the old one puts it back inside the premise, so a tangent query
  evaluated into a shadow sees the prescribed value rather than a silently
  zero field.
* **A scratch state gets defaults but no prescribed read.** The faces of an
  `empty` patch are addressed through the `polyPatch`, because the `fvPatch`
  has none, so there is no patch field to read a prescribed value from. They
  take their declared default. Those faces take no part in the discretisation,
  so nothing downstream sees the difference.

**Where a prescribed field is looked for, and in what order.** The file comes
first and the object registry only after it. The other way round reads more
naturally and is wrong: the legacy mechanical law registers a `sigma0` of its
own, read from the current time directory alone, so on a restart there is a
registered field of zeros standing in front of the file the case actually
wrote. Reading the file first means what the case says on disk is what the law
gets, and leaves the registry to serve fields that only ever exist in memory
because another model computes them.

File-first is right *for prescribed state* and would be wrong for a live
coupling input. A temperature or a pore pressure is solved for and changes
every iteration, so a stale file of the same name must never stand in front of
the computed field. That is the distinction between this mechanism and §8a.1:
prescribed state is supplied and then fixed, a coupling input is produced by
another equation. They travel by different routes for that reason, and a
future coupling input must not be routed through `prescribed` merely because
the plumbing looks similar.

The file is looked for in the current time directory and then in `0`. A
prescribed field describes the case rather than the state at a particular
time: it is written once, into `0`, and a restart from a later time would
otherwise not find it. The legacy Guccione law already does this for its fibre
directions, so the framework follows it rather than inventing a rule. This is
a deliberate divergence from legacy's `sigma0`, which looks only beside the
fields it restarts from and so loses the initial stress unless an earlier run
had already written it forward - a hole rather than a behaviour worth
reproducing. `cantilever2d` checks it by restarting and requiring the result
to match the continuous run exactly.

**Fibre fields: the mechanism is here, the verification is not.** A fibre
direction is the same shape of thing as an initial stress - a uniform value in
the law's dictionary when `uniformFibreField` is set, and a field of the same
name otherwise - so `prescribed` covers it as declared, with vector and tensor
types already supported and the `0` fallback matching what the legacy Guccione
law does for exactly these fields.

What is missing is a way to check a migrated fibre law against its legacy
counterpart over a real load path, which is the only check that has caught
anything. Every tutorial that uses one - `heartTissueBeam`, `ratCarotid`,
`idealisedVentricle` - runs through `coupledPressureDisplacementSolid`, which
is not one of the solid models wired to the framework and is foam-extend only,
and none of them carries a regression test. Migrating a fibre law before that
is settled would mean writing several hundred lines of anisotropic
constitutive code with no way to tell whether it reproduces the law it
replaces.

Two ways out, and the choice is open: port `coupledPressureDisplacementSolid`
to the framework, which is the direct route and inherits its foam-extend-only
restriction; or wait for a case that exercises a fibre law through a solid
model already on the framework. One further thing a fibre law will need that
`sigma0` did not: legacy *fatals* when a non-uniform fibre field is missing,
where a missing `sigma0` is simply zero, so the spec will need a way for a
declaration to say it is required. That is deliberately not added yet, because
adding it without a law that uses it would leave it untested.

What is deliberately *not* reproduced: legacy's point-field `correct` aborts
outright when `sigma0` is non-zero. The framework evaluates points like any
other integration point, so that case now runs. That is a capability legacy
lacked rather than a behaviour change to a case that worked before.

### 8a.9 The bi-material interface: a difference we accept

`layeredPipe` is the only multi-material tutorial, and legacy and framework
disagree there by 6.25e-10 in a displacement field of 1.6e-7, about 0.4%. The
difference is confined to the four cells immediately against the material
interface, all on the inner side; the case's own discretisation error against
the analytical solution is 1.9%, some five times larger.

Three candidate causes were tested and ruled out:

* Switching `grad(D)` to a material-aware least-squares scheme changed the
  disagreement by nothing at all - 6.252e-10 either way. Both arms take their
  `gradD` from the same `mechanical().grad()` call, which computes it per
  material on sub-meshes, so the two share a gradient by construction and no
  gradient scheme can separate them.
* Tightening `rTol` from 1e-6 to 1e-12, which takes the solve from 57 to 402
  iterations, left it at 6.15e-10. It is not an artefact of stopping early.
  Note in passing that `rTol` is the entry that governs this criterion;
  `solutionTolerance`, which appears in these dictionaries, does not.
* Making the two materials identical collapses it to 1.1e-12, which is
  round-off. The material contrast is the whole of it.

What remains is the interface treatment itself. Legacy carries a sub-mesh per
material and iteratively corrects the displacement at the interface faces
between them. That machinery is what the sub-meshes exist for, and it is not
being reproduced: the intended replacement is a material-aware least-squares
gradient, where a cell's stencil draws only on cells of its own material. That
is simpler, it is consistent, and it will not agree to the last bit with the
iterative correction it replaces.

So this difference is expected rather than a defect, and it is recorded here
rather than chased. What would be a defect is it appearing somewhere with no
material interface, or growing.

Two notes for whoever writes the material-aware scheme. It does not exist yet:
an earlier `cellZoneInterface` helper in `leastSquaresS4fGrad` skipped face
contributions across cell-zone interfaces, and it was removed by the change
that added the rank-deficient and simplex stencil widening, so nothing on
`development` is material aware today. And that widening is itself a second
reason a migrated case may not reproduce a legacy number exactly, independent
of anything here. Whether the new scheme keys on cell zones or looks up the
materials is open; cell zones were what the removed helper used.

### 8a.5 Why kinematics and inputs stay separate, and what goes in inputs

`kin` and `inputs` look like the same thing: both immutable, both per-call,
both integration-point views the law only reads. The question of whether to
merge them is fair, and the answer is no, for two reasons and one caveat.

**The formulation boundary lives in `kin`'s type.** There are two `evaluate`
overloads, one taking small-strain kinematics and one taking finite-strain
kinematics, and a finite-strain-only law such as `neoHookeanElastic` refuses a
small-strain evaluation simply by not implementing that overload. That is a
compile-time capability check with no runtime convention behind it. Merging
would either lose it or force `smallStrainContext` and `finiteStrainContext`
types that duplicate the input half in both - and temperature and pressure have
no business participating in formulation selection.

**The finite-difference tangent reconstructs the kinematics object.** It builds
a fresh one per Voigt component, forwarding by hand everything that is not
being perturbed:

```c++
    smallStrainMechanicalConstitutiveLawKinematics kinPert
    (
        gradDPertView, kin.gradD0(), kin.dt()
    );
```

Anything inside `kin` must be forwarded at every reconstruction, and a
forgotten one is silent. Every new input would add a term to that constructor
and another way to drop a value without any test noticing. Passing the same
`inputs` object straight through each perturbation cannot be got wrong:

```c++
    evaluate(kinPert, inputs, shadow, respPert);
```

**The caveat, which is the part the question gets right.** `dt` is in the wrong
object today. It is not a deformation measure, and in the framework the only
references to `kin.dt()` are the two FD reconstructions forwarding it. It is a
per-call input, and a real one: the legacy `viscousHookeanElastic` needs the
time increment for its Maxwell update. So `dt` moves into `inputs` as a fixed
global scalar - not a one-element integration-point field - in the same
signature change, and the FD reconstruction gets shorter rather than longer.

The instinct that the boundary was arbitrary was therefore correct; it was the
placement of `dt` that was wrong, not the existence of two objects.

### 8a.6 What inputs will actually be needed

Taken from what the legacy laws look up today, not from imagination.

<!-- markdownlint-disable MD013 -->

| Input | Wanted by | Category |
|---|---|---|
| Mixed-formulation pressure `p`, `pf` | `neoHookeanElastic`, `GuccioneElastic`, `HolzapfelGasserOgdenElastic` when `pressureDisplacement_` | live |
| Temperature `T` | `thermoMechanicalLaw`, `viscousHookeanElastic` | live |
| Pore pressure | `poroMechanicalLaw` | live |
| Active tension `Ta` | `electroMechanicalLaw` | live, optional |
| `dt` | `viscousHookeanElastic` | global scalar |
| Reference stress via `sigmaf` | `neoHookeanElastic`, `linearElasticMisesPlastic` | prescribed - §4 covers it |
| Fibre and sheet directions | `GuccioneElastic`, `HolzapfelGasserOgdenElastic` | prescribed - §4 covers it |
| `DEqnA`, the momentum matrix diagonal | only `updateSigmaHyd` | solver internal - see below |

<!-- markdownlint-enable MD013 -->

Four things follow.

**The mixed-formulation pressure is not a future requirement.** `solvePressure()`
already exists in `linGeomTotalDispSolid` and both `nonLinGeom` solvers, and
three laws already look `p` up. An `inputs` object is needed for the *next* law
migrated, not for some later coupling.

**Two different physical quantities are both called `p`.** The mixed-formulation
pressure is the solid model's own unknown; the pore pressure is a fluid
pressure that may come from another region entirely. They must be given
distinct semantic names in the declaration, or a case will one day supply one
where the law wanted the other and the answer will merely be wrong.

**The smoothing decision removes the worst coupling rather than relocating it.**
`DEqnA` is a law reaching into the linear system's diagonal, and it appears
only inside `updateSigmaHyd`. Since that moves to the solid model, `DEqnA` does
not need to become an input at all.

**Do not design for contact and cohesive laws.** Face normals, gaps and
characteristic lengths are wanted by contact and cohesive-zone models, but
those are a separate hierarchy - cohesive laws live under `dynamicFvMesh`, and
contact models take normals and gaps explicitly - and they are interface laws,
not bulk constitutive laws. Pushing that data through every cell law's loop to
serve them would be the obvious overreach here. They get their own interface
kinematics if and when they are migrated.

Plausible but not present in the repository, and therefore not to be designed
for yet: concentration or chemical potential for swelling, a damage phase
field, crystal orientation. Of these, only the phase field would be a live
input; orientation is prescribed material data that §4 already covers, and a
phase-field history variable is persistent state, not an input.

### 8a.7 Shape of the inputs object, and the response side

**Declared and name-keyed, not fixed members.** Fixed members would be pleasant
for `dt` alone and would fail on the first new coupling, because the type is
instantiated at every evaluation site and in every manager call path, internal
and boundary. A law declares what it needs - semantic name, type, dimensions,
entity, required or optional, and the mapping policy for shared and boundary
points - and reads it with typed getters.

**Absence must not read as zero.** `Ta` absent means "use the law's own
prescribed activation"; `Ta = 0` means explicitly inactive. A required `T` or
`p` that was never supplied must fail before evaluation rather than silently
behave as zero. That argues for a strict `get` and a separate `find` for
genuinely optional inputs, which is the distinction the state object already
makes.

**The response side needs the same treatment eventually, but not yet.** A
phase-field law must return a driving energy; a dissipative law may need to
report dissipated work. Those are outputs of the current evaluation and must
not be pushed into the state, where they would inherit exactly the shadow and
rollover hazards §8a.4 rejects. The response object is already the right
container and is already passed to every law, so adding declared outputs later
does not force a second signature migration - this can wait for a consuming
solver.

One rule is worth recording now, because it is the same trap in a new place:
**a finite-difference perturbation must never write auxiliary outputs into the
real coupled fields.** It evaluates the law repeatedly at perturbed kinematics,
so any such output must go to disposable storage or be explicitly suppressed,
exactly as the shadow state does for history.

### 8a.4 The rejected design, and why it is worth recording

The first draft proposed carrying `T` and `p` in the state object: the solid
model would push them in through `manager.setCouplingField("T", T)` before each
update, the law would read them like any other state field, and `evaluate`
would keep its signature. A "was this set during this update?" flag would catch
a solid model that forgot.

It cannot work, for a reason that is specific and checkable rather than
aesthetic. **A shadow state aliases the parent's old-time fields and owns its
own, fresh, current-time fields.** That is what makes a tangent query safe. But
it means a coupling value written into current-time state is simply *not there*
inside the shadow: every tangent query, and every finite-difference
perturbation, would evaluate with a missing input. Putting it in old-time state
instead is worse - it disguises a current input as history, and `storeOldTime()`
would then commit it.

The staleness flag was the tell. It was compensating for the wrong ownership
model: the design needed a protocol enforced by convention, outside the type
system, to keep a mutable manager-held copy correct. The inputs object removes
the need for the protocol instead of policing it. **Where a design needs a rule
that says "remember to call this first", the ownership is usually wrong.**

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

## 11. Writing state: two purposes that are not the same purpose

The first draft of section 9 treated writing state as one problem and, when
quadrature-point topologies would not fit, proposed refusing them. That was
wrong, and refusing is at best an interim. There are two reasons to write a
state field and they want different things:

1. **Restart.** The value read back must be the value written, exactly. An
   interpolated field restarts a *different* calculation that happens to look
   similar.
2. **Visualisation.** ParaView can display a cell field and a point field. It
   cannot display a surface field, and it certainly cannot display a value at a
   quadrature point. Here an approximation is not merely acceptable, it is the
   only thing on offer.

Conflating them is what produced the refusal. Separated, both are solvable, and
quadrature stops being a special case that has to be excluded.

### 11.1 Restart: write the state in the shape the state is already in

State is stored as a flat `Field<T>`, one entry per integration point. That is
already exactly what an `IOField<T>` holds, so for any topology whatsoever the
state can be written verbatim and read back in the same order, entry for
entry, and to full precision.

Not bit-for-bit by default, which is worth being accurate about: written as
ASCII a value is rounded to `writePrecision`, and at the default of 6 that
rounding alone moved a continued run by 6e-6, a thousand times the 1e-8 the
restart otherwise reaches. So the state is written at 17 significant digits
whatever the case asks of its other output. `writePrecision` is a choice about
numbers a person is going to read, and nobody reads this file.

What is left at the default is the ordinary fields, `D` among them, still
written at the case's precision - and that is shared with the legacy model,
which lands on the same figure to the same last digit. It is a property of
OpenFOAM's ASCII output rather than of either model. Quadrature points are
not awkward here; nothing in an `IOField` needs the entries to correspond to a
mesh entity.

So the fallback is total, and there is no topology this design cannot restart.

It was worth asking whether, where the integration points *do* correspond
one-to-one with a mesh entity, the natural `GeometricField` should be preferred
instead, for one concrete reason: `decomposePar` and `reconstructPar` know how
to map it, and an `IOField` they will not touch. Counting points against mesh
entities:

| topology    | points          | exact field  | ParaView | decomposes |
|-------------|-----------------|--------------|----------|------------|
| cellCentred | `nCells`        | vol          | yes      | yes        |
| pointCentred| `nPoints`       | point        | yes      | yes        |
| faceCentred | `nFaces`        | surface      | **no**   | yes        |
| compactCell | per-cell count  | none         | no       | no         |
| compactFace | per-face count  | none         | no       | no         |
| dualFace    | internal duals  | none         | no       | no         |

The first draft of this section made a rule of that - the natural geometric
field where the correspondence is exact, an `IOField` otherwise - and what is
built is **not** that. Every topology, `cellCentred` included, is written as an
`IOField`. The reasons are in 11.1a, and the rule above is left standing only as
the question it answered.

What the one path costs is that state reads back only on the decomposition that
wrote it, for every topology rather than the three that would have mapped. That
is a real limitation. It fails loudly on a mismatch rather than reading a
wrong-length list, because a restart that silently differs is worse than one
that refuses, and 13 is about removing it.

### 11.1a Why one path, and not the geometric field where it fits

Three things, the first of which is decisive and was not known when the rule
above was written.

**A `calculated` boundary patch does not write its values.** It writes its type
and nothing else, so boundary state parked there would be lost on the way out,
before decomposition is ever reached. `zeroGradient` is worse: it re-extrapolates
from the internal cells on every evaluation. Carrying boundary state in a
`GeometricField` boundary therefore needs a value-bearing patch type chosen and
justified per state variable, which is a decision the law would have to make and
no law has any business making.

**Per-law boundary state is a subset of a patch, not a patch.**
`lawBoundaryFaces_[lawI][patchI]` holds only the faces whose owner cell belongs
to that law. A `GeometricField`'s boundary is dense - one entry per patch face -
so a multi-material case would have to pad every patch with values belonging to
no law and then be careful never to read them. The same padding problem recurs
in the interior, where a per-law volField means something in that law's cells
and nothing in the rest, in a field any generic tool is free to average.

**No law declares a dimensionSet.** A `GeometricField` needs one. Adding it to
every declaration to satisfy a storage format is the format dictating the
interface.

Against that, the gain would have been a redecomposition path for three
topologies out of six, of which one is in use. It is not worth two read paths,
two failure modes, and a boundary convention per variable - especially as 13
gets the redecomposition for all six without any of it.

### 11.2 Visualisation: a projection, written, never read

Four of the six topologies produce nothing a user can look at. For those, and
only for looking at, project the state onto a `volField` - cell average of the
points belonging to each cell - and write that under a name that says what it
is, so nobody mistakes it for the state itself.

Two rules keep this honest. It is **written and never read**: the restart path
must not fall back to it, or the exactness above is lost at the first
inconvenience. And it is **opt-in**, because for a law with several state
variables on a fine mesh it is a lot of output that most runs do not want.

This also covers `faceCentred`, which restarts exactly through a surface field
but still cannot be viewed through one.

### 11.3 What this changes about the staging

Section 9 staged cell-centred first because it was the easy one. That still
holds - it is the default topology and the one every current tutorial uses - but
it is now the first case of a general rule rather than the only case that works,
and the quadrature topologies are scheduled rather than excluded.

## 12. What building it changed

### 12.1 Legacy restart does work, at least here

The plan was written on the understanding that restart was never properly
considered in the legacy approach and could be assumed not to work, and that
this freed the design from having to reproduce it. Measured on
`perforatedPlate`, stopping at t = 10 and continuing, the legacy mechanical
model reproduces the uninterrupted run to 1e-12. It works.

So the freedom was not there to take, and it is as well that nothing was
built on it. The framework has to restart at least as well as legacy does,
which is a stricter target than "better than nothing", and it is what the
comparison in the regression test now asserts.

### 12.2 The last error was not a state error at all

With the state written and read back, the restart still differed by 3e-4 in
max epsilonEq while every state field was byte-identical after the first
restarted step. The state was never the problem.

`impK` is formed once and kept. On a cold start it is formed before anything
has happened; on a restart it was being formed against the restored history,
and came out different. The solver stops on a residual measured relative to
its first one, so a different `impK` moves where the step stops - the run
converges to a slightly different point and, the material being path
dependent, keeps the difference.

The legacy `impK()` looks like it has the same dependence, and its comment
here used to say so. It does not: it scales by `1 - 2*mu*DLambda/magSTrial`,
and `DLambda` is `NO_READ`, so on a restart it is zero, the factor is exactly
one, and the result is the elastic tangent either way.

`frameworkImpK` now evaluates at zero gradient against a state built the way a
fresh run builds one, which is the elastic tangent and is the same whether the
run was continued or not. A cold start is unaffected, because there the state
already is that state, and no tutorial moves. Nothing is lost: this is a
scalar preconditioner for an approximate Jacobian, and the elastic value is in
practice as good as a scaled one.

The general point is worth keeping. Restoring state exactly is necessary and
was not sufficient: anything else frozen at construction is also part of the
restart, and a quantity that is *derived* from state inherits the state's
restart problem without ever appearing in the state.

### 12.3 Restart is a property of where a run started

The first version asked `timeIndex() > 0`. A topology is created the first
time something asks for it, which can be long after the run has advanced, so
an ordinary run became a "restart" the moment it took a step, and then refused
to find files it had never written. The unit tests caught it by advancing a
time step, which no amount of reading the code had.

It asks `startTimeIndex()` now, which is where the run began rather than where
it has got to.

### 12.4 What the tests assert, and why the middle one matters

An elastic restart is exact to 1e-13 with no state IO whatsoever, because the
history being restored is zero either way. A restart test that does not
establish that the material has yielded first is therefore vacuous, and the
test asserts the yielding rather than assuming it.

The negative control is that deleting one state file makes the run refuse.
Without it the test would only show that two runs agree, not that they agree
*because* the history came back.

## 13. Changing the decomposition

State reads back only on the decomposition that wrote it. `decomposePar` will
not copy these files into the processor directories, so the continued run finds
nothing and refuses. Honest, and not good enough: a user who redecomposes a
history-dependent case has to start it again.

### 13.1 The routes that do not work, and why

**Registering with the stock utilities.** There is no hook. The list of classes
`decomposePar` and `reconstructPar` map is hard-coded in the utilities, and the
processor `Time` databases they open are deliberately created with no function
objects and no extra libraries, so nothing this library does can be reached from
inside them. Ruled out as a fact about OpenFOAM rather than as a preference.

**Storing the state in a `GeometricField` so the utilities map it.** Covered in
11.1a. The boundary half does not survive its own write, and the per-law subset
does not fit a dense patch.

**Becoming a Lagrangian cloud.** This is the one thing in OpenFOAM that already
survives redecomposition, and it survives because a particle carries its own
identity: a position and a cell, from which `decomposePar` recomputes ownership
without any ordering convention. That is the right idea and the wrong vehicle -
using it means every integration point becomes a real particle with a valid
position and tet-location that this library then has to keep valid, and
foam-extend's cloud IO supports fewer field types than the others. The idea is
worth stealing; the machinery is not worth inheriting.

### 13.2 The route that does work

`decomposePar` writes `cellProcAddressing` into every processor directory, and
it is exactly the map wanted:

    localState[procCelli] = serialState[cellProcAddressing[procCelli]]

which is the same direction the official field reconstructor uses. So a
processor restarting can read the *serial* state file and permute it itself.
Nothing about `decomposePar` changes; the work happens in this library's read
path, which is where it belongs, and the case is decomposed with the stock
utility exactly as it is today.

This reaches further than the `GeometricField` route would have - but only
`cellCentred` is built. A cell maps through the addressing directly, and that
covers every tutorial there is. The rest is scoped and not written: a topology
whose points are grouped by cell, the compact ones included, can map the same
way once the per-cell counts are carried alongside, and the topology knows them
because it built the offsets; a point-based one needs `pointProcAddressing`
read and decoded. Those topologies write no locations at all today, so a
changed decomposition refuses rather than guessing, which is the right way for
an unbuilt thing to behave.

What has to be true, and has to be checked rather than assumed:

  - the serial state must still be there, and reachable from every processor;
  - the serial numbering must be the numbering the addressing file was written
    against, so a renumbered or refined mesh invalidates it;
  - the per-cell offsets must be reconstructible on the processor mesh, which
    is a stronger requirement than the flat length check that guards the
    same-decomposition case today.

Each of these should refuse rather than guess. The failure this removes is one
that currently refuses too, so nothing gets quieter: it gets possible.

### 13.3 One case that would be quiet, and should not be

If a case is redecomposed but stale processor time directories are left in
place, a processor can find a state file that exists, has the right header, and
happens to have the right length, and read another decomposition's values. The
length check cannot see it. That is the one path in this design where a wrong
restart would not announce itself, and it wants an identity in the file - the
decomposition it was written for - rather than a warning in a README.

## 14. The half restart

Restoring the constitutive state is not the whole of a restart, and the part
that is missing does not announce itself.

`stripFooting` - poroMechanicalLaw over linearElasticMohrCoulombPlastic, the
only shipped case where a composite's child keeps history of its own - came
back 46% wrong. Every state file was written and read correctly. The state was
never the problem.

The Mohr-Coulomb law is incremental: it builds this step's stress on the last
step's, over the strain increment `symm(gradD) - symm(gradD0)`. The solid model
writes `grad(D)` at old time only when asked, through a `restart` switch that
defaults to off, so `gradD0` came back as zero and one step's increment became
the whole run's strain - seventy-six times too large. The restored history was
then added to a trial stress that had already counted it, which is why
restoring the correct history looked *worse* than restoring none at all, and
why the bisection pointed at the state when the state was innocent.

Two reviews looked at this and both diagnosed it wrongly, one as an old-time
aliasing bug and one as a contaminated write, before an instrumented print of
what the law actually read settled it in one run. The lesson is not about the
reviewers. It is that a plausible mechanism is not evidence, and the cheap
measurement was available the whole time.

### 14.1 What guards it now

A restart that would restore constitutive history without the kinematic history
that history is measured against is refused, by name, with the switch to set.
The framework knows when it is doing the first half of a restart, so it is the
right place to notice that the second half is missing.

The solid model warns in the general case, where it cannot know whether the
material is incremental. A law written in total strain genuinely does not care,
and `perforatedPlate` restarted correctly for years' worth of steps without the
switch for exactly that reason, so refusing outright would stop runs that are
fine.

### 14.2 The switch itself

It is worth saying plainly that the switch is the problem, not the fix. It
guards nine fields, which is a real cost on a large case, so something like it
has to exist. But it defaults to off, twenty-four of a hundred and sixty-eight
tutorials set it, and getting it wrong costs tens of percent in silence. A
default of on, with an opt-out for cases that will never be restarted, would
put the cost where the choice is rather than the correctness.

That is a library-wide change and not this work's to make.

## 15. The other direction, and one numbering for both

13.2 solved half the problem: a case run in serial and then decomposed. The
commoner half is the reverse - run in parallel, `reconstructPar`, continue in
serial - and it was left refusing, because `reconstructPar` gathers these files
no more than `decomposePar` distributes them.

Both halves turn out to be the same problem said twice, and they are now the
same code. A run that finds no state of its own looks where the other shape of
the run would have left it: a decomposed run reads the undecomposed case and
distributes, a serial run reads the processor directories and gathers. Each
entry is matched by the piece of mesh it belongs to, so neither has to know how
the mesh was cut.

### 15.1 What made it one problem instead of two

The locations are written in the numbering of the *undecomposed* mesh, whichever
way the run that writes them is being run. A processor translates its own
through the addressing `decomposePar` left it, before writing. So a serial run's
locations are what a decomposed run needs to distribute the state, a decomposed
run's are what a reconstructed run needs to put it back together, and one list
serves both.

That forced a change to how a location is said. It used to be a cell index, or
a face index offset past every cell - and a processor holds a piece of the mesh,
so it does not know how many cells the whole one had. An encoding that needs a
number nobody has cannot be written from both ends. A face is now its index made
negative and offset by one, so that face zero is still negative, and nothing
needs a total.

### 15.2 The face that two processors both hold

A face `decomposePar` cut is a boundary on two processors and interior in the
undecomposed mesh. Both write a value for it; neither is ever asked for one,
because a cell-centred topology keeps no history on an interior face. So the
gather sets rather than inserts, and the duplicate is harmless by construction
rather than by luck.

Going the other way the same face is the one case with no answer, and it takes
its owner cell's history. Only a *processor* patch may: a cyclic is coupled too
and exists in the undecomposed mesh exactly as it does here, so a miss on one
is a real disagreement and has to be reported rather than smoothed over.

### 15.3 What it costs

Both directions reproduce the uninterrupted run to about 1e-8 on
`perforatedPlate`, which is where the solver's own restart sits. Neither
`decomposePar` nor `reconstructPar` is modified, and neither needs to be.

What is still missing is the same as before: only the cell-centred topology
writes locations, so only it can change shape. The others refuse, which remains
the right behaviour for something unbuilt.


### 15.4 What each direction refuses

A source being self-consistent does not make it this case, so both directions
check that it describes this mesh before believing it. A decomposition
preserves cells, which gives a count that has to agree either way: the
undecomposed file records the cells its mesh had, and a decomposed run holds
that many between its ranks; a set of processor directories records what each
rank held, and those have to add up to a reconstructed run's mesh.

Collated output is refused by name. It puts every rank's data in one
`processors<N>` directory rather than a `processorN` each, in a format this
does not read, and a run that simply found nothing there would look exactly
like one that had never been run in parallel.

Going from N processors straight to M, without reconstructing on the way, is
not supported and says so: the decomposed branch looks for an undecomposed
state that a never-serial case never wrote. Reconstruct and decompose again,
both of which now work.

## 16. What every law that keeps history is now tested for

Five laws keep history, and each has a restart test:

| law                            | case            | what it exercises            |
|--------------------------------|-----------------|------------------------------|
| linearElasticMisesPlastic      | perforatedPlate | plastic history; both        |
|                                |                 | decomposition directions     |
| linearElasticMohrCoulombPlastic| stripFooting    | history inside a composite   |
| poroMechanicalLaw              | stripFooting    | the composite's own history  |
| viscousHookeanElastic          | viscoTube       | a variable number of fields  |
| neoHookeanElasticMisesPlastic  | neckingBar      | finite strain, updated       |
|                                |                 | Lagrangian                   |

Every one of them checks three things and not one: that the state is written,
that deleting it stops the run, and that the answer comes back. The middle
check is what the other two rest on. Twice now a restart test has passed while
proving nothing - once because the material had not yielded before the restart
point, once because a viscous arm count was read with a glob that matched
nothing - and in both cases the assertion caught it rather than the comparison.

### 16.1 Two things these tests found that were not the framework's

`viscoTube` could not restart at all. `normalDisplacement` read its value under
that keyword and wrote it as `normalDisp`, so every time directory it produced
was one the case could not be started from: the field came back and the
condition it carried did not. Fixed, because a restart test cannot run on a
case that cannot restart - and then `normalDisplacementZeroShear` turned out to
have the identical defect, which a review found and this work's own search for
other readers had missed by filtering out the very keyword it should have been
looking for.

`neckingBar` does not reproduce an uninterrupted run exactly across a restart,
with or without the framework. On the axial force the legacy model differs by
1.3e-5 and the framework by 4.3e-6, and a force component that is otherwise
zero returns as -0.041 in both. That is the updated Lagrangian solid model's,
not this framework's, and it is left alone and written down rather than tuned
around. It is why that arm holds to 1e-4 where the small-strain arms hold to
1e-6, and why a negative control carries the weight there instead.

## 17. The fibre laws, and what porting them needs decided

`problem3` is now in the tree with a check of its own, and it is the only
tutorial whose material is anisotropic about a direction. It runs
`electroMechanicalLaw` over `GuccioneElastic`, neither of which is ported. This
records what the port needs decided, so the decisions are made rather than
discovered halfway through.

### 17.1 What the law holds

`GuccioneElastic` keeps four fields that are not history and are not read from
a dictionary: `f0`, `s0`, `n0` - the fibre, sheet and normal directions - and
`R_`, the rotation built from them. `f0` is produced by `setFibreField` before
the solver runs and read from the `0` directory.

That makes them **prescribed** state in this framework's terms, and the first
prescribed *vector* state anywhere: everything prescribed so far has been a
scalar or a symmetric tensor. The path exists and is untested for vectors,
which is a reason to port this law rather than an obstacle to it.

`R_` is the awkward one. It is derived from the three directions and constant
in time, so it is neither history nor input, and the choices are to recompute
it per evaluation, to hold it as a fourth prescribed field, or to give the
framework somewhere to put a per-material derived constant. The third is new
machinery and should not be invented for one law; the first is arithmetic in the
inner loop. It is probably the second, but it is a decision and not an obvious
one.

### 17.2 What the composite holds

`electroMechanicalLaw` wraps a passive law and adds an active tension along the
fibre direction, ramped over `rampTime`, or read from a `Ta` field on the
registry if one is there. The registry lookup is exactly the coupling input the
framework already has a path for, and the ramp is a function of time, which
`mechanicalConstitutiveLawInputs` already carries.

So it is the same shape as `poroMechanicalLaw`: a law with its own parameters
holding a child state for its sub-law. That one is ported and tested, including
its restart, so this has a worked example to follow rather than a blank page.

### 17.3 What has to be decided first

  - Whether `calculateStressInLocalCoordinateSystem` is ported at all. It
    doubles the law's evaluation path, and `problem3` sets it to `no`. A port
    that carries both has twice the surface and half the coverage.
  - Whether the mixed displacement-pressure variant is in scope. The tutorial
    ships both formulations and the legacy law branches on which is active.
  - Whether `R_` is prescribed state, recomputed, or something else.

None of these is a question the tests can answer, which is why they are written
down here rather than guessed at.

## 18. The mixed displacement-pressure formulation

A review, not a decision taken. Nothing here is implemented.

### 18.1 What the legacy laws do

`neoHookeanElastic` and `GuccioneElastic` each read a `pressureDisplacement`
switch from their own dictionary, and when it is set three things change:

  - the deviatoric stress is computed by a different formula,
    `mu*(b - I)/J` rather than `mu*dev(bEbar)`;
  - the volumetric part stops being computed and is taken from the solid
    model's pressure instead, reached through
    `mesh().lookupObject<volScalarField>("p")`;
  - `impK()` returns `mu` instead of `(4/3)mu + K`.

So a material dictionary entry changes which equations the law solves, where
half its answer comes from, and what tangent the solver is given.

### 18.2 Two things wrong with that, before any framework question

**The switch is in the wrong place.** Whether the pressure is a solved variable
is a property of the solution algorithm, not of the material. Steel does not
become a different steel because the solver chose a mixed formulation. In a two
material case the switch has to be set identically in every material's
sub-dictionary, and setting it in one and not the other is accepted silently
and produces a mixture of formulations.

**The stress conventions differ between the two branches, and only work by
accident of both being written the same way.** The ordinary path ends with

    sigma = (1/J)*(sigmaHyd*I + s)

so `s` is a Kirchhoff-like quantity and the division by `J` is what makes the
result Cauchy. The mixed path ends with

    sigma = -coeff*p*I + s        with s = mu*(b - I)/J

where `s` already carries its own `1/J` and the pressure is applied with none.
Both are Cauchy, so both are right, but nothing says so and the two lines do
not look like they agree.

That matters because of where the error would land. A law returning a Kirchhoff
deviatoric stress to a solid model subtracting a Cauchy pressure is wrong by a
factor of `J` - and a mixed formulation is used precisely when the material is
nearly incompressible, which is precisely when `J` is nearly one. The mistake
would be small, smooth, and entirely plausible. It is the kind that survives a
review and a plot.

### 18.2a Two things a fact-check changed

Both were wrong in the first draft of this section and both matter.

**What the mixed branch returns is not deviatoric.** `mu*(b - I)/J` has a
non-zero trace at finite deformation, whatever the comment above it says. So
"the law returns the deviatoric stress and the solid model adds the pressure"
is not a description of what happens now, and if it were implemented literally
the trace of that term would be counted twice - once by the law and once by the
pressure. Whatever the law returns has to be named precisely, and if the name is
"deviatoric" then something has to make it true.

**The pressure is not the pressure.** `p` is the mixed solver's variable, and it
enters the stress as `-pressureDisplacementCoeff*p*I`, where the coefficient
defaults to `3*nu/(1 + nu)` computed from the material's own `K` and `mu`. It is
one only in the incompressible limit; for `nu = 0.3` it is about 0.69.

That second fact is the awkward one for assembling at the solid model, because
the coefficient belongs to the material and the solid model would have to be
told it per material. It does not sink the idea, but it does mean the solid
model cannot simply subtract `p*I`.

### 18.3 What is proposed

There are two shapes this can take, and 18.2a is what decides between them.

**A: the solid model assembles.**

    sigma = s_dev(from the law) - alpha*p*I

The law returns a genuinely deviatoric stress; the solid model owns the
pressure and the assembly. It has to be given `alpha` per material, since that
is material data, which means the manager supplying a per-material coefficient
alongside the stress - new machinery, though small.

**B: the pressure is a coupling input and the law assembles.**

The solid model publishes `p`, the manager hands it to each law as a live
coupling input, and the law returns the whole stress as it does now. This is
exactly what `poroMechanicalLaw` already does with pore pressure, so it needs no
new machinery at all, and `alpha` stays where it is computed from `K` and `mu`.

B is the smaller change and fits what exists. A is the cleaner separation - the
material stops knowing that a pressure equation is being solved - but it needs
the coefficient plumbed and the deviatoric part made genuinely deviatoric.

I lean to B for the port, and A as the eventual shape, on the grounds that B
can be built and tested now while A needs the trace question in 18.2a settled
first. Either way the following three hold.

**The request comes from the solid model, not the material.** A solid model
running a mixed formulation says so once, and the manager asks every law for a
deviatoric stress. No material dictionary mentions it, and a case cannot be
half mixed.

**The convention is stated and not implied.** The framework's stress is Cauchy.
That is true today and written down in one comment; it should be part of the
law interface's contract, so that "deviatoric stress" has one meaning and a law
returning a Kirchhoff quantity is a violation rather than a coincidence waiting
to matter.

**The tangent follows the same split.** A mixed solid model needs `mu`, not
`(4/3)mu + K`, and that is the deviatoric tangent of the same request - not a
second switch to be set consistently with the first.

### 18.4 Why this has to be settled before the fibre port

`GuccioneElastic` carries the branch, and `problem3` ships both formulations.
Porting it as it stands means carrying `pressureDisplacement` into the
framework and giving the first ported law a material-dictionary entry that
selects a solution algorithm - which is the thing this framework exists to stop
laws doing.

The alternative is to port the displacement formulation only, leave the mixed
one to the legacy law until the above is built, and have `problem3`'s
regression test cover the displacement arm alone for now. That is the smaller
step and it does not foreclose anything.

### 18.5 What the reviews changed, and one thing they found in the code

**The convention error is not hypothetical. It is in the tree.**
`GuccioneElastic` computes, in its ordinary path,

    s = dev(symm(F & S & FT))/J          // Cauchy

and in its mixed path,

    s = dev(symm(F & S & FT))            // no /J, so Kirchhoff

under a comment that says "convert the second Piola-Kirchhoff stress to the
deviatoric Cauchy stress", and then forms `sigma = s - p*I`. That is a Cauchy
pressure subtracted from a Kirchhoff deviatoric stress: the very thing 18.2
described as the mistake that would be small, smooth and plausible, written
down and shipped, with a comment asserting the opposite.

How much it matters is a separate question from whether it is right. The error
is a factor of `J` on the deviatoric part, and a mixed formulation exists to
drive `J` to one, so for the nearly incompressible tissue this law is used on it
is small - which is exactly why it would survive. `neoHookeanElastic`'s mixed
branch keeps its `/J`, so the two laws do not agree.

**The two laws do not agree on the coefficient either.** `neoHookeanElastic`
scales the pressure by `3*nu/(1 + nu)`; `GuccioneElastic` has no coefficient at
all and subtracts `p` bare. So the "one meaning" a new contract has to fix does
not exist to be copied - it has to be chosen, and one of the two laws will
change when it is.

**The mechanism already exists, and it is not the one proposed.** The manager
injects `planeStress` and `solutionD` into every law's dictionary at
construction and refuses to start if a law's own sub-dictionary sets either.
That is precisely "the request comes from the solid model, not the material,
and a case cannot be half mixed", already built and already tested. A mixed
formulation flag belongs there, injected the same way, rather than in new
machinery.

And `tangentRequest::scalarDeviatoric` already exists, with a comment naming
mixed displacement-pressure formulations as its purpose. The tangent half of
this needs wiring, not designing.

**Keep the two signals apart.** Whether a law is in mixed mode is fixed for the
run; whether a tangent is wanted is decided per evaluation. Collapsing the first
into `tangentRequest` would make "deviatoric stress, no tangent" - an ordinary
residual iteration - inexpressible.

**Composites need the contract stated more carefully than "deviatoric".**
`poroMechanicalLaw` hands its sub-law a stress, gets a full one back, and
subtracts a pore pressure. Wrapping a mixed-formulation law would produce
something that is deviatoric except for a pore pressure it has already netted
out and a formulation pressure it has not, which is not a thing the word
"deviatoric" describes. And both default their field name to `p`, so a
poro-wrapped mixed law would find one field and count it twice.

**A trace check is not enough to catch it.** Scaling a traceless tensor by
`1/J` leaves it traceless, so asserting `tr(s) == 0` catches
`neoHookeanElastic`'s non-trace-free term and misses the `J` error entirely.
What would catch it is comparing a law's deviatoric evaluation against
`dev()` of its own ordinary evaluation, at a state deliberately away from
`J = 1`, as a unit test per law. That argues for keeping the ordinary path in
any law that gains a mixed one, as the oracle.

**The tangent is not `mu` outside isotropy.** `GuccioneElastic`'s mixed tangent
is an effective shear modulus computed by a dedicated routine, because the law
is anisotropic. Which is a further reason to do 18.4's smaller thing: port the
displacement formulation, leave the mixed one on the legacy law.

## 19. Where the mixed formulation actually belongs

18 asked whether the solid model or the law should own the assembly. The
answer is the solid model, for a reason 18 did not give: if the law chooses
whether to use the pressure, then a solid model can solve a pressure equation
and a law can quietly not use it, and nothing says so. The law should not have
the choice.

And it already works that way in one solid model. `linGeomTotalDispSolid`, when
its `solvePressure` switch is on, does

    sigma() = dev(sigma()) - p*I;

which takes whatever the law produced, discards its volumetric part and
substitutes the solved pressure. The law is not consulted, cannot opt out, and
needs no switch of its own. It also asks for `tangentRequest::scalarDeviatoric`
in that case, which is what that enum was declared for.

So the framework needs nothing. The two fibre laws just ported carry no mixed
branch and do not need one. The `dev()` disposes of the question of whether
what a law returns is trace-free, because the solid model makes it so.

### 19.1 What harmonising costs, which is not nothing

The two solid models do not solve for the same variable, and this is the part
that has to be got right rather than assumed.

Both solve for something that is approximately `-kappa*div(u)` rather than the
physical Cauchy pressure. `coupledPressureDisplacementSolid` then carries a
coefficient `alpha = lambda/kappa = 3*nu/(1 + nu)` into the momentum coupling,
as `nuCoeff_`, and `neoHookeanElastic`'s mixed branch scales the pressure by
the same `alpha`. They are a matched pair, not two laws disagreeing.
`GuccioneElastic` needs no coefficient because its mixed branch forces the bulk
modulus to `GREAT`, so `nu` is one half and `alpha` is one - which is also why
that branch can omit a `/J` that the ordinary branch keeps.

Converting `coupledPressureDisplacementSolid` to the replacement form is
therefore a change of variable, from `q` to `P = alpha*q`, and it has to be
carried through consistently: the pressure equation becomes `P/(alpha*kappa) +
div(u) = 0`, the momentum coupling loses its `alpha`, and the Rhie-Chow
coefficient transforms with it. Pressure values and pressure boundary
conditions change meaning. It is not a refactor.

### 19.2 It is already a finite-strain formulation

An earlier draft of this section said `solvePressure` was small strain only,
because `linGeomTotalDispSolid`'s constraint is `tr(gradD())`. That is true of
that model and it is not a limitation: it is the linear geometry model, and the
small-strain volume change is the right constraint there.

`nonLinGeomTotalLagTotalDispSolid` implements the same switch at finite strain,
with the same replacement and the constraint written on `J`:

    sigma() = dev(sigma()) - p()*I;
    ...
    - p*rKappa() + ... - 0.5*(J^2 - 1)/J

So both halves exist and agree, and putting a finite-strain law behind the
switch is a matter of using the model that already does it. There is no piece
of work here that has not been done.

### 19.3 The coefficient is not a per-material problem after all

18.2a worried that assembling at the solid model needs `alpha` plumbed through
per material, because it is computed from the material's own `K` and `mu`.
`coupledPressureDisplacementSolid` already answers that: it computes
`nuCoeff_ = 3*nu/(1 + nu)` itself, as a field, from the generic modulus
quantities every law already exposes through `impK()` and `rKappa()`. Nothing
is plumbed and multiple materials are handled by it being a field rather than a
number.

That removes the main objection to the solid model owning the assembly, and it
is worth noting that the objection was answered by code that already existed
rather than by an argument.

### 19.4 Two things that are actually wrong today

**The formulation pressure and the pore pressure are the same field.**
`solidModel::makeP()` registers its pressure as `"p"`, and
`poroMechanicalLaw`'s `pressureFieldName` defaults to `"p"` in both the legacy
and the ported law. A solid model with `solvePressure` on, running a material
that uses `poroMechanicalLaw` without overriding that name, has one field
serving as two different pressures: either the Biot term counts the formulation
pressure or the replacement counts the pore pressure. Renaming the solid
model's own field is the cheap fix, since the two concepts should not be able
to collide by accident.

**`dev()` is not the volumetric-deviatoric split for an anisotropic law.**
Removing the trace of the Cauchy stress is not the same operation as taking the
part of the strain energy that is independent of volume change, and for a fibre
law the two differ. That is why the legacy `GuccioneElastic` needed a dedicated
effective shear modulus for its mixed tangent rather than a scalar. The
solid-model-owned replacement is blind to this by construction: it deliberately
does not ask the law anything, so it cannot know when the split it performs is
the wrong one.

This is the real cost of the design, and it is worth stating plainly rather
than discovering later. It is not specific to finite strain - it would show up
at small strain too, with an anisotropic law. No shipped case exercises
`solvePressure` with an anisotropic law, so it is latent rather than observed.

So the order is: the design is settled and needs nothing from the laws; the
formulation pressure wants its own name; and converting the coupled solver is
the one real piece of work, with the change of variable in 19.1 done
deliberately rather than discovered.

## 20. What `dev()` needs from a law, and why it is not about anisotropy

19.4 said the solid-model replacement is blind to whether the split it performs
is the right one, and put that down to anisotropy. That is not quite the
criterion, and getting it right decides what the law has to provide.

### 20.1 The criterion is how the energy is written

For a law whose energy is written with the split,

    W = W_iso(Cbar) + U(J),      Cbar = J^(-2/3) C

the isochoric contribution to the Cauchy stress is trace-free by construction,
so `tr(sigma) = 3*U'(J)` and `dev(sigma)` recovers the isochoric stress
*exactly*. Replacing the rest with `-p*I` is then not an approximation at all.

For a law whose energy is written on the full deformation, `dev(sigma)` still
removes the spherical part - nothing is counted twice, and the first draft of
this section was wrong to say so. What happens instead is subtler and worse to
find: what is left is a deviatoric response that still depends on `J`, and
replacing the spherical part with an independent `p` gives *a different
material model* from the one the law describes. Not obviously wrong, not
double-counted, simply not the law that was asked for.

The sign is worth writing down, since the two conventions differ by one:
with `sigma = dev(sigma) - p*I`, the pressure is `p = -U'(J)`.

So the question is not isotropy. It is whether the energy was parameterised on
`Cbar` or on `C`. An anisotropic law written on `Cbar` is reproduced exactly;
an isotropic law written on `C` is replaced by a neighbouring model.

Among the ported laws, `neoHookeanElastic`, `MooneyRivlinElastic` and
`neoHookeanElasticMisesPlastic` are written on an isochoric measure and are
reproduced exactly. `StVenantKirchhoffElastic`, `OgdenElastic` and the
`GuccioneElastic` just ported are written on the full strain and are not -
Guccione's `Q` takes the whole Green-Lagrange strain, with the volumetric
response added separately as a penalty, and Ogden decomposes the full `C` and
applies `dev` only afterwards.

That last one is worth noticing: the legacy Guccione's mixed branch avoids the
problem by forcing the bulk modulus to `GREAT`, which is the incompressible
limit where `J` is one and the distinction vanishes. It does not solve the
problem; it stands in the formal incompressible limit, where the problem is
not. `GREAT` is a penalty rather than an identity, so `J` is near one rather
than one, and the distinction is small rather than absent.

### 20.2 So the law does have to say something - but not choose

The thing to avoid, and the reason the solid model owns the assembly, is a law
deciding whether to use a pressure the solid model solved for. That stays. What
the law has to provide is different in kind: not a choice, but a fact about how
it is written.

A law declares whether its energy is split. A solid model running a mixed
formulation asks the manager for a deviatoric stress, exactly as now, and a law
that has not declared the split is refused - by name, at construction, before
anything runs. It cannot opt out of the pressure, and it cannot silently be
used where `dev()` means something other than what the formulation needs.

A boolean is not the whole of it, though, and this is the part that needs
thought rather than a switch. The pressure equation carries its own
compressibility - `rKappa()` in the solvers - and that has to correspond to the
law's `U(J)`. A law can be perfectly split and still be paired with a pressure
equation assuming a different volumetric response, and nothing in a boolean
would notice. So what the law has to declare is not only *that* it is split but
*what its volumetric response is*, in a form the pressure equation can be built
from.

The mechanism for the request already exists: `planeStress` and `solutionD` are
injected into every law and a law that sets either itself is refused. A
`mixedFormulation` entry belongs there. The declaration is a new virtual on the
law, defaulting to false, which is the safe direction: a law that has not
thought about it does not silently qualify.

### 20.3 What that costs

`GuccioneElastic` would refuse to run under `solvePressure` until its `Q` is
written on the isochoric strain. That is a real reformulation and a change in
what the law computes, so it is the project's call rather than this work's -
but a refusal is the right behaviour meanwhile, because the alternative is the
quiet mis-split that 19.4 describes, in the one case anybody would actually
want it for.

The eventual and better answer is that a law provides its isochoric stress
directly, evaluated on `Cbar`, *and* the volumetric response that goes with it,
rather than the solid model taking `dev()` of a full stress and building a
pressure equation from a compressibility that may not match. Then nothing is
inferred and nothing has to be declared, because the interface says it. That is
more work and it is where this should end up.
