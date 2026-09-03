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

Two things this does not do.

**The shipped `rodAndSeabed` cannot run on the framework yet.** Its
effective-stress law is `anisotropicBiotElastic`, which has no framework
counterpart, and migrating it needs one new thing: its two-dimensional
behaviour comes from `mesh.solutionD()`, and a framework law is constructed
from a dictionary with no mesh. The manager already injects `planeStress` into
each law's dictionary for the same reason, so the mechanism exists; it just has
to be extended. The tutorial's comparison therefore substitutes `linearElastic`
in both arms, which puts the composite itself under test - the effective stress
it carries, the seeding, the pressure gather - and says so where it does it.

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
