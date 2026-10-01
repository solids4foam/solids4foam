---
sort: 27
---

# electroMechanicalLaw

A composite law that adds an active tension along a fibre direction to the
stress of a passive law. The runtime type is:

```text
electroMechanicalLaw
```

---

## User Guide

### What it computes

This is a **finite-strain** law: it implements the finite-strain update only,
so both it and its passive law are used with a nonlinear geometry solid
model.

The passive law, given in the `passiveMechanicalLaw` sub-dictionary, is
evaluated first. The active stress is then added:

```text
f0f0  = (f0/|f0|)*(f0/|f0|)
sigma = sigma_passive + symm(F & (Ta*f0f0) & F^T)/J
```

The active tension `Ta` is either a field or a ramped constant:

- with `activeTensionFromField yes`, `Ta` is the coupling input named by
  `activeTensionFieldName` (see [Coupling inputs](#coupling-inputs)), and
  `rampTime` is not applied;
- otherwise `Ta` is `activeTension`, ramped linearly from zero:

```text
Ta = (t/rampTime)*activeTension    for t < rampTime
Ta = activeTension                 otherwise, or when rampTime = 0
```

The density and the bulk modulus are those of the passive law. `rho` given
on this law is passed down to the passive law when its own dictionary has
none, so it can be given once, on the outer dictionary.

### Model options

| Entry | Required | Default | Description |
| --- | --- | --- | --- |
| `passiveMechanicalLaw` | yes | - | Sub-dictionary selecting the passive law |
| `activeTension` | yes | - | Constant active tension, pressure units |
| `rampTime` | yes | - | Ramp duration for `activeTension`, `>= 0` |
| `activeTensionFromField` | no | `no` | Read `Ta` as a coupling input field |
| `activeTensionFieldName` | no | `Ta` | Name of the active tension field |
| `uniformFibreField` | no | `no` | Take the default `f0` from this dictionary |
| `f0` | if uniform | - | Fibre vector; required if `uniformFibreField yes` |
| `rho` | see text | - | Density, passed to the passive law if absent there |
| `inputCaseDirectories` | no | - | Read `Ta` from another case directory |

`activeTension` and `rampTime` are required even when the tension comes from
a field. A negative `rampTime` is a fatal error.

The `passiveMechanicalLaw` sub-dictionary holds the passive law's `type` and
all of its entries. Any finite-strain law can be used; the tutorials use
[`GuccioneElastic`](../GuccioneElasticMechanicalConstitutiveLaw/README.md).

### Coupling inputs

With `activeTensionFromField yes` the law declares the active tension field
(`Ta` by default) as a required scalar input. The manager gathers it on every
evaluation, at cells, faces and boundary faces alike, and takes it from the
first of:

1. a `volScalarField` of that name registered by another model, e.g. an
   electrophysiology model coupled to the solid;
2. the case directory named for it in `inputCaseDirectories` (or in the
   older `TacaseDirectory` spelling), read at the current time; see
   [the framework README](../../README.md#scalar-inputs-from-another-case-directory);
3. a file of that name in the current time directory, or failing that in `0`.

If none exists the run stops with an error.

With `activeTensionFromField no` (the default) the field is not read, even
if one is registered. In that case the manager prints a warning, once, if a
`volScalarField` with the `activeTensionFieldName` name is registered, since
the case is then running on the constant rather than on that field. An
`inputCaseDirectories` entry for `Ta` while `activeTensionFromField` is `no`
is a fatal error, because the law does not read that input; a
`TacaseDirectory` entry is then ignored.

### Fibre direction

`f0` is a prescribed state variable. It is looked up as for `GuccioneElastic`:
a `volVectorField` `f0` in the current time directory, then in `0`, then a
registered field, and only if none exists the dictionary `f0` (with
`uniformFibreField yes`) or a zero vector. A zero-length fibre stops the run
at the first evaluation. `f0` is normalised before use.

This law and its passive law each hold their own `f0`. An `f0` field on disk
is read by both. `uniformFibreField` and `f0` are not passed down, however:
a passive `GuccioneElastic` with a uniform fibre needs them in its own
`passiveMechanicalLaw` dictionary as well as here.

### State variables

| Name | Type | Role | Description |
| --- | --- | --- | --- |
| `f0` | vector | prescribed | Reference fibre direction of the active stress |

The passive law's state lives in a child state named `passiveMechanicalLaw`.
Prescribed variables are neither written with the results nor restarted
from; a passive law with persistent history has it written as described in
that law's documentation, with restart files named
`<material>_<topology>_passiveMechanicalLaw_<variable>`.

### Tangents

- `scalar` and `scalarDeviatoric`: those of the passive law. The active
  tension is not included.
- `fourthOrderFiniteDifference`: recomputed by this law over the whole
  stress, active part included.
- `fourthOrder` (analytical): refused, since the passive law's tangent would
  miss the active term.

### Mixed displacement-pressure formulation

The law provides the volumetric split exactly when its passive law does. The
volumetric response is the passive law's `dU/dJ`; the active stress, whose
mean is not a function of `J`, stays in the isochoric part that the solved
pressure does not replace. So with a passive `GuccioneElastic` the law can be
used with a mixed formulation, as in `LandEtAl2015/problem3` run with
`./Allrun pressure`.

### Differences from the legacy law

- A registered `Ta` field is no longer picked up automatically. Set
  `activeTensionFromField yes` to use one; without it the constant is used
  and a warning is printed if a `Ta` field is registered. The field name can
  be changed with `activeTensionFieldName`.
- `Ta` can now also be read from a file in the time directory, or from
  another case directory with `inputCaseDirectories`.
- `f0` may be given in this dictionary with `uniformFibreField`, is
  normalised, and no longer needs a face field `f0f`.
- The legacy face-centred path multiplied the active term by `J` rather
  than dividing by it; every evaluation path now uses the relation above.
- The fourth-order finite-difference tangent includes the active stress; the
  legacy law returned the passive law's tangent.
- `regionName` is no longer read. `solvePressureEqn` is read by the solid
  model wherever it appears in the material's dictionary, including inside
  `passiveMechanicalLaw`, and applies to the whole mesh; see the
  [material models page](../../../materialModels/README.md).

### Example

From the `LandEtAl2015/problem3` tutorial (comments and unused entries
removed), where `f0` is written into `0` by `setFibreField`:

```text
planeStress     no;

mechanical
(
    heartTissue
    {
        type            electroMechanicalLaw;
        rho             rho [1 -3 0 0 0 0 0] 3000;
        activeTension   activeTension [1 -1 -2 0 0 0 0] 60.0e3;
        rampTime        1.0;

        passiveMechanicalLaw
        {
            type            GuccioneElastic;
            k               k [1 -1 -2 0 0 0 0] 2e3;
            cf              8.0;
            ct              2.0;
            cfs             4.0;
            bulkModulus     bulkModulus [1 -1 -2 0 0 0 0] 16e6;
        }
    }
);
```

To drive the tension from a field supplied by another model, or from files
`Ta` in the time directories, add:

```text
        activeTensionFromField  yes;
        // activeTensionFieldName  Ta;
```

---

## Tutorials

- [`solids/hyperelasticity/LandEtAl2015/problem3`](../../../../../tutorials/solids/hyperelasticity/LandEtAl2015/problem3)
