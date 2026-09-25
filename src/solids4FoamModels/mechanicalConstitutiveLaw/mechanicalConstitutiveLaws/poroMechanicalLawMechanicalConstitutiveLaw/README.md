---
sort: 7
---

# poroMechanicalLaw

A composite law giving the total stress of a porous medium from the effective
stress of a sub-law and a pore pressure. The runtime type is:

```text
poroMechanicalLaw
```

---

## User Guide

### What it computes

This is a **small-strain** law: it implements the small-strain update only,
so it and its sub-law are used with a linear geometry solid model, typically
[`poroLinearGeometry`](../../../solidModels/poroLinGeomSolid/README.md). A
nonlinear geometry solid model asks for the finite-strain update, and the run
stops with an error.

The sub-law, given in the `effectiveStressMechanicalLaw` sub-dictionary,
computes the effective stress, and the law returns

```text
sigma = sigmaEff - b*(p + p0)*I
```

where `b` is the Biot coefficient, `p` the pore pressure and `p0` a uniform
reference pressure.

The effective stress is kept as persistent state. At the first evaluation of
each integration point it is seeded from the total stress already standing
there:

```text
sigmaEff = sigma + b*(p + p0)*I
```

and that seeded value is handed to the sub-law as its incoming stress. After
each evaluation `sigmaEff` holds what the sub-law returned. None of the
current laws reads its incoming stress, so for them the seed does not change
the answer; it is there for a sub-law whose response depends on the stress
already standing.

The density and the bulk modulus are those of the sub-law. `rho` given on
this law is passed down to the sub-law when its own dictionary has none, so
it can be given once, on the outer dictionary.

### Model options

| Entry | Required | Default | Description |
| --- | --- | --- | --- |
| `effectiveStressMechanicalLaw` | yes | - | Sub-dictionary for the sub-law |
| `biotCoeff` | no | `1` | Biot coefficient `b`, dimensioned, dimensionless |
| `p0` | no | `0` | Reference pore pressure, `[1 -1 -2 0 0 0 0]` |
| `pressureFieldName` | no | `porePressure` | Name of the pore pressure field |
| `rho` | see text | - | Density, passed to the sub-law if absent there |
| `inputCaseDirectories` | no | - | Read the pore pressure from another case |

`biotCoeff` and `p0` are dimensioned scalars and are dimension-checked:
`biotCoeff` must be dimensionless and `p0` must have pressure dimensions.
The `effectiveStressMechanicalLaw` sub-dictionary holds the sub-law's `type`
and all of its entries. Any small-strain law can be used; the tutorials use
`linearElasticMohrCoulombPlastic` and `anisotropicBiotElastic`.

### Coupling inputs

The pore pressure is a required scalar input, named by `pressureFieldName`.
The manager gathers it on every evaluation, at cells, faces and boundary
faces alike, and takes it from the first of:

1. a `volScalarField` of that name registered by another model, such as the
   `porePressure` that `poroLinearGeometry` solves for;
2. the case directory given for it in `inputCaseDirectories`, or as
   `<name>caseDirectory` (e.g. `porePressurecaseDirectory`), read at the
   current time; see
   [the framework README](../../README.md#scalar-inputs-from-another-case-directory);
3. a file of that name in the current time directory, or failing that in `0`.

If none exists the run stops with an error. The field is looked up on the
solid mesh only.

### State variables

| Name | Type | Role | Description |
| --- | --- | --- | --- |
| `sigmaEff` | symmTensor | persistent | Effective stress |
| `sigmaEffSeeded` | scalar | persistent | `1` once `sigmaEff` has been seeded |

The sub-law's state lives in a child state named
`effectiveStressMechanicalLaw`.

Persistent state is written at each write time in two forms. For viewing,
each variable is written as a `volField` named after the variable
(`sigmaEff`, `sigmaEffSeeded`, and the sub-law's own, such as `activeYield`
for `linearElasticMohrCoulombPlastic`), covering every material that declares
it and zero elsewhere, and without dimensions; this is done for the
cell-centred topology only, and a name another model has registered is
skipped. For restart, the state is
written as it is held, in files named `<material>_<topology>_<variable>`,
e.g. `soil_cellCentredIntegrationPointTopology_sigmaEff`, with the sub-law's
variables as
`<material>_<topology>_effectiveStressMechanicalLaw_<variable>`. Only the
restart files are read back.

### Tangents

The pore pressure term does not depend on the strain, so the tangents are the
sub-law's, unchanged: whichever of `scalar`, `scalarDeviatoric`,
`fourthOrder` and `fourthOrderFiniteDifference` the sub-law supports.

### Mixed displacement-pressure formulation

The law does not declare a volumetric split, whatever the sub-law: it returns
the total stress only. A mixed displacement-pressure formulation therefore
refuses to start with it, and so does hydrostatic stress smoothing
(`solvePressureEqn`), wherever in the material's dictionary it is set.

### Differences from the legacy law

- `pressureFieldRegion` is no longer read: the pore pressure is taken from
  the solid mesh, and can also come from a file in the time directory or
  from another case directory. Unused legacy entries are ignored rather than
  reported.
- `p0` is a reference pressure added to `p`, as before, and still has no
  field form.
- `biotCoeff` and `p0` are dimension-checked.
- The effective stress is persistent state, seeded per integration point and
  restarted with the rest of the state; the legacy law held separate
  `sigmaEff` and `sigmaEfff` fields.
- `solvePressureEqn` and `pressureSmoothingScaleFactor` inside the sub-law
  dictionary no longer smooth the sub-law's hydrostatic stress: the request
  is found by the solid model and refused. `regionName` is no longer read.
- The sub-law's tangent is passed through; the legacy law had none.

### Example

From the `stripFooting` tutorial, with `poroLinearGeometry`:

```text
planeStress     no;

mechanical
(
    soil
    {
        type            poroMechanicalLaw;
        rho             rho [1 -3 0 0 0 0 0] 7854;
        biotCoeff       biotCoeff [0 0 0 0 0 0 0] 1.0;

        effectiveStressMechanicalLaw
        {
            type            linearElasticMohrCoulombPlastic;
            E               E [ 1 -1 -2 0 0 0 0 ] 20e6;
            nu              nu [0 0 0 0 0 0 0] 0.3;
            frictionAngle   frictionAngle [0 0 0 0 0 0 0] 30;
            dilationAngle   dilationAngle [0 0 0 0 0 0 0] 0;
            cohesion        cohesion [1 -1 -2 0 0 0 0] 1e5;
        }
    }
);
```

Optional entries, with their defaults:

```text
        p0                  p0 [1 -1 -2 0 0 0 0] 0;
        pressureFieldName   porePressure;
```

---

## Tutorials

- [`solids/poroelasticity/stripFooting`](../../../../../tutorials/solids/poroelasticity/stripFooting)
- [`solids/poroelasticity/suctionCaission`](../../../../../tutorials/solids/poroelasticity/suctionCaission)
- [`solids/poroelasticity/rodAndSeabed`](../../../../../tutorials/solids/poroelasticity/rodAndSeabed):
  with an `anisotropicBiotElastic` sub-law
