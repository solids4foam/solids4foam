---
sort: 8
---

# thermoMechanicalLaw

A composite law that adds an isotropic thermal stress to the stress of a
sub-law. The runtime type is:

```text
thermoMechanicalLaw
```

---

## User Guide

### What it computes

This is a **small-strain** law: it implements the small-strain update only,
so it and its sub-law are used with a linear geometry solid model. A
nonlinear geometry solid model asks for the finite-strain update, and the run
stops with an error.

The sub-law, given in the `mechanicalLaw` sub-dictionary, is evaluated first,
and the thermal stress is then subtracted:

```text
sigma = sigma_mechanical - 3*K*alpha*(T - T0)*I
```

where `K` is the sub-law's bulk modulus, `alpha` the linear expansion
coefficient, `T` the temperature and `T0` the stress-free reference
temperature. `K` is taken as a single value for the material, so the sub-law
must provide a bulk modulus; for a plastic sub-law it is the elastic one.

The density and the bulk modulus are those of the sub-law. `rho` given on
this law is passed down to the sub-law when its own dictionary has none, so
it can be given once, on the outer dictionary.

### Model options

| Entry | Required | Default | Description |
| --- | --- | --- | --- |
| `mechanicalLaw` | yes | - | Sub-dictionary selecting the sub-law |
| `alpha` | yes | - | Expansion coefficient, `[0 0 0 -1 0 0 0]` |
| `T0` | yes | - | Stress-free reference temperature, `[0 0 0 1 0 0 0]` |
| `TName` | no | `T` | Name of the temperature field |
| `rho` | see text | - | Density, passed to the sub-law if absent there |
| `TcaseDirectory` | no | - | Read `T` from another case directory |
| `inputCaseDirectories` | no | - | As `TcaseDirectory`, in a sub-dictionary |

`alpha` and `T0` are dimension-checked: any other dimensions are a fatal
error. The `mechanicalLaw` sub-dictionary holds the sub-law's `type` and all
of its entries. Any small-strain law can be used; the tutorials use
`linearElastic`.

### Coupling inputs

The temperature is a required scalar input, named by `TName`. The manager
gathers it on every evaluation, at cells, faces and boundary faces alike, and
takes it from the first of:

1. a `volScalarField` of that name registered by another model, such as the
   temperature that
   [`thermalLinearGeometry`](../../../solidModels/thermalLinGeomSolid/README.md)
   solves for;
2. the case directory given by `TcaseDirectory`, or for the input in
   `inputCaseDirectories`, read at the current time;
3. a file of that name in the current time directory, or failing that in `0`.

If none exists the run stops with an error.

The second is for a temperature computed beforehand. The path is relative to
the case (in parallel, to each `processorN` directory); the other case must
have the same mesh and patches, and its field is read once per time step at
the current time, keeping the last one read when a time directory has none.
The copy is registered under the input's name and written with the
results. Giving both
spellings for the same input is an error. See
[the framework README](../../README.md#scalar-inputs-from-another-case-directory)
for details.

### State variables

This law declares no state of its own. The sub-law's state lives in a child
state named `mechanicalLaw`; a sub-law with persistent history has it written
as described in that law's documentation, with restart files named
`<material>_<topology>_mechanicalLaw_<variable>`.

### Tangents

The thermal stress does not depend on the strain, so the tangents are the
sub-law's, unchanged: whichever of `scalar`, `scalarDeviatoric`,
`fourthOrder` and `fourthOrderFiniteDifference` the sub-law supports.

### Mixed displacement-pressure formulation

The law does not declare a volumetric split, whatever the sub-law: it
returns the total stress only. A mixed
displacement-pressure formulation therefore refuses to start with it, and so
does hydrostatic stress smoothing (`solvePressureEqn`), wherever in the
material's dictionary it is set.

### Differences from the legacy law

- The temperature can also come from a file in the current or `0` time
  directory; `TcaseDirectory` no longer defaults to `.`, and an
  `inputCaseDirectories` sub-dictionary may be used instead.
- `TName` is new, for a temperature field with another name.
- `alpha` and `T0` are dimension-checked.
- `solvePressureEqn` and `pressureSmoothingScaleFactor` inside
  `mechanicalLaw` no longer smooth the sub-law's hydrostatic stress: the
  request is found by the solid model and refused, since this law provides
  no volumetric split. `regionName` is no longer read.
- The sub-law's tangent is passed through, so solid models that ask for a
  fourth-order tangent can be used; the legacy law had none.

### Example

From the `hotSphere` tutorial, with
[`thermalLinearGeometry`](../../../solidModels/thermalLinGeomSolid/README.md):

```text
planeStress     no;

mechanical
(
    steel
    {
        type            thermoMechanicalLaw;
        rho             rho [1 -3 0 0 0 0 0] 7750;
        alpha           alpha [0 0 0 -1 0 0 0] 9.7e-06;
        T0              T0 [0 0 0 1 0 0 0] 300;

        mechanicalLaw
        {
            type            linearElastic;
            E               E [1 -1 -2 0 0 0 0] 190e+9;
            nu              nu [0 0 0 0 0 0 0] 0.305;
        }
    }
);
```

To read the temperature from a separate case instead, as
`hotCylinderPredefinedTFieldMultipleMaterials` does with
`linearGeometryTotalDisplacement`:

```text
        TcaseDirectory  "hotCylinderTemperatureField";
```

---

## Tutorials

- [`solids/thermoelasticity/hotSphere`](../../../../../tutorials/solids/thermoelasticity/hotSphere)
- [`solids/thermoelasticity/slabCooling`](../../../../../tutorials/solids/thermoelasticity/slabCooling)
- [`solids/thermoelasticity/hotCylinder/hotCylinder`](../../../../../tutorials/solids/thermoelasticity/hotCylinder/hotCylinder)
- [`solids/thermoelasticity/hotCylinder/hotCylinderPredefinedTFieldMultipleMaterials`](../../../../../tutorials/solids/thermoelasticity/hotCylinder/hotCylinderPredefinedTFieldMultipleMaterials):
  two materials, temperature read from another case directory
- [`thermoFluidSolidInteraction/flowOverHeatedPlate`](../../../../../tutorials/thermoFluidSolidInteraction/flowOverHeatedPlate)
- [`thermoFluidSolidInteraction/hotTJunction`](../../../../../tutorials/thermoFluidSolidInteraction/hotTJunction)
- [`thermoFluidSolidInteraction/thermalCavity`](../../../../../tutorials/thermoFluidSolidInteraction/thermalCavity)
