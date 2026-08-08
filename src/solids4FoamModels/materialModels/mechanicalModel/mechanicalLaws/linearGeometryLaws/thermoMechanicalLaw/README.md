---
sort: 8
---

# thermoMechanicalLaw

This law adds an isotropic thermal stress to the stress calculated by a nested
linear-geometry mechanical law. The runtime type is:

```text
thermoMechanicalLaw
```

---

## User Guide

### What it computes

The cell-centred stress is first calculated by the nested mechanical law and
then corrected using its bulk modulus:

```text
sigma = sigmaMech - 3*K*alpha*(T - T0)*I
```

For the face-centred overload, the temperature and nested-law bulk modulus are
interpolated from cells to faces:

```text
sigmaf = sigmaMech,f - 3*interpolate(K)*alpha*(interpolate(T) - T0)*I
```

The law implements `correct(volSymmTensorField&)` and
`correct(surfaceSymmTensorField&)`. It does not implement the
`pointSymmTensorField` overload, which aborts with `notImplemented`. On
OpenFOAM.com and OpenFOAM.org, the `CompactListList` overload is present in the
base-class interface but also aborts with `notImplemented`. That interface is
excluded from foam-extend by `#ifndef FOAMEXTEND`.

The nested law is selected from the linear-geometry runtime table. The
decoupled thermal correction may not be appropriate when the nested law is a
nonlinear function of volumetric stress, such as a plasticity law whose
response depends on volumetric stress.

### Model options

| Entry | Required | Description |
| --- | --- | --- |
| `rho` | yes | Density, `[1 -3 0 0 0 0 0]` |
| `alpha` | yes | Linear expansion coefficient, `[0 0 0 -1 0 0 0]` |
| `T0` | yes | Stress-free reference temperature, `[0 0 0 1 0 0 0]` |
| `mechanicalLaw` | yes | Nested linear-geometry law dictionary |
| `TcaseDirectory` | no | Relative temperature-case path; default `.` |
| `solvePressureEqn` | no | Outer base-class switch; default `no` |
| `pressureSmoothingScaleFactor` | no | Outer smoothing scale; default `100` |
| `regionName` | conditional | Base-mesh region name; detected if omitted |

The `mechanicalLaw` sub-dictionary must contain its own `type` and all entries
required by the selected law. The outer `solvePressureEqn` and
`pressureSmoothingScaleFactor` entries are read by the base constructor, but
this wrapper does not call `updateSigmaHyd()`, so they do not change its stress
calculation. Put those entries inside `mechanicalLaw` when the nested law uses
them.

If `regionName` is absent, the base class selects a registered mesh named
`solid`, then `region0`. Construction aborts if neither is available.

The law first uses an in-memory `T` field when one is available. Otherwise it
reads `T` from the current time of `TcaseDirectory`. A separate temperature
case must use the same base mesh point coordinates and the same numbers of
points, faces and cells as the mechanical case.

### Recommended dictionary setup

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

        // Optional: read T from a separate case relative to this case
        // TcaseDirectory  "steelTemperatureCase";

        mechanicalLaw
        {
            type            linearElastic;
            E               E [1 -1 -2 0 0 0 0] 190e9;
            nu              nu [0 0 0 0 0 0 0] 0.305;

            // Optional for the nested law
            // solvePressureEqn no;
        }
    }
);
```

### Field glossary

- `T`: temperature used in the thermal correction. It is either looked up
  from the mechanical mesh registry or read from disk. A disk-read field is
  set to `AUTO_WRITE` so it is available for restart and inspection.
- `Tf`: temporary face interpolation of `T` used by the surface correction.
- `Kf`: temporary face interpolation of the nested law's bulk modulus.
- The nested `mechanicalLaw` owns the strain, stress-history and other
  constitutive fields needed to calculate `sigmaMech`.

---

## Developer Notes

### Class role

`thermoMechanicalLaw` derives directly from `mechanicalLaw` and owns an
`autoPtr<mechanicalLaw>` for the nested law. It also stores the uniform
`dimensionedScalar` values `alpha_` and `T0_`, a lazily allocated temperature
field, a disk-read flag, the relative temperature-case path, lazily allocated
`Time` and `fvMesh` objects, and the last checked time index.

The class has no `#ifdef FOAMEXTEND` or `#ifndef FOAMEXTEND` block of its own.
It uses `OPENFOAM_COM` and `OPENFOAM_NOT_EXTEND` guards for constructor and
field-API differences. The inherited quadrature-point interface is guarded by
`#ifndef FOAMEXTEND`. The law is listed in both `Make/files.openfoam` and
`Make/files.foamextend`.

### Construction

The base constructor first reads `solvePressureEqn` and
`pressureSmoothingScaleFactor`, inserting defaults of `false` and `100.0`, and
then resolves `regionName`. The wrapper then:

1. Requires the `mechanicalLaw` sub-dictionary and its `type` entry.
2. Selects and constructs that law through `NewLinGeomMechLaw`.
3. Reads `alpha` and `T0` as dimensioned scalars.
4. Inserts `TcaseDirectory` with the default `.` when it is absent.
5. Initialises the temperature, time and mesh pointers and the time index.

Missing required entries cause dictionary lookup errors. An unknown nested
type, or a type from the nonlinear-geometry selection table, causes a fatal
selection error. The selected nested law can emit its own construction errors
and warnings. The wrapper constructor performs no range checks on `alpha` or
`T0` and emits no warning itself.

Temperature-case objects are created only when disk access is needed. Reading
aborts if the temperature mesh topology or point locations differ from the
base mesh, or if no `T` field is available in memory or on disk. Disk lookup is
attempted at most once per mechanical time step. Once a disk field has been
read, its last value remains available when a later time has no new `T` file.

### Key methods

- `impK()` delegates to the nested law. It is the diffusivity used for the
  solid model's implicit Laplacian term and affects outer-iteration convergence
  rate rather than the converged answer.
- `bulkModulus()` delegates to the nested law and supplies `K` for the thermal
  stress correction.
- `correct(volSymmTensorField&)` delegates the mechanical stress calculation,
  obtains `T`, and subtracts the cell-centred thermal stress.
- `correct(surfaceSymmTensorField&)` performs the same sequence with
  interpolated temperature and bulk modulus.
- `lookupTemperatureField()` prefers a solver-registered `T` until a
  temperature has been read from disk; after that, it uses the disk-backed
  field.
- `readTField()` synchronises the temperature-case time, reads at most once per
  time step, copies a single-material field directly, and maps a base-mesh
  field onto the current material mesh for a multi-material case.
- `materialTangent()` is not overridden. The base implementation aborts with
  `notImplemented`, so `materialTangentField()` also cannot be used.

### Extension points

A related additive wrapper can copy this class, select its underlying law from
the appropriate runtime table, and replace the correction terms in both
implemented `correct()` overloads. It must continue to delegate `impK()` and
any required modulus queries. Supporting vertex-centred or higher-order solid
models additionally requires point and `CompactListList` overloads instead of
the base-class `notImplemented` fallbacks. A new wrapper must also be added to
the applicable build lists and runtime selection table.

The source is at
[thermoMechanicalLaw.C](https://github.com/solids4foam/solids4foam/blob/master/src/solids4FoamModels/materialModels/mechanicalModel/mechanicalLaws/linearGeometryLaws/thermoMechanicalLaw/thermoMechanicalLaw.C).

---

## Tutorials

- `solids/thermoelasticity/slabCooling`
- `solids/thermoelasticity/hotCylinder/hotCylinder`
- `solids/thermoelasticity/hotCylinder/hotCylinderPredefinedTFieldMultipleMaterials`
- `solids/thermoelasticity/hotSphere`
- `thermoFluidSolidInteraction/flowOverHeatedPlate`
- `thermoFluidSolidInteraction/thermalCavity`
- `thermoFluidSolidInteraction/hotTJunction`
