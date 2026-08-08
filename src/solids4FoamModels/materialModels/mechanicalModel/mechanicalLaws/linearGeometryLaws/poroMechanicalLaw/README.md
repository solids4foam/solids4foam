---
sort: 7
---

# poroMechanicalLaw

This page documents the linear-geometry wrapper that obtains effective stress
from another run-time selectable mechanical law and subtracts a scaled pore
pressure. The runtime type is:

```text
poroMechanicalLaw
```

---

## User Guide

### What it computes

The law asks the nested mechanical law to update the effective stress and then
forms the total stress as

```text
sigma = sigmaEff - b*(p + p0)*I
```

where `b` is the Biot coefficient. For a face-centred stress field, the current
pressure and initial pressure are interpolated to the faces:

```text
pf    = interpolate(p)
p0f   = interpolate(p0)
sigma = sigmaEfff - b*(pf + p0f)*I
```

On the first call, the effective-stress field is initialised by adding the pore
term to the supplied total stress. Subsequent calls retain that field and pass
it to the nested law.

The law implements `correct(volSymmTensorField&)` and
`correct(surfaceSymmTensorField&)`. The point-centred overload is inherited
from `mechanicalLaw` and aborts with `notImplemented`. The
`CompactListList` quadrature-point overload is also inherited and aborts with
`notImplemented`; that overload is compiled only when `FOAMEXTEND` is not
defined. The nested effective-stress law must itself implement whichever of the
cell- or face-centred overloads the selected solid model calls.

### Model options

| Entry | Required | Description |
| --- | --- | --- |
| `effectiveStressMechanicalLaw` | yes | Nested mechanical-law dictionary |
| `biotCoeff` | no | Biot coefficient, dimensionless; default `1.0` |
| `pressureFieldName` | no | Pressure field name; default `p` |
| `pressureFieldRegion` | no | Pressure field region; default `region0` |
| `p0` | no | Uniform initial pore pressure; default `0`, `[1 -1 -2 0 0 0 0]` |
| `rho` | yes | Density, `[1 -3 0 0 0 0 0]` |
| `regionName` | no | Base solid-region name; otherwise detected |
| `solvePressureEqn` | no | Base switch; default `no`, unused by this wrapper |
| `pressureSmoothingScaleFactor` | no | Default `100`; unused by this wrapper |

The nested dictionary must contain its own `type` and all entries required by
that law. Density belongs to the outer dictionary because `rho()` reads it from
the wrapper's dictionary.

Pressure lookup first checks `pressureFieldRegion`. If that registry does not
exist, the law falls back to the `solid` registry. It does not fall back when
the selected registry exists but does not contain `pressureFieldName`.

The outer `solvePressureEqn` and `pressureSmoothingScaleFactor` entries are read
by the base constructor, but this wrapper never calls `updateSigmaHyd()`, so
they do not alter its pore-pressure term. A nested law may read and use its own
copies of those entries.

### Recommended dictionary setup

The following setup represents a saturated soil with an isotropic linear
elastic effective-stress response:

```text
planeStress     no;

mechanical
(
    saturatedSoil
    {
        type            poroMechanicalLaw;
        rho             rho [1 -3 0 0 0 0 0] 2000;

        effectiveStressMechanicalLaw
        {
            type        linearElastic;
            E           E [1 -1 -2 0 0 0 0] 20e6;
            nu          nu [0 0 0 0 0 0 0] 0.3;
        }

        // Optional
        // biotCoeff          biotCoeff [0 0 0 0 0 0 0] 1.0;
        // pressureFieldName  p;
        // pressureFieldRegion region0;
        // p0                 p0 [1 -1 -2 0 0 0 0] 0;
    }
);
```

### Field glossary

- `sigmaEff`: cell-centred effective stress, created on the first cell-centred
  update and written automatically.
- `sigmaEfff`: face-centred effective stress, created on the first face-centred
  update and written automatically.
- `p`: externally maintained cell-centred pressure field; its actual name is
  selected by `pressureFieldName`.
- `p0`: uniform cell-centred initial pore pressure, held by the law and neither
  read nor written as a field.
- `p0f`: demand-created face interpolation of `p0`.
- `biotCoeff`: temporary uniform cell field returned by `biotCoeff()` and not
  written.

---

## Developer Notes

### Class role

`poroMechanicalLaw` derives directly from `mechanicalLaw` and is registered in
the `linGeomMechLaw` selection table. It owns the nested law through an
`autoPtr<mechanicalLaw>`, keeps demand-created cell and face effective-stress
fields, and stores the Biot coefficient, pressure field name, pressure region,
initial pressure field and a pointer to its face interpolation.

The class contains no `FOAMEXTEND` preprocessor guards. Its source appears in
both `Make/files.openfoam` and `Make/files.foamextend`. The inherited
quadrature-point interface is absent from foam-extend because that interface is
guarded by `#ifndef FOAMEXTEND` in `mechanicalLaw`.

### Construction

The `mechanicalLaw` constructor first reads `solvePressureEqn` and
`pressureSmoothingScaleFactor`, then uses `regionName` or detects `solid` or
`region0`. It emits a fatal error if no solid region can be found.

The wrapper then requires the `effectiveStressMechanicalLaw` sub-dictionary,
reads its `type`, and constructs that linear-geometry law. It next reads
`biotCoeff`, `pressureFieldName`, `pressureFieldRegion` and `p0`, in that order,
using the defaults listed above. A non-zero `p0` produces an informational
message. The wrapper itself emits no construction warning. Missing required
dictionary entries cause dictionary lookup errors, and the nested law can emit
its own validation errors and warnings.

At stress correction, lookup aborts if the pressure field cannot be found in
the configured region or in `solid`. The demand-driven `p0f` constructor also
contains a defensive fatal error if its pointer is already set.

### Key methods

- `impK()` delegates directly to the nested effective-stress law. It is the
  diffusivity of the solid model's implicit Laplacian term and affects the
  outer-iteration convergence rate rather than the converged answer.
- `correct(volSymmTensorField&)` looks up `p`, creates `sigmaEff` if needed,
  updates it with the nested law, and subtracts the cell-centred pore term.
- `correct(surfaceSymmTensorField&)` interpolates `p`, creates `sigmaEfff` if
  needed, updates it with the nested law, and subtracts the face-centred pore
  term.
- `biotCoeff()` returns a uniform cell-centred field containing `b`.
- `materialTangent()` is not overridden or delegated. The inherited
  implementation aborts with `notImplemented`.

### Extension points

A related linear-geometry wrapper can copy this class, replace the added stress
contribution, and retain delegation of `impK()` and stress correction to a
nested law. Any additional dictionary entries should be read in the constructor.
Support for point or quadrature-point solid models requires explicit overloads
that also delegate to compatible nested-law interfaces. The new class must be
registered in the appropriate mechanical-law selection table and added to the
OpenFOAM and foam-extend build lists where supported.

The source is at
[poroMechanicalLaw.C](https://github.com/solids4foam/solids4foam/blob/master/src/solids4FoamModels/materialModels/mechanicalModel/mechanicalLaws/linearGeometryLaws/poroMechanicalLaw/poroMechanicalLaw.C).

---

## Tutorials

- `solids/poroelasticity/stripFooting`
- `solids/poroelasticity/suctionCaission`
- `solids/poroelasticity/rodAndSeabed`
