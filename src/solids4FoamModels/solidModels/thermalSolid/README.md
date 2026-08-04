---
sort: 13
---

# thermalSolid

This page documents the heat-conduction-only solid model. The runtime type is:

```text
thermalSolid
```

The model solves the transient heat equation for temperature `T` in a solid
region. It does **not** solve a momentum equation, and the displacement field
stays at its initial value. It exists to provide the solid side of conjugate
heat transfer, where a fluid model exchanges temperature and heat flux with a
rigid solid.

---

## User Guide

### What this model solves

`thermalSolid`:

- solves `rhoC*ddt(T) == laplacian(k, T)`, with no source terms;
- takes `C` and `k` from a `thermalModel`, configured in
  `constant/thermalProperties`;
- implements the `setTemperatureAndHeatFlux()`,
  `setEqInterHeatTransferCoeff()`, `faceZoneTemperature()`,
  `faceZoneHeatFlux()` and `faceZoneHeatTransferCoeff()` interface used by the
  thermal fluid-solid interfaces;
- returns a zero surface-normal gradient from `tractionBoundarySnGrad()`, so
  traction boundary conditions have no effect.

```warning
No mechanical equation is solved. If you need the deformation as well as the
temperature, use
[thermalLinGeomSolid](https://www.solids4foam.com/documentation/solid-models/thermalLinGeomSolid.html)
for a strongly coupled solve, or
[weakThermalLinGeomSolid](https://www.solids4foam.com/documentation/solid-models/weakThermalLinGeomSolid.html)
for a one-way coupled one.
```

### Supported solution algorithms

This model implements a single, segregated implicit loop for the heat
equation. It does not read `solutionAlgorithm`.

### Model options

| Entry | Default | Description |
| --- | --- | --- |
| `absoluteTemperatureTolerance` | `1e-06` | Convergence threshold on `T` |

`absoluteTemperatureTolerance` is an absolute change in degrees: once the
temperature stops moving by more than this between outer iterations, the
thermal loop is considered converged.

The relevant inherited `solidModel` entries are:

| Entry | Default | Relevance |
| --- | --- | --- |
| `nCorrectors` | `10000` | Maximum number of correctors |
| `solutionTolerance` | `1e-06` | Primary convergence tolerance |
| `alternativeTolerance` | `1e-07` | Secondary convergence tolerance |
| `infoFrequency` | `100` | Frequency for solver progress output |

The mechanical options inherited from `solidModel` are read but have no
effect, since no momentum equation is assembled.

### Required input files

- `T` and `D` in the time directory;
- `constant/solidProperties`;
- `constant/thermalProperties`, providing the `thermal` sub-dictionary;
- `constant/mechanicalProperties`, which is still needed because the density
  is taken from the mechanical law when forming `rhoC`;
- `constant/g`.

`T` is `MUST_READ`, so a case that omits it will fail at construction.

Example `constant/thermalProperties`:

```text
thermal
{
    type            constant;
    C               C [0 2 -2 -1 0 0 0] 434;
    k               k [1 1 -3 -1 0 0 0] 250;
}
```

### Recommended dictionary setup

Minimal example for `constant/solidProperties`, as used on the solid side of a
conjugate heat transfer case:

```text
solidModel        thermalSolid;

thermalSolidCoeffs
{
    nCorrectors                  10000;
    solutionTolerance            1e-06;
    alternativeTolerance         1e-07;
    absoluteTemperatureTolerance 1e-06;
    infoFrequency                100;
}
```

`fvSchemes` must provide a `laplacian(k,T)` entry, and `fvSolution` must
provide a solver for `T`.

### Boundary conditions

Temperature boundaries use the usual OpenFOAM scalar patch types, plus:

- `thermalRobin`, which carries both a temperature and a heat flux and is the
  type required on a thermal fluid-solid interface. Calling
  `setTemperatureAndHeatFlux()` on any other patch type is a fatal error.
- `thermalConvection`, for a convective heat-transfer boundary.

Displacement boundary conditions may be present but are not solved for.

### Field glossary

- `T`: temperature; the only field this model solves for.
- `grad(T)`: temperature gradient, used to report the heat flux.
- `heatFlux`: `-k*grad(T)`, constructed and written by `writeFields()`.
- `rhoC`: product of density and specific heat capacity.

Each write also prints the maximum and minimum temperature and the maximum
magnitude of the heat flux.

---

## Developer Notes

### Class role

`thermalSolid` inherits directly from `solidModel` and adds a `thermalModel`
member. It is the minimal solid model: it satisfies the `solidModel` interface
but leaves the mechanical part inert.

- `nonLinGeom()` returns `LINEAR_GEOMETRY`;
- `solutionD()` returns `D()`, but `D` is never solved for;
- `tractionBoundarySnGrad()` returns a zero field of patch size.

The header's `Description` block, inherited from `thermalLinGeomSolid`, still
mentions a strongly coupled momentum equation; the implementation solves the
heat equation only.

### Construction

The constructor calls the base `solidModel` constructor, constructs the
`thermalModel` from the mesh, builds `rhoC_` as
`thermal.C()*mechanical().rho()` and `k_` from `thermal.k()`, reads `T`
(`MUST_READ`), creates `grad(T)`, reads `absoluteTemperatureTolerance`, and
forces creation of `T.oldTime()`.

Note that `DisRequired()` is not called, consistent with the absence of a
momentum solve.

### `evolve()`

A single loop on the heat equation: store `T`, relax and solve
`rhoC*ddt(T) - laplacian(k, T)`, relax `T`, update `grad(T)`. On OpenFOAM.com
builds the solver-performance dictionary is cleared each iteration, to avoid an
expensive copy of the residual history when many correctors run. The loop
exits through the private `converged()` overload or when `iCorr` reaches
`nCorrectors`. `evolve()` always returns `true`.

### Conjugate heat transfer interface

`setTemperatureAndHeatFlux()` and `setEqInterHeatTransferCoeff()` write into a
`thermalRobin` patch field; `faceZoneTemperature()`, `faceZoneHeatFlux()` and
`faceZoneHeatTransferCoeff()` read the corresponding quantities back out on a
global face zone. These are the hooks the thermal fluid-solid interfaces call
each coupling iteration.

### Extension points

- `evolve()` if the heat equation needs source terms;
- the private `converged()` overload for a different convergence criterion.

If a mechanical solve is ever wanted here, prefer `thermalLinGeomSolid`
instead of extending this class.

---

## Tutorials

Cases that select `thermalSolid`:

- `thermoFluidSolidInteraction/flowOverHeatedPlate`
- `thermoFluidSolidInteraction/thermalCavity`

In both, `thermalSolid` is the solid region of a conjugate heat transfer
simulation, and the solid does not deform.
