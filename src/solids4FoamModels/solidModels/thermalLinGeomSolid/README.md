---
sort: 11
---

# thermalLinGeomSolid

This page documents the coupled thermal linear-geometry solid model. The
runtime type is:

```text
thermalLinearGeometry
```

The model assumes small strains and small rotations, and solves the energy
equation for temperature `T` and the linear momentum equation for total
displacement `D` inside one outer loop, so that the temperature-displacement
coupling is converged within each time step.

The stress is computed by the run-time selectable mechanical law.

---

## User Guide

### What this model solves

`thermalLinGeomSolid`:

- solves the transient heat equation
  `rhoC*ddt(T) == laplacian(k, T) + (sigma && grad(U))`, i.e. including the
  mechanical work term;
- solves the linear-geometry momentum equation for `D`;
- alternates the two equations inside a single outer loop, and requires both
  to converge before the time step ends;
- takes `C` and `k` from a `thermalModel`, configured in
  `constant/thermalProperties`;
- exposes the interface that thermal fluid-solid interaction uses to set
  patch temperature and heat flux.

```warning
The momentum equation does not include a thermal-expansion term by itself.
To get thermal stresses, select a mechanical law that consumes the temperature
field, such as `thermoMechanicalLaw`. Without one, the temperature solution
has no effect on the displacement.
```

### Supported solution algorithms

This model implements a single, segregated implicit algorithm. It does not
read `solutionAlgorithm`, and there is no PETSc SNES or explicit path.

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
| `nCorrectors` | `10000` | Maximum number of outer correctors |
| `solutionTolerance` | `1e-06` | Primary convergence tolerance, `T` and `D` |
| `alternativeTolerance` | `1e-07` | Secondary convergence tolerance |
| `materialTolerance` | `1e-05` | Mechanical-law convergence tolerance |
| `relaxationMethod` | `fixed` | Under-relaxation method (`fixed`, `aitken`) |
| `infoFrequency` | `100` | Frequency for solver progress output |
| `stabilisation` | auto-created if absent | `momentum` sub-dictionary is used |
| `cellDisplacements` | optional | Internal-cell displacement constraints |

Convergence requires the temperature _and_ the displacement residuals to be
satisfied, and at least two outer iterations are always performed.

### Required input files

- `T` and `D` in the time directory;
- `constant/solidProperties`;
- `constant/thermalProperties`, providing the `thermal` sub-dictionary;
- `constant/mechanicalProperties`;
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

Minimal example for `constant/solidProperties`:

```text
solidModel        thermalLinearGeometry;

thermalLinearGeometryCoeffs
{
    nCorrectors                  10000;
    solutionTolerance            1e-06;
    alternativeTolerance         1e-07;
    materialTolerance            1e-05;
    absoluteTemperatureTolerance 1e-06;
    infoFrequency                100;
}
```

`fvSchemes` must provide `laplacian(k,T)` and `laplacian(DD,D)` entries, and
`fvSolution` must provide solvers for `T` and `D`.

### Boundary conditions

Displacement boundaries use the standard solids4foam patch types. Temperature
boundaries use the usual OpenFOAM scalar patch types, plus:

- `thermalRobin`, which carries both a temperature and a heat flux and is the
  type required on a thermal fluid-solid interface. Calling
  `setTemperatureAndHeatFlux()` on any other patch type is a fatal error.
- `thermalConvection`, for a convective heat-transfer boundary.

### Field glossary

- `T`: temperature; primary thermal unknown, written each output time.
- `grad(T)`: temperature gradient, used to report the heat flux.
- `heatFlux`: `-k*grad(T)`, constructed and written by `writeFields()`.
- `D`, `DD`: total and incremental displacement.
- `pointD`, `pointDD`: displacement at mesh points.
- `grad(D)`, `grad(DD)`: displacement gradients.
- `sigma`: symmetric stress tensor.
- `U`: velocity field, computed as `ddt(D)`.
- `rhoC`: product of density and specific heat capacity.

Each write also prints the maximum and minimum temperature and the maximum
magnitude of the heat flux.

---

## Developer Notes

### Class role

`thermalLinGeomSolid` inherits directly from `solidModel` and adds a
`thermalModel` member. The key design choices are:

- `D` is the primary solution variable (`solutionD()` returns `D()`);
- `nonLinGeom()` returns `LINEAR_GEOMETRY`;
- the momentum equation is the same deferred-correction Laplacian form used by
  `linGeomTotalDispSolid`, with the momentum stabilisation term applied
  explicitly;
- the heat equation and the momentum equation share one outer loop, which is
  what distinguishes this model from
  [weakThermalLinGeomSolid](https://www.solids4foam.com/documentation/solid-models/weakThermalLinGeomSolid.html).

### Construction

The constructor calls the base `solidModel` constructor and `DisRequired()`,
constructs the `thermalModel` from the mesh, builds `rhoC_` as
`thermal.C()*mechanical().rho()` and `k_` from `thermal.k()`, reads `T`
(`MUST_READ`) and creates `grad(T)`, reads `absoluteTemperatureTolerance`, and
takes `impK_`, `impKf_` and `rImpK_` from the mechanical law. It then forces
creation of `T.oldTime()`.

### `evolve()`

Each outer iteration:

1. stores `T` for under-relaxation, assembles and solves the heat equation
   including the `sigma && grad(U)` work term, relaxes `T`, and updates
   `grad(T)`;
2. stores `D`, updates the momentum stabilisation, assembles

   ```text
   rho*d2dt2(D) == laplacian(impKf, D) - fvc::laplacian(impKf, D)
                 + fvc::div(sigma) + rho*g + impK*stabilisation
   ```

   then relaxes and solves it, and under-relaxes `D`;
3. updates `DD`, `U`, `grad(D)`, `grad(DD)` and `sigma`;
4. refreshes `impKf_` from the mechanical law every tenth iteration, to help
   convergence when the material stiffness changes. `impK_` and `rImpK_` are
   deliberately _not_ refreshed, because they are used by the traction
   boundary conditions.

The loop exits through the private `converged()` overload, which checks both
`T` and `D`, or when `iCorr` reaches `nCorrectors`. Afterwards `pointD`, `DD`
and `pointDD` are updated once.

### Thermal fluid-solid interaction interface

`setTemperatureAndHeatFlux()` and `setEqInterHeatTransferCoeff()` write into a
`thermalRobin` patch field, and `faceZoneTemperature()`,
`faceZoneHeatFlux()` and `faceZoneHeatTransferCoeff()` read the corresponding
quantities back out on a global face zone. These are the hooks used by the
thermal fluid-solid interfaces.

### Extension points

- `evolve()` for a different coupling strategy between the two equations;
- the private `converged()` overload if a different convergence criterion is
  needed;
- `writeFields()` for extra derived thermal output.

---

## Tutorials

Cases that select `thermalLinearGeometry`:

- `solids/thermoelasticity/hotCylinder`
- `solids/thermoelasticity/hotSphere`
- `solids/thermoelasticity/slabCooling`
- `thermoFluidSolidInteraction/hotTJunction`
