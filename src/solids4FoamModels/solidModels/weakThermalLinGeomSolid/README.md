---
sort: 12
---

# weakThermalLinGeomSolid

This page documents the weakly coupled thermal linear-geometry solid model.
The runtime type is:

```text
weakThermalLinearGeometry
```

The model assumes small strains and small rotations. It solves the energy
equation for temperature `T` to convergence, then solves the momentum equation
for total displacement `D` once, with no outer iterations between the two.
This one-way coupling is what "weak" refers to.

The stress is computed by the run-time selectable mechanical law.

---

## User Guide

### What this model solves

`weakThermalLinGeomSolid`:

- solves the transient heat equation `rhoC*ddt(T) == laplacian(k, T)`, with no
  mechanical work term;
- then hands over to
  [linGeomTotalDispSolid](https://www.solids4foam.com/documentation/solid-models/linGeomTotalDispSolid.html)
  for the momentum equation, from which it inherits;
- takes `C` and `k` from a `thermalModel`, configured in
  `constant/thermalProperties`.

Use this model when the temperature field drives the deformation but the
deformation does not appreciably change the temperature field, which covers
most thermoelastic problems. Use
[thermalLinGeomSolid](https://www.solids4foam.com/documentation/solid-models/thermalLinGeomSolid.html)
if the coupling needs to be converged within each time step.

```warning
As for the other thermal solid models, the momentum equation does not include
a thermal-expansion term by itself. To get thermal stresses, select a
mechanical law that consumes the temperature field, such as
`thermoMechanicalLaw`.
```

### Supported solution algorithms

The heat equation is always solved with a segregated implicit loop. The
momentum equation is solved by `linGeomTotalDispSolid`, so the
`solutionAlgorithm` entry applies to it as usual: `implicitSegregated`
(default), `PETScSNES` or `explicit`.

### Model options

| Entry | Default | Description |
| --- | --- | --- |
| `absoluteTemperatureTolerance` | `1e-06` | Convergence threshold on `T` |

`absoluteTemperatureTolerance` is an absolute change in degrees: once the
temperature stops moving by more than this between outer iterations, the
thermal loop is considered converged.

All of the `linGeomTotalDispSolid` and base `solidModel` options apply as
well; see that page for the full list. The entries most relevant to the
thermal loop are:

| Entry | Default | Relevance |
| --- | --- | --- |
| `nCorrectors` | `10000` | Maximum correctors, both loops |
| `solutionTolerance` | `1e-06` | Primary convergence tolerance |
| `alternativeTolerance` | `1e-07` | Secondary convergence tolerance |
| `infoFrequency` | `100` | Frequency for solver progress output |

The temperature loop always performs at least two iterations, and converges
when the absolute change in `T` falls below
`absoluteTemperatureTolerance`, or when the solver and relative residuals meet
`solutionTolerance` or `alternativeTolerance`.

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
solidModel        weakThermalLinearGeometry;

weakThermalLinearGeometryCoeffs
{
    solutionAlgorithm            implicitSegregated;

    nCorrectors                  10000;
    solutionTolerance            1e-06;
    alternativeTolerance         1e-07;
    materialTolerance            1e-05;
    absoluteTemperatureTolerance 1e-06;
    infoFrequency                100;
}
```

`fvSchemes` must provide `laplacian(k,T)` in addition to the entries
`linGeomTotalDispSolid` needs, and `fvSolution` must provide solvers for `T`
and `D`.

### Boundary conditions

Displacement boundaries use the standard solids4foam patch types. Temperature
boundaries use the usual OpenFOAM scalar patch types, plus `thermalConvection`
for a convective heat-transfer boundary.

Unlike `thermalLinGeomSolid` and `thermalSolid`, this model does not implement
the `setTemperatureAndHeatFlux()` interface, so it cannot be used as the solid
side of a thermal fluid-solid interface.

### Field glossary

- `T`: temperature; written each output time.
- `grad(T)`: temperature gradient, used to report the heat flux.
- `heatFlux`: `-k*grad(T)`, constructed and written by `writeFields()`.
- `rhoC`: product of density and specific heat capacity.
- All the fields of `linGeomTotalDispSolid`: `D`, `DD`, `pointD`, `pointDD`,
  `grad(D)`, `grad(DD)`, `U`, `sigma`.

Each write also prints the maximum and minimum temperature and the maximum
magnitude of the heat flux.

---

## Developer Notes

### Class role

`weakThermalLinGeomSolid` inherits from `linGeomTotalDispSolid` rather than
from `solidModel` directly, and adds a `thermalModel` member. That inheritance
is the whole design: the thermal part is a prefix to an otherwise unmodified
mechanical solve.

- `D` remains the primary solution variable;
- `nonLinGeom()` returns `LINEAR_GEOMETRY`;
- the heat equation is uncoupled from the stress state, so there is no
  `sigma && grad(U)` term.

### Construction

The constructor calls the `linGeomTotalDispSolid` constructor, constructs the
`thermalModel` from the mesh, builds `rhoC_` as
`thermal.C()*mechanical().rho()` and `k_` from `thermal.k()`, reads `T`
(`MUST_READ`), creates `grad(T)`, reads `absoluteTemperatureTolerance`, and
forces creation of `T.oldTime()`.

### `evolve()`

1. Loop on the heat equation `rhoC*ddt(T) == laplacian(k, T)`: store `T`,
   relax and solve, relax `T`, update `grad(T)`; repeat until the private
   `converged()` overload is satisfied or `nCorrectors` is reached.
2. Call `linGeomTotalDispSolid::evolve()` to solve the momentum equation.

There is no outer loop around the two, which is the difference from
`thermalLinGeomSolid`. `evolve()` always returns `true`.

### `writeFields()`

Reports the temperature range, constructs and writes the `heatFlux` field, and
then calls `linGeomTotalDispSolid::writeFields()`.

### Extension points

- `evolve()` if the heat equation needs source terms or a different coupling;
- the private `converged()` overload for a different thermal criterion.

Because the class derives from `linGeomTotalDispSolid`, any change to that
model's momentum solve is inherited here automatically.

---

## Tutorials

No tutorial selects `weakThermalLinearGeometry` by default. It is offered as a
commented alternative in the thermoelasticity cases, which are the natural
places to try it:

- `solids/thermoelasticity/hotCylinder`
- `solids/thermoelasticity/hotSphere`
- `solids/thermoelasticity/slabCooling`

Those cases run `thermalLinearGeometry` as shipped; switching the `solidModel`
entry to `weakThermalLinearGeometry` is a one-line change.
