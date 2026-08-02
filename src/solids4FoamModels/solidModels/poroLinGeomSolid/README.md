---
sort: 14
---

# poroLinGeomSolid

This page documents the poroelastic linear-geometry solid model. The runtime
type is:

```text
poroLinearGeometry
```

The model assumes small strains and small rotations. It solves a pore-pressure
equation for `p` and the linear momentum equation for total displacement `D`
inside one outer loop, so that the pressure-displacement coupling is converged
within each time step.

The approach follows `elastoPlasticBiotFoam` from Tian Tang's miniGeotechFoam
toolbox. See:

- T. Tang, O. Hededal and P. Cardiff (2014), _On finite volume method
  implementation of poro-elasto-plasticity soil model_, International Journal
  for Numerical and Analytical Methods in Geomechanics,
  [10.1002/nag.2361](https://doi.org/10.1002/nag.2361).
- T. Tang and O. Hededal (2014), _Simulation of pore pressure accumulation
  under cyclic loading using finite volume method_, Proceedings of NUMGE14,
  Volume 2, pages 1301-1306.

---

## User Guide

### What this model solves

`poroLinGeomSolid`:

- solves the pore-pressure equation
  `(porosity*rKprime)*ddt(p) + div(U) == (k/gammaWater)*laplacian(p)`,
  where `rKprime` is the reciprocal of the equivalent bulk modulus of the pore
  fluid;
- solves the linear-geometry momentum equation for `D`;
- alternates the two equations inside a single outer loop, and requires both
  to converge before the time step ends;
- couples the two through `div(U)`, the volumetric strain rate that appears as
  a source in the pressure equation.

```warning
The pore pressure enters the momentum equation only through the mechanical
law. Pair this model with `poroMechanicalLaw`, which subtracts the effective
pore pressure from the total stress; the momentum equation itself has no
explicit `p` term.
```

### Supported solution algorithms

This model implements a single, segregated implicit algorithm. It does not
read `solutionAlgorithm`, and there is no PETSc SNES or explicit path.

### Model options

These five entries are **required** — they are read with `lookup()`, so the
model fails at construction if any is missing:

| Entry | Dimensions | Description |
| --- | --- | --- |
| `hydraulicConductivity` | `[0 1 -1 0 0 0 0]` | Darcy hydraulic conductivity |
| `waterSpecificWeight` | `[1 -2 -2 0 0 0 0]` | Pore fluid specific weight |
| `porosity` | `[0 0 0 0 0 0 0]` | Porosity of the solid skeleton |
| `degreeOfSaturation` | `[0 0 0 0 0 0 0]` | Degree of saturation, 0 to 1 |
| `waterBulkModulus` | `[1 -1 -2 0 0 0 0]` | Bulk modulus of the pore fluid |

From the last two the model forms

```text
rKprime = saturation/KWater + (1 - saturation)/1e5 Pa
```

which is the compressibility of a partially saturated pore fluid; the
`1e5 Pa` is a hard-coded atmospheric pressure.

The relevant inherited `solidModel` entries are:

| Entry | Default | Relevance |
| --- | --- | --- |
| `nCorrectors` | `10000` | Maximum number of outer correctors |
| `solutionTolerance` | `1e-06` | Primary convergence tolerance, `p` and `D` |
| `alternativeTolerance` | `1e-07` | Secondary convergence tolerance |
| `materialTolerance` | `1e-05` | Mechanical-law convergence tolerance |
| `relaxationMethod` | `fixed` | Under-relaxation method (`fixed`, `aitken`) |
| `infoFrequency` | `100` | Frequency for solver progress output |
| `stabilisation` | auto-created if absent | `momentum` sub-dictionary is used |
| `cellDisplacements` | optional | Internal-cell displacement constraints |

### Required input files

- `p` and `D` in the time directory;
- `constant/solidProperties`;
- `constant/mechanicalProperties`, normally selecting `poroMechanicalLaw`;
- `constant/g`.

`p` is `MUST_READ`, so a case that omits it will fail at construction.

### Recommended dictionary setup

Minimal example for `constant/solidProperties`:

```text
solidModel     poroLinearGeometry;

poroLinearGeometryCoeffs
{
    nCorrectors          1000;
    solutionTolerance    1e-06;
    alternativeTolerance 1e-07;
    materialTolerance    1e-05;
    infoFrequency        100;

    hydraulicConductivity hydraulicConductivity [0 1 -1 0 0 0 0] 0.001;
    porosity              porosity [0 0 0 0 0 0 0] 0.2;
    waterSpecificWeight   waterSpecificWeight [1 -2 -2 0 0 0 0] 1e+04;
    degreeOfSaturation    degreeOfSaturation [0 0 0 0 0 0 0] 0.99;
    waterBulkModulus      waterBulkModulus [1 -1 -2 0 0 0 0] 2e+09;
}
```

`fvSchemes` must provide `laplacian(DD,D)` and a default `laplacian` entry for
the pressure equation, and `fvSolution` must provide solvers for `p` and `D`.

### Boundary conditions

Displacement boundaries use the standard solids4foam patch types. Pore
pressure uses the usual OpenFOAM scalar patch types: `fixedValue` for a
drained boundary, `zeroGradient` for an undrained one.

Note that on a traction boundary the `pressure` entry is the mechanical
traction pressure, not the pore pressure; the two are independent.

### Field glossary

- `p`: pore pressure; primary hydraulic unknown, written each output time.
- `grad(p)`: pore-pressure gradient.
- `D`, `DD`: total and incremental displacement.
- `pointD`, `pointDD`: displacement at mesh points.
- `grad(D)`, `grad(DD)`: displacement gradients.
- `sigma`: symmetric stress tensor; with `poroMechanicalLaw` this is the total
  stress, i.e. effective stress minus the pore-pressure contribution.
- `U`: velocity field, computed as `ddt(D)`; the `div(U)` term is what couples
  the momentum equation back into the pressure equation.

---

## Developer Notes

### Class role

`poroLinGeomSolid` inherits directly from `solidModel`. The key design choices
are:

- `D` is the primary solution variable (`solutionD()` returns `D()`);
- `nonLinGeom()` returns `LINEAR_GEOMETRY`;
- the momentum equation is the same deferred-correction Laplacian form used by
  `linGeomTotalDispSolid`, with the momentum stabilisation term applied
  explicitly;
- the poroelastic coupling is one-directional in each equation:
  displacement enters the pressure equation through `div(U)`, and pressure
  enters the momentum equation only through the mechanical law's stress.

### Construction

The constructor calls the base `solidModel` constructor, reads `p`
(`MUST_READ`), creates `grad(p)`, reads the five required poroelastic
properties, forms `rKprime_`, and forces creation of `p.oldTime()`. The
implicit stiffness fields `impK_`, `impKf_` and `rImpK_` come from the
mechanical law.

### `evolve()`

Each outer iteration:

1. stores `p`, assembles and solves the pressure equation, relaxes `p`, and
   updates `grad(p)`;
2. stores `D`, updates the momentum stabilisation, assembles

   ```text
   rho*d2dt2(D) == laplacian(impKf, D) - fvc::laplacian(impKf, D)
                 + fvc::div(sigma) + rho*g + impK*stabilisation
   ```

   then relaxes and solves it, and under-relaxes `D`;
3. updates `DD`, `grad(D)`, `grad(DD)`, then `U` — which the next pressure
   equation consumes — and `sigma`;
4. refreshes `impKf_` from the mechanical law every tenth iteration. `impK_`
   and `rImpK_` are deliberately _not_ refreshed, because they are used by the
   traction boundary conditions.

The loop exits through the private `converged()` overload, which checks both
`p` and `D`, or when `iCorr` reaches `nCorrectors`. Afterwards `pointD` and
`pointDD` are updated once. `evolve()` always returns `true`.

### Traction boundary treatment

`tractionBoundarySnGrad()` is the standard linear-geometry form,

```text
((traction - n*pressure) - (n & (sigma - impK*gradD)))*rImpK
```

evaluated with the cell-centred stress and gradient.

### Extension points

- `evolve()` for a different coupling strategy between the two equations;
- the private `converged()` overload if a different convergence criterion is
  needed.

The poroelastic properties are currently read as uniform `dimensionedScalar`
values from `solidProperties`. Spatially varying properties, or a
nonlinear permeability, would need those members to become fields.

---

## Tutorials

Cases that select `poroLinearGeometry`:

- `solids/poroelasticity/stripFooting`
- `solids/poroelasticity/rodAndSeabed`
- `solids/poroelasticity/suctionCaission`
