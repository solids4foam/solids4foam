---
sort: 3
---

# Material models

A material model describes how a material responds, as opposed to a solid
model, which describes how the governing equations are discretised and
solved. solids4foam splits this into two independent models:

- the **mechanical model**, the `mechanicalConstitutiveLaw` framework, which
  turns deformation into stress and provides the density and the implicit
  stiffness used by the momentum equation. It is read from
  `constant/mechanicalProperties`, and every solid model uses it.
- the **thermal model** (`thermalModel`), which provides the specific heat
  capacity and thermal conductivity used by the heat equation. It is read from
  `constant/thermalProperties` and is only created by the solid models that
  solve for temperature.

The two are deliberately separate. A thermal analysis therefore needs both
dictionaries: the density that multiplies the specific heat capacity is taken
from the mechanical model, not the thermal one.

The mechanical model's classes live in
[`src/solids4FoamModels/mechanicalConstitutiveLaw`](https://github.com/solids4foam/solids4foam/tree/development/src/solids4FoamModels/mechanicalConstitutiveLaw),
whose README describes the framework's design.

---

## User Guide

### The `constant/mechanicalProperties` dictionary

`mechanicalProperties` is read `MUST_READ`, so it must exist for every solid
case. It has these top-level entries:

| Entry | Required | Description |
| --- | --- | --- |
| `mechanical` | yes | List of named material sub-dictionaries |
| `planeStress` | no | Default `no`; plane stress rather than plane strain/3-D |

`planeStress` is given once, for all materials, and is passed to every law;
giving it inside a law's sub-dictionary is an error.

### Selecting a mechanical law

The `mechanical` entry is a **list**, not a dictionary. Each element is a
named sub-dictionary, and the `type` keyword inside it selects the mechanical
constitutive law at runtime:

```text
planeStress     no;

mechanical
(
    steel
    {
        type            linearElastic;
        rho             rho [1 -3 0 0 0 0 0] 7854;
        E               E [1 -1 -2 0 0 0 0] 200e9;
        nu              nu [0 0 0 0 0 0 0] 0.3;
    }
);
```

For a single-material case the name (`steel` above) is arbitrary. It becomes
significant as soon as there is more than one entry — see below.

A law implements a small-strain stress update, a finite-strain one, or both.
A linear geometry solid model asks for the small-strain update and a
nonlinear geometry solid model for the finite-strain one; a law that does not
implement the one asked for stops the run with an error.

### Available laws

Each law has a page describing its entries, its state and the tutorials that
use it:

| Law | Strain | Description |
| --- | --- | --- |
| [`linearElastic`](../mechanicalConstitutiveLaw/mechanicalConstitutiveLaws/linearElasticMechanicalConstitutiveLaw/README.md) | small | Hookean elasticity, with an optional residual stress |
| [`linearElasticMisesPlastic`](../mechanicalConstitutiveLaw/mechanicalConstitutiveLaws/linearElasticMisesPlasticMechanicalConstitutiveLaw/README.md) | small | Hookean elasticity with von Mises plasticity and hardening |
| [`linearElasticMohrCoulombPlastic`](../mechanicalConstitutiveLaw/mechanicalConstitutiveLaws/linearElasticMohrCoulombPlasticMechanicalConstitutiveLaw/README.md) | small | Hookean elasticity with Mohr-Coulomb plasticity |
| [`viscousHookeanElastic`](../mechanicalConstitutiveLaw/mechanicalConstitutiveLaws/viscousHookeanElasticMechanicalConstitutiveLaw/README.md) | small | Linear viscoelasticity (Prony series) |
| [`anisotropicBiotElastic`](../mechanicalConstitutiveLaw/mechanicalConstitutiveLaws/anisotropicBiotElasticMechanicalConstitutiveLaw/README.md) | small | Orthotropic elasticity for poroelastic soils |
| [`thermoMechanicalLaw`](../mechanicalConstitutiveLaw/mechanicalConstitutiveLaws/thermoMechanicalLawMechanicalConstitutiveLaw/README.md) | small | Thermal expansion over a sub-law, from a temperature field |
| [`poroMechanicalLaw`](../mechanicalConstitutiveLaw/mechanicalConstitutiveLaws/poroMechanicalLawMechanicalConstitutiveLaw/README.md) | small | Biot effective stress over a sub-law, from a pore pressure |
| [`StVenantKirchhoffElastic`](../mechanicalConstitutiveLaw/mechanicalConstitutiveLaws/StVenantKirchhoffElasticMechanicalConstitutiveLaw/README.md) | finite | St. Venant-Kirchhoff hyperelasticity |
| [`neoHookeanElastic`](../mechanicalConstitutiveLaw/mechanicalConstitutiveLaws/neoHookeanElasticMechanicalConstitutiveLaw/README.md) | finite | Neo-Hookean hyperelasticity |
| [`neoHookeanElasticMisesPlastic`](../mechanicalConstitutiveLaw/mechanicalConstitutiveLaws/neoHookeanElasticMisesPlasticMechanicalConstitutiveLaw/README.md) | finite | Neo-Hookean elasticity with von Mises plasticity |
| [`MooneyRivlinElastic`](../mechanicalConstitutiveLaw/mechanicalConstitutiveLaws/MooneyRivlinElasticMechanicalConstitutiveLaw/README.md) | finite | Mooney-Rivlin hyperelasticity |
| [`OgdenElastic`](../mechanicalConstitutiveLaw/mechanicalConstitutiveLaws/OgdenElasticMechanicalConstitutiveLaw/README.md) | finite | Three-term Ogden hyperelasticity |
| [`GuccioneElastic`](../mechanicalConstitutiveLaw/mechanicalConstitutiveLaws/GuccioneElasticMechanicalConstitutiveLaw/README.md) | finite | Transversely isotropic myocardium (Guccione) |
| [`HolzapfelGasserOgdenElastic`](../mechanicalConstitutiveLaw/mechanicalConstitutiveLaws/HolzapfelGasserOgdenElasticMechanicalConstitutiveLaw/README.md) | finite | Fibre-reinforced arterial tissue (HGO) |
| [`electroMechanicalLaw`](../mechanicalConstitutiveLaw/mechanicalConstitutiveLaws/electroMechanicalLawMechanicalConstitutiveLaw/README.md) | finite | Active fibre tension over a passive sub-law |

`rho`, the density, is read by every law. The optional `solvePressureEqn` and
`pressureSmoothingScaleFactor` (default `100`) entries of a law ask the solid
model to smooth the hydrostatic stress; see the solid model pages for which
solid models support it. It is supported for a single material only.

### Multi-material cases

Give the `mechanical` list more than one entry to assign different laws to
different parts of the mesh:

```text
planeStress     yes;

mechanical
(
    outer
    {
        type            linearElastic;
        rho             rho [1 -3 0 0 0 0 0] 1000;
        E               E [1 -1 -2 0 0 0 0] 200e+9;
        nu              nu [0 0 0 0 0 0 0] 0.3;
    }
    inner
    {
        type            linearElastic;
        rho             rho [1 -3 0 0 0 0 0] 1000;
        E               E [1 -1 -2 0 0 0 0] 20e+9;
        nu              nu [0 0 0 0 0 0 0] 0.35;
    }
);
```

The rules are:

1. **The name of each entry must be the name of a `cellZone`.** The keyword of
   each `mechanical` entry (`outer`, `inner` above) is taken verbatim as a
   cellZone name, and construction aborts if no cellZone of that name exists.
2. **Every cell must be in exactly one of those cellZones.** A cell in none of
   them, or in more than one, is an error.
3. **The displacement gradient scheme must be material-aware.** The stress is
   evaluated on the whole mesh, each cell with its own law, from one
   displacement gradient. Across an interface between materials the
   displacement is continuous but its gradient jumps, so a cell's gradient
   stencil must stay within its own material: set
   `grad(D) leastSquaresS4f;` (or `grad(DD)` for the incremental solid models)
   in `fvSchemes`. The cell-centred solid models refuse a multi-material case
   with any other scheme.

Cell zones are normally created after meshing with `setSet` followed by
`setsToZones` (or with `topoSet`). The `layeredPipe` tutorial does exactly
this, with a `batch.setSet` file that makes an `outer` and an `inner` cellSet
and subtracts one from the other so the two do not overlap.

```note
Because the materials are assigned by cellZone, a material region does not
need to be contiguous. It does, however, have to be a single mesh: separately
meshed parts must be merged first, as in the `punch` tutorial.
```

### Multi-material limitations

Some capabilities are single-material only, and fail with a clear error
rather than silently giving a wrong answer:

| Capability | Status with more than one material |
| --- | --- |
| Face-quadrature (high-order) gradient | Not implemented |
| Hydrostatic stress smoothing (`solvePressureEqn`) | Refused |
| `coupledPressureDisplacementSolid` | Refused |

The vertex-centred solid model supports more than one material: each dual
mesh face takes the law of the primary cell it lies in.

### Thermal properties

Solid models that solve a temperature equation additionally read
`constant/thermalProperties`, which selects a `thermalLaw` giving the specific
heat capacity and thermal conductivity. See the thermal model page in this
section.

---

## Tutorials

Cases with more than one entry in the `mechanical` list:

- `solids/multiMaterial/layeredPipe`: a bi-material thick-walled cylinder with
  cellZones built by `setSet`/`setsToZones` and an analytical solution to
  compare against.
- `solids/linearElasticity/punch`: two separately meshed parts merged into one
  mesh, then split into the `punch_top` and `punch_bottom` cellZones.
- `solids/thermoelasticity/hotCylinder/hotCylinderPredefinedTFieldMultipleMaterials`:
  `steel` and `aluminium` materials, each a `thermoMechanicalLaw` wrapping a
  `linearElastic` law.
