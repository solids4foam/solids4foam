---
sort: 2
---

# Linear geometry mechanical laws

This section documents the mechanical laws that assume the linear geometry
(small strain, small rotation) approximation.

A mechanical law is the constitutive part of a solids4foam case: the solid
model assembles and solves the momentum equation, and the mechanical law
answers the single question _given the current displacement field, what is the
stress?_ Linear geometry laws answer that question in terms of the small-strain
tensor

```text
epsilon = symm(grad(D))
```

or, for an incremental solid model, in terms of the increment
`symm(grad(DD))` accumulated onto the old-time strain. There is no distinction
between the reference and deformed configurations, no deformation gradient
`F`, and no distinction between Cauchy and Piola-Kirchhoff stress. Laws in this
directory are registered under the `linGeomMechLaw` run-time selection table
and are intended for use with the linear geometry solid models, for example
`linearGeometryTotalDisplacement`, `unsLinearGeometry` and
`poroLinearGeometry`.

```note
The nonlinear geometry counterparts live in `nonLinearGeometryLaws` and are
registered under a separate `nonLinGeomMechLaw` table. A law from one table
cannot be selected by a solid model that uses the other.
```

---

## Catalogue

For every law in this subsection the run-time selection name registered by
`addToRunTimeSelectionTable` is identical to the C++ class name, so the class
name is what you type in `constant/mechanicalProperties`.

| Runtime type | Purpose |
| --- | --- |
| `linearElastic` | Isotropic Hooke's law |
| `orthotropicLinearElastic` | Nine-parameter orthotropic Hooke's law |
| `linearElasticMisesPlastic` | Hooke's law with J2 (von Mises) plasticity |
| `linearElasticMohrCoulombPlastic` | Hooke's law with Mohr-Coulomb plasticity |
| `poroMechanicalLaw` | Wrapper adding a pore-pressure term |
| `thermoMechanicalLaw` | Wrapper adding a thermal expansion term |
| `viscousHookeanElastic` | Generalised Maxwell (Prony series) viscoelasticity |
| `anisotropicBiotElastic` | Orthotropic elasticity for soil skeletons |
| `diffusionElastic` | Hooke's law scaled by a mesh motion diffusivity |
| `linearElasticCt` | `E` from CT data; **not currently built** |
| `linearElasticFromFile` | Hooke's law with `E` read as a field from disk |

Three of these are _wrappers_ rather than stand-alone constitutive models:
`poroMechanicalLaw` and `thermoMechanicalLaw` each own a nested, run-time
selectable law and add a spherical stress contribution to whatever it returns.
They can therefore be combined with most of the other entries in the table.

`diffusionElastic` is not a physical material model; it exists so that a solid
model can be used as a mesh motion solver.

---

## Selecting a law

Mechanical laws are selected in `constant/mechanicalProperties`. The
`mechanical` entry is a _list_ of sub-dictionaries, one per material, so
multi-material cases are expressed by adding further entries:

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

The name of each sub-dictionary (`steel` above) is arbitrary; for
multi-material cases it must match the corresponding `cellZone`.

### Entries every law understands

These are read by the `mechanicalLaw` base class, so they are valid inside any
of the sub-dictionaries above:

| Entry | Default | Description |
| --- | --- | --- |
| `rho` | required | Density, `[1 -3 0 0 0 0 0]` |
| `solvePressureEqn` | `no` | Solve a Laplacian for hydrostatic pressure |
| `pressureSmoothingScaleFactor` | `100` | Scale factor for that equation |
| `regionName` | auto | Object registry holding the solid mesh |

`rho` is looked up on demand rather than at construction, so a missing `rho`
is reported the first time the solid model asks for the density.
`solvePressureEqn` is intended for nearly incompressible materials; not every
law honours it, and the individual pages say which do.

The top-level `planeStress` switch sits outside the `mechanical` list and is
read by the base class through the `mechanicalProperties` dictionary itself.
Several laws in this subsection abort if `planeStress` is `yes` — see the
individual pages.
