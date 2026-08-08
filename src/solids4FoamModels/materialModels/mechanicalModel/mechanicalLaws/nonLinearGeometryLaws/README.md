---
sort: 14
---

# Nonlinear geometry mechanical laws

This section documents the mechanical laws that make no small-strain or
small-rotation assumption.

A mechanical law is the constitutive part of a solids4foam case: the solid
model assembles and solves the momentum equation, and the mechanical law
answers the single question _given the current displacement field, what is the
stress?_ Nonlinear geometry laws answer that question in terms of the
deformation gradient

```text
F = I + grad(D)^T
J = det(F)
```

and quantities derived from it, such as the left Cauchy-Green tensor
`B = F & F.T()`, the Green strain `G = 0.5*(F.T() & F - I)`, and the isochoric
split `isoB = J^(-2/3)*B`. Because the reference and deformed configurations
are distinct, these laws must also be explicit about which stress measure they
return — the Cauchy stress `sigma` for updated-Lagrangian laws, the second
Piola-Kirchhoff stress `S` for total-Lagrangian ones.

Laws in this directory are registered under the `nonLinGeomMechLaw` run-time
selection table and are intended for the nonlinear geometry solid models, for
example `nonLinearGeometryUpdatedLagrangian`,
`nonLinearGeometryTotalLagrangianTotalDisplacement` and
`unsNonLinearGeometryTotalLagrangian`.

```note
The linear geometry counterparts live in `linearGeometryLaws` and are
registered under a separate `linGeomMechLaw` table. A law from one table
cannot be selected by a solid model that uses the other.
```

---

## Catalogue

For every law in this subsection the run-time selection name registered by
`addToRunTimeSelectionTable` is identical to the C++ class name, so the class
name is what you type in `constant/mechanicalProperties`.

| Runtime type | Purpose |
| --- | --- |
| `StVenantKirchhoffElastic` | Hookean law in Green strain and 2nd PK stress |
| `StVenantKirchhoffOrthotropicElastic` | Nine-parameter orthotropic St Venant |
| `neoHookeanElastic` | Compressible neo-Hookean, Simo and Hughes Eqn 9.2.6 |
| `MooneyRivlinElastic` | Three-parameter Mooney-Rivlin |
| `YeohElastic` | Three-parameter Yeoh |
| `OgdenElastic` | Third-order Ogden, in principal stretches |
| `GentElastic` | Gent-Flory law for the swelling behaviour of microgels |
| `isotropicFungElastic` | Two-parameter exponential Fung-like model |
| `GuccioneElastic` | Exponential Guccione law for myocardium |
| `HolzapfelGasserOgdenElastic` | Two-fibre-family anisotropic HGO law |
| `neoHookeanElasticMisesPlastic` | Neo-Hookean with J2 plasticity |
| `viscoNeoHookeanElastic` | Neo-Hookean generalised Maxwell viscoelasticity |
| `electroMechanicalLaw` | Wrapper adding an active electro-mechanical stress |
| `diffusionHyperElastic` | Neo-Hookean scaled by a mesh motion diffusivity |

`neoHookeanElastic`, `MooneyRivlinElastic`, `YeohElastic`, `OgdenElastic`,
`GentElastic` and `isotropicFungElastic` are isotropic and purely elastic, and
differ only in their strain energy function; `StVenantKirchhoffElastic` is the
degenerate case that keeps a linear stress-strain relation while accounting for
finite rotation. The remaining anisotropic laws each need directional data in
addition to their material constants: `StVenantKirchhoffOrthotropicElastic`
takes a local coordinate system, while `GuccioneElastic` and
`HolzapfelGasserOgdenElastic` take a fibre direction field. The latter two were
added for cardiovascular modelling.

`neoHookeanElasticMisesPlastic` and `viscoNeoHookeanElastic` are the two
inelastic laws in this subsection: they carry internal state between time
steps, so they require the solid model to be run with a sensible time step
rather than as a steady-state solve.

`electroMechanicalLaw` is a _wrapper_. It owns a nested, run-time selectable
passive hyperelastic law and adds an active stress contribution to whatever
that law returns, so it can be combined with most of the other entries.

`diffusionHyperElastic` is not a physical material model; like
`diffusionElastic` in the linear geometry subsection, it exists so that a solid
model can be used as a mesh motion solver.

```warning
`HolzapfelGasserOgdenElastic` is written for the mixed pressure-displacement
formulation and is intended for `coupledPressureDisplacementSolid`. It is not
a drop-in replacement for the other laws here.
```

---

## Selecting a law

Mechanical laws are selected in `constant/mechanicalProperties`. The
`mechanical` entry is a _list_ of sub-dictionaries, one per material, so
multi-material cases are expressed by adding further entries:

```text
planeStress     no;

mechanical
(
    rubber
    {
        type            neoHookeanElastic;
        rho             rho [1 -3 0 0 0 0 0] 1100;
        E               E [1 -1 -2 0 0 0 0] 6e6;
        nu              nu [0 0 0 0 0 0 0] 0.45;
    }
);
```

The name of each sub-dictionary (`rubber` above) is arbitrary; for
multi-material cases it must match the corresponding `cellZone`.

Most of these laws accept their stiffness as either the `E`/`nu` pair or the
`mu`/`K` pair, in the same way as `linearElastic`, and then add their own
law-specific constants. The individual pages give the exact entries.

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

The top-level `planeStress` switch sits outside the `mechanical` list and is
read by the base class through the `mechanicalProperties` dictionary itself.
Several laws here consult it when converting `E` and `nu` into the constants
they actually use: with `planeStress yes` the second Lame parameter becomes
`nu*E/((1 + nu)*(1 - nu))` rather than `nu*E/((1 + nu)*(1 - 2*nu))`. It changes
only that conversion — it does not impose a plane stress constraint on the
solution — so a genuine plane stress analysis still has to be set up as a
one-cell-thick three-dimensional model with a traction-free front and back.
