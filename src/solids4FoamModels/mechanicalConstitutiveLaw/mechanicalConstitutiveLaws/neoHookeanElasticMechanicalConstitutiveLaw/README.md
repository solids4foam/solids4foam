---
sort: 17
---

# neoHookeanElastic

Compressible neo-Hookean hyperelasticity, with the volumetric-isochoric split
of Simo and Hughes (1998). The runtime type is:

```text
neoHookeanElastic
```

It is the usual choice for large-strain rubber-like materials in solids4foam,
and the elastic part of `neoHookeanElasticMisesPlastic`. The dictionary layout
common to all laws is described in the
[material models page](../../../materialModels/README.md).

---

## User Guide

### What it computes

At each integration point, from the deformation gradient `F` and `J = det(F)`
supplied by the solid model:

```text
bEbar    = J^(-2/3)*symm(F & F.T())
sigmaHyd = 0.5*K*(J^2 - 1)
sigma    = (mu/J)*dev(bEbar) + (sigmaHyd/J)*I
```

which is the Cauchy stress of the stored energy

```text
Psi = 0.5*mu*(tr(bEbar) - 3) + U(J),   U(J) = 0.25*K*(J^2 - 1 - 2*ln(J))
```

The law is finite-strain only. A linear geometry solid model stops with an
error, since the law has no small-strain evaluation; use `linearElastic`
there. A point with `J <= sqrt(SMALL)` is a fatal error. The law holds no
history.

### Model options

| Entry | Required | Default | Description |
| --- | --- | --- | --- |
| `rho` | yes | - | Density, `[1 -3 0 0 0 0 0]` |
| `E` | with `nu` | - | Young's modulus, `[1 -1 -2 0 0 0 0]` |
| `nu` | with `E` | - | Poisson's ratio, dimensionless |
| `mu` | with `K` | - | Shear modulus, `[1 -1 -2 0 0 0 0]` |
| `K` | with `mu` | - | Bulk modulus, `[1 -1 -2 0 0 0 0]` |
| `solvePressureEqn` | no | `no` | Hydrostatic smoothing, see below |
| `pressureSmoothingScaleFactor` | no | `100` | Scale of that smoothing |

Exactly one of the pairs `E`/`nu` and `mu`/`K` must be given; both, or
neither, is a fatal error. Their dimensions are checked. A `mu`/`K` pair is
converted to `E` and `nu` with the 3-D relations, and then

```text
mu     = E/(2*(1 + nu))
lambda = E*nu/((1 + nu)*(1 - 2*nu))    planeStress no
lambda = E*nu/((1 + nu)*(1 - nu))      planeStress yes
K      = lambda + (2/3)*mu
```

so under `planeStress yes` the bulk modulus used is the plane-stress one, not
a `K` given in the dictionary. `nu` must satisfy `-1 < nu < 0.5`, or equal
`0.5` (below). Under `planeStress yes` `nu = 0.5` is accepted and gives a
finite bulk modulus.

### Fully incompressible material

With `nu 0.5` (to within `SMALL`) and `planeStress no` the law is fully
incompressible: `K` is held at `GREAT` and the law reports itself as
incompressible. The framework then refuses, naming the material, any
evaluation that returns the total stress or a tangent that includes the bulk
stiffness. Only a mixed displacement-pressure formulation, which replaces the
volumetric response with a solved pressure, can use it.

### Tangents and volumetric split

- Scalar tangent: `(4/3)*mu + K`; the deviatoric scalar tangent is
  `(4/3)*mu`.
- Fourth-order tangent: `fourthOrderFiniteDifference` only. Asking for the
  analytical `fourthOrder` tangent is a fatal error.
- Volumetric split: yes. The energy is written on `bEbar`, so the isochoric
  stress is `(mu/J)*dev(bEbar)` and the volumetric response is
  `dU/dJ = sigmaHyd/J`. The isochoric stress is invariant under a superposed
  dilation.

The split lets the law be used with a mixed formulation, such as
[`coupledPressureDisplacementSolid`](../../../solidModels/coupledPressureDisplacementSolid/README.md)
or
[`nonLinearGeometryTotalLagrangianTotalDisplacement`](../../../solidModels/nonLinGeomTotalLagTotalDispSolid/README.md)
with `solvePressure true`, where the solved pressure replaces `dU/dJ`. The
pressure equation of the latter assumes the `U(J)` above, which this law
satisfies. The mixed formulation is chosen by the solid model, not by this
dictionary.

### Hydrostatic stress smoothing

`solvePressureEqn yes` in the law's entry asks the solid model to replace the
explicit hydrostatic stress `J*dU/dJ = 0.5*K*(J^2 - 1)` with a smoothed field,
`sigmaHyd`, which is written at each write time. The smoothing diffusivity is
`pressureSmoothingScaleFactor` times the interpolate of `impK/DEqnA`. It is
done by the solid model, not the law, and is available only:

- with `nonLinearGeometryUpdatedLagrangian`,
  `nonLinearGeometryTotalLagrangianTotalDisplacement` or
  `linearGeometryTotalDisplacement`;
- with the `implicitSegregated` solution algorithm;
- for a single material, and not together with `solvePressure`.

Any other combination is refused at start-up.

### State variables

None. The law writes no restart files and no state fields.

### Example

From the `blockPunch` tutorial, with the moduli given as `mu` and `K`:

```text
planeStress     no;

mechanical
(
    steel
    {
        type            neoHookeanElastic;
        rho             rho [1 -3 0 0 0 0 0] 1000;
        mu              mu [1 -1 -2 0 0 0 0] 92.5;
        K               K [1 -1 -2 0 0 0 0] 462.416666667;
    }
);
```

The `rubberSealing` tutorial gives `E` and `nu` instead, with
`planeStress yes`, `solvePressureEqn yes` and
`pressureSmoothingScaleFactor 1e4`.

### Differences from the legacy law

- `pressureDisplacement`, `pressureDisplacementCoeff`,
  `alternatePressureDefinition`, `tangentEps` and `regionName` are no longer
  read. They are not reported as errors, so remove them from migrated cases.
- The pressure-displacement mode is gone. The stress is always built on
  `bEbar`; a mixed formulation now asks the law for its isochoric stress and
  volumetric response separately, and never reads a `p` field. The legacy
  mode used the full `b` in the deviatoric stress and `mu` as the implicit
  stiffness.
- `nu 0.5` now makes the material incompressible directly, where the legacy
  law set `K` to `GREAT` only in pressure-displacement mode. `nu` is
  validated.
- The hydrostatic smoothing is done by the solid model rather than the law,
  and is refused where it is not supported rather than ignored.
- A fourth-order tangent is available by finite differences
  (`fourthOrderFiniteDifference`); the legacy law built a face-based
  numerical tangent that needed `tangentEps`.
- `F` is not stored by the law; the solid model owns the kinematics.

---

## Tutorials

- [blockPunch](../../../../../tutorials/solids/hyperelasticity/blockPunch/README.md)
- [cantileverVibration](../../../../../tutorials/solids/hyperelasticity/cantileverVibration/README.md)
- [compressedSpheres](../../../../../tutorials/solids/hyperelasticity/compressedSpheres/README.md)
- [cooksMembrane](../../../../../tutorials/solids/hyperelasticity/cooksMembrane/README.md)
- [rubberSealing](../../../../../tutorials/solids/hyperelasticity/rubberSealing/README.md)
- [shallowIroning](../../../../../tutorials/solids/hyperelasticity/shallowIroning/README.md)
- [twistingHemisphere](../../../../../tutorials/solids/hyperelasticity/twistingHemisphere/README.md),
  with two materials
- [cylindricalPressureVessel](../../../../../tutorials/solids/hyperelasticity/cylindricalPressureVessel/README.md),
  in the three `caseOptions/pressureDisplacement` variants, which use
  `nu 0.5` with `coupledPressureDisplacementSolid`
- [plateHole](../../../../../tutorials/solids/linearElasticity/plateHole/README.md),
  in both `caseOptions/pressureDisplacement` variants
- [HronTurekFsi3](../../../../../tutorials/fluidSolidInteraction/HronTurekFsi3/README.md)
- [fillingElasticContainer](../../../../../tutorials/fluidSolidInteraction/fillingElasticContainer/README.md)
- [flexibleDamBreak](../../../../../tutorials/fluidSolidInteraction/flexibleDamBreak/README.md)

---

## References

- J.C. Simo and T.J.R. Hughes, _Computational Inelasticity_, Springer, 1998,
  [10.1007/b98904](https://doi.org/10.1007/b98904).
