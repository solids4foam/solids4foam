---
sort: 18
---

# MooneyRivlinElastic

Three-parameter Mooney-Rivlin (Rivlin-Saunders) hyperelasticity, with a
volumetric-isochoric split. The runtime type is:

```text
MooneyRivlinElastic
```

The model adds a dependence on the second invariant to the neo-Hookean form.
The dictionary layout common to all laws is described in the
[material models page](../../../materialModels/README.md).

---

## User Guide

### What it computes

At each integration point, from `F` and `J = det(F)`:

```text
isoB = J^(-2/3)*symm(F & F.T())
I1   = tr(isoB)
I2   = 0.5*(I1^2 - tr(isoB & isoB))
s    = 2*(c10 + c11*(I2 - 3))*isoB - 2*(c01 + c11*(I1 - 3))*inv(isoB)
p    = 0.5*K*(J^2 - 1)
sigma = (dev(s) + p*I)/J
```

which is the Cauchy stress of

```text
Psi = c10*(I1 - 3) + c01*(I2 - 3) + c11*(I1 - 3)*(I2 - 3) + U(J)
U(J) = 0.25*K*(J^2 - 1 - 2*ln(J))
```

`c11 = 0` gives the two-parameter Mooney-Rivlin model, and `c01 = c11 = 0` a
neo-Hookean model with `mu = 2*c10`. The law is finite-strain only; a linear
geometry solid model stops with an error. A point with `J <= sqrt(SMALL)` is
a fatal error. The law holds no history.

### Model options

| Entry | Required | Default | Description |
| --- | --- | --- | --- |
| `rho` | yes | - | Density, `[1 -3 0 0 0 0 0]` |
| `c10` | yes | - | First coefficient, `[1 -1 -2 0 0 0 0]` |
| `c01` | yes | - | Second coefficient, `[1 -1 -2 0 0 0 0]` |
| `c11` | yes | - | Coupled coefficient, `[1 -1 -2 0 0 0 0]` |
| `K` | `K` or `nu` | - | Bulk modulus, `[1 -1 -2 0 0 0 0]` |
| `nu` | `K` or `nu` | - | Poisson's ratio, used if `K` is absent |
| `solvePressureEqn` | no | `no` | Hydrostatic smoothing, see below |
| `pressureSmoothingScaleFactor` | no | `100` | Scale of that smoothing |

All three coefficients are required; give `c11 c11 [1 -1 -2 0 0 0 0] 0;` if
the coupled term is not wanted. Their dimensions, and those of `K` and `nu`,
are checked.

Giving neither `K` nor `nu` is a fatal error. If both are given, `K` is used
and a warning is printed. When only `nu` is given, `K` follows from the
linearised Young's modulus:

```text
E = 6*(c10 + c01),   K = E/(3*(1 - 2*nu))
```

and `nu >= 0.5 - SMALL` is a fatal error; give `K` directly for a nearly
incompressible material. The linearised shear modulus is always
`mu = 2*(c10 + c01)`. The top-level `planeStress` entry is not used by this
law: there is no plane-stress reduction.

### Tangents and volumetric split

- Scalar tangent: `(4/3)*mu + K`; the deviatoric scalar tangent is
  `(4/3)*mu`. These are constant, from the linearised moduli.
- Fourth-order tangent: `fourthOrderFiniteDifference` only. Asking for the
  analytical `fourthOrder` tangent is a fatal error.
- Volumetric split: yes. The isochoric stress is `dev(s)/J` and the
  volumetric response is `p/J = dU/dJ`. The isochoric stress is invariant
  under a superposed dilation.

The split lets the law be used with `solvePressureEqn` and with a mixed
displacement-pressure formulation, as for
[`neoHookeanElastic`](../neoHookeanElasticMechanicalConstitutiveLaw/README.md).
The `U(J)` above is the one the pressure equation of
[`nonLinearGeometryTotalLagrangianTotalDisplacement`](../../../solidModels/nonLinGeomTotalLagTotalDispSolid/README.md)
assumes.

### Hydrostatic stress smoothing

`solvePressureEqn yes` asks the solid model to replace the explicit
hydrostatic stress `p = 0.5*K*(J^2 - 1)` with a smoothed field,
`sigmaHyd`, which is written. It is done by the solid model, not the law, and
only with `nonLinearGeometryUpdatedLagrangian`,
`nonLinearGeometryTotalLagrangianTotalDisplacement` or
`linearGeometryTotalDisplacement`, the `implicitSegregated` algorithm, a
single material, and without `solvePressure`. Any other combination is
refused at start-up.

### State variables

None. The law writes no restart files and no state fields.

### Example

From the `cylinderCrush` hyperelasticity tutorial:

```text
planeStress     no;

mechanical
(
    rubber
    {
        type            MooneyRivlinElastic;
        rho             rho [1 -3 0 0 0 0 0] 1000;
        K               K [1 -1 -2 0 0 0 0] 1.410e+9;
        c10             c10 [1 -1 -2 0 0 0 0] 0.293e+6;
        c01             c01 [1 -1 -2 0 0 0 0] 0.177e+6;
        c11             c11 [1 -1 -2 0 0 0 0] 0;
        solvePressureEqn yes;
        pressureSmoothingScaleFactor 100;
    }
);
```

### Differences from the legacy law

- The coefficients are uniform. `c10`, `c01` and `c11` fields in the start
  time directory are no longer read, and are no longer written.
- No initial stress is carried: a `sigma0` field is not read, and
  `symm(F & sigma0 & F.T())` is no longer added to the stress.
- `regionName` is no longer read.
- Giving both `K` and `nu` now prints a warning; `K` is still used. A `nu`
  at or above `0.5 - SMALL` is now a fatal error.
- The hydrostatic smoothing is done by the solid model rather than the law,
  and is refused where it is not supported rather than ignored.
- A fourth-order tangent is available by finite differences.

---

## Tutorials

- [cylinderCrush](../../../../../tutorials/solids/hyperelasticity/cylinderCrush/README.md)
- [cylindricalPressureVessel](../../../../../tutorials/solids/hyperelasticity/cylindricalPressureVessel/README.md),
  given `nu` rather than `K`, in the base case and in
  `caseOptions/displacement`
- [longWall](../../../../../tutorials/solids/hyperelasticity/longWall/README.md)

---

## References

- M. Mooney, _A theory of large elastic deformation_, Journal of Applied
  Physics, 11, 582-592, 1940,
  [10.1063/1.1712836](https://doi.org/10.1063/1.1712836).
- R.S. Rivlin and D.W. Saunders, _Large elastic deformations of isotropic
  materials VII. Experiments on the deformation of rubber_, Philosophical
  Transactions of the Royal Society of London A, 243, 251-288, 1951,
  [10.1098/rsta.1951.0004](https://doi.org/10.1098/rsta.1951.0004).
