---
sort: 22
---

# isotropicFungElastic

An isotropic Fung-type exponential hyperelastic law, intended for soft tissue
whose response does not depend significantly on direction. The runtime type
is:

```text
isotropicFungElastic
```

No fibre or material direction is needed. The strain energy is the
exponential of a quadratic form of the isochoric Green-Lagrange strain, in the
form introduced by Y. C. Fung for soft tissue, restricted to isotropy.

---

## User Guide

### What it computes

This is a **finite-strain** law: it implements the finite-strain update only,
so it is used with a nonlinear geometry solid model. A linear geometry solid
model asks for the small-strain update, and the run stops with an error.

The strain energy is split into an isochoric and a volumetric part:

```text
W = 0.5*c*(exp(Q) - 1) + U(J)

Fbar = J^(-1/3)*F
Ebar = 0.5*(Fbar^T & Fbar - I)        (isochoric Green-Lagrange strain)

Q    = a*tr(Ebar & Ebar)
```

`Ebar` equals the Green-Lagrange strain when `J = 1`. The stress follows as

```text
S     = dW/dEbar = c*a*exp(Q)*Ebar
sigma = dev(symm(Fbar & S & Fbar^T))/J + dU/dJ*I

U(J)  = 0.25*K*(J^2 - 1 - 2*ln(J))
dU/dJ = 0.5*K*(J^2 - 1)/J
```

where `K` is `bulkModulus`, a penalty enforcing near incompressibility. The
energy and the stress are zero in the undeformed state.

`c` scales the stress linearly at all strains; `a`, inside the exponent, sets
how quickly the material stiffens with strain. At small strain only their
product matters, and the law reduces to linear elasticity with shear modulus

```text
mu = c*a/2
```

`c` is therefore not the small-strain modulus, and its value depends on the
form of `Q` it was fitted with.

A deformation gradient with `J <= sqrt(SMALL)` stops the run.

### Model options

| Entry | Required | Default | Description |
| --- | --- | --- | --- |
| `rho` | yes | - | Density, `[1 -3 0 0 0 0 0]` |
| `c` | yes | - | Stress scale multiplying the exponential, pressure units |
| `a` | yes | - | Coefficient in the exponent `Q`, dimensionless |
| `bulkModulus` | yes | - | Volumetric penalty `K`, pressure units, `> 0` |

A value of `c`, `a` or `bulkModulus` that is zero or negative is a fatal
error.

### State variables

The law has no state variables and no persistent history.

### Tangents

- `scalar`: `(4/3)*mu + K`, with `mu = c*a/2`. It is a constant stiffness
  estimate for the segregated solver, not a tangent of the exponential
  energy.
- `scalarDeviatoric`: `(4/3)*mu`, for mixed displacement-pressure
  formulations.
- `fourthOrderFiniteDifference`: supported, through the base class.
- `fourthOrder` (analytical): not implemented; requesting it is a fatal
  error.

### Mixed displacement-pressure formulation

The law provides the volumetric split: because `Q` is written on the
isochoric strain, it can return the isochoric stress and `dU/dJ` separately.
It can therefore be used with a mixed formulation, such as
[`coupledPressureDisplacementSolid`](../../../solidModels/coupledPressureDisplacementSolid/README.md)
or
[`nonLinearGeometryTotalLagrangianTotalDisplacement`](../../../solidModels/nonLinGeomTotalLagTotalDispSolid/README.md)
with `solvePressure true`, where the solved pressure replaces `dU/dJ`.
`coupledPressureDisplacementSolid` uses `1/bulkModulus` as the
compressibility of its pressure equation, so a very large `bulkModulus` gives
a fully incompressible material.

The mixed formulation is chosen by the solid model, not by this dictionary.

### Relation to GuccioneElastic

This law is identical to
[`GuccioneElastic`](../GuccioneElasticMechanicalConstitutiveLaw/README.md)
with

```text
k = c,    cf = ct = cfs = a
```

for which the Guccione exponent reduces to `a*tr(Ebar & Ebar)` and the fibre
direction has no effect. The two laws share no strain-energy code, so this
equivalence is a check of either.

### Using published Fung coefficients

Fung coefficients fitted to planar biaxial (membrane) data cannot be used
directly. A membrane fit gives `Q` in the two in-plane strains only, with the
through-thickness strain eliminated by incompressibility; that strain enters
a three-dimensional `Q` nonlinearly, so no choice of `a` reproduces the
membrane `Q`. Calibrate `c` and `a` together against stress-stretch curves
from the published fit, over the strain range of interest.

An isotropic incompressible material also implies, at small strain, that the
membrane cross coefficient is half the normal coefficients. A published fit
far from that is anisotropic, and this law will fit it less well.

### Differences from the legacy law

The legacy `isotropicFungElastic` law in `materialModels` is a different
material: its energy, `(c1/c2)*(exp(0.5*c2*(I1bar - 3)) - 1)`, is the Demiray
law, whose exponent is linear in strain rather than quadratic.

- Its entries `c1`, `c2`, `K` and `nu` are not read. A legacy dictionary
  stops with an error for the missing `c`, rather than being reinterpreted.
- `bulkModulus` is always required; the bulk modulus cannot be derived from
  `nu`.
- A spatially varying `c1` field and an initial stress `sigma0` are not
  supported.

### Example

The material of the `LandEtAl2015/problem2` tutorial (`k = 10e3`,
`cf = ct = cfs = 1`), written for this law:

```text
planeStress no;

mechanical
(
    tissue
    {
        type        isotropicFungElastic;
        rho         rho [ 1 -3 0 0 0 0 0 ] 3000;
        c           c [ 1 -1 -2 0 0 0 0 ] 10e3;
        a           1.0;
        bulkModulus bulkModulus [ 1 -1 -2 0 0 0 0 ] 1e6;
    }
);
```

This gives `mu = 5e3`, so `bulkModulus` is 200 times the shear modulus.

---

## Tutorials

- [`solids/hyperelasticity/LandEtAl2015/problem2`](../../../../../tutorials/solids/hyperelasticity/LandEtAl2015/problem2),
  with `./Allrun isotropicFung`, which runs the material given above in place
  of the equivalent `GuccioneElastic` law of the default `petsc` option

---

## References

- Y. C. Fung, _Biomechanics: Mechanical Properties of Living Tissues_, 2nd
  ed., Springer, 1993.
- C. J. Chuong and Y. C. Fung, _Three-dimensional stress distribution in
  arteries_, Journal of Biomechanical Engineering, 105 (1983) 268-274.
- S. Land et al., _Verification of cardiac mechanics software: benchmark
  problems and solutions for testing active and passive material behaviour_,
  Proceedings of the Royal Society A, 471 (2015) 20150641.
- G. A. Holzapfel, _Nonlinear Solid Mechanics: A Continuum Approach for
  Engineering_, Wiley, 2000.
