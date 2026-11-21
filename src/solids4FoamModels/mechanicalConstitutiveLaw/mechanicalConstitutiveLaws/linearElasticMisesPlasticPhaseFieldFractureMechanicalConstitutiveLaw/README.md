---
sort: 31
---

# linearElasticMisesPlasticPhaseFieldFracture

Small-strain J2 plasticity with isotropic hardening coupled to the ductile
phase field fracture model of Borden et al. (2016). The runtime type is:

```text
linearElasticMisesPlasticPhaseFieldFracture
```

This is the phase field model of Whelan et al. (2025). For finite strains,
use it as the sub-law of
[`logStrain`](../logStrainMechanicalConstitutiveLaw/README.md), which is the
formulation of that work. With a yield stress out of reach, it is the brittle
phase field model of Miehe et al. (2010).

---

## User Guide

### What it computes

The phase field `d`, between 0 for intact material and 1 for a crack,
degrades the stress by

```text
g(d) = (1 - k)(1 - d)^2 + k
```

where `k` is the residual stiffness. Its volumetric part is degraded only in
tension, so that a crack closes under compression (the split of Amor et al.,
2009):

```text
sigma = g 2 mu dev(epsilonE) + g_v K tr(epsilonE) I
g_v   = g if tr(epsilonE) >= 0, and 1 otherwise
```

The yield stress is degraded alike, so the plastic flow is that of
`linearElasticMisesPlastic`, and the stress is the intact one scaled.

The phase field is the solution of

```text
(1 + A) d - div((2 l)^2 grad(d)) = A
A = 4 l (1 - k)(H + <Wp - w0>)/Gc
```

with zero flux on the boundary, where `l` is the characteristic length,
`Gc` the critical fracture energy, `H` the history of the tensile elastic
energy,

```text
psi+ = mu dev(epsilonE):dev(epsilonE) + K <tr(epsilonE)>^2/2
H    = max(H0, psi+)
```

`Wp` the plastic work, which grows by `sigmaY dEpsilonPEq`, and `w0` the
plastic work below which the plastic work does not drive the crack. The
energies are those of the intact material, so `A` follows from the
kinematics and the old-time state alone. The plastic work is that of Whelan
et al., with the intact yield stress; Borden et al. write it with the
degraded stress. For a uniform field the equation
gives `d = A/(1 + A)`.

### The phase field equation

The solid model solves the phase field equation over the mesh, once per
outer iteration, after the stress update; the law reads the solution as the
input `phaseField`, lagged by one outer iteration, which the converged step
balances. It is an implicit-gradient field, of source `phaseFieldSource`,
`A`, and coefficient `phaseFieldCoefficient`, `1 + A`, over the length
`2 l`; see the
[framework README](../../README.md#implicit-gradient-fields).

The field `phaseField` is written, read on a restart, and needs a linear
solver entry in `system/fvSolution`:

```text
phaseField
{
    solver          PCG;
    preconditioner  FDIC;
    tolerance       1e-09;
    relTol          0;
}
```

The equation is solved by the updated Lagrangian, total Lagrangian and
linear geometry total displacement solid models with the cell-centred stress
and the `implicitSegregated` algorithm; others refuse the law. The
irreversibility of the crack is that of its driving force: `H` and `Wp`
never decrease.

### Model options

| Entry | Required | Default | Description |
| --- | --- | --- | --- |
| `rho` | yes | - | Density |
| `E`, `nu` | yes | - | Young's modulus and Poisson's ratio |
| `file` or `fileName` | yes | - | Table of yield stress against `epsilonPEq` |
| `outOfBounds` | yes | - | `clamp`, `warn`, `error` or `repeat` |
| `Gc` | yes | - | Critical fracture energy, per unit area |
| `characteristicLength` | yes | - | Length `l` |
| `w0` | no | `0` | Plastic work threshold |
| `residualStiffness` | no | `0` | Residual stiffness `k`, in `[0, 1)` |

`Gc` has dimensions of energy per unit area, `[1 0 -2 0 0 0 0]`, and `w0` of
stress. A table of two points gives linear hardening, returned in closed
form; more points give nonlinear hardening, returned by a Newton loop whose
tolerance and iteration limit are `NewtonLoopTol` (`1e-8`) and
`NewtonMaxIter` (`200`). Plane stress is not supported.

### State variables

| Variable | Description |
| --- | --- |
| `epsilonP` | Plastic strain (small-strain update only) |
| `epsilonPEq` | Equivalent plastic strain |
| `sigmaY` | Yield stress |
| `elasticEnergyHistory` | History `H` of the tensile elastic energy |
| `plasticWork` | Plastic work `Wp` |
| `phaseFieldSource` | Crack driving force `A` |
| `phaseFieldCoefficient` | `1 + A` |

Under `logStrain`, `epsilonP` stays zero: the plastic deformation is carried
by the elastic history of `logStrain`.

### Tangents

The scalar tangents are those of `linearElasticMisesPlastic` scaled by the
degradation. No analytical fourth-order tangent is provided;
`fourthOrderFiniteDifference` is available. No volumetric split is
provided, so mixed displacement-pressure formulations refuse the law.

### Diagnostics

At the end of each step, with the law's debug switch on (the default), the
yielding points, the maximum plastic increment and the largest local phase
field `A/(1 + A)` are reported.

### Example

The flat notched bar of Whelan et al. (2025), at finite strain:

```text
mechanical
(
    aluminium
    {
        type            logStrain;
        rho             rho [1 -3 0 0 0 0 0] 2700;

        mechanicalLaw
        {
            type            linearElasticMisesPlasticPhaseFieldFracture;
            E               E [1 -1 -2 0 0 0 0] 68.8e9;
            nu              nu [0 0 0 0 0 0 0] 0.33;
            file            "$FOAM_CASE/constant/plasticStrainVsYieldStress";
            outOfBounds     clamp;
            Gc              Gc [1 0 -2 0 0 0 0] 60e3;
            characteristicLength characteristicLength [0 1 0 0 0 0 0] 3.226e-4;
            w0              w0 [1 -1 -2 0 0 0 0] 1e7;
        }
    }
);
```

The tutorial `tutorials/solids/fracture/flatNotchedBarPhaseField` runs it
and compares the stress with the results of Borden et al. and of Eldahshan
et al.

### Verification

`applications/test/phaseFieldFractureMechanicalConstitutiveLaw` drives the
law at a material point, with the phase field supplied as the input, and
checks that:

- with the phase field zero, it is `linearElasticMisesPlastic`;
- the stress is degraded by `g(d)`, plastic flow included, except for its
  volumetric part in compression, and a crack keeps the residual stiffness;
- on a path that yields and unloads, `H`, `Wp`, `A` and `1 + A` are those
  recomputed independently;
- for the brittle model, with the phase field iterated to `A/(1 + A)`, the
  uniaxial strain response is the analytical one, peaking at `A = 1/3` at
  9/16 of the intact stress;
- under `logStrain`, on a path whose principal directions stay fixed, the
  Kirchhoff stress equals the small-strain law's at the logarithmic strain.

### References

- M.J. Borden, T.J.R. Hughes, C.M. Landis, A. Anvari and I.J. Lee, A
  phase-field formulation for fracture in ductile materials: finite
  deformation balance law derivation, plastic degradation, and stress
  triaxiality effects, Comput. Methods Appl. Mech. Eng. 312 (2016) 130-166.
- C. Miehe, M. Hofacker and F. Welschinger, A phase field model for
  rate-independent crack propagation, Comput. Methods Appl. Mech. Eng. 199
  (2010) 2765-2778.
- H. Amor, J.-J. Marigo and C. Maurini, Regularized formulation of the
  variational brittle fracture with unilateral contact, J. Mech. Phys.
  Solids 57 (2009) 1209-1229.
- A. Whelan, T. Tang, V. Pakrashi and P. Cardiff, A finite volume framework
  for damage and fracture prediction in wire drawing, Int. J. Numer. Meth.
  Eng. (2025), <https://doi.org/10.1002/nme.7640>.
