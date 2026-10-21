---
sort: 28
---

# viscoActiveHolzapfelOgden

The orthotropic Holzapfel-Ogden law for myocardium, with viscosity and a
tabulated active stress along the fibres, as defined by the cardiac
elastodynamics benchmark of Aróstica et al. (2025). The runtime type is:

```text
viscoActiveHolzapfelOgden
```

The response is anisotropic about a fibre direction `f0` and a sheet
direction `s0`, given either as fields or as single vectors.

Created by Claude Opus 5.5, supervised by Chanikya Valeti, UCD.

---

## User Guide

### What it computes

This is a **finite-strain** law: it is used with a nonlinear geometry solid
model. The second Piola-Kirchhoff stress is

```text
S = dPsi/dE + eta*Edot + tau(t)*f0*f0
```

with the passive energy

```text
Psi = a/(2b)*(exp(b*(I1bar - 3)) - 1)
    + sum_{i = f, s} a_i/(2b_i)*chi(I4i)*(exp(b_i*(I4i - 1)^2) - 1)
    + a_fs/(2b_fs)*(exp(b_fs*I8fs^2) - 1)
    + kappa/4*(J^2 - 1 - 2*ln(J))
```

where `I1bar = J^(-2/3)*tr(C)` is isochoric, while `I4f = f0 & C & f0`,
`I4s = s0 & C & s0` and `I8fs = f0 & C & s0` are formed on the full right
Cauchy-Green tensor `C`, as in the benchmark. `chi` switches off the fibre
and sheet families in compression:

```text
chi(x) = 1/(1 + exp(-k*(x - 1)))   fibreSwitch smooth (default)
chi(x) = 1 if x > 1, else 0        fibreSwitch sharp
chi(x) = 1                         fibreSwitch none
```

The derivative of the smooth switch is included in the stress, as in the
reference implementations, which differentiate the energy automatically.

- **Viscosity**: `S_visco = eta*Edot`, with the Green-Lagrange strain rate
  taken over the time step from the old-time deformation gradient,
  `Edot = (E(F) - E(F0))/dt`. It is skipped before the first time step.
- **Active stress**: `tau(t)` is uniform in space and read from an
  interpolation table (`activeStressSeries`), for example the solution of the
  benchmark's ODE, eq. (5), tabulated once in time. Without a table the law
  is passive.

The exponentials are evaluated with their argument clamped at 500, so that an
extreme trial state in a Newton line search gives a large, finite stress
rather than `inf`. A deformation gradient with `J <= sqrt(SMALL)` stops the
run.

### Model options

| Entry | Required | Default | Description |
| --- | --- | --- | --- |
| `rho` | yes | - | Density, `[1 -3 0 0 0 0 0]` |
| `a`, `af`, `as`, `afs` | yes | - | Stiffnesses, `[1 -1 -2 0 0 0 0]` |
| `b`, `bf`, `bs`, `bfs` | yes | - | Exponents, dimensionless, `> 0` |
| `bulkModulus` | yes | - | Penalty `kappa`, `[1 -1 -2 0 0 0 0]`, `> 0` |
| `eta` | no | `0` | Viscosity, `[1 -1 -1 0 0 0 0]`, `>= 0` |
| `fibreSwitch` | no | `smooth` | `smooth`, `sharp` or `none` |
| `k` | no | `100` | Steepness of the smooth switch |
| `tangentShearModulus` | no | see below | Scalar tangent shear modulus, Pa |
| `activeStressSeries` | no | - | Interpolation table of `tau(t)`, Pa |
| `uniformFibreField` | no | `no` | Take `f0` and `s0` from this dictionary |
| `f0`, `s0` | if uniform | - | Fibre and sheet vectors |

The default `tangentShearModulus` is `a + (af + as + afs)/6`. `f0` and `s0`
need not be of unit length: they are normalised before use.

### Fibre and sheet directions

`f0` and `s0` are prescribed state variables, found as for
[`GuccioneElastic`](../GuccioneElasticMechanicalConstitutiveLaw/README.md):
a `volVectorField` of that name in the start time directory, in `0`, or
registered by another model, or else the dictionary vectors. A zero-length
fibre or sheet stops the run at the first evaluation. The
[`setFibreField`](../../../../../applications/utilities/README.md) utility,
with `method ellipsoid`, writes both for the benchmark's truncated
ellipsoid.

### Tangents and the stabilisation scale

- `scalar`: `(4/3)*mu + kappa`, and `(4/3)*mu` for `scalarDeviatoric`, with
  `mu = tangentShearModulus`. It is not a tangent of the exponential energy.
  The default is the linearised passive shear modulus, without the viscous
  `eta/dt`.
- `fourthOrderFiniteDifference`: supported, through the base class.
- `fourthOrder` (analytical): not implemented; requesting it is a fatal
  error.

The solid model also scales its finite-volume stabilisation (e.g.
`diffStencilLaplacian`) by the scalar tangent, so `tangentShearModulus`
affects the answer, not only the convergence: the relaxed myocardium is very
soft (`a = 59 Pa`), and an estimate far above the physical stiffness makes the
wall artificially stiff. Including `eta/dt` in it made the benchmark far too
stiff.

### Mixed displacement-pressure formulation

The law provides the volumetric split: it returns `dU/dJ =
0.5*kappa*(J^2 - 1)/J` separately, so it can be used with a mixed formulation
(e.g. `nonLinearGeometryTotalLagrangianTotalDisplacement` with
`solvePressure yes`), whose solved pressure satisfies the same penalty
relation. The rest of the stress is written on the full deformation, as in
the benchmark, so it is not dilation invariant, and the law says so.

### State variables

| Name | Type | Role | Description |
| --- | --- | --- | --- |
| `f0` | vector | prescribed | Reference fibre direction |
| `s0` | vector | prescribed | Reference sheet direction |

The law has no persistent history: the viscous term uses the solid model's
old-time deformation gradient.

### Example

From the `ArosticaEtAl2025/benchmark1` tutorial (Tables 1-3 of the paper):

```text
mechanical
(
    heartTissue
    {
        type            viscoActiveHolzapfelOgden;
        rho             rho [1 -3 0 0 0 0 0] 1000;
        a               a [1 -1 -2 0 0 0 0] 59;
        b               8.023;
        af              af [1 -1 -2 0 0 0 0] 18472;
        bf              16.026;
        as              as [1 -1 -2 0 0 0 0] 2481;
        bs              11.12;
        afs             afs [1 -1 -2 0 0 0 0] 216;
        bfs             11.436;
        bulkModulus     bulkModulus [1 -1 -2 0 0 0 0] 1e6;
        eta             eta [1 -1 -1 0 0 0 0] 100;
        fibreSwitch     smooth;
        k               100;
        activeStressSeries
        {
            file        "$FOAM_CASE/constant/activeStress";
            outOfBounds clamp;
        }
    }
);
```

---

## Tutorials

- [`solids/hyperelasticity/ArosticaEtAl2025/benchmark1`](../../../../../tutorials/solids/hyperelasticity/ArosticaEtAl2025/benchmark1/README.md)
- [`solids/hyperelasticity/ArosticaEtAl2025/benchmark2`](../../../../../tutorials/solids/hyperelasticity/ArosticaEtAl2025/benchmark2/README.md)

## Reference

R. Aróstica et al., "A software benchmark for cardiac elastodynamics",
Comput. Methods Appl. Mech. Engrg. 435 (2025) 117485.
<https://doi.org/10.1016/j.cma.2024.117485>
