---
sort: 28
---

# logStrain

A composite law that uses a small-strain law at finite strain, by driving it
with the elastic logarithmic (Hencky) strain. The runtime type is:

```text
logStrain
```

---

## User Guide

### What it computes

This is a **finite-strain** law: it implements the finite-strain update
only, so it is used with a nonlinear geometry solid model, total or updated
Lagrangian.

The sub-law, given in the `mechanicalLaw` sub-dictionary, is a small-strain
law. Over each step, with `f = F & inv(F0)` the relative deformation gradient
and `epsilonE0` the elastic strain at the start of the step, this law forms
the trial elastic logarithmic strain

```text
bETrial       = f & exp(2*epsilonE0) & f^T
epsilonETrial = 0.5*log(bETrial)
```

The sub-law performs its own return mapping on `epsilonETrial`, exactly as it
would on a small strain, and gives back the elastic strain `epsilonE` and the
stress `tau` conjugate to it, which is the Kirchhoff stress. The Cauchy
stress is `sigma = tau/J`.

This is the exponential-map return mapping of Eterovic and Bathe (1990) and
Simo (1992), as used by Whelan et al. (2025). The multiplicative split
`F = Fe & Fp` is exact, not approximated, for an isotropic sub-law with
scalar history, and the step is unconditionally objective: a rigid rotation
only rotates the stress. Under `logStrain`:

- `linearElastic` is Hencky's hyperelastic law,
  `tau = 2*mu*dev(ln V) + K*tr(ln V)*I`;
- `linearElasticMisesPlastic` is finite-strain J2 plasticity, with the yield
  function on the Kirchhoff stress and the hardening table read as Kirchhoff
  yield stress against equivalent logarithmic plastic strain. For a metal,
  where `J` stays close to one, the Kirchhoff and Cauchy stresses barely
  differ.

Only a sub-law that can be driven by a trial elastic strain can be used:
`linearElastic`, `linearElasticMisesPlastic` and
[`linearElasticMisesPlasticLemaitreDamage`](../linearElasticMisesPlasticLemaitreDamageMechanicalConstitutiveLaw/README.md).
Any other is refused when
the case starts, naming it. A law with kinematic hardening or a preferred
material direction does not qualify, since the exponential map is exact only
when the return mapping keeps the elastic strain coaxial with its trial
value.

The density and the bulk modulus are those of the sub-law. `rho` given on
this law is passed down to the sub-law when its own dictionary has none.

### Model options

| Entry | Required | Default | Description |
| --- | --- | --- | --- |
| `mechanicalLaw` | yes | - | Sub-dictionary selecting the small-strain law |
| `rho` | see text | - | Density, passed to the sub-law if absent there |

The `mechanicalLaw` sub-dictionary holds the sub-law's `type` and all of its
entries, as that law documents them.

### Restrictions

- Plane stress is refused. Solve the case in 3-D instead, with a symmetry
  plane at the back and a traction-free front.
- Prestress is refused: the elastic history starts at zero, not at the
  prestressed state.
- A non-positive Jacobian `J` is a fatal error, since an inverted point has
  no logarithmic strain.

### State variables

| Variable | Description |
| --- | --- |
| `bE` | Elastic left Cauchy-Green tensor, the history carried between steps |
| `epsilonE` | Elastic logarithmic strain, `0.5*log(bE)` |

The sub-law's state lives in a child state named `mechanicalLaw`, with
restart files named `<material>_<topology>_mechanicalLaw_<variable>`. Under
`logStrain`, `linearElasticMisesPlastic` updates `epsilonPEq` and `sigmaY`
but not `epsilonP`, which stays zero: the plastic deformation is carried by
`bE`.

### Tangents

The scalar tangents are the sub-law's small-strain ones, which steer the
segregated solid models but are not consistent tangents. No analytical
fourth-order tangent is provided; `fourthOrderFiniteDifference` is
available.

### Mixed displacement-pressure formulation

The law does not declare a volumetric split, so a mixed
displacement-pressure formulation refuses to start with it.

### Example

```text
planeStress     no;

mechanical
(
    steel
    {
        type            logStrain;
        rho             rho [1 -3 0 0 0 0 0] 7800;

        mechanicalLaw
        {
            type            linearElasticMisesPlastic;
            E               E [1 -1 -2 0 0 0 0] 200e+9;
            nu              nu [0 0 0 0 0 0 0] 0.3;
            file            "$FOAM_CASE/constant/plasticStrainVsYieldStress";
            outOfBounds     clamp;
        }
    }
);
```

### Verification

`applications/test/logStrainMechanicalConstitutiveLaw` drives the law over
multi-step paths at a material point and checks that:

- with `linearElastic`, it is Hencky's law after a path that stretches and
  rotates;
- with `linearElasticMisesPlastic`, on a path whose principal directions
  stay fixed, its Kirchhoff stress equals the small-strain law's at the
  logarithmic strain, step by step, through yield and a change of loading
  direction;
- a rigid rotation after yield only rotates the stress;
- in simple shear to a shear of 5 with perfect plasticity, the shear stress
  stays within 0.1% of `sigmaY/sqrt(3)` and of
  `neoHookeanElasticMisesPlastic`.

### References

- A.L. Eterovic and K.-J. Bathe, A hyperelastic-based large strain
  elasto-plastic constitutive formulation with combined isotropic-kinematic
  hardening using the logarithmic stress and strain measures, Int. J. Numer.
  Meth. Eng. 30 (1990) 1099-1114.
- J.C. Simo, Algorithms for static and dynamic multiplicative plasticity that
  preserve the classical return mapping schemes of the infinitesimal theory,
  Comput. Methods Appl. Mech. Eng. 99 (1992) 61-112.
- A. Whelan, P. Cardiff et al., Int. J. Numer. Meth. Eng. (2025).
