---
sort: 29
---

# linearElasticMisesPlasticLemaitreDamage

Small-strain J2 plasticity with isotropic hardening, coupled to Lemaitre
ductile damage. The runtime type is:

```text
linearElasticMisesPlasticLemaitreDamage
```

This is the damage model of Whelan et al. (2025), in its local form or
regularised by a non-local damage. For finite strains, use it as the sub-law
of [`logStrain`](../logStrainMechanicalConstitutiveLaw/README.md), which is
the formulation of that work.

---

## User Guide

### What it computes

The damage `D`, between 0 for virgin material and 1, degrades the stress,
and the yield function and flow rule are written on the effective stress
`sigma/(1 - D)`:

```text
sigma     = (1 - D) C : epsilonE
Phi       = sigmaEq/(1 - D) - sigmaY(r)
dEpsilonP = (3/2) dr/(1 - D) dev(sigma)/sigmaEq
```

where `r` is the accumulated plastic multiplier, which drives the hardening
and is stored as `epsilonPEq`. The damage evolves with the energy release
rate `Y`:

```text
dD = dr/(1 - D) (-Y/S)^b     if r > epsilonD and eta > -1/3
```

where `eta = tr(sigma)/(3 sigmaEq)` is the stress triaxiality. The
triaxiality cut-off (`triaxialityCutOff yes`, the default) stops damage under
compressive triaxiality, as in wire drawing. The denominator `S` is selected
by `lemaitreType`:

| `lemaitreType` | `S` |
| --- | --- |
| `classic` | `S0` |
| `Malcher` | `S0Tension/(3 abs(eta) + (S0Tension/S0Shear) (1 - xi^2))` |
| `triaxiality` | `2 S0/(1 + 3 eta)` |

where `xi` is the normalised third invariant of the deviatoric stress (the
Lode parameter). `Malcher` follows Malcher et al. (2012), with `S0Tension`
and `S0Shear` the values in uniaxial tension and pure shear. `triaxiality`
is the form proposed by Whelan et al. (2025), equation 55; it is singular at
`eta = -1/3` and so requires the cut-off.

The energy release rate is the classic

```text
Y = -[(1 + nu) sigma:sigma - nu tr(sigma)^2]/(2 E (1 - D)^2)
```

unless a crack closure parameter `h < 1` is given, when the compressive
principal stresses contribute only in proportion to `h`:

```text
Y = -[(1 + nu) sigma+:sigma+ - nu <tr(sigma)>^2]/(2 E (1 - D)^2)
  - h [(1 + nu) sigma-:sigma- - nu <-tr(sigma)>^2]/(2 E (1 - h D)^2)
```

`h = 0.2` is common for steels.

When the damage would pass the critical damage `Dc`, the point fails: its
damage is set to 0.99, and stays there, so that it carries almost no stress.

### Non-local forms

The local damage softens the material without a length scale, so once it
localises the solution depends on the mesh: the finer the mesh, the narrower
the damaged band and the earlier the failure. With `nonLocal classic` or
`nonLocal effective`, the stress is degraded instead by the implicit
non-local damage `Dbar` of Peerlings et al. (1996), the solution of

```text
Dbar - div(l^2 grad(Dbar)) = D
```

where `l` is `characteristicLength`, with zero flux on the boundary. The
solid model solves this equation over the mesh, once per outer iteration,
after the stress update; the law reads the solution as the input
`damageNonLocal`, lagged by one outer iteration, which the converged step
balances. At each point, with `Deff` the damage the stress is degraded by:

```text
Phi   = sigmaEq/(1 - Deff) - sigmaY(r)
sigma = (1 - Deff) C : epsilonE
D     = D0 + dr/(1 - Deff) (-Y/S)^b
```

with `Y` evaluated at `Deff`; the local damage `D` is explicit, since the
stress does not depend on it.

- `classic`: `Deff = Dbar`, and a local damage past `Dc` is set to 0.99.
- `effective`: the effective non-local damage of Whelan et al. (2025).
  `Deff = Dbar` until `Dbar` reaches `Dc`, and 0.99 from then on; and the
  local damage stops growing once `Deff` has reached `Dc`, so that a failed
  region neither recovers nor feeds the damage of its neighbours. The local
  damage is then not bounded by 0.99.

With `characteristicLength` zero the equation gives `Dbar = D`, and at a
converged step the non-local forms are the local one.

The field `damageNonLocal` is written, read on a restart, and needs a
linear solver entry in `system/fvSolution`:

```text
damageNonLocal
{
    solver          PCG;
    preconditioner  FDIC;
    tolerance       1e-09;
    relTol          0;
}
```

The equation is solved by the updated Lagrangian, total Lagrangian and
linear geometry total displacement solid models with the cell-centred stress
and the `implicitSegregated` algorithm; others refuse the law. With more
than one material, the source and length are zero in the cells of a
material that does not read the field, and no flux crosses into them.

### Time integration

In the local form, the return mapping and the damage update are solved
together at each point, implicitly: for a trial damage the plastic
multiplier is found by the radial return, and the damage is the root of
`D - D0 - dD(D) = 0` between the damage at the start of the step and `Dc`,
found by regula falsi. Whelan et al. instead lag the damage by one outer
iteration of the solid model; that
scheme converges to the same solution, but this one does not depend on how
many times the law has been evaluated, as a law in this framework must not.

### Model options

| Entry | Required | Default | Description |
| --- | --- | --- | --- |
| `rho` | yes | - | Density |
| `E`, `nu` | yes | - | Young's modulus and Poisson's ratio |
| `file` or `fileName` | yes | - | Table of yield stress against `r` |
| `outOfBounds` | yes | - | `clamp`, `warn`, `error` or `repeat` |
| `lemaitreType` | no | `classic` | `classic`, `Malcher` or `triaxiality` |
| `S0` | see text | - | Damage denominator, for `classic` and `triaxiality` |
| `S0Tension` | see text | - | Malcher denominator in tension |
| `S0Shear` | see text | - | Malcher denominator in shear |
| `b` | no | `1` | Damage exponent |
| `epsilonD` | no | `0` | Value of `r` below which there is no damage |
| `Dc` | no | `0.99` | Critical damage, at most 0.99 |
| `crackClosure` | no | `1` | Crack closure parameter `h`, in `[0, 1]` |
| `triaxialityCutOff` | no | `yes` | No damage for `eta <= -1/3` |
| `nonLocal` | no | `none` | `none`, `classic` or `effective` |
| `characteristicLength` | non-local | - | Length `l` of the non-local damage |

`S0`, `S0Tension` and `S0Shear` have dimensions of stress. A table of two
points gives linear hardening, returned in closed form; more points give
nonlinear hardening, returned by a Newton loop whose tolerance and iteration
limit are `NewtonLoopTol` (`1e-8`) and `NewtonMaxIter` (`200`). The coupled
damage solution stops at `damageTol` (`1e-12`) or `damageMaxIter` (`100`).
Plane stress is not supported.

### State variables

| Variable | Description |
| --- | --- |
| `epsilonP` | Plastic strain (small-strain update only) |
| `epsilonPEq` | Accumulated plastic multiplier `r` |
| `sigmaY` | Yield stress |
| `damage` | Damage `D`; the local damage in the non-local forms |

Under `logStrain`, `epsilonP` stays zero: the plastic deformation is carried
by the elastic history of `logStrain`.

### Tangents

The scalar tangents are those of `linearElasticMisesPlastic` scaled by
`1 - D`, or `1 - Deff` in the non-local forms. No analytical fourth-order
tangent is provided; `fourthOrderFiniteDifference` is available. No
volumetric split is provided, so mixed displacement-pressure formulations
refuse the law.

### Diagnostics

At the end of each step, with the law's debug switch on (the default), the
yielding points, the maximum plastic increment, the maximum damage and the
number of failed points are reported.

### Example

The notched round bar of Whelan et al. (2025), at finite strain:

```text
mechanical
(
    aluminium
    {
        type            logStrain;
        rho             rho [1 -3 0 0 0 0 0] 2700;

        mechanicalLaw
        {
            type            linearElasticMisesPlasticLemaitreDamage;
            E               E [1 -1 -2 0 0 0 0] 69e9;
            nu              nu [0 0 0 0 0 0 0] 0.3;
            file            "$FOAM_CASE/constant/plasticStrainVsYieldStress";
            outOfBounds     clamp;
            S0              S0 [1 -1 -2 0 0 0 0] 1.1e6;
            b               1;
        }
    }
);
```

This is the local form, which is mesh dependent once the damage localises.
The tutorial `tutorials/solids/fracture/notchedBarLemaitre` adds

```text
            nonLocal        classic;
            characteristicLength characteristicLength [0 1 0 0 0 0 0] 6.325e-4;
```

and compares the reaction force with the Abaqus result of Whelan et al.

### Verification

`applications/test/LemaitreDamageMechanicalConstitutiveLaw` drives the law at
a material point and checks that:

- with the damage threshold out of reach, it is `linearElasticMisesPlastic`;
- for each `lemaitreType`, the returned state satisfies the discrete yield
  and damage equations, recomputed from the stress with `Y` in the
  triaxiality form of Whelan et al., equation 52;
- under uniaxial stress, the damage converges at first order in the step to
  the analytical solution of Doghri, to 0.3% with 200 steps;
- crack closure has no effect when every principal stress is tensile and
  reduces the damage when one is compressive;
- the cut-off keeps the damage at zero under compressive triaxiality;
- a point whose damage passes `Dc` fails and stays failed;
- under `logStrain`, on a path whose principal directions stay fixed, the
  Kirchhoff stress equals the small-strain law's at the logarithmic strain;
- in the non-local forms, with the input given: with it zero the stress is
  `linearElasticMisesPlastic`'s and the local damage follows the explicit
  update; with it held the state satisfies the discrete equations degraded
  by it; iterated to the local damage within each step, the converged outer
  iteration of a vanishing length, both forms give the local form to 2e-14;
  and the effective form fails a point once the input passes `Dc` and stops
  its local damage.

### References

- J. Lemaitre, A continuous damage mechanics model for ductile fracture,
  J. Eng. Mater. Technol. 107 (1985) 83-89.
- E.A. de Souza Neto, D. Peric and D.R.J. Owen, Computational Methods for
  Plasticity, Wiley, 2008.
- L. Malcher, F.M. Andrade Pires and J.M.A. Cesar de Sa, An assessment of
  isotropic constitutive models for ductile fracture under high and low
  stress triaxiality, Int. J. Plasticity 30-31 (2012) 81-115.
- I. Doghri, Mechanics of Deformable Solids, Springer, 2000.
- R.H.J. Peerlings, R. de Borst, W.A.M. Brekelmans and J.H.P. de Vree,
  Gradient enhanced damage for quasi-brittle materials, Int. J. Numer. Meth.
  Eng. 39 (1996) 3391-3403.
- A. Whelan, T. Tang, V. Pakrashi and P. Cardiff, A finite volume framework
  for damage and fracture prediction in wire drawing, Int. J. Numer. Meth.
  Eng. (2025), <https://doi.org/10.1002/nme.7640>.
