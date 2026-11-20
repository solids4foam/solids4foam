---
sort: 30
---

# linearElasticGTNPlastic

Small-strain Gurson-Tvergaard-Needleman (GTN) plasticity: a porous plastic
material whose void volume fraction grows with the plastic dilatation and
by nucleation. The runtime type is:

```text
linearElasticGTNPlastic
```

This is the GTN model of Whelan et al. (2025), in its local form or
regularised by a non-local void fraction. For finite strains, use it as the
sub-law of [`logStrain`](../logStrainMechanicalConstitutiveLaw/README.md),
which is the formulation of that work.

---

## User Guide

### What it computes

With `p = tr(sigma)/3` the mean stress, `q` the von Mises equivalent stress,
`sigmaY` the yield stress of the matrix and `f*` the effective void
fraction, the yield function is

```text
Phi = (q/sigmaY)^2 + 2 q1 f* cosh(3 q2 p/(2 sigmaY)) - 1 - q3 f*^2
```

with associated flow. The equivalent plastic strain of the matrix,
`epsilonPEq`, which drives the hardening, follows from the equivalence of
plastic work:

```text
(1 - f) sigmaY dEpsilonPEq = sigma : dEpsilonP
```

The void fraction `f` grows with the plastic dilatation, by nucleation, and
optionally in shear:

```text
df = (1 - f) tr(dEpsilonP) + A dEpsilonPEq + kw f omega dEpsilonQ
A  = fN/(sN sqrt(2 pi)) exp(-((epsilonPEq - epsilonN)/sN)^2/2)   if p > 0
```

where the nucleation is the strain-controlled form of Chu and Needleman,
under a tensile mean stress only, and the last term is the shear term of
Nahshon and Hutchinson, with `omega = 1 - (27 J3/(2 q^3))^2` and
`dEpsilonQ` the equivalent deviatoric plastic strain increment. The shear
term is written with `f`, as Nahshon and Hutchinson write it; the
implementation of Whelan used `f*`. Under a compressive mean stress the
voids close.

Without coalescence `f* = f`. With `coalescence yes`,

```text
f* = f                                  for f <= fC
f* = fC + (fu - fC)(f - fC)/(fF - fC)   for f >  fC
```

where `fu` is the effective void fraction at which the yield surface shrinks
to a point, `1/q1` for the usual `q3 = q1^2`.

When the void fraction would take `f*` past `0.99 fu`, the point fails: its
void fraction is held there, where, for `q3 = q1^2`, the material keeps
about 1% of its strength, and it stays failed. The failed void fraction is
reported when the case starts.

### Non-local form

The local void fraction softens the material without a length scale, so
once it localises the solution depends on the mesh. With `nonLocal classic`,
the yield function, the work equivalence and the stress are written instead
with the implicit non-local void fraction `fbar` of Whelan et al., the
solution of

```text
fbar - div(l^2 grad(fbar)) = f
```

where `l` is `characteristicLength`, with zero flux on the boundary. As for
the non-local forms of
[`linearElasticMisesPlasticLemaitreDamage`](../linearElasticMisesPlasticLemaitreDamageMechanicalConstitutiveLaw/README.md),
the solid model solves this equation once per outer iteration, after the
stress update, and the law reads its solution as the input
`voidFractionNonLocal`, lagged by one outer iteration. The local void
fraction is then updated explicitly, and is held at the failed void fraction
once it passes it. With `characteristicLength` zero, at a converged step,
the non-local form is the local one.

The field `voidFractionNonLocal` is written, read on a restart, and needs a
linear solver entry in `system/fvSolution`:

```text
voidFractionNonLocal
{
    solver          PCG;
    preconditioner  FDIC;
    tolerance       1e-09;
    relTol          0;
}
```

The equation is solved by the updated Lagrangian, total Lagrangian and
linear geometry total displacement solid models with the cell-centred stress
and the `implicitSegregated` algorithm; others refuse the law.

### Time integration

The return mapping and the void fraction update are solved together at each
point, implicitly. For a trial void fraction, the yield condition, the
normality of the flow and the work equivalence are solved for the
deviatoric, volumetric and matrix plastic strain increments by Newton's
method, with an analytical Jacobian and a backtracking line search; the void
fraction is then the root of the void growth equation, found by regula falsi
between the void fraction at the start of the step and the failed one, or
zero when the voids close. Whelan et al. instead lag the void fraction by
one outer iteration of the solid model; that scheme converges to the same
solution, but this one does not depend on how many times the law has been
evaluated, as a law in this framework must not.

### Model options

| Entry | Required | Default | Description |
| --- | --- | --- | --- |
| `rho` | yes | - | Density |
| `E`, `nu` | yes | - | Young's modulus and Poisson's ratio |
| `file` or `fileName` | yes | - | Matrix yield stress against `epsilonPEq` |
| `outOfBounds` | yes | - | `clamp`, `warn`, `error` or `repeat` |
| `q1`, `q2` | no | `1.5`, `1` | Tvergaard parameters |
| `q3` | no | `q1^2` | Tvergaard parameter, in `(0, q1^2]` |
| `f0` | yes | - | Initial void fraction |
| `fN` | no | `0` | Volume fraction of void-nucleating particles |
| `sN`, `epsilonN` | with `fN` | - | Deviation and mean nucleation strain |
| `kw` | no | `0` | Shear damage parameter |
| `coalescence` | no | `no` | Coalescence of voids |
| `fC`, `fF` | coalescence | - | Void fraction at its start and at failure |
| `nonLocal` | no | `none` | `none` or `classic` |
| `characteristicLength` | non-local | - | Length `l` of `fbar` |

The Newton loop stops when the update of the plastic strain increments,
relative to the largest strain, and the yield function both fall below
`NewtonLoopTol` (`1e-10`), or after `NewtonMaxIter` (`100`) iterations. The
coupled solution stops at `voidFractionTol` (`1e-12`) or
`voidFractionMaxIter` (`100`). Plane stress is not supported.

### State variables

| Variable | Description |
| --- | --- |
| `epsilonP` | Plastic strain (small-strain update only) |
| `epsilonPEq` | Equivalent plastic strain of the matrix |
| `sigmaY` | Yield stress of the matrix |
| `voidFraction` | Void fraction `f`; the local one in the non-local form |

Under `logStrain`, `epsilonP` stays zero: the plastic deformation is carried
by the elastic history of `logStrain`.

### Tangents

The scalar tangent is the elastic stiffness with its deviatoric part scaled
by the return, as for `linearElasticMisesPlastic`. It is stiffer than the
true tangent, which softens with the plastic dilatation, and serves only the
convergence of the solid model. No analytical fourth-order tangent is
provided; `fourthOrderFiniteDifference` is available. No volumetric split is
provided, so mixed displacement-pressure formulations refuse the law.

### Diagnostics

At the end of each step, with the law's debug switch on (the default), the
yielding points, the maximum plastic increment, the maximum void fraction and
the number of failed points are reported.

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
            type            linearElasticGTNPlastic;
            E               E [1 -1 -2 0 0 0 0] 69e9;
            nu              nu [0 0 0 0 0 0 0] 0.3;
            file            "$FOAM_CASE/constant/plasticStrainVsYieldStress";
            outOfBounds     clamp;
            q1              1.5;
            q2              1;
            q3              2.25;
            f0              0.002;
            fN              0.2;
            sN              0.08;
            epsilonN        0.15;
        }
    }
);
```

The tutorial `tutorials/solids/fracture/notchedBarGTN` runs it and compares
the reaction force with the Abaqus result of Whelan et al.

### Verification

`applications/test/GTNPlasticMechanicalConstitutiveLaw` drives the law at a
material point and checks that:

- with no voids, it is `linearElasticMisesPlastic`;
- on a path of tensile mean stress, the returned state satisfies the
  discrete yield condition, the normality of the flow, the work equivalence
  and the void growth with nucleation, recomputed from the stress, with
  linear and nonlinear hardening, with the shear term and with coalescence;
- under hydrostatic tension the voids grow, and under compression they
  close, the equations holding;
- a point whose void fraction passes the failed one fails, stays failed and
  keeps about 1% of its strength, while a step whose void growth equation
  has a root near the start of the step returns it, even where the equation
  changes sign again short of failure;
- the void fraction converges at first order in the step;
- under `logStrain`, on a path whose principal directions stay fixed, the
  Kirchhoff stress equals the small-strain law's at the logarithmic strain;
- in the non-local form, with the input zero the stress is
  `linearElasticMisesPlastic`'s while the local voids nucleate, and with the
  input iterated to the local void fraction within each step, the converged
  outer iteration of a vanishing length, the local form is recovered.

### References

- A.L. Gurson, Continuum theory of ductile rupture by void nucleation and
  growth: Part I, J. Eng. Mater. Technol. 99 (1977) 2-15.
- V. Tvergaard and A. Needleman, Analysis of the cup-cone fracture in a
  round tensile bar, Acta Metall. 32 (1984) 157-169.
- C.C. Chu and A. Needleman, Void nucleation effects in biaxially stretched
  sheets, J. Eng. Mater. Technol. 102 (1980) 249-256.
- K. Nahshon and J.W. Hutchinson, Modification of the Gurson model for
  shear failure, Eur. J. Mech. A Solids 27 (2008) 1-17.
- N. Aravas, On the numerical integration of a class of pressure-dependent
  plasticity models, Int. J. Numer. Meth. Eng. 24 (1987) 1395-1416.
- A. Whelan, T. Tang, V. Pakrashi and P. Cardiff, A finite volume framework
  for damage and fracture prediction in wire drawing, Int. J. Numer. Meth.
  Eng. (2025), <https://doi.org/10.1002/nme.7640>.
