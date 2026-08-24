---
sort: 26
---

# viscoNeoHookeanElastic

This law combines a finite-strain nonlinear deviatoric response with a
generalised Maxwell relaxation model and an elastic volumetric response. The
runtime type is:

```text
viscoNeoHookeanElastic
```

---

## User Guide

### What it computes

The law first forms the total deformation gradient, its determinant and its
volume-preserving part:

```text
J    = det(F)
Fbar = J^(-1/3)*F
C    = F^T*F
lambdaBar_a = J^(-1/3)*sqrt(eigenvalue_a(C))
```

For `m = 0, 1, 2`, the three supplied modulus/exponent pairs define a diagonal
tensor `D` and the current deviatoric response `S`:

```text
D_aa = sum_m mu_m*lambdaBar_a^(alpha_m - 1)
S    = dev(2*Fbar*D*Fbar^T)
s    = inv(Fbar)*S*inv(Fbar)^T
```

For each Maxwell branch `i`, the internal variables are updated from their
old-time values using the current time-step `deltaT`:

```text
H_i = exp(-deltaT/tau_i)*h_i.oldTime()
    - exp(-deltaT/(2*tau_i))*s.oldTime()

h_i = H_i + exp(-deltaT/(2*tau_i))*s
```

The relative moduli, relaxation factor, pressure and stress assigned by the
code are:

```text
E0       = EInfinity + sum_i E_i
gammaInf = EInfinity/E0
gamma_i  = E_i/E0
g        = gammaInf + sum_i gamma_i*exp(-deltaT/(2*tau_i))
p        = k/(beta*J)*(1 - J^(-beta))

sigma = J*p*I + (1 + g)*S
      + sum_i gamma_i*dev(Fbar*H_i*Fbar^T)
```

The class implements both `correct(volSymmTensorField&)` and
`correct(surfaceSymmTensorField&)`. The inherited point-centred overload aborts
with `notImplemented`. On OpenFOAM.com and OpenFOAM.org, the inherited
`CompactListList` quadrature-point overload also aborts with `notImplemented`;
that overload is not compiled for foam-extend.

### Model options

| Entry | Required | Description |
| --- | --- | --- |
| `EInfinity` | yes | Long-term modulus, `[1 -1 -2 0 0 0 0]` |
| `E` | yes | Maxwell moduli; scalar list used to form relative weights |
| `relaxationTimes` | yes | Positive scalar list, one value per `E` value |
| `nu` | yes | Long-term Poisson's ratio, dimensionless |
| `mu0` | yes | First nonlinear modulus, `[1 -1 -2 0 0 0 0]` |
| `mu1` | yes | Second nonlinear modulus, `[1 -1 -2 0 0 0 0]` |
| `mu2` | yes | Third nonlinear modulus, `[1 -1 -2 0 0 0 0]` |
| `alpha` | yes | Dimensionless exponent list; entries 0 through 2 are used |
| `beta` | yes | Dimensionless volumetric exponent |
| `k` | yes | Bulk modulus used by the stress law, `[1 -1 -2 0 0 0 0]` |
| `rho` | yes | Density, `[1 -3 0 0 0 0 0]` |
| `solvePressureEqn` | no | Base switch, default `no`; unused by `correct()` |
| `pressureSmoothingScaleFactor` | no | Default `100`; unused here |
| `regionName` | no | Base mesh region when it cannot be detected |

The `E` and `relaxationTimes` lists must have the same non-zero length. Every
entry of each list must be positive. The values in `E` are plain scalars, so
they must use the same numerical modulus unit as `EInfinity`.

The implementation indexes the first three `alpha` values without checking the
list length. Supply at least three values; later values are not used. The
constructor requires `k`, warns that it is read directly, and uses it instead
of the bulk modulus derived from `EInfinity` and `nu`.

`planeStress` is read from the enclosing `mechanicalProperties` dictionary. It
changes the derived `lambda`, but the supplied `k` still replaces the derived
bulk modulus. The pressure-equation entries are constructed by the base class,
but this law does not call `updateSigmaHyd()`.

### Recommended dictionary setup

The following values define an illustrative three-branch silicone-like
material:

```text
planeStress     no;

mechanical
(
    siliconeLike
    {
        type            viscoNeoHookeanElastic;
        rho             rho [1 -3 0 0 0 0 0] 1100;
        EInfinity       EInfinity [1 -1 -2 0 0 0 0] 1.0e6;
        E               (5.0e5 2.0e5 1.0e5);
        relaxationTimes (0.1 1 10);
        nu              nu [0 0 0 0 0 0 0] 0.49;
        mu0             mu0 [1 -1 -2 0 0 0 0] 3.0e5;
        mu1             mu1 [1 -1 -2 0 0 0 0] 1.0e5;
        mu2             mu2 [1 -1 -2 0 0 0 0] 5.0e4;
        alpha           (2 4 6);
        beta            beta [0 0 0 0 0 0 0] 2;
        k               k [1 -1 -2 0 0 0 0] 5.0e7;

        // Optional base-class entries
        // solvePressureEqn             no;
        // pressureSmoothingScaleFactor 100;
        // regionName                   region0;
    }
);
```

### Field glossary

- `F_`, `Ff_`: total deformation gradients at cell centres and faces, created
  by the base class and read on restart when present.
- `relF_`, `relFf_`: relative deformation gradients maintained by the base
  class while updating `F_` and `Ff_`.
- `s`, `sf`: pull-backs of the current deviatoric response at cells and faces;
  their old-time values enter the relaxation update.
- `h<i>`, `hf<i>`: internal stress variables for Maxwell branch `i`; their
  old-time values are stored.
- `H<i>`, `Hf<i>`: intermediate branch-history fields calculated during each
  stress correction.
- `transformH<i>`, `transformHf<i>`: transformed branch-history fields used in
  the final stress sum.
- `transformNeeded`, `transformNeededf`: diagonal nonlinear response tensors
  assembled from the modified principal stretches.
- `transformFbar`, `transformFbarf`: transformed nonlinear response tensors
  whose deviatoric parts give `S`.
- `impK`: implicit stiffness supplied to the solid model's Laplacian term.

---

## Developer Notes

### Class role

`viscoNeoHookeanElastic` derives directly from `mechanicalLaw` and is
registered in the nonlinear-geometry mechanical-law selection table. It stores
the long-term modulus and Poisson's ratio, Maxwell moduli and times, relative
weights, three nonlinear moduli, the exponent list, the volumetric exponent,
the bulk modulus, and cell- and face-centred history fields.

The class is listed in both `Make/files.openfoam` and `Make/files.foamextend`.
The header includes `surfaceFields.H` only under `OPENFOAM_NOT_EXTEND`. The
source uses `realEigenValues` under `OPENFOAM_COM` and fork-specific field
accessors under `OPENFOAM_NOT_EXTEND`; the constitutive algorithm is otherwise
present for all three forks.

### Construction

The base constructor first reads or inserts `solvePressureEqn` and
`pressureSmoothingScaleFactor`, then selects `regionName`. The derived
constructor reads `EInfinity`, `E`, `relaxationTimes`, `nu`, `mu0`, `mu1`,
`mu2`, `alpha`, `beta` and `k`, in that order, and creates the cell and face
history fields.

It then performs the following operations:

- a `nu` outside `[-1, 0.5]` causes a fatal error;
- `nu == 0.5` prints an incompressibility message, sets `lambda` and `k` to
  `GREAT`, and then encounters the fatal independent-`k` check;
- otherwise it derives `lambda` and a bulk modulus, warns that `k` is read
  directly, and replaces the derived bulk modulus with the supplied `k`;
- unequal `E` and `relaxationTimes` list lengths cause a fatal error;
- it computes `gammaInf` and each `gamma_i` from the total modulus;
- a relaxation time below `SMALL` causes a fatal error;
- an `E` value below `SMALL` causes a fatal error;
- it prints the relative moduli, allocates one set of history fields per
  Maxwell branch, and stores the old-time `h`, `hf`, `s` and `sf` fields.

The constructor does not validate `EInfinity`, `beta`, the total modulus, or
the length of `alpha`.

### Key methods

- `impK()`: returns `2*mu` when `nu == 0.5` or a cell field named `p` exists;
  otherwise it returns `2*mu + 4*k/3`. It is the diffusivity of the solid
  model's Laplacian term and affects the outer-iteration convergence rate
  rather than the converged answer.
- `K()`: returns a uniform cell field containing the stored `k` value.
- `correct()`: updates `F` or `Ff`, evaluates the nonlinear deviatoric
  response from the three modulus/exponent pairs, advances the Maxwell
  histories, and assigns the stress shown above. If the solid model enforces
  linear material behavior, the base-class `updateF()` assigns an elastic
  stress using `mu` and `k` and the method returns immediately.
- `materialTangent()`: is not overridden; the inherited implementation aborts
  with `notImplemented`.
- `setRestart()`: sets `F` and `Ff` to `AUTO_WRITE`.

### Extension points

A related law can copy this class and replace the nonlinear response and
history updates in the two `correct()` methods. New branch properties should
be read and validated before allocating matching cell and face history lists.
Implement the point and quadrature-point overloads, and `materialTangent()`, if
the new law must support vertex-centred or tangent-based solid models. Keep the
source in both build lists when all three OpenFOAM forks are supported.

The source is at
[viscoNeoHookeanElastic.C](https://github.com/solids4foam/solids4foam/blob/master/src/solids4FoamModels/materialModels/mechanicalModel/mechanicalLaws/nonLinearGeometryLaws/viscoNeoHookeanElastic/viscoNeoHookeanElastic.C).

---

## Tutorials

No tutorial currently uses `viscoNeoHookeanElastic`.
