---
sort: 20
---

# OgdenElastic

Three-term Ogden hyperelasticity, written on the principal stretches, with a
bulk-modulus penalty on volume change. The runtime type is:

```text
OgdenElastic
```

The dictionary layout common to all laws is described in the
[material models page](../../../materialModels/README.md).

---

## User Guide

### What it computes

At each integration point, from `F` and `J = det(F)`, the right Cauchy-Green
tensor is decomposed into its eigenvalues `c_i` and eigenvectors `N_i`:

```text
C        = F.T() & F
lambda_i = max(sqrt(c_i), VSMALL)
p_i      = mu1*lambda_i^alpha1 + mu2*lambda_i^alpha2 + mu3*lambda_i^alpha3
s        = sum_i p_i*N_i*N_i
sigmaHyd = 0.5*K*(J^2 - 1)
sigma    = (dev(s) + sigmaHyd*I + symm(F & sigma0 & F.T()))/J
```

The stretches are the full principal stretches, not isochoric ones, and `s`
is formed in the eigenvectors of `C`, the principal directions of the
reference configuration. `sigma0` is an optional initial stress (see State
variables), pushed forward with `F`.

The law is finite-strain only; a linear geometry solid model stops with an
error. It holds no history.

### Model options

| Entry | Required | Default | Description |
| --- | --- | --- | --- |
| `rho` | yes | - | Density, `[1 -3 0 0 0 0 0]` |
| `mu1` | yes | - | First modulus, `[1 -1 -2 0 0 0 0]` |
| `mu2` | yes | - | Second modulus, `[1 -1 -2 0 0 0 0]` |
| `mu3` | yes | - | Third modulus, `[1 -1 -2 0 0 0 0]` |
| `alpha1` | yes | - | First exponent, dimensionless |
| `alpha2` | yes | - | Second exponent, dimensionless |
| `alpha3` | yes | - | Third exponent, dimensionless |
| `K` | yes | - | Bulk modulus, `[1 -1 -2 0 0 0 0]` |

The entries are read as dimensioned scalars, but their dimensions and values
are not checked, and there is no check on `J`. The parameters are printed at
construction. The top-level `planeStress` entry is not used by this law:
there is no plane-stress reduction.

### Tangents and volumetric split

- Scalar tangent: `(4/3)*(mu1 + mu2 + mu3) + K`, returned for both the
  `scalar` and the `scalarDeviatoric` requests.
- Fourth-order tangent: `fourthOrderFiniteDifference` only. Asking for the
  analytical `fourthOrder` tangent is a fatal error.
- Volumetric split: no. The energy is not written on an isochoric measure, so
  the law cannot separate its isochoric stress from its volumetric response.
  A mixed displacement-pressure formulation refuses it, naming the material,
  and `solvePressureEqn` is refused at start-up.

### State variables

| Name | Type | Role | Description |
| --- | --- | --- | --- |
| `sigma0` | symmTensor | prescribed | Initial stress, default zero |

`sigma0` is read from a `volSymmTensorField` named `sigma0` in the current
time directory, or in `0`, or from a field of that name registered by another
model; otherwise it is zero. It is the user's input rather than history, so
the law writes no restart files and no state fields.

### Example

The equivalent Ogden parameters given, commented out, in the `cylinderCrush`
hyperelasticity tutorial:

```text
planeStress     no;

mechanical
(
    rubber
    {
        type            OgdenElastic;
        rho             rho [1 -3 0 0 0 0 0] 1000;
        K               K [1 -1 -2 0 0 0 0] 1.410e+9;
        mu1             mu1 [1 -1 -2 0 0 0 0] 0.746e+6;
        mu2             mu2 [1 -1 -2 0 0 0 0] -0.306e+6;
        mu3             mu3 [1 -1 -2 0 0 0 0] 6.609e-5;
        alpha1          alpha1 [0 0 0 0 0 0 0] 1.748;
        alpha2          alpha2 [0 0 0 0 0 0 0] -1.656;
        alpha3          alpha3 [0 0 0 0 0 0 0] 7.671;
    }
);
```

The tutorial itself runs `MooneyRivlinElastic` with `solvePressureEqn yes`,
which this law does not accept.

### Differences from the legacy law

- `solvePressureEqn` and `pressureSmoothingScaleFactor` are now refused: the
  law has no volumetric split, so the solid model cannot smooth its
  hydrostatic stress. The legacy law smoothed it itself.
- `regionName` is no longer read.
- The law is evaluated at whatever integration points the solid model uses;
  the legacy law was cell-centred only, and its face form aborted.
- The legacy fallback to a Hookean stress when the solid model enforced
  linearity is gone.
- A fourth-order tangent is available by finite differences.
- `sigma0` is still optional and still pushed forward with `F`; it is now
  read as a prescribed state field rather than by the base class.

---

## Tutorials

No tutorial selects this law. The
[cylinderCrush](../../../../../tutorials/solids/hyperelasticity/cylinderCrush/README.md)
and
[rubberSealing](../../../../../tutorials/solids/hyperelasticity/rubberSealing/README.md)
hyperelasticity tutorials keep Ogden parameters as a commented alternative;
`rubberSealing` notes that the law supports neither plane stress nor the mixed
formulation.
