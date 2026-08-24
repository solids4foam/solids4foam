---
sort: 24
---

# HolzapfelGasserOgdenElastic

This page documents the Holzapfel-Gasser-Ogden (HGO) anisotropic hyperelastic
law for fibre-reinforced soft tissue, in particular arterial wall. The runtime
type is:

```text
HolzapfelGasserOgdenElastic
```

The law is **incompressible only**: the bulk modulus is always set to `GREAT`
and the hydrostatic stress is taken from a solved pressure field. It is
therefore intended for `coupledPressureDisplacementSolid`.

The constitutive model is due to G. A. Holzapfel, T. C. Gasser and
R. W. Ogden, _A new constitutive framework for arterial wall mechanics and a
comparative study of material models_, Journal of Elasticity, 61 (2000) 1-48.

---

## User Guide

### Constitutive relation

The response is an isotropic neo-Hookean ground matrix plus two exponential
fibre families, symmetric about the circumferential direction. With
incompressibility (`J = 1`) enforced through the pressure field:

```text
b     = symm(F & F.T)
s     = mu*(b - I)

N4    =  cos(fibreAngle)*Ec + sin(fibreAngle)*Ea
N6    = -cos(fibreAngle)*Ec + sin(fibreAngle)*Ea
n4    = F & N4
n6    = F & N6

I4    = C && symm(N4*N4)
I6    = C && symm(N6*N6)

sigma = -p*I + s
      + 2*k1*(I4 - 1)*exp(k2*sqr(I4 - 1))*symm(n4*n4)
      + 2*k1*(I6 - 1)*exp(k2*sqr(I6 - 1))*symm(n6*n6)
```

This corresponds to the strain energy

```text
Psi = 0.5*mu*(I1 - 3)
    + (k1/(2*k2))*sum_i (exp(k2*sqr(I_i - 1)) - 1)
```

with `i` running over the two fibre families.

```warning
The fibre terms are active in compression as well as tension. The standard HGO
model switches them off when `I4 < 1`; this implementation does not, so the
fibres carry compressive load.
```

### Local material directions

The fibre directions are built from a local cylindrical triad that must already
exist as fields:

- `Ec`, `Ea`, `Er` as `volVectorField`s (circumferential, axial, radial);
- `Ecf`, `Eaf`, `Erf` as `surfaceVectorField`s.

All six are read `MUST_READ`, so the case must provide them before the solver
starts. In the `ratCarotid` tutorial they are produced by a case-local utility,
`calcLocCoordinates`, which the `Allrun` script compiles and runs first. Note
that `Er` is multiplied by zero when forming the fibre vectors, so only `Ec` and
`Ea` affect the result.

### Model options

Required entries:

| Entry | Type | Description |
| --- | --- | --- |
| `rho` | `dimensionedScalar` | Density, `[1 -3 0 0 0 0 0]` |
| `fibreAngle` | `dimensionedScalar` | Fibre angle in **degrees** from `Ec` |
| `k1` | `dimensionedScalar` | Fibre stress scale, `[1 -1 -2 0 0 0 0]` |
| `k2` | `dimensionedScalar` | Dimensionless fibre exponent |

The shear modulus is given either as `E` together with `nu`, or as `mu` alone;
supplying both forms, or neither, is a fatal error. When `E` and `nu` are used,
`mu = E/(2*(1 + nu))` and `nu` must be `0.5` — anything smaller aborts with
"This is incompressible solid model (nu=0.5)".

Optional entries:

| Entry | Default | Description |
| --- | --- | --- |
| `impKcoeff` | `1.0` | Scales the implicit stiffness `impK` |

`fibreAngle` is converted from degrees to radians in the constructor, so give
it in degrees.

### Compatibility

Both `correct(volSymmTensorField&)` and `correct(surfaceSymmTensorField&)` are
implemented, and both look up a pressure field — `p` at cells and `pf` at
faces. Only `coupledPressureDisplacementSolid` provides these, so that is
effectively the only compatible solid model. Total and updated Lagrangian
kinematics are both handled, because the deformation gradient is updated
through the base-class `updateF()`.

The law is registered for both OpenFOAM and foam-extend builds. The
`ratCarotid` tutorial is foam-extend-only for reasons unrelated to the law
itself.

### Recommended dictionary setup

```text
mechanical
(
    arteryWall
    {
        type            HolzapfelGasserOgdenElastic;
        rho             rho [1 -3 0 0 0 0 0] 1200;

        E               E [1 -1 -2 0 0 0 0] 132.69e3;
        nu              nu [0 0 0 0 0 0 0] 0.5;
        // or, equivalently:
        // mu           mu [1 -1 -2 0 0 0 0] 44.23e3;

        impKcoeff       1;

        // Fibre angle measured from the circumferential direction
        fibreAngle      fibreAngle [0 0 0 0 0 0 0] 39.76;
        k1              k1 [1 -1 -2 0 0 0 0] 0.206e3;
        k2              k2 [0 0 0 0 0 0 0] 1.465;
    }
);
```

`constant/solidProperties` must select `coupledPressureDisplacementSolid`, and
the start time directory must contain `Ec`, `Ea`, `Er`, `Ecf`, `Eaf` and `Erf`.

### Field glossary

- `Ec`, `Ea`, `Er`: local circumferential, axial and radial unit vectors at
  cell centres; `Ecf`, `Eaf`, `Erf` are the face equivalents. All are written
  with `AUTO_WRITE`.
- `muEff`: effective shear modulus used by `impK`, updated once per time step.
- `p`, `pf`: pressure at cells and faces, owned by the solid model.

---

## Developer Notes

### Class role

`HolzapfelGasserOgdenElastic` derives from `mechanicalLaw`. Unlike most laws it
overrides `rho()`, returning the locally read `rho_` as a uniform field with
`calculated` patch types.

### Construction

The constructor reads `rho`, `fibreAngle`, `k1`, `k2` and the six direction
fields, then resolves `mu_` from either `E`/`nu` or `mu` and pins `K_` to
`GREAT`. It converts `fibreAngle_` to radians, and finally calls
`calcInitialShearModulus()` followed by `calcEffectiveShearModulus()` to
populate `muEff_`.

### `correct()` walkthrough

`correct(volSymmTensorField&)` shadows the member `K_` with a **local**
zero-valued `dimensionedScalar` of the same name before calling `updateF()`:

```text
dimensionedScalar K_("K", mu_.dimensions(), 0);
if (updateF(sigma, mu_, K_)) { return; }
```

That local `K_` is what the enforced-linear fallback inside `updateF()` uses,
so if the solver ever enforces linearised elasticity the fallback behaves as a
zero-bulk-modulus material. The same shadowing appears in
`correct(surfaceSymmTensorField&)`.

After `updateF()`, the isotropic part is `mu*(b - I)` with `b = symm(F & F.T)`,
the pressure term is subtracted, and the two fibre contributions are added. The
Jacobian is never computed: the code assumes `J == 1`.

The face variant delegates to the private `correctF()`, which is the same
algebra written for `surface*` field types using `Ecf`, `Eaf` and `Erf`.

### Implicit stiffness `impK()`

```text
impK = impKcoeff*muEff
```

`muEff` is a field, not a constant. `calcEffectiveShearModulus()` perturbs the
displacement gradient by `0.001` in each of the three shear planes of the
principal deviatoric stress frame, recomputes the deviatoric Cauchy stress via
`calcDevCauchy()`, and stores the largest secant shear stiffness per cell. It
evaluates both the old-time and current deformation gradient and keeps the
maximum of the two. `updateTotalFields()` calls it at the end of each time step,
so `impK` tracks the stiffening fibre response.

`calcInitialShearModulus()` performs the same three shear perturbations once at
construction, starting from `F = I` and using the direction vectors of cell 0
only, and sets `mu_` to the largest result. This overwrites the value derived
from `E` and `nu`.

### Extension points

- Add the `I4 < 1` / `I6 < 1` switch to disable fibres in compression.
- Introduce the Gasser dispersion parameter `kappa`, which the full HGO model
  uses to blend `I1` and `I4`; this implementation has no dispersion.
- Generalise the two fibre families to an arbitrary list, or read fibre
  directions directly rather than deriving them from `Ec` and `Ea`.

---

## Tutorials

Cases that select `HolzapfelGasserOgdenElastic`:

- `solids/hyperelasticity/ratCarotid`
