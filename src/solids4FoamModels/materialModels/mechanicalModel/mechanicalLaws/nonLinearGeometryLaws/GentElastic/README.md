---
sort: 21
---

# GentElastic

This page documents the Gent-Flory elastic law for microgel behaviour. The
runtime type is:

```text
GentElastic
```

---

## User Guide

### What it computes

The law updates the deformation gradient and evaluates the following stress at
cell centres or faces:

```text
J     = det(F)
B     = symm(F & F.T())
C     = symm(F.T() & F)
omega = Vs/Na
Nv    = N/omega

sigma = kb*temperature*(Nv/J)
      * ((ilim/(ilim - tr(B) + 3))*C - I)
```

Here `Na` is fixed at `6.022141e23` and dimensionless, while `kb` is fixed at
`1.38064852e-23` with dimensions `[1 2 -2 -1 0 0 0]`. The nonlinear stress
does not use `E`, `nu`, `mu` or `K`; those quantities are used by `updateF()`
when the solid model enforces a linear response and by `impK()`.

The `volSymmTensorField` and `surfaceSymmTensorField` `correct()` overloads are
implemented. The inherited `pointSymmTensorField` overload aborts with
`notImplemented`. On OpenFOAM.com and OpenFOAM.org, the inherited
`CompactListList` overload also exists and aborts with `notImplemented`; that
interface is excluded from foam-extend by `#ifndef FOAMEXTEND`.

### Model options

| Entry | Required | Description |
| --- | --- | --- |
| `E` | yes | Young's modulus, `[1 -1 -2 0 0 0 0]` |
| `nu` | yes | Poisson's ratio, dimensionless |
| `temperature` | yes | Temperature, `[0 0 0 1 0 0 0]` |
| `ilim` | yes | Limiting invariant parameter, dimensionless |
| `Vs` | yes | Volume divided by `Na` to form `omega`, `[0 3 0 0 0 0 0]` |
| `N` | yes | Crosslink parameter divided by `omega`, dimensionless |
| `rho` | yes | Density, `[1 -3 0 0 0 0 0]` |
| `solvePressureEqn` | no | Base switch; default `no` |
| `pressureSmoothingScaleFactor` | no | Base scale factor; default `100` |
| `regionName` | no | Base mesh region; otherwise detected automatically |

The top-level `planeStress` switch is also required. It changes the `K` used
for linear enforcement and `impK()`, but it does not change the nonlinear
stress expression. Although the base constructor reads `solvePressureEqn` and
`pressureSmoothingScaleFactor`, this law never calls `updateSigmaHyd()`, so
neither option changes its computed stress.

No range or positivity checks are applied to the material entries. In
particular, the denominator `ilim - tr(B) + 3` is not guarded against zero.

### Recommended dictionary setup

Illustrative setup for a water-swollen microgel:

```text
planeStress     no;

mechanical
(
    waterSwollenMicrogel
    {
        type            GentElastic;
        rho             rho [1 -3 0 0 0 0 0] 1000;
        E               E [1 -1 -2 0 0 0 0] 3e4;
        nu              nu [0 0 0 0 0 0 0] 0.49;
        temperature     temperature [0 0 0 1 0 0 0] 298;
        ilim            ilim [0 0 0 0 0 0 0] 100;
        Vs              Vs [0 3 0 0 0 0 0] 1.8e-5;
        N               N [0 0 0 0 0 0 0] 1e-4;

        // Optional base-class entries
        // regionName                    region0;
        // solvePressureEqn              no;
        // pressureSmoothingScaleFactor  100;
    }
);
```

### Field glossary

- `F_`, `Ff_`: total deformation gradients at cell centres and faces. They
  are read if present and written for restart.
- `relF_`, `relFf_`: relative deformation gradients maintained by `updateF()`.
- `mu_<material>`, `K_<material>`: cell fields populated from the uniform
  derived moduli when the volume correction is used.
- `muf_<material>`, `Kf_<material>`: corresponding face fields populated by
  the surface correction.
- `impK`: implicit stiffness field returned to the solid model, with value
  `(4/3)*mu + K` unless an `impK` field is read from the current time.
- `sigma`: stress field supplied by the solid model and overwritten by
  `correct()`.

---

## Developer Notes

### Class role

`GentElastic` derives directly from `mechanicalLaw` and is registered in the
`nonLinGeomMechLaw` runtime selection table. It stores ten constant
`dimensionedScalar` members: `E_`, `nu_`, `mu_`, `temperature_`, `ilim_`,
`K_`, `Na_`, `kb_`, `omega_` and `N_`.

The class itself has no fork guards. The inherited quadrature-point interface
is guarded by `#ifndef FOAMEXTEND`. `GentElastic.C` appears in both
`Make/files.openfoam` and `Make/files.foamextend`.

### Construction

The base constructor first reads or inserts `solvePressureEqn` and
`pressureSmoothingScaleFactor`, then reads `regionName` or detects `solid` or
`region0`. Failure to find a base mesh region is fatal.

The derived constructor then performs these operations in order:

1. Reads `E` and `nu`, and derives `mu = E/(2*(1 + nu))`.
2. Reads `temperature` and `ilim`.
3. Reads the top-level `planeStress` switch and derives `K`:

   ```text
   K = nu*E/((1 + nu)*(1 - nu))     + (2/3)*mu  // plane stress
   K = nu*E/((1 + nu)*(1 - 2*nu))   + (2/3)*mu  // otherwise
   ```

4. Sets the fixed values of `Na` and `kb`.
5. Reads `Vs` and calculates `omega = Vs/Na`.
6. Reads `N` and replaces it internally by `N/omega`.

Missing or dimensionally inconsistent required entries cause dictionary or
dimension errors. The constructor contains no explicit material validation
and emits no law-specific fatal error or warning. If `solvePressureEqn` is
enabled, the base constructor only prints an informational message.

### Key methods

- `impK()` returns a cell field equal to `(4/3)*mu + K`. It is the
  diffusivity of the solid model's implicit Laplacian term and affects the
  outer-iteration convergence rate rather than the converged answer.
- Both implemented `correct()` methods call `updateF(sigma, mu, K)` first. If
  it returns `true`, the linearised stress supplied by `updateF()` is kept.
  Otherwise they evaluate the Gent-Flory expression above using `F` or `Ff`.
  `updateF()` aborts for non-incremental updated-Lagrangian use or an unknown
  nonlinear-geometry type, and warns when linearity is enforced for stability.
- `materialTangent()` is not overridden. The inherited implementation aborts
  with `notImplemented`, so solid models that request a material tangent
  cannot use this law.
- `setRestart()` sets both `F` and `Ff` to `AUTO_WRITE`, constructing either
  field if it has not already been used.

### Extension points

A related law can copy this class, replace both `correct()` expressions, and
update the `dimensionedScalar` inputs and `impK()` linearisation. It must use a
new `TypeName`, register in the `nonLinGeomMechLaw` table, and be added to both
build lists when it supports all forks. Implement the point, quadrature-point
and material-tangent interfaces if the target solid models require them.

The source is at
[GentElastic.C](https://github.com/solids4foam/solids4foam/blob/master/src/solids4FoamModels/materialModels/mechanicalModel/mechanicalLaws/nonLinearGeometryLaws/GentElastic/GentElastic.C).

---

## Tutorials

No tutorial currently selects the `GentElastic` runtime type.
