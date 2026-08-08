---
sort: 18
---

# MooneyRivlinElastic

This page documents the three-parameter Mooney-Rivlin hyperelastic law. The
runtime type is:

```text
MooneyRivlinElastic
```

The model extends the neo-Hookean form with a dependence on the second
invariant, which lets it reproduce rubber behaviour over a wider range of
deformation modes. It supports heterogeneous material coefficients and an
initial stress field.

References:

- M. Mooney (1940), _A theory of large elastic deformation_, Journal of
  Applied Physics 11, 582-592,
  [10.1063/1.1712836](https://doi.org/10.1063/1.1712836).
- R. S. Rivlin and D. W. Saunders (1951), _Large elastic deformations of
  isotropic materials VII. Experiments on the deformation of rubber_,
  Philosophical Transactions of the Royal Society of London A 243, 251-288,
  [10.1098/rsta.1951.0004](https://doi.org/10.1098/rsta.1951.0004).

---

## User Guide

### Strain energy function

The isochoric part uses the three-term Rivlin-Saunders form:

```text
Psi_iso = c10*(I1Bar - 3) + c01*(I2Bar - 3)
        + c11*(I1Bar - 3)*(I2Bar - 3)
```

where the isochoric left Cauchy-Green tensor and its invariants are

```text
bBar  = J^(-2/3)*(F & F.T()),   J = det(F)
I1Bar = tr(bBar)
I2Bar = 0.5*(I1Bar^2 - tr(bBar & bBar))
```

The volumetric part is the same as in `neoHookeanElastic`,
`U(J) = 0.5*K*(0.5*(J^2 - 1) - ln(J))`. The Cauchy stress evaluated by the
code is

```text
J*sigma = 0.5*K*(J^2 - 1)*I
        + dev( 2*(c10 + c11*(I2Bar - 3))*bBar
             - 2*(c01 + c11*(I1Bar - 3))*inv(bBar) )
        + F & sigma0 & F.T()
```

Setting `c11` to zero recovers the classical two-parameter Mooney-Rivlin
model, and additionally setting `c01` to zero recovers a neo-Hookean model
with `mu = 2*c10`.

The law is registered in the nonlinear-geometry table, so it works with both
total-Lagrangian and updated-Lagrangian solid models. It is not available to
linear-geometry solid models.

### Model options

| Entry | Requirement | Description |
| --- | --- | --- |
| `c10` | required | First coefficient, `[1 -1 -2 0 0 0 0]` |
| `c01` | required | Second coefficient, `[1 -1 -2 0 0 0 0]` |
| `c11` | required | Coupled coefficient, `[1 -1 -2 0 0 0 0]` |
| `rho` | required | Density, `[1 -3 0 0 0 0 0]` |

All three coefficients are read with `lookup()`, so all three must be present;
give `c11 c11 [1 -1 -2 0 0 0 0] 0;` if the coupled term is not wanted.

The bulk modulus is set by **either** `K` **or** `nu`; if neither is present
the constructor aborts. If both are present, `K` wins.

| Entry | Requirement | Description |
| --- | --- | --- |
| `K` | one of two | Bulk modulus, `[1 -1 -2 0 0 0 0]` |
| `nu` | one of two | Poisson's ratio, used if `K` is absent |

When `nu` is supplied, `K` is derived from the linearised Young's modulus,
`E = 6*(c10 + c01)` and `K = E/(3*(1 - 2*nu))`. The linearised shear modulus
is always `mu = 2*(c10 + c01)`.

Optional base-class entries:

| Entry | Default | Description |
| --- | --- | --- |
| `solvePressureEqn` | `false` | Smooth `sigmaHyd` with a Laplacian |
| `pressureSmoothingScaleFactor` | `100` | Scale factor for that smoothing |
| `regionName` | auto | Base mesh region name |

### Heterogeneous coefficients and pre-stress

`c10`, `c01` and `c11` are stored as `volScalarField`s created with
`READ_IF_PRESENT` and `AUTO_WRITE`. Placing `c10`, `c01` or `c11` fields in
the start-time directory therefore overrides the uniform dictionary value with
a spatially varying one, and the fields are written at every output time.

The base-class `sigma0` initial-stress field is pushed forward with `F` and
added to the Cauchy stress. Note that the uniform values in the dictionary
are still required even when field files are supplied.

### Recommended dictionary setup

Minimal example for `constant/mechanicalProperties`:

```text
planeStress     no;

mechanical
(
    rubber
    {
        type            MooneyRivlinElastic;
        rho             rho [1 -3 0 0 0 0 0] 1000;

        K               K [1 -1 -2 0 0 0 0] 1.410e+09;
        c10             c10 [1 -1 -2 0 0 0 0] 0.293e+06;
        c01             c01 [1 -1 -2 0 0 0 0] 0.177e+06;
        c11             c11 [1 -1 -2 0 0 0 0] 0;

        solvePressureEqn             yes;
        pressureSmoothingScaleFactor 100;
    }
);
```

### Field glossary

- `c10`, `c01`, `c11`: material coefficient fields, written at output times.
- `F`, `Ff`: total deformation gradient at cell centres and faces.
- `sigma0`: initial stress field, pushed forward and added to `sigma`.
- `sigmaHyd`: hydrostatic Cauchy stress, `-p`.
- `sigma`: Cauchy stress tensor returned to the solid model.
- `impK`: implicit stiffness used as the Laplacian diffusivity.

---

## Developer Notes

### Class role

`MooneyRivlinElastic` inherits from `mechanicalLaw` and is added to the
`nonLinGeomMechLaw` runtime selection table. Unlike `neoHookeanElastic`, it
uses the base class's field-valued `mu()`, `K()`, `muf()` and `Kf()`
accessors, which is what allows spatially varying coefficients.

### Construction

The three coefficient fields are constructed with `READ_IF_PRESENT` from the
dictionary values. The body then sets `mu()` from `2*(c10 + c01)`, resolves
`K()` from `K` or `nu`, prints the coefficient ranges, and finally
interpolates `mu()` and `K()` onto faces with `muf()` and `Kf()`.

Note that the class does not define a destructor and does not call
`F().storeOldTime()`; the base class handles old-time creation lazily when
required.

### `correct()` walkthrough

`correct(volSymmTensorField&)`:

1. calls `updateF(sigma, mu(), K())` and returns early if the solid model has
   enforced linearity for stability;
2. forms `J`, `bBar` and `bBar & bBar`, then `I1Bar` and `I2Bar`;
3. builds the deviatoric stress from the two coefficient groups;
4. calls `updateSigmaHyd(0.5*K()*(J^2 - 1), (4/3)*mu() + K())`;
5. assembles
   `sigma = (1/J)*(dev(s) + sigmaHyd*I + symm(F & sigma0 & F.T()))`.

The surface overload mirrors this using `Ff()`, `c10f_`, `c01f_`, `c11f_` and
`Kf()`. It deliberately bypasses `updateSigmaHyd()` and computes
`0.5*Kf()*(J^2 - 1)` locally, so the pressure smoothing enabled by
`solvePressureEqn` applies to the cell-centred stress only.

### Implicit stiffness

`impK()` returns `(4/3)*mu() + K()`, evaluated from the linearised moduli. It
is a genuine field, so the segregated diffusivity follows any spatial
variation in the coefficients. It does not stiffen as the material stiffens
under large strain, so heavily deformed cases may need more outer correctors.

### Extension points

- Additional Rivlin terms only require extending the deviatoric-stress
  expression in both `correct()` overloads and updating the linearised `mu`.
- The face coefficients `c10f_`, `c01f_` and `c11f_` are interpolated once at
  construction; a law with evolving coefficients would need to refresh them.

---

## Tutorials

Cases that select `MooneyRivlinElastic`:

- `solids/hyperelasticity/cylinderCrush`
- `solids/hyperelasticity/cylindricalPressureVessel`
- `solids/hyperelasticity/cylindricalPressureVessel/caseOptions/displacement`
- `solids/hyperelasticity/longWall`
