---
sort: 16
---

# StVenantKirchhoffOrthotropicElastic

This page documents the orthotropic St. Venant-Kirchhoff hyperelastic law. The
runtime type is:

```text
StVenantKirchhoffOrthotropicElastic
```

The law is the orthotropic generalisation of
`StVenantKirchhoffElastic`: the second Piola-Kirchhoff stress is a linear
function of the Green-Lagrange strain through a fourth-order orthotropic
stiffness tensor, expressed in a user-defined material coordinate system.

```warning
The entire implementation file is wrapped in `#ifdef FOAMEXTEND`. Although the
source is listed in both `Make/files.openfoam` and `Make/files.foamextend`, it
compiles to an empty translation unit on OpenFOAM.com and OpenFOAM.org
builds, so the law is **only selectable in foam-extend**.
```

---

## User Guide

### Constitutive relation

```text
S = C : E,   E = 0.5*(F.T() & F - I),   sigma = (1/J)*F & S & F.T()
```

`C` is the orthotropic stiffness in the material frame, built from nine
independent constants. With

```text
Jc = (1 - nu12*nu21 - nu23*nu32 - nu31*nu13 - 2*nu21*nu32*nu13)/(E1*E2*E3)
```

the normal-stress block is

```text
C11 = (1 - nu23*nu32)/(Jc*E2*E3)
C22 = (1 - nu13*nu31)/(Jc*E1*E3)
C33 = (1 - nu21*nu12)/(Jc*E2*E1)
C12 = (nu12 + nu32*nu13)/(Jc*E1*E3)
C31 = (nu31 + nu21*nu32)/(Jc*E2*E3)
C23 = (nu23 + nu21*nu13)/(Jc*E1*E2)
```

and the shear block is `C44 = 2*G12`, `C55 = 2*G23`, `C66 = 2*G31`. The
reciprocal Poisson's ratios are derived from the symmetry conditions
`nu21 = nu12*E2/E1`, `nu32 = nu23*E3/E2` and `nu13 = nu31*E1/E3`; do not
specify them.

`C` is assembled in the local material frame and then rotated into the global
frame with the `materialDirections` tensor field.

Registered in the nonlinear-geometry table, the law works with both
total-Lagrangian and updated-Lagrangian solid models. Plane stress is
explicitly rejected at construction.

### Model options

All nine elastic constants are required and are read with `lookup()`:

| Entry | Dimensions | Description |
| --- | --- | --- |
| `E1`, `E2`, `E3` | `[1 -1 -2 0 0 0 0]` | Young's moduli, material axes |
| `nu12`, `nu23`, `nu31` | `[0 0 0 0 0 0 0]` | Independent Poisson's ratios |
| `G12`, `G23`, `G31` | `[1 -1 -2 0 0 0 0]` | Shear moduli |
| `rho` | `[1 -3 0 0 0 0 0]` | Density |

Two optional entries set the material axes:

| Entry | Default | Description |
| --- | --- | --- |
| `materialDirection1` | `(1 0 0)` | First material axis |
| `materialDirection2` | `(0 1 0)` | Second material axis |

Both are read with `lookupOrAddDefault`, so a missing entry is written back
into the dictionary.

```warning
There is no working third material axis. The constructor reads
`materialDirection2` a second time, with a default of `(0 0 1)`, in the slot
where `materialDirection3` is intended. The third row of the rotation tensor
therefore always ends up equal to the second, which makes the tensor singular
and the resulting global stiffness wrong. A `materialDirection3` entry in the
dictionary is silently ignored. Until this is fixed in the source, supply the
orientation through a `materialDirections` field instead.
```

The two entries only build the uniform default of the `materialDirections`
tensor field, whose rows are the three material axes. A spatially varying, and
correct, orientation can be supplied by placing a `materialDirections`
`volTensorField` in the time directory; it is read `READ_IF_PRESENT` and takes
precedence over the uniform default. The rows are normalised after reading,
and zero-length rows are a fatal error.

### Physical-admissibility checks

The constructor aborts if:

- any of `E1`, `E2`, `E3`, `G12`, `G23`, `G31` is negative;
- `mag(nu_ij)` is not less than `sqrt(E_i/E_j)` for the three supplied ratios;
- the determinant term
  `1 - nu12*nu21 - nu23*nu32 - nu31*nu13 - 2*nu21*nu32*nu13` is not positive;
- `planeStress` is `yes`.

### Recommended dictionary setup

Minimal example for `constant/mechanicalProperties`:

```text
planeStress     no;

mechanical
(
    composite
    {
        type            StVenantKirchhoffOrthotropicElastic;
        rho             rho [1 -3 0 0 0 0 0] 1600;

        E1              E1 [1 -1 -2 0 0 0 0] 130e+09;
        E2              E2 [1 -1 -2 0 0 0 0] 10e+09;
        E3              E3 [1 -1 -2 0 0 0 0] 10e+09;

        nu12            nu12 [0 0 0 0 0 0 0] 0.3;
        nu23            nu23 [0 0 0 0 0 0 0] 0.4;
        nu31            nu31 [0 0 0 0 0 0 0] 0.023;

        G12             G12 [1 -1 -2 0 0 0 0] 5e+09;
        G23             G23 [1 -1 -2 0 0 0 0] 3.6e+09;
        G31             G31 [1 -1 -2 0 0 0 0] 5e+09;

        materialDirection1 (1 0 0);
        materialDirection2 (0 1 0);
    }
);
```

### Field glossary

- `materialDirections`: tensor field whose rows are the normalised material
  axes; read from the time directory if present.
- `elasticC`, `elasticCf`: fourth-order stiffness in the global frame, at cell
  centres and faces; demand-driven and not written.
- `F`, `Ff`: total deformation gradient.
- `sigma`: Cauchy stress tensor returned to the solid model.
- `impK`: implicit stiffness used as the Laplacian diffusivity.

---

## Developer Notes

### Class role

`StVenantKirchhoffOrthotropicElastic` inherits from `mechanicalLaw` and is
added to the `nonLinGeomMechLaw` runtime selection table. It relies on the
`symmTensor4thOrder` type and its `transform` overloads, which exist only in
foam-extend; this is the reason for the `FOAMEXTEND` guard.

### Construction

The nine constants and the three reciprocal ratios are set in the initialiser
list, together with the `materialDirections` field. The body then runs the
admissibility checks, normalises the direction vectors, and initialises the
scalar `mu_` and `K_` used by the base class:

```text
mu = (G12 + G23 + G31)/3
K  = Ebar*mu/(3*(3*mu - Ebar)),   Ebar = (E1 + E2 + E3)/3
```

These two scalars are only passed to `updateF()`, where they define the
isotropic Hooke's law fallback used when the solid model enforces linearity.

### Demand-driven stiffness

`elasticC()` and `elasticCf()` build their fields on first access through
`makeElasticC()` and `makeElasticCf()`, which duplicate the same nine
coefficient expressions. The cell-centred tensor is rotated by `matDir_` and
the face tensor by `fvc::interpolate(matDir_)`. Both pointers are released in
the destructor via `deleteDemandDrivenData`.

### `correct()` walkthrough

Both overloads:

1. call `updateF(sigma, mu_, K_)` and return early if enforced linearity is
   active;
2. materialise `FT` before contracting, as required by the note in
   `StVenantKirchhoffElastic.C`;
3. form `C = symm(FT & F)` and the Green strain `E = 0.5*(C - I)`;
4. evaluate `S = elasticC() && E`;
5. push forward with `sigma = (1/J)*transform(F, S)`.

### Implicit stiffness

`impK()` returns the average of the three normal diagonal components of
`elasticC`, `(C_XXXX + C_YYYY + C_ZZZZ)/3`. A `diagTensor` diffusivity would
be more faithful, but the source comment notes that the scalar average gives
similar convergence in practice.

### Extension points

- The coefficient assembly is duplicated between `makeElasticC()` and
  `makeElasticCf()`; factor it out before changing the stiffness definition.
- Fixing the third material direction requires changing the second
  `materialDirection2` lookup in the `matDir_` initialiser to
  `materialDirection3`.
