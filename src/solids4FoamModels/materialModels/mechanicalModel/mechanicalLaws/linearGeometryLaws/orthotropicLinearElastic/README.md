---
sort: 4
---

# orthotropicLinearElastic

This page documents the orthotropic Hookean linear elastic mechanical law. The
runtime type is:

```text
orthotropicLinearElastic
```

The law relates stress to strain through nine independent elastic constants:
three Young's moduli, three shear moduli and three Poisson's ratios. The
constants are given in a local material coordinate system, which under
foam-extend can then be rotated into global Cartesian coordinates using
per-cell material direction fields.

The formulation follows P. Cardiff, A. Karač and A. Ivanković (2014), _A large
strain finite volume method for orthotropic bodies with general material
orientations_, Computer Methods in Applied Mechanics and Engineering, 268(1),
318-335,
[10.1016/j.cma.2013.09.008](https://doi.org/10.1016/j.cma.2013.09.008).

---

## User Guide

### What it computes

The law builds a fourth-order stiffness tensor `C` from the nine constants and
evaluates

```text
sigma = C && epsilon
```

at cell centres and at faces. Only these two `correct` overloads are
implemented, so the law works with the cell-centred and face-based solid
models but not with the vertex-centred ones.

### Model options

All nine entries below are read with `lookup()` and are therefore
**required**; the law aborts at construction if any is missing.

| Entry | Dimensions | Description |
| --- | --- | --- |
| `E1` | `[1 -1 -2 0 0 0 0]` | Young's modulus, local 1 direction |
| `E2` | `[1 -1 -2 0 0 0 0]` | Young's modulus, local 2 direction |
| `E3` | `[1 -1 -2 0 0 0 0]` | Young's modulus, local 3 direction |
| `nu12` | dimensionless | Poisson's ratio, 1-2 plane |
| `nu23` | dimensionless | Poisson's ratio, 2-3 plane |
| `nu31` | dimensionless | Poisson's ratio, 3-1 plane |
| `G12` | `[1 -1 -2 0 0 0 0]` | Shear modulus, 1-2 plane |
| `G23` | `[1 -1 -2 0 0 0 0]` | Shear modulus, 2-3 plane |
| `G31` | `[1 -1 -2 0 0 0 0]` | Shear modulus, 3-1 plane |

The reciprocal ratios are derived rather than read:

```text
nu21 = nu12*E2/E1
nu32 = nu23*E3/E2
nu13 = nu31*E1/E3
```

`rho` is also required, as for every mechanical law.

### Material directions (foam-extend only)

Under foam-extend, three optional vector entries set the local material axes:

| Entry | Default | Description |
| --- | --- | --- |
| `materialDirectionX` | `(1 0 0)` | Local 1 axis |
| `materialDirectionY` | `(0 1 0)` | Local 2 axis |
| `materialDirectionZ` | `(0 0 1)` | Local 3 axis |

Each is used as the uniform value of a corresponding `materialDirectionsX`,
`materialDirectionsY` or `materialDirectionsZ` field, which is itself
`READ_IF_PRESENT`. Supplying those fields in the time directory therefore
gives a spatially varying material orientation, which is the point of the law.
The vectors are normalised at construction and must be locally orthogonal and
of non-zero length.

```warning
Under OpenFOAM.com and OpenFOAM.org these fields are constructed `NO_READ`
with hard-coded global axes, and the rotation of `C` into global coordinates
is skipped entirely. With those forks the law is restricted to materials whose
principal axes coincide with the global x, y and z axes; the
`materialDirection*` entries are ignored.
```

### Physical constraints

The constructor rejects a configuration if any of the following holds:

- `planeStress` is `yes` — the law is not implemented for plane stress;
- any of `E1`, `E2`, `E3`, `G12`, `G23`, `G31` is negative;
- `mag(nu_ij) >= sqrt(E_i/E_j)` for any of the three supplied ratios;
- `1 - nu12*nu21 - nu23*nu32 - nu31*nu13 - 2*nu21*nu32*nu13 <= 0`, which is
  the positive-definiteness condition on the compliance matrix.

### Recommended dictionary setup

```text
planeStress     no;

mechanical
(
    composite
    {
        type            orthotropicLinearElastic;
        rho             rho [1 -3 0 0 0 0 0] 1600;

        E1              E1 [1 -1 -2 0 0 0 0] 140e9;
        E2              E2 [1 -1 -2 0 0 0 0] 10e9;
        E3              E3 [1 -1 -2 0 0 0 0] 10e9;

        nu12            nu12 [0 0 0 0 0 0 0] 0.3;
        nu23            nu23 [0 0 0 0 0 0 0] 0.4;
        nu31            nu31 [0 0 0 0 0 0 0] 0.021;

        G12             G12 [1 -1 -2 0 0 0 0] 5e9;
        G23             G23 [1 -1 -2 0 0 0 0] 3.6e9;
        G31             G31 [1 -1 -2 0 0 0 0] 5e9;

        // foam-extend only
        // materialDirectionX (1 0 0);
        // materialDirectionY (0 1 0);
        // materialDirectionZ (0 0 1);
    }
);
```

### Field glossary

- `epsilon`, `epsilonf`: small-strain tensor at cells and faces, with old-time
  values stored.
- `elasticC`, `elasticCf`: the stiffness tensor at cells and faces, built on
  first use.
- `materialDirectionsX/Y/Z`: unit vector fields defining the local axes
  (foam-extend only).
- `impK`: implicit stiffness, the mean of the three normal-direction diagonal
  stiffnesses.

---

## Developer Notes

### Class role

`orthotropicLinearElastic` derives from `mechanicalLaw`. Its distinguishing
feature is a fork-dependent storage type for the stiffness tensor:

- under foam-extend, `volSymmTensor4thOrderField` / its surface counterpart,
  which support `transform()` and the `&&` double-inner-product directly;
- under OpenFOAM, a plain `volTensorField` is used to hold the same nine
  components, with the packing documented in the header:
  `XXXX, XXYY, XXZZ, YYYY, YYZZ, ZZZZ, XYXY, YZYZ, ZXZX` mapped onto the
  tensor components `XX, XY, XZ, YX, YY, YZ, ZX, ZY, ZZ`. The double product
  is then performed by the free function in `doubleDotProduct.H`.

The law is listed in both `Make/files.openfoam` and `Make/files.foamextend`.

### Construction

Elastic constants are read directly in the initialiser list; the derived
ratios, the direction fields and the physical checks follow. `elasticC` and
`elasticCf` are _not_ built in the constructor: `makeElasticC()` and
`makeElasticCf()` are called lazily on first access through `elasticC()` and
`elasticCf()`, each guarded against double initialisation.

### Key methods

- `impK()`: returns `(C_XXXX + C_YYYY + C_ZZZZ)/3`. The code notes that a
  `diagTensor` implicit stiffness would be more faithful, but that the
  averaged scalar achieves comparable convergence.
- `correct(volSymmTensorField&)`, `correct(surfaceSymmTensorField&)`: update
  the strain then apply the double product.

Both `makeElasticC()` and `makeElasticCf()` recompute the same `A11` through
`A66` coefficients from scratch; the duplication is harmless because each runs
at most once.

### Extension points

Adding the missing global rotation for the OpenFOAM forks is the obvious
improvement: it requires a `transform()` equivalent for the nine-component
tensor packing, after which the `#ifdef FOAMEXTEND` guards around the
direction fields can be removed. Implementing `materialTangent()` would also
let the law be used by the block-coupled and PETSc SNES solid models.

The source is at
[orthotropicLinearElastic.C](https://github.com/solids4foam/solids4foam/blob/master/src/solids4FoamModels/materialModels/mechanicalModel/mechanicalLaws/linearGeometryLaws/orthotropicLinearElastic/orthotropicLinearElastic.C).
