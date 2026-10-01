---
sort: 10
---

# anisotropicBiotElastic

Small-strain orthotropic linear elasticity, with the stiffness given as
Young's moduli, Poisson's ratios and shear moduli in the coordinate
directions. Pore pressure is included by nesting it inside
`poroMechanicalLaw`. The runtime type is:

```text
anisotropicBiotElastic
```

The dictionary layout common to all laws is described in the
[material models page](../../../materialModels/README.md).

---

## User Guide

### What it computes

With `epsilon = symm(grad(D))` from the solid model, the stress is
`sigma = A : epsilon`, component by component. The form is chosen by the
mesh, not by the dictionary.

On a 3-D mesh:

```text
sigma_xx = A11*epsilon_xx + A12*epsilon_yy + A31*epsilon_zz
sigma_yy = A12*epsilon_xx + A22*epsilon_yy + A23*epsilon_zz
sigma_zz = A31*epsilon_xx + A23*epsilon_yy + A33*epsilon_zz
sigma_xy = A44*epsilon_xy
sigma_yz = A55*epsilon_yz
sigma_xz = A66*epsilon_xz
```

with

```text
nuyx = nuxy*Ey/Ex
nuxz = nuzx*Ex/Ez
nuzy = nuyz*Ez/Ey
J    = (1 - nuxy*nuyx - nuyz*nuzy - nuzx*nuxz
        - 2*nuyx*nuzy*nuxz)/(Ex*Ey*Ez)
A11  = (1 - nuyz*nuzy)/(J*Ey*Ez)
A22  = (1 - nuxz*nuzx)/(J*Ex*Ez)
A33  = (1 - nuyx*nuxy)/(J*Ey*Ex)
A12  = (nuxy + nuzy*nuxz)/(J*Ex*Ez)
A31  = (nuzx + nuyx*nuzy)/(J*Ey*Ez)
A23  = (nuyz + nuyx*nuxz)/(J*Ex*Ey)
A44  = 2*Gxy,  A55 = 2*Gyz,  A66 = 2*Gzx
```

On a 2-D mesh whose empty direction is z, the plane stress reduction is used:

```text
sigma_xx = A11*epsilon_xx + A12*epsilon_yy
sigma_yy = A21*epsilon_xx + A22*epsilon_yy
sigma_xy = A44*epsilon_xy
sigma_zz = sigma_yz = sigma_xz = 0

nuyx = nuxy*Ey/Ex
J    = 1/(1 - nuxy*nuyx)
A11  = J*Ex,  A22 = J*Ey,  A12 = J*nuyx*Ex,  A21 = J*nuxy*Ey
A44  = 2*Gxy
```

The law is small-strain only; a nonlinear geometry solid model stops with an
error. It keeps no history.

### Model options

The elastic constants are read as plain scalars, without dimension sets; use
one consistent pressure unit for all moduli.

| Entry | Required | Default | Description |
| --- | --- | --- | --- |
| `rho` | yes | - | Density, `[1 -3 0 0 0 0 0]` |
| `Ex` | yes | - | Young's modulus in x |
| `Ey` | yes | - | Young's modulus in y |
| `Ez` | 3-D mesh | - | Young's modulus in z |
| `nuxy` | yes | - | Poisson's ratio `nuxy` |
| `nuyz` | 3-D mesh | - | Poisson's ratio `nuyz` |
| `nuzx` | 3-D mesh | - | Poisson's ratio `nuzx` |
| `Gxy` | yes | - | Shear modulus in the xy plane |
| `Gyz` | 3-D mesh | - | Shear modulus in the yz plane |
| `Gzx` | 3-D mesh | - | Shear modulus in the zx plane |

On a 2-D mesh the z entries are not read. When the law is nested inside
`poroMechanicalLaw`, `rho` is taken from the enclosing law's dictionary and
need not be repeated.

The following are fatal errors at construction:

- a 2-D mesh without `planeStress yes;` in `mechanicalProperties`: the 2-D
  reduction is plane stress, and there is no plane strain form;
- a 2-D mesh whose empty direction is x or y.

`solutionD` is supplied to the law by the framework from the mesh and must not
be given in the dictionary. The moduli, Poisson's ratios and the determinant
`J` are not checked.

### Tangents and volumetric split

- Scalar tangent: `max(A11, A22, A33)` for both the scalar and the deviatoric
  scalar tangent; the same value is returned as the law's bulk modulus.
- Fourth-order tangent: `fourthOrderFiniteDifference` only. Asking for the
  analytical `fourthOrder` tangent is a fatal error.
- Volumetric split: no. The law cannot be used with `solvePressureEqn` or with
  the mixed displacement-pressure formulations.

### State variables

None. The law declares no state, so it writes no state fields or restart
files.

### Example

From the `rodAndSeabed` tutorial, a 3-D case:

```text
planeStress     no;

mechanical
(
    steel
    {
        type            poroMechanicalLaw;
        rho             rho [1 -3 0 0 0 0 0] 2650;
        biotCoeff       biotCoeff [0 0 0 0 0 0 0] 1.0;
        effectiveStressMechanicalLaw
        {
            type            anisotropicBiotElastic;

            // Young's moduli (in Pa)
            Ex              1.2e7;
            Ey              1.2e7;
            Ez              2e7;

            // Poisson's ratios (dimensionless)
            nuxy            0.2;
            nuyz            0.24;
            nuzx            0.4;

            // Shear moduli (in Pa)
            Gxy             0.5e7;
            Gyz             1.2e7;
            Gzx             1.2e7;
        }
    }
);
```

### Migrating from the legacy law

- The 2-D/3-D selection is fixed
  ([issue #334](https://github.com/solids4foam/solids4foam/issues/334)). The
  legacy law took its reduced branch on 3-D meshes, silently ignoring `Ez`,
  `nuyz`, `nuzx`, `Gyz` and `Gzx` and leaving the out-of-plane stresses at
  zero. A 3-D case now uses all nine constants, so its results change.
- A 2-D case (empty z) now needs `planeStress yes;` and only the four in-plane
  constants; the legacy law required all nine there.
- In the 2-D form the out-of-plane stress components are now set to zero
  rather than left unassigned.
- `rho` is read at construction; `regionName` is no longer read. The legacy
  law silently ignored `solvePressureEqn`; the solid model now refuses it for
  this law.
- The law now serves every integration-point location, and a fourth-order
  tangent is available by finite differences.

---

## Tutorials

- [rodAndSeabed](../../../../../tutorials/solids/poroelasticity/rodAndSeabed/README.md),
  inside the `effectiveStressMechanicalLaw` sub-dictionary of
  `poroMechanicalLaw`
