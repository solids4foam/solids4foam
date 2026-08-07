---
sort: 10
---

# anisotropicBiotElastic

This law calculates small-strain orthotropic effective stress. Pore pressure
can be included by nesting it inside `poroMechanicalLaw`. The runtime type is:

```text
anisotropicBiotElastic
```

---

## User Guide

### What it computes

The cell-centred `correct` obtains the total strain from the displacement
gradient:

```text
epsilon = symm(grad(D))                         // total formulation
epsilon = epsilon.oldTime() + symm(grad(DD))   // incremental formulation
```

In the branch labelled as 2-D, the implemented stress relation is:

```text
sigma_xx = A11*epsilon_xx + A12*epsilon_yy
sigma_yy = A21*epsilon_xx + A22*epsilon_yy
sigma_xy = A44*epsilon_xy
```

The other stress components are not assigned in this branch. In the 3-D
branch, the relation is:

```text
sigma_xx = A11*epsilon_xx + A12*epsilon_yy + A31*epsilon_zz
sigma_yy = A12*epsilon_xx + A22*epsilon_yy + A23*epsilon_zz
sigma_zz = A31*epsilon_xx + A23*epsilon_yy + A33*epsilon_zz
sigma_xy = A44*epsilon_xy
sigma_yz = A55*epsilon_yz
sigma_xz = A66*epsilon_xz
```

The same relations are applied to the boundary values. The law implements
only `correct(volSymmTensorField&)`. Its surface overload aborts with
`notImplemented`. The inherited point overload also aborts. The inherited
`CompactListList` overload is unavailable on foam-extend and aborts with
`notImplemented` on the other forks.

### Model options

The elastic properties are read as plain scalars, so the dictionary does not
attach dimension sets to them. Use one consistent pressure unit for all
Young's and shear moduli.

| Entry | Required | Description |
| --- | --- | --- |
| `Ex` | yes | Young's modulus in the x direction |
| `Ey` | yes | Young's modulus in the y direction |
| `Ez` | in the 3-D branch | Young's modulus in the z direction |
| `nuxy` | yes | Poisson's ratio `nuxy` |
| `nuyz` | in the 3-D branch | Poisson's ratio `nuyz` |
| `nuzx` | in the 3-D branch | Poisson's ratio `nuzx` |
| `Gxy` | yes | Shear modulus in the xy plane |
| `Gyz` | in the 3-D branch | Shear modulus in the yz plane |
| `Gzx` | in the 3-D branch | Shear modulus in the zx plane |
| `rho` | for direct use | Density, `[1 -3 0 0 0 0 0]` |
| `regionName` | no | Base mesh region; otherwise `solid` or `region0` |
| `solvePressureEqn` | no | Accepted switch; default `no`, unused here |
| `pressureSmoothingScaleFactor` | no | Default `100`, unused here |

Which entries are required depends on a branch flag that is currently
inverted, so read the following carefully before setting up a case.

```warning
The law selects its 2-D branch when `mesh.solutionD()[vector::Z] > 0`, which
is true when z is a **solved** direction — that is, for a **3-D** mesh. The
sense of the test is the opposite of the one intended, so today a 3-D case
takes the 2-D branch and reads only `Ex`, `Ey`, `nuxy` and `Gxy`, silently
ignoring `Ez`, `nuyz`, `nuzx`, `Gyz` and `Gzx` and leaving the out-of-plane
stress components at zero. A genuinely 2-D case, with z empty, takes the 3-D
branch and aborts with `keyword Ez is undefined` unless all nine constants are
supplied. This is tracked as
[issue #334](https://github.com/solids4foam/solids4foam/issues/334); supply all
nine constants until it is fixed.
```

`rho` is read through the base-class density functions. When the law is nested
inside `poroMechanicalLaw`, density belongs to the outer law, as in the setup
below. The pressure-equation entries are stored by `mechanicalLaw`, but this
law never calls `updateSigmaHyd()`, so they do not change its stress update.

### Recommended dictionary setup

The following seabed-soil values follow the available poroelastic tutorial:

```text
planeStress     no;

mechanical
(
    seabedSoil
    {
        type            poroMechanicalLaw;
        rho             rho [1 -3 0 0 0 0 0] 2650;
        biotCoeff       biotCoeff [0 0 0 0 0 0 0] 1.0;

        effectiveStressMechanicalLaw
        {
            type        anisotropicBiotElastic;
            Ex          1.2e7;
            Ey          1.2e7;
            Ez          2e7;
            nuxy        0.2;
            nuyz        0.24;
            nuzx        0.4;
            Gxy         0.5e7;
            Gyz         1.2e7;
            Gzx         1.2e7;

            // Optional for a non-standard base mesh name
            // regionName solid;
        }
    }
);
```

### Field glossary

- `epsilon`: cell-centred total small-strain tensor. It is created with
  `NO_READ` and `NO_WRITE` and accumulated from its old-time value for an
  incremental solid model.
- `grad(D)`, `grad(DD)`: displacement-gradient fields looked up for total and
  incremental solid models respectively; the law does not own them.
- `sigma`: caller-owned cell-centred stress updated by `correct`.
- `impK`: temporary uniform pressure-dimensioned field returned to the solid
  model.

---

## Developer Notes

### Class role

`anisotropicBiotElastic` derives directly from `mechanicalLaw`. It stores the
branch flag `model2d_`, ten scalar stiffness coefficients and the
cell-centred `epsilon_` field. It is registered in the `linGeomMechLaw`
selection table under `anisotropicBiotElastic`.

The source uses `#ifdef OPENFOAM_NOT_EXTEND` to select writable internal and
boundary field accessors. The class itself has no `FOAMEXTEND`-only interface
guards. Its source is listed in both `Make/files.openfoam` and
`Make/files.foamextend`.

### Construction

The base constructor first reads `solvePressureEqn` with default `false` and
`pressureSmoothingScaleFactor` with default `100.0`. It then reads
`regionName`, or selects an existing `solid` or `region0` mesh. If none is
available, it raises a fatal error asking for `regionName`.

The derived constructor sets `model2d_` from whether the z solution direction
is greater than zero, zeroes every stiffness coefficient and creates
`epsilon_`. In the 2-D branch it raises a fatal error if x or y is an empty
direction. It then reads `Ex`, `Ey`, `nuxy` and `Gxy`, and derives:

```text
nuyx = nuxy*Ey/Ex
J    = 1/(1 - nuxy*nuyx)
A11  = J*Ex
A22  = J*Ey
A12  = J*nuyx*Ex
A21  = J*nuxy*Ey
A44  = 2*Gxy
```

In the 3-D branch it reads all nine elastic entries and derives:

```text
nuyx = nuxy*Ey/Ex
nuxz = nuzx*Ex/Ez
nuzy = nuyz*Ez/Ey
J = (1 - nuxy*nuyx - nuyz*nuzy - nuzx*nuxz
     - 2*nuyx*nuzy*nuxz)/(Ex*Ey*Ez)
A11 = (1 - nuyz*nuzy)/(J*Ey*Ez)
A22 = (1 - nuxz*nuzx)/(J*Ex*Ez)
A33 = (1 - nuyx*nuxy)/(J*Ey*Ex)
A12 = (nuxy + nuzy*nuxz)/(J*Ex*Ez)
A31 = (nuzx + nuyx*nuzy)/(J*Ey*Ez)
A23 = (nuyz + nuyx*nuxz)/(J*Ex*Ey)
A44 = 2*Gxy
A55 = 2*Gyz
A66 = 2*Gzx
```

The constructor prints the supplied properties and derived coefficients. It
does not validate moduli, Poisson's ratios or the denominators, and emits no
law-specific warnings.

### Key methods

- `impK()` returns a uniform field containing
  `max(A11, max(A22, A33))`. In the 2-D branch `A33` remains zero. This is the
  diffusivity of the solid model's implicit Laplacian term and affects the
  outer-iteration convergence rate rather than the converged answer.
- `correct(volSymmTensorField&)` updates `epsilon_` from `grad(D)` or
  `grad(DD)`, then applies the component-wise constitutive relation to cell
  and boundary values.
- `correct(surfaceSymmTensorField&)` always aborts with `notImplemented`.
- `materialTangent()` is not overridden. The inherited implementation aborts
  with `notImplemented`, so this law does not supply a Newton material
  tangent.

### Extension points

A related orthotropic law can copy this class, replace the property reads and
stiffness derivation, and update the component equations in `correct`. It must
be registered in the appropriate mechanical-law selection table and added to
both build lists if all forks are supported. Implement the surface, point and
quadrature overloads, and `materialTangent()`, when the target solid models
require those interfaces.

The source is at
[anisotropicBiotElastic.C][source].

[source]: https://github.com/solids4foam/solids4foam/blob/master/src/solids4FoamModels/materialModels/mechanicalModel/mechanicalLaws/linearGeometryLaws/anisotropicBiotElastic/anisotropicBiotElastic.C

---

## Tutorials

- `solids/poroelasticity/rodAndSeabed`
