---
sort: 4
---

# Uniform heating of a plane stress plate: `heatedPlatePlaneStress`

---

Prepared by Philip Cardiff

---

## Tutorial Aims

- Demonstrate a thermo-elastic analysis in plane stress;
- Verify the thermal stress of `thermoMechanicalLaw` in plane stress against
  the analytical solution.

---

## Case Overview

A $$1 \times 1$$ m plate is heated uniformly from the stress-free reference
temperature $$T_0 = 300$$ K to $$T = 400$$ K. The plate is linear elastic,
with $$E = 200$$ GPa and $$\nu = 0.3$$, and its coefficient of linear thermal
expansion is $$\alpha = 10^{-5}$$ K$$^{-1}$$. The case is two-dimensional and
solved under the plane stress assumption (`planeStress yes;` in
`constant/mechanicalProperties`), on a mesh of $$10 \times 10$$ cells.

The `thermalLinearGeometry` solid model is used. The temperature is
$$T = 400$$ K initially and on the right and top edges, so the steady
temperature field is uniform. Thermal inertia and gravity are neglected, and
one loading increment is used.

The left and bottom edges are symmetry planes, which remove the rigid-body
motion, and the right and top edges are traction-free, so the plate is free to
expand.

In plane stress, the thermal stress of an isotropic linear elastic material is

$$
\boldsymbol{\sigma} = \frac{E}{1 - \nu^2}\left[(1 - \nu)\boldsymbol{\varepsilon}
+ \nu\,\text{tr}(\boldsymbol{\varepsilon})\boldsymbol{I}\right]
- \frac{E\alpha(T - T_0)}{1 - \nu}\boldsymbol{I}
$$

for the in-plane components. A free plate is therefore stress-free, with the
displacement $$\boldsymbol{u} = \alpha (T - T_0)(x, y)$$, and a clamped plate
has $$\sigma_{xx} = \sigma_{yy} = -E\alpha(T - T_0)/(1 - \nu)$$.

---

## Running the Case

The tutorial case is located at
`solids4foam/tutorials/solids/thermoelasticity/heatedPlatePlaneStress`. The
case can be run using the included `Allrun` script, i.e. `> ./Allrun`. The
`Allrun` creates the mesh using `blockMesh` (`> blockMesh`), after which the
case is run with the `solids4Foam` solver (`> solids4Foam`).

---

## Expected Results

The plate expands without stress, and the displacement at the top-right
corner is $$\alpha (T - T_0)(1, 1) = (10^{-3}, 10^{-3})$$ m.

In plane stress, only the in-plane components of the stress and strain are
meaningful: the out-of-plane components are not post-processed, so they are
not the zero stress and the free out-of-plane strain, and should be ignored.
For this reason, the `Max sigmaEq` printed in the log is not zero: it comes
from $$\sigma_{zz} = -E\alpha(T - T_0)/(1 + \nu)$$, about 154 MPa, while the
in-plane stress is zero.

---

## Regression Testing

The `regressionTest.sh` script runs the case as it is and checks the
displacement against the analytical solution. It then runs it again with the
right and top edges clamped, and checks the stress against
$$-E\alpha(T - T_0)/(1 - \nu)$$. Both checks failed before the plane stress
thermal stress was corrected
([#137](https://github.com/solids4foam/solids4foam/issues/137)), when it was
$$(1 + 2\nu)/(1 + \nu)$$ times too large.

```bash
./regressionTest.sh
```
