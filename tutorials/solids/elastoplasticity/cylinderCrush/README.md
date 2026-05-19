---
sort: x
---

# Pressing a Long Cylinder: `cylinderCrush`

Prepared by Philip Cardiff and Ivan Batistić

---

## Tutorial Aims

- Demonstrate a large-strain elastoplastic analysis with frictionless solid
  contact.

```note
A hyperelastic version of this tutorial is also available in solids4foam at
`tutorials/solids/hyperelasticity/cylinderCrush`.
```

---

## Case Overview

This case considers a long solid cylinder crushed between two rigid plates, as
shown in Figure 1. The cylinder diameter is $$0.4$$ m, and the cylinder is
compressed by half of its diameter. Because the problem is symmetric, only one
quarter of the cylinder is modelled. The rigid plate in the computational
domain is displaced by $$0.1$$ m over 30 equal loading increments.

The plate-cylinder contact is modelled using the `solidContact` boundary
condition with a frictionless contact law. Gravity and transient effects are
neglected, so the case is solved as a quasi-static 2-D plane-strain problem. The
updated Lagrangian solid model is selected in `constant/solidProperties` using
`nonLinearGeometryUpdatedLagrangian`.

![Figure 1: Solid cylinder compressed between two rigid
plates](./images/cylinderCrush-geometry.png)

Figure 1 - Solid cylinder compressed between two rigid plates

The cylinder is modelled using the `neoHookeanElasticMisesPlastic` mechanical
law. The mechanical properties are given in Table 1, and the tabulated
hardening curve is given in Table 2. The hardening data is the same as in the
`elastoplasticity/pipeCrush` tutorial and was taken from Cardiff et al. [1].

Table 1 - Mechanical properties

| Parameter            | Symbol       | Value           |
| -------------------- | ------------ | --------------- |
| Young's modulus      | $$E$$        | $$186$$ GPa     |
| Poisson's ratio      | $$\nu$$      | $$0.3$$         |
| Initial yield stress | $$\sigma_Y$$ | $$241.32$$ MPa  |
| Initial density      | $$\rho$$     | $$7800$$ kg/m^3 |

Table 2 - Hardening behaviour

| Plastic strain | Yield stress (in MPa) |
| -------------- | --------------------- |
| $$0.000$$      | $$241.32$$            |
| $$0.0035$$     | $$275.79$$            |
| $$0.0083$$     | $$301.65$$            |
| $$0.0133$$     | $$318.88$$            |
| $$0.0182$$     | $$344.74$$            |
| $$0.0281$$     | $$361.98$$            |
| $$0.0380$$     | $$379.21$$            |
| $$1.000$$      | $$1482.37$$           |

---

## Expected Results

Figure 2 shows the initial mesh and the deformed mesh coloured by the von Mises
stress field at an intermediate loading step and at the final loading step. 

![Figure 2: Initial and deformed meshes coloured by von Mises
stress](./images/cylinderCrush-results.jpeg)

Figure 2 - Initial and deformed meshes coloured by von Mises stress

The `solidForces` and `solidPointDisplacement` function objects are added in
`system/controlDict`. They record the contact force on the `cylinderContact`
patch and the displacement of the point at `(0 0.2 0)`, respectively. These
outputs can be used to plot the contact force as a function of the plate
displacement. Reference force-displacement data is not available in the
literature for direct comparison.

---

## Running the Case

The tutorial case is located at
`solids4foam/tutorials/solids/elastoplasticity/cylinderCrush`. The case can be
run using the included `Allrun` script, i.e. `> ./Allrun`.

The `Allrun` script creates the mesh using `blockMesh` and then runs the
`solids4Foam` solver.

---

### References

[1]
[P. Cardiff, Ž. Tuković, P. D. Jaeger, M. Clancy, and A. Ivanković, "A
Lagrangian cell-centred finite volume method for metal forming simulation,"
International Journal for Numerical Methods in Engineering, vol. 109, no. 13,
pp. 1777-1803, 2017.](https://onlinelibrary.wiley.com/doi/abs/10.1002/nme.5345)
