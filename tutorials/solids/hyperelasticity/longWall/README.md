---
sort: 1
---

# Tutorial: `longWall`

---

## Tutorial Aims

- Demonstrate how to perform a solid-only analysis in solids4foam;
- Exemplify the use of a hyperelastic mechanical law with large deformations.


## Case overview

<img src="images/long-wall.png" width="400" />

**Figure 1: Deformation of a long wall under two forces**

This benchmark consists in a wall with *1 m* x *2 m* (width x height) and
assumed infinitely long in length, i.e. along the z-axis (see Figure 1). Due to
symmetry, it can be assumed to be in a biaxial state of deformation. A
compressing load of *50 MPa* is applied on the right surface and a traction of
*100 MPa* on the upper surface. The left and bottom surfaces are free to slide
along its tangential, but constrained to move along their normal. The material
was assumed isotropic and incompressible and charaterized by the Mooney-Rivlin
hyperelastic law, with material parameters *c<sub>10</sub> = 80 MPa*,
*c<sub>01</sub> = 20 MPa*, and *c<sub>11</sub> = 0.0 MPa*, under plane-strain
conditions (see the file `constant/mechanicalProperties` to change these
values).

The default mesh in the tutorial had 20 x 20 volumes. You can use both total
and updated Lagrangian approaches in this tutorial because both yielded the
same results in terms of accuracy and convergence speed. For the results
presented here, we employed the `nonLinearGeometryTotalLagrangian` model, one
of the total Lagrangian approaches available, that solves for the *increment*
of the displacement, which is why you hav to use the `DD` file in the directory
`0/`.

The boundary conditions were applied in 100 equal incremental steps and note
the use of the `fixedDisplacementZeroShear` boundary condition applied to the
bottom and left surfaces of the wall that must be allowed to slip freely but
consttrained in their normal direction. A normalized residual tolerance for the
momentum equation of *10<sup>-8</sup>* was used.

---

## Expected results

The results are compared in Table 1 againt two references: one numerical and
another analytical provided by [I. Bijelonja, I. Demirdžić, and S. Muzaferija,
“A finite volume method for large strain analysis of incompressible
hyperelastic materials,” International Journal for Numerical Methods in
Engineering, vol.  64, pp. 1594–1609, Nov. 2005, doi: 10.1002/nme.1413.
](https://hrcak.srce.hr/206941). The *σ<sub>xx</sub>* and *σ<sub>yy</sub>* are
the Cauchy stress components along the *x* and *y* directions, respectively,
while *u<sub>x</sub>* and *u<sub>y</sub>* are the total displacements along the
x and y directions.

**Table 1: COmparison of Cauchy stress and the displacement of the wall in the x and y directions.**

| Source      | *σ<sub>xx</sub>* (MPa)| *σ<sub>yy</sub>* (MPa)|  *u<sub>x</sub>* (m)| *u<sub>y</sub>* (m)|
|-------------|-----------------------|-----------------------|---------------------|--------------------|
| solids4foam | -49.99                | 100.0                 | -0.1636             | 0.4010             |
| Reference   | -49.00                | 100.0                 | -0.1676             | 0.4022             |
| Analytical  | -50.00                | 100.0                 | -0.1675             | 0.4025             |

