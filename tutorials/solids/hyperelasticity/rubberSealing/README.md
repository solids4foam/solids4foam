---
sort: 16
---

# Compression of a Rubber Seal: `rubberSealing`

Prepared by Ivan Batistić

---

## Tutorial Aims

- Demonstrate how to perform a finite-strain hyperelastic solid-only analysis
  of a nearly incompressible rubber component in solids4foam.
- Demonstrate the updated Lagrangian solid model
  (`nonLinearGeometryUpdatedLagrangian`) in a displacement-controlled
  compression problem involving large deformations and rotations.
- Demonstrate the use of the `planeStress` option, the smoothed hydrostatic
  pressure equation (`solvePressureEqn`), and the `solidForces` and
  `solidPointDisplacement` function objects.

## Case Overview

This case examines the compression of a rubber seal, a benchmark problem
considered by Brink and Stein [1], Angoshtari et al. [2] and Pascon [3]. The
seal is a hollow trapezoid with 2 mm thick upper and lower flanges and
inclined side walls, as shown in Figure 1: the base is 17.47 mm wide, the
seal is 9 mm tall, and the trapezoidal cavity is 5 mm tall with widths of
11 mm (bottom) and 6 mm (top). To exploit symmetry, only half of the seal is
modelled, and the problem is solved as 2-D. Gravitational effects are
neglected and there are no body forces.

The bottom surface is fixed (`fixedDisplacement` of zero), while the upper
surface is given a prescribed displacement of 2.2 mm downwards, with its
horizontal displacement constrained. The displacement is applied in 100 equal
increments, using the displacement series in `constant/timeVsDisp`. The
remaining surfaces, including the surfaces of the cavity, are traction-free,
and the `solidSymmetry` condition is applied on the symmetry plane.

In the literature [1, 2], the rubber is described by a compressible
third-order Ogden model with $$\mu_1 = 0.63$$ MPa, $$\mu_2 = 0.0012$$ MPa,
$$\mu_3 = -0.01$$ MPa, $$\alpha_1 = 1.3$$, $$\alpha_2 = 5$$,
$$\alpha_3 = -2$$ and a bulk modulus of $$\kappa = 1000$$ MPa, under plane
strain conditions. This tutorial departs from that setup in two places:

1. **A neo-Hookean law is used instead of the Ogden law.** The material is
   modelled with `neoHookeanElastic`, with $$E = 1.2673215$$ MPa and
   $$\nu = 0.4997887797$$. These values give the initial shear modulus of the
   Ogden model, $$\mu = \frac{1}{2}\sum_p \mu_p \alpha_p = 0.4225$$ MPa, and a
   3-D bulk modulus of $$\kappa = 1000$$ MPa. The Ogden law (`OgdenElastic`)
   is not used because it currently supports neither the plane stress
   approximation nor the mixed pressure-displacement formulation; its
   parameters are kept, commented out, in `constant/mechanicalProperties`.
2. **Plane stress is assumed instead of plane strain** (`planeStress yes`).
   With a nearly incompressible material under plane strain conditions, the
   solution fails at about 50% of the stroke, independently of the
   formulation (displacement or mixed pressure-displacement), the load step,
   the relaxation, the pressure smoothing and the mesh. The iteration count
   grows smoothly towards this point whatever the step size, indicating a
   plane-strain structural instability rather than a numerical one. The
   reference results of Pascon [3] are also for plane stress.

```note
The `planeStress yes` option in solids4foam does not solve a finite-strain
plane stress problem. Instead, `neoHookeanElastic` replaces the bulk modulus
by its linear plane stress equivalent,
$$K = \frac{\nu E}{(1 + \nu)(1 - \nu)} + \frac{2}{3}\mu$$,
which is exact for small strains and an approximation for finite strains.
Here, this gives an effective in-plane bulk modulus of approximately
1.13 MPa, so the 3-D bulk modulus of 1000 MPa does not appear directly in the
computation. In addition, the out-of-plane stress is not constrained to be
zero, so it contributes to the reported equivalent stress.
```

The Ogden parameters are also subject to a known ambiguity. In the Ogden
convention used by solids4foam, $$W = \sum_p \frac{\mu_p}{\alpha_p}
(\lambda_1^{\alpha_p} + \lambda_2^{\alpha_p} + \lambda_3^{\alpha_p} - 3)$$,
they give an initial shear modulus of 0.4225 MPa; in the other common
convention, $$W = \sum_p \frac{2\mu_p}{\alpha_p^2}(\ldots)$$, the same
parameters give $$\sum_p \mu_p = 0.621$$ MPa, which is close to Pascon's
value of $$2c_1 = 0.625$$ MPa, so the literature parameters may have been
intended for that convention. This tutorial follows the solids4foam
`OgdenElastic` convention and uses 0.4225 MPa.

The hydrostatic stress is smoothed by solving a pressure equation
(`solvePressureEqn yes`) with `pressureSmoothingScaleFactor 1e4`; the default
value of 100 fails to converge at about 90% of the stroke. The problem is
solved as static using the segregated updated Lagrangian
`nonLinearGeometryUpdatedLagrangian` solid model, with the
`solutionTolerance` and `alternativeTolerance` settings given in
`constant/solidProperties`. The mesh, created with `blockMesh`, consists of
258 hexahedral cells in five blocks.

![Figure 1: Problem geometry [2]](./images/rubberSeal-geometry.png)

Figure 1: Problem geometry [2]

---

## Running the Case

The tutorial case is located at
`solids4foam/tutorials/solids/hyperelasticity/rubberSealing`. The case can be
run using the included `Allrun` script, i.e. `./Allrun`. The `Allrun` script
first creates the mesh using `blockMesh` and then runs the `solids4Foam`
solver. The case takes approximately 25 seconds to run in serial.

```warning
This case currently only runs with foam-extend. When run with OpenFOAM.com or
OpenFOAM.org, the `Allrun` script exits without running the case.
```

---

## Expected Results

As the upper surface is pushed down, the inclined side wall of the seal
rotates and bends, and the cavity closes progressively. The equivalent (von
Mises) stress is highest in the bending inclined wall, while the flanges away
from the wall remain lightly stressed. As noted by Pascon [3], the complex
displacement and stress fields near the corners of the seal can lead to
numerical instabilities.

The reaction force on the upper surface is written by the `solidForces`
function object to `postProcessing/0/solidForcestop.dat`, and the displacement
of the outer upper corner of the inclined wall, at (5.235, 7) mm in the
undeformed configuration, is written by the `solidPointDisplacement` function
object to `postProcessing/0/solidPointDisplacement_pointDisp.dat`:

```c++
functions
{
    forceTop
    {
        type            solidForces;
        historyPatch    top;
    }

    pointDisp
    {
        type            solidPointDisplacement;
        point           (0.005235 0.007 0);
    }
}
```

At the end of the loading (2.2 mm), the tutorial predicts (using
foam-extend-4.1) a vertical force of 0.171 N on the upper surface of the half
model (per 1 mm depth), a displacement magnitude of 0.975 mm at the monitored
point, and a maximum cell equivalent stress of approximately 904 kPa. Halving
the load step changes the force and the peak equivalent stress by less than
1%.

The peak equivalent stress is within approximately 2.5% of the value of
882 kPa reported by Pascon [3]. This is encouraging agreement rather than a
like-for-like validation: Pascon [3] uses an incompressible Yeoh-type model
with a different initial shear modulus ($$2c_1 = 0.625$$ MPa, compared with
0.4225 MPa here) and a true finite-strain plane stress formulation, whereas
this tutorial uses the plane stress approximation described above.

Figure 2 shows the deformed seal coloured by the equivalent stress, where the
left image is from Pascon [3] and the right image is a solids4foam result
from the original benchmark study, whose settings were not recorded (its peak
equivalent stress is 836 kPa). These images are included for a qualitative
comparison of the deformed shape and the stress distribution only; no
quantitative reference data are provided with this tutorial.

![Figure 2: Equivalent stress (in kPa) in the deformed seal: Pascon [3] (left) and solids4foam from the original benchmark study (right)](./images/rubberSeal-sigmaEq.jpeg)

Figure 2: Equivalent stress (in kPa) in the deformed seal: Pascon [3] (left)
and solids4foam from the original benchmark study, with unrecorded settings
(right)

---

## References

[1]
[U. Brink and E. Stein, A posteriori error estimation in large-strain
elasticity using equilibrated local Neumann problems. _Computer Methods in
Applied Mechanics and Engineering_, 161(1-2), 77-101,
1998.](https://www.sciencedirect.com/science/article/abs/pii/S0045782597003101)

[2]
[A. Angoshtari, M. Faghih Shojaei and A. Yavari, Compatible-strain mixed
finite element methods for 2D compressible nonlinear elasticity. _Computer
Methods in Applied Mechanics and Engineering_, 313, 596-631,
2017.](https://www.sciencedirect.com/science/article/abs/pii/S0045782516312798)

[3]
[J. P. Pascon, Large deformation analysis of plane-stress hyperelastic
problems via triangular membrane finite elements. _International Journal of
Advanced Structural Engineering_, 11(3), 331-350,
2019.](https://link.springer.com/article/10.1007/s40091-019-00234-w)
