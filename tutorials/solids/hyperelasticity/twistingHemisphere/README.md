---
sort: 12
---

# Twisting contact between a hemisphere and a block: `twistingHemisphere`

---

Prepared by Ivan Batistić

---

## Tutorial Aims

- Demonstrate how to perform a solid analysis with a hyperelastic material,
  large deformations and large sliding frictional contact between two
  deformable bodies in 3-D;
- Show how to prescribe a combined translation and rotation of a boundary with
  the `fixedRotation` condition;
- Compare the vertical contact force and the twisting torque with published
  solutions.

---

## Case Overview

A thick-walled hollow hemisphere is pressed into a deformable block and then
twisted about its vertical axis (Figure 1). The block is a cube of side
$$L = 2$$ fixed on its bottom face, and the hemisphere has an outer radius of
$$R_o = L/2 = 1$$. Both bodies are described by the `neoHookeanElastic` law with
Poisson's ratio $$\nu = 0.3$$; the hemisphere has Young's modulus $$E = 5$$ and
the block $$E = 1$$ (see `constant/mechanicalProperties`). Coulomb friction with
a high coefficient of friction, $$\mu = 0.5$$, acts between the bodies, and the
`standardPenalty` normal and friction contact models are used. Gravitational and
inertial forces are neglected, and the updated Lagrangian formulation
(`nonLinearGeometryUpdatedLagrangian`) is used.

The top surface of the hemisphere (patch `sphere-displacement`) is loaded with
the `fixedRotation` condition in two phases:

1. indentation: a downward displacement of $$u_y = 1$$
   (`constant/timeVsDisp`) is applied from $$t = 0$$ to $$t = 4$$ in 40 equal
   increments;
2. twisting: with the indentation held, a rigid rotation of $$180^{\circ}$$
   about the $$y$$ axis (`constant/rotVsDisp`) is applied from $$t = 4$$ to
   $$t = 13$$ in 90 equal increments, i.e. $$2^{\circ}$$ per increment.

The rotation angle is therefore $$\theta = 20(t - 4)$$ degrees during the
twisting phase.

```note
The original benchmark prescribes the indentation in 10 increments. With the
current solids4foam version, an increment of 0.1 makes the first momentum
corrector of the first time step produce a non-positive Jacobian next to the
loaded patch, and the solver fails; the tutorial therefore uses 40 indentation
increments. In the original benchmark, and in [2] and [3], the indentation is
also frictionless; this option is not available in solids4foam, so friction is
active from the start of the simulation.
```

The tutorial ships the coarse mesh of the original benchmark (2 349 cells: 1 620
in the hemisphere and 729 in the block) as the compressed archive
`constant/polyMesh_coarse.tar.xz`. The archived mesh is stored with
$$L = 0.002$$, so `Allrun` scales it by 1000 with `transformPoints`. On this
mesh, the block patch `block-contact` is the master contact patch; with the
hemisphere as master, as in the original benchmark, the momentum loop stalls
from about half-way through the indentation. For the same reason, the friction
penalty scale factor is reduced from 1.5 to 0.5. The finer meshes used in [1]
are available in the `solid-benchmarks` repository
(`hyperElasticity/twistingHemisphere/meshes`).

![Figure 1: Problem geometry and computational mesh [1]](./images/twistingHemisphere-geometry.png)

Figure 1: a) Problem geometry and material parameters, in the notation of
[2]; b) the 35 424-cell mesh of [1] (block 27 744 cells, hemisphere 7 680
cells). The material parameters used in the tutorial are those given above
and in `constant/mechanicalProperties`.

```warning
This case currently only runs with foam-extend.
```

---

## Running the Case

The tutorial case is located at
`solids4foam/tutorials/solids/hyperelasticity/twistingHemisphere`. The case can
be run using the included `Allrun` script, i.e. `> ./Allrun`. The `Allrun`
script extracts the mesh from `constant/polyMesh_coarse.tar.xz`, scales it with
`transformPoints` and runs the `solids4Foam` solver. Optionally, if `gnuplot` is
installed, the vertical force and the twisting torque on the
`sphere-displacement` patch are plotted against the reference data in the
`force.png` and `torque.png` files.

With foam-extend-4.1, the complete run (130 increments) takes about 13 minutes
on one core. In 20 increments, between $$t = 3.8$$ and $$t = 5.7$$ (the end of
the indentation and the start of the twisting, when the contact changes from
stick to slip), the momentum loop reaches its limit of 2 000 correctors before
meeting the tolerance; all other increments converge.

The `regressionTest.sh` script runs a shortened version of the case to keep the
test short: it stops at $$t = 6$$, i.e. after the full indentation and the first
$$40^{\circ}$$ of the twist (about 6 minutes). Its checks are made only at
increments whose momentum loop converged, because values from the
non-converged stick-to-slip increments ($$t = 3.8$$ to $$5.7$$) depend on the
iteration path and are not reliable regression targets. The script confirms
from the solver log that the increments at $$t = 3.5$$ and $$t = 6$$ converged,
and checks:

- the vertical force at $$t = 3.5$$ (indentation) and the vertical force and
  twisting torque at $$t = 6$$ ($$40^{\circ}$$) against their calibrated
  values;
- the vertical force and twisting torque at $$40^{\circ}$$ against the solution
  of [3], interpolated from `deLorenzisForce.dat` and `deLorenzisMoment.dat`
  (to within 3% and 25%, respectively; the calibration run gives deviations of
  0.8% and 18%).

The twisting phase beyond $$40^{\circ}$$ is not part of the regression test.

---

## Expected Results

The high coefficient of friction initially makes the contact stick, so the
block is twisted with the hemisphere (Figure 2). As the rotation increases, the
contact area gradually changes to slip, after which an almost constant torque is
transmitted to the block, and the block deformation no longer changes.

![Figure 2: Deformed mesh after compression and at 30, 60 and 180 degrees of rotation [1]](./images/twistingHemisphere-geometryEvolution.png)

Figure 2: Deformed mesh after compression and at 30, 60 and 180 degrees of
rotation, on the 35 424-cell mesh of [1]

Table 1 compares the vertical force and the twisting torque on the
`sphere-displacement` patch computed with the tutorial (coarse mesh,
foam-extend-4.1) with the digitised solutions of [3] (`deLorenzisForce.dat`,
`deLorenzisMoment.dat`) and [2] (`febioMoment.dat`, divided by 10 to use the
same scale).

**Table 1: Vertical force and twisting torque on the coarse mesh.**

| Angle | Force, solids4foam | Force [3] | Torque, solids4foam | Torque [3] | Torque [2] |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 0° | 1.282 | 1.230 | 0.000 | 0.002 | 0.003 |
| 20° | 1.241 | 1.234 | 0.139 | 0.146 | 0.144 |
| 40° | 1.235 | 1.245 | 0.219 | 0.267 | 0.273 |
| 80° | 1.241 | 1.273 | 0.313 | 0.395 | 0.395 |
| 120° | 1.249 | 1.282 | 0.346 | 0.403 | 0.410 |
| 180° | 1.254 | 1.289 | 0.361 | 0.401 | 0.405 |

Over the twisting phase, the vertical force agrees with [3] to 2.4% RMS (4.3% at
most, at $$0^{\circ}$$). The torque is higher than the references at small
angles (16% at $$10^{\circ}$$) and agrees with them to about 4% at
$$20^{\circ}$$, but on the coarse mesh the transition from stick to slip is
more gradual, so the torque is up to 22% lower around $$60^{\circ}$$ and is
10% lower than [3] at $$180^{\circ}$$, where it is still slowly increasing.

Figures 3 and 4 show the corresponding results obtained in [1] on the
35 424-cell mesh. On that mesh the twisting torque agrees well with [2] and [3]
during the transition from stick to slip, and is slightly overpredicted in the
slip regime because of the higher normal contact force. Note that in [2] and [3]
the indentation is frictionless, which affects the results. The plots in
Figures 3 and 4 show the force and torque multiplied by 10.

![Figure 3: Evolution of the vertical contact force [1]](./images/twistingHemisphere-force.png)

Figure 3: Evolution of the vertical contact force on the 35 424-cell mesh [1]

![Figure 4: Evolution of the twisting moment [1]](./images/twistingHemisphere-torque.png)

Figure 4: Evolution of the twisting moment on the 35 424-cell mesh [1]

The results of [2] and [3] were digitised using
[WebPlotDigitizer](https://apps.automeris.io/wpd/).

---

### References

[1]
[I. Batistić, Segment-to-Segment Algorithm for Finite Volume Mechanical Contact
Simulations, PhD thesis, University of Zagreb,
2022.](https://www.sciencedirect.com/science/article/abs/pii/S0307904X21004248)

[2]
[B. K. Zimmerman and G. A. Ateshian, “A surface-to-surface finite element
algorithm for large deformation frictional contact in FEBio,” Journal of
Biomechanical Engineering, vol. 140, no. 8,
2018.](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6056201/)

[3]
[R. A. Sauer and L. De Lorenzis, “An unbiased computational contact formulation
for 3D friction,” International Journal for Numerical Methods in Engineering,
vol. 101, no. 4, pp. 251–280, 2015.](https://doi.org/10.1002/nme.4794)
