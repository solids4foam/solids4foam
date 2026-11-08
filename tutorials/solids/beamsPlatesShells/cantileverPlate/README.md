---
sort: 2
---

# Cantilever plate with free edges: `cantileverPlate`

---

Prepared by Ivan Batistić

---

## Tutorial Aims

- Demonstrate the free-edge boundary conditions of the finite area
  Kirchhoff-Love plate model on a plate that is clamped along one edge and
  free along the other three;
- Show how the predicted tip deflection converges with mesh refinement and
  measure the order of convergence against a converged finite element
  solution.

---

## Case Overview

The plate is the same as in the `squarePlate` tutorial: $$L \times L \times h =
10$$ m $$\times~10$$ m $$\times~0.1$$ m, loaded by a uniform transverse pressure
$$p = 1000$$ Pa, with Young's modulus $$E = 200$$ GPa and Poisson's ratio
$$\nu = 0.3$$. The `left` edge ($$x = 0$$) is clamped (zero deflection and zero
rotation) and the `right`, `top` and `bottom` edges are free (zero bending
moment and zero Kirchhoff shear force). The plate weight is ignored and the
density is set to zero in `constant/mechanicalProperties`, so the single time
step gives the static solution: with a non-zero density the solver's inertia
term would add $$\rho h w / \Delta t^2$$ to the load, about 2% of $$p$$ here.

The plate model solves the two coupled second-order equations

$$
\nabla^2 M = -p, \qquad M = -D \nabla^2 w,
$$

where $$w$$ is the transverse displacement, $$M$$ the moment sum and $$D$$ the
bending stiffness. On a free edge with outward normal $$n$$ and tangent $$s$$
the bending moment $$M_n$$ and the Kirchhoff shear force $$V_n = Q_n + \partial
M_{ns}/\partial s$$ vanish, which the model imposes as a pair of conditions on
the two unknowns:

- `freeEdgeDisplacement` (on `w`): the normal curvature is
  $$\partial^2 w / \partial n^2 = \dfrac{\nu}{1 - \nu} \dfrac{M}{D}$$, imposed
  as a normal gradient of $$w$$ extrapolated from the adjacent cell;
- `freeEdgeMoment` (on `M`): the normal gradient of the moment sum is
  $$\partial M / \partial n = D (1 - \nu)\, \partial^2 (\partial w / \partial n)
  / \partial s^2$$, which couples the moment to the second derivative of the
  normal slope along the edge. At a corner shared by two free edges the
  twisting moment vanishes, so the normal slope is mirrored about the corner;
  at a corner with a clamped edge both slopes vanish there.

The clamped edge uses `fixedValue` for `w` and `clampedMoment` for `M`, which
sets the edge moment from the normal curvature estimated as the slope of the
adjacent cell centre over its distance $$d$$ to the edge, where the slope
vanishes: $$\partial^2 w / \partial n^2 = -(n \cdot \theta_P) / d$$ with
$$\theta_P = -\nabla w$$ at the cell centre. The one-sided combination
$$6 w_P / d^2 - 2 (n \cdot \theta_P) / d$$, formally of higher order, was found
to converge only at first order, while this estimate converges at second order.

The conditions are coupled explicitly, so the `relaxationFactor` of each
condition and the field under-relaxation of `w` and `M` in
`system/finite-area/faSolution` are needed for the corrector loop to converge.
The loop gain of `clampedMoment` grows with mesh refinement, so finer meshes
need a smaller factor than the `0.01` used here (the verification study uses
`0.002` on 100 x 100 cells and `0.001` on 200 x 200 cells).

---

## Expected Results

There is no closed-form Kirchhoff solution for a plate with three free edges.
The reference is a converged Kirchhoff plate finite element solution of the
same problem (12-DOF ACM rectangular element, Richardson-extrapolated from
64 x 64 and 128 x 128 elements), which gives the deflection at the midpoint
of the free tip edge as $$w = 0.070475$$ m, the mean deflection along that edge
as $$0.070121$$ m, and the corner deflection as $$0.069470$$ m. For comparison,
the Euler-Bernoulli strip estimate $$p L^4 / (8 D)$$ is $$0.06825$$ m: the plate
deflects more than a beam strip because the free side edges relax the
anticlastic curvature that a strip would suppress.

The tutorial mesh of 50 x 50 cells predicts a tip midpoint deflection of
$$0.07057$$ m, 0.13% above the reference; the overprediction falls by a factor
of two to three with every mesh refinement (see the verification study
below). The bending moment along
the free side edges and its normal gradient also agree with the finite element
solution. Figure 1 shows the deflection field.

![Figure 1: Transverse deflection of the cantilever
plate](./images/cantileverPlate-deflection.png)

Figure 1: Transverse deflection of the cantilever plate (50 x 50 cells).

---

## Running the Case

The tutorial case is located at
`solids4foam/tutorials/solids/beamsPlatesShells/cantileverPlate`. The case can
be run using the included `Allrun` script, i.e. `> ./Allrun`, which creates
the mesh using `blockMesh` (`> blockMesh`), creates the finite area mesh using
`makeFaMesh` (`> makeFaMesh`), and runs the case with the `solids4foam` solver
(`> solids4Foam`). The case runs with OpenFOAM.com v2412 or newer and takes
about 20 s on the 50 x 50 mesh.

The segregated algorithm couples the clamped and free edge conditions
explicitly and needs thousands of under-relaxed correctors. Setting
`coupled yes;` in `constant/solidProperties` instead assembles the `M` and `w`
equations with their edge conditions into one linear system and solves it with
PETSc (a direct solve with MUMPS by default; a `petscOptions` file in the case
directory can select another solver). This gives the same solution in two
outer iterations and a fraction of a second, and requires solids4foam to be
compiled with PETSc.

---

## Verification and Convergence Study

The opt-in [`verification/`](verification/) directory adds the mesh-convergence
sweep. It scales the two in-plane divisions by powers of two from $$25 \times
25$$ to $$200 \times 200$$ cells, reads the deflection along the free tip edge
from the written finite area field, compares it with the finite element
reference and reports the observed order of convergence:

```bash
cd verification
./Allverify
```

The tip deflection converges towards the finite element reference at an
observed order between one and two (relative error against the reference):

| Cells | Tip midpoint | Error | Tip mean | Error |
|:---:|---:|---:|---:|---:|
| 25 x 25 | 0.070734 m | +0.37% | 0.070387 m | +0.38% |
| 50 x 50 | 0.070567 m | +0.13% | 0.070217 m | +0.14% |
| 100 x 100 | 0.070511 m | +0.05% | 0.070160 m | +0.06% |
| 200 x 200 | 0.070490 m | +0.02% | 0.070138 m | +0.02% |

The net order of the tip midpoint deflection is 1.35 (successive orders 1.50,
1.34 and 1.22: the clamped edge treatment is second order while the one-sided
free edge estimates are first order, so the observed order falls towards one
as the free-edge error takes over), and Richardson extrapolation of the three
finest meshes gives an order of 1.43 and a limit of $$0.070478$$ m, 0.005%
above the finite element reference, so the scheme converges to the right
answer. The whole sweep takes about an hour in
serial on an Apple M1 Ultra with the segregated solver, most of it on the
finest mesh, and under a minute with `./Allverify --coupled`.

The study is separate from `regressionTest.sh` and is not run by the normal
tutorial tests. See the verification README for the mesh levels, options and
acceptance criteria.

![Figure 2: Relative error of the tip deflection against the mesh
spacing](./images/cantileverPlate-convergence.png)

Figure 2: Relative error of the tip deflection against the mesh spacing.

---

### References

[1] [Timoshenko, S., & Woinowsky-Krieger, S., Theory of plates and shells.
1959.](https://www.cap-recifal.com/ccs_files/articles/cuveaqua1_denisio/Timoshenko_-_Theory_of_plates_and_shells.pdf)

[2] Adini, A., & Clough, R. W., Analysis of plate bending by the finite
element method. NSF Report G7337, 1960 (the ACM rectangular plate element
used for the reference solution).
