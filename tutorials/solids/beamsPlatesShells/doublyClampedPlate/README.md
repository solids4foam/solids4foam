---
sort: 3
---

# Plate clamped on two opposite edges: `doublyClampedPlate`

---

Prepared by Ivan Batistić

---

## Tutorial Aims

- Demonstrate the free-edge boundary conditions of the finite area
  Kirchhoff-Love plate model on a plate that is clamped along two opposite
  edges and free along the other two;
- Show how the predicted deflection converges with mesh refinement and measure
  the order of convergence against a converged finite element solution.

---

## Case Overview

The plate is the same as in the `squarePlate` and `cantileverPlate` tutorials:
$$L \times L \times h = 10$$ m $$\times~10$$ m $$\times~0.1$$ m, loaded by a
uniform transverse pressure $$p = 1000$$ Pa, with Young's modulus $$E = 200$$
GPa and Poisson's ratio $$\nu = 0.3$$. The `left` ($$x = 0$$) and `right`
($$x = 10$$ m) edges are clamped and the `top` and `bottom` edges are free. The
density is set to zero so that the single time step gives the static solution
(see the `cantileverPlate` tutorial for the reason).

The boundary conditions are the same as in the `cantileverPlate` tutorial:
`fixedValue` and `clampedMoment` on the clamped edges, `freeEdgeDisplacement`
and `freeEdgeMoment` on the free edges. Here every corner joins a clamped and
a free edge, so the free-free corner treatment of the conditions is not
exercised, but the clamped moment condition is now applied on two edges that
both drive the explicit corrector loop.

---

## Expected Results

The plate behaves like a wide clamped-clamped beam whose free side edges
deflect more than its centre line. The reference is a converged Kirchhoff plate
finite element solution of the same problem (12-DOF ACM rectangular element,
Richardson-extrapolated from 64 x 64 and 128 x 128 elements), which gives the
deflection at the midpoint of a free edge as $$w = 0.0015882$$ m and at the
plate centre as $$0.0013976$$ m. For comparison, the Euler-Bernoulli
clamped-clamped strip estimate $$p L^4 / (384 D)$$ is $$0.0014219$$ m, which
lies between the two because the strip neither relaxes the anticlastic
curvature at the free edges nor carries the full plate stiffness at the
centre.

The tutorial mesh of 50 x 50 cells predicts a free-edge midpoint deflection of
$$0.001605$$ m, 1.1% above the reference, and a centre deflection of
$$0.001402$$ m, 0.3% above; both overpredictions fall by a factor of two to
three with every mesh refinement (see the verification study below).

---

## Running the Case

The tutorial case is located at
`solids4foam/tutorials/solids/beamsPlatesShells/doublyClampedPlate`. The case
can be run using the included `Allrun` script, i.e. `> ./Allrun`, which
creates the mesh using `blockMesh` (`> blockMesh`), creates the finite area
mesh using `makeFaMesh` (`> makeFaMesh`), and runs the case with the
`solids4foam` solver (`> solids4Foam`). The case runs with OpenFOAM.com v2412
or newer and takes a few seconds on the 50 x 50 mesh.

Setting `coupled yes;` in `constant/solidProperties` solves the `M` and `w`
equations with their edge conditions as one linear system with PETSc instead
of the segregated corrector loop; see the `cantileverPlate` tutorial. It gives
the same solution in two outer iterations and requires solids4foam to be
compiled with PETSc.

---

## Verification and Convergence Study

The opt-in [`verification/`](verification/) directory adds the mesh-convergence
sweep. It scales the two in-plane divisions by powers of two from $$25 \times
25$$ to $$200 \times 200$$ cells, reads the deflection at the midpoint of the
free `top` edge and at the plate centre from the written finite area field,
compares them with the finite element reference and reports the observed order
of convergence:

```bash
cd verification
./Allverify
```

The deflections converge towards the finite element reference at an observed
order between one and two (relative error against the reference):

| Cells | Free-edge midpoint | Error | Centre | Error |
|:---:|---:|---:|---:|---:|
| 25 x 25 | 0.0016397 m | +3.2% | 0.0014197 m | +1.6% |
| 50 x 50 | 0.0016051 m | +1.1% | 0.0014017 m | +0.3% |
| 100 x 100 | 0.0015951 m | +0.43% | 0.0013985 m | +0.06% |
| 200 x 200 | 0.0015912 m | +0.19% | 0.0013978 m | +0.01% |

The net order of the free-edge midpoint deflection is 1.36 (successive orders
1.61, 1.30 and 1.18, falling towards the first order of the one-sided free
edge estimates), and Richardson extrapolation of the three finest meshes
gives an order of 1.39 and a limit of $$0.0015889$$ m, 0.04% above the finite
element reference, so the scheme converges to the right answer. The whole
sweep takes about an hour in serial on an Apple M1 Ultra with the segregated
solver, almost all of it on the finest mesh, and a few seconds with
`./Allverify --coupled`.

The study is separate from `regressionTest.sh` and is not run by the normal
tutorial tests. See the verification README for the mesh levels, options and
acceptance criteria.

![Figure 1: Relative error of the deflection against the mesh
spacing](./images/doublyClampedPlate-convergence.png)

Figure 1: Relative error of the deflection against the mesh spacing.

---

### References

[1] [Timoshenko, S., & Woinowsky-Krieger, S., Theory of plates and shells.
1959.](https://www.cap-recifal.com/ccs_files/articles/cuveaqua1_denisio/Timoshenko_-_Theory_of_plates_and_shells.pdf)

[2] Adini, A., & Clough, R. W., Analysis of plate bending by the finite
element method. NSF Report G7337, 1960 (the ACM rectangular plate element
used for the reference solution).
