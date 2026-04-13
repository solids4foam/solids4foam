# Flow over a rectangular building with a flexible membrane roof: `membraneRoof`

Prepared by Ivan Batistić, Philip Cardiff

## Tutorial Aims

- Demonstrate a three-dimensional fluid-structure interaction benchmark with a
  very slender membrane roof.
- Show the partitioned Dirichlet-Neumann coupling and the monolithic fluid-solid interaction
  solver.

## Case Overview

This tutorial studies a rectangular building covered by a flexible
$10 \times 10$ m membrane roof and exposed to an incoming three-dimensional flow; see Figure 1.
The setup follows the membrane-roof example described in Section 5.1 of the
paper by von Scheven and Ramm [1], which extends the earlier two-dimensional benchmark
from [2]. With a roof length of $10$ m, the length-to-thickness ratio is $1000$, so the structure is
extremely slender.

![Fluid and solid geometry together with the boundary-condition layout](./images/membraneRoof-geometry.png)

**Figure 1:** Geometry and boundary conditions for the membrane-roof benchmark.

The physical parameters of the problem are:

- **Geometry**
  - Building base: $10 \times 10$ m
  - Building height: $5$ m
  - Overall domain size ($L \times W \times H$): $150 \times 100 \times 75$ m
  - Membrane thickness: $0.01$ m
  - Gravity ($g$): $0$ $\mathrm{m/s}^2$
- **Solid**
  - St. Venant-Kirchhoff material model
  - Density ($\rho_s$): $1000$ $\mathrm{kg/m}^3$
  - Young’s modulus ($E_s$): $1\cdot10^9$ Pa
  - Poisson’s ratio ($\nu_s$): $0$
- **Fluid**
  - Laminar flow
  - Density ($\rho_f$): $1.25$ $\mathrm{kg/m}^3$
  - Kinematic viscosity ($\nu_f$): $0.08$ $\mathrm{m}^2$/s

As in the reference paper [1], the case is solved using $600$ uniform time steps of $\Delta t = 0.02$ s.
The published study reports that the maximum inlet speed is $71.26$ m/s, reached at $z = 75$ m and
$t = 5$ s, corresponding to $Re \approx 8900$ when the roof length is used as the characteristic length.

At the inlet boundary ($x = 0$), the paper prescribes the streamwise velocity $u_x(z, t) = 100 \hat{u}_t(t) \hat{u}_z(z).$

The spatial and temporal factors are given by (see Figure 2):

- $\hat{u}_z(z) = (z/350)^{0.22}$
- $\hat{u}_t(t) = 0.5 [ sin(\pi (t/5 - 0.5)) + 1 ]$ for $t < 5$ s and $\hat{u}_t(t) = 1$ for $t >= 5$ s

<img src="./images/membraneRoof-inlet.png" alt="Spatial (left) and temporal (right) parts of the inlet velocity profile" style="zoom:25%;" />

**Figure 2:** Spatial and temporal variations of the inflow boundary condition.

The inlet condition is implemented as a tutorial-local library in
`src/membraneRoofVelocityFvPatchVectorField`
and loaded through `system/controlDict`. This avoids `codedFixedValue`, which
is not portable across all supported OpenFOAM forks.

The library is compiled automatically by `./Allrun` through `src/Allwmake`.

In this tutorial mesh, the coordinate corresponding to the paper's vertical
`z` direction is the case `y` direction. For this reason, the compiled inlet
boundary condition evaluates the profile using the patch-face `y` coordinate,
rather than the case `z` coordinate.



## Mesh Generation

The tutorial uses a coarser mesh intended for routine testing and demonstration.
The fluid mesh consists of $28 028$ cells, and the solid mesh consists of $1 536$ cells.

For the fluid region, the mesh is generated with `blockMesh` and then refined
locally above and around the roof using `setSet` and `refineMesh`. The `setSet`
command creates the `cellsToRefine` cell set from the selection specified in
`setSet.batch`, while `refineMesh`, using `system/fluid/refineMeshDict`,
uniformly splits each selected cell into eight cells; see Figure 4.
Figure 3 shows the block-structured layout used to generate the base fluid mesh.

![Base blockMesh layout used for the fluid region](./images/membraneRoof-blockMeshDict.png)

**Figure 3:** Base `blockMesh` layout for the fluid region.

![Fluid and solid geometry together with the boundary-condition layout](./images/membraneRoof-refinement.png)

**Figure 4:** Base `blockMesh` mesh (left) and mesh after refinement using `refineMesh` utility (right).

## Running the Case

From the case directory, run

```bash
./Allrun                      # partitioned (Aitken), serial
./Allrun monolithic           # monolithic (Newton/PETSc), serial
./Allrun monolithic parallel  # monolithic (Newton/PETSc), parallel
./Allrun parallel             # partitioned (Aitken), parallel
```

```note
The monolithic approach requires OpenFOAM.com (ESI) and a PETSc installation
with `PETSC_DIR` set. It does not run with foam-extend or OpenFOAM.org.
```

The `Allrun` script will:

1. Compile the local inlet boundary-condition library.
2. Link the partitioned or monolithic case files.
3. Generate the fluid and solid meshes.
4. Refine the fluid mesh.
5. Run `solids4Foam`.
6. Generate deflection and force plots with gnuplot (if available).

To clean the generated files, use

```bash
./Allclean
```

```note
For more information about the monolithic approach, as well as the Aitken and IQNILS approaches, see the other FSI tutorials.
```



## Results

The solution shows vortex shedding around the building and a transition of the membrane roof shape from concave to convex.

One useful quantity of interest is the membrane displacement at the roof centre point.
In this case it is recorded by the `solidPointDisplacement` function object in
`system/controlDict`.

```c++
    pointDisp
    {
        type                solidPointDisplacement;
        point               (50 5 0);
    }
```



![Fluid and solid geometry together with the boundary-condition layout](./images/membraneRoof-deflection.png)

**Figure 5:** Deflection of the roof centre point $(50, 5, 0)$.

Other useful outputs to check are the interface forces from the `forces` function object
and the number of fluid-structure coupling iterations per time step, which are also reported
using function objects in `system/controlDict`. These are plotted using a gnuplot script that
generates `deflection.pdf`, `force.pdf`, and `iterations.pdf`.

The regression script (`regressionTest.sh`) tracks the tip displacement in
`postProcessing/0/solidPointDisplacement_displacement.dat` and the total fluid
force in `postProcessing/fluid/forces/0/force.dat` at the final time ($t = 12$ s).

## References

[1] [M. von Scheven and E. Ramm, "Strong coupling schemes for interaction of
thin-walled structures and incompressible flows," *International Journal for
Numerical Methods in Engineering*, 87(1-5), 2011, pp. 214-231](
https://doi.org/10.1002/nme.3033)

[2] [P. Le Tallec and J. Mouro, "Fluid structure interaction with large
structural displacements," *Computer Methods in Applied Mechanics and
Engineering*, 190(24-25), 2001, pp. 3039-3067.](https://www.sciencedirect.com/science/article/abs/pii/S0045782500003819)
