---
sort: 4
---

# My fourth tutorial: `3dTube`

---

## Case overview

This case consists of a flow pulse moving in an elastic thick tube. It is used
to exemplify the application of different formulations and models that are
currently implemented in solids4foam to solve the fluid-solid interaction (FSI)
of the flow pulse with the tube deformation. A schematic figure of the case can
be seen in the figure below.

![](images/3dTube.png)

The fluid is assumed incompressible and Newtonian, with density 1000.0 kg/m3
and kinematic viscosity 3e-6 m2/s, flowing due to a pressure wave, with peak
1333.3 Pa, applied at the tube inlet during 3e-3s. Hence, the Naview-Stokes
equations must be solved, and they take the form:

The tube wall is assumed an isotropic elastic body under the small-strains
regime, modeled with Hooke's law (only for consistency with
the original publication that proposed this benchmark, although a linear
elastic law may suffice) with density 1200.0 kg/m3, Young's modulus 3e5 Pa and
Poisson's ratio 0.3.

The case exemplify a situation where the coupling between the fluid and the
solid is considered to be strong and this is mainly caused by the close value
between the density of the fluid and solid sub-domains. Rigorously, this
situation occurs when the ratio 'solid density/fluid density' is near to 1.0.
When this situation occurs and you are using a *partitioned strategy* to solve
the FSI coupling, i.e. where the fluid and solid sub-problems are solved
separately and the case in solids4foam, the so-called 'added-mass operator'
render the problem difficult to solve numerically due to instabilities in the
course of the solution. 

The solution to improve the convergence is to apply a coupling algorithm that
may estabilize the solution. The typical strongly-coupled Dirichlet-Neumann
coupling algorithm  is of the form:

```pseudocode
for all time-steps
    do
        solve fluid domain
        pass fluid interface forces to the solid interface
        solve solid domain
        pass solid interface velocities to the fluid
        interface using under-relaxation
        update the fluid mesh
    while not converged
end
```

For this tutorial, you can use two approaches that are currently implemented in
solids4foam, independently of the OpenFOAM versions you are using, and an extra
one in case you use foam-extend-4.1.

The firt two approaches use the typical *Dirichlet-Neumann formulation* by
which the condition on the FSI interface is modeled as follows: the solid
receives the traction force from the flow  (Neumann boundary condition) and the
fluid domain receives the velocity of the wall (Dirichlet boundary condition).
The difference on the two approaches used lies on the fluid model:

1) Using a purely incompressible formulation for the fluid flow: we can achieve
that with the pimpleFluid model, the 'newMovingWallVelocity' BC in the 'wall'
boundary for the flow velocity, and the 'zeroGradient' BC on the 'wall'
boundary for the pressure.  This configuration is the typical way of handling
FSI problems, but may be unstable even when iterating several times over the
fluid and solid solution, so we need to use a strong form of coupling
algorithm. Two options available are: 

(a) Aitken's dynamic relaxation and;
(b) the IQN-ILS algorithm. 
    
We label this alternative DNF for reference.

2) Using a 'weakly compressible' fluid model: the convergence of the FSI
coupling can be improved by using a compressible law to model the fluid and
controlling the compressibility by the bulk modulus included in the model.
Thus, the incompressible behavior could be achieve by increasing the bulk
modulus sufficiently, although this value is application-dependent. We can
achieve that in solids4Foam with the 'sonicLiquidFluid' fluid model and the
same BCs as used in item (1). Both Aitken's dynamic relaxation or the IQN-ILS
algorithms can be used for the FSI coupling, but as we will see below, the
convergence improvement when using this weakly compressible model is
impressive. In the tutorial, we employed the bulk modulus of the fluid of 1
MPa.

More information on this weakly compressible law and its effect in FSI
computations can be found at:
[E. Tandis and A. Ashrafizadeh, “A numerical study on the fluid compressibility effects in strongly coupled fluid–solid interaction problems,” Engineering with Computers, 2019, doi: 10.1007/s00366-019-00880-4.
](https://doi.org/10.1007/s00366-019-00880-4)

We label this alternative WCM for reference (from 'weakly compressible model').

3) The third alternative, that currently is only implemented when compiling
with foam-extend, uses a *Robin-Neumann formulation*, consisting in a different
way of writing the mathematical formulation of the FSI problem by taking
directly into account the added-mass operator. It leads to a Robin boundary
condition for the pressure field of the flow on the FSI interface.  Support for
it is implemented in the pimpleFluid model, but we need to apply special BCs
for this case to work: they are the 'elasticWallPressure' for the pressure
field and 'elasticWallVelocity' for the velocity field of the fluid flow
domain. 

We label this alternative RNF for reference (from 'Robin-Neumann formulation').

In all approaches, the solid domain setup is exactly the same. Hence, the
momentum and continuity equations are solved for the fluid flow and the
momentum equation for the solid wall, in this case, with an updated Lagrangian
formulation. Due to symmetry, only a quarter of the tube's cross section is
considered and the test will run for 0.02 s.

---

## Running the case

The tutorial case can be run using the included `Allrun` scripts:
`Allrun.pimpleFluid` for the approaches 1 and 3, i.e. `> ./Allrun.pimpleFluid`,
and `> ./Allrun.sonicLiquidFluid` for the option 2. The script will update the
case with links to the correct files to be used by each approach.

Additionally, if you will test the option (1), you can use the Aitken or
IQN-ILS algorithm, which should be changed in the `constant/fsiProperties` file
prior to running the script. 

If running the option (3), i.e. with the Robin boundary condition, you should
further change the boundary conditions in the `0/fluid/p` and `0/fluid/U` files
to the `elastic*` versions that are commented in those files. 

```tip
Remember that a tutorial case can be cleaned and reset using the included
`Allrun` script, i.e. `> ./Allclean`.
```

---

## Expected results

This case has been proposed as a benchmark for FSI problems and the solution
for the axial and radial displacement of the point A (see figure above) is
available in the following reference:

[A. Lozovskiy, M. A. Olshanskii, and Y. V. Vassilevski, “Analysis and assessment of a monolithic FSI finite element method,” Computers and Fluids, vol. 179, pp. 277–288, 2019, doi: 10.1016/j.compfluid.2018.11.004.](https://doi.org/10.1016/j.compfluid.2018.11.004)

We compare the solution for the different approaches in the two figures shown
below. We ran tests with the option (1) with both the Aitken and the IQN-ILS
algorithms.

![](./images/axial-displacement.png)
![](./images/radial-displacement.png)

Note that, overall, all solutions agree very well with the reference solution.
Only the RNF solution that varies according to the time-step used (dt in the
legend). This behavior is expected due to the first-order approach used for the
temporal discretization, as detailed in the following paper:

Ž. Tuković, M. Bukač, P. Cardiff, H. Jasak, and A. Ivanković, “Added Mass Partitioned Fluid–Structure Interaction Solver Based on a Robin Boundary Condition for Pressure,” OpenFOAM, pp. 1–22, 2019, doi: 10.1007/978-3-319-60846-4_1.

In the figure below, we can see the performance comparison among the different strategies in terms of the number of FSI coupling iterations.

![](./images/coupling-iterations.png)

---

## Tips for fluid-solid analyses

All the tips provided in the beamInCrossFlow tutorial also applies for this
case. But note from this case that different coupling algorithms may have a
much better performance, although they may not always help or you will have to
tweak with the setup to find a converging/stable solution.
