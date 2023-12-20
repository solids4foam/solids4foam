# The cooling of a solid slab: `slabCooling`

---

Prepared by Philip Cardiff, Duo Huang and Ivan Batistić

---

## Tutorial Aims

- Demonstrate the solution of a thermal stress problem;
- Check the accuracy of a solid model against the expected solution.

---

## Case Overview

This case demonstrates the process of cooling the slab from the initial $$T_o$$ to the reference temperature, $$T_{ref}$$. The slab is $$6$$ m long, $$2$$ m wide, and $$1$$ m deep (Figure 1). The material of the slab is linear elastic with a Young modulus of $$E = 69$$ GPa and a Poisson's ratio of $$\nu = 0.33$$. Material density is $$\rho=7854$$ kg/m$$^3$$, thermal conductivity is $$k=250$$ W/mK, specific heat capacity is $$C = 434$$ J/kgK and coefficient of linear thermal expansion is $$\alpha = 2.3\cdot 10^{-05}$$ 1/K. The initial slab temperature is set to $$800$$ K, and it is cooling to the reference temperature of $$T_{ref} = 300$$ K. Gravitation is neglected, and the case is solved under the plane strain assumption. All slab faces are traction-free. The slab is discretised with $$30 \times 10 \times 5$$ CVs.

A steady-state solution is required, meaning the case can be solved without temporal discretisation (using one loading increment). Note that in this case, the value of density, thermal conductivity and specific heat capacity are not important for the resulting stress and displacement fields.

<div style="text-align: center;">
  <img src="./images/slabCooling-geometry.png" alt="Image" width="450">
    <figcaption>
     <strong>Figure 1: Problem geometry; computational mesh</strong>
    </figcaption>
</div>

---

## Expected Results

- All slab faces are traction-free, meaning that the slab can expand/contract freely in all directions. Accordingly, the reached solution should be stress-free. This is shown in Figure 3 where one can see that the resulting stress field is in magnitude of a few hundred pascals which is negligible and is caused by the usage of iterative procedure and chosen tolerances;
- Herein, the case considers slab cooling, meaning that the slab will contract towards its centre, as shown in Figure 2;

<div style="text-align: center;">
  <img src="./images/slabCooling-D.png" alt="Image" width="400">
    <figcaption>
     <strong>Figure 2: Displacement magnitude field</strong>
    </figcaption>
</div>

<div style="text-align: center;">
  <img src="./images/slabCooling-sigmaEq.png" alt="Image" width="400">
    <figcaption>
     <strong>Figure 3: Equivalent stress field</strong>
    </figcaption>
</div>

Since all slab faces have specified the traction-free boundary condition,  a displacement value is not specified at any boundary point. This means that the body is unconstrained in all three directions and can freely float in space. To avoid this, a displacement value should be specified at some location. This can be done in the following ways:
- set the displacement of one internal cell; 
- set the displacement of one boundary face; 
- if possible, take advantage of symmetry planes and use them when modelling. 

The first two points restrict rigid translations; however, they can not restrict rigid rotation, what to be careful about. Symmetry plane remove rigid translation in the direction normal to it.

In this case, zero displacement of internal cell is specified in `solidProperties` dict using following lines:

```c++
// Set displacements of internal cells
cellDisplacements
{
    // Set displacement to zero for the cell closest to the centre of the block
    cellDisp1
    {
        approximateCoordinate    (3 1 0.5);
    	displacement             (0 0 0);
	}
}
```

---

## Running the Case

The tutorial case is located at `solids4foam/tutorials/solids/thermoelasticity/slabCooling`. The case can be run using the included `Allrun` script, i.e. `> ./Allrun`.  The `Allrun` creates the mesh using `blockMesh` (`> blockMesh`) after which the case is run with the `solids4Foam` solver (`> solids4Foam`). 


---

### References

[1] [Link to forum discussion where initial case is provided](https://www.cfd-online.com/Forums/openfoam-community-contributions/126706-support-thread-solid-mechanics-solvers-added-openfoam-extend-22.html#post726239)
