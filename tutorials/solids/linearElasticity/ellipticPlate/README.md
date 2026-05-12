# Out-of-plane bending of an elliptic plate: `ellipticPlate`

---

Prepared by Philip Cardiff and Ivan Batistić

---

## Tutorial Aims

- Demonstrate how to perform a 3D linear-static stress analysis in solids4foam
- Compare the performance of three variants of the same tutorial case:
  `petscSnes`, which uses the `linearGeometryTotalDisplacement` solid model
  with a PETSc SNES solution algorithm; `segregated`, which uses the same
  solid model with an implicit segregated solution algorithm; and
  `unsCoupled`, which uses the original
  `coupledUnsLinearGeometryLinearElastic` solid model.
- Demonstrates using the solids4foam solid solver based on the PETSc SNES
  nonlinear solver

---

## Case Overview

This 3-D test case consists of a thick elliptic plate with a centred elliptic hole (Fig. 1); a constant pressure of $$1$$ MPa is applied to the upper surface, and the outer surface is fully clamped. The case has been described by the National Agency for Finite Element Methods and Standards (NAFEMS) [1], and in the context of the FV method has been examined and benchmarked in detail by Demirdžić et al. [2]. The assumed mechanical properties are Young’s modulus = $$210$$ GPa and Poisson’s ratio $$\nu=0.3$$. Due to a double symmetry, only a quarter of the plate is analysed. The thickness of the plate is $$0.6$$ m, and the inner and outer ellipses are given by:
$$
\left(\dfrac{x}{2}\right)^2 + \left(\dfrac{y}{1}\right)^2 = 1 \qquad \text{inner elipse}
$$

$$
\left(\dfrac{x}{3.25}\right)^2 + \left(\dfrac{y}{2.75}\right)^2 = 1 \qquad \text{outer elipse}
$$

The problem is solved as quasi-static with one-loading increment and neglecting gravitational effects. Mesh is generated using blockMesh uditlity and consists of 12 layers of cell in plane depth,  24 in circumferential direction and 16 in radial direction, figure 2. Mesh density corresponds to mesh 3 used in [3]. 

```Note
The target stress value given in the NAFEMS benchmark cannot be directly used for comparison as the fixed displacement boundary conditions used in the current study mimic those used in Demirdžić et al. [2] and do not correspond exactly to those in the NAFEMS benchmark.
```

![](images/ellipticPlate-geometry.png)



**Figure 1 - Problem geometry and loading [2]**



![Geometry](images/ellipticPlate-mesh.png)

**Figure 2 - Hexahedral mesh (mesh 3 in [2])**

## Expected Results

Figure 3 shows the stress contours one the $$z=0.3$$ m $$x$$–$$y$$ plane for the finest mesh.
The predictions agree closely with the results of Demirdzic et al. [2]. 
Further comparisons can be found in [3]:
[Cardiff et al., 2016, A block-coupled Finite Volume methodology for linear
elasticity and unstructured meshes, Comp. and
Struct.](https://www.sciencedirect.com/science/article/pii/S0045794916306046)

![](images/ellipticPlate-stress.jpeg)

**Figure 3 - Stress component distributions on the plane $$z=0.3$$ m. Second image is from [2] and third is generated using finite elements with Abaqus [3].**

The wall-clock times and memory requirements for each run are given in Table 1,
where the results for the `unsCoupled` and `segregated` variants are compared.
The `unsCoupled` solver used a bi-conjugate gradient stabilised linear solver
with ILU(0) preconditioner, while the `segregated` variant used an implicit
segregated solve with the same convergence tolerance of
$$1 \times 10^{-6}$$. In this case, the `unsCoupled` variant is approximately
four times faster than the `segregated` variant but requires about four times
more memory. From Table 1, it can be seen that the `unsCoupled` is faster than the `segregated `
method by a factor of 2.5–6 times; as expected, the memory requirements are greater, by approximately 4.5 times in the largest mesh case. Interestingly, when compared to the FE solution, the `unsCoupled` is faster and considerably more memory efficient in all cases; in the finest mesh case the `unsCoupled` is almost 6 times faster and uses 8 times less memory. 

```note
The wall-clock times given in Table 1 were recorded in 2015 using one core of a
2.4 GHz Intel Ivy Bridge CPU. Better performance can be expected using a new
machine.
```

**Table 1: Wall-clock time (in s) and maximum memory usage (in MB) [2].** 

| Mesh    | unsCoupled |            | segregated |            | Abaqus   |            |
| ------- | ---------- | ---------- | ---------- | ---------- | -------- | ---------- |
|         | **Time**   | **Memory** | **Time**   | **Memory** | **Time** | **Memory** |
| 500     | 0.08       | 11         | 58         | 20         | 3        | 40         |
| 4 500   | 0.5        | 50         | 384        | 27         | 4        | 113        |
| 12 500  | 1.4        | 140        | 1 387      | 43         | 5        | 197        |
| 50 000  | 6          | 570        | 4 737      | 112        | 16       | 881        |
| 200 000 | 36         | 2 500      | -          | -          | 73       | 1800       |

```warning
The `coupledUnsLinearGeometryLinearElastic` solid model currently does not run
in parallel. For a coupled solid model that *does* run in parallel, use the
`vertexCentredLinGeomSolid` solid model.
```

## Running the Case

The tutorial case is located at
`solids4foam/tutorials/solids/linearElasticity/ellipticPlate`. The case can be
run using the included `Allrun` script. The `Allrun` script optionally takes an
argument which specifies the solution approach:

```bash
./Allrun                # Defaults to the petscSnes variant
./Allrun segregated     # Segregated variant
./Allrun unsCoupled     # Original foam-extend-only variant
```

In all cases, the `Allrun` script simply creates the mesh using `blockMesh`
and then runs the solids4foam solver.

The default `petscSnes` variant uses the `linearGeometryTotalDisplacement`
solid model with PETSc SNES. The original `unsCoupled` variant is also
provided through the `.unsCoupled` file suffixes, and uses the
`coupledUnsLinearGeometryLinearElastic` solid model together with the original
`block*` boundary conditions and coupled discretisation settings. The
`unsCoupled` variant can currently only be run using solids4foam built on
foam-extend.

The `segregated` variant uses the same `linearGeometryTotalDisplacement`
solid model as `petscSnes`, but with `solutionAlgorithm implicitSegregated;`
in `constant/solidProperties`. This variant is intended for comparison with the
PETSc SNES path using the same geometry and boundary conditions.

---

## References

[1] National Agency for Finite Element Methods and Standards (U.K.). The standard NAFEMS benchmarks. NAFEMS; 1990.

[2] [Demirdžić, I., Muzaferija, S., Perić, M.: Benchmark solutions of some structural analysis problems using finite-volume method and multigrid acceleration. International journal for numerical methods in engineering 40(10), 1893–1908 (1997)](https://onlinelibrary.wiley.com/doi/10.1002/(SICI)1097-0207(19970530)40:10%3C1893::AID-NME146%3E3.0.CO;2-L)

[3] [P. Cardiff, Ž. Tuković, H. Jasak, A. Ivanković, A block-coupled finite volume
 methodology for linear elasticity and unstructured meshes. *Computers and Structures*,
 175, 2016, 100–122, 10.1016/j.compstruc.2016.07.004.](https://doi.org/10.1016/j.compstruc.2016.07.004)

[4] [P. Cardiff, D. Armfield, Ž. Tuković, I. Batistić, A Jacobian-free
Newton-Krylov method for cell-centred finite volume solid mechanics.
_International Journal for Numerical Methods in Engineering_, 127, e70268,
2026, 10.1002/nme.70268.](https://doi.org/10.1002/nme.70268)