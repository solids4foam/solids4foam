---
sort: 14
---

# Wobbly Newton: `wobblyNewton`

---

Prepared by Philip Cardiff and Ivan Batistić

---

## Tutorial Aims

- Demonstrate a transient small-strain linear elastic analysis on a complex
  geometry with an unstructured mesh obtained using `cfMesh`;
- Demonstrate the use of `surfaceToPatch` for extracting patches corresponding
  to the provided `.stl` files.

---

## Case Overview

This case considers a flexible Newton-shaped body attached to a fixed base.
The geometry is supplied by the `newtonDecimated.stl` surface file.
The volume mesh is generated using the `cfMesh` `cartesianMesh`
utility, with the meshing settings prescribed in `system/meshDict`.

The body is modelled as linear elastic rubber using the small-strain
assumption, with density $$\rho = 1500$$ kg/m$$^3$$, Young's modulus
$$E = 1$$ MPa, and Poisson's ratio $$\nu = 0.3$$. The flat base of the statue
has a fixed displacement boundary condition, and the remaining surfaces are
traction-free. The case is solved as a transient problem from $$t = 0$$ s to
$$t = 0.36$$ s using a time-step of $$\Delta t = 0.01$$ s. The body force is
set through the `g` field as $$(0, 9.81, 0)$$ m/s$$^2$$.

![Figure 1: Input geometry and computational
mesh](images/wobblyNewton-geometry-mesh.png)**Figure 1: Input geometry (left) and computational mesh (right)**

```note
The mesh density does not capture all geometry features, which is acceptable
for this demonstration. For real cases, a finer decimated surface and a finer
volume mesh should be used.
```

---

## Expected Results

This case does not have an analytical solution for direct quantitative
comparison. The expected response is qualitative: the fixed base constrains the
lower part of the geometry while the flexible Newton-shaped body undergoes a
transient oscillating displacement under the applied gravitational (body) force.

![Figure 2: Deformed geometry coloured by displacement
magnitude](images/wobblyNewton-results.png)**Figure 2: Deformed geometry coloured by displacement magnitude; deformation is
scaled**

The `solidPointDisplacement` function object records the displacement of a
point close to $$(0.0408, 0.0048, 0.0929)$$. The output is written to
`postProcessing/0/solidPointDisplacement_pointDisp.dat`. The monitored point
displacement can be used to check the transient response.

---

## Running the Case

The tutorial case is located at
`solids4foam/tutorials/solids/linearElasticity/wobblyNewton`.

The case is solved with the `PETScSNES` solution algorithm. By default, the
PETSc options are read from `petscOptions.mf.hypre`, which uses a matrix-free
Newton method with a HYPRE preconditioner.

Initially, this case was a demonstration of vertex-centred discretisation;
therefore, `system/fvSchemes` also contains time discretisation schemes used by
the vertex-centred solid model.

This case requires OpenFOAM.com, PETSc, and the `cartesianMesh` utility. It is
skipped by the `Allrun` script if `PETSC_DIR` is not set or if `cartesianMesh`
is not available.

The case can be run in serial using:

```bash
./Allrun
```

The script converts the case format if required, creates the mesh using
`cartesianMesh`, creates the `base` patch using `surfaceToPatch`, and then runs
the `solids4Foam` solver.

The case can also be run in parallel using:

```bash
./Allrun parallel
```

In parallel, the case is decomposed using `decomposePar`, solved in parallel,
and reconstructed using `reconstructPar`.
