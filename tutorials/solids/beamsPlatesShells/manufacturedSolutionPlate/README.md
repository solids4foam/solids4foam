---
sort: 5
---

# Manufactured solution for the Kirchhoff plate: `manufacturedSolutionPlate`

---

Prepared by Ivan Batistić

---

## Tutorial Aims

- Verify the finite area Kirchhoff-Love plate model, including its clamped,
  free and symmetry edge conditions, with the method of manufactured
  solutions, where the exact deflection, moment sum and rotation are known
  everywhere;
- Measure the order of accuracy of the deflection, the moment sum and the
  rotation in the $$L_2$$ and $$L_\infty$$ norms, and along the edges;
- Show how the discretisation behaves on a uniform quadrilateral mesh, on a
  randomly perturbed quadrilateral mesh and on an unstructured triangular
  mesh.

---

## Case Overview

The plate is the $$L \times W \times h = 10$$ m $$\times~10$$ m $$\times~0.1$$
m plate of the other plate tutorials, with Young's modulus $$E = 200$$ GPa and
Poisson's ratio $$\nu = 0.3$$. The `left` edge ($$x = 0$$) is clamped
(`fixedValue` for `w`, `clampedMoment` for `M`), the `right` edge ($$x = L$$)
is free (`freeEdgeDisplacement`, `freeEdgeMoment`) and the `top` and `bottom`
edges are symmetry lines (`zeroGradient` for both fields). The density is
zero, so the single time step gives the static solution.

The manufactured deflection is

$$
w = w_0 \left[ a \left(\frac{x}{L}\right)^2 + b \left(\frac{x}{L}\right)^3
+ \left(\frac{x}{L}\right)^4 \right] \cos\left(\frac{\pi y}{W}\right),
$$

with $$w_0 = 0.1039$$ m, $$a = 1.8657$$ and $$b = -2.3845$$. The quartic
vanishes with its slope at the clamped edge, and $$a$$ and $$b$$ are chosen
so that the bending moment $$M_n = -D (\partial^2 w/\partial n^2 + \nu\,
\partial^2 w/\partial s^2)$$ and the Kirchhoff shear force $$V_n$$ vanish on
the free edge; the cosine makes the normal derivatives of $$w$$ and $$M$$
vanish on the top and bottom edges. The largest deflection is $$0.05$$ m on
the free edge. The load that produces this deflection,

$$
p = D \nabla^4 w,
$$

is written into `0/finite-area/p` by `writeManufacturedSolution.py` from the
face centres of the finite area mesh, together with the exact deflection and
moment sum as the fields `wExact` and `MExact` for comparison in ParaView.
The script holds the plate data and the coefficients, and the verification
study imports it, so the manufactured solution is defined in one place.

Three meshes with the same nominal spacing are available:

- `hex`: a uniform $$50 \times 50$$ quadrilateral mesh from `blockMesh`
  (default);
- `distHex`: the same mesh with every interior point moved by a random
  fraction of up to 15% of the local spacing, by `perturbPlateMesh.py`,
  which moves the two surfaces of the one-cell-thick mesh together so that the
  cells stay straight prisms (the `perturbMeshPoints` utility perturbs the two
  surfaces independently and its mesh checks then fail on the twisted side
  faces);
- `tri`: an unstructured triangular mesh from Gmsh with the same target
  spacing, extruded to one layer of prisms (`gmsh/tri-unstructured.geo`).

---

## Running the Case

The tutorial case is located at
`solids4foam/tutorials/solids/beamsPlatesShells/manufacturedSolutionPlate`.
It is run with the included `Allrun` script, which takes the mesh type as an
optional argument:

```bash
./Allrun            # uniform quadrilaterals (hex)
./Allrun distHex    # randomly perturbed quadrilaterals
./Allrun tri        # unstructured triangles, requires gmsh
```

The script creates the mesh, creates the finite area mesh with `makeFaMesh`,
writes the manufactured load with `writeManufacturedSolution.py` (python3)
and runs `solids4Foam`. The case runs with OpenFOAM.com v2412 or newer and
takes about 10 s on the quadrilateral meshes and 90 s on the triangles with
the segregated algorithm.

Setting `coupled yes;` in `constant/solidProperties` assembles the `M` and
`w` equations with their edge conditions into one linear system that is
solved with PETSc (a direct solve with MUMPS by default). The non-orthogonal
correction of the Laplacians is part of the matrix, so one solve gives the
solution on any mesh, in a fraction of a second; this requires solids4foam to
be compiled with PETSc.

---

## Expected Results

The errors are measured against the exact solution at the cell centres. On
the $$50 \times 50$$ `hex` mesh the relative $$L_2$$ error of the deflection
is 0.57% and that of the moment sum 0.59%; the largest deflection is
$$0.04984$$ m at the cells next to the free edge, where the exact value at
those cell centres is $$0.04873$$ m. The perturbed mesh gives errors of 2.9%
and 1.4%, and the triangular mesh of 44% and 52%, see the verification study
below for the reason. Figure 1 shows the convergence of the errors on the
three meshes.

![Figure 1: Convergence of the relative L2 errors of w and M on the three
meshes](./images/manufacturedSolutionPlate-convergence.png)

Figure 1: Relative $$L_2$$ errors of the deflection (solid) and the moment sum
(dashed) against the mesh spacing for the three mesh types.

---

## Verification and Convergence Study

The opt-in [`verification/`](verification/) directory refines each mesh type
from 25 to 200 cells per side and reports the errors of `w`, `M` and the
rotation `theta` in the $$L_2$$ and $$L_\infty$$ norms, the errors of the
boundary values on the free and clamped edges, and the observed orders:

```bash
cd verification
./Allverify              # segregated solver
./Allverify --coupled    # block-coupled solver, seconds per sweep
```

Both solvers give the same discrete solution to better than $$10^{-5}$$
relative on every mesh. On the uniform mesh the errors converge at first
order, which is the order of the one-sided estimates of the free edge
conditions:

| Cells | $$L_2$$ error $$w$$ | $$L_2$$ error $$M$$ | $$L_2$$ error $$\theta$$ |
|:---:|---:|---:|---:|
| 25 x 25 | 1.22% | 1.07% | 1.68% |
| 50 x 50 | 0.57% | 0.59% | 0.86% |
| 100 x 100 | 0.28% | 0.31% | 0.44% |
| 200 x 200 | 0.14% | 0.16% | 0.22% |

with net orders of 1.06, 0.91 and 0.98. On the perturbed and triangular
meshes the errors do not converge: they stay between 0.5% and 3% on the
perturbed quadrilaterals and at about 40% on the triangles. The study
isolates the cause by replacing each edge condition in turn by the exact
boundary values: the interior discretisation, the clamped condition, the
symmetry condition and the `freeEdgeMoment` condition all keep their
accuracy on the skewed meshes, whereas the `freeEdgeDisplacement` condition
does not. On a cell whose centre is offset tangentially from its boundary
edge, the extrapolation of the edge displacement from the cell gradient
excites a smooth spurious bending mode of the plate, equivalent to a shear
force along the free edge, whose amplitude is proportional to the skewness
and independent of the mesh spacing. On the uniform mesh the truncation
errors are symmetric and this mode is never excited. The free edge
conditions are therefore only verified on meshes whose boundary cells are
orthogonal to the plate edges, and a form of the displacement condition that
imposes the vanishing bending moment through the edge value of `M` is the
subject of further work.

---

## References

P. J. Roache, Code verification by the method of manufactured solutions,
Journal of Fluids Engineering 124 (2002) 4-10.

M. Torlak, A finite-volume method for coupled numerical analysis of
incompressible fluid flow and linear deformation of elastic structures, PhD
thesis, TU Hamburg-Harburg, 2006 (the rotation-free plate formulation).
