---
sort: 4
---

# Annular plate with a curved free edge: `annularPlate`

---

Prepared by Ivan Batistić

---

## Tutorial Aims

- Demonstrate the free-edge boundary conditions of the finite area
  Kirchhoff-Love plate model on a curved edge, where the edge curvature enters
  the Kirchhoff shear force condition;
- Verify the predictions against closed-form solutions of the annular plate
  and measure the order of convergence.

---

## Case Overview

An annular plate with inner radius $$a = 1$$ m, outer radius $$b = 2$$ m and
thickness $$h = 0.1$$ m is loaded by a transverse pressure of amplitude $$p_0 =
1000$$ Pa, with Young's modulus $$E = 200$$ GPa and Poisson's ratio $$\nu =
0.3$$. The density is set to zero so that the single time step gives the static
solution. Three configurations are provided and selected by an argument of the
`Allrun` script:

| Configuration | Clamped edge | Free edge | Load |
| --- | :---: | :---: | --- |
| `clampedInner` (default) | $$r = a$$ | $$r = b$$ | uniform $$p_0$$ |
| `clampedInnerCosineLoad` | $$r = a$$ | $$r = b$$ | $$p_0 \cos\varphi$$ |
| `clampedOuter` | $$r = b$$ | $$r = a$$ | uniform $$p_0$$ |

The mesh has one block per quadrant with the circumferential edges as arcs;
the tutorial mesh has 16 radial and 152 circumferential cells. The fields `M`
and `w` in `0/finite-area` are symbolic links to the `*.clampedInner` or
`*.clampedOuter` versions, which `Allrun` sets. The cosine load varies over the
finite area faces, and since there is no field-expression utility for finite
area fields, `Allrun` writes it with the small Python script
`writeCosineLoad.py` after the finite area mesh is created.

On a curved free edge with curvature $$\kappa$$ (positive on the convex outer
edge, negative on the concave inner edge) the Kirchhoff shear force condition
contains a curvature term, so the `freeEdgeMoment` condition imposes

$$
\frac{\partial M}{\partial n} = D (1 - \nu) \frac{\partial}{\partial s}
\left( \frac{\partial}{\partial s} \frac{\partial w}{\partial n}
- \kappa \frac{\partial w}{\partial s} \right),
$$

where the curvature is evaluated from the normals of the neighbouring edges.
For the axisymmetric configurations both terms vanish on the circular edge; the
cosine load exercises them, and the curvature term is about 40% of the
straight-edge term there. The `freeEdgeDisplacement` condition includes the
non-orthogonal correction of the normal slope, which matters on a curved edge.
The clamped edge uses `fixedValue` and `clampedMoment` as in the
`cantileverPlate` tutorial. On this geometry the explicit coupling of the
conditions is driven entirely by the clamped edge; a curved free edge alone
(for example with a simply supported inner edge) does not converge with the
present segregated scheme.

---

## Expected Results

The axisymmetric configurations have the classical solution [1]

$$
w(r) = C_1 + C_2 r^2 + C_3 \ln r + C_4 r^2 \ln r + \frac{p_0 r^4}{64 D},
$$

where $$C_4 = -p_0 r_f^2 / (8 D)$$ follows from the vanishing shear force at
the free edge $$r_f$$, and $$C_1$$ to $$C_3$$ from $$w = \partial w /
\partial r = 0$$ at the clamped edge and the vanishing radial moment $$w'' +
\nu w' / r = 0$$ at the free edge. The cosine load has the solution $$w =
f(r) \cos\varphi$$ with

$$
f(r) = A r + B r^3 + \frac{C}{r} + E\, r \ln r + \frac{p_0 r^4}{45 D},
$$

where the four constants follow from $$f = f' = 0$$ at the clamped edge and the
vanishing radial moment and Kirchhoff shear force at the free edge. The
verification script evaluates both solutions.

The deflection of the free edge is $$7.503 \times 10^{-6}$$ m for
`clampedInner`, $$7.214 \times 10^{-6}$$ m (amplitude of the $$\cos\varphi$$
mode) for `clampedInnerCosineLoad` and $$4.603 \times 10^{-6}$$ m for
`clampedOuter`. The tutorial mesh overpredicts them by 8.0%, 8.5% and 8.6%,
and the overprediction halves with every mesh refinement (see the verification
study below). The moment sum $$M$$ is predicted to within 1% in the whole plate
for the inner-clamped configurations, and for the cosine load the normal
gradient of $$M$$ on the free edge, which the curvature term determines, is
within 4.4% of the exact value on the tutorial mesh.

---

## Running the Case

The tutorial case is located at
`solids4foam/tutorials/solids/beamsPlatesShells/annularPlate`. The case can be
run using the included `Allrun` script, i.e. `> ./Allrun` for the default
configuration or `> ./Allrun clampedInnerCosineLoad` and `> ./Allrun
clampedOuter` for the other two. The script links the boundary condition
files, creates the mesh using `blockMesh` (`> blockMesh`), creates the finite
area mesh using `makeFaMesh` (`> makeFaMesh`), writes the cosine load when
requested, and runs the case with the `solids4foam` solver (`> solids4Foam`).
The case runs with OpenFOAM.com v2412 or newer and takes a few seconds on the
tutorial mesh.

---

## Verification and Convergence Study

The opt-in [`verification/`](verification/) directory adds the mesh-convergence
sweep for all three configurations. It scales the radial and circumferential
divisions by powers of two from $$8 \times 76$$ to $$64 \times 608$$ cells,
measures the error of the free-edge deflection, the relative $$L_2$$ errors of
$$w$$ and $$M$$ over the plate and, for the cosine load, the error of the
normal gradient of $$M$$ on the free edge against the closed-form solutions,
and reports the observed orders of convergence:

```bash
cd verification
./Allverify
```

The free-edge deflection converges at first order in all three configurations
(relative error against the exact solution):

| Mesh | Inner clamped | Inner clamped, cosine | Outer clamped |
|:---:|---:|---:|---:|
| 8 x 76 | 13.7% | 14.4% | 15.5% |
| 16 x 152 | 8.0% | 8.5% | 8.6% |
| 32 x 304 | 4.3% | 4.6% | 4.5% |
| 64 x 608 | 2.2% | 2.4% | 2.3% |
| Net order | 0.88 | 0.87 | 0.92 |

The relative $$L_2$$ error of $$w$$ follows the same first-order trend (16.4% to
2.6% for `clampedInner`), and the error of the free-edge normal gradient of
$$M$$ for the cosine load falls from 7.2% to 1.2% (net order 0.84). The whole
sweep takes about 25 minutes in serial on an Apple M1 Ultra, almost all of it
on the finest meshes.

The study is separate from `regressionTest.sh` and is not run by the normal
tutorial tests. See the verification README for the mesh levels, options and
acceptance criteria.

![Figure 1: Relative errors against the radial mesh
spacing](./images/annularPlate-convergence.png)

Figure 1: Relative errors of the free-edge deflection and of $$w$$ against the
radial mesh spacing for the three configurations.

---

### References

[1] [Timoshenko, S., & Woinowsky-Krieger, S., Theory of plates and shells.
1959.](https://www.cap-recifal.com/ccs_files/articles/cuveaqua1_denisio/Timoshenko_-_Theory_of_plates_and_shells.pdf)
