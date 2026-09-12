---
sort: 35
---

# Method of Manufactured Solutions: `manufacturedSolution`

Prepared by Philip Cardiff and Federico Mazzanti

## Tutorial aims

- Demonstrate verification with the method of manufactured solutions (MMS).
- Measure displacement and stress errors against a smooth three-dimensional
  analytical field.
- Provide a fast default case and opt-in mesh-convergence studies.

## Case overview

The domain is a 0.2 m cube with Young's modulus 200 GPa and Poisson's ratio
0.3. The manufactured displacement is

$$
\boldsymbol{u} =
\begin{bmatrix}a_x \\ a_y \\ a_z\end{bmatrix}
\sin(4\pi x)\sin(2\pi y)\sin(\pi z),
$$

where $a_x=2\times10^{-6}$ m, $a_y=4\times10^{-6}$ m, and
$a_z=6\times10^{-6}$ m. The tutorial-local library supplies the corresponding
body force, displacement boundary condition, analytical stress, and error
function object.

The default case uses the segregated solution procedure and a regular
$5\times5\times5$ hexahedral mesh, so it runs quickly:

```bash
./Allrun
```

The run script also accepts a solution approach and mesh type:

```bash
./Allrun segregated hex
./Allrun petscSnes hex
./Allrun segregated distHex
./Allrun petscSnes tet
./Allrun petscSnes poly
```

The `petscSnes` approach requires a PETSc-enabled solids4foam build. The `tet`
and `poly` meshes require Gmsh. This tutorial currently supports OpenFOAM.com.

## Verification and regression

The default case has a regression test:

```bash
./regressionTest.sh
```

The opt-in [`verification/`](verification/) directory migrates the mesh and
solver variants from `solid-benchmarks/linearElasticity/manufacturedSolution`.
It is separate from the normal tutorial regression suite because it runs many
cases. See the verification README for variants, commands, and acceptance
criteria.
