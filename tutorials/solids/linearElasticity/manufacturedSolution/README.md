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
./Allrun highOrder-movingLeastSquares hex
./Allrun highOrder-kExactLeastSquares hex
```

The `petscSnes` and high-order approaches require a PETSc-enabled solids4foam
build. The `tet` and `poly` meshes require Gmsh. The tutorial supports OpenFOAM.com,
OpenFOAM.org, and foam-extend.

The `tet` and `poly` runs use `gmsh/tet-structured.geo` by default. Set
`GMSH_MESH=tet-unstructured` in `Allrun` to use the unstructured tetrahedral
alternative. A standalone `gmsh/hex-structured.geo` is also provided; the
standard `hex` run uses `blockMesh`. All Gmsh scripts read `gmsh/meshSpacing.geo`.

The high-order approaches use cubic displacement reconstruction, face
quadrature, and volume integration of the manufactured body force. The
integrated body force is divided by cell volume before insertion as an
`fvOptions` source density on OpenFOAM.com. On OpenFOAM.org and foam-extend,
`Allrun` selects the tutorial-local `manufacturedSolutionSolid` model, which
adds the same source through the linear solid model's `fvOptionsSource()` function.
The source coefficients remain in `constant/fvOptions` on all versions. The
boundary condition evaluates the analytical displacement at face quadrature
points.

Rebuild solids4foam after updating the core source interface. `Allrun`
builds the tutorial library, but does not rebuild the core libraries or solver.

`movingLeastSquares` stores point values, so displacement errors use the
analytical solution at cell centres. `kExactLeastSquares` stores cell averages,
so its displacement errors use the volume-averaged analytical solution.
Analytical stress is evaluated at cell centres for both approaches.

## Verification and regression

The regression test checks `segregated`, `petscSnes`, and both high-order
approaches on coarse hex and tet meshes. Both high-order reconstructions share
the same displacement and stress L2 tolerances:

```bash
./regressionTest.sh
```

PETSc regression runs are skipped when `PETSC_DIR` is unset; tet runs are
skipped when Gmsh is unavailable. Logs are retained under `regressionTests/`;
use `./regressionTest.sh --check-only` to check existing results.

The opt-in [`verification/`](verification/) directory migrates the mesh and
solver variants from `solid-benchmarks/linearElasticity/manufacturedSolution`.
It is separate from the normal tutorial regression suite because it runs many
cases. See the verification README for variants, commands, and acceptance
criteria.
