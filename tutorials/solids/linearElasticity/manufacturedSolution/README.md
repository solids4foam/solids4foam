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
`GMSH_MESH` when running `Allrun` to use another Gmsh script, for example
`GMSH_MESH=tet-unstructured ./Allrun petscSnes tet` for the unstructured
tetrahedral alternative. `gmsh/hex-structured.geo` is provided for manual use
in the same way (`GMSH_MESH=hex-structured`); no script selects it, and the
standard `hex` run uses `blockMesh`. All Gmsh scripts read
`gmsh/meshSpacing.geo`.

The high-order approaches use cubic displacement reconstruction, face
quadrature, and volume integration of the manufactured body force. The
integrated body force is divided by cell volume before insertion as an
`fvOptions` source density on OpenFOAM.com. On OpenFOAM.org and foam-extend,
`Allrun` selects the tutorial-local `manufacturedSolutionSolid` model, which
adds the same source through the linear solid model's `fvOptionsSource()` function.
The source coefficients remain in `constant/fvOptions` on all versions. The
boundary condition evaluates the analytical displacement at face quadrature
points.

`movingLeastSquares` stores point values, so displacement errors use the
analytical solution at cell centres. `kExactLeastSquares` stores cell averages,
so its displacement errors use the volume-averaged analytical solution.
Analytical stress is evaluated at cell centres for both approaches.

## Mixed displacement-pressure approach

`petscSnesPressure` solves the same problem with the mixed
displacement-pressure formulation (`solvePressure yes`) and the PETSc SNES
solution procedure:

```bash
./Allrun petscSnesPressure
```

The pressure is compression positive, $\boldsymbol\sigma =
\tilde{\boldsymbol\sigma} - p\mathbf{I}$, and the pressure equation is
$-p/K - \nabla\cdot\boldsymbol{u} = 0$ plus stabilisation. For this linear
elastic solution the analytical pressure is

$$
p = -\tfrac{1}{3}\operatorname{tr}\boldsymbol\sigma
  = -\left(\lambda + \tfrac{2}{3}\mu\right)\nabla\cdot\boldsymbol{u}
  = -K\,\nabla\cdot\boldsymbol{u},
$$

evaluated at the cell centres. The case uses one stabilisation configuration:
momentum `diffStencilLaplacian` with scale factor 0.1, and pressure
`generalisedEvenOrderLaplacian` with `laplacianPower` 1 (the bi-Laplacian),
scale factor 0.125 and Jacobian scale factor 1.0, unnormalised. The PETSc
options use a Schur-complement field split with BoomerAMG on both blocks.
These settings verify the current implementation; they are not intended to
reproduce published results.

The error function object reports the pressure error whenever the solid model
solves for the pressure (`cellPressure` may be set to `yes` or `no` in
`system/controlDict` to override this; `yes` without `solvePressure` is an
error). It writes `analyticalP` and `pDifference` and logs

```text
Writing pDifference field
    Pressure error norms: mean L1, cell RMS, LInf, volume-weighted RMS:
    Magnitude: 183305 230706 908771 230706
```

As for the displacement and stress (whose second column is labelled
"mean L2"), the cell RMS is the root mean square over cells, unweighted. The
volume-weighted RMS, $\sqrt{\sum_P V_P e_P^2 / \sum_P V_P}$, is the
continuous $L^2$ norm of the cell-wise constant error divided by the square
root of the domain volume; on this uniform mesh the two coincide. The values above are those of the default
$5^3$ mesh. Displacement-only approaches are unaffected.

The case scales the mixed PETSc system (`scaleMixedPETScFields yes`): the
pressure unknown is stored as $p/2\mu$ and the pressure rows are multiplied
by $2\mu$, the volume average of the shear stiffness the law declares. The
converged solution is unchanged. Without it, the first Krylov solve on the
$20^3$ and finer meshes stopped at `ksp_max_it`, because the matrix-free
differencing step, scaled by the unknowns, is too small for pressures of
order $10^6$ Pa. With the scaling no linear solve fails, so the case keeps
PETSc's default limit of one failed linear solve per Newton solve. The tests
require nonlinear convergence and report any failed linear solves with their
PETSc reasons (`DIVERGED_ITS` is the iteration limit).

## Verification and regression

The regression test checks `segregated`, `petscSnes`, and both high-order
approaches on coarse hex and tet meshes, and `petscSnesPressure` on the coarse
hex mesh. Both high-order reconstructions share the same displacement and
stress L2 tolerances. The mixed check also requires nonlinear convergence and a
pressure error below its tolerance; a missing pressure error fails the test.
Approach names and `--mesh` select a subset:

```bash
./regressionTest.sh
./regressionTest.sh --mesh hex petscSnesPressure
```

The opt-in [`verification/`](verification/) directory holds mesh and solver
variants of this case. It is separate from the normal tutorial regression suite
because it runs many cases. See the verification README for variants, commands,
and acceptance criteria.
