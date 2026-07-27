# RBF Mesh Motion Solver

`RBFMeshMotionSolver` moves internal mesh points by radial basis function
interpolation from selected boundary control points. It is intended for cases
where the moving boundary displacement is known, for example the fluid mesh in
a fluid-solid interaction case.

## Example: `perpendicularFlap`

The `perpendicularFlap` tutorial has a small two-dimensional fluid mesh and is a
convenient starting point for testing the solver. The default tutorial keeps the
finite-volume `velocityLaplacian` mesh motion solver enabled, but the following
settings can be used in
`tutorials/fluidSolidInteraction/perpendicularFlap/constant/fluid/dynamicMeshDict`
with OpenFOAM.com:

```cpp
dynamicFvMesh dynamicMotionSolverFvMesh;

"solver|motionSolver" RBFMeshMotionSolver;
motionSolverLibs ("libRBFMeshMotionSolver.so");

RBFMeshMotionSolverCoeffs
{
    staticPatches    (upperWall lowerWall);
    movingPatches    (flap);
    fixedPatches     (inlet outlet);

    interpolation
    {
        function     TPS;
        polynomial   no;
        cpu          no;
        fullCPU      no;
    }

    coarsening
    {
        enabled                 yes;
        tol                     0.05;
        minPoints               20;
        maxPoints               200;
        livePointSelection      yes;
        tolLivePointSelection   0.05;
        exportSelectedPoints    no;
        twoPointSelection       no;
        surfaceCorrection       no;
    }
}
```

When switching this tutorial to `RBFMeshMotionSolver`, remove or rename
`0/fluid/pointMotionU`. If that field is present, the FSI mesh-motion code treats
the case as a finite-volume mesh-motion case and does not pass the interface
motion to the RBF solver.

For foam-extend builds, the entries inside `RBFMeshMotionSolverCoeffs` are read
directly from `dynamicMeshDict`, so omit the wrapper dictionary.

These settings were checked with the `perpendicularFlap` case on OpenFOAM v2312
for a short transient run. They are not intended to be universal optimum values;
increase `maxPoints`, tighten the tolerances, and inspect mesh quality for larger
or more distorted meshes.

## Patch Lists

`movingPatches` are the boundary patches whose displacement drives the mesh
motion. In the flap case this is the fluid-side FSI patch, `flap`.

`staticPatches` are included as zero-motion control points. Use these for walls
or other boundaries that must anchor the interpolation.

`fixedPatches` are enforced back to zero after interpolation. Use these for
boundaries such as inlet and outlet patches where the mesh must remain fixed.

Do not put the same patch in more than one of these lists. Empty, wedge, and
processor patches should normally be omitted.

## Interpolation Options

`function` selects the radial basis function. Valid values are `TPS`,
`WendlandC0`, `WendlandC2`, `WendlandC4`, and `WendlandC6`. The Wendland
functions require a positive `radius` entry in the `interpolation` dictionary.
`TPS` does not use `radius`.

`polynomial` adds a constant and linear polynomial term to the interpolation
system. This can improve rigid-body-like motion reproduction, but it increases
the dense system size.

`cpu` stores the factorisation of the control-point system and evaluates the
point interpolation matrix at each interpolation. With `cpu no`, the solver
stores the full point interpolation matrix. `fullCPU yes` requires `cpu yes` and
keeps more of the interpolation work out of the precomputed matrix path.

## Coarsening Options

Coarsening reduces the number of boundary control points used by the dense RBF
solve. This is useful for large boundary meshes, where using every moving and
static face centre can become expensive.

`enabled` switches coarsening on or off. With `enabled no`, all moving and static
control points are used.

`tol` is the relative error tolerance used by the greedy point selection. The
selection stops when both the global 2-norm error and maximum point error are
below this value, after at least `minPoints` have been selected, or when
`maxPoints` is reached.

`minPoints` and `maxPoints` bound the selected control-point count. For a larger
mesh, raise `maxPoints` if the reported coarsening error remains high.

`livePointSelection` recomputes the selected control points from the accumulated
motion during the run. This is usually preferable for FSI cases because the
important control points may change as the interface deforms.

`tolLivePointSelection` controls when live point selection is reused. If the
current selected points interpolate the accumulated boundary motion below this
relative tolerance, the solver reuses them; otherwise it selects a new set.

`exportSelectedPoints` writes the selected point coordinates to
`rbf-coarsening-greedy-selection-*.txt` and `.csv` files for inspection.

`twoPointSelection` may add a second point in the opposite error direction during
each greedy selection step. This can improve selection for motions with opposing
directions, at the cost of selecting points faster.

`surfaceCorrection` applies an additional WendlandC2 correction based on the
coarse interpolation error at the boundary. `ratioRadiusError` controls the
correction radius as a multiple of the maximum boundary interpolation error and
defaults to `10.0` when `surfaceCorrection` is enabled but `ratioRadiusError` is
not specified.
