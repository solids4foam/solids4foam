---
sort: 8
---

# Immersed Turek-Hron FSI2 benchmark: `immersedHronTurekFsi2`

You can find the files for this tutorial under
[`tutorials/fluids/immersedBoundary/immersedHronTurekFsi2`](https://github.com/solids4foam/solids4foam/tree/master/tutorials/fluids/immersedBoundary/immersedHronTurekFsi2).

---

## Tutorial Aims

- Demonstrate a fluid-solid interaction in which the solid interface patch
  drives an immersed body of the `immersedBoundaryForce` finite volume
  option, on a fixed fluid mesh, rather than a fluid patch on a moving mesh;
- Compare the immersed solution of the Turek and Hron [1] FSI2 benchmark with
  the body-fitted solution of the `HronTurekFsi3` tutorial set up with the
  FSI2 parameters, and with the benchmark values.

## Case Overview

The geometry is that of the `HronTurekFsi3` tutorial: a channel of 2.5 m by
0.41 m with a rigid cylinder of radius 0.05 m centred at (0.2, 0.2) m and an
elastic flag of 0.35 m by 0.02 m attached to its downstream side, with the
tip point A at (0.6, 0.2) m. The FSI2 parameters are used: fluid density
1000 kg/m$$^3$$ and kinematic viscosity 0.001 m$$^2$$/s, mean inlet velocity
1 m/s (parabolic profile, applied without a ramp), solid density
10000 kg/m$$^3$$, Young's modulus 1.4 MPa and Poisson's ratio 0.4
(Saint Venant-Kirchhoff, plane strain). The coupling is enabled at
$$t = 2$$ s, once the flow has developed around the rigid flag, as in the
body-fitted tutorial.

The solid is the flag of the body-fitted tutorial, meshed with 105 by 6
cells (`nonLinearGeometryTotalLagrangianTotalDisplacement` solid model,
PETSc SNES solution algorithm). The fluid is solved with `pimpleFluid` on a
fixed background mesh of the whole channel, uniform with 5 mm cells (four
cells across the flag) in $$0.1 < x < 0.7$$ m, $$0.1 < y < 0.3$$ m, and
coarsening towards the outlet, for `MESH_LEVEL=1`; each level halves the cell
size. The cylinder and the flag are immersed bodies of the
`immersedBoundaryForce` option (`cutLink` method) in
`constant/fluid/fvOptions`: the cylinder is a static body given by
`constant/fluid/triSurface/cylinder.stl` (written by `makeCylinderStl.py`),
and the flag is a body with the `fsiDriven` motion and no surface file:

```c++
bodies
{
    cylinder
    {
        surface     "cylinder.stl";
    }

    flag
    {
        motion
        {
            type        fsiDriven;
        }
    }
}
```

The fluid-solid interface (`constant/fsiProperties`, IQN-ILS coupling)
pairs the solid interface patch `plate` with `fluidPatch none;` and names
the body in the `immersedInterfaces` dictionary, with the clamped root of
the flag, `plateFix`, as a closure patch:

```c++
solidPatch      plate;
fluidPatch      none;

immersedInterfaces
{
    plate
    {
        body            flag;
        closurePatches  (plateFix);
    }
}
```

The surface of the flag body is built by the interface from the `plate` and
`plateFix` patches of the solid; as the fluid mesh is two-dimensional, the
surface is extended through the mesh and capped in the empty direction, so
that the patches of that direction are not needed. In each coupling
iteration the interface moves the surface to the relaxed (or IQN-ILS
accelerated) solid displacement, the bodies are re-positioned in the fluid
mesh, and the surface traction of the immersed boundary, averaged over the
quadrature points of each face of the `plate` patch, is applied to the solid.
The log prints the total traction force on the flag and, for comparison, the
momentum exchange between the fluid and the flag body, which also includes
the inertia of the fluid inside the body and the root of the flag.

The displacement of the tip point A is written to
`postProcessing/0/solidPointDisplacement_pointDisp.dat`, and the forces on
the bodies to `postProcessing/fluid/immersedBoundary/0/<body>.dat` (columns
2-4 the momentum exchange, 11-13 the surface traction).

## Expected Results

RESULTS_PLACEHOLDER

## Running the Case

The case is run with `./Allrun` (`./Allrun parallel` for four processors),
which creates the solid and fluid meshes with `blockMesh` and runs
`solids4Foam`. The fluid mesh density is set by the `MESH_LEVEL` environment
variable (default 1), e.g. `MESH_LEVEL=2 ./Allrun parallel`. The case
requires OpenFOAM.com (the `immersedBoundary` library) and solids4foam built
with PETSc.

## References

[1]
[Turek, S., Hron, J. (2006). Proposal for Numerical Benchmarking of Fluid-Structure Interaction between an Elastic Object and Laminar Incompressible Flow. In: Bungartz, HJ., Schäfer, M. (eds) Fluid-Structure Interaction. Lecture Notes in Computational Science and Engineering, vol 53. Springer, Berlin, Heidelberg.](https://doi.org/10.1007/3-540-34596-5_15)
