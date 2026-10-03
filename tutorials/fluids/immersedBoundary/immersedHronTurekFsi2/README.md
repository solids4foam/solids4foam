---
sort: 8
---

# Immersed Turek-Hron FSI2 benchmark: `immersedHronTurekFsi2`

## Tutorial Aims

- Demonstrate a fluid-solid interaction in which the solid interface patch
  drives an immersed body of the `immersedBoundaryForce` finite volume
  option, on a fixed fluid mesh, rather than a fluid patch on a moving mesh;
- Compare the immersed solution of the Turek and Hron [1] FSI2 benchmark with
  the benchmark values, and with the body-fitted solution of the
  `HronTurek` tutorial with the FSI2 parameters (`./Allrun fsi2`) while the
  oscillation grows.

## Case Overview

The geometry is that of the `HronTurek` tutorial: a channel of 2.5 m by
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

The oscillation of the flag grows from the start of the coupling at
$$t = 2$$ s and saturates at about $$t = 8$$ s. The tables give the mean and
amplitude (half the difference between the extrema) over 2 s of the
saturated oscillation, $$13 < t < 15$$ s for `MESH_LEVEL=1` and
$$8 < t < 10$$ s for `MESH_LEVEL=2`, with the forces smoothed over 10 ms to
remove their fluctuations as the surface crosses the cells, and the
frequency from the crossings of the mean.

Tip displacement of the flag and frequency of $$u_y$$:

| | Tip $$u_y$$ (mm) | Tip $$u_x$$ (mm) | Frequency (Hz) |
| --- | --- | --- | --- |
| `MESH_LEVEL=1` | 1.17 ± 75.3 | -12.8 ± 11.7 | 2.07 |
| `MESH_LEVEL=2` | 1.23 ± 76.8 | -13.3 ± 12.0 | 2.05 |
| Turek and Hron [1] | 1.23 ± 80.6 | -14.58 ± 12.44 | 2.0 |

Drag and lift (N/m):

| | Drag | Lift |
| --- | --- | --- |
| `MESH_LEVEL=1` | 204.8 ± 73.6 | -2.5 ± 236.7 |
| `MESH_LEVEL=2` | 209.2 ± 77.9 | -0.9 ± 255.0 |
| Turek and Hron [1] | 208.83 ± 73.75 | 0.88 ± 234.2 |

The drag and lift are those of the cylinder and the flag, from the surface
traction, per unit depth. The amplitude of the tip displacement is 7% and 5%
low with `MESH_LEVEL=1` and 2, and the frequency 3% and 2% high. The mean
drag is within 2% of the benchmark on both meshes; the amplitudes of the drag
and lift are within 1% of it with `MESH_LEVEL=1` but 6% and 9% high with
`MESH_LEVEL=2`, and so do not yet converge with the mesh (the smoothed extrema
are sensitive to the fluctuations of the traction as the surface crosses the
cells).

While the oscillation grows, the tip displacement follows that of the
body-fitted solution closely (for $$3.5 < t < 4$$ s, $$u_y$$ from -0.7 to
4.5 mm immersed and from -0.9 to 4.6 mm body-fitted). The body-fitted
solution with the FSI2 parameters and the mesh of the `HronTurekFsi3`
tutorial (now `HronTurek`) diverged at about $$t = 6.5$$ s, before the
oscillation saturates (solids4foam issue #489), whereas the immersed solution,
with the same solid mesh and solid model, runs to $$t = 15$$ s.

The regression test runs `MESH_LEVEL=1` with the coupling started at
$$t = 0.5$$ s to $$t = 0.6$$ s, and checks the tip displacement and the
force on the flag.

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
