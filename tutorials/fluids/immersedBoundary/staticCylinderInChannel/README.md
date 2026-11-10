---
sort: 1
---

# Immersed cylinder in a channel: `staticCylinderInChannel`

Prepared by Sairam Pamulaparthi Venkata and Philip Cardiff

## Tutorial Aims

- Demonstrate how to represent a body with the immersed boundary finite volume
  option `immersedBoundaryForce`, instead of meshing it;
- Compare the drag and lift coefficients with the laminar flow around a
  cylinder benchmark of Schäfer and Turek (1996).

## Case Overview

This is the steady 2-D benchmark case 2D-1 of Schäfer and Turek (1996): a
cylinder of diameter $$D = 0.1$$ m, centred at $$(0.2, 0.2)$$ m, in a channel
of height $$H = 0.41$$ m and length $$2.2$$ m. The inlet velocity profile is
parabolic,

$$
U_x = 4 U_m \frac{y (H - y)}{H^2}
$$

with the maximum velocity $$U_m = 0.3$$ m/s, so the mean velocity is
$$\bar{U} = 0.2$$ m/s. The outlet pressure is zero, and the channel walls are
no-slip walls. The fluid is Newtonian with a kinematic viscosity
$$\nu = 10^{-3}$$ m$$^2$$/s and density $$\rho = 1$$ kg/m$$^3$$, giving a
Reynolds number $$\bar{U} D/\nu = 20$$.

The cylinder is not meshed: the channel mesh (`system/blockMeshDict`) is
uniform upstream of $$x = 0.6$$ m and the cylinder is given as a closed surface,
`constant/triSurface/cylinder.stl`, which extends beyond the mesh in the
$$z$$ direction. The immersed boundary library is loaded in
`system/controlDict`:

```c++
libs            (immersedBoundary);
```

and the cylinder is added to the momentum equation of the `pimpleFluid` fluid
model in `constant/fvOptions`:

```c++
immersedBoundary
{
    type            immersedBoundaryForce;

    method          cutLink;

    bodies
    {
        cylinder
        {
            surface     "cylinder.stl";
        }
    }
}
```

With the `cutLink` method, the option penalises, implicitly in the momentum
equation, the velocity of the cells whose centre is inside the cylinder
towards the velocity of the cylinder, here zero, and that of each fluid cell
next to them towards the velocity of the cylinder at the surface, with a rate
from the distance to the surface along the lines between the cell centres, so
that the viscous flux through these faces is that of a linear profile through
the velocity of the cylinder on the surface. The method is described in
`src/immersedBoundary/README.md`. The force on the cylinder is written every
time step to `postProcessing/immersedBoundary/0/cylinder.dat`, and the solid
volume fraction, target velocity and forcing fields are written as
`immersedBoundary:lambda`, `immersedBoundary:Ui` and `immersedBoundary:f`.

The drag and lift coefficients are

$$
C_d = \frac{2 F_x}{\rho \bar{U}^2 D L_z}, \qquad
C_l = \frac{2 F_y}{\rho \bar{U}^2 D L_z}
$$

where $$L_z = 0.1$$ m is the thickness of the mesh, i.e. $$C_d = 5000 F_x$$.

## Running the Case

The `Allrun` script optionally takes the forcing method:

```bash
./Allrun              # The method in constant/fvOptions, cutLink by default
./Allrun cutLink      # Sharp interface cut-link penalisation
./Allrun ghostCell    # Sharp interface ghost cell method
./Allrun penalty      # Diffuse interface volume penalisation
./Allrun incremental  # Direct forcing of pimpleHFDIBFoam
./Allrun --help
```

With a method, `constant/fvOptions` is linked to `constant/fvOptions.<method>`,
which holds the settings of that method, and `nOuterCorrectors` in
`system/fvSolution` is set to 5 for `incremental` and 1 for the other methods,
as in the verification study (`verification/`); without one, the case is run
with the current `constant/fvOptions` and `system/fvSolution`. `./Allclean`
links `constant/fvOptions` back to `fvOptions.cutLink` and sets
`nOuterCorrectors` back to 1.

The mesh density is set by the `MESH_LEVEL` environment variable (default 1),
where each level halves the cell size, e.g.

```bash
MESH_LEVEL=2 ./Allrun
MESH_LEVEL=2 ./Allrun ghostCell
```

The time step is adjusted to give a maximum Courant number of 0.5. The flow is
steady after about 5 s. If gnuplot is installed, `Allrun` plots the force
coefficients in `forceCoeffs.pdf`.

## Expected Results

| `MESH_LEVEL` | Cells | Cells across $$D$$ | $$C_d$$ | $$C_l$$ |
| ------------ | ----- | ------------------ | ------- | ------- |
| 1 | 4 510 | 10 | 5.61 | 0.007 |
| 2 | 18 040 | 20 | 5.58 | 0.010 |
| 3 | 72 160 | 40 | 5.58 | 0.011 |
| Schäfer and Turek (1996) | | | 5.57–5.59 | 0.0104–0.0110 |

The drag is within 0.6% of the reference from 10 cells across the cylinder,
and within the reference bounds from 20 cells. The results depend little on
the time step, the number of pressure correctors or the number of processors:
with `MESH_LEVEL=1`, fixed time steps of 0.002 s and 0.0005 s give
$$C_d = 5.607$$ and $$C_d = 5.600$$, one pressure corrector gives the same
$$C_d$$, and with `MESH_LEVEL=2`, 1 and 4 processors give the same $$C_d$$
to the written precision. Columns 11-13 of the force file hold the force from
the surface traction, which, with the default `tractionQuadrature clustered`,
gives $$C_d$$ = 5.35, 5.46 and 5.49 and $$C_l$$ = 0.0026, 0.0078 and 0.0084:
the reconstructed wall traction varies on the scale of the cells, and the
lift remains low on these meshes. With `tractionQuadrature greedy;` (the
default of earlier versions), the lift is -0.060, -0.063 and -0.048, with a
sign set by the order of the facets of the STL file.

The other forcing methods of `immersedBoundaryForce` (see
`src/immersedBoundary/README.md`) give:

- `./Allrun ghostCell`: $$C_d$$ = 5.42, 5.51 and 5.55 for `MESH_LEVEL` 1, 2
  and 3;
- `./Allrun penalty`: $$C_d$$ = 5.24, 5.38 and 5.48;
- `./Allrun penalty`, with `weighting occupancy;` and
  `occupancy vertexFraction;` added to `constant/fvOptions.penalty`:
  $$C_d$$ = 5.45, 5.61 and 5.61;
- `./Allrun incremental` (`couplingCoeff 0.2;`, `occupancy vertexFraction;`
  and `nOuterCorrectors 5;`, the settings of the verification study):
  $$C_d$$ = 5.52, 5.70 and 5.72.

The `incremental` method is the direct forcing of the `pimpleHFDIBFoam`
solver from which the option is derived. Unlike the other methods, its results
also depend on the coupling coefficient, the number of outer correctors and
the time step; the verification study (`verification/README.md`) quantifies
this sensitivity.

## References

- Schäfer, M., Turek, S. (1996). Benchmark computations of laminar flow around
  a cylinder. In: Hirschel, E.H. (ed.) Flow Simulation with High-Performance
  Computers II. Notes on Numerical Fluid Mechanics, 52, 547-566.
