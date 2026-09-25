---
sort: 1
---

# Immersed cylinder in a channel: `staticCylinderInChannel`

You can find the files for this tutorial under
[`tutorials/fluids/immersedBoundary/staticCylinderInChannel`](https://github.com/solids4foam/solids4foam/tree/master/tutorials/fluids/immersedBoundary/staticCylinderInChannel).

---

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

    method          penalty;

    bodies
    {
        cylinder
        {
            surface     "cylinder.stl";
        }
    }
}
```

The option penalises the difference between the velocity of the cells covered
by the cylinder and the velocity of the cylinder, here zero, implicitly in the
momentum equation. The method is described in
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

```bash
./Allrun
```

The mesh density is set by the `MESH_LEVEL` environment variable (default 1),
where each level halves the cell size, e.g.

```bash
MESH_LEVEL=2 ./Allrun
```

The time step is adjusted to give a maximum Courant number of 0.5. The flow is
steady after about 5 s. If gnuplot is installed, `Allrun` plots the force
coefficients in `forceCoeffs.pdf`.

## Expected Results

| `MESH_LEVEL` | Cells | Cells across $$D$$ | $$C_d$$ | $$C_l$$ |
| ------------ | ----- | ------------------ | ------- | ------- |
| 1 | 4 510 | 10 | 5.24 | 0.004 |
| 2 | 18 040 | 20 | 5.38 | 0.009 |
| 3 | 72 160 | 40 | 5.48 | 0.010 |
| Schäfer and Turek (1996) | | | 5.57–5.59 | 0.0104–0.0110 |

The drag converges towards the reference at first order in the cell size.
The results do not depend on the penalty coefficient, the number of pressure
correctors or the time step: with `MESH_LEVEL=1`, fixed time steps of 0.01 s
and 0.0025 s give $$C_d = 5.236$$ and $$C_d = 5.239$$.

The other forcing methods of `immersedBoundaryForce` can be compared by
editing `constant/fvOptions` (see `src/immersedBoundary/README.md`):

- `weighting occupancy;` and `occupancy vertexFraction;`: $$C_d$$ = 5.45,
  5.61 and 5.61 for `MESH_LEVEL` 1, 2 and 3;
- `method incremental;`, `couplingCoeff 0.8;` and
  `occupancy vertexFraction;`: $$C_d$$ = 6.04, 6.42 and 6.60.

The `incremental` method is the direct forcing of the `pimpleHFDIBFoam`
solver from which the option is derived. Its results also depend on the
coupling coefficient, the time step and the number of pressure correctors: for
example, with `MESH_LEVEL=1`, a coupling coefficient of 0.2 gives
$$C_d = 5.99$$, and a fixed time step of 0.00125 s gives $$C_d = 5.79$$.

## References

- Schäfer, M., Turek, S. (1996). Benchmark computations of laminar flow around
  a cylinder. In: Hirschel, E.H. (ed.) Flow Simulation with High-Performance
  Computers II. Notes on Numerical Fluid Mechanics, 52, 547-566.
