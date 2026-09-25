---
sort: 2
---

# Immersed oscillating cylinder in a channel: `oscillatingCylinderInChannel`

You can find the files for this tutorial under
[`tutorials/fluids/immersedBoundary/oscillatingCylinderInChannel`](https://github.com/solids4foam/solids4foam/tree/master/tutorials/fluids/immersedBoundary/oscillatingCylinderInChannel).

---

## Tutorial Aims

- Demonstrate how to prescribe the motion of an immersed body with the
  immersed boundary finite volume option `immersedBoundaryForce`, on a static
  mesh;
- Compare the drag and lift coefficients with those of Wan and Turek (2006).

## Case Overview

Wan and Turek (2006) computed the 2-D laminar flow around a cylinder that
oscillates horizontally in a closed channel filled with fluid initially at
rest. The channel is $$[0, 2.2] \times [0, 0.41]$$ m, and the cylinder, of
diameter $$D = 0.1$$ m, is centred at

$$
\mathbf{x}_c(t) = \left(1.1 + A \sin(2 \pi t/T),\ 0.2\right)
$$

with the amplitude $$A = 0.25$$ m and the period $$T = 4$$ s. The fluid is
Newtonian with a kinematic viscosity $$\nu = 10^{-3}$$ m$$^2$$/s and density
$$\rho = 1$$ kg/m$$^3$$. All the channel walls are no-slip walls. The same
problem is solved with a moving body-fitted mesh in the
`oscillatingCylinderInChannel` case of
[fluid-benchmarks](https://github.com/solids4foam/fluid-benchmarks).

Here the mesh does not move: the cylinder is given as a closed surface,
`constant/triSurface/cylinder.stl`, which moves through a mesh that is uniform
in the band $$0.78 < x < 1.42$$ m that the cylinder sweeps. The motion is
prescribed in `constant/fvOptions`:

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

            motion
            {
                type        sinusoidalTranslation;
                amplitude   0.25;
                period      4;
                direction   (1 0 0);
            }
        }
    }
}
```

and the immersed boundary library is loaded in `system/controlDict` with
`libs (immersedBoundary);`. The method is described in
`src/immersedBoundary/README.md`.

The force on the cylinder is written every time step to
`postProcessing/immersedBoundary/0/cylinder.dat`: columns 2-4 are the force
exerted by the immersed boundary forcing, and columns 8-10 the inertia of the
fluid inside the cylinder, whose sum is the hydrodynamic force on the
cylinder. The drag and lift coefficients are

$$
C_d = \frac{2 F_x}{\rho U_{ref}^2 D L_z}, \qquad
C_l = \frac{2 F_y}{\rho U_{ref}^2 D L_z}
$$

where the reference velocity is the maximum velocity of the cylinder,
$$U_{ref} = 2 \pi A/T = 0.3927$$ m/s, and $$L_z = 0.1$$ m is the thickness of
the mesh, i.e. $$C_d = 1296.9 F_x$$. The coefficients of Wan and Turek (2006)
are in `verificationData`.

## Running the Case

```bash
./Allrun
```

The mesh density is set by the `MESH_LEVEL` environment variable (default 1),
where each level halves the cell size, e.g.

```bash
MESH_LEVEL=2 ./Allrun
```

The case runs for two periods. The time step is at most 0.0025 s, and is
reduced if needed to give a maximum Courant number of 0.5. If gnuplot is
installed, `Allrun` plots the force coefficients against those of Wan and
Turek (2006) in `forceCoeffs.pdf`.

## Expected Results

The root mean square difference between the drag coefficient and that of Wan
and Turek (2006), over the two periods ($$0.25 < t < 8$$ s), whose root mean
square is 2.16:

| `MESH_LEVEL` | Cells | Cells across $$D$$ | Forcing | With inertia |
| ------------ | ----- | ------------------ | ------- | ------------ |
| 1 | 5 084 | 10 | 0.48 | 0.12 |
| 2 | 20 336 | 20 | 0.40 | 0.07 |
| 3 | 81 344 | 40 | 0.38 | 0.08 |

Including the inertia of the fluid inside the cylinder, the drag agrees with
the reference to within about 5% of its root mean square from the coarsest
mesh. As the occupancy varies continuously with the position of the cylinder,
the force is smooth as the cylinder crosses the cells.

With the `incremental` forcing method of `pimpleHFDIBFoam` (`method
incremental;` and `occupancy vertexFraction;` in `constant/fvOptions`), the
differences including the inertia are 1.36, 0.53 and 0.49 for levels 1-3;
the inertia is noisy, as it jumps when cells enter or leave the cylinder.

## References

- Wan, D., Turek, S. (2006). Fictitious boundary and moving mesh methods for
  the numerical simulation of rigid particulate flows. Journal of
  Computational Physics, 222, 28-56.
  <https://doi.org/10.1016/j.jcp.2006.06.002>
