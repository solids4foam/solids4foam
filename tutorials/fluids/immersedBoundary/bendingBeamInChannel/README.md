---
sort: 6
---

# Immersed beam bending in a channel: `bendingBeamInChannel`

You can find the files for this tutorial under
[`tutorials/fluids/immersedBoundary/bendingBeamInChannel`](https://github.com/solids4foam/solids4foam/tree/master/tutorials/fluids/immersedBoundary/bendingBeamInChannel).

---

## Tutorial Aims

- Verify the force on a deforming immersed body against a body-fitted
  solution with a deforming mesh;
- Demonstrate the deforming body motions of the `immersedBoundaryForce`
  finite volume option.

## Case Overview

A beam of width 0.2 m and height $$H = 2.0295$$ m stands on the bottom wall
of a channel of length 24 m and height 4.1 m, at $$x_c = 12$$ m, in a laminar
flow with the parabolic inlet velocity profile of mean $$\bar{U} = 1$$ m/s and
the kinematic viscosity $$\nu = 0.101475$$ m$$^2$$/s, i.e. the Reynolds number
$$\bar{U} H/\nu = 20$$. The beam bends periodically about its base, with the
prescribed displacement of a point at the reference position
$$(x_0, y_0)$$ (the `customProfileBend` motion)

$$
\Delta x = A_x g(\xi) \sin \omega t, \qquad
\Delta y = -(x_0 - x_c) \frac{A_y}{H} g'(\xi) \sin \omega t
$$

with $$\xi = y_0/H$$, the quintic profile $$g(\xi) = 10 \xi^3 - 15 \xi^4 +
6 \xi^5$$, $$A_x = 0.4$$ m, $$A_y = 0.8$$ m and $$\omega = 2 \pi$$ rad/s: the
tip moves 0.4 m either side of the axis, and the sections rotate with the
slope of the axis. The motion is prescribed in `constant/fvOptions`:

```c++
immersedBoundary
{
    type            immersedBoundaryForce;

    forceEstimator  surfaceTraction;

    bodies
    {
        beam
        {
            surface     "beam.stl";
            CofR        (12 1.01475 0);

            motion
            {
                type        customProfileBend;
                amplitudeX  0.4;
                amplitudeY  0.8;
                period      1;
                yMin        0;
                yMax        2.0295;
                xCenter     12;
            }
        }
    }
}
```

The surface of the beam, `constant/triSurface/beam.stl`, written by
`makeBeamStl.py`, has 80 segments along its height, so that it can bend; the
velocity of the body at a point is interpolated from the vertices of the
nearest triangle. The mesh is uniform around the beam, $$8.5 < x < 15.5$$ m,
with about three cells across the beam for `MESH_LEVEL=1`.

The drag and lift coefficients are $$C_d = 2 F_x/(\rho \bar{U}^2 H L_z)$$ and
$$C_l = 2 F_y/(\rho \bar{U}^2 H L_z)$$, with the mesh thickness
$$L_z = 0.1$$ m, i.e. $$C_d = 9.855 F_x$$. The force on the beam is written
every time step to `postProcessing/immersedBoundary/0/beam.dat`: columns 2-4
are the force from the surface traction, and columns 11-13 that from the
momentum exchange.

The reference, in `verificationData/CdBodyFitted.dat`, is the solution of
`pimpleFoam` on a body-fitted mesh deformed with the beam by a Laplacian
motion solver, from the embedded beam case of Sairam Pamulaparthi Venkata
(solids4foam/cardiacFoam#20), with twice as many cells in each direction and
the time step 0.0025 s; the solutions with 1 and 1.5 times as many cells
differ from it by 1% of the root mean square of the drag coefficient.

## Running the Case

```bash
./Allrun
```

with `MESH_LEVEL=2 ./Allrun` etc. for the finer meshes. If gnuplot is
installed, `Allrun` plots the force coefficients against those of the
body-fitted solution in `forceCoeffs.pdf`.

## Expected Results

The root mean square difference of the drag coefficient from the body-fitted
solution over $$1 < t < 8$$ s, relative to the root mean square of the drag
coefficient (53.9):

| `MESH_LEVEL` | Cells across the beam | Surface traction | Momentum exchange |
| ------------ | --------------------- | ---------------- | ----------------- |
| 1 | 3 | 8.9% | 7.4% |
| 2 | 6 | 5.8% | 4.6% |
| 3 | 11 | 4.7% | 3.5% |

The difference is largest when the beam is bent furthest upstream, where the
immersed drag is less negative than the body-fitted drag. The sections of the
beam stretch on one side and shorten on the other, so the flux of the body
velocity out of the cells inside the beam is not zero; the `pimpleFluid`
fluid model does not impose continuity in these cells, otherwise this flux
leaks into the flow, and the difference no longer decreases with the mesh
size (10.2%, 8.0% and 7.4% for the surface traction). Without the rotation of
the sections (`amplitudeY 0`), the body velocity is solenoidal, and the
differences from the corresponding body-fitted solution are 8.5%, 4.8% and
4.3%.
