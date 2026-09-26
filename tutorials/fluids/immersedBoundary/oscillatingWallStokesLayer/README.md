---
sort: 4
---

# Stokes layer on an immersed oscillating wall: `oscillatingWallStokesLayer`

You can find the files for this tutorial under
[`tutorials/fluids/immersedBoundary/oscillatingWallStokesLayer`](https://github.com/solids4foam/solids4foam/tree/master/tutorials/fluids/immersedBoundary/oscillatingWallStokesLayer).

---

## Tutorial Aims

- Verify the velocity and the wall shear stress next to a moving immersed
  wall that is not aligned with the cells, against the exact solution of the
  Stokes second problem;
- Demonstrate the `cutLink` method of the `immersedBoundaryForce` finite
  volume option.

## Case Overview

A slab below $$y_w = 0.1234$$ m, given as a closed surface
(`constant/triSurface/slab.stl`) that extends beyond the mesh in the $$x$$
and $$z$$ directions, oscillates in its plane with the velocity

$$
U_x = U_0 \sin(\omega t)
$$

with $$U_0 = 0.1$$ m/s and $$\omega = 2 \pi$$ rad/s, in a fluid of kinematic
viscosity $$\nu = 7.854 \times 10^{-3}$$ m$$^2$$/s, so that the thickness of
the Stokes layer is $$\delta = \sqrt{2 \nu/\omega} = 0.05$$ m. The mesh is
uniform, with $$N = 20 \times 2^{L - 1}$$ cells across the height of 0.5 m
for `MESH_LEVEL` $$L$$, periodic in the $$x$$ direction, with slip walls at
the bottom and top. The wall of the slab is not at a cell face.

Once the start-up transient has decayed, the velocity is that of the Stokes
second problem,

$$
U_x = U_0 e^{-k s} \sin(\omega t - k s), \qquad k = 1/\delta
$$

where $$s = y - y_w$$, and the wall shear stress is

$$
\tau_w = -\nu k U_0 \left( \sin(\omega t) + \cos(\omega t) \right)
$$

for a unit density. The time step is 0.001 s with the second order backward
scheme, and the case runs for four periods.

## Running the Case

```bash
./Allrun
```

with `MESH_LEVEL=2 ./Allrun` etc. for the finer meshes. If gnuplot is
installed, `Allrun` plots the wall shear stress, from the force on the slab
divided by the wall area in the mesh, against the periodic solution in
`wallShearStress.pdf`.

## Expected Results

The root mean square differences from a body-fitted solution over the fourth
period, relative to $$U_0$$ for the velocity in the fluid cells within
$$6 \delta$$ of the wall, and to the root mean square wall shear stress,
from the momentum exchange and from the surface traction:

| Level | $$\delta/h$$ | Velocity | Shear (exchange) | Shear (traction) |
| ----- | ------------ | -------- | ---------------- | ---------------- |
| 1 | 2 | 9.4e-3 | 9.5% | 39% |
| 2 | 4 | 1.8e-3 | 5.5% | 5.1% |
| 3 | 8 | 5.8e-4 | 4.5% | 3.8% |
| 4 | 16 | 7.0e-5 | 4.2% | 0.7% |

The velocity converges at second order or better. The wall shear stress of
the momentum exchange, which is the force that the fluid receives, is that of
a one-sided difference between the wall and the first fluid cell, and
converges at first order; that of the surface traction, which is fitted to
the fluid velocity over several cells, is more accurate on the finer meshes.
With the `penalty` method, the wall shear stress differs by 28%, 14%, 6% and
4% for levels 1-4, and the velocity by 5.4e-2, 2.7e-2, 1.2e-2 and 7.6e-3:
the effective wall is displaced into the slab.
