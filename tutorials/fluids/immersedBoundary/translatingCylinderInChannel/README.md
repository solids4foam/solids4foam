---
sort: 3
---

# Immersed cylinder translating in a channel: `translatingCylinderInChannel`

You can find the files for this tutorial under
[`tutorials/fluids/immersedBoundary/translatingCylinderInChannel`](https://github.com/solids4foam/solids4foam/tree/master/tutorials/fluids/immersedBoundary/translatingCylinderInChannel).

---

## Tutorial Aims

- Verify the force on an immersed body that moves through the mesh, against
  the same flow computed with the body at rest;
- Demonstrate the `cutLink` method of the `immersedBoundaryForce` finite
  volume option for moving bodies.

## Case Overview

This is the laminar flow around a cylinder benchmark 2D-1 of Schäfer and Turek
(1996), see the `staticCylinderInChannel` tutorial, seen from a frame moving
upstream at 0.05 m/s: the cylinder of diameter $$D = 0.1$$ m translates
downstream at 0.05 m/s, from $$(0.2, 0.2)$$ m at $$t = 0$$, and the channel
walls move at 0.05 m/s, with the inlet velocity profile

$$
U_x = 0.05 + 4 U_m \frac{y (H - y)}{H^2}
$$

with $$U_m = 0.3$$ m/s and $$H = 0.41$$ m. The flow relative to the cylinder
is then that of the benchmark with the cylinder at the position
$$x = 0.2 + 0.05 t$$ m, so that, once the flow has developed, the forces on
the translating cylinder equal those on a static cylinder at the same
position. These are in `verificationData/CdStatic.dat`, from static cylinders
at 0.05 m intervals on the finest mesh: the drag coefficient is 5.58 with the
cylinder at $$x = 0.2$$ m, as in the benchmark, and falls to 5.29 as the
cylinder moves away from the inlet, beyond $$x = 0.4$$ m ($$t = 4$$ s).

The cylinder moves through a mesh that is uniform upstream of $$x = 1.2$$ m,
with 10 cells across the cylinder for `MESH_LEVEL=1`, so the cells that it
covers change as it moves. Its motion is prescribed in `constant/fvOptions`:

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

            motion
            {
                type        uniformTranslation;
                velocity    (0.05 0 0);
            }
        }
    }
}
```

The drag and lift coefficients are $$C_d = 2 F_x/(\rho \bar{U}^2 D L_z)$$ and
$$C_l = 2 F_y/(\rho \bar{U}^2 D L_z)$$, with the mean velocity relative to the
cylinder $$\bar{U} = 0.2$$ m/s and the mesh thickness $$L_z = 0.1$$ m, i.e.
$$C_d = 5000 F_x$$. The force on the cylinder is written every time step to
`postProcessing/immersedBoundary/0/cylinder.dat`: with the `cutLink` method,
columns 2-4 are the force from the momentum exchange, which includes the
inertia of the fluid inside the cylinder, and columns 11-13 that from the
surface traction.

## Running the Case

```bash
./Allrun
```

with `MESH_LEVEL=2 ./Allrun` etc. for the finer meshes. If gnuplot is
installed, `Allrun` plots the force coefficients against those of the static
cylinder in `forceCoeffs.pdf`.

## Expected Results

The mean drag coefficient over $$4 < t < 6$$ s, where the reference is 5.29,
and the root mean square of its fluctuations about a linear fit, which come
from the cells that the cylinder crosses:

| `MESH_LEVEL` | Cells across $$D$$ | `cutLink` | `penalty` | `ghostCell` |
| ------------ | ------------------ | --------- | --------- | ----------- |
| 1 | 10 | 5.298 ± 0.053 | 5.034 ± 0.058 | 5.128 ± 0.62 |
| 2 | 20 | 5.294 ± 0.021 | 5.127 ± 0.025 | 5.187 ± 0.27 |
| 3 | 40 | 5.299 ± 0.014 | 5.210 ± 0.010 | 5.210 ± 0.07 |

The `cutLink` drag is within 0.2% of the static reference on all the meshes,
whereas the `penalty` drag is low by 5%, 3% and 1.5%, as the effective wall is
displaced into the body, and the `ghostCell` drag is disturbed as the cells
that the cylinder leaves become fluid. The fluctuations of the `cutLink` drag
decrease by about a factor of three for each refinement, and do not grow as
the time step decreases.

## References

- Schäfer, M., Turek, S. (1996). Benchmark computations of laminar flow around
  a cylinder. In: Hirschel, E.H. (ed.) Flow Simulation with High-Performance
  Computers II. Notes on Numerical Fluid Mechanics, 52, 547-566.
