---
sort: 5
---

# Immersed Taylor-Couette flow: `immersedTaylorCouette`

You can find the files for this tutorial under
[`tutorials/fluids/immersedBoundary/immersedTaylorCouette`](https://github.com/solids4foam/solids4foam/tree/master/tutorials/fluids/immersedBoundary/immersedTaylorCouette).

---

## Tutorial Aims

- Verify the velocity and the torque for a rotating immersed body, against
  the exact solution of the circular Couette flow;
- Demonstrate the `solidBodyRotation` body motion and several immersed
  bodies with the `immersedBoundaryForce` finite volume option.

## Case Overview

The fluid fills the gap between two immersed bodies on a Cartesian mesh of
the square $$[-0.16, 0.16]^2$$ m: a cylinder of radius $$R_1 = 0.05$$ m
(`constant/triSurface/rotor.stl`), which rotates about its axis at
$$\omega = 1$$ rad/s with the `solidBodyRotation` motion, and a fixed
annulus between $$R_2 = 0.1$$ m and 0.25 m (`stator.stl`), which covers the
corners of the mesh. Both surfaces cut the cells arbitrarily. The kinematic
viscosity is $$\nu = 10^{-3}$$ m$$^2$$/s, so that the Reynolds number
$$\omega R_1 (R_2 - R_1)/\nu$$ is 2.5, and the flow is steady after about
5 s.

The azimuthal velocity is

$$
u_\theta = A r + \frac{B}{r}, \qquad
A = -\frac{\omega R_1^2}{R_2^2 - R_1^2}, \qquad
B = \frac{\omega R_1^2 R_2^2}{R_2^2 - R_1^2}
$$

and the torque on the rotor, over the mesh thickness $$L_z = 0.1$$ m, is

$$
T = -\frac{4 \pi \rho \nu \omega R_1^2 R_2^2}{R_2^2 - R_1^2} L_z
  = -4.189 \times 10^{-6} \text{ N m}
$$

for a unit density. The torque about the centre of rotation (`CofR`) of
each body is written in columns 5-7 of
`postProcessing/immersedBoundary/0/<body>.dat`.

## Running the Case

```bash
./Allrun
```

with `MESH_LEVEL=2 ./Allrun` etc. for the finer meshes, which halve the cell
size (32, 64 and 128 cells across the mesh, i.e. 10, 20 and 40 across the
rotor).

## Expected Results

The difference of the torque on the rotor from the exact torque, and the root
mean square difference of the azimuthal velocity in the gap from the exact
velocity, relative to $$\omega R_1$$:

| `MESH_LEVEL` | Torque | Velocity | Torque (`penalty`) | Velocity (`penalty`) |
| ------------ | ------ | -------- | ------------------ | -------------------- |
| 1 | +0.8% | 7.2e-3 | -25% | 9.4e-2 |
| 2 | +0.3% | 2.3e-3 | -16% | 5.9e-2 |
| 3 | +0.2% | 9.8e-4 | -10% | 3.6e-2 |

With the `cutLink` method, the walls of both bodies are on their surfaces, and
the torque and velocity converge at about first and second order from 1%
on the coarsest mesh. With the `penalty` method (`method penalty;` in
`constant/fvOptions`), the effective walls are displaced into the bodies,
which widens the gap and lowers the torque.
