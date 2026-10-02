---
sort: 3
---

# Pressure-driven laminar channel flow: `poiseuilleChannel`

You can find the files for this tutorial under
[`tutorials/fluids/poiseuilleChannel`](https://github.com/solids4foam/solids4foam/tree/master/tutorials/fluids/poiseuilleChannel).

---

## Tutorial Aims

- Demonstrate how to add finite volume options (`fvOptions` in OpenFOAM.com,
  `fvModels` and `fvConstraints` in OpenFOAM.org) to a `pimpleFluid` analysis;
- Verify the `pimpleFluid` fluid model against the analytical solution for
  plane Poiseuille flow.

## Case Overview

A Newtonian fluid (kinematic viscosity $$\nu = 0.1$$ m$$^2$$/s) flows between
two parallel no-slip walls a distance $$h = 1$$ m apart. The domain is 1 m long
and periodic (`cyclic`) in the flow direction, so there is no inlet or outlet.
Instead, the flow is driven by the `meanVelocityForce` option in
`constant/fvOptions`, which adds the uniform streamwise pressure gradient that
gives a prescribed mean velocity $$\bar{U} = 1$$ m/s. The fluid starts from
rest and the case is run until the flow is steady.

The steady solution is the parabolic velocity profile

$$
u(y) = 6 \bar{U} \frac{y}{h} \left( 1 - \frac{y}{h} \right)
$$

with the maximum velocity $$1.5 \bar{U}$$ on the channel centreline, driven
by the kinematic pressure gradient

$$
-\frac{1}{\rho}\frac{\partial p}{\partial x} = \frac{12 \nu \bar{U}}{h^2} = 1.2
\text{ m/s}^2
$$

The meanVelocityForce option reports the pressure gradient in the solver log
on every time step. On the 10 $$\times$$ 40 cell mesh, the steady pressure
gradient is 1.1985 m/s$$^2$$ and the maximum velocity is 1.4981 m/s, within
0.13% of the analytical values.

## Running the Case

```bash
./Allrun
```

The `regressionTest.sh` script runs the case and checks the final pressure
gradient, mean velocity and maximum velocity against the analytical solution.

OpenFOAM.org reads `constant/fvOptions` as `fvConstraints`, with a warning, so
the same case runs with both OpenFOAM.com and OpenFOAM.org. The case does not
run with foam-extend, as `pimpleFluid` does not support finite volume options
there.
