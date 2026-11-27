---
sort: 4
---

# Unsteady stagnation-point flow: `stagnationPointFlow`

You can find the files for this tutorial under
[`tutorials/fluids/stagnationPointFlow`](https://github.com/solids4foam/solids4foam/tree/master/tutorials/fluids/stagnationPointFlow).

---

## Tutorial Aims

- Demonstrate the `fluxConsistentPatches` PIMPLE entry of the OpenFOAM.com
  `pimpleFluid` fluid model on an outlet with a fixed pressure and a
  prescribed (`fixedGradient`) velocity gradient;
- Verify `pimpleFluid` against an exact unsteady solution of the
  incompressible Navier–Stokes equations.

## Case Overview

The flow towards a stagnation point at the origin, with a strain rate
$$a(t)$$ that varies in time,

$$
u = a(t)\, x, \qquad v = -a(t)\, y, \qquad
a(t) = 1 + \tfrac{1}{2}\sin(2 \pi t),
$$

is an exact solution of the incompressible Navier–Stokes equations for any
viscosity, as the viscous term is zero. The kinematic pressure is

$$
p = -\tfrac{1}{2}\left(a^2 + \dot{a}\right) x^2
    - \tfrac{1}{2}\left(a^2 - \dot{a}\right) y^2 .
$$

The domain is the unit square, $$0 \le x, y \le 1$$, with
$$\nu = 0.01$$ m$$^2$$/s, meshed with 20 $$\times$$ 20 cells. The case
starts from the exact solution and runs one period, to $$t = 1$$ s, with
$$\Delta t = 0.01$$ s:

| Boundary | Velocity | Pressure |
| --- | --- | --- |
| `left`, `bottom` | `symmetryPlane` | `symmetryPlane` |
| `top`, inflow | exact value | `fixedFluxPressure` |
| `outlet` | exact gradient | exact value |

The `left` and `bottom` boundaries are at $$x = 0$$ and $$y = 0$$, `top` is
at $$y = 1$$ and `outlet` at $$x = 1$$. The exact values use
`uniformFixedValue` with an `expression` function, and the exact velocity
gradient, $$\partial \mathbf{u}/\partial x = (a, 0)$$, uses
`uniformFixedGradient`.

The `outlet` is listed in `fluxConsistentPatches` in `system/fvSolution`.
Without it, the outlet flux in the pressure equation is that of the
boundary-cell velocity, which differs from the flux of the boundary velocity
by $$(h/2)\, a$$ per unit face area. This O($$h$$) flux error becomes a
pressure error in the outlet cells. With the entry, the increment
$$(\mathbf{u}_b - \mathbf{u}_P) \cdot \mathbf{S}_f$$ of the gradient condition
is added to the outlet flux.

The case runs without a momentum predictor, with two outer correctors and
two pressure correctors per time step. The errors are written at the final
time as the fields `UError` $$= |\mathbf{u} - \mathbf{u}_{exact}|$$ and
`pError` $$= p - p_{exact}$$.

## Results

Root-mean-square errors at $$t = 1$$ s, on 20 $$\times$$ 20 and
40 $$\times$$ 40 cells with $$\Delta t \propto h$$:

| Predictor | `fluxConsistentPatches` | p, 20 | p, 40 | U, 20 | U, 40 |
| --- | --- | --- | --- | --- | --- |
| off | none | 7.4e-2 | 3.1e-2 | 1.9e-2 | 5.7e-3 |
| off | `(outlet)` | 1.5e-2 | 6.0e-3 | 9.1e-3 | 2.2e-3 |
| on | none | 8.7e-2 | 3.6e-2 | 1.9e-2 | 5.6e-3 |
| on | `(outlet)` | 7.8e-3 | 2.4e-3 | 8.7e-3 | 1.9e-3 |

With the outlet in `fluxConsistentPatches`, the pressure error falls by a
factor of 5–15 and the velocity error by a factor of 2–3.

## Running the Case

```bash
./Allrun
```

The `regressionTest.sh` script runs the case and compares the largest and the
mean velocity and pressure errors at $$t = 1$$ s with reference values. It
also checks that an `inletOutlet` outlet velocity, which neither fixes the
value nor is a `fixedGradient`, is rejected at start-up.

The case runs with OpenFOAM.com only: `fluxConsistentPatches` is implemented
in the OpenFOAM.com `pimpleFluid`, and the boundary conditions use the
OpenFOAM.com expression functions.
