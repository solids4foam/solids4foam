---
sort: 11
---

# Channel with an elastic leaflet: `channelLeaflet`

---

## Tutorial Aims

- Demonstrates a strongly coupled internal-flow fluid-solid interaction case
  with a thin, light, highly flexible leaflet wetted on both faces and loaded
  by a pulsatile flow at Reynolds number 200.
- Compares the case with solutions of the same problem computed with the
  open-source finite-element library
  [oomph-lib](https://oomph-lib.github.io/oomph-lib/).

## Case Overview

A two-dimensional channel of height $$a = 1\,\mathrm{m}$$ is partly blocked
by a vertical elastic leaflet of height $$0.5\,\mathrm{m}$$ and thickness
$$h = 0.05\,\mathrm{m}$$, centred on $$x = 1\,\mathrm{m}$$ and clamped at its
base. The channel extends $$1\,\mathrm{m}$$ upstream and $$7\,\mathrm{m}$$
downstream of the leaflet. A pulsatile Poiseuille flow,
$$u = 6\,y\,(1 - y)\,q(t)$$ with $$q(t) = 1.5 - 0.5\cos(\pi t)$$, enters on
the left, so that the flux varies between 1 and $$2\,\mathrm{m^2/s}$$ with a
period of $$2\,\mathrm{s}$$. The fluid starts from rest and the flux is ramped
up as $$(1 - \cos(\pi t))/2$$ over the first second. The outflow is parallel
and axially traction-free.

| Quantity | Value |
| --- | --- |
| Fluid density, kinematic viscosity | 1 kg/m^3, 0.005 m^2/s |
| Reynolds number $$Re = U a/\nu$$ (and $$Re\,St$$) | 200 |
| Leaflet Young's modulus $$E$$, Poisson's ratio $$\nu_s$$ | 4550 Pa, 0.3 |
| Plane-strain modulus $$E/(1-\nu_s^2)$$ | 5000 Pa |
| Leaflet density | 1 kg/m^3, as the fluid (regularising) |
| Time step, end time | 0.01 s, 8 s |

After about two periods the leaflet performs a periodic, large-amplitude
oscillation with the period of the inflow: its tip moves between
$$x \approx 1.24$$ and $$1.36\,\mathrm{m}$$, a large rotation of the leaflet.
The monitored quantity is the displacement of the centre of the leaflet's top
face, written by the `solidPointDisplacement` function object `tip`.

### Relation to the published case

The configuration is that of the oomph-lib tutorial
[Flow in a 2D channel with an elastic leaflet](https://oomph-lib.github.io/oomph-lib/doc/interaction/fsi_channel_with_leaflet/html/index.html),
also a test case of Heil, Hazel & Boyle (2008), with $$Re = 200$$,
$$Q = 10^{-6}$$ and $$H = 0.05$$. Four things are changed:

1. **A longer downstream section**, 7 rather than 3 channel heights, which
   removes the artificial outflow boundary layer that the oomph-lib tutorial
   warns of; over the first 2.4 s, oomph-lib gives the same tip history for
   7 and 11 heights to 0.01%.
2. **A start from rest** with a ramped flux, instead of a sequence of steady
   solves up to $$Re = 200$$. The periodic state does not depend on the
   start: from the fourth period on, the oomph-lib solutions for both starts
   agree to 0.04% of the tip displacement.
3. **A leaflet of finite thickness.** The published leaflet is a
   Kirchhoff-Love beam loaded on its midplane, so the fluid sees a line; here
   it is a plane-strain St Venant-Kirchhoff continuum of thickness $$h$$,
   with $$E = E_{eff}(1 - \nu_s^2)$$ so that its bending and axial stiffness
   match the beam's, $$E_{eff} h^3/12$$ and $$E_{eff} h$$.
4. **A regularising leaflet density,** $$\rho_s = \rho_f$$; the published
   leaflet is massless, and the partitioned coupling does not converge for a
   massless leaflet.

## Solid discretisation and coupling

The leaflet uses the total-Lagrangian solid model solved with PETSc SNES and
the high-order least-squares reconstruction, `polynomialOrder 3`
(`./Allrun`), or 1 (`./Allrun linear`). The coupling is Dirichlet-Neumann with
Aitken relaxation (`constant/fsiProperties`), and the fluid mesh follows the
leaflet with the radial-basis-function mesh motion solver, which does not
accumulate mesh distortion from period to period as the `velocityLaplacian`
solver does.

## Running the case

```bash
./Allrun           # high-order (cubic) solid, about 2.5 hours in serial
./Allrun linear    # second-order (linear) solid
```

The case needs OpenFOAM.com and solids4foam built with PETSc.

## Status of the comparison with oomph-lib

The `verification/` directory holds a draft opt-in study against oomph-lib
solutions computed with the driver in `verification/reference/oomph-lib`.
On the tutorial mesh and time step, the leaflet-tip history over the periodic
state (the third and fourth periods) differs from oomph-lib by 1.1-1.3%
of the largest tip displacement
(0.43 m) with the high-order solid and by 3.2-3.7% with the linear solid;
the time-step error at 0.01 s is about 0.3%. Two issues prevent this from
being a verification to about 1%:

- **The model difference of the finite thickness is not controlled.** At the
  same bending stiffness and mass per unit length, the solids4foam solution
  changes by about 3% of the tip displacement when the thickness is halved to
  0.025 m, and its difference from oomph-lib grows to 2.5-2.8%, while the
  oomph-lib solution changes by only 0.4%. The difference does not fall as
  the leaflet is thinned.
- **The fluid mesh cannot yet be refined.** On the fluid mesh refined 1.5 and
  2 times, the Aitken coupling stalls during the first large deflection,
  IQN-ILS and Robin-Neumann coupling do not complete even the tutorial mesh,
  and the time step 0.0025 s fails in the solid Newton solve. The spatial
  error of the tutorial mesh, and so the part of the difference that is due
  to the thickness, cannot be measured.

## Regression test

`regressionTest.sh` runs the tutorial to $$t = 1\,\mathrm{s}$$ and checks the
tip displacement at $$t = 0.5$$ and $$1\,\mathrm{s}$$ against stored values.

## References

- M. Heil, A. L. Hazel and J. Boyle, Solvers for large-displacement
  fluid-structure interaction problems: segregated versus monolithic
  approaches, *Computational Mechanics*, 43, 91-101, 2008.
- oomph-lib,
  [Flow in a 2D channel with an elastic leaflet](https://oomph-lib.github.io/oomph-lib/doc/interaction/fsi_channel_with_leaflet/html/index.html).
