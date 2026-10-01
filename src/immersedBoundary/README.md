# immersedBoundary

Hybrid fictitious domain-immersed boundary (HFDIB) method for bodies with a
prescribed motion immersed in an incompressible flow, as the finite volume
option `immersedBoundaryForce`. It is currently only available for
OpenFOAM.com, where `libimmersedBoundary` is built with solids4foam.

The option works with the solids4foam `pimpleFluid` fluid model and with the
standard `pimpleFoam` solver. See `tutorials/fluids/immersedBoundary`.

## Usage

In `system/controlDict`:

```c++
libs (immersedBoundary);
```

In `constant/fvOptions` (for a fluid-only solids4foam case), or
`constant/fluid/fvOptions` (for the fluid region of a fluid-solid interaction
case):

```c++
immersedBoundary
{
    type            immersedBoundaryForce;

    // Forcing method, optional, default penalty
    method          penalty;

    bodies
    {
        cylinder
        {
            // Closed surface, relative to constant/triSurface
            surface     "cylinder.stl";

            // Prescribed motion, optional, default static
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

The complete list of entries is in the header of
`fvOptions/immersedBoundaryForce/immersedBoundaryForce.H`.

## Method

Each body is a closed surface that covers some cells of the fluid mesh with an
occupancy (solid volume fraction) `lambda`. By default, `lambda` is found from
the signed distance of the cell centres to the surface, relative to the width
of the cells normal to the surface, so that it varies continuously as the
body moves (`occupancy signedDistance`). Alternatively, as in openHFDIB-DEM,
it is half the fraction of the cell vertices inside the surface, plus a half
if the cell centre is inside the surface (`occupancy vertexFraction`).

Three forcing methods are available:

- `penalty` (the default): the implicit volume penalisation
  `kappa*(Ui - U)`, where `Ui` is the velocity of the body, is added to the
  momentum equation, with `-kappa*U` in the matrix. The penalisation rate is
  `K/deltaT` in the fully covered cells, where `K` is `penaltyCoeff` (default
  1e3). In the partially covered cells, `kappa*deltaT = lambda/(1 - lambda)`
  (`weighting volumeFraction`, the default), so that the penalised velocity is
  `lambda*Ui + (1 - lambda)*U`; or `kappa = K*lambda/deltaT` (`weighting
  occupancy`), which blocks every covered cell. With the `volumeFraction`
  weighting, the rate in the partially covered cells is also limited to
  `lambda/(1 - lambda)*C*(nu/w^2 + max(|Ui|, |U|)/w)`, where `w` is the
  cell width and `C` is `surfaceRateCoeff` (default 3). The maximum preserves
  the body-motion scale and includes fluid advection for static bodies.
  This numerical velocity scale is relative to the fixed mesh, and the cap
  has no explicit time-step dependence. It is recomputed from the current
  fluid velocity at each momentum assembly; timestep and iteration
  independence must therefore be checked for each case. Local fluid speed
  is itself damped by the penalty, so this is not an estimate of incident
  flow speed or a guarantee of a Reynolds-independent boundary location.
  With `surfaceRateCoeff 0`, the partially covered cells are relaxed towards
  `Ui` once per time step.
- `ghostCell`: a sharp interface variant of `penalty`, for static bodies.
  The cells whose centre is inside a body are penalised with the rate
  `K/deltaT`. In the cells next to the fluid (the ghost cells), the target
  velocity is extrapolated linearly along the surface normal, through the
  body velocity at the nearest surface point and the velocity at an image
  point in the fluid, `imageDistance` (default 1.5) cell widths beyond the
  surface. The image point velocity is interpolated in its donor cell
  (`imageInterpolation cellPoint`, the default) or is the donor cell value
  (`imageInterpolation cell`). The donor may be on another processor. The
  fluid velocity then equals the body velocity on the surface, so the
  effective wall is not displaced into the body as it is with `penalty`.
  The results do not depend on `K`, the time step, the number of pressure
  correctors or the number of processors. For moving bodies, the cells that
  the body leaves (fresh cells) disturb the force: `penalty` is recommended.
- `incremental`: the direct forcing of the openHFDIB-DEM `pimpleHFDIBFoam`
  solver. An explicit forcing `f` is added to the momentum equation. After
  each pressure corrector, `f` is increased by `couplingCoeff*(Ui - U)/deltaT`
  in the covered cells, and the updated forcing is used by the next pressure
  corrector; at the start of each time step, `f` is multiplied by the new
  occupancy. The results depend on `couplingCoeff`, the time step and the
  number of pressure correctors. This method is kept for comparison: with
  `occupancy vertexFraction` it reproduces `pimpleHFDIBFoam`.
  `couplingCoeff` controls how much of the velocity error is corrected after
  each pressure corrector. The default `0.8` gives strong enforcement without
  the oscillation that full correction can produce. Reduce it when forcing
  or pressure corrections oscillate or diverge; a smaller value enforces
  body velocity more slowly and can require more pressure correctors or a
  smaller timestep for the same accuracy.

The force `-rho*sum(f*V)` and torque on each body, where `f` is the forcing
(an acceleration) exerted on the fluid, are written every time step to
`postProcessing/<option name>/<start time>/<body>.dat`, together with the
inertia of the fluid inside the body, `rho*d/dt(sum(lambda*U*V))`. For a
moving body, the hydrodynamic force is the sum of the two (Uhlmann, 2005).

The fields `<option name>:lambda`, `<option name>:Ui`, `<option name>:kappa`
and `<option name>:f`
are written at write times; with the `incremental` method,
`<option name>:f` is read on restart.

### Accuracy

Before adding fluid speed to the surface-rate cap, the two tutorials
with 10, 20 and 40 cells across the cylinder gave the following results. The
static cylinder drag coefficient (reference 5.57-5.59), and the root mean
square difference of the oscillating cylinder drag coefficient, including the
inertia of the fluid inside the cylinder, from a moving body-fitted mesh
solution and from Wan and Turek (2006), over 0.25 < t < 7.5 s (where the root
mean square of the drag coefficient is 2.05), are:

| Settings | Static `Cd` | Oscillating: body-fitted | Wan and Turek |
| -------- | ----------- | ------------------------ | ------------- |
| A | 5.24, 5.38, 5.48 | 0.11, 0.05, 0.02 | 0.13, 0.09, 0.08 |
| B | 5.45, 5.61, 5.61 | 2.00, 0.44, 0.29 | 1.42, 0.33, 0.20 |
| C | 6.04, 6.42, 6.60 | 1.73, 0.60, 0.53 | 1.32, 0.52, 0.48 |
| D | 5.42, 5.51, 5.55 | 0.80, 0.30, - | 0.59, 0.27, - |

where the settings are:

- A: `penalty`, `volumeFraction`, `signedDistance`, with the previous
  body-speed-only surface-rate cap;
- B: `penalty` with `weighting occupancy` and `occupancy vertexFraction`;
- C: `incremental` with `occupancy vertexFraction`;
- D: `ghostCell` (the setting of the static tutorial).

For an immersed plane Poiseuille flow between walls that are not aligned with
the cell faces, with 20, 40 and 80 cells across the channel, the flow rate
differs from the exact solution by 20%, 11% and 7% with setting A, and by
0.9%, 0.03% and 0.02% with setting D.

These results showed first-order accuracy in cell size for all the methods
(the ghost cell method was not run with 40 cells across the moving cylinder).
The differences from Wan and Turek (2006) stop decreasing at about 0.08, the
difference between the body-fitted mesh solution and Wan and Turek (2006),
whose coefficients lag the converged solutions by about 0.015 s (see the
`oscillatingCylinderInChannel` tutorial).

### Fluid-speed surface-rate check

With the maximum of body and fluid speed in the cap, OpenFOAM v2512 checks
with fixed timesteps gave static `Cd = 5.3027, 5.3944` at 10 and 20 cells
across the cylinder (at 10 s). On the coarse mesh, reducing the timestep
from 0.01 to 0.0025 s gave `Cd = 5.3029`, and three instead of one PIMPLE
outer correctors gave `Cd = 5.3027`.

For the oscillating cylinder, full 8 s runs at a fixed timestep of 0.0025 s
gave RMS drag differences from the body-fitted reference of `0.1039, 0.0454`
at 10 and 20 cells across the cylinder, over `0.25 < t < 7.5 s`, including
the inertia of the internal fluid. On the coarse mesh, a timestep of
0.000625 s gave `0.1363`, while three outer correctors gave `0.1018`.
These checks show measurable timestep sensitivity in the moving case;
absence of explicit timestep dependence in the cap does not imply timestep
independence of the solution. The finest mesh was not reverified with this
cap; the three-level table above records the earlier body-speed-only cap.

At one-tenth viscosity, the RMS fluid speed in the partially covered cells
at 10 s fell from 0.131 to 0.065 m/s compared with the body-speed-only cap.
This is a damping diagnostic, not a boundary-location or drag-accuracy test:
that flow is unsteady and has a different Reynolds number from the supplied
static reference. In a short zero-viscosity diagnostic, all 36 partially
covered cells had positive penalty rates, compared with zero for the old
cap. This does not establish inviscid no-slip accuracy.

## Provenance

This library is a rewrite, for prescribed motion and as a finite volume
option, of the immersed boundary library and `pimpleHFDIBFoam` solver
contributed to cardiacFoam in
[solids4foam/cardiacFoam#20](https://github.com/solids4foam/cardiacFoam/pull/20)
by Sairam Pamulaparthi Venkata (UCD). That code is based on:

- **openHFDIB-DEM**, <https://github.com/techMathGroup/openHFDIB-DEM>,
  developed mostly by members of the techMathGroup of the Institute of
  Thermomechanics of the Czech Academy of Sciences and the Monolith group of
  the University of Chemistry and Technology, Prague; its main contributors
  are Martin Isoz, Martin Kotouč Šourek, Ondřej Studeník and Petr Kočí. The
  imported version is the library in `src/HFDIBDEM` and the
  `pimpleHFDIBFoam` solver at commit
  [`b4ea321`](https://github.com/techMathGroup/openHFDIB-DEM/commit/b4ea321)
  (2024-08-19), on the `main` and `Ver2.6` branches.
- **openHFDIB**, <https://github.com/fmuni/openHFDIB>, by Federico Municchi,
  from which openHFDIB-DEM is derived.

The lineage of the code is:

1. openHFDIB-DEM `b4ea321`, imported into
   [xenosim-erc/immersedBoundaryRigidMotion](https://github.com/xenosim-erc/immersedBoundaryRigidMotion)
   by Sairam Pamulaparthi Venkata;
2. the discrete element (contact, virtual mesh) parts removed and the code
   ported to OpenFOAM v2512 by Philip Cardiff (February 2026, branch
   `sairam-leftHeart-IB-Verf`);
3. prescribed motion laws (sinusoidal translation, bending, heart valves) and
   benchmark cases added by Sairam Pamulaparthi Venkata (June 2026,
   `xenosim-erc/openHFDIB-Benchmarks-v2412`, and
   solids4foam/cardiacFoam#20 commit `251e887`);
4. rewritten here for solids4foam by Philip Cardiff: the discrete element
   method, free bodies, body addition, adaptive mesh refinement and restart
   files are removed; bodies are moved identically on every processor at the
   current time; and the solver is replaced by the finite volume option.

The occupancy (`immersedBody::addOccupancy`) follows the openHFDIB-DEM
`nonConvexBody`, and the forcing and force calculation
(`immersedBoundaryForce`, `immersedBody::force`) follow the openHFDIB-DEM
`pimpleHFDIBFoam` solver and `immersedBody`. The sinusoidal translation
follows the version of Sairam Pamulaparthi Venkata. The openHFDIB-DEM
interpolation of the velocity at the immersed boundary (`lineInt`,
`leastSquares`), which the benchmark cases do not use, is not included.

With the `incremental` method, the `vertexFraction` occupancy and the same
settings, `immersedBoundaryForce` reproduces the forces of `pimpleHFDIBFoam`
for the static cylinder in a channel to six significant figures. The
`penalty` method, the `signedDistance` occupancy and the `volumeFraction`
weighting are new in solids4foam.

## Licence

openHFDIB is distributed under the GNU Lesser General Public License version
3, and the openHFDIB-DEM repository under the GNU General Public License
version 3 (its source files state the GNU Lesser General Public License
version 3). The Lesser General Public License allows the code to be
distributed under the GNU General Public License version 3, under which this
library is distributed as part of solids4foam.

## References

- Municchi, F. openHFDIB. <https://doi.org/10.5281/zenodo.3940236>
- Isoz, M., Kotouč Šourek, M., Studeník, O., Kočí, P. (2022). Hybrid
  fictitious domain-immersed boundary solver coupled with discrete element
  method for simulations of flows laden with arbitrarily-shaped particles.
  Computers & Fluids, 244, 105538.
  <https://doi.org/10.1016/j.compfluid.2022.105538>
- Studeník, O., Isoz, M., Kotouč Šourek, M., Kočí, P. (2024). OpenHFDIB-DEM:
  An extension to OpenFOAM for CFD-DEM simulations with arbitrary particle
  shapes. SoftwareX, 27. <https://doi.org/10.1016/j.softx.2024.101871>
- Uhlmann, M. (2005). An immersed boundary method with direct forcing for the
  simulation of particulate flows. Journal of Computational Physics, 209,
  448-476.
