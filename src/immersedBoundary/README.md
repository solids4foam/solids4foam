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

Four forcing methods are available:

- `penalty` (the default): the implicit volume penalisation
  `kappa*(Ui - U)`, where `Ui` is the velocity of the body, is added to the
  momentum equation, with `-kappa*U` in the matrix. The penalisation rate is
  `K/deltaT` in the fully covered cells, where `K` is `penaltyCoeff` (default
  1e3). In the partially covered cells, `kappa*deltaT = lambda/(1 - lambda)`
  (`weighting volumeFraction`, the default), so that the penalised velocity is
  `lambda*Ui + (1 - lambda)*U`; or `kappa = K*lambda/deltaT` (`weighting
  occupancy`), which blocks every covered cell. With the `volumeFraction`
  weighting, the rate in the partially covered cells is also limited to
  `lambda/(1 - lambda)*C*(nu/w^2 + |Ui|/w)`, where `w` is the cell width and
  `C` is `surfaceRateCoeff` (default 3), so that it does not depend on the
  time step. The results do not depend on `K` for `K` of 1e3 or more, the
  number of pressure correctors or the time step: for the static cylinder
  with 10 cells across it, `Cd` is 5.236 and 5.239 for time steps of 0.01 s
  and 0.0025 s. With `surfaceRateCoeff 0`, the partially covered cells are
  relaxed towards `Ui` once per time step, and `Cd` is 5.37 and 5.49 for the
  same time steps.
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
  the body leaves (fresh cells) disturb the force: `cutLink` or `penalty`
  is recommended.
- `cutLink`: a sharp interface penalisation for static and moving bodies.
  The cells whose centre is inside a body are penalised towards the body
  velocity, with the rate `min(K/deltaT, C*(nu/w^2 + |Ub|/w))`, where `C` is
  `pinnedRateCoeff` (default 100). Each fluid cell next to them is penalised
  towards the body velocity at the surface with the rate
  `sum(nu*|Sf|*deltaCoeff*(|d|/phi - 1))/V` over its faces with a penalised
  neighbour, where `|d|` is the distance between the cell centres and `phi`
  the distance from the fluid cell centre to the surface along the line
  between them, from the intersection with the surface; with the explicit
  correction of the flux to the penalised cell (`linkCorrection`, default
  yes), the viscous flux through the face is that of a linear profile
  through the body velocity on the surface (Shortley-Weller), so that the
  wall is on the surface to second order. The rates do not depend on the
  time step and vary continuously as the body moves. With
  `apertureCoupling` (default yes), the option also registers the fluid
  fraction of the area of the faces cut by the bodies, from the geometric
  cutting of the faces (`cutFaceIso`), and the flux of the body velocity
  through their solid part; the `pimpleFluid` fluid model uses them in its
  pressure equation, so that it changes continuously as the body moves. The
  force is the momentum exchange (`forceEstimator momentumExchange`, the
  default), which includes the inertia of the fluid inside the body; the
  surface traction (`forceEstimator surfaceTraction`), from quadratic least
  squares fits of the pressure and velocity at quadrature points moving
  with the surface, is written in the columns after the inertia (or is the
  force, with the momentum exchange after the inertia). The momentum
  exchange assumes a laminar flow of constant viscosity on a static mesh.
- `incremental`: the direct forcing of the openHFDIB-DEM `pimpleHFDIBFoam`
  solver. An explicit forcing `f` is added to the momentum equation. After
  each pressure corrector, `f` is increased by `couplingCoeff*(Ui - U)/deltaT`
  in the covered cells, and the updated forcing is used by the next pressure
  corrector; at the start of each time step, `f` is multiplied by the new
  occupancy. The results depend on `couplingCoeff`, the time step and the
  number of pressure correctors. This method is kept for comparison: with
  `occupancy vertexFraction` it reproduces `pimpleHFDIBFoam`.

The force `-rho*sum(f*V)` and torque on each body, where `f` is the forcing
(an acceleration) exerted on the fluid, are written every time step to
`postProcessing/<option name>/<start time>/<body>.dat`, together with the
inertia of the fluid inside the body, `rho*d/dt(sum(lambda*U*V))`. For a
moving body, the hydrodynamic force is the sum of the two (Uhlmann, 2005).

The fields `<option name>:lambda`, `<option name>:Ui` and `<option name>:f`
are written at write times; with the `incremental` method,
`<option name>:f` is read on restart.

### Accuracy

For the two tutorials, with 10, 20 and 40 cells across the cylinder, the
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
| E | 5.61, 5.58, 5.58 | 0.19, 0.07, 0.06 | 0.17, 0.07, 0.06 |
| F | 5.36, 5.46, 5.48 | 0.11, 0.04, 0.02 | 0.12, 0.08, 0.08 |

where the settings are:

- A: the defaults (`penalty`, `volumeFraction`, `signedDistance`);
- B: `penalty` with `weighting occupancy` and `occupancy vertexFraction`;
- C: `incremental` with `occupancy vertexFraction`;
- D: `ghostCell`;
- E: `cutLink` (the setting of the static tutorial), with the force from the
  momentum exchange (the default);
- F: `cutLink`, with the force from the surface traction.

For an immersed plane Poiseuille flow between walls that are not aligned with
the cell faces, with 20, 40 and 80 cells across the channel, the flow rate
differs from the exact solution by 20%, 11% and 7% with setting A, by
0.9%, 0.03% and 0.02% with setting D, and by 0.25%, 0.23% and 0.05% with
setting E. For a cylinder translating at constant velocity through the mesh
(the `translatingCylinderInChannel` tutorial), `cutLink` gives the drag of the
static cylinder at the same position to within 0.2% on all the meshes, with
fluctuations of 1%, 0.4% and 0.3% as the cylinder crosses the cells, whereas
the `penalty` drag is low by 5%, 3% and 1.5%. For the immersed Stokes layer
(the `oscillatingWallStokesLayer` tutorial), the `cutLink` velocity converges
at second order. The momentum exchange is the force that the fluid receives,
and is the more accurate for static and steadily moving bodies; the surface
traction is the more accurate for accelerating bodies and for the wall shear
stress on fine meshes.

The cylinder forces converge at about first order in the cell size for all the
methods, except the static `cutLink` drag, which is within 0.6% of the
reference from 10 cells across the cylinder (the ghost cell method was not
run with 40 cells across the moving cylinder). The differences from Wan and
Turek (2006) stop decreasing at about 0.08, the difference between the
body-fitted mesh solution and Wan and Turek (2006), whose coefficients lag the
converged solutions by about 0.015 s (see the `oscillatingCylinderInChannel`
tutorial).

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
