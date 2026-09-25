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

    // Coupling coefficient, optional, default 0.8
    couplingCoeff   0.2;

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
occupancy `lambda`: half the fraction of the cell vertices inside the surface,
plus a half if the cell centre is inside the surface. The option adds a
forcing `f` (an acceleration) to the momentum equation. After each pressure
corrector, the forcing in the covered cells is increased by
`couplingCoeff*(Ui - U)/deltaT`, where `Ui` is the velocity of the body, and
the updated forcing is used by the next pressure corrector. At the start of
each time step, the bodies are moved to the new time and the forcing is
multiplied by the new occupancy.

The force `-rho*sum(f*V)` and torque on each body are written every time step
to `postProcessing/<option name>/<start time>/<body>.dat`, together with the
inertia of the fluid inside the body, `rho*d/dt(sum(lambda*U*V))`. For a
moving body, the hydrodynamic force is the sum of the two (Uhlmann, 2005).

The fields `<option name>:lambda`, `<option name>:Ui` and `<option name>:f`
are written at write times; `<option name>:f` is read on restart.

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

With the same settings, `immersedBoundaryForce` reproduces the forces of
`pimpleHFDIBFoam` for the static cylinder in a channel to six significant
figures.

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
