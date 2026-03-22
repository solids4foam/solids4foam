---
sort: 5
---

# Tutorial: `perpendicularFlap`

This tutorial demonstrates a two-dimensional fluid-solid interaction between a steady channel flow and a rigid flap that protrudes perpendicularly from the lower wall. The fluid is incompressible and Newtonian (`rho = 1`, `nu = 1`) with a uniform inflow of `(10, 0, 0)` m/s. The solid is a linear-elastic flap (`rho = 3000 kg/m³`, `E = 4e6 Pa`, `nu = 0.3`) that is free to rotate about its fixed root and lies on the `flap` interface patch. The `frontAndBack` patch is marked `empty`, so the simulation is purely two-dimensional.

The coupling interface uses `solidPatch interface` and `fluidPatch flap`, while the kinematic and dynamic conditions are enforced by the partitioned `solids4Foam` solver. Strong coupling occurs early in the unsteady response, which is why most benchmarking focuses on the first second of the simulation even though the control dictionary extends to `t = 10 s`.

## Coupling options

This case ships with a tuned IQNILS configuration as the default (`./Allrun`). The new `Allrun` script can also run the original Aitken under-relaxation scheme and control whether the solver runs serially or in parallel:

- `./Allrun` (default): IQNILS coupling.
- `./Allrun aitken`: original Aitken relaxation scheme.
- `./Allrun parallel`: run in parallel (can be combined with either coupling option).

The script uses `foamDictionary` to update `constant/fsiProperties` on the fly and resets the dictionary back to the IQNILS default after the run. This means you can easily reproduce either algorithm without manually editing dictionaries.

## Running the case

1. Source your OpenFOAM environment (so that `solids4Foam` and `foamDictionary` are on `PATH`) and ensure PETSc is available.
2. Optional: reduce `system/controlDict` → `endTime` to `1.0` if you only want to probe the strongly-coupled transient window.
3. Run `./Allrun [coupling] [parallel]`.

`solids4Foam` will build the solid and fluid meshes, execute the coupling loop, and, if requested, run in parallel with automatic decompose/reconstruct steps.

More background on the available coupling methods and `fluidSolidInterface` options can be found in `src/solids4FoamModels/fluidSolidInterfaces/README.md`.
