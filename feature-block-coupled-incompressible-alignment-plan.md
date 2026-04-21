# feature-block-coupled-incompressible Alignment Plan

This plan records the cleanup needed before merging PR #248 from
`feature-block-coupled-incompressible` into `development`.

## Goals

- Remove duplicated mechanical-law code where the pressure-displacement
  variants are forks of existing solids4foam laws.
- Keep the new coupled pressure-displacement capability, but express it through
  existing solids4foam/OpenFOAM design patterns.
- Avoid introducing new top-level tutorial families for cases that are variants
  of existing benchmark cases.
- Keep the cleanup split into reviewable patches.

## Mechanical Laws

- [x] Merge `pdNeoHookeanElastic` behavior into `neoHookeanElastic`.
  - [x] Add an explicit pressure-displacement mode to `neoHookeanElastic`.
  - [x] Preserve existing `neoHookeanElastic` behavior by default.
  - [x] Update tutorials from `pdNeoHookeanElastic` to
    `neoHookeanElastic`.
  - [x] Remove `pdNeoHookeanElastic` from build lists once no references
    remain.
- [x] Align `pdGuccioneElastic` with `GuccioneElastic`.
  - [x] Reuse the existing fibre-field creation logic where possible.
  - [x] Avoid duplicating the Guccione invariant/stress assembly.
  - [x] Fix the duplicated `R_`/`Rf_` component assignments.
  - [x] Decide whether the pressure-displacement behavior is a mode of
    `GuccioneElastic` or a thin derived/sibling law.
- [ ] Align and rename `HolzapfelGasserOgdenElastic`.
  - [x] Rename without the `pd` prefix.
  - [x] Keep `TypeName` and file/class names consistent.
  - [ ] Document required direction fields and dictionary entries.
  - [ ] Factor duplicated volume/surface stress code into shared helpers where
    practical.
  - [x] Keep direct `p`/`pf` coupling explicit if required by the coupled
    pressure-displacement solver.

## Plate-Hole Boundary Conditions and Validation

- [ ] Replace or shrink `setPlateHoleBC`.
  - Treat it as a case-specific boundary-condition mutator and validation
    helper, not as a general function object.
  - Reuse `analyticalPlateHoleTraction` for boundary tractions.
  - Reuse or extend `plateHoleAnalyticalSolution` for analytical error fields
    and sampled outputs.
  - Avoid duplicate Kirsch/plate-hole analytical kernels in separate classes.

## Tutorials

- [ ] Move `tutorials/solids/coupledPressureDisplacement` cases into existing
  `tutorials/solids` families.
  - Put pressure-vessel/cylinder cases under the existing cylinder or
    hyperelasticity/linear-elasticity families as appropriate.
  - Add pressure-displacement as a solution-algorithm variant when the physical
    case already exists.
  - Merge `plateHole` additions into the existing `plateHole` benchmark layout
    instead of maintaining a separate duplicate tree.
  - Keep genuinely new biomechanics cases under an appropriate existing
    hyperelasticity/biomechanics-style grouping, adding a README if a new group
    is necessary.
- [ ] Update `Allrun`, `Allclean`, `README.md`, and regression hooks after
  moving cases.

## Verification

- [x] Build with OpenFOAM-v2512:
  `. ~/OpenFOAM/OpenFOAM-v2512/etc/bashrc && ./Allwmake -j`
- [ ] Decide supported OpenFOAM targets for `coupledPressureDisplacementSolid`.
  It currently builds only through `files.foamextend`; adding it to
  `files.openfoam` exposes foam-extend-only dependencies such as
  `fvBlockMatrix`, GGI patches, and older point-patch APIs.
- [ ] Run targeted tutorial checks for the changed pressure-displacement cases.
  A v2512 plate-hole copy currently stops before law construction because
  `coupledPressureDisplacementSolid` is not registered in the OpenFOAM build.
- [ ] Make clean `/tmp` copies for FE41 and OF9 builds where supported.
- [ ] Run relevant regression tests before finalizing.
