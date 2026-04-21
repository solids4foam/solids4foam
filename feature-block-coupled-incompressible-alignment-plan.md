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
- [x] Align and rename `HolzapfelGasserOgdenElastic`.
  - [x] Rename without the `pd` prefix.
  - [x] Keep `TypeName` and file/class names consistent.
  - [x] Document required direction fields and dictionary entries.
  - [x] Factor duplicated volume/surface stress code into shared helpers where
    practical. (Added `correctF()` private helper; surface `correct()` now
    delegates to it, reducing duplication by ~65 lines.)
  - [x] Keep direct `p`/`pf` coupling explicit if required by the coupled
    pressure-displacement solver.

## Plate-Hole Boundary Conditions and Validation

- [~] Replace or shrink `setPlateHoleBC`.
  - Treat it as a case-specific boundary-condition mutator and validation
    helper, not as a general function object. ✓ (it already is; no code removed)
  - Reuse `analyticalPlateHoleTraction` for boundary tractions. DEFERRED —
    requires extending that BC to support
    `tractionPressureDisplacementFvPatchVectorField`, which is out of scope.
  - Reuse or extend `plateHoleAnalyticalSolution` for analytical error fields
    and sampled outputs. DEFERRED (same reason).
  - Avoid duplicate Kirsch/plate-hole analytical kernels in separate classes.
    Added a comment in `setPlateHoleBC.C` documenting the duplication as a
    future TODO.

## Tutorials

- [x] Move `tutorials/solids/coupledPressureDisplacement` cases into existing
  `tutorials/solids` families.
  - [x] cylinder, cylinderLin, cylinderUnsteady →
    `hyperelasticity/cylindricalPressureVessel/pressureDisplacement/`
  - [x] plateHole →
    `linearElasticity/plateHole/pressureDisplacement/`
    (compressible/ and incompressible/ subtrees preserved)
  - [x] heartTissueBeam, ventricleSymm →
    `biomechanics/cardiac/`
  - [x] HGO/ratCarotid →
    `biomechanics/vascular/HGO/ratCarotid/`
  - [x] `tutorials/solids/coupledPressureDisplacement/` removed.
  - [x] New `tutorials/solids/biomechanics/README.md` added.
- [x] Update `Allrun`, `Allclean`, `README.md`, and regression hooks after
  moving cases.
  - Added `caseOnlyRunsWithFoamExtend` guard to `ventricleSymm/Allrun` and
    `heartTissueBeam/Allrun` (the others already had it).

## Verification

- [x] Build with OpenFOAM-v2512:
  `. ~/OpenFOAM/OpenFOAM-v2512/etc/bashrc && ./Allwmake -j`
  (re-verified after HGO `correctF` refactor — exits 0)
- [x] Decide supported OpenFOAM targets for `coupledPressureDisplacementSolid`.
  Decision: FE-only for this PR. Added explanatory comment to
  `src/solids4FoamModels/Make/files.foamextend`.
- [ ] Run targeted tutorial checks for the changed pressure-displacement cases.
- [ ] Make clean `/tmp` copies for FE41 and OF9 builds where supported.
- [ ] Run relevant regression tests before finalizing.
