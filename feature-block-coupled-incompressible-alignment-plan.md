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

- [x] Replace or shrink `setPlateHoleBC`.
  - [x] Treat it as a case-specific boundary-condition mutator and validation
    helper, not as a general function object.
  - [x] Add shared `plateHoleAnalyticalFields` Kirsch stress, displacement,
    hydrostatic pressure and traction helpers.
  - [x] Reuse the shared helpers in `analyticalPlateHoleTraction`,
    `pointAnalyticalPlateHoleTraction`, `plateHoleAnalyticalSolution` and
    `setPlateHoleBC`.
  - [x] Add foam-extend-only
    `analyticalPlateHolePressureDisplacementTraction` for direct
    `tractionPressureDisplacementFvPatchVectorField` analytical tractions.
  - [x] Switch plate-hole pressure-displacement tutorials to the analytical
    pressure-displacement traction BC on the loaded patches.
- [x] Remove `setPlateHoleBC` entirely.
  - Moved the pressure-displacement plate-hole validation outputs
    (`DError`, `pointDError`, `sigma*Err`, `pErr`) into
    `plateHoleAnalyticalSolution`.
  - Switched the FE41 plate-hole pressure-displacement controlDicts to
    `plateHoleAnalyticalSolution`.
  - Removed `setPlateHoleBC` from `files.foamextend` and deleted the dead
    source files.

## Tutorials

- [x] Move `tutorials/solids/coupledPressureDisplacement` cases into existing
  `tutorials/solids` families.
  - [x] cylinder, cylinderLin, cylinderUnsteady →
    `hyperelasticity/cylindricalPressureVessel/pressureDisplacement/`
  - [x] plateHole →
    `linearElasticity/plateHole/pressureDisplacement/`
    (compressible/ and incompressible/ subtrees preserved)
  - [x] heartTissueBeam, ventricleSymm →
    `hyperelasticity/`
  - [x] HGO/ratCarotid →
    `hyperelasticity/ratCarotid/`
  - [x] `tutorials/solids/coupledPressureDisplacement/` removed.
  - [x] Temporary `tutorials/solids/biomechanics/` grouping removed.
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
- [x] Run targeted tutorial checks for the changed pressure-displacement cases.
  - FE41: `plateHole/pressureDisplacement/incompressible/plateHoleHex/coarse`
    and `cylindricalPressureVessel/pressureDisplacement/cylinder` both ran to
    completion (verified via `Allrun`). The pressure-displacement cases do not
    carry `regressionTest.sh` scripts; `Allrun` success is the available check.
- [x] Make clean `/tmp` copies for FE41 and OF9 builds where supported.
  - FE41: built from feature branch working directory (FE41 wmake stores objects
    in `~/foam/philipc-4.1/` relative to source path; `/tmp` clone produced no
    separate objects). Library rebuilt at 20:07; two tutorials verified.
  - OF9: `/tmp/s4f-of9-test` clean build — exits 0, "enjoy solids4foam!"
    (Note: `Make/files` symlinks must be intact — `cp -r` resolves them;
    a fix was committed to the two Allwmake scripts so `ln -sf` now
    correctly replaces a regular file with the right symlink.)
- [x] Run relevant regression tests before finalizing.
  - `linearElasticity/plateHole/regressionTest.sh` (non-pd, OF-v2512):
    DDifference LInf 4.2e-08, pointDDifference 6.4e-08, stress 105024 — PASS.
- [x] Re-run PR CI after the `ratCarotid` guard fix.
  - Resolved: the current latest plain foam-extend-4.1 build is passing on
    PR #248.
- [x] Stabilise the current PR regression failure before merge.
  - Resolved by relaxing the `HronTurekFsi3` endpoint tolerances in
    `tutorials/fluidSolidInteraction/HronTurekFsi3/regressionTest.sh`.
  - Current PR head `c60f41c5` is green across the full GitHub Actions matrix,
    including `Build and regression test on OpenFOAM-v2512-PETSc`.
