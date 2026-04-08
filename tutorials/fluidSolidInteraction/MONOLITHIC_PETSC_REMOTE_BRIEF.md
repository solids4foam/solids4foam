# Monolithic PETSc Remote Run Brief

Date: 2026-04-08

## Target Tree

- branch: `feature-petsc-snes-quasi-monolithic`
- sync the latest pushed tip of this branch before benchmarking
- this PETSc campaign currently depends on the branch tip including these files:
  - `src/solids4FoamModels/fluidModels/newtonIcoFluid/newtonIcoFluid.C`
  - `src/solids4FoamModels/fluidModels/newtonIcoFluid/newtonIcoFluid.H`
  - `src/solids4FoamModels/fluidSolidInterfaces/newtonQuasiMonolithicCouplingInterface/newtonQuasiMonolithicCouplingInterface.C`
  - `tutorials/fluidSolidInteraction/cavityFlexibleBottom/system/decomposeParDict.monolithic`
  - `tutorials/fluidSolidInteraction/blobInTreacle/system/fvSolution`
  - `tutorials/fluidSolidInteraction/cavityFlexibleBottom/system/fvSolution.monolithic`
  - `tutorials/fluidSolidInteraction/beamInCrossFlow/system/fvSolution.monolithic`
  - `tutorials/fluidSolidInteraction/foilInWind/system/fvSolution.monolithic`
  - `tutorials/fluidSolidInteraction/foilInWind/petscOptions.active`
  - `tutorials/fluidSolidInteraction/membraneRoof/system/fvSolution.monolithic`
  - `tutorials/fluidSolidInteraction/membraneRoof/petscOptions.active`
  - `tutorials/fluidSolidInteraction/compareAsciiFieldL2.sh`
  - `tutorials/fluidSolidInteraction/MONOLITHIC_PETSC_TODO.md`
  - `tutorials/fluidSolidInteraction/MONOLITHIC_PETSC_RESULTS.md`
  - `tutorials/fluidSolidInteraction/MONOLITHIC_PETSC_REMOTE_BRIEF.md`

## Current Local Winner

Use this as the working general preset candidate:

- top-level `fluid/solid` split
- nested fluid `U/p` Schur split
- Schur settings:
  `pc_fieldsplit_schur_fact_type lower`,
  `pc_fieldsplit_schur_precondition selfp`
- fluid velocity block:
  `bjacobi + ilu(0)`
- fluid pressure block:
  `hypre boomeramg` with:
  `max_iter 1`, `strong_threshold 0.6`, `grid_sweeps_up/down 1`,
  `agg_nl 1`, `agg_num_paths 1`, `max_levels 25`, `coarsen_type HMIS`,
  `interp_type ext+i`, `P_max 1`, `truncfactor 0.3`
- solid block:
  `bjacobi + lu`
- tolerances:
  `snes_rtol 1e-4`, `snes_stol 1e-4`, `ksp_rtol 1e-4`,
  `ksp_max_it 1000`, `ksp_type lgmres`, `ksp_gmres_restart 200`
- dimension-aware fluid block split:
  - 2-D: `block_size 3`, fields `0,1` / `2`
  - 3-D: `block_size 4`, fields `0,1,2` / `3`

## Local Validation Status

This preset is now the best local general candidate.

- `blobInTreacle`, `cavityFlexibleBottom`, `beamInCrossFlow`:
  passed serial and MPI screening
- `foilInWind`:
  passed the short fully coupled acceptance run at essentially the same cost as
  the previous foil-style `0.7` preset
- tightened field validation:
  see `MONOLITHIC_PETSC_RESULTS.md`
- field-difference helper:
  `tutorials/fluidSolidInteraction/compareAsciiFieldL2.sh`

## Requested Remote Work

1. Reproduce the local winning preset on a clean scratch workflow for:
   `blobInTreacle`, `cavityFlexibleBottom`, `beamInCrossFlow`, `foilInWind`.
2. Compare against these references:
   the checked-in monolithic baseline on `blobInTreacle` and
   `cavityFlexibleBottom`, the monolithic `bjacobi+lu` reference on
   `beamInCrossFlow`, and the current foil-style `0.7` split on `foilInWind`.
3. Run MPI scaling for the winning preset at the largest sensible rank counts
   for the remote machine.
   Suggested starting points:
   - `blobInTreacle`: `2, 4, 8`
   - `cavityFlexibleBottom`: `2, 4, 8`
   - `beamInCrossFlow`: `2, 4, 8, 16`
   - `foilInWind`: `2, 4, 8, 16, 32`
4. Record for every run:
   Newton iterations, total linear iterations, wall time, failure mode if any,
   and whether memory pressure or setup cost becomes the real limiter.
5. On at least one serial and one representative MPI run per case, compare
   `fluid/U`, `fluid/p`, `solid/U`, and `solid/sigmaEq` against the chosen
   reference using `compareAsciiFieldL2.sh`.
6. Report one scalar observable per case:
   final fluid force vector, and `Max sigmaEq` from the solver log.

## Deliverable Format

Return a markdown summary with one table row per run containing:

- case
- ranks
- preset
- Newton iterations
- linear iterations
- wall time
- status
- notes

Then add a short conclusion section stating:

- whether `strong_threshold 0.6` still looks like the best general preset
- whether a different preset is better at higher rank counts
- whether any case shows a real scaling or robustness regression
