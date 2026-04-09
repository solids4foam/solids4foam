# Monolithic PETSc TODO

Date: 2026-04-08

## Scope

This is the active to-do list for finding a robust and efficient PETSc solver
strategy for the monolithic FSI solver on the
`feature-petsc-snes-quasi-monolithic` branch.

Primary goal:

- find a solver/preconditioner approach that is robust across most or all
  monolithic PETSc tutorial cases
- keep the work grounded in the current `v2512` branch reality
- update documentation as soon as the branch state changes materially

## Current Status

- [x] Reviewed the previous Codex and Claude session summaries against the
  current branch state
- [x] Confirmed that the current branch already contains the important
  split-lifecycle, serial `MatConvert`, and parallel solid direct-solver guard
  fixes in
  `newtonQuasiMonolithicCouplingInterface`
- [x] Identified that the old `foilInWind.codex/PETSC_TUNING_SUMMARY.md` note
  is historically useful but stale in several places
- [x] Identified that the checked-in `foilInWind` variants are not appropriate
  as the primary fast debugging case
- [x] Restored `pressureScaleFactor` handling from the local Codex branch into
  the current working tree
- [x] Re-gated per-block residual diagnostics behind
  `-s4f_petsc_report_block_residuals`
- [x] Built the branch successfully after the local code recoveries
- [x] Ran a short serial `blobInTreacle` scratch case with
  `pressureScaleFactor 3.0` to confirm that the recovered pressure scaling
  path executes and that block residual diagnostics remain quiet by default
- [x] Completed the first short serial baseline screen on
  `blobInTreacle`, `cavityFlexibleBottom`, and `beamInCrossFlow`
- [x] Completed short parallel baseline checks for
  `beamInCrossFlow` and `cavityFlexibleBottom`
- [x] Found and fixed a real monolithic parallel wrapper bug in
  `cavityFlexibleBottom/system/decomposeParDict.monolithic`
- [x] Completed a first serial `fluid/solid` split comparison on
  `blobInTreacle`, `cavityFlexibleBottom`, and `beamInCrossFlow`
- [x] Identified that the first `beamInCrossFlow` and `foilInWind` nested
  fluid Schur transplants into `blobInTreacle` and `cavityFlexibleBottom`
  were misconfigured because they kept the 3-D fluid `block_size 4` split on
  2-D cases
- [x] Re-ran the `foilInWind` split family with a dimension-aware 2-D nested
  fluid split (`block_size 3`, velocity fields `0,1`, pressure field `2`) on
  `blobInTreacle` and `cavityFlexibleBottom`
- [x] Found that the corrected dimension-aware `foilInWind` split preset is
  now competitive on the small 2-D cases and credible across
  `blobInTreacle`, `cavityFlexibleBottom`, `beamInCrossFlow`, and
  short-run `foilInWind`
- [x] Tuned the fluid-side pressure-AMG path inside the top-level
  `fluid/solid` family; the best current balance is the full foil-style
  BoomerAMG block with `strong_threshold 0.6`
- [x] Checked a fluid-block `ilu(0)` alternative inside the top-level
  `fluid/solid` family; it removes the `beamInCrossFlow` `lu` cost spike but
  degrades `blobInTreacle` and `cavityFlexibleBottom` too much to advance as
  the next cross-case candidate
- [x] Confirmed that the in-sandbox PRTE socket-bind failures were not
  solids4foam case bugs by reproducing them with `mpirun -np 2 hostname` and
  clearing them outside the sandbox
- [x] Revalidated the short parallel baseline matrix on
  `blobInTreacle`, `cavityFlexibleBottom`, and `beamInCrossFlow` outside the
  sandbox
- [x] Confirmed that the plain monolithic `bjacobi+lu` setup also converges on
  `beamInCrossFlow` in serial and parallel, but remains substantially slower
  than the tuned beam split baseline
- [x] Ran the slow acceptance case on short fully coupled `foilInWind`
  scratch copies
- [x] Confirmed that the current tuned split preset completes a one-step
  `foilInWind` acceptance run, while the plain monolithic `bjacobi+lu`
  reference fails its first linear solve with a PETSc out-of-memory error
- [x] Decided that the plain monolithic `bjacobi+lu` setup is a useful
  small-case robustness reference, but not a universal preset for the larger
  beam/foil-like cases
- [x] Compared the corrected dimension-aware `beamInCrossFlow` and
  `foilInWind` split variants directly on `blobInTreacle`,
  `cavityFlexibleBottom`, `beamInCrossFlow`, and short-run `foilInWind`
- [x] Chose the dimension-aware `foilInWind` split preset as the better
  current general split-family baseline: it remains competitive on the small
  cases, viable on `beamInCrossFlow`, and faster than the simpler beam-style
  pressure-AMG variant on short-run `foilInWind`
- [x] Tested hybrid pressure-AMG variants between the foil-style and beam-style
  settings rather than revisiting the direct-LU monolithic baseline
- [x] Screened two off-axis hybrids on the fast serial matrix:
  `strong_threshold 0.5` with the full foil-style AMG block and
  `strong_threshold 0.7` with the simpler beam-style AMG block
- [x] Found that the full foil-style AMG block with
  `strong_threshold 0.6` is the best current local balance: it improves
  `blobInTreacle`, stays strong on `cavityFlexibleBottom`, clearly improves
  `beamInCrossFlow`, survives MPI on all three fast cases, and matches the
  existing short-run `foilInWind` acceptance cost
- [x] Validated the `strong_threshold 0.6` hybrid against write-enabled
  tightened references on `blobInTreacle`, `cavityFlexibleBottom`, and
  `beamInCrossFlow`, plus the current `foilInWind` split acceptance run
- [x] Promoted the full foil-style pressure-AMG block with
  `strong_threshold 0.6` as the current working general monolithic PETSc
  preset candidate
- [x] Updated the shipped monolithic tutorial defaults to use the current
  working preset in `blobInTreacle`, `cavityFlexibleBottom`,
  `beamInCrossFlow`, `foilInWind`, and `membraneRoof`
- [x] Completed the larger-machine remote scaling pass on xenosim for
  `blobInTreacle`, `cavityFlexibleBottom`, `beamInCrossFlow`, and
  `foilInWind`; confirmed clean MPI behaviour on the first three cases and
  found a real foil MPI solid-response regression at 4+ ranks
- [x] Completed the larger-machine MPI scaling assessment of the current
  preset candidate
- [ ] Resolve the `foilInWind` MPI solid-response regression before treating
  the current preset as a universal default recommendation
- [ ] Test whether tighter outer Newton tolerances reduce the foil MPI
  `sigmaEq` drift that now appears with both the `strong_threshold 0.6` and
  `0.7` presets

## Ground Rules

- Use clean scratch copies for benchmarking
- Use `coupled yes;` and `couplingStartTime 0;` for debugging and tuning runs
- Shorten `endTime` and output where possible to keep inner-loop runs cheap
- Prefer breadth first across fast direct-map-compatible monolithic cases
  before returning to `foilInWind`
- Defer large-machine scaling work until a candidate survives the local 2-D
  case matrix
- On this machine, run MPI cases outside the sandbox; in-sandbox PRTE socket
  binds fail before solver startup

## Benchmark Order

Phase 1: local debugging and solver screening

- [x] `blobInTreacle`
- [x] `cavityFlexibleBottom`
- [x] `beamInCrossFlow`

Phase 2: slower acceptance case

- [x] `foilInWind`

Phase 3: later only, after the 2-D matrix is stable

- [ ] prepare a monolithic `3dTube` variant

## Solver Family Order

- [ ] establish a tightened baseline on each case using the safest current
  configuration
- [x] test top-level `fluid/solid` split variants
- [x] test nested fluid `U/p` Schur variants on top of the surviving
  `fluid/solid` family
- [ ] keep top-level `velocity/pressure/solid` split as exploratory only until
  the simpler families are stable across multiple cases
- [ ] defer MUMPS, PCLSC, GAMG, and high-rank scaling until a cross-case
  winner exists locally

## Parameter Sweep Order

- [ ] top-level split family
- [ ] fluid sub-PC
- [ ] solid sub-PC wrapper
- [ ] `fluidSystemScaleFactor`
- [ ] `solidSystemScaleFactor`
- [ ] `pressureScaleFactor`
- [ ] Schur details such as `upper` vs `lower` and `selfp`

## Validation Rules

For a candidate to advance:

- [x] it converges in serial on a short run
- [x] it converges in parallel on a short run
- [x] it does not rely on repeated linear failures
- [x] it stays close to a tightened reference on solid response
- [x] it remains competitive enough to justify additional testing

Compare at least:

- [x] `fluid/U`
- [x] `fluid/p`
- [x] `solid/U`
- [x] a solid stress quantity such as `sigmaEq`
- [x] one case-specific scalar observable such as force or tip displacement

## Documentation Tasks

- [x] Create a single active to-do list for the monolithic PETSc work
- [x] Replace stale foil-specific recommendations with branch-accurate status
- [x] Record benchmark results in a compact cross-case table instead of only in
  case-local narrative notes
- [x] Create a remote-run brief when a candidate is ready for a larger machine

## Exit Criteria For Local Phase

- [x] one solver family works cleanly on `blobInTreacle`,
  `cavityFlexibleBottom`, and `beamInCrossFlow`
- [x] the same family remains credible on short `foilInWind` runs
- [x] the branch documentation reflects the actual current state
- [x] the next large-machine request can be handed to another Codex instance as
  an explicit run brief
