# Monolithic PETSc Results

Date: 2026-04-08

This table is the compact cross-case record for the active monolithic PETSc
screening work on `feature-petsc-snes-quasi-monolithic`.

## Code Recovery Checks

| Case | Mode | Scratch case | Settings | Newton | Linear iters | Exec time | Notes |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `blobInTreacle` | serial | `/tmp/blobInTreacleCheck.jE3rNv` | current monolithic setup, `endTime 0.025`, `pressureScaleFactor 3.0` | 2 | 69 | 0.19 s | recovered `pressureScaleFactor` path exercised; no block residual diagnostics printed by default |

## Breadth-First Screening

| Case | Mode | Scratch case | Settings | Newton | Linear iters | Exec time | Status | Notes |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| `blobInTreacle` | serial | `/tmp/monolithicPetscScreen.kUAclp/blobInTreacle` | current checked-in monolithic baseline, `endTime 0.025` | 2 | 49 | 0.16 s | complete | simple `bjacobi+lu` reference case |
| `cavityFlexibleBottom` | serial | `/tmp/monolithicPetscScreen.kUAclp/cavityFlexibleBottom` | current checked-in monolithic baseline, `endTime 0.5` | 2 | 28 | 0.12 s | complete | solver log finished cleanly; `Allrun` wrapper exited nonzero after solver/post-processing tail |
| `beamInCrossFlow` | serial | `/tmp/monolithicPetscScreen.kUAclp/beamInCrossFlow` | current checked-in monolithic setup, `endTime 0.05` | 2 | 27 | 1.9 s | complete | already uses top-level `fluid/solid` split with nested fluid Schur |
| `blobInTreacle` | parallel (`2`) | `/tmp/monolithicPetscScreen.kUAclp/blobInTreacle-par-rerun` | current checked-in monolithic baseline, `endTime 0.025` | 2 | 86 | 0.1 s | complete | runs cleanly outside the sandbox; earlier PRTE socket-bind failure was environmental rather than case-specific |
| `cavityFlexibleBottom` | parallel (`4`) | `/tmp/monolithicPetscScreen.kUAclp/cavityFlexibleBottom-par4` | current checked-in monolithic baseline, `endTime 0.5` | 2 | 69 | 0.1 s | complete | revalidated outside the sandbox after making the top-level monolithic `decomposeParDict` consistent with the region decomposition |
| `beamInCrossFlow` | parallel (`2`) | `/tmp/monolithicPetscScreen.kUAclp/beamInCrossFlow-par` | current checked-in monolithic setup, `endTime 0.05` | 2 | 37 | 1.05 s | complete | revalidated outside the sandbox; useful checked-in beam-style reference |
| `blobInTreacle` | serial | `/tmp/monolithicPetscScreen.kUAclp/blobInTreacle-beamSplit` | initial `beamInCrossFlow`-style top-level `fluid/solid` split transplant, `endTime 0.025` | 2 | 139 | 0.27 s | superseded | used the wrong 3-D nested fluid split (`block_size 4`) on a 2-D case, so this row should not be used for final comparison |
| `cavityFlexibleBottom` | serial | `/tmp/monolithicPetscScreen.kUAclp/cavityFlexibleBottom-beamSplit` | initial `beamInCrossFlow`-style top-level `fluid/solid` split transplant, `endTime 0.5` | 2 | 122 | 0.27 s | superseded | used the wrong 3-D nested fluid split (`block_size 4`) on a 2-D case, so this row should not be used for final comparison |
| `blobInTreacle` | serial | `/tmp/monolithicPetscScreen.kUAclp/blobInTreacle-fsSplit` | plain top-level `fluid/solid` split, fluid `bjacobi+lu`, solid `bjacobi+lu`, `endTime 0.025` | 2 | 54 | 0.17 s | complete | close to baseline, but slightly worse than the 49-iteration `bjacobi+lu` monolith |
| `cavityFlexibleBottom` | serial | `/tmp/monolithicPetscScreen.kUAclp/cavityFlexibleBottom-fsSplit` | plain top-level `fluid/solid` split, fluid `bjacobi+lu`, solid `bjacobi+lu`, `endTime 0.5` | 2 | 28 | 0.11 s | complete | essentially matches the baseline iteration count and runtime |
| `beamInCrossFlow` | serial | `/tmp/monolithicPetscScreen.kUAclp/beamInCrossFlow-fsSplit` | plain top-level `fluid/solid` split, fluid `bjacobi+lu`, solid `bjacobi+lu`, `endTime 0.05` | 2 | 24 | 20.66 s | complete | fewer iterations than the nested Schur baseline, but far slower because direct `lu` on the full fluid block is too expensive |
| `blobInTreacle` | serial | `/tmp/monolithicPetscScreen.kUAclp/blobInTreacle-fsSplit-fluidILU0` | plain top-level `fluid/solid` split, fluid `bjacobi+ilu(0)`, solid `bjacobi+lu`, `endTime 0.025` | 2 | 127 | 0.23 s | complete | removes direct fluid `lu`, but is much worse than both the baseline and the fluid-`lu` split |
| `cavityFlexibleBottom` | serial | `/tmp/monolithicPetscScreen.kUAclp/cavityFlexibleBottom-fsSplit-fluidILU0` | plain top-level `fluid/solid` split, fluid `bjacobi+ilu(0)`, solid `bjacobi+lu`, `endTime 0.5` | 2 | 111 | 0.2 s | complete | much worse than both the baseline and the fluid-`lu` split |
| `beamInCrossFlow` | serial | `/tmp/monolithicPetscScreen.kUAclp/beamInCrossFlow-fsSplit-fluidILU0` | plain top-level `fluid/solid` split, fluid `bjacobi+ilu(0)`, solid `bjacobi+lu`, `endTime 0.05` | 2 | 94 | 2.2 s | complete | far cheaper than fluid `lu` on the beam case, but still slower and higher-iteration than the checked-in nested-Schur beam baseline |
| `beamInCrossFlow` | serial | `/tmp/monolithicPetscScreen.kUAclp/beamInCrossFlow-monoLU` | plain monolithic `bjacobi+lu` baseline from the smaller cases, `endTime 0.05` | 2 | 36 | 11.26 s | complete | common conservative baseline does converge on beam, but with a large runtime penalty versus the tuned beam split |
| `beamInCrossFlow` | parallel (`2`) | `/tmp/monolithicPetscScreen.kUAclp/beamInCrossFlow-monoLU` | plain monolithic `bjacobi+lu` baseline from the smaller cases, `endTime 0.05` | 2 | 69 | 4.59 s | complete | common conservative baseline also converges in MPI, but remains clearly slower than the tuned beam split baseline |
| `blobInTreacle` | serial | `/tmp/monolithicPetscScreen.kUAclp/blobInTreacle-foilSplit-serial2d` | dimension-aware `foilInWind` split preset, 2-D nested fluid split (`block_size 3`, fields `0,1` / `2`), `endTime 0.025` | 2 | 55 | 0.13 s | complete | near the 49-iteration baseline and far better than the misconfigured 3-D transplant |
| `blobInTreacle` | parallel (`2`) | `/tmp/monolithicPetscScreen.kUAclp/blobInTreacle-foilSplit-par2d` | dimension-aware `foilInWind` split preset, 2-D nested fluid split (`block_size 3`, fields `0,1` / `2`), `endTime 0.025` | 2 | 76 | 0.12 s | complete | corrected split now runs cleanly in MPI and improves on the 86-iteration plain baseline |
| `cavityFlexibleBottom` | serial | `/tmp/monolithicPetscScreen.kUAclp/cavityFlexibleBottom-foilSplit-serial2d` | dimension-aware `foilInWind` split preset, 2-D nested fluid split (`block_size 3`, fields `0,1` / `2`), `endTime 0.5` | 2 | 30 | 0.12 s | complete | close to the 28-iteration small-case baseline |
| `cavityFlexibleBottom` | parallel (`4`) | `/tmp/monolithicPetscScreen.kUAclp/cavityFlexibleBottom-foilSplit-par2d` | dimension-aware `foilInWind` split preset, 2-D nested fluid split (`block_size 3`, fields `0,1` / `2`), `endTime 0.5` | 2 | 38 | 0.09 s | complete | substantially better than the 69-iteration plain baseline in MPI |
| `beamInCrossFlow` | serial | `/tmp/monolithicPetscScreen.kUAclp/beamInCrossFlow-foilSplit` | `foilInWind` split preset on the 3-D beam case, `endTime 0.05` | 2 | 36 | 1.7 s | complete | still slower in iterations than the checked-in beam split, but same order of runtime and clearly viable |
| `beamInCrossFlow` | parallel (`2`) | `/tmp/monolithicPetscScreen.kUAclp/beamInCrossFlow-foilSplit` | `foilInWind` split preset on the 3-D beam case, `endTime 0.05` | 2 | 46 | 1.0 s | complete | viable in MPI, though the checked-in beam split still has fewer linear iterations |
| `blobInTreacle` | serial | `/tmp/monolithicPetscScreen.kUAclp/blobInTreacle-beamSplit2d-serial` | dimension-aware `beamInCrossFlow` split variant, 2-D nested fluid split (`block_size 3`, fields `0,1` / `2`), `endTime 0.025` | 3 | 77 | 0.17 s | complete | much worse than the foil-style variant on this case |
| `blobInTreacle` | parallel (`2`) | `/tmp/monolithicPetscScreen.kUAclp/blobInTreacle-beamSplit2d-par` | dimension-aware `beamInCrossFlow` split variant, 2-D nested fluid split (`block_size 3`, fields `0,1` / `2`), `endTime 0.025` | 3 | 106 | 0.16 s | complete | substantially worse than the foil-style variant in MPI too |
| `cavityFlexibleBottom` | serial | `/tmp/monolithicPetscScreen.kUAclp/cavityFlexibleBottom-beamSplit2d-serial` | dimension-aware `beamInCrossFlow` split variant, 2-D nested fluid split (`block_size 3`, fields `0,1` / `2`), `endTime 0.5` | 2 | 14 | 0.12 s | complete | excellent on this case; clearly better than the foil-style variant here |
| `cavityFlexibleBottom` | parallel (`4`) | `/tmp/monolithicPetscScreen.kUAclp/cavityFlexibleBottom-beamSplit2d-par` | dimension-aware `beamInCrossFlow` split variant, 2-D nested fluid split (`block_size 3`, fields `0,1` / `2`), `endTime 0.5` | 2 | 22 | 0.07 s | complete | excellent MPI result on this case; still a case-specific winner rather than the best balanced preset |
| `blobInTreacle` | serial | `/tmp/monolithicPetscScreen.kUAclp/blobInTreacle-hybrid05full-serial` | hybrid pressure AMG: full foil-style BoomerAMG block with `strong_threshold 0.5`, `endTime 0.025` | 3 | 46 | 0.19 s | complete | improves iteration count over the foil-style `0.7` preset, but does not beat the best balanced hybrid overall |
| `cavityFlexibleBottom` | serial | `/tmp/monolithicPetscScreen.kUAclp/cavityFlexibleBottom-hybrid05full-serial` | hybrid pressure AMG: full foil-style BoomerAMG block with `strong_threshold 0.5`, `endTime 0.5` | 2 | 21 | 0.13 s | complete | better than the foil-style `0.7` preset here, but not clearly the best cavity choice |
| `beamInCrossFlow` | serial | `/tmp/monolithicPetscScreen.kUAclp/beamInCrossFlow-hybrid05full-serial` | hybrid pressure AMG: full foil-style BoomerAMG block with `strong_threshold 0.5`, `endTime 0.05` | 2 | 23 | 1.78 s | complete | clear beam improvement over the foil-style `0.7` preset, but slightly weaker than the `0.6` full-AMG hybrid |
| `blobInTreacle` | serial | `/tmp/monolithicPetscScreen.kUAclp/blobInTreacle-hybrid07min-serial` | hybrid pressure AMG: simpler beam-style BoomerAMG block with `strong_threshold 0.7`, `endTime 0.025` | 3 | 47 | 0.2 s | complete | strong iteration count, but slower than the `0.6` full-AMG hybrid |
| `cavityFlexibleBottom` | serial | `/tmp/monolithicPetscScreen.kUAclp/cavityFlexibleBottom-hybrid07min-serial` | hybrid pressure AMG: simpler beam-style BoomerAMG block with `strong_threshold 0.7`, `endTime 0.5` | 2 | 12 | 0.12 s | complete | essentially matches the cavity-focused beam-style behaviour, but the beam runtime penalty kept it from advancing |
| `beamInCrossFlow` | serial | `/tmp/monolithicPetscScreen.kUAclp/beamInCrossFlow-hybrid07min-serial` | hybrid pressure AMG: simpler beam-style BoomerAMG block with `strong_threshold 0.7`, `endTime 0.05` | 2 | 21 | 2.06 s | complete | excellent iteration count, but slower than both the checked-in beam reference and the `0.6` full-AMG hybrid |
| `blobInTreacle` | serial | `/tmp/monolithicPetscScreen.kUAclp/blobInTreacle-hybrid06full-serial` | hybrid pressure AMG: full foil-style BoomerAMG block with `strong_threshold 0.6`, `endTime 0.025` | 3 | 42 | 0.18 s | complete | best blob iteration count so far in this split family; runtime remains close enough to advance |
| `blobInTreacle` | parallel (`2`) | `/tmp/monolithicPetscScreen.kUAclp/blobInTreacle-hybrid06full-par` | hybrid pressure AMG: full foil-style BoomerAMG block with `strong_threshold 0.6`, `endTime 0.025` | 3 | 64 | 0.14 s | complete | MPI stays clean and improves on the foil-style `0.7` split |
| `cavityFlexibleBottom` | serial | `/tmp/monolithicPetscScreen.kUAclp/cavityFlexibleBottom-hybrid06full-serial` | hybrid pressure AMG: full foil-style BoomerAMG block with `strong_threshold 0.6`, `endTime 0.5` | 2 | 20 | 0.12 s | complete | strong middle ground between the foil-style and cavity-specific beam-style variants |
| `cavityFlexibleBottom` | parallel (`4`) | `/tmp/monolithicPetscScreen.kUAclp/cavityFlexibleBottom-hybrid06full-par` | hybrid pressure AMG: full foil-style BoomerAMG block with `strong_threshold 0.6`, `endTime 0.5` | 2 | 23 | 0.08 s | complete | nearly matches the best cavity MPI result while staying better balanced elsewhere |
| `beamInCrossFlow` | serial | `/tmp/monolithicPetscScreen.kUAclp/beamInCrossFlow-hybrid06full-serial` | hybrid pressure AMG: full foil-style BoomerAMG block with `strong_threshold 0.6`, `endTime 0.05` | 2 | 22 | 1.64 s | complete | best local beam balance so far: fewer iterations than the checked-in beam split and lower wall time than the foil-style `0.7` split |
| `beamInCrossFlow` | parallel (`2`) | `/tmp/monolithicPetscScreen.kUAclp/beamInCrossFlow-hybrid06full-par` | hybrid pressure AMG: full foil-style BoomerAMG block with `strong_threshold 0.6`, `endTime 0.05` | 2 | 29 | 0.98 s | complete | strongest beam MPI result so far among the general split-family candidates |

## Slow Acceptance

| Case | Mode | Scratch case | Settings | Newton | Linear iters | Exec time | Status | Notes |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| `foilInWind` | serial | `/tmp/monolithicPetscScreen.kUAclp/foilInWind-currentSplit` | current tuned foil split preset, `endTime 0.01` | 5 | 378 | 94.84 s | complete | one-step fully coupled foil acceptance run completed; this is expensive but viable as an acceptance case |
| `foilInWind` | serial | `/tmp/monolithicPetscScreen.kUAclp/foilInWind-beamSplit` | simpler beam-style pressure-AMG variant on the foil split family, `endTime 0.01` | 5 | 308 | 117.64 s | complete | fewer linear iterations than the foil-style preset, but slower overall wall time on the foil acceptance case |
| `foilInWind` | serial | `/tmp/monolithicPetscScreen.kUAclp/foilInWind-hybrid06full` | hybrid pressure AMG: full foil-style BoomerAMG block with `strong_threshold 0.6`, `endTime 0.01` | 5 | 378 | 94.7 s | complete | matches the current foil-style acceptance cost while materially improving the faster screening cases |
| `foilInWind` | serial | `/tmp/monolithicPetscScreen.kUAclp/foilInWind-monoLU` | plain monolithic `bjacobi+lu` baseline from the smaller cases, `endTime 0.01` |  |  |  | failed | PETSc ran out of memory during the first linear solve (`PCSetUp_LU` / error code 55); this rules out plain monolithic `bjacobi+lu` as a universal foil-scale preset |

## Tightened Validation

| Case | Reference | Candidate | Relative L2 (`fluid/U`, `fluid/p`, `solid/U`, `solid/sigmaEq`) | Scalar check | Outcome |
| --- | --- | --- | --- | --- | --- |
| `blobInTreacle` | `/private/tmp/monolithicPetscValidate.kUAclp/blobInTreacle-ref` | `/private/tmp/monolithicPetscValidate.kUAclp/blobInTreacle-hybrid06full` | `1.90e-06`, `1.31e-06`, `1.15e-05`, `2.26e-05` | final total force stayed at `5.80851e-03, -1.91318e-02, ~0`; `Max sigmaEq` `8.10239e-03 -> 8.10229e-03` | passes cleanly against the tighter monolithic reference |
| `cavityFlexibleBottom` | `/private/tmp/monolithicPetscValidate.kUAclp/cavityFlexibleBottom-ref` | `/private/tmp/monolithicPetscValidate.kUAclp/cavityFlexibleBottom-hybrid06full` | `4.42e-07`, `4.69e-07`, `4.20e-05`, `8.57e-05` | final total force stayed at `3.70699e-05, -1.94747e-02, ~0`; `Max sigmaEq` `2.06917e-03 -> 2.06902e-03` | passes cleanly against the tighter monolithic reference |
| `beamInCrossFlow` | `/private/tmp/monolithicPetscValidate.kUAclp/beamInCrossFlow-ref` | `/private/tmp/monolithicPetscValidate.kUAclp/beamInCrossFlow-hybrid06full` | `5.87e-07`, `5.15e-07`, `5.08e-06`, `5.14e-06` | final total force stayed at `1.29019e-01, -2.55603e-01, 2.59067e-01`; `Max sigmaEq` stayed `4.42764` | passes cleanly against the tighter monolithic reference |
| `foilInWind` | `/private/tmp/monolithicPetscScreen.kUAclp/foilInWind-currentSplit` | `/private/tmp/monolithicPetscScreen.kUAclp/foilInWind-hybrid06full` | `1.46e-06`, `4.74e-06`, `4.38e-03`, `1.37e-04` | final total force stayed at `7.72069e-05, ~-1e-10, ~1e-11`; `Max sigmaEq` `4.99774 -> 4.99808` | acceptance remains effectively unchanged while the fast cases improve |

## Notes

- Use clean scratch copies only.
- Use `coupled yes;` and `couplingStartTime 0;`.
- For initial local screening, shorten to one coupled step where possible.
- On this machine, MPI runs must be launched outside the sandbox; in-sandbox
  PRTE socket binds fail before solver startup, and the issue reproduces even
  with `mpirun -np 2 hostname`.
- The nested fluid `U/p` split must be dimension-aware:
  use `block_size 3` with fields `0,1` / `2` on 2-D cases and `block_size 4`
  with fields `0,1,2` / `3` on 3-D cases.
- On the short screening cases, the PETSc log reports cumulative linear
  iterations inside the SNES solve, so the table uses the last reported value.
  On `foilInWind`, the reported linear iterations are per Newton step, so the
  table uses the summed total across the five coupled Newton solves.
- The current working recommendation after local validation is the full
  foil-style `fluid/solid` split with the dimension-aware nested fluid
  `U/p` Schur block and pressure AMG `strong_threshold 0.6`.
- The checked-in monolithic tutorial defaults have been updated to use this
  preset in `blobInTreacle`, `cavityFlexibleBottom`, `beamInCrossFlow`,
  `foilInWind`, and `membraneRoof`. `membraneRoof` was aligned as the 3-D
  analogue of the validated beam/foil setup and still needs its own dedicated
  run.
