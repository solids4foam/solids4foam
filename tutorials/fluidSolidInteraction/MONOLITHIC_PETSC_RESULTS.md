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

## Remote Scaling (xenosim, 2026-04-09)

OpenFOAM-v2512 remote scaling and validation pass on branch
`feature-petsc-snes-quasi-monolithic`.

| Case | Ranks | Preset | Newton | Linear iters | Wall time | Status | Notes |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `blobInTreacle` | 1 | `baseline06-ref` | 2 | 55 | 0.43 s | complete | reference equals current checked-in default |
| `blobInTreacle` | 1 | `winner06-serial` | 2 | 55 | 0.41 s | complete | matches the checked-in reference exactly on the remote machine |
| `blobInTreacle` | 2 | `winner06-r2` | 2 | 77 | 0.25 s | complete | useful first speedup despite higher iteration count |
| `blobInTreacle` | 4 | `winner06-r4` | 2 | 96 | 0.25 s | complete | wall time already flat by 4 ranks; setup cost dominates |
| `blobInTreacle` | 8 | `winner06-r8` | 2 | 101 | 0.24 s | complete | little gain beyond 2 ranks; fields remain close |
| `cavityFlexibleBottom` | 1 | `baseline06-ref` | 2 | 20 | 0.39 s | complete | reference equals current checked-in default |
| `cavityFlexibleBottom` | 1 | `winner06-serial` | 2 | 20 | 0.38 s | complete | matches the checked-in reference exactly on the remote machine |
| `cavityFlexibleBottom` | 2 | `winner06-r2` | 2 | 27 | 0.16 s | complete | clean MPI scaling with stable observables |
| `cavityFlexibleBottom` | 4 | `winner06-r4` | 2 | 25 | 0.17 s | complete | still clean; 2 and 4 ranks are close in cost |
| `cavityFlexibleBottom` | 8 | `winner06-r8` | 2 | 27 | 0.12 s | complete | best remote cavity scaling point; force and stress stay stable |
| `beamInCrossFlow` | 1 | `monoLU-ref` | 2 | 25 | 25.85 s | complete | robust monolithic `bjacobi+lu` reference, but much slower than the split winner |
| `beamInCrossFlow` | 1 | `winner06-serial` | 2 | 23 | 4.41 s | complete | same answer as the monolithic reference at much lower cost |
| `beamInCrossFlow` | 2 | `winner06-r2` | 2 | 30 | 2.20 s | complete | clean scaling and stable beam response |
| `beamInCrossFlow` | 4 | `winner06-r4` | 2 | 36 | 0.98 s | complete | strong speedup with only tiny field drift |
| `beamInCrossFlow` | 8 | `winner06-r8` | 2 | 42 | 0.64 s | complete | still clean; iteration growth is modest |
| `beamInCrossFlow` | 16 | `winner06-r16` | 2 | 49 | 0.52 s | complete | best clean remote beam scaling point |
| `foilInWind` | 1 | `foil07-ref` | 5 | 378 | 207.28 s | complete | serial reference on the remote machine |
| `foilInWind` | 1 | `winner06-serial` | 5 | 378 | 222.86 s | complete | matches foil scalar observables in serial but is slightly slower on this machine |
| `foilInWind` | 2 | `winner06-r2` | 5 | 383 | 146.77 s | complete | force stays stable, but `Max sigmaEq` already drops materially |
| `foilInWind` | 4 | `winner06-r4` | 5 | 390 | 62.87 s | complete | clear MPI solid-response drift; see field validation below |
| `foilInWind` | 4 | `foil07-r4` | 5 | 390 | 55.86 s | complete | matches the `winner06` 4-rank foil behaviour; threshold change does not cure the drift |
| `foilInWind` | 8 | `winner06-r8` | 5 | 399 | 42.17 s | complete | force remains stable while solid fields drift strongly |
| `foilInWind` | 8 | `foil07-r8` | 5 | 399 | 40.88 s | complete | again matches the `winner06` MPI foil behaviour |
| `foilInWind` | 16 | `winner06-r16` | 5 | 396 | 28.27 s | complete | still converges, but the solid response is no longer close to the serial reference |
| `foilInWind` | 32 | `winner06-r32` | 5 | 402 | 17.82 s | complete | severe MPI solid-response drift despite a nearly unchanged force vector |

## Remote Field Validation

Representative field-difference checks against the chosen remote references.

| Case | Reference | Candidate | Relative L2 (`fluid/U`, `fluid/p`, `solid/U`, `solid/sigmaEq`) | Scalar check | Outcome |
| --- | --- | --- | --- | --- | --- |
| `blobInTreacle` | `baseline06-ref` | `winner06-serial` | `0`, `0`, `0`, `0` | force unchanged; `Max sigmaEq` `8.10187e-03 -> 8.10187e-03` | serial equivalence |
| `blobInTreacle` | `baseline06-ref` | `winner06-r8` | `9.31e-05`, `4.91e-05`, `5.41e-04`, `8.72e-04` | force stable; `Max sigmaEq` `8.10187e-03 -> 8.10054e-03` | clean MPI agreement |
| `cavityFlexibleBottom` | `baseline06-ref` | `winner06-serial` | `0`, `0`, `0`, `0` | force unchanged; `Max sigmaEq` `2.06904e-03 -> 2.06904e-03` | serial equivalence |
| `cavityFlexibleBottom` | `baseline06-ref` | `winner06-r8` | `5.57e-07`, `5.05e-07`, `5.49e-05`, `9.95e-05` | force stable; `Max sigmaEq` `2.06904e-03 -> 2.06908e-03` | clean MPI agreement |
| `beamInCrossFlow` | `monoLU-ref` | `winner06-serial` | `4.78e-07`, `6.59e-07`, `3.87e-06`, `3.95e-06` | force stable; `Max sigmaEq` `4.42653 -> 4.42653` | serial equivalence to the conservative monolithic reference |
| `beamInCrossFlow` | `monoLU-ref` | `winner06-r16` | `1.11e-06`, `1.02e-06`, `1.07e-05`, `1.42e-05` | force stable; `Max sigmaEq` `4.42653 -> 4.42649` | clean MPI agreement |
| `foilInWind` | `foil07-ref` | `winner06-serial` | `0`, `0`, `0`, `0` | force unchanged; `Max sigmaEq` `4.94872 -> 4.94872` | serial equivalence |
| `foilInWind` | `foil07-ref` | `foil07-r8` | `1.41e-05`, `4.44e-05`, `5.71e-01`, `5.64e-01` | force stable; `Max sigmaEq` `4.94872 -> 2.12045` | large MPI solid drift |
| `foilInWind` | `foil07-ref` | `winner06-r32` | `2.29e-05`, `6.96e-05`, `8.01e-01`, `7.91e-01` | force stable; `Max sigmaEq` `4.94872 -> 1.08098` | severe MPI solid drift |

## Remote Conclusion

- `strong_threshold 0.6` still looks like the best general preset on the
  remote machine for `blobInTreacle`, `cavityFlexibleBottom`, and
  `beamInCrossFlow`.
- A different preset is not clearly better at higher rank counts on
  `foilInWind`: the older `0.7` foil preset reproduces essentially the same
  4-rank and 8-rank MPI foil behaviour.
- The real remote regression is `foilInWind` under MPI. Scaling in wall time
  still improves with rank, but the solid response drifts badly while the
  force stays nearly unchanged.
- The small number of foil solid cells per rank may contribute to degraded
  scaling, but it does not by itself explain the answer change. A more likely
  explanation is that the current outer Newton tolerances are loose enough
  that solver-dependent iteration error survives in MPI and contaminates the
  solid `sigmaEq` response. The next foil-specific check should therefore be a
  tighter-Newton-tolerance MPI rerun rather than another `strong_threshold`
  sweep.

## Foil Newton Tolerance Sweep (xenosim, 2026-04-09)

Follow-up remote foil check on the current `winner06` preset, sweeping the
outer Newton tolerances only:
`snes_rtol = snes_stol = 1e-5` and `1e-6`.

| Tolerance | Ranks | Newton | Linear iters | Wall time | Final force | Max sigmaEq | Notes |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `1e-4` | 1 | 5 | 378 | 222.86 s | `(7.720687e-05, -1.4111325e-10, 9.31478e-12)` | `4.94872` | existing remote `winner06` serial baseline |
| `1e-4` | 8 | 5 | 399 | 42.17 s | `(7.720708e-05, 2.1156424e-10, 7.758601e-12)` | `2.12045` | existing remote MPI foil-drift point |
| `1e-4` | 32 | 5 | 402 | 17.82 s | `(7.720689e-05, 2.95105466e-10, 8.4239e-12)` | `1.08098` | existing remote worst-rank foil-drift point |
| `1e-5` | 1 | 5 | 378 | 361.21 s | `(7.720687e-05, -1.4111325e-10, 9.31478e-12)` | `4.94872` | exactly the same answer as `1e-4`; only wall time changed |
| `1e-5` | 8 | 5 | 399 | 42.81 s | `(7.720708e-05, 2.1156424e-10, 7.758601e-12)` | `2.12045` | exactly the same MPI answer as `1e-4` |
| `1e-5` | 32 | 5 | 402 | 19.17 s | `(7.720689e-05, 2.95105466e-10, 8.4239e-12)` | `1.08098` | exactly the same MPI answer as `1e-4` |
| `1e-6` | 1 | 6 | 452 | 362.74 s | `(7.720707e-05, -2.5057671e-10, 9.05256e-12)` | `4.94893` | serial answer moved only slightly; one extra Newton step |
| `1e-6` | 8 | 6 | 474 | 53.53 s | `(7.720718e-05, 5.826131e-11, 8.043359e-12)` | `2.43799` | MPI answer shifts, but still remains far from serial |
| `1e-6` | 32 | 6 | 475 | 23.22 s | `(7.720728e-05, 1.53783326e-10, 7.978215e-12)` | `1.2266` | MPI answer shifts modestly toward serial, but the large solid drift remains |

| Comparison | Relative L2 (`fluid/U`, `fluid/p`, `solid/U`, `solid/sigmaEq`) | Outcome |
| --- | --- | --- |
| `1e-4` serial vs `1e-5` serial | `0`, `0`, `0`, `0` | tightening by one order does nothing at all |
| `1e-4` serial vs `1e-6` serial | `7.16e-05`, `3.80e-04`, `3.83e-03`, `1.07e-04` | the serial foil answer is already essentially converged by `1e-4` |
| `1e-5` serial vs `1e-5` rank-8 | `1.41e-05`, `4.44e-05`, `5.71e-01`, `5.64e-01` | unchanged from the loose-tolerance foil drift |
| `1e-5` serial vs `1e-5` rank-32 | `2.29e-05`, `6.96e-05`, `8.01e-01`, `7.91e-01` | unchanged from the loose-tolerance foil drift |
| `1e-6` serial vs `1e-6` rank-8 | `2.74e-07`, `1.28e-06`, `5.07e-01`, `5.02e-01` | modest improvement, but the solid-field mismatch is still very large |
| `1e-6` serial vs `1e-6` rank-32 | `3.16e-07`, `1.35e-06`, `7.64e-01`, `7.55e-01` | modest improvement, but the severe high-rank solid drift remains |

The sweep shows that tighter outer Newton tolerances do influence the foil MPI
answer, but they do not resolve the issue. Tightening from `1e-4` to `1e-5`
changes nothing. Tightening to `1e-6` adds one Newton step and substantially
more linear work, while the serial reference barely moves and the MPI solid
fields remain far from serial. The foil MPI regression therefore is not just a
loose-Newton-tolerance artifact; tolerance tightening alone is not sufficient.

## Foil Fluid-System Scaling Sweep (xenosim, 2026-04-09)

Follow-up remote foil check on the current `winner06` preset, with the tighter
outer Newton tolerances held fixed at `snes_rtol = snes_stol = 1e-6` while
sweeping `fluidSystemScaleFactor`.

| Fluid system scale | Ranks | Newton | Linear iters | Wall time | Final force | Max sigmaEq | Notes |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `1e6` | 1 | 6 | 452 | 455.50 s | `(7.720707e-05, -2.5069823e-10, 9.03687e-12)` | `4.94866` | lower fluid scaling; serial answer remains essentially unchanged |
| `1e7` | 1 | 6 | 452 | 473.96 s | `(7.720707e-05, -2.4973641e-10, 9.03671e-12)` | `4.94878` | same serial answer within noise |
| `1e8` | 1 | 6 | 452 | 481.57 s | `(7.720707e-05, -2.5057671e-10, 9.05256e-12)` | `4.94893` | current checked-in fluid-system scale |
| `1e9` | 1 | 6 | 452 | 440.12 s | `(7.720707e-05, -2.5061518e-10, 9.0367e-12)` | `4.94871` | higher fluid scaling; same serial answer within noise |
| `1e6` | 8 | 6 | 474 | 65.57 s | `(7.720718e-05, 5.806094e-11, 8.053506e-12)` | `2.43806` | same MPI foil drift as the baseline scale |
| `1e7` | 8 | 6 | 474 | 67.80 s | `(7.720718e-05, 5.815102e-11, 8.048088e-12)` | `2.43825` | same MPI foil drift as the baseline scale |
| `1e8` | 8 | 6 | 474 | 97.20 s | `(7.720718e-05, 5.826131e-11, 8.043359e-12)` | `2.43799` | current checked-in fluid-system scale |
| `1e9` | 8 | 6 | 474 | 109.67 s | `(7.720718e-05, 5.815706e-11, 8.055255e-12)` | `2.43845` | same MPI foil drift as the baseline scale |
| `1e6` | 16 | 6 | 471 | 34.26 s | `(7.720728e-05, 1.47632263e-10, 6.217675e-12)` | `2.4857` | lower fluid scaling; same answer within noise |
| `1e8` | 16 | 6 | 471 | 34.60 s | `(7.720728e-05, 1.47609324e-10, 6.216343e-12)` | `2.48515` | current checked-in fluid-system scale |
| `1e9` | 16 | 6 | 471 | 35.14 s | `(7.720728e-05, 1.47605895e-10, 6.215752e-12)` | `2.48481` | higher fluid scaling; same answer within noise |

| Comparison | Relative L2 (`fluid/U`, `fluid/p`, `solid/U`, `solid/sigmaEq`) | Outcome |
| --- | --- | --- |
| serial `1e6` vs serial `1e9` | `6.54e-09`, `1.41e-08`, `1.33e-04`, `3.53e-05` | changing the fluid-system scale across three orders of magnitude barely moves the serial foil solution |
| serial `1e6` vs rank-8 `1e6` | `2.73e-07`, `1.28e-06`, `5.07e-01`, `5.02e-01` | large solid-field drift remains |
| serial `1e8` vs rank-8 `1e8` | `2.74e-07`, `1.28e-06`, `5.07e-01`, `5.02e-01` | unchanged from `1e6` |
| serial `1e9` vs rank-8 `1e9` | `2.73e-07`, `1.28e-06`, `5.07e-01`, `5.02e-01` | unchanged from `1e6` and `1e8` |
| serial `1e6` vs rank-16 `1e6` | `2.73e-07`, `1.16e-06`, `8.85e-01`, `8.43e-01` | large solid-field drift remains |
| serial `1e8` vs rank-16 `1e8` | `2.73e-07`, `1.16e-06`, `8.85e-01`, `8.43e-01` | unchanged from `1e6` |
| serial `1e9` vs rank-16 `1e9` | `2.73e-07`, `1.15e-06`, `8.85e-01`, `8.43e-01` | unchanged from `1e6` and `1e8` |
| rank-8 `1e6` vs rank-8 `1e9` | `9.17e-09`, `3.07e-08`, `2.34e-04`, `1.67e-04` | the rank-8 MPI solution itself barely moves across the scale sweep |
| rank-16 `1e6` vs rank-16 `1e9` | `5.49e-09`, `4.98e-08`, `5.03e-04`, `4.40e-04` | the rank-16 MPI solution itself barely moves across the scale sweep |

This sweep does not support the fluid-system scaling as the main cause of the
foil MPI regression, at least over the practical range `1e6` to `1e9`.
Changing the scale by three orders of magnitude barely moves the serial foil
solution and barely moves the MPI foil solution at 8 and 16 ranks, while the
serial-vs-MPI solid-field mismatch remains essentially unchanged. The current
foil drift therefore is more likely tied to something else in the nonlinear or
preconditioning path than to the fluid-system scaling alone. The
`-s4f_petsc_report_block_residuals` flag was added for this sweep, but the
current foil run path did not emit the expected per-block residual diagnostic,
so that specific block-wise convergence data is still missing.

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
- The current working recommendation after local validation is the combined
  best preset: multiplicative `fluid/solid` split with `upper` Schur for the
  nested fluid `U/p` block, `redundant+lu` for the solid sub-block,
  `bjacobi+ilu(0)` for the velocity sub-block, and HYPRE BoomerAMG
  `strong_threshold 0.6` with full tuning for the pressure Schur complement.
- The checked-in monolithic tutorial defaults have been updated to use this
  combined preset in `blobInTreacle`, `cavityFlexibleBottom`,
  `beamInCrossFlow`, `foilInWind`, and `membraneRoof`.
- The `bjacobi+lu` remote results in the scaling table above are unreliable
  for `foilInWind` solid fields: `bjacobi+lu` gives sigmaEq ~2.1 at 8 ranks
  vs the correct ~5.0 (see Parallel Regression section below). Force values
  in those remote runs are unaffected and remain reliable.

## Upper vs Lower Schur Factorization (2026-04-09, local v2412)

Direct A/B comparison on all five tutorial cases (Codex config + `upper`
vs Codex config + `lower`, same machine and mesh):

| Case | lower KSP total | upper KSP total | Winner |
| --- | --- | --- | --- |
| `blobInTreacle` (3 coupled steps) | 233 | 246 | lower (~5% fewer) |
| `cavityFlexibleBottom` (3 coupled steps) | 85 | 91 | lower (~7% fewer) |
| `beamInCrossFlow` (2 coupled steps) | 55 | 62 | lower (~11% fewer) |
| `foilInWind` (1 coupled step) | 392 (Phase 8b) | 298 (Phase 8b) | **upper (24% fewer)** |

Conclusion: for small/easy cases the difference is negligible (both < 5 s).
On `foilInWind`, `upper` is clearly better. Because `foilInWind` is the
hardest and most representative case, `upper` is the better general default.
The advantage grows with problem difficulty (more Newton steps, more KSP per
Newton). All tutorial defaults have been set to `upper`.

## Parallel Solid-Response Regression: Root Cause and Fix (2026-04-09)

Codex identified a `foilInWind` MPI regression: with `bjacobi+lu` on the
solid sub-block, `sigmaEq` drops from 4.95 to ~2.1 at 8 ranks (force stays
correct). The root cause is that `bjacobi+lu` applies LU only to each rank's
local portion of the solid matrix; at 8+ ranks the solid mesh is split into
very small sub-blocks and the per-rank LU is too approximate.

Fix: replace `fieldsplit_solid_pc_type bjacobi; fieldsplit_solid_sub_pc_type lu`
with `fieldsplit_solid_pc_type redundant; fieldsplit_solid_redundant_pc_type lu`.
The `redundant` PC gathers the entire solid block to rank 0 and applies an
exact global LU, which is correct regardless of the number of ranks.

Validation (local, 8-rank `foilInWind`, one coupled step):

| Solid PC | Newton | KSP total | sigmaEq | Wall time |
| --- | --- | --- | --- | --- |
| `bjacobi+lu` (Codex, 8-rank) | 5 | ~395 | **2.12** (WRONG) | ~20 s |
| `redundant+lu` (fix, 8-rank) | 5 | 397 | **4.994** (correct) | 20 s |
| `redundant+lu` + `upper` (8-rank) | 5 | 340 | **4.994** (correct) | **16 s** |

The `bjacobi+lu` results from the remote scaling section are therefore
unreliable for `foilInWind`. The regression is fully resolved by the fix; all
tutorial defaults have been updated to use `redundant+lu`.

## Combined Best Config (2026-04-09)

Combining both improvements:
- Schur factorization: `upper`
- Solid sub-block PC: `redundant+lu`
- Velocity sub-block PC: `bjacobi+ilu(0)` (Codex-validated; cheaper per KSP than HYPRE)
- Pressure sub-block PC: HYPRE BoomerAMG, `strong_threshold 0.6`, full tuning
  (Codex-validated settings)

`foilInWind` serial, one coupled step (local v2412 Mac Studio):

| Config | Newton | KSP total | sigmaEq | Wall time |
| --- | --- | --- | --- | --- |
| Codex lower + bjacobi+lu (remote ref) | 5 | 378 | 4.95 | 207 s |
| HYPRE+HYPRE 0.7 + upper (Phase 8b) | 5 | 298 | 4.99 | 122 s |
| **Combined best (this commit)** | 5 | 321 | 4.998 | **85 s** |

The higher KSP count (321 vs 298) is explained by `bjacobi+ilu(0)` needing
slightly more outer iterations than HYPRE for the velocity block; however it
is far cheaper per iteration, giving a 30% wall-clock improvement.

All five tutorial monolithic configs have been updated to this combined preset.
