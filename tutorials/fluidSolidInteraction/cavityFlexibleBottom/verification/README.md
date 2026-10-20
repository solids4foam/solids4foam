# cavityFlexibleBottom verification study

This opt-in study is a mesh study. It uniformly refines the tutorial fluid and
solid meshes, runs each level to a steady FSI response, and compares the steady
vertical displacement at `(4 -1 0.5)` and the steady vertical interface force
with the mesh study of Tuković et al. (2018).

It is deliberately separate from `regressionTest.sh`: the regression test
checks that the tutorial remains numerically stable, whereas this study checks
convergence towards published steady values. Nothing here is run by
`tutorials/Alltest` or `tutorials/Alltest-regression`.

## Running

Source an OpenFOAM environment with solids4foam built, and run:

```bash
cd tutorials/fluidSolidInteraction/cavityFlexibleBottom/verification
./Allverify
```

Useful options:

```bash
./Allverify --quick              # levels 1 and 2 only, smoke test
./Allverify --variants iqnils    # IQN-ILS instead of Aitken coupling
./Allverify --variants monolithic # quasi-monolithic Newton (OpenFOAM.com)
./Allverify --end-time 400       # longer pseudo-transient
./Allverify --delta-t 0.1        # smaller time step, needed by level 4
./Allverify --reuse              # resume a sweep without re-running cases
```

Each level is a complete copy of the tutorial under the ignored
`verification/work/` directory, so the tutorial itself and its regression test
are never modified. Results are written to the ignored
`verification/postProcessing/` directory as `mesh_convergence.csv`,
`verification_summary.md`, and, when `gnuplot` is available,
`cavityFlexibleBottom_meshConvergence.pdf`.

The study is serial, like the campaign it is migrated from.

## Mesh levels

Each level doubles the in-plane block divisions of both the fluid and the solid
mesh, leaving the single spanwise cell alone. The resulting spacings match the
four published mesh-study points:

| Level | ΔX (m) | Fluid cells | Solid cells |
|---:|---:|---:|---:|
| 1 | 0.1 | 1 688 | 160 |
| 2 | 0.05 | 6 752 | 640 |
| 3 | 0.025 | 27 008 | 2 560 |
| 4 | 0.0125 | 108 032 | 10 240 |

Level 1 is the mesh shipped with the tutorial. Levels 1 to 3 form the default
sweep.

Level 4 is reachable with `--levels 1,2,3,4`, but **not at the tutorial time
step**: partitioned coupling on that mesh diverges at `deltaT = 0.5` with both
the Aitken and the IQN-ILS variant, because the added-mass effect grows as the
interface cells shrink. Run it with a smaller `--delta-t`, and expect a
correspondingly long serial run. The driver reports a diverged level as a
failure rather than reading its blown-up final value as a steady result.

## Steady state

The tutorial is a pseudo-transient route to a steady response, and the
published quantities are steady values. The driver runs each level to a fixed
end time, then checks that the monitored displacement and force have actually
plateaued: the spread of each quantity over the closing tenth of the run must
be within `1e-4` of its final value. A level that has not settled fails rather
than being compared with the published steady value.

The default end time is the tutorial's own `200`. Use `--end-time` when a finer
mesh needs longer.

## Acceptance criteria

- Every level must reach a steady response, as defined above.
- The change in the steady displacement between successive levels must decrease
  monotonically with a positive net order.
- On the finest mesh, the steady displacement must be within 10% and the steady
  interface force within 5% of the published value at the same spacing.

The displacement tolerance is looser than the force tolerance because the two
quantities converge very differently here: the interface force agrees with the
published value to a fraction of a percent on every mesh, whereas the
displacement carries a visible discretisation error on the coarser meshes and
approaches the published curve as the mesh is refined.

The published force values are quoted for an out-of-plane thickness of
`0.05 m`; the tutorial is one metre thick, so the driver compares against
twenty times the published values.

A `--quick` run only exercises the two coarsest meshes, which are not expected
to meet the accuracy tolerances, so it checks only that both levels reach a
steady response.

`verification/Allclean` removes the generated `work/` and `postProcessing/`
directories, as does `./Allclean all` in the case directory; the plain case
`./Allclean` keeps them.

## Reference results

Recorded with OpenFOAM v2512 and the default `aitken` variant on an Apple M1
Ultra. The whole default sweep took about ten minutes in serial.

| Level | ΔX (m) | uy (m) | Error | Fy (N) | Error | Δuy (m) |
|---:|---:|---:|---:|---:|---:|---:|
| 1 | 0.1 | -0.20246 | 13.2% | -5.15312 | 0.13% | – |
| 2 | 0.05 | -0.22883 | 6.88% | -5.23894 | 0.19% | 2.64e-02 |
| 3 | 0.025 | -0.23563 | 5.59% | -5.25632 | 0.14% | 6.79e-03 |

The interface force agrees with the published mesh study to within 0.2% on
every level. The displacement converges cleanly at a net order of 1.96, but
towards a value about 5% away from the published curve rather than onto it.
That residual offset is consistent with the difference already noted in the
tutorial README between the solids4foam prediction and the published
displacement at the coarsest spacing, and it is why the displacement tolerance
is looser than the force tolerance.

## Quasi-monolithic coupling

`--variants monolithic` runs the quasi-monolithic Newton coupling
(`./Allrun monolithic`, OpenFOAM.com only). Its inlet profile is ramped in over
the first 20 s (`inlet_ramp_s` in the reference JSON) with
`transitionalParabolicVelocity`, which gives the same parabola: the monolithic
coupling resolves the pressure impulse of the impulsive start, which inverted
the fluid mesh on level 3, and the ramp does not change the steady response.
Recorded with OpenFOAM v2512 on macOS, serial, at `deltaT = 0.5`, with the
partitioned Aitken coupling on levels 1 to 3 for comparison (level 3 at
`deltaT = 0.25`, at which its coupling converges; mono and part are the
monolithic and partitioned couplings). **The monolithic values in this
table are pre-#597** (before the `newtonIcoFluid` continuity fix); the rerun
is in "Rerun with the newtonIcoFluid continuity fix" below and supersedes them:

| Level | ΔX (m) | Uy mono (m) | Uy part (m) | Fy mono (N) | Fy part (N) |
|---:|---:|---:|---:|---:|---:|
| 1 | 0.1 | -0.201557 | -0.202495 | -5.1225 | -5.1527 |
| 2 | 0.05 | -0.228817 | -0.228799 | -5.2381 | -5.2378 |
| 3 | 0.025 | -0.235676 | -0.235628 | -5.2576 | -5.2564 |
| 4 | 0.0125 | -0.237336 | – | -5.2606 | – |

The two couplings agree to 0.02% in Uy and Fy on level 3. The monolithic Uy
converges at observed orders of 1.99 and 2.05 to an extrapolated
-0.23787 m. The published mesh-study values are about 5% larger in magnitude
(-0.2496 m at ΔX = 0.025 m), as is an independent ANSYS calculation (about
-0.25 m); with both couplings converged to the same value, the difference is
not a discretisation error of either coupling, and the definition of the
published case (monitor point, boundary conditions, material) is being
checked.

## References

Ž. Tuković, A. Karač, P. Cardiff, H. Jasak, A. Ivanković, OpenFOAM finite
volume solver for fluid-solid interaction, *Transactions of FAMENA*, 42(3),
1–31, 2018.

## Partitioned four-level sequence and the monitor point (recorded 17 Oct 2026)

Data: `results/four_level_partitioned_monolithic.csv` and `.json`
(`results/sample_pointD.py` reads the solid `pointD` at a mesh node).

### The tutorial monitor point is on the plate underside

The solid mesh of this tutorial occupies `-0.4 <= y <= 0` (blockMeshDict), with
the wetted interface at `y = 0`; the geometry figure puts point A on the top
(wetted) surface at `(4, 0)`. The `solidPointDisplacement` function object uses
the mesh point nearest to `(4 -1 0.5)`, which is `(4, -0.4)`: the **underside**
of the plate. The `uy` recorded above (and in `mesh_convergence.csv`) is
therefore the underside value, not the wetted-surface value. The values below
are sampled from the final-time `pointD` at the mesh nodes `(4, y)` (all three
`y` are nodes on every level; `z` averaged over the two spanwise nodes). The
underside value reproduces the recorded `solidPointDisplacement` value, and
`postProcess -region solid` probes of `D` with `cellPoint` interpolation agree
with the nodal `pointD` to <= 0.3 % (mesh 1) and <= 2e-5 m (meshes 3, 4).

Quasi-monolithic `uy` (m) (**pre-#597** values; see the rerun section below), `y` in tutorial coordinates (top `y = 0`, mid
`y = -0.2`, underside `y = -0.4`); the mid-thickness point corresponds to
(4, -1.2) in ANSYS coordinates if its plate top is at `y = -1` (assumption,
not confirmed by the ANSYS author):

| ΔX (m) | top | mid | underside (recorded) | ANSYS mid (Tandis) | Tukovic A |
|---:|---:|---:|---:|---:|---:|
| 0.1 | -0.20361 | -0.20504 | -0.20156 | -0.2502 | -0.2333 |
| 0.05 | -0.23143 | -0.23268 | -0.22882 | -0.2535 | -0.2457 |
| 0.025 | -0.23845 | -0.23965 | -0.23568 | -0.2525 (not steady) | -0.2496 |
| 0.0125 | -0.24014 | -0.24134 | -0.23734 | – | -0.2506 |

Orders 1.99 and 2.05; Richardson (p = 2): top -0.2407, mid -0.2419,
underside -0.2379. The top-to-mid difference is 0.5 % (0.0012 m), the
mid-to-underside 1.7 %, so the choice of surface changes `uy` by at most 2 %.
Against ANSYS at the same (mid) point we are lower in magnitude by 18 %, 8.2 %
and 5.1 % on the three meshes; Tukovic's finest value (-0.2506) is 4.2 % above
our top value. The surface-versus-mid choice cannot explain a 4 to 5 % offset,
so the discrepancy with the published and ANSYS values is not resolved by the
monitor-point correction. The ANSYS values were reported to hardly change with
mesh (-0.2502, -0.2535), ours change by 13 % over the same meshes. `ux` at all
three points is below 1.1e-4 m (ANSYS `ux` of +0.0325 m reported at its first
point is not reproduced).

### Partitioned versus quasi-monolithic, four levels

`uy` is the underside (recorded) value; `Fy` is the total interface force
(N per metre of span).

| ΔX (m) | uy partitioned | uy monolithic (pre-#597) | Fy partitioned | Fy monolithic (pre-#597) |
|---:|---:|---:|---:|---:|
| 0.1 | -0.202495 | -0.201557 | -5.15266 | -5.12249 |
| 0.05 | -0.228799 | -0.228817 | -5.23784 | -5.23811 |
| 0.025 | -0.235628 | -0.235676 | -5.25637 | -5.25764 |
| 0.0125 | -0.237307 | -0.237336 | -5.25993 | -5.26056 |
| order | 1.95, 2.02 | 1.99, 2.05 | 2.20, 2.38 | 2.57, 2.74 |
| Richardson (p = 2) | -0.23787 | -0.23789 | -5.2611 | -5.2615 |

The `uy` orders are compatible (about 2) so the p = 2 extrapolation is
meaningful; the `Fy` orders are higher and not constant, so the Fy
extrapolations (using the observed last order: -5.2608 and -5.2611) are only
indicative. Partitioned and monolithic agree to 0.01 to 0.02 % in `uy` and
`Fy` on the three finest meshes.

### Rerun with the newtonIcoFluid continuity fix (solids4foam PR #597, recorded 20 Oct 2026)

`newtonIcoFluid` built its continuity residual from `fvc::div(U)` rather than
from the flux that convects momentum; #597 fixes that. Everything else is as
for the pre-#597 monolithic runs (same case files from the tutorial variants,
inlet ramped over 20 s, `deltaT = 0.5`, SNES and linear tolerances unchanged).
Only the code version changed. Data: `monolithic_597` rows of the CSV/JSON;
the earlier monolithic rows are now labelled `monolithic_pre597`.

Build: solids4foam branch `feature-petsc-snes-quasi-monolithic` at 257b2e564
(contains 597cc6317), OpenFOAM v2512 (OpenCFD) on xenosim (Ubuntu 22.04),
gcc 11.4.0, wmake `-O3` (default `Opt`), PETSc 3.24.1 with MUMPS and HYPRE;
private lib/app directories; the `solids4Foam` solver was linked and run.
Levels 1 to 3 ran serially; level 4 ran on 4 MPI ranks (`Allrun monolithic
parallel`, scotch) because serial took about 55 to 80 s per step, with the
fluid mesh-motion PETSc preconditioner changed from `bjacobi`+LU to
`redundant`+LU (exact LU as in serial: a 4-rank `bjacobi` run failed its
mesh-motion linear solve at t = 15.5 s). The first three levels at the same
settings were not rerun in parallel, so a serial/parallel difference on level 4
is not separately quantified. End time 300 s (not 200 s): at 200 s the level 2
`uy` still oscillated with peak-to-peak 2.3e-5 m (decaying), above the 1e-5 m
target; at 300 s the change in `uy` over the last 10 % is at most 5e-6 m on every
level and `Fy` changes by at most 2e-7 N. Levels 1 and 2 at 200 s gave
-0.201552 and -0.228473, within 6e-6 m of the 300 s values.

`uy` underside (recorded point (4, -0.4)); top and mid are sampled from the
final `pointD` at the nodes (4, 0) and (4, -0.2):

| ΔX (m) | uy pre-#597 | uy mono (#597) | uy partitioned | Fy pre-#597 | Fy mono (#597) | Fy partitioned |
|---:|---:|---:|---:|---:|---:|---:|
| 0.1 | -0.201557 | -0.201546 | -0.202495 | -5.12249 | -5.12258 | -5.15266 |
| 0.05 | -0.228817 | -0.228468 | -0.228799 | -5.23811 | -5.22750 | -5.23784 |
| 0.025 | -0.235676 | -0.235448 | -0.235628 | -5.25764 | -5.25072 | -5.25637 |
| 0.0125 | -0.237336 | -0.237187 | -0.237307 | -5.26056 | -5.25611 | -5.25993 |
| order | 1.99, 2.05 | 1.95, 2.01 | 1.95, 2.02 | 2.57, 2.74 | 2.18, 2.11 | 2.20, 2.38 |
| Richardson (p = 2, levels 3-4) | -0.23789 | -0.23777 | -0.23787 | -5.2615 | -5.2579 | -5.2611 |

Top (y = 0) `uy`: -0.203597, -0.231074, -0.238211, -0.239990; mid (y = -0.2):
-0.205030, -0.232320, -0.239415, -0.241184 (pre-#597 top -0.203607,
-0.231431, -0.238445, -0.240142). `ux` at all three points is below 1.1e-4 m.
The monitor point in the tutorial is now (4 0 0.5), the top surface (PR
#584); all three points are recorded regardless.

Change relative to pre-#597 (`uy` / `Fy`): +0.005 % / -0.002 % (level 1),
+0.15 % / +0.20 % (level 2), +0.10 % / +0.13 % (level 3), +0.06 % / +0.08 %
(level 4), where positive means a smaller magnitude. The fix lowers the
magnitude of both `uy` and `Fy`, by up to 0.15 to 0.2 % at ΔX = 0.05 m. The
monolithic values are now below the partitioned ones by 0.14 % (level 2),
0.08 % (level 3) and 0.05 % (level 4) in `uy` (0.20, 0.11 and 0.07 % in `Fy`),
where pre-#597 they agreed to 0.01 to 0.02 % on those levels (level 1: 0.46 %
in `uy` and 0.58 % in `Fy` before and after; the level-1 offset is the
partitioned one). The p = 2 extrapolated `uy` values still agree to 0.04 %
(-0.23777 monolithic against -0.23787 partitioned); the `Fy` extrapolations
differ by 0.06 % (-5.2579 against -5.2611, with `Fy` orders around 2.1 to 2.4
so these are indicative only). The monolithic-minus-partitioned `uy` difference
shrinks with refinement (3.3e-4, 1.8e-4, 1.2e-4 m on levels 2 to 4, a ratio of
only about 1.5 to 1.8 per halving, so not at second order); its origin is not
established. The
quasi-monolithic solver shares most of its code with the partitioned one, so
this agreement is a consistency check between formulations, not an independent
reference. The offset to the published values (about 5 %) is unchanged.

### Partitioned level 4 (ΔX = 0.0125 m)

Platform and build: macOS arm64 (Mac Studio), Apple Clang, OpenFOAM v2512
(build 87ed40d256-20251219), solids4foam v2.4 development library built from
branch `feature-ring-monolithic-verification` at 0698d55cd (HEAD when the
`ring3` library used by the recorded partitioned runs was built on 15 Oct
2026; the working tree at build time is not independently verified), Aitken
coupling, `relaxationFactor 0.1`, `outerCorrTolerance 1e-5`, coupling from
t = 10 s, end time 200 s. Level 4 ran on 4 MPI ranks (scotch).

* At the tutorial pseudo-time step (0.5 s) level 4 diverged earlier; the
  fluid-only start-up (the coupling starts at 10 s) also diverges at
  ΔT = 0.125 s (mean Courant number 8, blow-up at t = 9) and with a 20 s
  ramped inlet at ΔT = 0.125 and 0.25 s. ΔT = 0.0625 s (mean Courant number 4)
  is stable.
* With the default FSI test (`min(r1, r2)`) at ΔT = 0.0625 s the run reaches
  t = 200 s and `uy` = -0.237312, but steps are accepted after one iteration
  with `r1` of order 1 (median 0.76; issue #552) and `Fy` fluctuates by 0.1 N
  from step to step (-5.20 N at the end): the coupling is not converged and
  that run is not used.
* `requireAllResidualMeasures yes` with the default `nOuterCorr 20` fails at
  t = 56.75 s (residual floor about 1.6e-4 from the 1e-6 fluid linear solver
  tolerance); with `nOuterCorr 40` it fails at t = 142.6 s (2.1e-5). The
  accepted run uses `requireAllResidualMeasures yes`, `nOuterCorr 40` and
  the fluid linear solver tolerance 1e-8: `r1` at acceptance in the last 10 %
  of the run has a maximum of 6.0e-6 (below the 1e-5 tolerance, 3.8 outer
  iterations on average), r2 below 5e-13.
* Steadiness: the change in `uy` over the last 10 % of the run is 1.3e-5 m
  (target 1e-5; the monolithic runs are at 2.3e-5 m, oscillating with a
  decaying amplitude, extrapolated remaining about 5e-6 m).

The time-step reduction does not change the steady state on the coarser
meshes (strict test, smaller ΔT, compared with the recorded runs): level 2
ΔT = 0.25 s: `uy` -0.228830 against -0.228799 (0.014 %), `Fy` -5.23881 against
-5.23784 (0.02 %); level 3 ΔT = 0.125 s: -0.235636 against -0.235628 (0.003 %),
-5.25661 against -5.25637 (0.005 %); level 1 at ΔT = 0.5 s with the strict
test: -0.202499 against -0.202495. The recorded levels 1 to 3 used the default
FSI test with `r1` at acceptance between 1e-4 and 6e-2 in the closing 10 % of
the run (50 to 92 % of steps above the tolerance); the strict reruns
(`partitioned_strict` rows) have `r1` below 1e-6, and agree with them to
0.02 % or better. The `partitioned_strict` rows are the four-level set with
the acceptance check passed; the `partitioned_legacy_test` rows are the
recorded runs (level 4 shown for the default test, not usable).
