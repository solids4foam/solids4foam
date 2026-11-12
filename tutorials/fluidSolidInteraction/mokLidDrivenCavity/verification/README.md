# mokLidDrivenCavity verification study

This opt-in study checks the `mokLidDrivenCavity` tutorial against the
published midpoint displacement histories of the lid-driven cavity with a
flexible bottom. It is separate from `regressionTest.sh`: the regression test
checks that the tutorial is numerically stable, whereas this study checks mesh
and time-step convergence, compares the high-order and the standard solid
discretisations, and compares the converged histories with the literature.
Nothing here is run by `tutorials/Alltest` or `tutorials/Alltest-regression`.

## Running

Source an OpenFOAM environment with solids4foam built with PETSc, and run:

```bash
cd tutorials/fluidSolidInteraction/mokLidDrivenCavity/verification
./Allverify                  # all three studies
./Allverify --study mesh     # one study: mesh, timestep or solid
./Allverify --quick          # smoke test: two members per study, to t = 10 s
./Allverify --reuse          # reuse members completed with the same settings
./Allverify --cores 4        # run every member on 4 MPI ranks
```

Each member is a complete copy of the tutorial under the ignored
`verification/work/` directory, so the tutorial and its regression test are
never modified. Results go to the ignored `verification/postProcessing/`
directory: one `<study>_study.csv` table and one `<study>_histories.csv` file
per study, `verification_summary.md`, and, when `gnuplot` is available, a
`<study>_histories.png` plot against the published curves. `Allverify` returns
zero only when every check passes. Each completed member records its settings,
rank count and a fingerprint of the tutorial inputs in
`verification_member.json`; `--reuse` reuses a member only when these match
and its solver log ends normally, and reruns it otherwise.

## Benchmark definition and reference quality

The benchmark was introduced by Wall (1999) and Mok (2001) and has since been
used in many studies of partitioned coupling. All sources agree on the cavity
(1 m x 1 m), the lid velocity $$1 - \cos(2 \pi t / 5)$$ m/s, the fluid
($$\rho_f = 1$$ kg/m$$^3$$, $$\nu_f = 0.01$$ m$$^2$$/s), the membrane
($$h = 0.002$$ m, $$E = 250$$ Pa, $$\nu_s = 0$$, $$\rho_s = 500$$ kg/m$$^3$$)
and the monitored midpoint displacement. They disagree on the openings, and
the published curves fall into two groups (see the figure below).

- **Group A**: Mok (2001, Bild 6.6), Gerbeau and Vidrascu (2003, Fig. 4) and
  Wall (1999, Bild 7.16). Mok's thesis defines the openings only through the
  mesh: its figures show in- and outflow arrows at the top corners, and
  Küttler and Wall (2008) state that "at each side there are two unconstrained
  nodes". Neither Mok nor Wall states the opening size, the inflow profile or
  how the fluid pressure level is fixed; only the external pressure on the
  underside, $$P_{ext} = 0$$, is given. The opening therefore depends on the
  mesh, and the group is not self-consistent: Wall's 64 x 64 solution has a
  peak of 0.177 m, Mok's 32 x 32 solution 0.215 m. Mok's thesis does not give a
  better-defined alternative to group B.
- **Group B**: Valdés (2007, Fig. 5.12), Tiba et al. (2026, Table 3) and an
  independent Kratos Multiphysics solution with the openings fixed (Taran
  2026, below). They prescribe a linear inflow
  $$v_x = \bar{v}(8y - 7)$$ over a fixed 0.125 m opening on the left and
  $$p = 0$$ on the right opening, which fully specifies the problem. Valdés
  notes that his solution differs from Mok's because not all of Mok's boundary
  conditions are available. This study uses the group B definition.
- The Kratos Multiphysics example (Zorrilla) was previously counted in group B.
  It fixes only the horizontal mesh displacement of the side walls, so its
  inlet and outlet nodes slide up the walls with the mesh, and at the peak
  both openings are about 15% smaller than 0.125 m. It therefore does not solve
  the group B problem, and its peak lies 9-11% below the group B solutions at
  every level. It is kept for context only.
- Kassiotis et al. (2011, Fig. 7a) use OpenFOAM with two unconstrained faces on
  each side wall and $$p = 0$$ at the outlet; their curve lies between the
  groups. They also note that imposing the pressure at both the inlet and the
  outlet gives a negative mean displacement.

Two pilot runs on the 32 x 32 mesh with the standard solid confirm which
definition each group corresponds to:

| Openings | Peak (m) | Trough (m) | Mean (m) | Reproduces |
|---|---|---|---|---|
| Linear inflow, $$p = 0$$ outlet (group B) | 0.289 | 0.199 | 0.243 | group B |
| Both openings traction free | -0.087 | -0.164 | -0.126 | neither group |

With traction-free openings on both sides the membrane moves down, out of the
cavity, as Kassiotis et al. observed, and matches neither group.

![Published midpoint histories and the two pilot runs](reference/references_and_pilots.png)

The curves are stored in `reference/`, each with a header recording its
source, figure and extraction method. The Valdés, Mok, Wall, Gerbeau and
Vidrascu, and Kassiotis curves are extracted exactly from the vector paths of
the PDF figures and calibrated against each figure's axis ticks or labels. The
Kratos example curve is digitised from the example's 2400 x 1800 PNG, to
about 0.5% of the peak. The Valdés and Mok extractions agree with the independent
digitisations shipped with the CoCoNuT example to 1 mm and 2.3 mm (RMS). Over
the periodic window $$t = 20$$-$$70$$ s the references give:

| Reference | Group | Peak (m) | Trough (m) | Mean (m) |
|---|---|---|---|---|
| Valdés (2007) | B | 0.2695 | 0.1967 | 0.2321 |
| Kratos, openings fixed (Taran 2026), 128 x 128 | B | 0.2835 | 0.2057 | 0.2439 |
| Tiba et al. (2026), Table 3 only | B | 0.2769 | - | 0.2466 |
| Kratos example (Zorrilla), moving openings | - | 0.2619 | 0.2033 | 0.2316 |
| Kassiotis et al. (2011) | - | 0.2433 | 0.2066 | 0.2255 |
| Mok (2001) | A | 0.2149 | 0.1426 | 0.1773 |
| Gerbeau and Vidrascu (2003) | A | 0.2132 | 0.1545 | 0.1828 |
| Wall (1999), to 50 s | A | 0.1770 | 0.1474 | 0.1610 |

Within group B the peak spreads by 5.2% (Valdés to Kratos), the trough by
4.6% and the mean by 6.2%. None of the published sources reports a mesh or
time-step study. The periods agree (5 s), and in every curve the membrane
rises into the cavity.

The fixed-opening Kratos solution is the only group B reference with its own
convergence evidence. It was computed with Kratos Multiphysics 10.3.0, which
shares no code with solids4foam: a P1/P1 VMS-stabilised fluid with Bossak time
integration, a Q1 total Lagrangian St Venant-Kirchhoff solid, and
Dirichlet-Neumann coupling with MVQN at a tolerance of $$10^{-6}$$ m, met in
every step. The side-wall mesh nodes with $$y \ge 0.875$$ are fixed in both
directions, so the openings keep their size. On fluid meshes from 16 x 16 to
128 x 128 at $$\Delta t = 0.0125$$ s the trough converges at order 1.9 and the
mean at 1.35, but the peak, monotone only from 32 x 32, at order 0.4: the
peak extrapolates to about 0.280 m, below the 128 x 128 value of 0.2835 m.
The time step, the coupling tolerance and the solid resolution through the
thickness each change the response by less than $$10^{-5}$$ m or 0.03%. The
one model-form difference from solids4foam is the outlet: Kratos leaves the
outlet velocity to its natural condition, whereas solids4foam imposes a zero
gradient, and imposing $$u_y = 0$$ at the Kratos outlet changes the peak,
trough and mean by 0.5%, 1.7% and 1.0%. The stored curve is the 128 x 128
solution, written every fourth time step; the runs, inputs and scripts are in
the `paper-fsi-verification` repository,
`paper-data/kratos_mok_lid_driven_cavity/`.

## Studies

All members use the IQN-ILS coupling with the interface predictor, the FSI
tolerance $$10^{-5}$$, `pimpleFluid` with the second-order backward scheme, and
a uniform mesh-motion diffusivity, as the tutorial does. The comparison window
is $$t = 20$$-$$70$$ s, after the start-up transient. For each member the driver
records the mean peak and trough per lid period, the time-mean, their distance
outside the group B envelope, and the largest distance of the history from the
band between the Valdés and fixed-opening Kratos curves. In the tables below,
"FSI its" is the mean/maximum number of coupling iterations per time step, and
"Time" the wall-clock time in seconds.

- `mesh`: the fluid and solid meshes are refined together (16, 32, 64 and 96
  cells across the cavity and along the membrane, eight cells through the
  membrane), with the time step reduced with the cell size (0.2, 0.1, 0.05 and
  0.025 s), on the high-order solid.
- `timestep`: $$\Delta t = 0.1$$, 0.05 and 0.025 s on the 64 mesh, high-order
  solid.
- `solid`: the standard solid with two and eight cells through the thickness
  against the high-order solid with eight, on the 64 mesh at
  $$\Delta t = 0.1$$ s.

### Mesh study (high-order solid, eight cells through the membrane, 8 ranks)

| Mesh | $$\Delta t$$ (s) | Peak (m) | Trough (m) | Mean (m) | FSI its | Time |
|---|---|---|---|---|---|---|
| 16 x 16 | 0.2 | 0.2934 | 0.1943 | 0.2425 | 5.3/26 | 123 |
| 32 x 32 | 0.1 | 0.2889 | 0.1986 | 0.2426 | 2.7/14 | 227 |
| 64 x 64 | 0.05 | 0.2870 | 0.2022 | 0.2437 | 2.3/6 | 26096 |
| 96 x 96 | 0.025 | 0.2861 | 0.2044 | 0.2445 | 2.0/6 | 39963 |

The 16 x 16 and 32 x 32 levels were run serially: on eight ranks these small
meshes leave each rank so few cells that the high-order solid gives wrong
results (solids4foam issue #696). The clock times are explained under
"Recorded results".

The response changes by 1.58%, 1.25% and 0.77% of the peak between successive
levels. The mean is mesh-independent to 0.8%, and the peak converges at an
observed order of about 1.2 (from the 16, 32 and 64 levels) towards about
0.285 m. The trough converges more slowly, at an observed order below 0.5, and
is still rising by 0.8% of the peak at the finest level; it is uncertain by
about 1-2% of the peak. The time step matters much less than the mesh (see
below), so halving it with each mesh level only removes a small contribution.

A 128 x 128 level was attempted and could not be completed. With both the
standard and the high-order solid, with $$\Delta t = 0.1$$ and 0.025 s, with an
initial relaxation factor of 0.02 and with eight PIMPLE outer correctors, the
coupled iteration diverges within a single time step and the solid PETSc SNES
solve then fails: at $$t = 4.3$$ s with $$\Delta t = 0.1$$ s, and at
$$t = 13.925$$ s in all three runs with $$\Delta t = 0.025$$ s. The 96 x 96
level therefore stands in for it. The cause was not isolated, and the
divergence is recorded as a known limitation.

### Time-step study (64 mesh, high-order solid, 4 ranks)

| $$\Delta t$$ (s) | Peak (m) | Trough (m) | Mean (m) | FSI its | Time |
|---|---|---|---|---|---|
| 0.1 | 0.2870 | 0.2026 | 0.2440 | 3.1/25 | 16558 |
| 0.05 | 0.2870 | 0.2022 | 0.2437 | 2.3/7 | 18020 |
| 0.025 | 0.2870 | 0.2020 | 0.2436 | 1.9/6 | 21885 |

The response changes by 0.16% and then 0.06% of the peak: the benchmark time
step of 0.1 s is already converged to well within the reference spread.

### Solid discretisation study (64 mesh, $$\Delta t = 0.1$$ s, 2 ranks)

| Solid | Cells in $$h$$ | Peak (m) | Trough (m) | Mean (m) | FSI its | Time |
|---|---|---|---|---|---|---|
| standard | 2 | 0.2869 | 0.2027 | 0.2440 | 3.2/26 | 2791 |
| standard | 8 | 0.2870 | 0.2027 | 0.2441 | 3.5/32 | 3323 |
| high-order | 8 | 0.2870 | 0.2026 | 0.2440 | 3.4/34 | 5316 |

The standard and the high-order solids give the same response to 0.03% of the
peak, and the standard solid needs only two cells through the thickness. The
membrane carries the pressure almost entirely through tension (its bending
stiffness is $$h^2 / (12 L^2) \approx 3 \times 10^{-7}$$ of $$EA L^2$$), so the
second-order finite volume solid does not suffer the slow convergence seen for
bending-dominated slender structures. The high-order solid gives no accuracy
benefit here and costs about 1.6 times as much. It is also more demanding of the
mesh: on the 32 mesh, with two or four cells through the thickness its PETSc
SNES solve stalls within the first second, whereas eight cells work, so the
tutorial and the studies use eight. The high-order solid is the tutorial default
because it was requested; `./Allrun standard` selects the cheaper standard solid
with the same result.

### FSI iterations

With IQN-ILS and the predictor, a time step needs 2.7 coupling iterations on
average on the 32 mesh at $$\Delta t = 0.1$$ s (at most 14), and 3.1 on the 64
mesh at $$\Delta t = 0.1$$ s; halving the time step reduces this to 2.3 and
1.9. Aitken relaxation (`./Allrun aitken`) gave the same response with 9.9
iterations per step on the 32 mesh in an earlier run with OpenFOAM v2412. For
comparison, at $$\Delta t = 0.1$$ s Mok (2001) reports 15 iterations in the
first step without relaxation and 8 with the best fixed
relaxation factor 0.825, with Aitken and steepest-descent relaxation 1-2
iterations cheaper than that over the run; Küttler and Wall (2008, Fig. 4,
tolerance $$10^{-7}$$, constant predictor) report about 5-7 iterations for
Aitken and 2-4 for their Newton-Krylov variants; Kassiotis et al. (2011) report
a mean of 17 for Aitken at tolerance $$10^{-7}$$; Valdés (2007) reports 6-8
for Aitken at $$\Delta t = 0.01$$ s. Bogaers et al. (2014) do not use this
case. The residual definitions and tolerances differ, so the comparison is
indicative only: IQN-ILS needs fewer iterations than the published
relaxation methods, whereas the solids4foam Aitken implementation needs more.

## Acceptance criteria

`Allverify` returns zero only when every check passes. The numbers are in
`reference/mokLidDrivenCavity_verification_references.json`.

- `mesh` and `timestep`: the change in the periodic response (the larger of
  the peak and trough changes) between successive members decreases, or is
  already below 0.1% of the peak, and the last change is below 1% of the peak.
- `solid`: every member agrees with the high-order solid to 1% of the peak.
- For the finest `mesh` and `timestep` members and every `solid` member, the
  result is compared with the group B references rather than with any single
  one of them:
  - the peak, trough and mean must lie within 5% of the envelope spanned by
    Valdés, the fixed-opening Kratos solution and Tiba et al. (Tiba gives
    only the peak and the mean): at most 5% above its maximum and at most 5%
    below its minimum. The envelopes are 0.2695-0.2835 m for the peak,
    0.1967-0.2057 m for the trough and 0.2321-0.2466 m for the mean;
  - after $$t = 20$$ s, the history must stay within 10% of the peak of the
    band between the Valdés and fixed-opening Kratos curves, which at each
    time runs from the lower to the higher of the two. The peak used is the
    smaller of the two curve peaks.

The published group B values are single-mesh results, so the checks accept any
result within the spread of the independent group B codes plus a 5% allowance,
and reject results that leave it. The fixed-opening Kratos solution, with its
own convergence study, is the closest of them to the converged answer.

A `--quick` run exercises the two coarsest members of each study to
$$t = 10$$ s, before the periodic state, and checks only that they complete.

## Interpretation

- None of the published references demonstrates mesh or time-step
  insensitivity: each gives one solution on one mesh (32 x 32 fluid elements
  in Mok, Valdés, the Kratos example and Kassiotis et al., 64 x 64 in Wall) at
  one time step. The fixed-opening Kratos solution is the exception.
- The solids4foam solution is mesh- and time-step-converged: the peak changes
  by 0.3% of itself between the 64 and 96 meshes (the whole response by
  0.8% of the peak), and by 0.06% of the peak between
  $$\Delta t = 0.05$$ and 0.025 s. The standard and the high-order solids
  agree to 0.03%.
- The fixed-opening Kratos solution agrees with solids4foam. On the same
  combined sequence (16 x 16 at 0.2 s to 96 x 96 at 0.025 s) the peak differs
  by -2.0%, -0.25%, -0.70% and -0.77%, the trough by +5.7%, +4.6%, +2.0% and
  +0.71%, and the mean by +1.3%, +1.9%, +0.48% and -0.13% (Kratos relative to
  solids4foam). The trough and mean differences shrink with refinement, the
  peak difference stays below 1% from 32 x 32 on, and the largest difference
  between the 96 x 96 histories is 2.4 mm, 0.8% of the peak. This is within
  the resolution of about 1% set by the different outlet conditions.
- Both codes therefore place the converged peak and trough above the
  published group B values of Valdés and Tiba et al., and the mean inside
  their range. The published solutions use numerically dissipative time
  integrators (generalised-alpha and Bossak schemes) on single, coarse meshes,
  and appear to be under-resolved.
- The peak converges slowly in both codes (observed order about 1.2 in
  solids4foam and 0.4 in Kratos at fixed time step), so the remaining
  difference of the 96 x 96 solids4foam peak from the extrapolated Kratos peak
  (0.9-2.2%) is mostly the solids4foam spatial error, not a disagreement
  between the codes.

## Recorded results

Recorded with OpenFOAM v2512 on Linux, on the UCD Sonic cluster (Intel Xeon
Gold 6542Y nodes), with solids4foam `development` at `206437983` plus the
outlet condition of PR #676 and this branch's verification scripts. The
studies ran concurrently in two Slurm jobs: `mesh` on eight ranks, with the
16 x 16 and 32 x 32 members rerun serially (`--cores 1`) because of
solids4foam issue #696; `timestep` on four ranks; and `solid` on two. The
combined summary was then written by `./Allverify --reuse`, which passes every
check. The eight-rank `mesh` pass on its own reports FAIL only because its two coarsest
members cannot complete there; the serial pass completes them.

The margins of the finest members to the acceptance limits are:

| Member | Check | Value | Outside envelope | Margin to limit |
|---|---|---|---|---|
| mesh96 | peak | 0.2861 m | +0.92% | 4.08% |
| mesh96 | trough | 0.2044 m | inside | 5.67% |
| mesh96 | mean | 0.2445 m | inside | 5.86% |
| mesh96 | history band | 1.11% of peak | - | 8.89% |
| dt0.025 | peak | 0.2870 m | +1.23% | 3.77% |
| dt0.025 | trough | 0.2020 m | inside | 6.82% |
| dt0.025 | mean | 0.2436 m | inside | 6.20% |
| dt0.025 | history band | 1.42% of peak | - | 8.58% |
| highOrder_ny8 | peak | 0.2870 m | +1.25% | 3.75% |
| highOrder_ny8 | trough | 0.2026 m | inside | 6.52% |
| highOrder_ny8 | mean | 0.2440 m | inside | 6.06% |
| highOrder_ny8 | history band | 1.37% of peak | - | 8.63% |

Every result lies inside, or within 1.3% of, the group B envelope; the peak
is the binding check. The clock times in the study tables are not
representative of the cost: the parallel members shared their nodes and were
slowed by MPI overhead on Sonic. The serial 16 x 16 and 32 x 32 members took
123 s and 227 s; in an earlier run on an Apple Silicon Mac Studio with two
ranks (OpenFOAM v2412), the 64 x 64 and 96 x 96 members took 827 s and
3729 s.

With `--quick` the driver runs the two coarsest members of each study to
$$t = 10$$ s, which took about 7 minutes on two ranks.

## Known limitations

- A 128 x 128 mesh level could not be completed (see the mesh study), so the
  finest mesh level is 96 x 96.
- The study was run on OpenFOAM.com only. The tutorial also runs on
  OpenFOAM.org and foam-extend, but their results differ (see the tutorial
  README): on foam-extend the midpoint displacement at t = 5 s is 3% smaller
  and the vertical force on the membrane about a tenth.
- Two group B references give the full history: Valdés, extracted from
  vector data, and the fixed-opening Kratos solution, as computed. Tiba et al.
  give scalar values only. The Kratos and solids4foam outlet conditions
  differ, which limits their comparison to about 1%.

`verification/Allclean` removes the generated `work/` and `postProcessing/`
directories, as does `./Allclean all` in the case directory; the plain case
`./Allclean` keeps them.

## References

- Wall, W.A. (1999). Fluid-Struktur-Interaktion mit stabilisierten Finiten
  Elementen. PhD thesis, Universität Stuttgart.
- Mok, D.P. (2001). Partitionierte Lösungsansätze in der Strukturdynamik und
  der Fluid-Struktur-Interaktion. PhD thesis, Universität Stuttgart.
- Gerbeau, J.-F., Vidrascu, M. (2003). A quasi-Newton algorithm based on a
  reduced model for fluid-structure interaction problems in blood flows.
  INRIA Research Report 4691.
- Förster, C., Wall, W.A., Ramm, E. (2007). Artificial added mass instabilities
  in sequential staggered coupling of nonlinear structures and incompressible
  viscous flows. Comput. Methods Appl. Mech. Engrg. 196:1278-1293.
- Valdés Vázquez, J.G. (2007). Nonlinear analysis of orthotropic membrane and
  shell structures including fluid-structure interaction. PhD thesis,
  Universitat Politècnica de Catalunya.
- Küttler, U., Wall, W.A. (2008). Fixed-point fluid-structure interaction
  solvers with dynamic relaxation. Comput. Mech. 43:61-72.
- Kassiotis, C., Ibrahimbegovic, A., Niekamp, R., Matthies, H.G. (2011).
  Nonlinear fluid-structure interaction problem. Part I: implicit partitioned
  algorithm, nonlinear stability proof and validation examples. Comput. Mech.
  47:305-323.
- Kratos Multiphysics Examples, FSI lid driven cavity (R. Zorrilla).
- Taran, A. (2026). Kratos Multiphysics 10.3.0 solution of the group B lid
  driven cavity with fixed openings. paper-fsi-verification,
  `paper-data/kratos_mok_lid_driven_cavity/`.
- Tiba, A. et al. (2026). Online adaptive non-intrusive model reduction via
  manifold interpolation and subspace updates: application to FSI convergence
  acceleration. arXiv:2609.16876.
