# 3dTube verification studies

This directory contains opt-in verification studies for the `3dTube` tutorial:
a pressure pulse travelling along a thick-walled elastic tube, computed with
the partitioned Robin-Neumann coupling (the tutorial default) and, for
comparison, the Dirichlet-Neumann IQN-ILS and Aitken couplings. The studies
are deliberately separate from `regressionTest.sh`: the regression test checks
that the tutorial remains numerically stable, whereas these studies check
convergence and agreement with published results. Nothing here is run by
`tutorials/Alltest` or `tutorials/Alltest-regression`.

**Quality of the references.** No publication tabulates results for this
case. The published point-A histories used here are read from figures: the
Tuković et al. (2018) curve exactly from the vector graphics of the PDF, and
the Lozovskiy et al. (2019) and Eken (2016) curves from raster figures, to
about ±3e-6 m and ±0.05 ms. The two independent monolithic references (finite
elements in Lozovskiy et al.; a side-centred finite volume fluid with a finite
element wall in Eken) use first-order implicit Euler with `Δt = 1e-4 s`, which
damps the pulse peak by about a quarter, so they are compared only with a run
that uses the same time discretisation; the backward, small-`Δt` results are
compared with Tuković et al. (2018), an earlier version of the same code. The
analytical pulse-wave speed is a thin- or thick-wall long-wave estimate and is
only an approximate check. See [Reference quality](#reference-quality) and,
for the three-level results and their limitations,
[Reference results](#reference-results).

## Running

Source a supported OpenFOAM environment, build solids4foam with PETSc, and run
from this directory:

```bash
cd tutorials/fluidSolidInteraction/3dTube/verification
./Allverify                        # mesh study, levels 1 and 2 (Robin-Neumann)
./Allverify --levels 1,2,3         # add the 4x mesh (very expensive)
./Allverify --study timestep       # time-step study on the tutorial mesh
./Allverify --study literature     # the published Euler, dt = 1e-4 s setup
./Allverify --study coupling       # Robin-Neumann vs IQN-ILS on the tutorial mesh
./Allverify --study coupling --couplings robin,iqnils,aitken
./Allverify --coupling iqnils      # any sweep with another coupling
./Allverify --quick                # smoke run of level 1 to t = 10 ms, no checks
./Allverify --reuse                # re-evaluate completed runs
./Allverify --cores 6              # MPI ranks per run (default: auto)
./Allverify --levels 1,2,3 --delta-t 2.5e-5   # refine the mesh at a fixed time step
./Allverify --probe-z-shift 7.8125e-5         # axis probes off the cell faces (c_p)
./Allverify --tight-tolerances     # diagnostic: bound the iterative error
```

The driver requires Python 3.8 or newer, `blockMesh` and `solids4Foam`;
`gnuplot` (5.0 or newer) is optional and is used for the history plots. Each run
is a complete copy of the tutorial under `verification/work/`, prepared by
editing the copy and then running the tutorial's own `Allrun` (with
`dirichletNeumann` for IQN-ILS and Aitken), so the tutorial itself and its
regression test are not modified. Results are written to
`verification/postProcessing/` as CSV files, `verification_summary.md`, and PNG
and PDF history plots. Both directories are ignored by Git and are retained to
make a failed run diagnosable. `Allverify` returns zero only when every
acceptance check of the selected study passes.

A completed run is reused by `--reuse` only when the settings recorded in its
`verification_settings.json` match the requested ones: coupling, refinement,
time step, end time, time scheme, solid preconditioner, monitors, MPI ranks,
OpenFOAM version and a hash of the tutorial's `0`, `constant`, `system` and
`Allrun`. Otherwise it is run again. Every monitored history must also have
one complete, finite row per time step up to the end time, the residual
evaluated for each step must be that of its last FSI iteration, and the solver
log must contain `End`; a run that fails these checks is reported as a failure
rather than evaluated.

## Case and verification copies

The tutorial reproduces the case as defined by Fernández and Moubachir (2005)
and used by Küttler and Wall (2008), Degroote et al. (2009, 2010), Tuković et
al. (2018), Lozovskiy et al. (2019) and others: a tube of length 5 cm, inner
radius 0.5 cm and wall thickness 0.1 cm, clamped at both ends; fluid density
1000 kg/m³ and viscosity 3e-6 m²/s; wall density 1200 kg/m³, `E = 3e5 Pa` and
`ν = 0.3`; an inlet pressure of 1333.2 Pa for 3 ms and zero outlet pressure. A
quarter of the tube is modelled. Point A is on the inner wall at mid-length,
`(0, 0.005, 0.025) m`, where `u_y` is the radial and `u_z` the axial
displacement.

The verification copies differ from the tutorial as follows:

- the fluid and solid use the second-order `backward` scheme (`--time-scheme
  Euler` restores the tutorial's first-order scheme); the `literature` study
  uses `Euler`, as the published finite element results do;
- the solid uses an algebraic multigrid (`hypre`) preconditioner instead of the
  block-Jacobi LU of the tutorial, whose factorisation does not scale to the
  finer meshes (`--solid-preconditioner lu` restores it); the converged
  solution does not depend on it;
- for IQN-ILS, the fluid-solid interface `predictor` is switched on, which
  avoids a first-iterate added-mass spike (issue #489); Aitken predicts the
  solid itself (`predictSolid`, on by default), and the Robin-Neumann setup is
  used as the tutorial has it;
- wall-displacement monitors at `z = 0.5, 1.0, ..., 4.5 cm` on the inner wall
  and pressure probes near the axis at the same positions are added, and the
  fields are written only at the end time;
- for more than one MPI rank, the copy's `Allrun` decomposes both regions and
  runs the solver in parallel (the tutorial `Allrun` is serial).

Two differences from the published definitions are not changed: the wall is
Hookean small-strain elastic, as in the linear elasticity of Fernández and
Moubachir (2005) and Deparis et al. (2015), whereas the other studies use St
Venant-Kirchhoff (the peak hoop strain is about 3%); and the inlet pressure is
ramped to zero between 3.0 and 3.1 ms rather than switched off, as in the
preCICE `elastic-tube-3d` tutorial.

## Monitored quantities

- `u_r,max(A)` and its time: the peak radial displacement at A, refined with
  a parabola through the largest sample and its neighbours;
- `u_z,min(A)`: the incident trough of the axial displacement at A, the
  minimum over `t < 6.5 ms` (at about 4.7 ms). The window matters: over the
  whole run the axial displacement has a second, reflected trough of similar
  depth at about 15.5 ms, and a global minimum switches between the two
  troughs from run to run. An earlier draft of the driver used the global
  minimum, which gave -0.0938 mm on level 1 (the reflected trough) against
  -0.0882 mm in the 10 ms quick run and a spurious 7.4% change between mesh
  levels; the incident trough is the well-defined quantity, and the one
  compared with the published values;
- `u_r,min(A), t > 14 ms`: the trough of the radial displacement after the
  reflection from the outlet (at about 17-18 ms);
- `t_arr(A)`: the arrival of the wave front at A, the time at which the radial
  displacement first reaches half of its first peak. It is used for the
  timing convergence checks instead of the time of `u_r,max(A)`: the peak is
  flat, so its time moves by tens of microseconds for a 0.01% change in the
  shape of the curve (7.204 ms at `Δt = 2.5e-5 s` against 7.242 ms at
  `Δt = 1.25e-5 s` for the same peak value to 0.02%);
- `c_p` and `c_wall`: the pulse-wave speed, the least-squares slope of the
  front arrival time against position for the stations `z = 1-4 cm`. For
  `c_p`, the arrival is the time at which the axis pressure first reaches half
  the inlet pressure (a fixed level, because the probes sample the containing
  cell and step slightly as the fluid mesh moves); for `c_wall`, the time at
  which the wall radial displacement reaches half of its first peak. `c_p` is
  the checked quantity: the wall signal also carries the faster axial stress
  wave of the wall. The stations `z = 5, 10, ..., 45 mm` lie on cell faces of
  every mesh level, so the cell that contains a default probe is a tie that
  different builds break differently (`c_p` = 4.767 m/s with OpenFOAM v2412
  and 4.660 m/s with v2512 on the same level-1 case). Use
  `--probe-z-shift 7.8125e-5` (half a level-3 cell) for `c_p`; even then the
  probes return the value of the containing cell, and the axial wall motion
  moves them between cells, so `c_p` carries a sampling noise of about one
  axial cell per station (`Δz/c`, 0.03 ms on level 3, about 0.5-1% in `c_p`).

## Studies and acceptance criteria

The reference values and tolerances are in
`reference/3dTube_verification_references.json`. Every study requires all its
runs to reach the end time (`t = 20 ms`) and every time step to satisfy the
coupling's convergence criteria: for Robin-Neumann, the
`robinConvergenceState` recorded in `postProcessing/fsiResiduals.dat` (the
displacement, pressure-change and leakage-flux residuals below their
tolerances, or stalled within the stall tolerance); for IQN-ILS and Aitken,
the final FSI residual below `outerCorrTolerance`.

### Mesh study (`--study mesh`, default)

Level 1 is the tutorial mesh and time step. Each further level doubles the
block divisions of both meshes in all three directions and halves the time
step:

| Level | Refinement | Fluid cells | Solid cells | Δt (s) | Default cores |
|---:|---:|---:|---:|---:|---:|
| 1 | 1x | 16 000 | 6 400 | 2.5e-5 | 4 |
| 2 | 2x | 128 000 | 51 200 | 1.25e-5 | 6 |
| 3 | 4x | 1 024 000 | 409 600 | 6.25e-6 | 32 |

Levels 1 and 2 form the default sweep; level 3 is opt-in and costs about
6 h on 32 ranks (720 core-hours). Because the time step is halved with the
mesh, the default sweep is a combined space-time refinement; `--delta-t`
holds the time step fixed so that the mesh alone is refined. The acceptance criteria are:

- the change between the two finest levels is at most 5% for `u_r,max(A)` and
  `u_z,min(A)`, 2% for `t_arr(A)` and 3% for `c_p`, and, with three or more
  levels, each successive change is smaller than the previous one or below
  0.2%, the resolution floor below which the order of the changes is noise;
- on the finest level, `u_r,max(A)` is within 5% of the Tuković et al. (2018)
  value of 0.15694 mm (backward, `Δt = 2.5e-5 s`) and its time within 3% of
  7.23 ms;
- on the finest level, `c_p` is within 10% of the thick-wall estimate with
  wall inertia, 4.81 m/s (see [Wave speed](#wave-speed)).

### Time-step study (`--study timestep`)

The tutorial mesh with `Δt = 1e-4, 5e-5, 2.5e-5` and `1.25e-5 s` (the
default members); a fifth member, `6.25e-6 s` (`--levels 5`), uses the
level-3 time step but is unstable on this mesh (see
[Small-time-step instability](#small-time-step-instability)). The change
between the two smallest time steps is at most 2% for `u_r,max(A)` and
`u_z,min(A)` and 1% for `t_arr(A)` and `c_p`; each successive change is
smaller than the previous one or below the 0.2% floor, and the observed order
from the three smallest time steps is reported where it is defined.

### Published-discretisation study (`--study literature`)

The tutorial mesh with first-order `Euler` and `Δt = 1e-4 s`, the time
discretisation of Lozovskiy et al. (2019) and Eken (2016). `u_r,max(A)` is
required within 8% of each published peak, its time within 5%, and
`u_z,min(A)` and `u_r,min(A), t > 14 ms` within 10%. The largest difference
from each published radial history is also reported. The tolerances cover the
spread between the two publications (4% in the peak) and the reading accuracy
of the raster figures.

### Coupling study (`--study coupling`)

The tutorial mesh and time step with Robin-Neumann and IQN-ILS coupling (and
Aitken on request). Both converge the same interface problem at every time
step, so the largest difference between the two radial histories at A,
normalised by the IQN-ILS peak, and the relative differences in
`u_r,max(A)`, `u_z,min(A)`, `t_arr(A)` and `c_p`, must each be within 1%. The
number of FSI iterations of each coupling is reported.

## Reference quality

- Tuković et al. (2018), Section 4.4, Fig. 25: finite volume, backward,
  `Δt = 2.5e-5` to `1e-4 s`, IQN-ILS, St Venant-Kirchhoff wall, full tube
  with 449 600 fluid and 288 000 solid cells, run to 10 ms; one mesh, with
  the peak given for three time steps (0.15203, 0.15631 and 0.15694 mm). Read
  exactly from the vector graphics of the PDF. An earlier version of this
  code (same lineage), so not an independent reference.
- Tuković et al. (2018), Eqs. (35)-(36): the analytical thick-wall wave speed,
  4.81 m/s. Printed.
- Tuković et al. (2018), Section 4.4: the simulated wave speed, 4.54 m/s.
  Printed, but obtained with another method; information only.
- Lozovskiy et al. (2019), Section 5.1, Fig. 3: monolithic finite elements
  (P2-P1 Taylor-Hood fluid, P2 displacement, Ani3D), St Venant-Kirchhoff wall,
  implicit Euler `Δt = 1e-4 s` with the geometry and advection linearly
  extrapolated from the previous steps (one linear solve per step), the
  pressure switched off instantaneously at 3 ms; three tetrahedral meshes
  (13 200/6 336, 29 202/11 904 and 89 232/38 016 fluid/solid cells), with a
  stated largest fine-to-finer history difference of 0.7% (axial) and 2.3%
  (radial); no time-step study. The `u_r(A)` and `u_z(A)` histories of the
  finer mesh to 20 ms, read from a raster figure, about ±3e-6 m. They state
  that their results are consistent with Eken and Sahin (2016).
- Eken (2016), Section 5.2, Fig. 5.6, and Eken and Sahin (2016), Section 3.2:
  monolithic ALE, side-centred (staggered) unstructured finite volume fluid
  with a Galerkin finite element St Venant-Kirchhoff wall, first-order implicit
  Euler in the fluid with `Δt = 1e-4 s` (chosen "to be consistent with Gee et
  al."), hexahedral meshes M1 (324 122 DOF) and M2 (2 557 571 DOF); no
  time-step study. The `u_r(A)` and `u_z(A)` histories of M2 to 20 ms, read
  from a raster figure, about ±3e-6 m.
- Moens-Korteweg and thick-wall formulas: long-wave analytical wave speeds of
  4.81-5.74 m/s. Approximate.

The histories are stored, in metres and resampled, in
`reference/Tukovic2018_fig25a_pointA.csv`,
`reference/Lozovskiy2019_fig3_pointA.csv` and
`reference/Eken2016_fig5.6_pointA.csv`; the JSON file records each value with
its figure and quality. Formaggia et al. (2001), Gerbeau and Vidrascu (2003),
Fernández and Moubachir (2005), Küttler and Wall (2008) and Degroote et al.
(2009, 2010) define or use the case but publish only iteration counts and
field snapshots. Formaggia et al. and Gerbeau and Vidrascu use a 5 ms pulse
and a shell wall; Degroote et al. use a shell wall. The preCICE
`elastic-tube-3d` tutorial uses the same parameters, but its reference plot
shows per-step displacement increments and is not a usable reference.

The implicit-Euler references and the backward finite volume reference differ
by about 21-24% in the peak radial displacement (0.119-0.124 mm against
0.157 mm). The split is between time discretisations rather than between
finite elements and finite volumes (the Eken fluid is itself a finite volume
discretisation). Much of it, and within the reading accuracy all of it for
`u_r,max(A)`, is reproduced by running this code with the published
`Euler`, `Δt = 1e-4 s` setup; see
[Implicit Euler and the published results](#implicit-euler-and-the-published-results).
The arrival of the Eken peak about 0.3 ms earlier is not explained by it.

### Wave speed

The Moens-Korteweg speed `c = sqrt(E h / (2 ρ_f R))` is 5.48 m/s with the inner
radius, 5.22 m/s with the mean radius and 5.74 m/s for an axially tethered
wall (`1/(1 - ν²)`). It assumes a thin wall, an inviscid fluid, a wavelength
much longer than the radius and no wall inertia. Here `h/R = 0.2`, and the
3 ms pulse is only about three radii long, so these values are upper bounds
rather than targets. The thick-wall (Lamé) estimate of Tuković et al. (2018)
Eq. (35) gives 5.07 m/s, and its correction for the axial stress waves in the
wall, Eq. (36), gives 4.81 m/s, which is used for the check with a 10%
tolerance.

## Reference results

Recorded in October 2026 with solids4foam branch `verification/3dtube-level3`
at commit `0ff7f462b`, which contains two fixes found by this study (see
[Platform dependence](#platform-dependence-and-the-two-fixes)), Robin-Neumann
coupling and the backward scheme unless stated. Levels 1 and 2 and the
level-1 studies: OpenFOAM v2512 (Ubuntu package) on a shared 192-core AMD
node; level 3: MeluXina, OpenFOAM v2412 (EasyBuild foss-2024a), 128 ranks.
The two platforms agree on levels 1 and 2 to 0.005% in every point-A
quantity. Displacements at point A are in mm, times in ms, speeds in m/s; `c_p`
from the probe-shifted runs (`--probe-z-shift 7.8125e-5`). The compact
evidence is in `results/` (`3dTube_evidence.json`, `3dTube_mesh_study.csv`,
`3dTube_runs.csv`, `3dTube_meluxina_runs.csv`), written by
`scripts/analyse_3dTube_evidence.py`.

**All values recorded before these fixes are superseded**, including the
earlier versions of this section and the Apple M1 values of the original
two-level study (which used the faulty point transfer on level 2).

### Mesh study: three levels

Time step halved with the mesh (`Iter.` is the mean number of FSI iterations
per step):

| Level | Fluid / solid cells | Δt (s) | Ranks | u_r,max | t(u_r,max) | u_z,min | t_arr | c_p | Iter. | Clock (s) |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | 16 000 / 6 400 | 2.5e-5 | 4 | 0.15978 | 7.201 | -0.08819 | 5.873 | 4.654 | 4.78 | 531 |
| 2 | 128 000 / 51 200 | 1.25e-5 | 8 | 0.15946 | 7.215 | -0.08659 | 5.843 | 4.635 | 3.87 | 3 744 |
| 3 | 1 024 000 / 409 600 | 6.25e-6 | 128 | 0.15851 | 7.252 | -0.08571 | 5.825 | 4.594 | 3.18 | 15 317 |

| Quantity | Change 1→2 | Change 2→3 | Observed order | Verdict |
|---|---:|---:|---:|---|
| u_r,max(A) | -0.20% | -0.60% | undefined | monotone, but the change grows: not asymptotic |
| u_z,min(A) | 1.87% | 1.02% | 0.87 | monotone, about first order (nominal 2) |
| t_arr(A) | -0.51% | -0.31% | 0.72 | monotone, sub-nominal |
| c_p | -0.41% | -0.88% | undefined | the change grows; of the size of the probe sampling noise |
| u_r,min(A), t > 14 ms | 6.3% | 2.0% | undefined | non-monotone |

With the time step fixed at `2.5e-5 s` on all three meshes (pure mesh
refinement; level 3 on MeluXina): `u_r,max(A)` 0.15978, 0.15921, 0.15835
(changes 0.36%, 0.54%; undefined); `u_z,min(A)` -0.08819, -0.08643, -0.08554
(2.0%, 1.03%; order 0.98); `t_arr(A)` 5.873, 5.837, 5.815 (0.61%, 0.38%; order
0.65).

Every step of every level met all three Robin criteria (worst final residuals
6.5e-7, 1.0e-5 and 9.1e-6 against tolerances 1e-6, 1e-5 and 1e-5; no stalled
steps); the iterations per step fall with refinement (4.78, 3.87, 3.18; at
most 9). Tightening every tolerance (`--tight-tolerances`) changed no point-A
quantity by more than 3e-6 (relative) on any level in the runs before the
fixes; the fixes do not change the solvers or the tolerances.

What this establishes:

- `u_r,max(A)` decreases monotonically but is not in an asymptotic range: the
  last change (0.6%, 0.5% at fixed time step) is larger than the first. No
  order can be given; the level-3 value carries an estimated discretisation
  uncertainty of about 0.6% (the last change; Roache's band with `Fs = 3`,
  `p = 2`) and plausibly more, since the changes are not yet decreasing.
- `u_z,min(A)` and `t_arr(A)` converge monotonically at about first order
  (0.87-0.98 and 0.65-0.72), well below the nominal second order. With
  `p = 1` the remaining error of `u_z,min(A)` on level 3 is about 1%
  (-0.0848 mm extrapolated), and of `t_arr(A)` about 0.3%.
- `c_p` stays within 4.59-4.66 m/s on every level and run, 3-5% below the
  thick-wall estimate of 4.81 m/s; its changes are of the size of the probe
  sampling noise and are not ranked.
- On level 3, `u_r,max(A)` is 1.0% above and its time 0.3% after the
  Tuković et al. (2018) value (0.15694 mm, 7.23 ms, `Δt = 2.5e-5 s`, an
  earlier version of this code on a full-tube mesh).

### Space-time separation

| Mesh | Δt change (s) | Δu_r,max | Δu_z,min | Δt_arr |
|---|---|---:|---:|---:|
| level 1 | 2.5e-5 → 1.25e-5 | +0.07% | -0.01% | -0.13% |
| level 2 | 2.5e-5 → 1.25e-5 | +0.16% | +0.18% | +0.10% |
| level 3 | 2.5e-5 → 6.25e-6 | +0.10% | +0.20% | +0.17% |

The time-step error at the level-2 and level-3 time steps is 0.1-0.2%,
against level-2-to-3 mesh changes at fixed time step of 0.4-1.0%: the scaled
sequence is dominated by the mesh, and the fixed-step sequence leads to the
same conclusions. On level 1 the backward time-step study gives successive
`u_r,max(A)` changes of 1.43%, 0.78% and 0.07% for `Δt = 1e-4` to
`1.25e-5 s`.

### Small-time-step instability

With the backward scheme, the time step must not be reduced far on a coarse
mesh. On level 1, `Δt = 6.25e-6 s` (the level-3 step) develops a growing
spurious mode from about 9 ms (the maximum Courant number, physically about
0.003, grows to 0.008 at 11 ms and 0.07 at 12 ms) and diverges at 13.9 ms,
although every coupling iteration converged; it is unchanged by the two fixes
below. At `Δt = 1.25e-5 s` the same mode is visible after 17 ms but stays
bounded to 20 ms, and with implicit Euler it does not appear. It depends on
`Δt` relative to the mesh, not on `Δt` alone: levels 2 and 3 at their scaled
time steps show no trace of it to 20 ms. The scaled sweep is not affected;
the cause (the Robin wall conditions with the backward scheme are a
candidate) has not been investigated.

### Platform dependence and the two fixes

The original three-level study gave a `u_z,min(A)` that differed between two
platforms by 1.5% on level 1 and 3% on level 3. Two defects were found and
fixed (details, evidence and reproducers in
[`platform/README.md`](platform/README.md)):

1. **Aliasing miscompilation of the fluid wall viscous force** (all levels,
   `-O3` builds only). `nf() & (-devReff().boundaryField())` writes the inner
   product into the storage of the tmp normal, which OpenFOAM's field loops
   access through `__restrict__` pointers; at `-O3` the axial component was
   formed from the overwritten normal. The axial wall shear was wrong by up
   to a factor of 2.5 on the first step. Builds with `-O2` (the EasyBuild
   OpenFOAM rules) were not affected. Fixed in every fluid model (commit
   `635ff464f`).
2. **Degenerate barycentric weights in the FSI point transfer** (levels 2 and
   3, all platforms). `triangle::pointToBarycentric` compares `4 A^2` (m^4)
   with `SMALL = 1e-15`, so the fan triangles of every interface face finer
   than about 0.3 mm got the weights (1/3, 1/3, 1/3): the solid displacement
   was transferred to the fluid as a smeared average on levels 2 and 3 but
   correctly on level 1, which broke the quarter-tube symmetry and made the
   result depend on decomposition and build. A scale-invariant weight
   (`triangleWeights.H`) is used instead (commit `0ff7f462b`).

With both fixes, the platforms agree on levels 1 and 2 to 0.005% in
`u_r,max(A)`, `u_z,min(A)` and `t_arr(A)` (the resolution of the reported
digits), and level 2 is independent of the number of ranks.

### Implicit Euler and the published results

| Run | u_r,max | t(u_r,max) | u_z,min | u_r,min, t > 14 ms |
|---|---:|---:|---:|---:|
| Euler, Δt 1e-4, level 1 | 0.12241 | 7.005 | -0.07480 | -0.07767 |
| Euler, Δt 1e-4, level 2 | 0.12100 | 7.035 | -0.07326 | -0.07628 |
| Euler, Δt 5e-5, level 1 | 0.13649 | 7.077 | -0.07969 | -0.08584 |
| Euler, Δt 2.5e-5, level 1 | 0.14651 | 7.121 | -0.08299 | -0.09048 |
| Euler, Δt 1.25e-5, level 1 | 0.15272 | 7.155 | -0.08517 | -0.09380 |
| Lozovskiy et al. (2019) | 0.1194 | 7.00 | -0.0740 | -0.0755 |
| Eken (2016) | 0.1242 | 6.69 | -0.0794 | -0.0763 |
| Backward, level 3 | 0.15851 | 7.252 | -0.08571 | -0.09603 |

With the published discretisation the peak is within +2.5%/-1.4% (level 1)
and +1.3%/-2.6% (level 2) of the two published peaks, against 33% and 28%
with the backward level-3 solution; the axial trough moves from 8-16% to
1-6% (level 1; 1-8% on level 2) of the published values and the late radial
trough from 26-27% to 0-3%. The damping by Euler at `Δt = 1e-4 s` (0.036-0.038 mm) amounts to
92-115% of the gap between the published peaks and the backward solution
(86-126% within the reading accuracy of the figures). The Euler results
converge slowly in time (successive `u_r,max(A)` changes of 11.5%, 7.3% and
4.2%, still 4.4% below the backward value at `Δt = 1.25e-5 s`), so neither
published finite element peak is converged in time. The time discretisation
therefore reproduces most, and within the reading accuracy all, of the
published difference in `u_r,max(A)`, and much of the difference in the
troughs; it does not explain the 0.3 ms earlier peak of Eken (2016), and the
published results also differ in their wall model (St Venant-Kirchhoff),
pressure switch-off (instantaneous) and meshes.

### Coupling study

Tutorial mesh, backward, `Δt = 2.5e-5 s`, 4 ranks:

| Coupling | u_r,max | u_z,min | t_arr | Iter. (mean / max) | Clock (s) |
|---|---:|---:|---:|---:|---:|
| Robin | 0.15978 | -0.08819 | 5.873 | 4.78 / 9 | 483 |
| IQN-ILS | 0.15969 | -0.08815 | 5.869 | 15.52 / 20 | 1 994 |

The radial histories at A differ by at most 0.26% of the peak and the
monitored quantities by at most 0.07%, below the level-2-to-3 changes.

## References

L. Formaggia, J.-F. Gerbeau, F. Nobile, A. Quarteroni, On the coupling of 3D
and 1D Navier-Stokes equations for flow problems in compliant vessels,
*Computer Methods in Applied Mechanics and Engineering* 191(6-7):561-582, 2001.

J.-F. Gerbeau, M. Vidrascu, A quasi-Newton algorithm based on a reduced model
for fluid-structure interaction problems in blood flows, *ESAIM: M2AN*
37(4):631-647, 2003.

M.A. Fernández, M. Moubachir, A Newton method using exact Jacobians for
solving fluid-structure coupling, *Computers & Structures* 83(2-3):127-142,
2005; INRIA research report RR-5085, 2004.

U. Küttler, W.A. Wall, Fixed-point fluid-structure interaction solvers with
dynamic relaxation, *Computational Mechanics* 43:61-72, 2008.

J. Degroote, K.-J. Bathe, J. Vierendeels, Performance of a new partitioned
procedure versus a monolithic procedure in fluid-structure interaction,
*Computers & Structures* 87:793-801, 2009.

J. Degroote, R. Haelterman, S. Annerel, P. Bruggeman, J. Vierendeels,
Performance of partitioned procedures in fluid-structure interaction,
*Computers & Structures* 88:446-457, 2010.

Ž. Tuković, A. Karač, P. Cardiff, H. Jasak, A. Ivanković, OpenFOAM finite
volume solver for fluid-solid interaction, *Transactions of FAMENA*
42(3):1-31, 2018, doi:10.21278/TOF.42301.

A. Lozovskiy, M.A. Olshanskii, Y.V. Vassilevski, Analysis and assessment of a
monolithic FSI finite element method, *Computers & Fluids* 179:277-288, 2019,
doi:10.1016/j.compfluid.2018.11.004.

A. Eken, *A parallel monolithic approach for the numerical simulation of
fluid-structure interaction problems*, PhD thesis, Istanbul Technical
University, 2016, Section 5.2; the study is also published as A. Eken, M.
Sahin, *International Journal for Numerical Methods in Fluids*
80(12):687-714, 2016.

S. Deparis, D. Forti, A. Quarteroni, A fluid-structure interaction algorithm
using radial basis function interpolation between non-conforming interfaces,
MATHICSE Technical Report 16.2015, EPFL, 2015.
