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
about ±3e-6 m and ±0.05 ms. The published finite element results use
first-order implicit Euler with `Δt = 1e-4 s`, which damps the pulse peak by
about a quarter, so they are compared only with a run that uses the same time
discretisation; the converged (backward, small `Δt`) results are compared with
Tuković et al. (2018). The analytical pulse-wave speed is a thin- or
thick-wall long-wave estimate and is only an approximate check. See
[Reference quality](#reference-quality).

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
```

The driver requires `python3`, `blockMesh` and `solids4Foam`; `gnuplot` is
optional and is used for the history plots. Each run is a complete copy of the
tutorial under `verification/work/`, prepared by editing the copy and then
running the tutorial's own `Allrun` (with `dirichletNeumann` for IQN-ILS and
Aitken), so the tutorial itself and its regression test are not modified.
Results are written to `verification/postProcessing/` as CSV files,
`verification_summary.md`, and PNG and PDF history plots. Both directories are
ignored by Git and are retained to make a failed run diagnosable. `Allverify`
returns zero only when every acceptance check of the selected study passes.

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
  wave of the wall.

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
| 3 | 4x | 1 024 000 | 409 600 | 6.25e-6 | 6 |

Levels 1 and 2 form the default sweep; level 3 is about sixteen times as
expensive as level 2 and is opt-in. The acceptance criteria are:

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

The tutorial mesh with `Δt = 1e-4, 5e-5, 2.5e-5` and `1.25e-5 s`. The change
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
  `Δt = 2.5e-5` to `1e-4 s`; the `u_r(A)` history to 10 ms and its peak for
  three time steps. Read exactly from the vector graphics of the PDF.
- Tuković et al. (2018), Eqs. (35)-(36): the analytical thick-wall wave speed,
  4.81 m/s. Printed.
- Tuković et al. (2018), Section 4.4: the simulated wave speed, 4.54 m/s.
  Printed, but obtained with another method; information only.
- Lozovskiy et al. (2019), Section 5.1, Fig. 3: monolithic finite elements,
  implicit Euler, `Δt = 1e-4 s`; the `u_r(A)` and `u_z(A)` histories to 20 ms.
  Read from a raster figure, about ±3e-6 m.
- Eken (2016), Section 5.2, Fig. 5.6: monolithic finite elements, implicit
  Euler, `Δt = 1e-4 s`; the `u_r(A)` and `u_z(A)` histories to 20 ms. Read
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

The finite element and finite volume references differ by about 25% in the
peak radial displacement (0.119-0.124 mm against 0.157 mm). This is the time
discretisation, not the benchmark: the `literature` study reproduces the
finite element values with the published `Euler`, `Δt = 1e-4 s` setup, and the
mesh and time-step studies converge towards the Tuković et al. (2018) value
with the backward scheme.

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

Recorded with OpenFOAM v2412 on an Apple M1 Ultra shared with other jobs, so
the clock times are indicative only. Every study passed. Radial and axial
displacements are at point A, in mm; times in ms; speeds in m/s.

### Mesh study

Robin-Neumann, backward, 4 MPI ranks for both levels (`Iter.` is the mean
number of FSI iterations per time step):

| Level | u_r,max | t(u_r,max) | u_z,min | t_arr | c_p | Iter. | Clock (s) |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | 0.15984 | 7.204 | -0.08823 | 5.874 | 4.771 | 4.74 | 425 |
| 2 | 0.15976 | 7.206 | -0.08660 | 5.854 | 4.695 | 4.16 | 7 197 |

Between the levels, `u_r,max(A)` changes by 0.05%, `u_z,min(A)` by 1.85%,
`t_arr(A)` by 0.34% and `c_p` by 1.58%. On level 2, `u_r,max(A)` is 1.8% above
the Tuković et al. (2018) value of 0.15694 mm and its time 0.3% before
7.23 ms; the largest difference from their radial history over 0-10 ms is 3.9%
of the peak. `c_p` is 2.4% below the thick-wall estimate of 4.81 m/s and 3.4%
above the 4.54 m/s they report from a different evaluation of their
simulation. The Moens-Korteweg values of 5.2-5.7 m/s are 11-22% above `c_p`,
as expected for a wall with `h/R = 0.2`.

### Time-step study

Robin-Neumann, backward, tutorial mesh, 2 MPI ranks (`Iter.` as above):

| Δt (s) | u_r,max | t(u_r,max) | u_z,min | t_arr | c_p | Iter. | Clock (s) |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 1e-4 | 0.15581 | 7.224 | -0.08558 | 5.721 | 4.752 | 4.55 | 246 |
| 5e-5 | 0.15853 | 7.207 | -0.08789 | 5.848 | 4.768 | 4.26 | 466 |
| 2.5e-5 | 0.15984 | 7.204 | -0.08823 | 5.874 | 4.771 | 4.74 | 942 |
| 1.25e-5 | 0.15987 | 7.242 | -0.08821 | 5.866 | 4.764 | 7.36 | 1894 |

The successive changes in `u_r,max(A)` are 1.74%, 0.83% and 0.016%, and in
`t_arr(A)` 2.23%, 0.44% and 0.14%: the tutorial time step of `2.5e-5 s` is
converged in time to about 0.1%. The mesh therefore dominates the remaining
level-1 error. The time-step dependence of the peak agrees with Tuković et al.
(2018), whose peak rises from 0.15203 to 0.15694 mm over the same range of
time steps. The mean number of Robin-Neumann iterations rises from 4.7 to 7.4
at the smallest time step (not investigated further here); every step still
met all its criteria.

### Published-discretisation study

Robin-Neumann, `Euler`, `Δt = 1e-4 s`, tutorial mesh, 2 MPI ranks, 232 s:

| Quantity | solids4foam | Lozovskiy et al. (2019) | Eken (2016) |
|---|---:|---:|---:|
| u_r,max (mm) | 0.1224 | 0.1194 (+2.5%) | 0.1242 (-1.4%) |
| t(u_r,max) (ms) | 7.005 | 7.00 (+0.1%) | 6.69 (+4.7%) |
| u_z,min (mm) | -0.0748 | -0.0740 (+1.1%) | -0.0794 (-5.8%) |
| u_r,min after 14 ms (mm) | -0.0777 | -0.0755 (+2.9%) | -0.0763 (+1.8%) |

The largest difference from the published radial histories over 0-20 ms is
4.7% (Lozovskiy et al.) and 11.5% (Eken) of the published peak; the latter is
mostly the 0.3 ms earlier peak of the Eken curve. With this discretisation
the solids4foam peak is 23% below its converged backward value, which is the
source of the difference between the two families of published results.

### Coupling study

Tutorial mesh, backward, `Δt = 2.5e-5 s`, 2 MPI ranks (`Iter.` is the total
number of FSI iterations over the 800 time steps):

| Coupling | u_r,max | u_z,min | t_arr | c_p | Iter. | Mean | Max | Clock (s) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| Robin | 0.15984 | -0.08823 | 5.874 | 4.771 | 3 791 | 4.74 | 9 | 758 |
| IQN-ILS | 0.15971 | -0.08818 | 5.870 | 4.771 | 12 369 | 15.46 | 20 | 3 077 |

The two radial histories at A differ by at most 0.27% of the peak, and the
monitored quantities by at most 0.08%. Every Robin-Neumann step met its
displacement, pressure and leakage-flux criteria without stalling (worst
final residuals 6.1e-7, 9.9e-6 and 9.2e-6). The Robin-Neumann coupling needs
3.3 times fewer FSI iterations than IQN-ILS and runs 4.1 times faster on this
case, for which it was designed.

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
