# hessenthalerFsi numerical verification

This directory contains the opt-in **numerical verification** of Phase I of
the `hessenthalerFsi` tutorial: how far the computed steady state is from the
exact solution of the mathematical model that the case defines. It is kept
separate from the **validation** in `../validation`, which compares the
model with the MRI measurements. The measurements are not used here, and
agreement with them is not evidence of numerical convergence.

Nothing here is run by `tutorials/Alltest` or `tutorials/Alltest-regression`.

## The verification problem

Every run solves the same mathematical problem; only numerical parameters
change. The problem is the tutorial's Phase I model with ONE fixed shear
modulus:

| Item | Definition |
| --- | --- |
| Geometry | the published fluid-domain surface (`geometry/fluidDomain.vtk.gz`) with the 250 mm outlet extension of `makeFluidSurface.py`; flap 11 x 2 x 65 mm, clamped at z = 0 |
| Fluid | incompressible Newtonian, ρ = 1163.3 kg/m³, μ_f = 12.50 mPa s, laminar, no gravity |
| Inlets | parabolic, peak 630 mm/s (upper) and 615 mm/s (lower), radius 10.95 mm, centres y = ±27.15 mm, ramp 3s² − 2s³ over 0.5 s |
| Outlet | p = 0, `inletOutlet` velocity |
| Walls, interface | no slip; moving-wall velocity on the flap |
| Solid | neo-Hookean (`neoHookeanElastic`), μ = μ* (below), ν = 0.45, K = 2μ(1+ν)/(3(1−2ν)), ρ_s = 1058.3 kg/m³ |
| Buoyancy | net body force on the solid, g_eff = (ρ_s − ρ_f)/ρ_s g = 0.973306 m/s² in +y, ramped with the inflow |
| Initial state | fluid at rest, straight flap |
| QoIs | steady state, as time averages over the last 20 % of the run |

μ* is held fixed for every mesh, time step and coupling setting. It is not
recalibrated per mesh: that would change the mathematical problem under
refinement. Its value and the evidence for it are in
[Fixed shear modulus](#fixed-shear-modulus).

The momentum stabilisation (`diffStencilLaplacian`, scale factor 0.01) is a
numerical term: a consistent face flux
`sf*(snGrad(D) − n·interp(grad D))` that vanishes as O(h²) for a smooth
displacement. It is therefore part of the discretisation, and its effect is
measured as part of the solid discretisation error.

## Running

Source OpenFOAM, build solids4foam with PETSc, and run from this directory:

```bash
./Allverify --study solid                      # zero-flow problem, S1-S4
./Allverify --study solid --levels S1,S2 --jobs 4
./Allverify --study meshes --threads 32        # build F1-F5, record quality
./Allverify --study run --spec F2:S2 --cores 48
./Allverify --study run --spec F2:S2:dt=0.001 --cores 48
./Allverify --study analyse                    # tables, orders, uncertainty
```

A coupled run is named by its spec: the fluid level `F1`-`F5`, the solid
level `S1`-`S4`, and optional `dt=`, `T=` (end time), `tol=` (coupling
tolerance), `coupling=robin`, `fluidtol=tight` (fluid linear solvers),
`pimple=` (PIMPLE outer correctors), `sf=` (stabilisation) and `mu=`.
Each run is a copy of the tutorial under `verification/work`; its evaluated
quantities are written to `verification/postProcessing/runs/<name>.json`.
Both directories are ignored by Git. `--reuse` re-evaluates a completed run
only if its settings, tutorial inputs and build match; `--evaluate-only`
re-evaluates an existing run directory.

The coupled runs are expensive (hours to a day on 24-512 cores each); the
campaign recorded below ran on MeluXina.

## Mesh families

**Solid** (nested hexahedra; S2 is the tutorial mesh):

| Level | Cells (x, y, z) | Cells | h_x, h_y, h_z (mm) | Through thickness |
| --- | --- | --- | --- | --- |
| S1 | 6 x 4 x 33 | 792 | 1.83, 0.50, 1.97 | 4 |
| S2 | 12 x 8 x 65 | 6 240 | 0.917, 0.25, 1.00 | 8 |
| S3 | 24 x 16 x 130 | 49 920 | 0.458, 0.125, 0.50 | 16 |
| S4 | 48 x 32 x 260 | 399 360 | 0.229, 0.0625, 0.25 | 32 |

S2 → S3 → S4 refine by exactly 2 in every direction. S1 → S2 is 2 in x and
y but 65/33 = 1.97 in z.

**Fluid** (cfMesh `cartesianMesh`; F1 is the tutorial's coarse mesh). Level
Fk scales every cell size of `system/meshDict.coarse`, including the octree
root size `maxCellSize`, by 2^(−(k−1)/2). Every region is refined
consistently by a nominal √2 per level, and F1 → F3 → F5 by 2. The meshes
are not nested (the octree is re-fitted on each level).

FLUID_MESH_TABLE

The tutorial's `medium` mesh (0.81 M cells) is not part of the family: it
refines only the wall and merging region and keeps the near-flap and flap
sizes of `coarse`.

## Quantities of interest

- tip displacement y and z of the flap tip centre (0, 0, 65 mm);
- the flap centreline: displacement at the material stations z₀ = 5, 10,
  …, 65 mm (exact nodes of S2-S4), compared between levels as the RMS and
  maximum of the displacement-vector difference;
- the fluid force on the flap (y and z), and the inlet pressures;
- the velocity on the MRI voxels of the planes z = 10 and 30 mm, averaged
  over each voxel (4 x 4 x 12 points) and over the time window at the FIXED
  voxel positions (sampled every 0.1 s); compared between levels as the mean
  and maximum vector difference normalised by 630 mm/s, plus the peak vz of
  each jet on each plane;
- FSI iterations, wall time, the residual transient and the deformed-mesh
  quality.

The centreline monitors are placed on mesh nodes, because
`solidPointDisplacement` reports the closest node. The run's final D field,
sampled at the stations, checks them.

## Steady state

The flap settles through a slowly decaying oscillation, so every QoI is a
time average over the final 20 % of the run, the end time is 15 s (not the
tutorial's 10 s), and each run estimates its own residual transient error.
The tip history is smoothed over one oscillation period and fitted with
a + b exp(−t/τ) over the second half of the run (and, as a sensitivity, its
last 40 % and 30 %). The distance between the window mean and the fitted
steady value a is the residual transient error. One run to 40 s tests the
estimator.

## Observed orders and uncertainty

For three levels with a constant nominal ratio r, the observed order is
computed from successive differences. It is reported only when the sequence
is monotone, the convergence ratio is in (0, 1), the order lies in
[0.5, 4], and the changes exceed the QoI's resolution (its residual
transient error). This is a screening criterion that is consistent with
entry into the asymptotic range; it does not prove it. Only then are a
Richardson estimate and a fine-level GCI (safety factor 1.25) reported.
Otherwise, a heuristic envelope of three times the largest successive
difference is given, and no order is claimed.

RESULTS
