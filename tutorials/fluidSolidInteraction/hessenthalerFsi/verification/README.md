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

- **Geometry:** the published fluid-domain surface
  (`geometry/fluidDomain.vtk.gz`) with the 250 mm outlet extension of
  `validation/scripts/makeFluidSurface.py` (committed as
  `constant/triSurface/fluidDomain.stl.gz`); flap 11 x 2 x 65 mm, clamped at
  z = 0.
- **Fluid:** incompressible Newtonian, ρ = 1163.3 kg/m³, μ_f = 12.50 mPa s,
  laminar, no gravity.
- **Inlets:** parabolic, peak 630 mm/s (upper) and 615 mm/s (lower), radius
  10.95 mm, centres y = ±27.15 mm, ramp 3s² − 2s³ over 0.5 s.
- **Outlet:** p = 0 and the backflow-stabilised traction
  σn = βρ min(u·n, 0) u with β = 0.2, as in the benchmark definition
  (Hessenthaler et al. 2017, Eq. 9; `backflowStabilisedOutletVelocity`). A
  zero-gradient velocity is a sensitivity check only.
- **Walls, interface:** no slip; moving-wall velocity on the flap.
- **Solid:** neo-Hookean (`neoHookeanElastic`), μ = μ* (below), ν = 0.45,
  K = 2μ(1+ν)/(3(1−2ν)), ρ_s = 1058.3 kg/m³.
- **Buoyancy:** net body force on the solid,
  g_eff = (ρ_s − ρ_f)/ρ_s g = 0.973306 m/s² in +y, ramped with the inflow.
- **Initial state:** fluid at rest, straight flap.
- **QoIs:** steady state, as time averages over the last 20 % of the run.

μ* is held fixed for every mesh, time step and coupling setting. It is not
recalibrated per mesh: that would change the mathematical problem under
refinement. Its value and the evidence for it are in
[Fixed shear modulus](#fixed-shear-modulus).

The momentum stabilisation (`diffStencilLaplacian`, scale factor 0.01) is a
numerical term: a consistent face flux
`sf*(snGrad(D) − n·interp(grad D))` that vanishes as O(h²) for a smooth
displacement. It is therefore part of the discretisation, and its effect is
measured as part of the solid discretisation error.

## Fixed shear modulus

The silicone's shear modulus is calibrated against the measured zero-flow
tip deflection (`./Allverify --study solid`): the flap alone under its net
buoyancy, damped to rest, on the nested hexahedral meshes S1-S4 (ratio 2 in
every direction), at two stabilisation scale factors. The calibrated value
converges at second order on the three finest levels:

| Level | μ_cal, sf 0.01 (kPa) | μ_cal, sf 0.001 (kPa) |
| --- | --- | --- |
| S1 | 44.81 | 50.07 |
| S2 | 58.42 | 60.68 |
| S3 | 62.58 | 63.23 |
| S4 | 63.61 | 63.79 |
| Richardson (S2-S4) | 63.95 (p = 2.01, GCI 0.42) | 63.94 (p = 2.21, GCI 0.19) |

Both scale factors give the same continuum limit, 63.9-64.0 kPa. μ* =
64.2 kPa was fixed from the S1-S3 limits (64.41 and 64.03 kPa) before S4
was available; it lies within the S2-S4 GCI of the limit. It is retained
because the steady coupled tip changes by only about -0.007 mm per kPa
(F1S2 with μ = 58.45 and 64.2 kPa), so the 0.25 kPa difference shifts the
tip by about 0.002 mm, far below the numerical uncertainties of the coupled
problem.

- The tutorial's 58.45 kPa is the S2-specific calibration (58.42 kPa at sf
  0.01): it compensates for the discretisation error of that mesh and is
  not used for verification.
- With ν = 0.45, the bending stiffness 2μ/(1 − ν) of a plate (2μ(1 + ν) of a
  beam) needs a higher μ than with ν ≈ 0.5, so 64 kPa is consistent with the
  ≈ 61 kPa reported by Hessenthaler et al. for a nearly incompressible model.
- The calibration uncertainty from the measured deflection is a model
  (validation) uncertainty. It is not part of the numerical uncertainty: μ*
  is a fixed parameter of the verified mathematical problem.

## Running

Source OpenFOAM, build solids4foam with PETSc, and run from this directory:

The `beta = 0.2` outlet needs the `backflowStabilisedOutletVelocity`
condition of PR #676; without it, use `outlet=zg` or `outlet=inletOutlet`.

```bash
./Allverify --study solid                      # zero-flow problem, S1-S4
./Allverify --study solid --levels S1,S2 --jobs 4
./Allverify --study meshes --threads 16        # build rF0-rF3, record quality
./Allverify --study run --spec rF1 --cores 24  # the campaign set-up on rF1
./Allverify --study run --spec rF0:outlet=zg --cores 24
./Allverify --study run --spec rF1 --restart --cores 24
./Allverify --study analyse                    # tables, orders, uncertainty
```

A coupled run is named by its spec: the fluid level (`rF0`-`rF3`; the
earlier family `F-1`-`F5` is kept for reference), the solid level `S1`-`S4`,
and the entries that differ from the campaign defaults
(`DEFAULT_SPEC` in the driver):

- fluid, solid (`rF1`, `S2`): the mesh levels.
- `dt`, `T` (0.002 s, 45 s): time step and end time; every run starts at
  t = 0.
- `outlet` (`b02`): `b02` (β = 0.2), `zg` (zero-gradient velocity) or
  `inletOutlet` (the tutorial's).
- `fsitest` (`strict`): `strict` requires both IQN residual measures (see the
  finding below); `default` is the solids4foam test.
- `tol` (1e-4): the coupling tolerance.
- `reuse` (2): IQN-ILS modes reused from previous steps (-1: the tutorial's
  5).
- `schemes` (`lust`): least-squares gradients, the full non-orthogonal
  correction and `Gauss LUST` convection; `lsq` keeps `linearUpwind`;
  `tutorial` is the limited set.
- `pimple`, `pimpletol` (20, 1e-5): PIMPLE outer correctors per coupling
  iteration, ended once the initial U and p residuals are below `pimpletol`.
- `sf`, `mu` (0.01, 64.2 kPa): the solid stabilisation and shear modulus.
- `fluidtol`, `solidtol` (`default`): `tight` linear-solver tolerances.
- `wall`: the wall cell size of level Fk (a diagnostic).
- `rep`: a replicate tag.

`--extend-from` continues a finished run to a later end time; `--restart`
continues an interrupted one from its latest written time.
Each run is a copy of the tutorial under `verification/work`; its evaluated
quantities are written to `verification/postProcessing/runs/<name>.json`.
Both directories are ignored by Git. `--reuse` re-evaluates a completed run
only if its settings, tutorial inputs and build match; `--evaluate-only`
re-evaluates an existing run directory.

The coupled runs are expensive: rF1 takes about a day on 24 cores, and rF3
needs a few hundred cores.

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

**Fluid** (cfMesh `cartesianMesh` from the feature-edge surface written by
`surfaceFeatureEdges -angle 45`). The campaign family rF0-rF3 takes the cell
sizes of the tutorial's `coarse` meshDict, scaled by 2^(−(k−1)/2) for level
rFk, so each level refines every region by a nominal √2, including the
octree root size `maxCellSize`. Each level adds:

- a 1.5 mm thick band of interface-size cells around the flap;
- a box at the interface size over the path of the flap tip;
- a shift of the octree, by one unused point appended to the surface, so
  that no flap face lies on an octree plane;
- `improveMeshQuality` after `cartesianMesh`.

Without feature edges, cartesianMesh chamfers the flap, and the transferred
load is then 5-150 % from the fluid force (#665). Without the shift, the flap
faces of some levels lie on octree planes, which leaves concave cells and a
0.02 mm bulge on the flap. Every level passes
`checkMesh -allGeometry -allTopology` with an exact flap; the driver records
the result as `strict_checkmesh_ok`. The meshes are not nested.

| Level | Cells | Interface faces | Cell sizes (mm) | Non-orth. |
| --- | --- | --- | --- | --- |
| rF0 | 122 435 | 3 784 | 0.714 / 1.43 / 2.86 | 39.5° |
| rF1 | 312 180 | 7 451 | 0.505 / 1.01 / 2.02 | 37.5° |
| rF2 | 832 552 | 14 654 | 0.357 / 0.714 / 1.43 | 37.3° |
| rF3 | 2 170 398 | 28 593 | 0.253 / 0.505 / 1.01 | 37.9° |

The cell sizes are those at the interface, in the flap sweep and at the
walls; Non-orth. is the maximum non-orthogonality. The effective linear
ratios (N_k/N_k−1)^(1/3) are 1.37, 1.39 and 1.38.

The earlier family F-1-F5 (`meshDict.coarse` scaled, without the band, box,
shift and quality pass) is kept in the driver for reference only. Its levels
F1 and F3 fail the strict check, and its results below are superseded.

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

With the strict FSI test the tip approaches its steady value geometrically (τ
≈ 4 s on F1S2 of the earlier family) after the start-up oscillation: 35-45 s
of simulated time are needed (the tutorial's 10 s leaves the F1S2 tip 0.25 mm
below its steady value). The steady QoIs are means over the last 4 s of runs
whose final drift is at most about 1e-4 mm/s; an unchanged restart of the
settled F1S2 stays constant to 5e-7 mm/s. The driver still reports its time
averages over the final 20 % of a run and an exponential-fit estimate of the
residual transient (a + b exp(−t/τ) on the smoothed history); with the
default convergence test these estimates were unreliable, because the
"transient" was an accumulating iterative error.

## Finding: the default FSI convergence test lets an iterative error accumulate

This finding was made on the earlier mesh family (F2S2, `inletOutlet`
outlet). It concerns the coupling test, not the mesh, and is why every run
uses `fsitest=strict`.

With the default IQN-ILS settings (`outerCorrTolerance 1e-4`), every coupled
run's tip kept rising long after the initial transient: +0.3-0.4 mm from 8 to
30 s, and still +0.003-0.005 mm/s at 60-80 s. Runs that differed only in rank
count or time step ended up 0.05-0.3 mm apart, so this drift exceeded every
spatial difference. Diagnosis:

1. Not the fluid alone: restarted at 30 s with the coupling off (the mesh is
   static), the flap force and a wake probe settle within ~2 s and stay
   constant to 0.3 % until 60 s.
2. Not the solid alone: the flap under buoyancy alone (production solver
   settings, velocity damping) is constant to 1e-15 mm.
3. The convergence test. A step is accepted when `min(r1, r2)` is below the
   tolerance. Here r1 is the interface residual relative to its largest value
   in the step. r2 is normalised by the TOTAL interface displacement, which is
   large for this flap (tip deflection ~17 mm). Late steps are therefore
   accepted as soon as r2 is small, whether or not the coupling iteration
   converged. The first iterations of every step are also under-relaxed
   (`relaxationFactor 0.05`). Per-step histories of the F2S2 run (0-80 s;
   iterations per step, the median residuals at the accepting iteration,
   and the share of steps accepted with r1 > 1e-2):

   | Window (s) | Iterations | r1 (median) | r1 > 1e-2 | r2 (median) |
   | --- | --- | --- | --- | --- |
   | 0-1 | 4.65 | 9.3e-4 | 4 % | 5.4e-5 |
   | 1-5 | 3.19 | 4.9e-2 | 87 % | 3.9e-5 |
   | 5-80 (each window) | 3.10-3.12 | 6.8e-2 | 90-91 % | 3.5e-5 |

   After the first second, every step is accepted with r1 above the
   tolerance: the coupling iteration reduces the step's interface mismatch
   by only ~15x. The fluid interface lags the solid (measured: ~0.03 mm in
   -y between the fluid interface and the deformed solid surface at 30, 55
   and 60 s), and the lag drives the slow drift.

4. With both measures required (`requireAllResidualMeasures yes`,
   `nOuterCorr 60`, `allowUnconvergedCoupling yes`; spec key
   `fsitest=strict`), from the F2S2 state at 80 s at the production time
   step: 10-11 iterations per step, median r1 5e-5 and r2 < 1e-8 at
   acceptance, and 0.0-0.4 % of steps above the tolerance. The tip drops
   within ~0.3 s from 17.031 to 17.0195 mm, then stays constant to
   +-0.0002 mm from 80 to 87.5 s. The same run with solid velocity damping
   (15 1/s) gives the same plateau (17.0193 mm).

Restarts from states built up under the default test keep part of that
history: a strict restart from another default-test state (the 48-rank
replicate at 60 s, 16.92 mm) settles at a different value (16.90 mm after
0.9 s, still being checked). The verification runs are therefore repeated
with the strict test from t = 0, which defines one trajectory and one end
state. The default-test runs are kept as the record of the finding.

The test is a solids4foam default, so other steady or slowly varying FSI
cases with large interface displacements may be affected in the same way.

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

## Results

The campaign on rF0-rF3 with the β = 0.2 outlet is in progress, and its
results will be added here when it is complete. The default series of
`--study analyse` are:

- fluid rF0-rF3 at S2;
- solid S1-S3 at rF1;
- Δt 0.004/0.002/0.001 s at rF0;
- the outlet (`zg`), coupling tolerance (1e-5), PIMPLE tolerance (1e-6),
  stabilisation (0.001) and a replicate, each paired with the rF0 run.

The earlier campaign (meshes F1-F3, `inletOutlet` outlet) is superseded and
its results have been removed from this directory. Those meshes chamfered the
flap (#665), and the outlet was not that of the benchmark definition.

Related issues: #552 (the default convergence test), #665 (the chamfered
flap) and #675 (the outlet condition).
