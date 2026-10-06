# hessenthalerFsi validation studies

This directory contains opt-in studies that validate the `hessenthalerFsi`
tutorial against the measurements of Phase I (steady inflow) of the
Hessenthaler et al. (2017) FSI experiment. Unlike a verification study, the
reference is a measurement, with its own uncertainty, not an exact or
converged numerical solution. Nothing here is run by `tutorials/Alltest` or
`tutorials/Alltest-regression`.

## Running

Source OpenFOAM, build solids4foam with PETSc, and run from this directory:

```bash
cd tutorials/fluidSolidInteraction/hessenthalerFsi/validation
./Allvalidate --study calibration              # solid alone, several meshes
./Allvalidate --study calibration --jobs 8     # run the cases concurrently
./Allvalidate --study phaseI --cores 32        # coarse and medium fluid meshes
./Allvalidate --study phaseI --levels medium --cores 32
./Allvalidate --study phaseI --case DIR        # evaluate a finished run
./Allvalidate --quick                          # short smoke runs, no checks
./Allvalidate --reuse                          # re-evaluate completed runs
```

The driver needs Python 3.8 or newer; `matplotlib` is optional and is used for
the figures. Each run is a copy of the tutorial under `validation/work`, and
the results are written to `validation/postProcessing` as CSV files, Markdown
summaries and PNG figures. Both directories are ignored by Git. The copies
never include result time directories of the tutorial.

A run is checked before it is used. Its solver log must end with `End` and
contain no fatal error. Every monitored history must be finite and cover every
time step from the first to the end time; restart segments are merged.
`--reuse` re-evaluates a completed run only if the settings stored with it
match the request: the study, mesh, material, time step, end time and ranks,
a hash of the tutorial inputs, a hash of this driver, and a fingerprint of
the OpenFOAM version and the solids4foam build. `Allvalidate` returns zero
only if every acceptance check passes.

## Calibration study

The silicone keeps curing, so Hessenthaler et al. recommend calibrating the
solid stiffness to the zero-flow deflection of the flap, 29.50 mm at the tip
for Phase I, rather than using the uniaxial test. Each calibration run is the
tutorial's solid region alone, under its net buoyancy, with μ = 61 kPa. The
buoyancy is raised in three plateaus of 0.9, 1.05 and 1.2 times its Phase I
value, and the flap is relaxed with damping to its static deflection on each.
The static deflection of a neo-Hookean solid with a fixed Poisson's ratio
depends only on ρg/μ, so load factor f at 61 kPa is the static state at
μ = 61/f kPa. Interpolation to 29.50 mm gives the calibrated μ.

The study runs the tutorial's total Lagrangian solid (stabilisation 0.01)
on three meshes, the tutorial mesh with stabilisation 0.05 and 0.001 and with
ν = 0.49, the updated Lagrangian solid, the cubic high-order solid with the
compact Jacobian and `faceStencilExtraCells 60` on the coarse mesh, and an
8-rank run. With ν = 0.45 the calibrated μ is 58.4 kPa on the tutorial mesh
(12 x 8 x 65) and 61.7 kPa on the finest mesh (18 x 12 x 98); the coarse mesh
(6 x 4 x 33) does not reach 29.50 mm within the tested range. The updated
Lagrangian solid gives 60.9 kPa and the high-order solid 62.0 kPa. The
8-rank run reproduces the serial run. The results and their discussion are
in the tutorial README.

## Phase I study

The coupled case is run to 10 s, with the buoyancy and the inflow ramped in
over the first 0.5 s from the straight flap. It is compared with:

- the measured centreline of the flap at x = 0: the root-mean-square (RMS)
  and maximum differences in y at the 59 measured points, and the tip
  position;
- the three velocity components on the planes z = 10 and 30 mm. The computed
  velocity is averaged over each MRI voxel (1.302 x 1.302 x 6 mm, sampled at
  4 x 4 x 12 points), and a voxel is compared only if at least half of it
  lies in the computed fluid. As in Hessenthaler, Röhrle and Nordsletten
  (2017), Eqs. 14 and 15, d̄ and d∞ are the mean and the maximum Euclidean
  distance between the computed and the measured velocity vectors,
  normalised by the peak inflow velocity of 630 mm/s.

A run passes if the tip is within 1 mm of the measured 16.41 mm, the
centreline RMS difference is below 1 mm, and the tip moves by less than
0.1 mm over the last tenth of the run.

Both fluid meshes pass with the tutorial's total Lagrangian solid: the tip
is at 16.54 mm (coarse) and 16.66 mm (medium) against the measured 16.41 mm,
the centreline RMS difference is 0.23 and 0.28 mm, and d̄ = 0.050 and 0.048
and d∞ = 0.34 and 0.32. See the tutorial README for the figures and the other
solids.

## Reference data and quality

- `phaseI_centreline.csv`, `phaseI_velocity.csv` and `uniaxial_test.csv`: the
  measured data, from the CC0 data set of the experiment,
  [doi:10.6084/m9.figshare.4141836.v1](https://doi.org/10.6084/m9.figshare.4141836.v1).
  The flap position is uncertain to about one MRI voxel, 0.977 mm, and the
  velocity to about 5 % of the encoding velocity: 6 mm/s for vx and vy and
  40 mm/s for vz.
- `hessenthalerFsi_validation_references.json`: the geometry and material
  data, the polynomial fits of the zero-flow and Phase I centrelines
  (Hessenthaler et al. 2017, Figure 3), and published tip positions.
- `Lozovskiy2019_fig9_upper_surface.csv`: the computed Phase I upper surface
  of the flap of Lozovskiy et al. (2019), digitised from their Figure 9 to
  about ±0.2 mm.
- The CHeart tip positions of Hessenthaler, Röhrle and Nordsletten (2017) are
  read from the inset of their Figure 5 to about ±0.1 mm. No published result
  is tabulated.

An independent COMSOL solution of Phase I has been requested and will be added
when it is available.
