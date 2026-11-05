---
sort: 15
---

# Prestressed pulmonary artery fluid-solid interaction: `pulmonaryArteries`

Prepared by Conor Moore

## Tutorial Aims

- Demonstrate a patient-specific vascular fluid-solid interaction simulation
  using `solids4foam`.
- Demonstrate the prestress functionality in `solids4foam`: computing the
  prestress of an imaged, loaded vessel wall with the `updatePrestress`
  utility, and applying it in a fluid-solid interaction simulation.
- Demonstrate a finite-strain arterial wall using the `isotropicFungElastic`
  mechanical law.
- Demonstrate fluid and arterial-wall mesh generation using `cartesianMesh` and
  `extrudeMesh`.

## Case Overview

This case models pulsatile blood flow and compliant-wall deformation in the
main pulmonary artery and its proximal branches. The geometry is based on
model `h08` [1], and has one inlet and eleven outlets.

> **The parameter values in this case are representative only.** The wall
> thickness, wall stiffness, inflow waveform, outflow (Windkessel) parameters
> and their distribution between the outlets, and the prestress pressure are
> taken from the literature or chosen for illustration; none is specific to
> the subject whose geometry is used. The case demonstrates the numerical
> workflow and is not intended as a clinically validated prediction.

### Blood

Blood is represented as an incompressible, shear-thinning fluid with a density
of $$1060\,\mathrm{kg\,m^{-3}}$$, using the Carreau-Yasuda model (`BirdCarreau`
in `constant/fluid/transportProperties`) with a zero-shear viscosity of
$$0.16\,\mathrm{Pa\,s}$$, an infinite-shear viscosity of
$$3.5\times10^{-3}\,\mathrm{Pa\,s}$$, $$\lambda = 8.2\,\mathrm{s}$$,
$$n = 0.2128$$ and $$a = 0.64$$ [2]. The flow is laminar.

### Inflow and outflow

A representative pulmonary artery flow-rate waveform,
`constant/fluid/inletFlow/humanFlow`, is prescribed at the inlet with
`flowRateInletVelocity`. Its period is $$0.86\,\mathrm{s}$$, its mean
$$3.95\,\mathrm{L\,min^{-1}}$$ and its peak
$$272\,\mathrm{mL\,s^{-1}}$$ at $$t = 0.13\,\mathrm{s}$$ [1].

Each outlet uses a three-element Windkessel condition (`windkesselPressure` in
`0/fluid/p`) to represent the downstream pulmonary vasculature. The total
proximal resistance $$R_p = 3.3\times10^{6}\,\mathrm{Pa\,s\,m^{-3}}$$, distal
resistance $$R_d = 2\times10^{7}\,\mathrm{Pa\,s\,m^{-3}}$$ and compliance
$$C = 3.5\times10^{-8}\,\mathrm{m^3\,Pa^{-1}}$$ are based on figures for a
healthy person [3]. These values are split between the
outlets by a flow fraction proportional to the outlet area, $$f_i = r_i^2 /
\sum_j r_j^2$$ [4]:

```text
Rch_i = Rp/f_i      R_i = Rd/f_i      C_i = C*f_i
```

so that a larger outlet has a lower resistance and a higher compliance, and
the outlets in parallel add back up to the totals. The flow fraction of each
outlet is noted beside its entry in `0/fluid/p`.

### Arterial wall

The arterial wall has a uniform thickness of $$1.5\,\mathrm{mm}$$
[8], discretised with three cells through the thickness, and a density of
$$1200\,\mathrm{kg\,m^{-3}}$$. It is represented by the finite-strain
`isotropicFungElastic` law,

```text
W = 0.5*c*(exp(a*tr(Ebar & Ebar)) - 1) + U(J)
```

with $$c = 7.87\,\mathrm{kPa}$$, $$a = 7.05$$ and a bulk modulus
$$K = 1\,\mathrm{MPa}$$. The small-strain shear modulus is
$$\mu = c\,a/2 \approx 27.7\,\mathrm{kPa}$$.

The values of $$c$$ and $$a$$ were obtained by fitting this isotropic,
three-dimensional law to a planar biaxial (membrane) Fung fit of pulmonary
artery tissue, $$c = 15\,\mathrm{kPa}$$,
$$c_{\theta\theta} = c_{zz} = 6$$, $$c_{\theta z} = 1.8$$ [8]. The membrane
constants cannot be used directly in a three-dimensional law, because the
membrane fit does not include the through-thickness strain, so the fit is made
under the conditions of the test: plane stress and incompressibility, over
biaxial stretch ratios of 1:1, 1:0.75 and 0.75:1, up to a stretch of 1.3. The
isotropic law matches the membrane stresses to within about 1% RMS over this
range.

The bulk modulus is a penalty that keeps the wall nearly incompressible
(small-strain Poisson's ratio of about 0.486). Larger values make the wall
closer to incompressible, but slow the linear solver considerably: with
$$K = 100\,\mathrm{MPa}$$ each load step of the prestress run needed over 3000
linear iterations.

The inlet and outlet rings of the wall are fixed (`fixedDisplacement`), and the
outer wall is traction-free. The surrounding tissue is not represented; it
could be added with `solidSpringDashpot` on `outerWall`, as in the
`cerebralAneurysm` tutorial.

## Case Structure

The tutorial directory is the fluid-solid interaction case, with the
solid-only prestress case in `prestress/`:

```text
pulmonaryArteries/
├── Allrun, Allclean, Allmesh     run, clean and mesh the tutorial
├── 0/, constant/, system/        the fluid-solid interaction case
├── plotInletVsPressure.gp        inlet flow and probe pressure plot
├── regressionTest.sh             shortened run checked against references
└── prestress/
    ├── Allrun, Allclean          prestress passes
    ├── petscOptions.hypre        PETSc options for the prestress runs
    └── 0/, constant/, system/    the solid-only prestress case
```

The wall material is defined once, in `constant/solid/mechanicalProperties`;
`prestress/Allrun` copies it into the prestress case, so that the prestress
is computed with the same material that the fluid-solid interaction run uses.

## Mesh Generation

The vascular surface, `constant/triSurface/h08_capped.fms.gz`, is in
millimetres. It is stored compressed; OpenFOAM reads the `.gz` file directly,
so `system/meshDict` names it `h08_capped.fms`.

`Allmesh` creates the meshes:

1. `cartesianMesh` generates the fluid mesh from `system/meshDict`, with a
   maximum cell size of $$1.4\,\mathrm{mm}$$, one boundary layer, and
   $$0.7\,\mathrm{mm}$$ refinement around the two smallest outlets.
2. The arterial wall is extruded outwards from the fluid `wall` patch by
   `extrudeMesh` and `system/extrudeMeshDict`, giving three cells across the
   $$1.5\,\mathrm{mm}$$ thickness. The exposed inner face is `innerWall`, and
   `createPatch` with `system/createPatchDict.solid` renames the outer face
   `outerWall`.
3. `transformPoints` scales both meshes from millimetres to metres.
4. The solid mesh is copied to `prestress/constant/polyMesh`.

The fluid mesh has about 87,000 cells and the solid mesh about 69,000.

Minimal root-level `fvSchemes` and `fvSolution` dictionaries are included
because `extrudeMesh` and `createPatch` require them; the simulation uses the
region dictionaries under `system/fluid` and `system/solid`.

## Prestress

Medical images show the artery while it is loaded by the blood pressure, so
the imaged geometry is not stress-free. If it is treated as stress-free, the
first part of the simulation inflates it beyond the imaged shape. The
prestress, `Fpre`, is a deformation gradient field which, applied to the
imaged geometry, puts the wall in equilibrium with the imaging pressure at
zero displacement. The deformation gradient seen by the mechanical law is
then $$\mathbf{F}\,\mathbf{F}_{pre}$$.

`Fpre` is found with a modified updated Lagrangian formulation (MULF) [5]. The
prestress case, `prestress/`, applies a pressure of $$10\,\mathrm{mmHg}$$
($$1333\,\mathrm{Pa}$$) to `innerWall` in a steady-state solid-only
simulation, and `prestress/Allrun` repeats the following pass:

1. Run `solids4Foam` to equilibrium under the pressure.
2. Run `updatePrestress`, which composes the deformation gradient of the run
   with the `Fpre` that the run used, $$\mathbf{F}_{pre} \leftarrow
   \mathbf{F}\,\mathbf{F}_{pre}$$, writes it to `0/Fpre` and removes the
   later time directories, so that the next pass starts again from the imaged
   geometry with zero displacement.

The first pass runs with `prestress no` in `mechanicalProperties`, and ramps
the pressure up from zero (`constant/timeVsPressure.ramp`). The later passes
run with `prestress yes`, and apply the full pressure from the start
(`constant/timeVsPressure.full`) in a single load step (`deltaT` is changed
from 0.25 to 1 in `system/controlDict`): the wall then already carries the
prestress, which balances the full pressure, and ramping the pressure again
would first deflate the prestressed wall, which the solver does not converge
through. `prestress/Allrun` restores `deltaT 0.25` when it exits.

Each pass appends a line to `prestress/prestressHistory`, with the convergence
measures of the run that used the previous `Fpre`:

- `maxMag(D)` and `rmsMag(D)`: the maximum and volume-weighted RMS
  displacement, i.e. how far the loaded wall still moves from the imaged
  geometry. Both tend to zero as `Fpre` converges.
- `maxMag(F-I)`: the largest deviation of the deformation gradient from the
  identity. This includes rigid rotation, so it can remain large where the
  wall turns but hardly strains.

The tutorial runs a fixed 15 passes. The convergence is gradual, roughly 15
to 20% reduction in `maxMag(D)` per pass for this case. For another case, use
these measures to choose the number of passes, or a tolerance on `maxMag(D)`,
to suit the accuracy required. `prestress/Allrun` finally copies `Fpre` to
`0/solid/Fpre`.

The prestress pressure of $$10\,\mathrm{mmHg}$$ is representative; it should
be the pressure at the time of imaging. The fluid-solid interaction case
starts from the same pressure: `0/fluid/p` initialises the internal field
and the Windkessel outlets at $$10\,\mathrm{mmHg}$$, so the prestressed wall
starts in equilibrium.

## Numerical Approach

- **Fluid:** incompressible non-Newtonian laminar flow using `pimpleFluid`,
  with mesh motion by `velocityLaplacian` and `quadratic inverseDistance`
  diffusivity from the wall.
- **Solid:** total Lagrangian, finite-strain formulation
  (`nonLinearGeometryTotalLagrangianTotalDisplacement`) with the
  `isotropicFungElastic` law, solved with PETSc SNES using a Jacobian-free
  Newton-Krylov (JFNK) approach [7], with Rhie-Chow stabilisation. The
  fluid-solid interaction run preconditions the linear solves with block
  Jacobi and LU (`system/solid/fvSolution`); the prestress runs use HYPRE
  BoomerAMG (`prestress/petscOptions.hypre`).
- **Coupling:** partitioned fluid-solid interaction using fixed relaxation with
  a relaxation factor of 1.0 and a Robin pressure boundary condition
  (`elasticWallPressure` with the `secant` `hsModel`), which provides the
  stability required by the strongly coupled blood-wall system [6]. The
  outer-loop tolerance is $$10^{-5}$$, the Robin pressure and flux tolerances
  are $$10^{-2}$$, and up to 30 FSI correctors are permitted per time step.
- **Interface:** `wall` in the fluid region and `innerWall` in the solid
  region.
- **Duration:** five cardiac cycles, `endTime 4.3`, with a fixed time step of
  $$2.5\times10^{-4}\,\mathrm{s}$$ (17,200 time steps). The inlet waveform
  uses `outOfBounds repeat`, so it repeats for each cycle. The Windkessel
  outlets take about two cycles to settle from their initial pressure (their
  time constant is $$R_iC_i = R_dC = 0.7\,\mathrm{s}$$), so the later cycles
  are the periodic, representative ones.

## Expected Results

> **Note:** reference images of the velocity, pressure, wall shear stress and
> wall displacement, and plots of the probe pressure, inner-wall displacement
> and prestress convergence, will be added here once the full five-cycle run of
> the case, which is currently running, has finished.

The pressure at a probe in the main pulmonary artery
($$(0.04, -0.008, 0.065)\,\mathrm{m}$$) is written every time step to
`postProcessing/pressureProbe/fluid`. When `gnuplot` is installed, `Allrun`
plots it against the inlet flow rate with `plotInletVsPressure.gp`, giving
`inletFlow_vs_probePressure.png`: both against time, and the pressure-flow
loop. The wall shear stress is written at each write time. The range of the
inner-wall displacement components is written every time step to
`postProcessing/0/solidDisplacementsinnerWall.dat`, and the number of FSI
iterations per time step is recorded by `fsiConvergenceData`.

## Running the Case

From the tutorial directory, run

```bash
./Allclean
./Allrun parallel
```

`Allrun` runs the three stages in order:

```text
Allmesh -> prestress/Allrun -> solids4Foam (fluid-solid interaction)
```

The case is configured for 16 subdomains (`system/decomposeParDict`, its
`fluid` and `solid` counterparts, and `prestress/system/decomposeParDict`);
change `numberOfSubdomains` in all four to run on a different number of cores.
Omit `parallel` to run in serial. The stages can also be run separately: for
example, `prestress/Allrun` can be rerun with a different number of passes
after `Allmesh`, without remeshing.

The case requires PETSc and `cartesianMesh`; if either is unavailable,
`Allrun` exits without running. The mesh-generation workflow is not supported
on foam-extend.

As a guide to run times, on 16 cores the first prestress pass took about one
minute and each later pass about 15 s of solver time. Three heartbeats (2.58 s)
of the FSI simulation time took 12 hours on the same 16 cores.

## Regression Test

The case includes a `regressionTest.sh` script:

```bash
./regressionTest.sh
```

Because the full case takes about a day, the script copies the case inputs to
`pulmonaryArteries/regressionTests/main/` and runs a shortened version there;
the tutorial directory itself is untouched. The copy keeps the tutorial meshes,
since coarser meshes degrade on this geometry, but runs three prestress passes
in place of 15, stops the fluid-solid interaction run after five time steps
($$1.25\times10^{-3}\,\mathrm{s}$$), and runs on four processes in place of
16. On four cores it takes about nine minutes.

The following are then checked, against stored reference values within loose
tolerances where a value is compared:

- `maxMag(D)` in `prestress/prestressHistory` falls with each pass, and
  matches the reference for the first and last passes,
- every FSI time step converges within `nOuterCorr` correctors,
- the prestressed wall stays near the imaged geometry: its largest
  displacement over the five steps is well below that of the unprestressed
  wall under the same pressure,
- the largest `innerWall` displacement component from
  `postProcessing/0/solidDisplacementsinnerWall.dat`,
- the probe pressure from `postProcessing/pressureProbe/fluid`.

Pass `--check-only` to re-check an existing regression run without re-running
the case. If PETSc or `cartesianMesh` is unavailable, `Allrun` exits silently
and the script skips the checks rather than reporting a failure. Like
`cerebralAneurysm`, the test is run by hand rather than from
`tutorials/Alltest-regression`.

## References

1. Yevtushenko, P., Kuehne, T., Bruening, J., & Goubergrits, L. (2025). A
   Simulation-Based Comparison of Human, Porcine and Ovine Pulmonary Artery
   Hemodynamics. Evaluating the Suitability of Large Animal Models for
   Endopulmonary Device Evaluation from a Hemodynamics Point of
   View. Cardiovascular Engineering and Technology, 16(6),
   688–702. <https://doi.org/10.1007/s13239-025-00803-z>
2. Whitty, J. P. M., Wain, R. a. J., Fsadni, A., Simo, J., & Francis,
   J. (2016). Computational Non-Newtonian Hemodynamics of Small Vessels. Journal
   of Bioinformatics Computational and Systems Biology, 1(1), 103.
3. Lungu, A., Wild, J. M., Capener, D., Kiely, D. G., Swift, A. J., & Hose,
   D. R. (2014). MRI model-based non-invasive differential diagnosis in
   pulmonary hypertension. Journal of Biomechanics, 47(12),
   2941–2947. <https://doi.org/10.1016/j.jbiomech.2014.07.024>
4. Yang, W., Dong, M., Rabinovitch, M., Chan, F. P., Marsden, A. L., &
   Feinstein, J. A. (2019). Evolution of hemodynamic forces in the pulmonary
   tree with progressively worsening pulmonary arterial hypertension in
   pediatric patients. Biomechanics and Modeling in Mechanobiology, 18(3),
   779–796. <https://doi.org/10.1007/s10237-018-01114-0>
5. M. W. Gee, C. Förster, and W. A. Wall, "A computational strategy for
   prestressing patient-specific biomechanical problems under finite
   deformation," *International Journal for Numerical Methods in Biomedical
   Engineering*, 26(1), 52-72 (2010). <https://doi.org/10.1002/cnm.1236>.
6. Ž. Tuković, M. Bukač, P. Cardiff, H. Jasak, and A. Ivanković (2019).
   "Added Mass Partitioned Fluid–Structure Interaction Solver Based on a
   Robin Boundary Condition for Pressure." In: J. Nóbrega and H. Jasak (eds),
   *OpenFOAM®*. Springer, Cham. <https://doi.org/10.1007/978-3-319-60846-4_1>.
7. P. Cardiff, D. Armfield, Ž. Tuković, and I. Batistić, "A Jacobian-Free
   Newton-Krylov Method for Cell-Centred Finite Volume Solid Mechanics,"
   *International Journal for Numerical Methods in Engineering*, 127(3),
   e70268 (2026). <https://doi.org/10.1002/nme.70268>.
8. Azadani AN, Chitsaz S, Matthews PB, Jaussaud N, Leung J, Ge L, Tseng
   EE. Regional mechanical properties of human pulmonary root used for the Ross
   operation. J Heart Valve Dis. 2012 Jul;21(4):527-34. PMID: 22953683.
