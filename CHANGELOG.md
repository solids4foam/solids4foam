# Changelog

This changelog highlights significant user-facing changes in each solids4foam
release. For complete commit-level details and contributor information, see the
[GitHub Releases](https://github.com/solids4foam/solids4foam/releases) page.

## [Unreleased]

### Added

- A documentation page for each `mechanicalConstitutiveLaw`, next to its
  class, giving its entries and defaults, its state, its tangents, whether it
  supports the mixed displacement-pressure formulation, how it differs from
  the removed legacy law, and the tutorials that use it, with an index in the
  material models page. The legacy laws' pages were removed with them.

- The `mechanicalConstitutiveLaw` history, such as the plastic strain, is
  written at each write time as volFields named after the variable
  (`epsilonPEq`, `epsilonP`, `sigmaY` and so on), so that it can be viewed in
  ParaView as the legacy laws' fields could. They cover every material that
  declares the variable and are zero in the others' cells, carry no
  dimensions, and are written, not read: the restart files remain the
  authoritative copy.

- Regression coverage of the `mechanicalConstitutiveLaw` framework in
  `squarePlate`, `cantilever2d`, `thermalCavity`, `curvedBeams` and `3dTube`,
  each held to the answer of the removed legacy model. `curvedBeams` and
  `3dTube` are the first coverage of the `impK` registry lookup that the
  contact penalty models, the cohesive zone models and `elasticWallPressure`
  all depend on, and `cantilever2d`'s `unsCoupled` case is the first
  regression coverage of `coupledUnsLinearGeometryLinearElastic`.
- `cantilever2d` regression cases for `vertexCentredLinearGeometry`, which no
  regression case ran before, on the PETSc SNES path and a short run of the
  explicit path, which no tutorial exercises.
- Added `tests/precice`, which runs solids4foam's preCICE coupling cases from
  the [preCICE tutorials](https://github.com/precice/tutorials) against the
  current source and checks them against stored reference values. The preCICE
  team's own system tests pin a released solids4foam image, so they cannot
  catch a regression introduced on a branch; this closes that gap. Run by the
  `preCICE coupling test` workflow for pull requests targeting `master`, for
  any pull request labelled `test-precice`, and on request.
- Added six hyperelastic tutorials ported from the `solid-benchmarks`
  repository, each with a README comparing it with its published reference and
  a `regressionTest.sh`: `cantileverVibration`, `cooksMembrane`,
  `compressedSpheres`, `twistingHemisphere`, `rubberSealing` and
  `shallowIroning` (all under `tutorials/solids/hyperelasticity`).
- Added opt-in mesh-convergence verification studies, migrated from
  `solid-benchmarks`, to the linear-elastic and elastoplastic Cook's membrane,
  `sphericalCavity`, `ellipticPlate`, `narrowTmember`, `plateHole`,
  `cantilever2d`, `curvedCantilever`, `curvedBeams`, `idealisedVentricle` and
  `cavityFlexibleBottom` tutorials. Each lives in the tutorial's
  `verification/` directory, is run with `./Allverify`, works on copies under
  the ignored `verification/work/` directory, and is not run by
  `tutorials/Alltest` or `tutorials/Alltest-regression`.

### Changed

- `OgdenElastic` assembles its Cauchy stress along the principal directions
  of the left Cauchy-Green tensor, the spatial ones. It used the right
  tensor's, the material ones, so under a stretch with a finite rotation the
  stress was rotated back by the rotation and the law was not objective. No
  tutorial selects it; results of a case with large rotations change.
- `linearElasticMohrCoulombPlastic` refuses a friction angle below 1e-3
  degrees in magnitude again, as the legacy law did: the apex of the
  Mohr-Coulomb surface moves to infinity in that limit and the law's derived
  parameters divide by zero. Its deviatoric scalar tangent is `(4/3)*mu`, as
  for `linearElastic`, and so is `OgdenElastic`'s.
- `linearElasticMisesPlastic` normalises its Newton residual for nonlinear
  hardening by the largest strain over every rank, as
  `neoHookeanElasticMisesPlastic` already did, so the plastic update no
  longer depends on the decomposition. No tutorial uses the nonlinear
  hardening branch.
- A law dictionary that still sets an entry only the removed legacy laws read,
  such as `pressureDisplacement`, `impKcoeff`, `tangentEps`, `regionName` or
  `pressureFieldRegion`, gives a warning naming it; the entries were ignored
  silently.

- The coupling inputs a `mechanicalConstitutiveLaw` reads, such as the
  active tension of `electroMechanicalLaw` or the temperature of
  `thermoMechanicalLaw`, are gathered on every evaluation path, where before
  only the small-strain cell paths gathered them and the others stopped. So
  `electroMechanicalLaw` with `activeTensionFromField yes` now runs under the
  finite-strain solid models, which are the only ones it supports. With the
  switch off, the law uses its ramped constant and the run warns once if a
  field of that name is registered: unlike the removed legacy law, a
  registered `Ta` is not picked up unasked.
- The restart files of the `mechanicalConstitutiveLaw` history are named
  `<material>_<topology>_<variable>`, where the development branch had
  colons, which Windows file systems, archive tools and CI artifact uploads
  refuse.

- `vertexCentredLinearGeometry` takes its whole constitutive response from the
  `mechanicalConstitutiveLaw` framework: the residual stress at the dual mesh
  faces, on the PETSc SNES and explicit paths alike, as well as the Jacobian
  tangent, the cell stress and the stiffness behind the default
  `fixedDofScale` and the explicit time step. More than one material is now
  supported; the removed legacy path stopped with "SubMesh not found when
  looking for a field in the base mesh".
- `kirchhoffPlate`, `coupledUnsLinearGeometryLinearElastic` and `thermalSolid`
  take their material from the `mechanicalConstitutiveLaw` framework. The
  first two read material constants from a single isotropic linear elastic law
  rather than asking it for a stress - a plate has one bending stiffness, and
  the block-coupled solver assembles one set of Lame constants - so both keep
  that restriction.
- On foam-extend, a run with more than one material now runs, in `linearGeometryTotalDisplacement`,
  `nonLinearGeometryTotalLagrangianTotalDisplacement`,
  `nonLinearGeometryUpdatedLagrangian`, `poroLinearGeometry` and
  `thermalLinearGeometry`. The displacement is interpolated to the points the
  way the framework does it on OpenFOAM.com and OpenFOAM.org: each cell around
  a point extrapolates its value with its own gradient, which the
  material-aware `leastSquaresS4f` scheme keeps to one material, and the
  extrapolated values are averaged with inverse-distance weights. foam-extend's
  least squares interpolation, which a stencil straddling a material interface
  smears, is not used for these runs, and nor are the legacy per-material
  sub-meshes. On `layeredPipe` the point displacement is within 1.5e-3 of the
  legacy model's, where the least squares interpolation is 4.9e-3 away.
  Single-material runs on foam-extend keep the least squares interpolation, as
  before.
- **Breaking:** the pore pressure field of `poroLinearGeometry` and the default
  pore pressure field name of `poroMechanicalLaw` are now `porePressure` rather
  than `p`. An existing case carrying `0/p` fails at construction with a
  `MUST_READ` error and must rename the field; the bundled tutorials are
  already migrated. The rename is needed because a solid model solving the
  mixed displacement-pressure formulation has its own `p`, which is a different
  quantity, and the two collided in the same registry. The old name is not
  accepted as a fallback: the two fields are not interchangeable, so silently
  reading one where the other was meant would be worse than a clear failure.
  This applies whether or not the `mechanicalConstitutiveLaw` framework is in
  use.
- `linearGeometryTotalDisplacement` with `solvePressure yes` now uses an
  implicit stiffness of `(4/3)*mu` rather than `2*mu`. The mixed
  displacement-pressure formulation solves `div(dev(sigma))`, whose scalar
  Laplacian surrogate is `mu*lap(D) + (1/3)*mu*grad(div(D))`, so `(4/3)*mu` is
  the consistent coefficient and `2*mu` had no derivation behind it. The
  implicit operator is the iteration path rather than the equation being
  solved, so converged results are unchanged within tolerance, but the
  iteration count and the path taken to get there move for every case that
  selects this option.

### Removed

- **Breaking:** the legacy `mechanicalModel` is removed - `mechanicalModel`,
  `dualMechanicalModel`, `solidSubMeshes`, `mechanicalLaw` and every legacy
  law - together with the `useMechanicalConstitutiveLawManager` switch that
  chose it. The `mechanicalConstitutiveLaw` framework is the only mechanical
  model, reading the same `constant/mechanicalProperties`, and results are
  those the switch's `yes` gave. A case that still sets
  `useMechanicalConstitutiveLawManager no;` stops with an error saying so,
  rather than silently changing its results; `yes` is accepted with a note
  that the entry can be removed. Entries only the legacy model read, such as
  a law's `pressureDisplacement` and `writeSubMeshes`, are ignored. The
  sources of the withdrawn `unsLinearGeometry`,
  `unsNonLinearGeometryTotalLagrangian`,
  `unsNonLinearGeometryUpdatedLagrangian`, `weakThermalLinearGeometry` and
  `vertexCentredNonLinTotalLagGeometry` solid models, which depended on the
  legacy model, are deleted.
- The material residual is removed from every solid model, with the `matRes`
  column of the residual output and of `residual.dat`, which now hold one
  column fewer. The legacy laws under-relaxed their plastic strain increment,
  so the constitutive update lagged the displacement between outer iterations
  and was tested against `materialTolerance`; a `mechanicalConstitutiveLaw` is
  a pure function of the kinematics and the old-time state, so there was
  nothing left to measure and the column was always zero. `materialTolerance`
  is ignored with a note that it can be removed, and the tutorials no longer
  set it.
- **Breaking:** `nonLinearGeometryUpdatedLagrangian` no longer supports
  `solvePressure`. The mixed formulation replaces a law's volumetric response
  with the solved pressure; this model takes `dev()` of the total stress to
  recover the rest, which is right only where everything the pressure must not
  replace is trace free - true of an isotropic hyperelastic law, false of an
  active tension along a fibre. It stops at construction rather than solving a
  different material. Use `nonLinearGeometryTotalLagrangianTotalDisplacement`,
  which asks the law for the two parts separately.
- **Breaking:** ten mechanical laws are removed: `diffusionElastic`,
  `linearElasticCt`, `linearElasticFromFile`, `orthotropicLinearElastic`,
  `GentElastic`, `StVenantKirchhoffOrthotropicElastic`, `YeohElastic`,
  `diffusionHyperElastic`, `isotropicFungElastic` and `viscoNeoHookeanElastic`.
  A case selecting any of these now fails to construct its mechanical model.
  None of the ten was ported to the `mechanicalConstitutiveLaw` framework, and
  none is selected by any tutorial, so none has ever had regression coverage on
  either path; three carry open correctness issues that nothing would have
  caught. Keeping unported, untested laws alive is what the framework exists to
  get away from, and `linearElasticCt` was in neither build list already, so it
  could not be selected at run time in any case. Each may return, but only
  together with a tutorial and a regression case that pins its numbers; see
  #457, and #458 for the two orthotropic laws in particular.
- **Breaking:** the `abaqusUMATs` library, its `plateHoleTotalDispUMAT`
  tutorial and the `S4F_USE_GFORTRAN` build hook are removed. It was a proof of
  concept that never covered the full Abaqus UMAT interface. The
  `abaqusMeshToFoam` and `foamMeshToAbaqus` mesh utilities are unaffected.
- **Breaking:** `weakThermalLinearGeometry` is no longer compiled, on any fork,
  following the precedent set for `vertexCentredNonLinTotalLagGeometry`. No
  tutorial selects it, so nothing exercises it, and it derives from
  `linearGeometryTotalDisplacement` and so inherits the
  `mechanicalConstitutiveLaw` framework path without anything having tested
  that it works there. `thermalLinearGeometry` may in any case cover every
  problem it was written for. Both build lists carry the reason in a comment,
  and the solver's sources are untouched, so restoring it is a matter of
  uncommenting the entry in both files - together with a tutorial. See #459,
  which also asks whether `thermalLinearGeometry` supersedes it entirely.
- **Breaking:** the public virtual `solidModel::newDeltaT()`, which forwarded to
  `mechanicalModel::newDeltaT()`, is removed. Its purpose was to let a
  constitutive law ask for a smaller time step, but nothing in solids4foam
  called it: the `solids4Foam` solver never consulted it, so no law could
  actually influence the time step through it. The `mechanicalConstitutiveLaw`
  framework provides no equivalent, and adding one is left until there is a
  caller to justify its shape. An external driver calling
  `solidModel::newDeltaT()` no longer compiles, and adaptive time stepping
  driven by material state will need a new interface rather than this one.
- `vertexCentredNonLinTotalLagGeometry` is no longer compiled, on any fork. The
  solver does not currently run: no tutorial selected it, so nothing exercised
  it, an attempt to give it one failed inside the PETSc solve, and its
  `useGeometricStiffness`, `compactImplicitStencil` and `tangentEps` entries
  have no defaults and are set by no case. Selecting that runtime type now
  fails to construct the solid model rather than failing inside it. Both build
  lists carry the reason in a comment, and the solver's sources and `README.md`
  are untouched. It is withdrawn rather than removed: the intention is to
  revisit it once the mechanical constitutive law framework has landed, at
  which point restoring it is a matter of uncommenting the entry in both
  files.
- Removed `tutorials/fluidSolidInteraction-preCICE`, which held standalone
  `3dTube` and `flexibleOversetCylinder` preCICE cases that were not covered by
  any test. solids4foam's preCICE cases are now maintained upstream in the
  preCICE tutorials and tested by `tests/precice`. The removed cases remain
  available as an archive from the solids4foam website.

## [v2.4] - 2026-08-24

### Added in v2.4

- Added support for OpenFOAM-v2606, extending compatibility to OpenFOAM-v2312
  through OpenFOAM-v2606, alongside the existing OpenFOAM-9 and
  foam-extend-4.1 support.
- Added a configurable stabilisation framework, including combined
  stabilisation models, volumetric-strain-rate stabilisation, and PETSc SNES
  Jacobian support for JST and even-order schemes.
- Added high-order finite-volume solid mechanics, with high-order
  configurations for `cantilever2d`, Cook's membrane, `pressurisedCylinder`,
  `plateHole`, and `sphericalCavity`
  (<https://doi.org/10.1016/j.jcp.2026.115056>).
- Added a block-coupled incompressible solid formulation
  (<https://doi.org/10.3390/app152312660>).
- Added `decomposeParMonolithic` for consistently decomposing coupled meshes.
- Added adaptive time stepping to `newtonIcoFluid`.
- Added the `solidPressureMinMax` function object and utilities for detecting
  PETSc build state and summarising solver logs.
- Added an optional unity build for `libsolids4FoamModels`, reducing
  compilation time.
- Added validation of the PETSc configuration at build time, continuous
  integration across supported OpenFOAM versions, and a v2512 Docker release
  image.
- Added further tutorial cases and per-tutorial README.md files, including the
  poroelasticity, elastoplasticity, `ellipticPlate`, and one-way cavity cases,
  plus expanded website documentation.

### Changed in v2.4

- Tutorial cases are now stored in the OpenFOAM.com format and converted to
  foam-extend or OpenFOAM.org at run time, rather than the reverse. All three
  variants continue to work; see the upgrade notes below.
- Incremental (updated Lagrangian) solid models now stop with a clear message
  when an initial `D` field is present outside a restart; see the upgrade notes
  below.
- Improved IQN-ILS coupling and simplified FSI region selection.
- Extended `electroMechanicalLaw` with field-based active tension and
  independent fibre tensors.
- Made point-patch enforcement optional in the `fixedDisplacement` boundary
  condition, and made the `fluidModel` fluid-property accessors public.
- Unsupported foam-extend `cyclicGgi` patches now fail with a clear error
  rather than segfaulting during volume-to-point interpolation.
- Expanded regression coverage for solid, contact, FSI, and least-squares
  gradient cases, and added a `Test-leastSquaresS4fGrad` test application.

### Fixed in v2.4

- Fixed a potential segmentation fault when PETSc SNES backs up its solution
  before the first solve.
- Fixed PETSc Jacobian preallocation by creating the Jacobian as a block AIJ
  matrix, resolving segmentation faults with recent PETSc versions.
- Fixed old-time deformation gradient storage in `StVenantKirchhoffElastic`.
- Fixed least-squares gradient evaluation at boundary faces, and ported
  `enhancedVolPointInterpolation` to OpenFOAM-9, correcting rigid-rotation
  cases solved with the updated Lagrangian approach.
- Fixed the sign convention in the even-order Laplacian for the `m == 0` case.
- Fixed boundary condition enforcement in the `diffusionElastic` mechanical law
  and the tabulated acceleration fields in `sonicLiquidFluid`.
- Fixed the `RBFMeshMotionSolver` setup for fluid-solid interaction, and
  documented a verified configuration.
- Fixed `tmp<fvMatrix>` assembly on OpenFOAM.org.
- Fixed the restart behaviour of the `elasticWallVelocity` boundary condition,
  which previously discarded its face-centre history and so restarted with a
  zero interface velocity. The history is now written to, and read from, the
  time directories, and is mapped by `decomposePar` and `reconstructPar`.
- Fixed the restart behaviour of `fluxCorrectedVelocity` boundaries in
  `fluidModel`. As the condition is derived from `zeroGradient`, it does not
  read the `value` entry that it writes, so the normal component of the
  velocity was replaced by a zero-gradient extrapolation on restart. The
  current boundary values of `U` are now re-derived from `phi`, and those of
  its old-time levels, for which `phi` is not written, are read back from the
  time directories.

### Removed in v2.4

- Removed tutorial-specific analytical solution function objects from
  `src/solids4FoamModels`; these now live in case-local libraries within the
  corresponding tutorials.

### Upgrade notes for v2.4

- **Tutorial case format**: tutorial cases are stored in the OpenFOAM.com
  format. Cases copied from the repository run directly on OpenFOAM.com, and
  the `Allrun` scripts convert them for foam-extend and OpenFOAM.org. Local
  copies of older cases are unaffected, but any scripts that assumed the
  foam-extend layout (for example `constant/polyMesh/blockMeshDict`) should be
  updated.
- **Initial `D` field with incremental solid models**: incremental updated
  Lagrangian solid models solve for `DD`, and an initial `D` field alongside it
  can make boundary conditions pick up an inconsistent displacement history.
  Such cases now fail at startup unless the run is an explicit restart. Remove
  `0/D` from affected cases.
- **Tutorial analytical solutions**: if a case or user library referenced an
  analytical solution function object from `libsolids4FoamModels`, link against
  the case-local library in the corresponding tutorial instead.
- **Interface-to-interface mapping and the undeformed mesh**: all
  interface-to-interface mappings (`AMI`, `GGI`, `RBF` and `directMap`) now
  build the interface correspondence from the undeformed mesh points in the
  `constant` instance, rather than from the current interface. This keeps the
  correspondence independent of when it is first constructed, now that the
  cached interface geometry follows the mesh motion and the solid deformation.
  Cases whose mesh changes topology during the run, or whose points are not
  available in the `constant` instance, now fail with an explicit error instead
  of silently building the correspondence from whatever configuration the
  interface happened to be in. Previously only `AMI` read the `constant`
  points, so this affects a wider set of cases than before.

### Related to v2.4

- [`beamFoam`](https://github.com/solids4foam/beamFoam) was released alongside
  the v2.4 update (<https://doi.org/10.51560/ofj.v5.170>). It is developed and
  maintained in a separate repository and is not part of solids4foam itself.

## [v2.3] - 2026-02-04

### Added in v2.3

- Added a Newton-Krylov solid solver, offering improved robustness and
  efficiency compared to traditional segregated solvers.
- Added a PETSc SNES interface for nonlinear solution procedures.
- Added Robin-Neumann coupling support for two-phase flows via the
  `interFluid` fluid model.
- Added the `cavityFlexibleBottom` fluid-solid interaction tutorial and new
  regression coverage for major FSI cases.
- Added the option to install solids4foam via the OpenFOAM package manager
  styro.
- Added README.md files across the tutorials and expanded website
  documentation.

### Changed in v2.3

- Extended compatibility to OpenFOAM-v2312 through OpenFOAM-v2512.
- Improved nonlinear vertex-centred solvers and several FSI tutorials.
- Made optional OpenFOAM source-file fixes opt-in by default.

## [v2.2] - 2025-04-02

### Added in v2.2

- Added support for OpenFOAM-v2406.
- Updated the preCICE tutorial configurations for preCICE v3.
- Added an analytical solution for the square-plate tutorial.
- Added a correction procedure to `perturbMeshPoints` for avoiding poor-quality
  cells.

### Changed in v2.2

- Improved least-squares interpolation at symmetry boundaries and corrected
  the plane-stress bulk modulus.
- Added citation, licensing, contribution, and pull-request metadata, and
  introduced automated Markdown checks.

## [v2.1] - 2024-06-22

### Added in v2.1

- Added conjugate heat-transfer and thermo-fluid-solid interaction support.
- Added segment-to-segment contact to the `solidContact` boundary condition.
- Added unified linear and nonlinear vertex-centred solid models with
  block-coupled, segregated, and explicit solution algorithms.
- Added the general `electroMechanicalLaw` and `poroMechLaw` wrappers and the
  `GuccioneElastic` mechanical law.
- Added further documented tutorials, including Cook's membrane, contact,
  curved-beam, and thermo-fluid-solid interaction cases.

### Changed in v2.1

- Extended support through OpenFOAM-v2312 while retaining OpenFOAM.org and
  foam-extend compatibility.
- Enabled Robin-Neumann fluid-solid coupling across OpenFOAM variants.
- Added the option to write case dictionaries with their default values.

## [v2.0] - 2022-12-20

### Added in v2.0

- Added support for coupling solids4foam to preCICE, including dedicated FSI
  tutorials.
- Added Mooney-Rivlin, Yeoh, isotropic Fung, and Ogden hyperelastic laws.
- Added PETSc-based vertex-centred solid mechanics capabilities.
- Added automated GitHub build testing and published Docker configurations for
  supported OpenFOAM variants.

### Changed in v2.0

- Moved project development from Bitbucket to GitHub and introduced the
  solids4foam website.
- Expanded OpenFOAM.com and OpenFOAM.org support for major features, including
  multi-material solids and solid contact.
- Redesigned the build and tutorial-test scripts for consistent behaviour
  across OpenFOAM variants.

## [v2.0-alpha] - 2022-09-16

### Added in v2.0-alpha

- Published the first preview of the cross-version compatibility work for
  v2.0.
- Added initial GitHub Actions and Docker build coverage for OpenFOAM-v2012,
  OpenFOAM-9, and foam-extend-4.1.

## [v1.1] - 2022-01-26

### Added in v1.1

- Added the `unsIcoFluid` fluid model and the `poroAnisotropicBiotElastic`
  solid model.
- Added Kirchhoff plate capabilities and the `solidTorque` function object.
- Added pressure and traction input from point-cloud data to the solid-traction
  boundary condition.
- Added the `abaqusMeshToFoam` utility for hexahedral Abaqus meshes.

### Changed in v1.1

- Added OpenFOAM-v1912 and Clang 12 compatibility.
- Improved parallel least-squares volume-to-point interpolation and Aitken FSI
  relaxation controls.
- Made Fortran-based Abaqus UMAT support optional.

## [v1.0] - 2021-07-22

- Published the initial solids4foam release for finite-volume solid mechanics
  and fluid-solid interaction simulations.
- Included linear, nonlinear, thermal, poromechanical, contact, and coupled
  solid models with a broad tutorial collection.
- Supported foam-extend-4.0 and foam-extend-4.1, with initial support for
  OpenFOAM-7 and OpenFOAM-v1812.

## [v0.1] - 2016-07-24

- Added the initial solids4foam codebase (then developed on Bitbucket).
- Introduced three working `solidFoam` solid-mechanics solvers.
- Established the foundation for further solids4foam development.

[v2.4]: https://github.com/solids4foam/solids4foam/releases/tag/v2.4
[v2.3]: https://github.com/solids4foam/solids4foam/releases/tag/v2.3
[v2.2]: https://github.com/solids4foam/solids4foam/releases/tag/v2.2
[v2.1]: https://github.com/solids4foam/solids4foam/releases/tag/v2.1
[v2.0]: https://github.com/solids4foam/solids4foam/releases/tag/v2.0
[v2.0-alpha]: https://github.com/solids4foam/solids4foam/releases/tag/v2.0-alpha
[v1.1]: https://github.com/solids4foam/solids4foam/releases/tag/v1.1
[v1.0]: https://github.com/solids4foam/solids4foam/releases/tag/v1.0
[v0.1]: https://github.com/solids4foam/solids4foam/commit/ccc2e752d6620c18f6cd42b38dddfe35c8a168b4
