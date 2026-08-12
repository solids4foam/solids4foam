# Changelog

This changelog highlights significant user-facing changes in each solids4foam
release. For complete commit-level details and contributor information, see the
[GitHub Releases](https://github.com/solids4foam/solids4foam/releases) page.

## [v2.4] - Unreleased

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
  interface-to-interface mappings (`ami`, `ggi`, `RBF` and `directMap`) now
  build the interface correspondence from the undeformed mesh points in the
  `constant` instance, rather than from the current interface. This keeps the
  correspondence independent of when it is first constructed, now that the
  cached interface geometry follows the mesh motion and the solid deformation.
  Cases whose mesh changes topology during the run, or whose points are not
  available in the `constant` instance, now fail with an explicit error instead
  of silently building the correspondence from whatever configuration the
  interface happened to be in. Previously only `ami` read the `constant`
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

[v2.4]: https://github.com/solids4foam/solids4foam/compare/v2.3...development
[v2.3]: https://github.com/solids4foam/solids4foam/releases/tag/v2.3
[v2.2]: https://github.com/solids4foam/solids4foam/releases/tag/v2.2
[v2.1]: https://github.com/solids4foam/solids4foam/releases/tag/v2.1
[v2.0]: https://github.com/solids4foam/solids4foam/releases/tag/v2.0
[v2.0-alpha]: https://github.com/solids4foam/solids4foam/releases/tag/v2.0-alpha
[v1.1]: https://github.com/solids4foam/solids4foam/releases/tag/v1.1
[v1.0]: https://github.com/solids4foam/solids4foam/releases/tag/v1.0
[v0.1]: https://github.com/solids4foam/solids4foam/commit/ccc2e752d6620c18f6cd42b38dddfe35c8a168b4
