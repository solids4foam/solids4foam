# Changelog

This changelog highlights significant user-facing changes in each solids4foam
release. For complete commit-level details and contributor information, see the
[GitHub Releases](https://github.com/solids4foam/solids4foam/releases) page.

## [v2.4] - Unreleased

### Added in v2.4

- Added a configurable stabilisation framework, including combined
  stabilisation models, volumetric-strain-rate stabilisation, and PETSc SNES
  Jacobian support for JST and even-order schemes.
- Added a block-coupled incompressible solid formulation.
- Added `decomposeParMonolithic` for consistently decomposing coupled meshes.
- Added adaptive time stepping to `newtonIcoFluid`.
- Added the `solidPressureMinMax` function object and utilities for detecting
  PETSc build state and summarising solver logs.
-Added new tutorial cases and per-tutorial README.md files, plus expanded website documentation.

### Changed in v2.4

- Expanded high-order solid mechanics and provided high-order configurations
  for `cantilever2d`, Cook's membrane, `pressurisedCylinder`, `plateHole`, and
  `sphericalCavity`.
- Improved IQN-ILS coupling and simplified FSI region selection.
- Extended `electroMechanicalLaw` with field-based active tension and
  independent fibre tensors.
- Expanded regression coverage for solid, contact, FSI, and least-squares
  gradient cases.

## [v2.3] - 2026-02-04
### Added in v2.3
- Added a Newton-Krylov solid solver, offering improved robustness and
  efficiency compared to traditional segregated solvers.
- Added high-order finite-volume capabilities for solid mechanics.
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

[v2.4]: https://github.com/solids4foam/solids4foam/compare/v2.3...development
[v2.3]: https://github.com/solids4foam/solids4foam/releases/tag/v2.3
[v2.2]: https://github.com/solids4foam/solids4foam/releases/tag/v2.2
[v2.1]: https://github.com/solids4foam/solids4foam/releases/tag/v2.1
[v2.0]: https://github.com/solids4foam/solids4foam/releases/tag/v2.0
[v2.0-alpha]: https://github.com/solids4foam/solids4foam/releases/tag/v2.0-alpha
[v1.1]: https://github.com/solids4foam/solids4foam/releases/tag/v1.1
[v1.0]: https://github.com/solids4foam/solids4foam/releases/tag/v1.0
