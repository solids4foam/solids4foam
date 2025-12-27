# solids4foam Tutorials

This directory contains tutorial cases for the **solids4foam** toolbox, covering
solid mechanics, fluid mechanics, and fluid–solid interaction problems. The
tutorials serve two main purposes:

1. **User examples** demonstrating how to set up and run typical solids4foam
   simulations.
2. **Automated tests** used to verify correct compilation and behaviour of
   solids4foam across different OpenFOAM versions and build configurations.

The tutorials are organised into the following subdirectories:

- `solids/` – solid mechanics cases
- `fluids/` – fluid mechanics cases
- `fluidSolidInteraction/` – monolithic FSI cases
- `fluidSolidInteraction-preCICE/` – partitioned FSI cases using preCICE
- `thermoFluidSolidInteraction/` – thermo-mechanical and coupled problems

---

## Prerequisites

All tutorial scripts assume that:

- OpenFOAM has already been installed and sourced (i.e. the OpenFOAM environment
  is active).
- solids4foam has been compiled successfully for the active OpenFOAM version.
- A standard OpenFOAM directory layout is being used.

No additional environment setup is performed by the tutorial scripts themselves.

---

## Tutorial Control Scripts

Several top-level helper scripts are provided to manage and test the tutorials.

### `Allrun`

Runs the tutorial cases.
Depending on the subdirectory, this may involve mesh generation and solver
execution. Individual tutorial directories may also provide their own `Allrun`
scripts with additional options.

This script is primarily intended for **interactive use** by users.

---

### `Allclean`

Cleans the tutorial directories by removing generated meshes, time directories,
and log files, returning cases to a clean state.

---

### `Alltest` (Smoke Tests)

Runs a **fast smoke test** across the tutorials.
The purpose of this script is to verify that:

- The tutorials start correctly.
- The solvers run and complete at least one iteration or time step.

No attempt is made to validate the numerical correctness of the results.

This script is:

- **Fast**
- Used routinely in **GitHub Actions CI**
- Intended to catch build, linking, and obvious runtime failures

---

### `Alltest-regression` (Regression Tests)

Runs a small number of **physics-based regression tests** on selected tutorial
cases. These tests go beyond simply checking that a case runs, and instead verify
that key numerical quantities remain within expected bounds.

For example, some regression tests:

- Extract scalar values from solver logs (e.g. stresses or residuals)
- Check that these values are below or close to known reference thresholds

Regression tests are:

- **More expensive** than smoke tests
- Intended primarily for **developers**
- Useful when modifying solvers, material models, or coupling algorithms

Not all tutorials are covered by regression tests; coverage is expected to grow
incrementally over time.

---

## Continuous Integration (CI)

In GitHub Actions, the tutorials are used in two stages:

1. **Smoke testing** via `Alltest` to ensure that solids4foam builds and runs
   correctly across supported OpenFOAM versions and configurations.
2. **Regression testing** via `Alltest-regression` for selected cases, providing
   additional confidence that numerical behaviour remains correct.

Contributors are encouraged to run the relevant tutorial tests locally when
making significant changes to solver logic or numerical methods.

---

## Notes

- Some tutorials may be computationally more expensive than others.
- Regression tests are intentionally limited in scope to keep CI runtimes
  reasonable.
- Individual tutorial directories may include additional documentation describing
  the physical problem being solved.

For detailed descriptions of specific tutorial cases, see the documentation in
the corresponding subdirectories or on
<https://www.solids4foam.com>.
