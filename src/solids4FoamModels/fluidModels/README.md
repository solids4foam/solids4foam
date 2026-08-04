---
sort: 5
---

# Fluid models

---

## Section Aims

This section aims to:

- Describe the fluid models available within solids4foam;
- Provide insight into the different solver types and their intended use cases.

---

## What is a Fluid Model?

In solids4foam, fluid models are classes which solve the governing equations in
a fluid domain, primarily to calculate the velocity and pressure fields. These
fluid models (base class `fluidModel`) can be interpreted as modular solvers,
packaged into a class structure so that they can be selected at runtime and
combined with solid models in fluid-solid interaction analyses.

Each fluid model manages its own mesh, time-stepping and boundary conditions.
The base class provides the velocity (`U`), pressure (`p`), volumetric flux
(`phi`), and related fields that coupling interfaces use to transfer loads and
displacements between the fluid and solid domains.

---

## Fluid Models Available in solids4foam

Here, we list the fluid models currently available in solids4foam, and briefly
summarise their characteristic features.

- **Incompressible single-phase models**

  - `pimpleFluid`: transient incompressible Newtonian fluid solver based on the
    PIMPLE (merged PISO-SIMPLE) algorithm; the most commonly used fluid model
    in solids4foam FSI analyses.
  - `newtonIcoFluid`: incompressible laminar solver using a Jacobian-free
    Newton-Krylov approach via the PETSc SNES interface; suitable for cases
    where a fully implicit nonlinear solve is preferred.

- **Incompressible multi-phase models**

  - `interFluid`: two-phase flow solver with volume-of-fluid (VOF) interface
    capturing, based on the `interFoam` family of solvers; supports dynamic
    meshes for FSI with a free surface.

- **Compressible models**

  - `sonicLiquidFluid`: transient compressible solver for weakly compressible
    liquids (sonic liquid), based on the `sonicLiquidFoam` solver.

- **Specialised models**

  - `pimpleOversetFluid`: a variant of `pimpleFluid` with dynamic mesh and
    overset (Chimera) mesh support; _currently only available with foam-extend_.
  - `oneWayFsiFluid`: a pseudo-fluid model that reads pre-computed velocity and
    pressure fields from a previously run fluid case; no fluid equations are
    solved. Used for one-way FSI where the fluid loading is known in advance.

---

## Notes

- All fluid models support dynamic meshes for use in FSI analyses, where the
  fluid domain deforms in response to solid motion.
- The `pimpleFluid` model is available for all supported OpenFOAM variants
  (OpenFOAM.com, OpenFOAM.org, and foam-extend), whereas some models are
  variant-specific.
- For one-way FSI (fluid drives solid, but solid does not affect fluid), use
  `oneWayFsiFluid` together with the `oneWayCoupling` fluid-solid interface.
