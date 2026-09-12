---
sort: 5
---

# Idealised Ventricle Inflation: `idealisedVentricle`

---

## Tutorial Aims

- Demonstrate the simulation of large-strain incompressible hyperelastic
  inflation of an idealised left ventricle.
- Show two solution approaches available in solids4foam:
  - PETSc SNES (default), as used in [1].
  - Block-coupled pressure-displacement (foam-extend), as used in [2].

---

## Case Overview

The case models the quasi-static inflation of an idealised, axisymmetric left
ventricle under an internal endocardial pressure that ramps linearly from
0 to 10 kPa over a unit pseudo-time. The myocardium is described with a
nearly-incompressible Guccione passive law; for the pressure-displacement
approach, an isotropic Guccione form is used (`cf = ct = cfs = 1`).

Two solution approaches are provided:

- `petsc` (default): solved with the
  `nonLinearGeometryTotalLagrangianTotalDisplacement` solid model and the
  PETSc SNES nonlinear solver. The mesh is generated with `blockMesh` and
  rotationally extruded with `extrudeMesh` to produce a full-revolution
  ventricle.
- `pressureDisplacement`: solved with the foam-extend
  `coupledPressureDisplacementSolid` block-coupled solid model on a
  symmetric (one-quarter) tetrahedral mesh imported from a Fluent `.msh`
  file.

---

## Running the Case

The tutorial case is located at
`solids4foam/tutorials/solids/hyperelasticity/idealisedVentricle`. The case
can be run using the included `Allrun` script:

```bash
./Allrun
```

The default option uses the PETSc SNES approach. The script also supports
the pressure-displacement solution option:

```bash
./Allrun petsc
./Allrun pressureDisplacement
```

The default `petsc` case requires OpenFOAM.com with PETSc enabled. The
`pressureDisplacement` option is foam-extend-only.

The default approach can also be run in parallel using:

```bash
./Allrun parallel
```

## Verification and Convergence Study

The opt-in [`verification/`](verification/) directory migrates the mesh
convergence study from `solid-benchmarks/hyperElasticity/idealisedVentricle`
into this tutorial. It refines the `blockMesh` divisions and the rotational
extrusion layers together, samples the deformed mid-wall line of Problem 2 of
Land et al. [3], and checks self-convergence of that line:

```bash
cd verification
./Allverify
```

The default sweep runs the 1 620-, 12 960-, and 103 680-cell levels; an
829 440-cell level is available explicitly. The study is separate from
`regressionTest.sh` and is not run by the normal tutorial test suites. See the
verification README for the mesh levels, options, and acceptance criteria.

---

## References

[1]
[P. Cardiff, D. Armfield, Ž. Tuković, and I. Batistić, "A Jacobian-Free
Newton-Krylov Method for Cell-Centred Finite Volume Solid Mechanics",
International Journal for Numerical Methods in Engineering, 127(3),
e70268, 2026.](https://doi.org/10.1002/nme.70268)

[2]
[A. Horvat, P. Milović, I. Karšaj, and Ž. Tuković, "A Block-Coupled
Finite Volume Method for Incompressible Hyperelastic Solids", Applied
Sciences, 15(23), 12660, 2025.](https://doi.org/10.3390/app152312660)

[3]
[S. Land, V. Gurev, S. Arens, et al., "Verification of cardiac mechanics
software: benchmark problems and solutions for testing active and passive
material behaviour", Proceedings of the Royal Society A, 471(2184), 20150641,
2015.](https://doi.org/10.1098/rspa.2015.0641)
