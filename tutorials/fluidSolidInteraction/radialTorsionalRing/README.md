---
sort: 12
---

# Radial-torsional FSI of an elastic ring in a viscous annulus: `radialTorsionalRing`

---

Prepared by Sairam Pamulaparthi Venkata and Philip Cardiff

---

## Tutorial Aims

- Provides a finite-amplitude fluid-solid interaction benchmark (neo-Hookean
  solid, convective fluid) with a semi-analytical reference solution.
- Couples radial inflation and torsion through the nonlinear solid response,
  and exercises swirl, wall torque and angular-momentum transport, which the
  other FSI verification cases do not.
- Shows an exact open outer boundary for a swirling, breathing flow (coded
  boundary conditions) and conservation diagnostics written by coded function
  objects.

---

## Case Overview

A neo-Hookean ring, $$A \le R \le B$$ with $$A = 0.5$$ m and $$B = 1$$ m, is
surrounded by a viscous fluid in the annulus $$B \le r \le c$$, $$c = 2$$ m
(Figure 1). The full 2-D annulus is modelled in plane strain, one cell deep,
with a conformal interface at $$r = B$$. The inner surface of the ring is
rotated by $$\alpha_A(t)$$ (`fixedRotation`) and is free radially; the outer
fluid boundary carries a pressure $$p_c(t)$$ and an exact open condition for
the velocity. Both drives start with a $$C^\infty$$ ramp, so the solution is
smooth in time.

![Figure 1: geometry and loading (schematic, interface deformation
exaggerated)](images/geometry.png)

**Figure 1:** Geometry and loading (schematic; interface deformation
exaggerated).

The fields depend on $$r$$ and $$t$$ only (the rotationally invariant
$$n = 0$$ branch). The reference solution (`verification/scripts/
rt_reference.py`) solves the reduced 1-D problem with Gauss-Lobatto-Legendre
modes in the solid and Chebyshev collocation in the fluid. It has been checked
against an independently coded coupled solver, symbolic derivations and
several limiting cases (solids4foam pull request #637).

### Table 1: Problem Parameters

<!-- markdownlint-disable MD013 -->

| Parameter | Value | Units |
| :-- | :-: | :-: |
| Solid inner and outer radius $$A$$, $$B$$ | 0.5, 1 | m |
| Fluid outer radius $$c$$ | 2 | m |
| Solid law | neo-Hookean, plane strain | |
| Young's modulus $$E$$, Poisson's ratio $$\nu$$ | 2e4, 0.3 | Pa, - |
| Solid and fluid density | 1000 | kg/m$$^3$$ |
| Fluid dynamic viscosity $$\mu$$ | 300 | Pa s |
| Torsion $$\alpha_A = \theta_0 S(t/T_r)\sin\Omega t$$ | $$\theta_0 = 0.1$$ rad, $$\Omega = 3$$ rad/s | |
| Outer pressure $$p_c = P_0 S(t/T_r)$$ | $$P_0 = -3000$$ | Pa |
| Ramp $$S$$ | $$C^\infty$$, $$T_r = 1$$ s | |
| Time step, end time | 0.0025, 2 | s |

<!-- markdownlint-enable MD013 -->

Three drives are defined: torsion only ($$P_0 = 0$$), radial only
($$\theta_0 = 0$$) and combined (both). The stored case is the combined drive
on mesh M1; `verification/Allverify` generates the others.

### Table 2: Meshes

| Mesh | $$n_\theta \times n_{rs} \times n_{rf}$$ | Solid cells | Fluid cells |
| :--: | :--------------------------------------: | :---------: | :---------: |
|  M1  |                64 x 8 x 16               |     512     |    1024     |
|  M2  |               128 x 16 x 32              |    2048     |    4096     |
|  M3  |               256 x 32 x 64              |    8192     |    16384    |

---

## Numerical Set-up

- **Coupling:** `IQNILS` with `directMap` transfer, initial relaxation 0.05,
  `couplingReuse 5`, `outerCorrTolerance 1e-6`. The $$C^\infty$$ ramp starts
  at round-off level, where the relative residual cannot reach its tolerance,
  so `outerCorrAbsoluteTolerance` is set explicitly (1e-12 m on M1 and M2,
  1e-10 m on M3).
- **Fluid:** `pimpleFluid`, laminar, `newMovingWallVelocity` on the interface
  (with its default continuity-based wall-normal gradient, #671),
  `displacementLaplacian` mesh motion, `backward` in time, 2 outer and 3
  pressure correctors with `nNonOrthogonalCorrectors 2`.
- **Solid:** `nonLinearGeometryTotalLagrangianTotalDisplacement` solved with
  PETSc SNES, `diffStencilLaplacian` stabilisation.
- **Outer boundary:** `outerPressureSeeded` (`codedFixedValue` for `p`) reads
  $$p_c(t)/\rho$$ from `constant/pcTable`; `openRadialSwirlFsi`
  (`codedMixed` for `U`) imposes the exact open condition
  $$\partial \mathbf{U}/\partial n = (\mathbf{U} - 2(\mathbf{n}\cdot
  \mathbf{U})\mathbf{n})/r$$. The drive tables are tabulated at every time
  step.
- **Diagnostics:** the coded function objects `rtDiagFluid` and `rtDiagSolid`
  in `system/controlDict` write, every time step, the interface traction,
  wall torque and power, the angular momentum and energies of both regions,
  the interface radius and rotation, $$n = 2, 3, 4$$ Fourier amplitudes, and
  the range of $$J$$, to `rtDiagFluid.dat` and `rtDiagSolid.dat`.

---

## Running the Case

```bash
./Allrun            # M1, combined drive, serial (a few minutes)
./Allrun parallel   # the same on the 4 ranks of decomposeParDict
```

The case needs solids4foam built with PETSc on OpenFOAM.com; the coded
conditions are not available on foam-extend or OpenFOAM.org. Where `python3`
with numpy and scipy is available, `Allrun` also writes the semi-analytical
reference (`reference.json`) and the errors, balances and quantities of
interest (`analysis.json`, `series.csv`); without them only the solver runs.

The mesh and drive study runs separately:

```bash
cd verification
./Allverify                    # M1, all three drives
./Allverify torsion M1 M2 M3   # torsion convergence study
```

See `verification/README.md`.

---

## Expected Results

Errors against the reference are max-norm over 0-2 s, as fractions of
max|reference|. For the **torsion** drive the rotational quantities converge
at close to second order:

<!-- markdownlint-disable MD013 -->

| Quantity | M1 | M2 | M3 | Observed orders |
| :-- | :-: | :-: | :-: | :-: |
| interface rotation $$\alpha_B$$ | 1.17e-2 | 2.97e-3 | 7.6e-4 | 1.98, 1.98 |
| wall torque $$T_B$$ | 2.84e-2 | 7.24e-3 | 1.9e-3 | 1.97, 1.94 |
| shear traction $$\sigma_{r\theta}$$ | 2.76e-2 | 7.04e-3 | 1.9e-3 | 1.97, 1.92 |
| drive torque $$T_A$$ | 2.83e-2 | 6.75e-3 | 1.6e-3 | 2.07, 2.04 |
| fluid angular momentum $$H_f$$ | 2.19e-2 | 5.86e-3 | 1.8e-3 | 1.90, 1.72 |

<!-- markdownlint-enable MD013 -->

These values were obtained before #671. The torsion drive does not stretch the
interface, so it is not expected to depend on the change, but the table has not
been regenerated.

In the **combined** drive the rotational errors stopped decreasing between M2
and M3 before #671, at about 3-7e-3 of max|reference|. Halving or quartering
the time step and tightening the coupling tolerance by three orders did not
change this floor, whereas a Gauss linear instead of `leastSquares` pressure
gradient removed it (M3 $$\alpha_B$$ error 7.6e-4 instead of 3.1e-3), which
points to the pressure-gradient discretisation.

Before #671 the viscous normal wall stress was transmitted at about one third of
the reference, without converging with the mesh. With the continuity-based
wall-normal gradient of #671, the stored M1 combined case gives:

| Quantity (M1, combined) | Before #671 | With #671 |
| :-- | :-: | :-: |
| transmitted viscous normal stress / reference | 0.29 | 0.90 |
| wall pressure $$p_b$$ error | 2.62e-2 | 1.74e-2 |
| normal traction $$\sigma_{rr}$$ error | 1.48e-2 | 1.38e-2 |
| interface displacement $$u_B$$ error | 1.78e-2 | 1.70e-2 |
| rotation $$\alpha_B$$ error | 1.27e-2 | 1.05e-2 |
| wall torque $$T_B$$ error | 2.28e-2 | 2.41e-2 |

The M2 and M3 combined and radial runs have not yet been repeated with #671.
