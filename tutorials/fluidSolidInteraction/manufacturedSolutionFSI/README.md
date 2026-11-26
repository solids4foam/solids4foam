---
sort: 16
---

# Method of Manufactured Solutions for FSI: `manufacturedSolutionFSI`

Prepared by Philip Cardiff, from the solids4foam verification campaign

## Tutorial aims

- Verify a transient fluid-solid interaction (FSI) solution on a moving
  (arbitrary Lagrangian-Eulerian, ALE) fluid mesh with the method of
  manufactured solutions (MMS).
- Measure the fluid, solid and interface errors against a closed-form
  solution, and the observed orders of accuracy under refinement.
- Provide a fast coarse case and an opt-in convergence study.

## Case overview

A two-dimensional channel flow is coupled to a compliant wall:

- the fluid occupies $$0 < x < L$$, $$0 < y < h(x, t)$$, with an inlet at
  $$x = 0$$, an outlet at $$x = L$$, a no-slip wall at $$y = 0$$ and the
  fluid-solid interface at $$y = h(x, t)$$;
- the solid occupies $$0 < X < L$$, $$H < Y < H + T$$ in its initial
  configuration, with the interface at $$Y = H$$, and is clamped on its top
  and ends.

The parameters are $$L = 2$$, $$H = 1$$, $$T = 1/2$$; fluid density
$$\rho_F = 1$$ and dynamic viscosity $$\mu_F = 1/100$$; a compressible
neo-Hookean solid (`neoHookeanElastic`, plane strain) with density
$$\rho_S = 1$$, $$E_S = 10$$ and $$\nu_S = 0.3$$. The density ratio is one,
so added-mass effects are strong. The period of the motion is one.

The manufactured solution is built so that both interface conditions, the
continuity of velocity and of traction, hold exactly. No interface source is
then needed, and the coupling itself is tested:

1. The solid displacement is chosen in closed form. Its trace on the
   interface gives the interface motion: a tangential Moebius map of
   $$[0, L]$$ onto itself, which stretches the interface by up to 10%, and a
   normal motion prescribed as the channel height
   $$h(x, t) = H + \frac{1}{10} \sin^2(\pi x / L)\, r(t)
   \left[1 + \frac{1}{2} \sin(2 \pi x / L - \omega t)\right]$$,
   with $$r(t) = (1 - \cos \omega t)/2$$, a pulsing wave with a travelling
   part. The interface tilts, travels and stretches.
2. The fluid velocity comes from a stream function, so it is divergence free.
   Hermite polynomials in $$y/h$$ impose no slip on the bottom wall and the
   solid velocity on the deformed interface.
3. A term that changes only the second normal derivative of the stream
   function on the interface makes the fluid shear traction equal to the
   solid shear traction, and the interface pressure matches the normal
   traction.
4. The fluid and solid momentum sources follow from the residuals of the
   governing equations.

The solution starts from rest: the displacements, velocities and pressure are
zero at $$t = 0$$. At $$t = 0.4$$ the fluid speed reaches 2.15, the solid
Green-Lagrange strain 0.26, and the solid Jacobian lies between 0.69
and 1.04.

The tutorial-local library `src/` (`libmmsFsi`) provides:

- `mmsFsiSource`, an `fvOption` with the fluid momentum source (per unit
  mass, at the current fluid cell centres) and the solid body force (per unit
  reference volume);
- the exact fluid boundary conditions `mmsFsiVelocity` (inlet),
  `mmsFsiVelocityGradient` (outlet, the exact normal velocity gradient) and
  `mmsFsiPressure` (outlet);
- the `mmsFsiErrors` function object, which writes the L1, L2 and LInf errors
  every time step to `postProcessing/mmsFsiErrors/errors.dat`: the fluid
  velocity `U` and pressure `p` on the moving mesh, the solid displacement
  `D` and velocity `Us`, the interface displacement `intD`, the traction
  applied to the solid interface `intTrac`, and the fluid interface pressure
  `intP`. The fluid error fields `UError` and `pError` are written at the
  write times.

The exact solution is evaluated by the generated C++ header
`src/mmsFsiExact.H`; see [Regenerating the exact solution](#regenerating-the-exact-solution).

## Running the case

The case uses the partitioned `pimpleFluid` fluid model, the
`nonLinearGeometryTotalLagrangianTotalDisplacement` solid model with the
PETSc SNES solution algorithm, and the `displacementLaplacian` mesh motion. It
requires OpenFOAM.com and a PETSc-enabled solids4foam build; `Allrun` skips
the case otherwise.

```bash
./Allrun          # Dirichlet-Neumann coupling with Aitken relaxation
./Allrun robin    # Robin-Neumann coupling
```

The default mesh has 16 cells along the interface: 16 x 8 square cells in
the fluid and 16 x 4 in the solid. The time step is 0.02, and the case runs
to $$t = 0.4$$ in a few seconds. The mesh and time step are set in
`system/caseParameters`. `Allrun` compiles the library, links
`constant/fsiProperties` and the fluid interface conditions
(`0/fluid/U` and `0/fluid/p`) for the chosen coupling, creates the meshes and
runs `solids4Foam`.

- `aitken`: `newMovingWallVelocity` and `fixedFluxPressure` on the fluid
  interface, Aitken relaxation (initial factor 0.05);
- `robin`: `elasticWallVelocity` and `elasticWallPressure` on the fluid
  interface (`fixedRelaxation` with factor 1).

Both use an outer coupling tolerance of 1e-6. The inlet velocity, the outlet
velocity gradient and the outlet pressure are exact; the inlet and wall
pressures are `fixedFluxPressure`.

## Verification and regression

`./regressionTest.sh` runs the coarse case with both couplings and checks the
L2 errors of `U`, `p`, `D` and `intTrac` at $$t = 0.4$$ against the reference
errors, with a relative tolerance of 2%. It is part of
`tutorials/Alltest-regression`.

The opt-in [`verification/`](verification/) study refines the mesh and the
time step together, $$\Delta t = 0.32/n$$ for $$n$$ = 16, 32 and 64 cells
along the interface, and reports the observed orders:

```bash
cd verification
./Allverify                          # Aitken, n = 16, 32, 64
./Allverify --couplings aitken,robin
./Allverify --quick                  # n = 16 and 32 only
```

The partitioned fluid model is first order in this combined refinement. The
expected L2 orders between n = 32 and 64 with Aitken relaxation are about 1.3
for `U`, 1.0 for `p`, 0.9 to 2 for `D` and `intD`, and 0.8 to 1.4 for
`intTrac` and `intP` (at $$t = 0.2$$ and 0.4); see the verification README for
the tables.

## Regenerating the exact solution

The generated header `src/mmsFsiExact.H` is committed, so the case needs
neither Python nor sympy. The construction is in `derivation/`:

- `mms_fsi.py` derives the solution with sympy and symengine, runs 11
  symbolic and 40-digit checks (divergence-free velocity, interface velocity
  and traction continuity, boundary conditions, zero initial state, and
  others), writes 40-digit reference values and generates the C++ header;
- `codegen.py` writes the header with common-subexpression elimination;
- `check_generated.C` checks the generated header independently: the fluid
  momentum residual and the solid neo-Hookean residual with finite
  differences, the interface conditions, and the 40-digit reference values.

```bash
cd derivation
./run_checks.sh     # about six minutes; Python 3 with sympy and symengine
```

The default call, the full variant, writes `derivation/generated/` with the
same solution as `src/mmsFsiExact.H` and `derivation/reference_values.txt`.
The files need not be byte-identical, as the order of the common
subexpressions and the last of the 40 digits can change from run to run;
`check_generated.C` compares the header with the 40-digit values. Its
arguments select other variants, which the case does not use: for example
`./run_checks.sh generated_flat 1 1 flat` gives the `flat` variant, in which
the interface stays flat and moves only vertically, and the density and
viscosity arguments change the fluid and solid parameters.

## Not included

The verification paper also runs this case with the quasi-monolithic
coupling of `newtonIcoFluid` and the velocity-form total Lagrangian solid,
and with the `flat` variant. These need solver changes that are not on the
`development` branch, so they are not part of this tutorial. IQN-ILS with
the default settings of this case does not converge on the 32-cell mesh.

## Reference

- The solids4foam fluid-solid interaction verification paper (in
  preparation), section "Manufactured solution for transient ALE FSI with a
  finite-strain solid" and its appendix, which give the construction in full.
