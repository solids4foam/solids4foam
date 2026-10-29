---
sort: 12
---

# Radial-torsional FSI of an elastic ring in a viscous annulus: `radialTorsionalRing`

---

Derived from a separate benchmark study of a radially breathing and torsionally twisting
neo-Hookean ring coupled to a viscous fluid. The verification is scoped: read the Scope and
Claim Boundaries before using any result. The provenance of every setting is recorded
outside the tutorial.

---

## Tutorial Aims

- Provides a finite-amplitude (nonlinear solid, convective fluid) fluid-solid
  interaction benchmark whose reference is a semi-analytical reduced model, with
  three drives: torsion only, radial (pressure) only, and both.
- Shows an exact open outer boundary for a swirling, breathing flow (coded
  boundary conditions), a conformal 2-D annulus mesh at three resolutions, and
  IQN-ILS coupling with an explicit per-mesh absolute coupling tolerance.
- Provides a generator, a reference-solution script, a quantity-of-interest
  extraction, a convergence and mixed-term analysis, and conservation
  diagnostics (angular momentum and energy balances) as function objects.

---

## Scope and Claim Boundaries

What the tutorial supports, and no more:

- **Torsion.** 0-2 s three-mesh verification at $$\Delta t = 0.0025$$ s and
  `nNonOrthogonalCorrectors 2`: near-second-order spatial convergence is **observed** for the primary
  rotational quantities. The measured orders (Table 5) are 1.72-2.07 in the max-norm and
  1.71-2.04 in L2-in-time over the two windows (0-2 s and 0-1.8 s); the value 1.68 that appears in
  the audits is the dt-corrected $$H_f$$ order, not a plain one. They are not all in 1.8-2.0 and are never
  rounded to it. In the max-norm, outside 1.8-2.0 are the fluid angular momentum $$H_f$$ at 1.72 (second pair,
  both windows, L2 1.71; its first pair 1.90 is inside) and the drive torque $$T_A$$ at 2.07 and 2.04
  (0-2 s) and 2.02 (0-1.8 s); several L2 values are 2.01-2.04 (Table 5). Inside 1.8-2.0, for example, are $$T_B$$ (1.97, 1.94 over 0-2 s;
  1.94, 1.85 over 0-1.8 s) and $$\sigma_{r\theta}$$ (1.97, 1.92; 1.94, 1.83). An asymptotic order is
  **not established** (three meshes), and there is **no temporal order claim**. At the documented
  numerical noise the second-pair orders are uncertain, each quantity with its own interval (Table 6, 0-2 s, max-norm):
  $$T_B$$ 1.94 is 1.75-2.12 ($$\pm$$0.19), $$\sigma_{r\theta}$$ 1.92 is 1.73-2.11 ($$\pm$$0.19), $$T_A$$ 2.04 is 1.82-2.25 ($$\pm$$0.22),
  $$H_s$$ 1.95 is 1.64-2.26 and $$b\omega_B$$ 1.82 is 1.51-2.12 (both about $$\pm$$0.31, wider than the $$\pm$$0.19-0.22 of the torque-type
  quantities), $$P_{wall}$$ and $$KE_s$$ $$\pm$$0.21 and $$\pm$$0.18, while $$\alpha_B$$ and $$H_f$$ are sharp (zero width to two digits). Intervals wider than $$\pm$$0.15 are
  labelled noise-limited by the driver. The intervals qualify the orders, they do not turn an observation into a proof.
- **Combined.** Full-duration three-mesh results with nominal
  self-convergence orders near 2 (sharply resolved only for $$\alpha_B$$ and $$H_f$$; the others noise-limited, see below), with a disclosed reference-error floor of
  about 2.5-5.2e-3 of max|ref| **for $$\alpha_B$$, $$T_B$$, $$\sigma_{r\theta}$$, $$T_A$$ and $$H_f$$ only**
  (measured on M3: 3.06e-3 for $$H_f$$ to 5.15e-3 for $$\sigma_{r\theta}$$; 2.5e-3 is the lower edge of the
  band, not a measured value). The floor does not describe the other quantities: on M3 the error is
  6.9e-3 for $$H_s$$, 6.3e-3 for $$b\omega_B$$ (wall angular velocity), 1.06e-2 for $$KE_s$$ and 4.8e-2 for the interface
  power $$P_{wall}$$ (Table 7). This is **not** clean second-order convergence to the reference, and the floor is
  not a verification PASS. The self-convergence orders (0-2 s) are near 2 for $$\alpha_B$$ (1.97), $$T_B$$ (2.00),
  $$\sigma_{r\theta}$$ (1.97), $$T_A$$ (2.09) and $$H_f$$ (1.94), **but only $$\alpha_B$$ and $$H_f$$ are robustly resolved**
  (intervals $$\pm$$0.01: 1.96-1.97 and 1.93-1.94). $$T_B$$ (2.00, interval 1.39-2.60), $$\sigma_{r\theta}$$ (1.97, 1.39-2.54) and $$T_A$$
  (2.09, 1.63-2.55) have nominal orders near 2 but noise-limited intervals, so they do **not** independently demonstrate a precise
  second-order rate. There is no second-order self-convergence for $$H_s$$ (1.63), $$b\omega_B$$ (1.66) and $$P_{wall}$$ (0.70), and $$KE_s$$
  (2.17, noise-limited interval 1.42-2.92) has a reference error above the floor. The origin of the floor is only **partially explained**:
  about one third of the rotation excess follows from the radial wall-traction (viscous normal stress) discrepancy, which is not a proven unique mechanism.
- **Mixed terms.** $$\Delta\alpha$$, $$\Delta T$$ and $$\Delta\sigma_{r\theta}$$
  (combined minus torsion-only) are reproduced to about 1.5-2.5% of the signal on M3
  over 0-2 s ($$\Delta\sigma_{r\theta}$$ up to about 3.4% over 0-1.8 s). Their orders (against the
  reference and self) lie between 0.2 and 1.4; no order is claimed for them.
- **Radial only: verification INCOMPLETE.** At `nNonOrthogonalCorrectors 2` only
  M1 completed; M2 stalled at t = 0.11 s and M3 was not run. No M1-M3 radial
  convergence is claimed or implied. The radial configuration is provided as an
  input set, and its M1 run is reported only as a single-mesh result.
- **Limitations.**
    - Temporal contamination: 5-25% of the M3 error, with no temporal order.
    - Iteration and partition noise of a max-norm error (one statement, with provenance in
      `verification/thresholds.json`): about 2.4e-4 of max|ref| for the torsion runs (restart and
      attempt-to-attempt differences of the M3 torsion torque 1.8e-4 to 2.4e-4) and about 1.3e-3 of
      max|ref| for the combined runs (serial against 4-rank M2 combined: torque 1.26e-3, shear 1.21e-3,
      $$T_A$$ 8.9e-4; the absolute-tolerance change gives 8.8e-4, 8.6e-4, 5.8e-4). The "about 1e-3
      torque and shear band" of the audits is the combined-run value, not the torsion one. Rotation and
      angular momentum are three orders smaller. Scotch partitions are not reproducible between invocations.
    - The absolute coupling tolerance differs between meshes (1e-12 m on M1 and M2,
      1e-10 m on M3). Its effect was quantified on M2 only; on M3 it is not bounded.
    - Normal-traction discrepancy on M3: $$u_B$$ about 1.1%, $$\sigma_{rr}$$ about
      0.85%, $$p_b$$ about 2.4%. These quantities are excluded from every order claim.
    - The wall viscous normal stress is transmitted at about 0.29-0.32 of the reference value
      (ratio 0.286, 0.307, 0.317 on M1, M2, M3 combined; 0.295 on radial M1) without
      mesh convergence (see Known Caveats). It is excluded from every order claim.
    - The ALE/GCL consecutive-save checks are incomplete.
    - The reference is only numerically resolved. No second, independent nonlinear
      implementation has validated it.

---

## Case Overview

A solid ring $$A \le R \le B$$ ($$A = 0.5$$ m, $$B = 1$$ m) surrounds the axis and
is surrounded by a viscous fluid $$B \le r \le c$$ ($$c = 2$$ m). The full 2-D
annulus is modelled (plane strain, one cell in depth) with a conformal interface at
$$r = B$$. This is the planar, rotationally invariant ($$n = 0$$) branch: the fields depend on $$r$$ and
$$t$$ only, and $$n \ne 0$$ modes are monitored through the Fourier amplitudes of the
diagnostics, not driven. The reference uses Gauss-Lobatto-Legendre (GLL) solid modes and Chebyshev
fluid modes (collocation). The inner solid surface is driven in rotation, $$\alpha_A(t)$$, with
`fixedRotation`; the solid is free of radial displacement there. The outer fluid
boundary carries an imposed pressure $$p_c(t)$$ and an exact open condition for the
velocity.

**Torque units.** $$T_B$$, $$T_A$$, $$T_c$$ and the angular momenta $$H_f$$, $$H_s$$ are per radian and per
metre of depth (the `rtDiag` column divided by $$2\pi\,dz$$: N m per rad per m of depth, i.e. N). The torque of the
whole annulus per metre of depth is $$2\pi$$ times that value; for example the signal of the mixed term
$$\Delta T$$, 74.4, is 468 N m per m of depth for the whole annulus.

### Table 1: Problem Parameters

|               Parameter               |      Value       |   Units    |
| :-----------------------------------: | :--------------: | :--------: |
|   Solid inner / outer radius $$A$$, $$B$$ |   0.5 / 1     |     m      |
|   Fluid outer radius $$c$$            |        2         |     m      |
|   Solid law                           | neo-Hookean, plane strain |   |
|   Young's modulus $$E$$               |   2e4            |    Pa      |
|   Poisson's ratio $$\nu$$             |       0.3        |            |
|   Solid density $$\rho_s$$            |      1000        | kg/m$$^3$$ |
|   Fluid density $$\rho$$              |      1000        | kg/m$$^3$$ |
|   Fluid dynamic viscosity $$\mu$$     |      300         |   Pa s     |
|   Fluid convection                    |        on        |            |
|   Torsion $$\alpha_A = \theta_0 S(t/T_r)\sin\Omega t$$ | $$\theta_0 = 0.1$$ rad, $$\Omega = 3$$ rad/s | |
|   Outer pressure $$p_c = P_0 S(t/T_r)$$ (ramp-hold)  | $$P_0 = -3000$$ Pa | Pa |
|   Ramp $$S$$                          | C-infinity, $$T_r = 1$$ s | |
|   Time step, end time                 | 0.0025, 2.0      |     s      |

The three configurations are torsion only ($$\theta_0 = 0.1$$, $$P_0 = 0$$), radial
only ($$\theta_0 = 0$$, $$P_0 = -3000$$ Pa) and combined (both).

### Table 2: Meshes ($$n_\theta \times n_{rs} \times n_{rf}$$)

| Mesh | $$n_\theta \times n_{rs} \times n_{rf}$$ | Solid cells | Fluid cells | Decomposition (M3) |
| :--: | :--------------------------------------: | :---------: | :---------: | :----------------: |
|  M1  |                64 x 8 x 16               |     512     |    1024     |         -          |
|  M2  |               128 x 16 x 32              |    2048     |    4096     |         -          |
|  M3  |               256 x 32 x 64              |    8192     |   16384     |  4 ranks, scotch   |

---

## Numerical Set-up

- **Coupling.** `IQNILS` with `directMap` on the interface, initial relaxation
  0.05, `couplingReuse 5`, `nOuterCorr 100`, `outerCorrTolerance 1e-6`, with the
  predictors on. `allowUnconvergedCoupling` is **not** set, so an unconverged
  step aborts the run.
- **Absolute coupling tolerance (explicit, per mesh).** `outerCorrAbsoluteTolerance`
  is 1e-12 m for M1 and M2 and 1e-10 m for M3, as in the accepted runs. It is
  needed because the C-infinity ramp starts at round-off level, where the
  relative residual cannot reach 1e-6 (a 1e-12 m attempt on M3 stalled). Override
  it with `./Allrun absTol=X`. The setting is written by `makeCase.py` into
  `constant/fsiProperties`; the stored dictionary is the M1 value.
- **Convergence test (the solver's actual acceptance logic).** An FSI step is accepted when
  `min(residual1, residual2) <= outerCorrTolerance`, that is when **either** relative residual measure
  satisfies its tolerance (`requireAllResidualMeasures` is unset in every accepted run), or when the
  absolute point mismatch is below `outerCorrAbsoluteTolerance` (`accepted by outerCorrAbsoluteTolerance`,
  printed only if both relative measures are above tolerance). With `allowUnconvergedCoupling` unset an
  unconverged step aborts the run. The driver records the criterion of the **last iteration** of each step
  (Table 4); this definition counts 68 and 70 absolute acceptances for the M1 and M2 torsion runs, where the
  audits' printed-line counts are 76 and 83 (they include acceptances printed at the solid-predictor evaluation
  before iteration 1, which is not a loop iteration). The counts of the other cases agree.
  A stricter test, both measures at 1e-6, is **not** what these runs satisfy. Two different
  quantities must not be mixed (Table 4): **(A)** the steps failing the both-relative-measures test (the absolute-tolerance steps count as failing), 93.0-97.5% of the
  steps of the six torsion and combined runs; and **(B)** the steps the solver's own strict mode (`requireAllResidualMeasures`, which still accepts
  absolute-tolerance steps) would reject at the accepted iterate, 82.1-92.25% over the eight runs (83.4-92.9% of the steps are accepted on a single relative measure;
  5-14 of those per run also have an absolute residual below the tolerance and would be accepted by strict mode too) (A: 93.0-97.5% of the six runs, 97.88% for the absolute-tolerance test case and 93.88% for radial M1; 17-56 of 800 steps meet both measures). The runs are therefore **not**
  labelled converged under that stricter criterion. Proxy evidence only: the radial M2 pilot at relative tolerance
  1e-6 and 1e-8 (800 steps each; at 1e-8, 334 of 701 late steps meet both measures at 1e-6) differs by
  u_B 8.5e-5, sigma_rr 4.0e-4, p_b 3.9e-4, P_wall 4.3e-4 and KE_s 2.5e-4 of max|ref|, which is 0.7-4% of the
  M2 errors. The rotational quantities were never run at a tighter relative tolerance; the disclosed
  combined-run noise of about 1.3e-3 of max|ref| (torque and shear) covers the proxy effect.
- **Pressure-solver iteration cap.** In every radial and combined run (M1, M2, M3, F3B and the stalled radial M2) three
  start-up steps (t = 0.0275-0.0325 s) have 51 pressure solves that stop at `maxIter` 1000 with a final residual
  above 1e-12 (2e-9 to 7e-9) because the pressure is at round-off level (|p| of 1e-16 to 1e-13 of the maximum). This costs
  0.8% (M3 combined) to 1.7% (M1 combined) of the run time and is immaterial according to the audit.
  `maxIter 1000` is kept so that the stored regression values and the accepted runs are reproduced bit for bit.
- **Fluid.** `pimpleFluid`, laminar, `newMovingWallVelocity` on the interface,
  `displacementLaplacian` mesh motion, `backward` time scheme,
  `nOuterCorrectors 2`, `nCorrectors 3`, `momentumPredictor on`, `p` and `U`
  solved to 1e-12 and 1e-13 (GAMG and PBiCG/DILU, `relTol 0`),
  `fluxConsistentPatches (outer)`.
- **Pressure correction: `nNonOrthogonalCorrectors 2`.** This is the verification
  setting. With 0 non-orthogonal correctors the M3 torsion run fails near 1.81 s
  (the fluid PIMPLE iteration loses contractivity; no unique root cause was identified, and a different
  predictor setting or half the time step also removes it), so M3 results exist only
  with 2. M1 and M2 were also run with 2, so that the family is consistent.
- **Outer boundary.** `outer` uses two coded conditions, both identical in all
  accepted runs:
    - `outerPressureSeeded` (`codedFixedValue` for `p`): reads $$p_c(t)/\rho$$ from
      `constant/pcTable`. A seeded $$n = 2$$ pulse (`seed`, 0 here) is part of the
      code but inactive.
    - `openRadialSwirlFsi` (`codedMixed` for `U`): the exact open condition
      $$\partial U / \partial n = (U - 2 (n\cdot U) n)/r$$ with zero value fraction.
- **Drive tables.** `constant/timeVsAngle` (degrees) and `constant/pcTable`
  ($$p_c/\rho$$) are tabulated at every time step, so the interpolation is exact
  at step times.
- **Solid.** `nonLinearGeometryTotalLagrangianTotalDisplacement` with PETSc SNES
  (Newton, matrix-free operator, `lgmres`, block-Jacobi LU), `diffStencilLaplacian`
  stabilisation of 0.5, `backward` $$\partial^2/\partial t^2$$.
- **Conservation diagnostics.** The function objects `rtDiagFluid` and
  `rtDiagSolid` (coded, in `system/controlDict`) write one row per time step to
  `rtDiagFluid.dat` and `rtDiagSolid.dat`: interface and outer loads, wall torque,
  angular momentum, kinetic and strain energy, dissipation, radius and rotation of
  both solid surfaces, $$n = 2, 3, 4$$ Fourier amplitudes, and the solid
  principal stretch and $$J$$ range. `verification/scripts/analyse_case.py` forms the
  balances from them:
    - fluid angular momentum, $$T_{bal} = -dH_f/dt + T_c - F_c$$ (against the wall torque);
    - solid angular momentum, $$T_A = dH_s/dt - T_B$$ (against the drive torque);
    - fluid energy residual, $$dKE_f/dt + P_{wall} - P_c + D_f$$;
    - solid energy balance, $$d(KE_s + W_s)/dt - P_{wall}$$ against the reference;
    - each reported per window 0-0.5, 0.5-1, 1-1.5, 1.5-1.8, 1.8-2 and 0-2 s, relative to the reference
      torque or interface power in the window.

---

## Running the Case

```bash
./Allrun                          # M1, combined, serial
./Allrun M2 torsion               # serial
./Allrun M3 combined parallel     # 4 ranks
./Allrun M3 torsion parallel np=4
./Allrun M1 radial endTime=0.5
```

Options, in any order: `M1|M2|M3`, `torsion|radial|combined`, `serial|parallel`,
`np=N`, `nNonOrth=N`, `absTol=X`, `endTime=T`, `regenerate`, `reconstruct`.

**Stored or regenerated dictionaries.** A plain `./Allrun` (or `./Allrun serial`, or `./Allrun parallel`) uses the
**stored** dictionaries, the M1, combined set with 4-rank decomposition, without rewriting any of them, and only writes
the reference solution; any script that edits the dictionaries beforehand (such as `tutorials/Alltest`) is respected.
Any of `M1|M2|M3`, `torsion|radial|combined`, `nNonOrth=`, `absTol=`, `endTime=`, `regenerate`, or `np=N` with N other than 4
in parallel runs `verification/scripts/makeCase.py`, which rewrites the meshes, tables, settings and reference for that
selection, including `controlDict`. After a regenerating run, `./Allclean` restores the stored set; a default `./Allrun`
refuses to run on a stored set that is not the M1 one.

`Allrun` then runs `blockMesh`, `checkMesh`, optionally `decomposePar`, the solver, and `analyse_case.py` (skipped
for a shortened run with fewer than 10 diagnostic rows). The stored set also runs by hand without Python after
`blockMesh -region solid` and `blockMesh -region fluid`. Run `./Allclean` before switching to another variant,
as `Allrun` does not rerun a solver whose log exists.

The tutorial needs solids4foam built with PETSc and OpenFOAM.com (coded
conditions and function objects are not available in foam-extend or
OpenFOAM.org). `python3` with numpy and scipy is needed for the generator and
the reference (tested with Python 3.12.3, numpy 1.26.4, scipy 1.11.4;
OpenFOAM v2512).

### Table 3: Cost (dt = 0.0025 s, 800 steps; wall time from the `time.txt` of each accepted run)

| Mesh | Configuration | Ranks | Wall time |
| :--: | :-----------: | :---: | :-------: |
|  M1  | torsion       |   1   | 90.2 s (90.19 s) |
|  M1  | radial        |   1   | 91.1 s (91.07 s) |
|  M1  | combined      |   1   | 98.6 s (98.56 s) |
|  M2  | torsion       |   1   | 328.4 s (328.37 s) |
|  M2  | combined      |   1   | 408.4 s (408.38 s); 417.9 s with `absTol=1e-10` |
|  M3  | torsion       |   4   | 807.7 s (807.73 s) |
|  M3  | combined      |   4   | 1178.4 s (1178.40 s) |

The wall times are those of the accepted runs of the study, measured on one desktop
machine; M1 and M2 ran in serial, M3 on 4 ranks. The first run in a case directory also compiles the coded conditions (about 20 s on M1).

---

## Verification Workflow

The verification driver keeps two layers apart. They have separate statuses, and a passing
numerical layer is never a scientific verification.

```bash
verification/Allverify                   # M1: torsion, radial, combined (about 5 min)
verification/Allverify torsion M1 M2
verification/Allverify combined M1 M2 M3 np=4
verification/Allverify analyse           # re-analyse the cases in verification/work
python3 verification/scripts/verify.py WORK_DIR [--out DIR] [--log LOG_OR_CASE ...]
python3 verification/test_verify.py      # tests without CFD (synthetic logs, Allrun paths, family logic)
```

`Allverify` copies the tutorial into `verification/work/<config>_<mesh>`, runs it and calls
`scripts/verify.py`, which reads the cases (never writing into them) and writes
`verification_report/` with `summary.txt`, `verification_report.json`, `layerA_runs.csv`,
`layerA_steps_<run>.csv` (the criterion that accepted every time step), and `layerB_errors.csv`. The exit
code is 1 if any run is FAIL, INCOMPLETE, INVALID-INPUT or UNKNOWN in layer A (and if a run could not be started). All
thresholds, each with its provenance, are in `verification/thresholds.json`, in separate sections: `layerA`,
`layerB_deterministic` (bands from the accepted results), `numerical_noise` (partition, restart and iteration
noise allowances, from existing evidence only) and `family_consistency`. `Allclean all` removes
`verification/work`. The reference alone: `python3 verification/scripts/make_reference.py CASE --config combined --mesh M2` (about 2 s).

**Layer A: numerical acceptance (per run, from the solver log and diagnostics).**
Statuses PASS, CONDITIONAL, FAIL, INCOMPLETE; NOT-VERIFIED for checks that need field data. One ordered severity definition is shared by
all modules (`verification/scripts/status.py`: PASS < INFO < NOT-VERIFIED < CONDITIONAL = DISCLOSED = ABOVE-FLOOR < UNKNOWN < INVALID-INPUT < INCOMPLETE < FAIL; a combination of statuses is never PASS unless all are PASS, and the status of nothing is NOT-VERIFIED); an unrecognised
status is reported as UNKNOWN (non-PASS) and blocks the exit code like FAIL and INCOMPLETE.

- Termination: `End`, last time equals `endTime`, no `FOAM FATAL`, no NaN, no SNES failure. A log
  that stops without `End` and without an error (killed, truncated) is **INCOMPLETE**, never PASS; so is
  a log with `End` that stops short of `endTime`.
- Time steps monotone, unique, equal to `deltaT`.
- Every FSI step is accepted by an explicit criterion recorded per step from its last iteration:
  `both` relative measures at or below `outerCorrTolerance`, `r1-only`, `r2-only`, or `absolute`
  (see Numerical Set-up). The counts per criterion and the absolute-accepted times are reported.
  Single-measure and absolute acceptances are flagged explicitly (the other residual is above tolerance), so
  a run with any of them is CONDITIONAL. A step with no accepting criterion, an absolute acceptance that is
  not below the absolute tolerance, or a step beyond `nOuterCorr` is FAIL. An absolute acceptance after
  t = 0.35 s is CONDITIONAL, not a failure.
- Iterations: at most 30 for t >= 0.25 s (FAIL above 30); before 0.25 s (start-up) larger counts are only
  reported. The M2 combined run has one 34-iteration step at t = 0.0925 s, inside the exception.
- 3 pressure solves per continuity line (`nNonOrthogonalCorrectors` + 1).
- Linear finals (p 1e-12, U 1e-13): exceedances in the start-up (t < 0.25 s) are CONDITIONAL, later ones FAIL
  (see the pressure-solver cap in Numerical Set-up).
- Continuity error (informational, with the radial-flow rebase); contraction (geometric mean per window and runs of
  steps above 1; start-up excluded to 0.5 s; windows above 0.90 between 0.5 and 1.5 s are CONDITIONAL, from 1.5 s
  FAIL, as in the audits); solid minimum $$J$$ (> 0.99); diagnostics row consistency.
- Angular-momentum and energy balances are informational; they exceed 1e-2 in most accepted runs and decrease
  with mesh refinement in torsion.
- NOT-VERIFIED: minimum cell volume, swept-volume GCL, ALE, $$n \ne 0$$ mode growth (field scripts are not in the
  tutorial).

**Layer B: scientific reference-error expectations (per family of three meshes).** Gates, in this order:
(1) layer A of every member is PASS or CONDITIONAL, otherwise **INVALID-INPUT**;
(2) **family consistency**, read from the actual dictionaries of every case (never from metadata) and compared semantically (comments,
whitespace, number formats and the order of entries do not matter). What is checked: physical parameters (solid and fluid $$E$$, $$\nu$$,
$$\rho$$), $$\theta_0$$ and $$P_0$$ (the drive tables are recomputed with the reference and compared), the mesh of **both** regions (cell
counts and refinement ratio 2, all vertex coordinates and arc points, that is the inner radius $$A$$, the interface radius $$B$$, the outer radius
$$c$$ and the depth, the block topology, the grading and the patch definitions of the `blockMeshDict`s), the whole `controlDict` (every top-level
entry **before and after** the `functions` block, for example `deltaT`, `adjustTimeStep`, `writeControl`, `startFrom`, except `endTime` and
`writeInterval`, the parameters of the function objects, and the normalised code of `rtDiagFluid` and `rtDiagSolid` against the documented
current and older versions), PIMPLE settings (`nNonOrthogonalCorrectors`, `nOuterCorrectors`, `nCorrectors`, `momentumPredictor`,
`fluxConsistentPatches`) and the linear solvers, the coupling algorithm (`IQNILS`), `interfaceTransferMethod`, the predictors and the FSI tolerances,
the solid model, solution algorithm, stabilisation (type and `scaleFactor`) and material law, the fluid model, the dynamic mesh motion (solver and
diffusivity), the boundary conditions of `0/fluid/U`, `p`, `pointDisplacement` and `0/solid/D` (patch types and the normalised text of the coded
conditions, case paths and comments ignored), and a semantic fingerprint of every small dictionary against the documented benchmark set and against the other
meshes (so that a changed entry without a named check is still found). Not checked: the solver binary and library versions, the points of the generated `polyMesh`
(the `blockMeshDict`s that generate them are checked), function-object code other than the two documented versions, and the log files (layer A).
A missing run, or a missing, empty, unreadable or garbled dictionary, or a truncated drive table gives **INCOMPLETE** or **INVALID-INPUT**, never a crash and never
PASS. An inconsistent setting between the meshes, or a uniform departure from the documented benchmark set (for example Aitken instead of IQN-ILS, a different outer
radius, a changed `deltaT` anywhere in the `controlDict`), is **INVALID-INPUT** and stops the assessment. Documented-equivalent differences only make the family
**CONDITIONAL**: the absolute tolerance 1e-12 m on M1 and M2 against 1e-10 m on M3 (bounded on M2 only), the older diagnostics function objects and the older
seeded-pressure code with seed 0 (both on M2), and different rank counts and scotch partitions. **A CONDITIONAL family never reaches a layer-B PASS** (its claims are capped at
CONDITIONAL). The accepted torsion and combined families are CONDITIONAL for this reason;
(3) errors against the reference (max-norm and L2-in-time, 0-2 s and 0-1.8 s), observed orders with their first-order noise intervals,
self-difference orders and mixed terms, each compared with a deterministic band from the accepted results (PASS inside the band
shrunk by the numerical-noise allowance where one applies, CONDITIONAL inside the band but within the allowance of
an edge or within 20% outside, FAIL beyond; an order interval that leaves the band, or that is unresolved because the noise is not below the error, downgrades and never upgrades; intervals wider than $$\pm$$0.15 are labelled noise-limited in the output). The mixed-term
allowance of each term is noise $$\times$$ max|reference quantity of the combined run| / mixed-term signal, from the actual magnitudes
(for example $$\Delta T$$: 1.3e-3 $$\times$$ 454 / 74.4 = 0.79% of the signal; $$\Delta\sigma_{r\theta}$$: 1.21e-3 $$\times$$ 418 / 90.3 = 0.56%).
Marginal values are therefore CONDITIONAL: the combined $$\sigma_{r\theta}$$ error of 5.15e-3 on M3, $$\Delta T$$ (2.47% over 0-2 s,
2.65% over 0-1.8 s) and $$\Delta\sigma_{r\theta}$$ of 3.36% over 0-1.8 s. The combined reference-error floor is judged **by
quantity** (five quantities inside the floor statement, `ABOVE-FLOOR` for the others) and is never a PASS: its best status is **DISCLOSED**.

Every report ends with the list of NOT-VERIFIED items: the radial-only M1-M3 family (incomplete), the
three-way radial term $$d_{u3}$$, asymptotic convergence, temporal order, the origin of the combined floor,
the absolute-tolerance effect on M3, the field-level checks, and the independent validation of the reference.

Conventions of the analysis (identical to the study):

- Error: $$E = \max_t |q - q_{ref}| / \max_t |q_{ref}|$$ over $$t > 2\Delta t$$ in
  the window, on the reference at every step. Observed
  order $$\log_2(E_1/E_2)$$, only if $$E_1 > E_2 > E_3$$.
- Self-convergence: $$d_{12} = \max|q_1 - q_2|$$ and $$d_{23} = \max|q_2 - q_3|$$
  (M2 reference scale), $$p_{self} = \log_2(d_{12}/d_{23})$$.
- Primary rotational quantities (order claims): `alpha_B`, `w_wall`, `T_B`,
  `sig_rt`, `T_A`, `H_s`, `H_f`, `P_wall`, `KE_s`. Secondary quantities
  (`sig_rr`, `p_b`, `u_B`) are reported but excluded from order claims.
- Mixed terms: $$\Delta\alpha = \alpha_C - \alpha_T$$, $$\Delta T = T_{B,C} - T_{B,T}$$,
  $$\Delta\sigma_{r\theta} = \sigma_{r\theta,C} - \sigma_{r\theta,T}$$ (the three-way term
  needs the radial runs and is NOT-VERIFIED), as the error of the mixed term over the maximum of the
  reference mixed term.

---

## Expected Results

### Numerical acceptance criteria (does the run behave?)

These do not depend on the reference. They are the criteria of the study's audits and are
implemented as layer A of `verification/scripts/verify.py` (see the Verification Workflow).

| Criterion | Threshold | Observed in the accepted runs |
| :--- | :---: | :--- |
| Solver ends at t = 2 s, exit 0, 800 steps, no `FOAM FATAL`, no `nan` | required | all accepted runs |
| Coupling converged on every step (`allowUnconvergedCoupling` unset) | required | 800/800 |
| Maximum outer iterations per step | <= 30 for t >= 0.25 s (start-up exception) | M1: 7 (torsion), 14 (combined), 27 (radial); M2: 22 (torsion), 34 (combined: one step at t = 0.0925 s, inside the start-up exception, and 21 after 0.25 s); M3: 16 (torsion), 25 (combined) |
| Steps accepted by the absolute tolerance | start-up (t < 0.35 s); a later one is CONDITIONAL, not a failure | M1 68 (torsion), 31 (combined), 27 (radial); M2 70, 32; M3 77, 65, all between t = 0.03 and 0.34 s (Table 4) |
| SNES failures | 0 | 0 |
| Final linear residuals | p <= 1e-12, U <= 1e-13 | met for t >= 0.25 s; the combined runs have 51 p solves stopped at maxIter 1000 (final 2e-9 to 7e-9) in three start-up steps (0.0275-0.0325 s), a CONDITIONAL item |
| Continuity error, sum local | <= 1e-16 (torsion); radial and combined <= 2e-16 | 8e-17 (torsion); 1.3e-16 (M3 combined) |
| Minimum cell volume, inverted cells | > 0, none | none |
| Maximum non-orthogonality | <= 18 degrees | 16.5-16.9 |
| Contraction, geometric mean per window | <= 0.90 from t = 1.5 s (CONDITIONAL between 0.5 and 1.5 s) | 0.83-0.90 over the accepted runs (torsion M3 0.835-0.888, combined M3 0.859-0.895); no run of 3 or more steps above 1 |
| Solid $$J_{min}$$ | > 0.99 | 0.998 (torsion M3), 1.000 (combined M3) |
| $$n = 2$$ solid amplitude growth | no exponential growth | maximum amplitude 1.6e-7 m (torsion M3) and 2.1e-7 m (combined M3), 3e-6 of $$u_B$$; log-slope from 1.5 s +3.0 /s and +4.9 /s at every-step sampling, against about +410 /s in the unstable `nNonOrthogonalCorrectors 0` run (which grows exponentially) |
| Angular-momentum and fluid-energy residual per window (relative to the window maximum of the reference torque or interface power) | <= 1e-2 | Met on M3 torsion (6.6e-3 and 7.9e-3). **Exceeded** elsewhere: torsion M1 2.1e-2 and 5.0e-2, M2 1.0e-2 and 1.9e-2; combined M1 4.2e-2 and 7.3e-3, M2 2.9e-2 and 1.06e-2, M3 2.4e-2 and 1.13e-2. The residuals fall with refinement in torsion |

`analyse_case.py` reports the coupling statistics, the balance windows, $$J$$, the principal stretch
and the $$n = 2$$ amplitude; `verify.py` evaluates the criteria above. The mesh-quality, swept-volume and mode-growth
criteria were evaluated by project-level audit scripts on the field output and are not shipped here
(NOT-VERIFIED). FAIL: any row above violated except those marked CONDITIONAL; INCOMPLETE: a truncated or short run, or
a missing `rtDiag*.dat`. The balance-residual row is informational: the 1e-2 threshold is exceeded in most
accepted runs, which is disclosed here rather than hidden. Every accepted run is CONDITIONAL, not PASS, because of
the single-measure and absolute acceptances.

**Table 4: steps per accepting criterion (driver definition, last iteration of each step; 800 steps per run; recomputed from the raw solver logs).** (A) fails the both-measures test = not in `both`; (B) strict-mode reject = both relative measures not satisfied by one measure alone **and** absolute residual not below the absolute tolerance (max(r1, r2) above 1e-6 and no absolute rescue).

| Case | both | r1-only | r2-only | absolute | absolute-accepted times (s) | single-measure share | (A) fail both-measures test | (B) strict-mode reject |
| :-- | :-: | :-: | :-: | :-: | :-: | :-: | :-: | :-: |
| torsion M1 | 32 | 0 | 700 | 68 | 0.0625-0.24 | 87.5% | 96.0% (768) | 86.62% (693) |
| torsion M2 | 47 | 1 | 682 | 70 | 0.1-0.2875 | 85.4% | 94.12% (753) | 84.5% (676) |
| torsion M3 | 56 | 0 | 667 | 77 | 0.14-0.34 | 83.4% | 93.0% (744) | 82.12% (657) |
| combined M1 | 26 | 0 | 743 | 31 | 0.0325-0.11 | 92.9% | 96.75% (774) | 92.25% (738) |
| combined M2 | 26 | 0 | 742 | 32 | 0.0325-0.1175 | 92.8% | 96.75% (774) | 92.12% (737) |
| combined M3 | 20 | 0 | 715 | 65 | 0.0325-0.2 | 89.4% | 97.5% (780) | 88.0% (704) |
| combined M2, `absTol` 1e-10 (F3B) | 17 | 0 | 720 | 63 | 0.0325-0.2 | 90.0% | 97.88% (783) | 88.25% (706) |
| radial M1 | 49 | 0 | 724 | 27 | 0.0325-0.0975 | 90.5% | 93.88% (751) | 89.75% (718) |

`r2-only`: residual2 at or below 1e-6 and residual1 above. `r1-only`: the reverse (one step in torsion M2). No step ended
without a criterion. Test (A), both measures at 1e-6, rejects every step that is not in the `both` column: 93.0-97.5% of the steps of the six torsion and combined runs (97.88% F3B, 93.88% radial M1); only 17-56 of 800 steps meet both measures. Strict mode (B) is less severe because it still accepts the absolute-tolerance steps and the few single-measure steps with an absolute residual below the tolerance: 82.1-92.25%. Counts of the audits that differ (76 and 83 for torsion M1 and M2) are printed-line counts, as
explained in Numerical Set-up.

### Scientific reference-error expectations (how close is the solution to the reference?)

These are errors against the semi-analytical reference, max-norm over 0-2 s, as
fractions of max|ref|, with $$\Delta t = 0.0025$$ s and
`nNonOrthogonalCorrectors 2`. They are expectations, not acceptance limits.

**Torsion** (M1, M2, M3; orders are $$p_{12}, p_{23}$$):

| QoI | M1 | M2 | M3 | $$p_{12}$$, $$p_{23}$$ |
| :-- | :-: | :-: | :-: | :-: |
| rotation $$\alpha_B$$ | 1.17e-2 | 2.97e-3 | 7.6e-4 | 1.98, 1.98 |
| wall torque $$T_B$$ | 2.84e-2 | 7.24e-3 | 1.9e-3 | 1.97, 1.94 |
| shear traction $$\sigma_{r\theta}$$ | 2.76e-2 | 7.04e-3 | 1.9e-3 | 1.97, 1.92 |
| drive torque $$T_A$$ | 2.83e-2 | 6.75e-3 | 1.6e-3 | 2.07, 2.04 |
| fluid angular momentum $$H_f$$ | 2.19e-2 | 5.86e-3 | 1.8e-3 | 1.90, 1.72 |

**Table 5: observed orders $$p_{12}$$, $$p_{23}$$ of the nine primary rotational quantities** (torsion; max-norm and
L2-in-time, 0-2 s and 0-1.8 s; each cell is $$p_{12}$$, $$p_{23}$$):

| QoI | max-norm, 0-2 s | L2, 0-2 s | max-norm, 0-1.8 s | L2, 0-1.8 s |
| :-- | :-: | :-: | :-: | :-: |
| $$\alpha_B$$ | 1.98, 1.98 | 2.01, 2.03 | 1.98, 1.98 | 1.98, 1.98 |
| $$b\omega_B$$ (`w_wall`) | 1.95, 1.82 | 1.92, 1.81 | 1.95, 1.89 | 1.96, 1.89 |
| $$T_B$$ | 1.97, 1.94 | 2.00, 2.00 | 1.94, 1.85 | 2.01, 2.00 |
| $$\sigma_{r\theta}$$ | 1.97, 1.92 | 2.00, 2.00 | 1.94, 1.83 | 2.01, 2.00 |
| $$T_A$$ | 2.07, 2.04 | 2.00, 2.04 | 1.99, 2.02 | 1.99, 1.99 |
| $$H_s$$ | 1.95, 1.95 | 1.95, 1.89 | 1.95, 1.95 | 1.97, 1.97 |
| $$H_f$$ | 1.90, 1.72 | 1.91, 1.71 | 1.90, 1.72 | 1.90, 1.71 |
| $$P_{wall}$$ | 1.96, 1.93 | 1.99, 2.00 | 1.96, 1.93 | 2.00, 2.01 |
| $$KE_s$$ | 1.93, 1.93 | 1.97, 2.01 | 1.93, 1.93 | 1.97, 2.01 |

The max-norm range is 1.72-2.07 and the L2 range 1.71-2.04 (1.68 is the dt-corrected $$H_f$$, a different quantity).
In the max-norm, the values outside 1.8-2.0 are $$H_f$$ 1.72 (both windows), $$T_A$$ 2.07 and 2.04 (0-2 s) and 2.02
(0-1.8 s); $$b\omega_B$$ (1.82), $$T_B$$ (1.85) and $$\sigma_{r\theta}$$ (1.83) are inside. These are observed orders;
an asymptotic order is not established. The secondary quantities converge more slowly ($$u_B$$ 1.1 and 0.7,
$$\sigma_{rr}$$ 1.2 and 0.8, $$p_b$$ 0.8 and 0.5) and are excluded.

**Table 6: observed max-norm orders with the first-order uncertainty interval at the documented noise** (nominal and, in brackets, nominal $$\pm$$ the noise propagated into the error ratio; noise from the `numerical_noise` section of `verification/thresholds.json`: 2.4e-4 of max|ref| for the torque-type torsion quantities and far less for $$\alpha_B$$, $$H_f$$; intervals with half-width above 0.15 are marked noise-limited, those below 0.02 sharp; an interval is unresolved when the noise is not below the error, which the driver treats as a downgrade; L2-in-time intervals are not quantified because the noise is documented only in the max-norm). Torsion:

| QoI | noise (of max\|ref\|) | 0-2 s: $$p_{12}$$ | $$p_{23}$$ | 0-1.8 s: $$p_{12}$$ | $$p_{23}$$ |
| :-- | :-: | :-: | :-: | :-: | :-: |
| $$\alpha_B$$ | 6.9e-07 | 1.98 (1.98-1.98 sharp) | 1.98 (1.97-1.98 sharp) | 1.98 (1.98-1.98 sharp) | 1.98 (1.97-1.98 sharp) |
| $$b\omega_B$$ | 0.00024 | 1.95 (1.86-2.03) | 1.82 (1.51-2.12 noise-limited) | 1.95 (1.83-2.07) | 1.89 (1.45-2.33 noise-limited) |
| $$T_B$$ | 0.00024 | 1.97 (1.92-2.02) | 1.94 (1.75-2.12 noise-limited) | 1.94 (1.88-2.00) | 1.85 (1.63-2.07 noise-limited) |
| $$\sigma_{r\theta}$$ | 0.00024 | 1.97 (1.92-2.02) | 1.92 (1.73-2.11 noise-limited) | 1.94 (1.87-2.00) | 1.83 (1.60-2.07 noise-limited) |
| $$T_A$$ | 0.00024 | 2.07 (2.02-2.12) | 2.04 (1.82-2.25 noise-limited) | 1.99 (1.94-2.05) | 2.02 (1.79-2.24 noise-limited) |
| $$H_s$$ | 0.00024 | 1.95 (1.87-2.03) | 1.95 (1.64-2.26 noise-limited) | 1.95 (1.87-2.03) | 1.95 (1.64-2.26 noise-limited) |
| $$H_f$$ | 8.6e-07 | 1.90 (1.90-1.90 sharp) | 1.72 (1.72-1.72 sharp) | 1.90 (1.90-1.90 sharp) | 1.72 (1.72-1.72 sharp) |
| $$P_{wall}$$ | 0.00024 | 1.96 (1.90-2.01) | 1.93 (1.72-2.14 noise-limited) | 1.96 (1.90-2.01) | 1.93 (1.72-2.14 noise-limited) |
| $$KE_s$$ | 0.00024 | 1.93 (1.88-1.97) | 1.93 (1.76-2.11 noise-limited) | 1.93 (1.88-1.97) | 1.93 (1.76-2.11 noise-limited) |

Combined (noise 1.3e-3 for the torque-type quantities):

| QoI | noise (of max\|ref\|) | 0-2 s: $$p_{12}$$ | $$p_{23}$$ | 0-1.8 s: $$p_{12}$$ | $$p_{23}$$ |
| :-- | :-: | :-: | :-: | :-: | :-: |
| $$\alpha_B$$ | 5.7e-06 | 1.30 (1.30-1.30 sharp) | 0.72 (0.72-0.73 sharp) | 1.30 (1.30-1.30 sharp) | 0.72 (0.72-0.73 sharp) |
| $$b\omega_B$$ | 0.00011 | 1.11 (1.09-1.13 sharp) | 0.64 (0.61-0.67) | 1.07 (1.05-1.10) | 0.61 (0.57-0.65) |
| $$T_B$$ | 0.0013 | 1.85 (1.54-2.16 noise-limited) | 0.57 (0.04-1.10 noise-limited) | 1.76 (1.45-2.07 noise-limited) | 0.57 (0.04-1.10 noise-limited) |
| $$\sigma_{r\theta}$$ | 0.0012 | 1.71 (1.45-1.97 noise-limited) | 0.43 (0.01-0.85 noise-limited) | 1.60 (1.33-1.86 noise-limited) | 0.43 (0.01-0.85 noise-limited) |
| $$T_A$$ | 0.00089 | 1.82 (1.60-2.03 noise-limited) | 0.99 (0.52-1.46 noise-limited) | 1.73 (1.50-1.95 noise-limited) | 0.95 (0.48-1.42 noise-limited) |
| $$H_s$$ | 1.2e-05 | 1.06 (1.06-1.07 sharp) | 0.62 (0.61-0.62 sharp) | 1.13 (1.13-1.13 sharp) | 0.71 (0.71-0.71 sharp) |
| $$H_f$$ | 1.1e-05 | 1.63 (1.63-1.64 sharp) | 1.11 (1.11-1.12 sharp) | 1.63 (1.63-1.64 sharp) | 1.11 (1.11-1.12 sharp) |
| $$P_{wall}$$ | 0.00064 | 0.40 (0.38-0.42) | 0.18 (0.16-0.21) | 0.43 (0.41-0.45) | 0.30 (0.27-0.33) |
| $$KE_s$$ | 0.0013 | 1.16 (1.00-1.33 noise-limited) | 0.28 (0.05-0.51 noise-limited) | 1.16 (1.00-1.33 noise-limited) | 0.28 (0.05-0.51 noise-limited) |

The intervals qualify the orders. They can lower the confidence in an order, never turn the observation into a proof of an asymptotic order, and the driver never converts them into a PASS.

**Combined** (sharply resolved near-second-order self-convergence only for $$\alpha_B$$ and $$H_f$$; $$T_B$$, $$\sigma_{r\theta}$$ and $$T_A$$ have nominal orders near 2 with broad noise-limited intervals; the floor statement covers the same five quantities):

| QoI | M1 | M2 | M3 | self order |
| :-- | :-: | :-: | :-: | :-: |
| rotation $$\alpha_B$$ | 1.27e-2 | 5.14e-3 | 3.1e-3 | 1.97 |
| wall torque $$T_B$$ | 2.28e-2 | 6.33e-3 | 4.3e-3 | 2.00 |
| shear traction $$\sigma_{r\theta}$$ | 2.28e-2 | 6.95e-3 | 5.2e-3 | 1.97 |
| drive torque $$T_A$$ | 2.16e-2 | 6.14e-3 | 3.1e-3 | 2.09 |

The M2-to-M3 reduction is far below a factor of four: for these four quantities (and $$H_f$$) the errors at M3 sit on a
floor of about 2.5-5.2e-3 of max|ref|, so the observed order against the reference
is 0.43-1.11 there (0.43, 0.57, 0.72, 0.99, 1.11 for $$\sigma_{r\theta}$$, $$T_B$$, $$\alpha_B$$, $$T_A$$, $$H_f$$). This is the disclosed floor, not second-order convergence to the
reference. The normal-traction quantities are worse and non-monotone in order
($$u_B$$ 1.1%, $$\sigma_{rr}$$ 0.85%, $$p_b$$ 2.4% at M3).

**Table 7: combined configuration by quantity (0-2 s, max-norm errors as fractions of max|ref|).**

| QoI | M1 | M2 | M3 | self-order | floor statement |
| :-- | :-: | :-: | :-: | :-: | :-- |
| $$\alpha_B$$ | 1.27e-02 | 5.14e-03 | 3.11e-03 | 1.97 | inside (2.5-5.2e-3 floor) |
| $$b\omega_B$$ | 2.10e-02 | 9.75e-03 | 6.27e-03 | 1.66 | outside: above the floor |
| $$T_B$$ | 2.28e-02 | 6.33e-03 | 4.26e-03 | 2.00 | inside (2.5-5.2e-3 floor) |
| $$\sigma_{r\theta}$$ | 2.28e-02 | 6.95e-03 | 5.15e-03 | 1.97 | inside (2.5-5.2e-3 floor) |
| $$T_A$$ | 2.16e-02 | 6.14e-03 | 3.09e-03 | 2.09 | inside (2.5-5.2e-3 floor) |
| $$H_s$$ | 2.21e-02 | 1.06e-02 | 6.91e-03 | 1.63 | outside: above the floor |
| $$H_f$$ | 2.05e-02 | 6.61e-03 | 3.06e-03 | 1.94 | inside (2.5-5.2e-3 floor) |
| $$P_{wall}$$ | 7.17e-02 | 5.43e-02 | 4.78e-02 | 0.70 | outside: above the floor |
| $$KE_s$$ | 2.87e-02 | 1.28e-02 | 1.06e-02 | 2.17 | outside: above the floor |
| $$u_B$$ (excluded) | 1.78e-02 | 1.27e-02 | 1.06e-02 | 0.78 | outside: above the floor |
| $$\sigma_{rr}$$ (excluded) | 1.48e-02 | 1.07e-02 | 8.50e-03 | 0.62 | outside: above the floor |
| $$p_b$$ (excluded) | 2.62e-02 | 2.49e-02 | 2.43e-02 | 0.67 | outside: above the floor |

Sharply resolved near-second-order self-convergence belongs only to $$\alpha_B$$ (1.97, interval 1.96-1.97) and $$H_f$$ (1.94, 1.93-1.94). $$T_B$$ (2.00, 1.39-2.60), $$\sigma_{r\theta}$$ (1.97, 1.39-2.54) and $$T_A$$ (2.09, 1.63-2.55) have nominal orders near 2 but broad noise-propagated intervals, and the floor statement covers these five quantities. $$KE_s$$ (2.17, 1.42-2.92) is also noise-limited and its reference error (1.06e-2) is above the floor; $$H_s$$, $$b\omega_B$$ and $$P_{wall}$$ do not self-converge at second order. The driver prints this table and a status per quantity (`CONDITIONAL` or `DISCLOSED` inside the statement, `ABOVE-FLOOR` outside; never PASS).

**Mixed terms** (error over signal; M1, M2, M3, 0-2 s): $$\Delta\alpha$$ 2.2%,
1.8%, 1.5%; $$\Delta T$$ 8.1%, 3.9%, 2.5%; $$\Delta\sigma_{r\theta}$$ 7.2%, 3.1%,
2.2% (3.4% at M3 over 0-1.8 s).

**Radial only, M1** (single mesh, no convergence claim): $$u_B$$ 1.70e-2,
$$\sigma_{rr}$$ 1.37e-2, $$p_b$$ 2.58e-2 of max|ref| over 0-2 s.

### Regression test

`regressionTest.sh` runs `./Allrun M1 combined serial endTime=0.3` in a copy of the tutorial (an `endTime` override, so
the dictionaries are regenerated with the shortened horizon): 120 time-steps, about 35 s including the compilation of the
coded conditions. It checks that the run completes, that the interface displacement $$u_B$$ lies within 3e-5 m of the
reference (the M1 mesh gives 8e-6 m) and that $$u_B$$ and $$\alpha_B$$ equal the stored values to 1e-7 m and 2e-8 rad. The
stored values reproduce the study's accepted M1 combined run to round-off. It checks launch and reproducibility only and
not any verification claim.

**Interaction with `tutorials/Alltest`.** `Alltest` cuts every case to one step by editing `controlDict`, `fvSolution`
and `fsiProperties` and then runs `./Allrun`. The default `./Allrun` therefore uses the **stored** dictionaries and
writes only the reference solution; it regenerates them only on request (`regenerate`, a mesh, a configuration,
`nNonOrth=`, `absTol=`, `endTime=`, or `np=N` other than 4 in parallel). The solver entries of
`system/fluid/fvSolution` are written one per line because `Alltest` replaces `maxIter.*;` greedily. Tested: the
default `./Allrun` after the `Alltest` edits ends after one step (`End`, 21 s), and an overridden run (`endTime=0.05`)
runs the 20 steps requested; `verification/test_verify.py` also checks the option logic without a solver. A run
with fewer than 10 diagnostic rows skips the analysis.

---

## Known Caveats

- M1 and M2 results are from serial runs, M3 results from 4-rank runs. The scotch partition is not reproducible
  between `decomposePar` invocations (the accepted M3 runs used different partitions, and so did two fresh
  decompositions of the same mesh), so parallel runs agree with each other only at the level of the iteration and
  partition noise, not bitwise. From existing evidence this is about 2e-4 of max|ref| for the M3 torsion torque
  (restart and attempt differences 1.8e-4 and 2.4e-4) and about 1.3e-3 for the combined torque and shear (serial
  against 4 ranks 1.26e-3 and 1.21e-3). Serial M1 runs reproduce bitwise. The numerical-noise section of
  `verification/thresholds.json` records these values; they cap marginal results at CONDITIONAL and are
  never used to widen a scientific band.
- M3 was run on 4 ranks with `absTol=1e-10`. A 1e-12 m absolute tolerance stalled at start-up.
- **Persistent viscous normal wall-stress discrepancy.** The CFD wall viscous normal stress (`sig_nv`, from the
  `pimpleFluid` deviatoric stress at the moving wall) is 0.286, 0.307 and 0.317 of the reference value on the M1, M2 and
  M3 combined meshes (0.295 on radial M1), a ratio near one third that does not converge with the mesh. `analyse_case.py`
  therefore prints a `sig_nv` error of about 0.70, and the note in its output says so. It is a documented discrepancy of the
  fluid treatment of the moving wall, not a bug of the analysis; the normal-traction quantities ($$u_B$$, $$\sigma_{rr}$$,
  $$p_b$$) inherit it and are excluded from every order claim. A reference run with a transmitted fraction of
  about one third reproduces about one third of the rotation excess of the combined case; this is partial
  evidence, not a proven unique mechanism, and the combined reference-error floor stays disclosed.
- The reference-error floor in the combined configuration is not removed by mesh refinement at this time step.
- The earlier diagnostics function objects (serial-only, without the $$n = 3, 4$$ columns) were used by the
  M2 accepted cases only (`F2_M2T`, `F3_M2C`, `F3B`, `F3_M2R`); the M1 and M3 cases used the version stored here. The shared
  columns are the same quantities.
- The ALE/GCL consecutive-save checks, the mode-growth checks and any temporal-order study are outside this
  tutorial.

---

## Deferred Diagnostics (designs only; nothing here was run)

These designs address the origin of the combined reference-error floor. None was run, and none changes a claim above.

- **(a) Fluid-only, prescribed wall motion.** The wall trajectory of each reference (or the coupled CFD trajectory) is imposed
  on the fluid alone, with the same meshes, schemes and boundary conditions, and the fluid error at given motion is measured with no
  feedback, then pushed through the reduced coupled model as an interface-load perturbation (closure). It would show
  whether the floor is born in the fluid (discretisation, outer open boundary, wall traction) or in the rest of the coupled path.
  It would **not** separate the coupling error from the solid-side error, does not test the fluid-to-solid traction transfer,
  the coupled temporal error or the M3 absolute-tolerance effect, cannot reach the rotation mixed term except through the
  closure, and cannot resolve effects below about 2e-4 (smoothed) or 1.3e-3 (max-norm torque of the combined runs). One flaw of the design: its smallest-N
  criterion for the PIMPLE passes is not computable from N = 2, 12, 24, so N = 2, 4, 8, 16, 32 or a residual control is needed.
- **(b) Solid only, combined, driven by the reference fluid tractions.** It would measure the solid-side discretisation error of the
  combined case on its own and so split the closure remainder of (a) into solid and coupling parts. It would **not**
  test the fluid, the transfer or the iteration.
- **(c) Larger outer radius c = 3.** It would test the sensitivity of the floor to the outer open boundary and
  to the truncation of the domain at r = c (a pressure or normal effect at r = c, which the tangential part of (a) cannot
  address). It needs a new reference for the new geometry and would **not** identify a fluid discretisation error.

(The design document labels solid-only as test (c) and c = 3 as test (b).) Suggested order (design audit): (b) solid-only first, because it is cheap and splits what (a) cannot; then (a) with its
closure; (c) last. The partial radial-traction explanation (about one third of the rotation excess, from the
transmitted fraction of the viscous normal stress) is partial evidence and is not presented as a proven unique mechanism.
