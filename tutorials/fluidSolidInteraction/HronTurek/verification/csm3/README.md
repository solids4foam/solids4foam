# CSM3 structural verification of the Hron-Turek flag

Structural-only study that isolates the spatial convergence of the solids4foam
solid used in the Hron-Turek FSI3 case, using the Turek-Hron CSM benchmarks
(Featflow, "CSM tests"). It is separate from the FSI3 driver and is not run by
any regression suite.

| | CSM2 | CSM3 |
|---|---|---|
| geometry | the FSI flag: x in [0.24899, 0.6] m, y in [0.19, 0.21] m, 2D plane strain | same |
| support | `plateFix` (x = 0.24899 m) clamped; other faces traction free | same |
| material | St. Venant-Kirchhoff, rho = 1000, nu = 0.4, E = 5.6e6 (mu = 2e6, **the FSI3 plate**) | same, E = 1.4e6 (mu = 0.5e6) |
| load | gravity g = (0, -2, 0) m/s^2 on the flag only | same |
| definition | steady | transient from rest, undeformed start |
| formulation | `nonLinearGeometryTotalLagrangianTotalDisplacement`, PETSc SNES | same, `backward` ddt/d2dt2 |
| QoIs | `u_x`, `u_y` of A = (0.6, 0.2) | mean, amplitude, frequency of `u_x`, `u_y` of A |
| Featflow reference | `-0.469000`, `-16.9739` mm (level 5+1) | `u_x -14.305 +- 14.305`, `u_y -63.607 +- 65.160` mm, f = 1.0995 Hz (dt = 0.005 s, level 4+0) |

The CSM3 reference depends on its own time step (`u_x` mean -14.40, -14.65,
-14.30 mm for dt = 0.02, 0.01, 0.005 s; frequency 1.0956 to 1.0995 Hz), so it
is a 1-2% reference. The steady CSM1/CSM2 values are converged to about 5
digits. CSM1 was not run: the single-step Newton solve stalls at 2x and above
for the soft material.

## Files

- `case/` template (CSM3: transient, E = 1.4e6). `overlay/` holds the variants
  (steady, E = 5.6e6, stabilisation factor, high-order reconstruction, tight
  SNES tolerance) that `run_case.sh` copies over it.
- `run_case.sh <name> <nx> <ny> <dt> <endTime> [overlay]`, `submit.sh` (Slurm).
- `analyse.py` statistics of the tip history; `summarise.py` builds `results/`;
  `pointcheck.py` compares the monitored point with the written `pointD`.
- `results/` compact CSV/JSON: `csm3_levels.csv`, `csm3_observed_order.csv`,
  `csm3_dt_check.csv`, `csm2_static_levels.csv`,
  `csm2_static_observed_order.csv`, `csm2_static_anisotropic.csv`,
  `csm3_structural_convergence.json`.

Mesh level `kx` is 105k x 6k cells in a single block, so h = 3.343/k mm.
