# Figures of the tutorial README

<!-- markdownlint-disable MD013 -->

| File | Content | Source table | Script |
| :-- | :-- | :-- | :-- |
| `fig1_geometry.png` / `.svg` | geometry and loading (schematic), drive curves | `../constant/timeVsAngle`, `../constant/pcTable` | `fig1_geometry.py` |
| `fig2_rotational_convergence.png` | alpha_B and T_B errors against mesh (A, B), time step and strict coupling | `data/fig2_convergence.csv`, `data/fig2_timestep_strict.csv` | `fig2_rotational_convergence.py` |
| `fig3_normal_traction.png` | normal wall traction: reference against CFD, pressure and viscous error budget | `data/fig3_traction_timeseries.csv` | `fig3_normal_traction.py` |

<!-- markdownlint-enable MD013 -->

Regenerate the figures from the shipped tables (matplotlib and numpy only, no
CFD):

```bash
cd figures
python3 fig1_geometry.py .. fig1_geometry
python3 fig2_rotational_convergence.py data fig2_rotational_convergence
python3 fig3_normal_traction.py data fig3_normal_traction
```

The tables are extracted from run diagnostics with `make_figure_data.py DIAG_DIR
../verification/scripts .`, where each `DIAG_DIR/<case>/` holds
`rtDiagFluid.dat`, `rtDiagSolid.dat` and `reference.json` (per-rank files
recombined for parallel runs). The runs are not bundled.

## Formulations and data provenance

* **A** is the tutorial as stored: default `leastSquares` for `grad(p)` and
  `grad(U)`, stock `newMovingWallVelocity`. Runs: torsion F2_M1T, F2_M2T,
  F1_M3T; combined F3_M1C, F3_M2C, F3_M3C (4 ranks on M3).
* **B** is A with `grad(p) Gauss linear` only (optional variant): combined
  C1_M1C, C2_M2C, C3_M3C.
* New controls (formulation A, combined M2, two ranks): P0 (2.5 ms), P1 (1.25
  ms), P2 (0.625 ms), P3 (strict coupling, `requireAllResidualMeasures yes`,
  `outerCorrTolerance 1e-7`; aborted at t = 1.825 s, first 1.8225 s used).
* Errors are max-norm over 0-2 s divided by max|reference|; finer-step runs are
  compared on the common 2.5 ms grid. The bars in Figure 2 are the documented
  partition/iteration noise of `verification/thresholds.json`.

## Limitations

A single time step per mesh in panels (a, b): no temporal order is claimed, and
order lines are drawn only where the observed order is resolved. Formulation B
is not clean in numerical acceptance at M3 (Layer A FAIL) and is not the
tutorial default. Panel (c) concerns combined M2 under A only. Figure 3 uses one
mesh (M2) and does not use any research wall condition.
