---
sort: 36
---

# Pressure-Stabilisation Fourier Verification: `pressureStabilisationFourier`

---

Prepared by Aaron Mullen-Hales

---

## Tutorial Aims

- Verify the pressure-stabilisation operators of the mixed
  displacement-pressure solid models against their analytical Fourier
  symbols, with and without spectral normalisation.
- Check the parameter validation of the spectral normalisation.

## Case Overview

`stabilisationFourierCheck` constructs each runtime-selected
`stabilisationModel` from `system/stabilisationFourierCheckDict` and applies it,
through the same calls a solid model makes (`updateScalar(p, &gradP)` and
`cellScalar(&rAUf, true)`), to discrete Fourier modes of the pressure on a
uniform periodic mesh: a unit square with cyclic x and y boundaries and one
empty z cell. `rAUf` is 1. The measured eigenvalue of each mode is compared
with the analytical symbol, written independently in
`applications/utilities/stabilisationFourierCheck/fourierSymbols.H`.

The models are `laplacian`, `JamesonSchmidtTurkel`,
`generalisedEvenOrderLaplacian` with m = `laplacianPower` = 0, 1 and 2,
`diffStencilLaplacian` and `RhieChow`, all with `scaleFactor 0.37`. The order
m is counted from zero: m = 0 is the Laplacian, m = 1 the bi-Laplacian and
m = 2 the sixth-order operator. Each check runs twice: with the legacy (raw)
operators, and with `spectralNormalisation true` and
`referenceNyquistDirections 2`. The
pairs that must coincide (`laplacian` and m = 0, `JamesonSchmidtTurkel` and
m = 1, `diffStencilLaplacian` and `RhieChow`) are also compared. The pressure
gradient is `Gauss linear`, face interpolation `linear` and normal gradients
`uncorrected`; these choices are part of the `diffStencilLaplacian`/`RhieChow`
symbol.

A 3-D level repeats the checks on an 8^3 mesh of the unit cube, cyclic in x, y
and z (`system/blockMeshDict.3D`), with `referenceNyquistDirections 3` and
modes in all three directions (`system/stabilisationFourierCheckDict.3D`).

The ten `system/reject*` dictionaries are invalid configurations, each of
which must be refused with its message: `referenceNyquistDirections` zero,
missing, or larger than the number of geometric mesh directions (2 here);
normalisation of models that do not support it (`alpha`, `gradDiv`,
`volStrainRate`, `multiple`); a negative `laplacianPower`; and both
`spectralNormalisation` and its deprecated alias `normalise` in one model.
The output of each refused run is kept in `expectedError.<dictionary>`.
The two `system/accept*` dictionaries must be accepted: `spectralNormalisation`
without a warning, and `normalise` with a deprecation warning.

## Running the Case

```bash
./Allrun           # 16 cells per axis, then the 3-D level
./Allrun refine    # 16, 32 and 64 cells per axis, then the 3-D level
```

`Allrun` stops if any eigenvalue or pair check fails (the utility exits
nonzero) or if an invalid dictionary is accepted. The results are written to
`postProcessing/stabilisationFourierCheck`:

- `results_NxNx1.csv` (`results_8x8x8.csv` for the 3-D level): measured and
  analytical eigenvalues per model and mode, for the raw
  (`spectralNormalisation` 0) and normalised (1) campaigns. The
  `paperOrder` column is the number of repeated Laplacian factors, m + 1, and
  the `spectralRelation` labels (`legacyM1`, `coupledRepeatedLaplacian`,
  `splitM2Shape`) also count the order from one;
- `crossChecks_NxNx1.csv` and `gammaChecks_NxNx1.csv`: pair and constant-gamma
  checks;
- `symbols_NxNx1.csv` and `symbolMaps_NxNx1.csv`: analytical cuts and maps
  (`-symbolSweep`).

With gnuplot, `plotSymbols.gnuplot` writes `symbols_N16.png`: the analytical
symbols along theta_y = 0 (lines) and the measured eigenvalues (points), each
scaled by its value at the Nyquist mode.

## Expected Results

All eigenvalue and pair checks pass to the tolerances in
`system/stabilisationFourierCheckDict` (relative eigenvalue error 1e-9), and all
ten invalid configurations are refused, and both spellings of the switch are
accepted as described above, on the 2-D levels and on the 3-D level.
`regressionTest.sh` checks all of these.

## Notes

The case is skipped on foam-extend. The utility needs a uniform, axis-aligned
mesh with cyclic and empty patches only.
