# womersleyTube: interface-traction diagnostic (optional)

This opt-in diagnostic measures each part of the fluid–solid interface traction
in the womersleyTube FSI tutorial and compares it with the exact linear
solution. It also splits the tutorial's attenuation error into temporal,
spatial, amplitude-dependent and residual parts.

- It is separate from the standard verification (`../Allverify`), which it does
  not change, and it does not change the solver or the tutorial defaults.
- Its function object (`scripts/interfaceTractionDiagnostic`) is not part of
  the tutorial's `system/`; it is included only in the study copies it
  creates.

## Requirements

- OpenFOAM.com with solids4foam built with PETSc (MUMPS and HYPRE), as for
  the tutorial.
- python3 with numpy; matplotlib is optional (plots are skipped without it).
- The `coded` function object is compiled at run time (the usual OpenFOAM
  `dynamicCode` requirements). It is serial only and stops with a fatal error
  in a decomposed run.

## Usage

```bash
./AllverifyTraction --quick          # one m1/50 run + consistency checks
./AllverifyTraction --cores N        # full study (about 4 CPU-hours)
./AllverifyTraction --analysis-only  # re-analyse the runs in work/
./AllverifyTraction --force ...      # discard stale study copies and runs
```

Runs are made in copies of the tutorial under `work/`; results are written to
`postProcessing/`. Both are git-ignored and removed by `../Allclean` (and so
`./Allclean all`).

- A study copy records hashes of the files the solver reads, the analysis
  files and its settings. A copy that no longer matches is prepared again; if
  it holds runs, the script stops with exit status 3 unless `--force` is
  given.
- A completed run is reused when its stored key matches: case, settings,
  tutorial files and the content hashes of the `solids4Foam` executable and
  libraries. Editing a README or the analysis scripts does not invalidate
  runs.
- Before any analysis (also with `--analysis-only`) every required run is
  validated; a missing or invalid run stops the script with exit status 1.
  `--skip-build-check` omits only the solver-build comparison.

Outputs: `interfaceTraction_cases.csv`, `traction_cases_all.csv`,
`traction_convergence.{json,md}`, `findings_summary.{json,txt}`,
`amplitude_series_m1.json`, `traction_checks.txt` and plots.

## What is measured

`scripts/interfaceTractionDiagnostic` is a `coded` function object (region
`fluid`) that records, every time step over the last two periods:

<!-- markdownlint-disable MD013 -->
| Field | Definition | Purpose |
|---|---|---|
| `p_b` | boundary pressure on `interface` | pressure traction, as used by `patchPressureForce` and the momentum equation |
| `tTr` | `n & -devReff_b` from the registered turbulence model | the viscous traction `pimpleFluid::patchViscousForce` transfers (before ρ) |
| `tLap`, `tExp` | `gradientInternalCoeffs*U_P + gradientBoundaryCoeffs`, `n & dev2(T(grad(U)_b))` | the momentum equation's own viscous wall flux |
| `trGradU` | `tr(grad(U)_b)` | continuity of the boundary velocity gradient |
| `tApplied` | `traction - n pressure` of the solid `solidTraction` patch | the traction applied to the solid |
| `tSigma` | `n & sigma_b` (solid) | the reconstructed boundary stress, for reference only |
<!-- markdownlint-enable MD013 -->

`scripts/womersley_interface_traction.py` takes the first temporal Fourier
coefficient of each face value and projects it on the exact travelling wave
over the central half of the tube, `â = Σ a(x) e^{ikx} w / Σ w`, using the
undeformed solid face positions and face areas. It compares the result with
the exact `p`, `2μ ∂u_r/∂r` and `σ_rx` at `r = R_i` from
`../scripts/womersley_exact.py`. Integrated (modal) quantities are used
because the finite-volume wall flux is first order pointwise but second order
in its integral ([#583][i583]).

## Study

<!-- markdownlint-disable MD013 -->
| Copy | Runs | Purpose |
|---|---|---|
| `mu5` | standard study (m1/m2/m4 at 200 steps/period; m2 at 50/100/200/400; IQN-ILS and Robin) plus m4/400 | spatial and temporal error |
| `amp05`, `amp2`, `amp4` | m2/200 (`amp4` also m1/200) | amplitude-dependent part and the mesh dependence of its slope |
| `tol` | m2/200, coupling tolerances × 0.01 | coupling-iteration error |
| `ampm1`, `p10`, `p10a4` (optional) | m2/200: a = −1; 10 periods at a = 1 and 4 | parity of the amplitude dependence; does it decay with time |
<!-- markdownlint-enable MD013 -->

At 4× amplitude the standard driver's 5 % check on the profile sample points
fails because the wall moves further; `amplitude_series.py` relaxes that check
only, to 50 %.

The optional copies are not run by `AllverifyTraction`. Prepare them with
`scripts/prepare_study.py --ensure ../.. work/ampm1/womersleyTube
--amplitude=-1` (and `--periods 10`, `--periods 10 --amplitude 4` for
`p10`, `p10a4`), run them with `scripts/run_single.py
work/<copy>/womersleyTube robin:2:200`, and analyse them with
`scripts/amplitude_series.py work OUT.json` and
`scripts/periods_slope.py work OUT.json`.

## Pass/fail checks

These cover consistency only:

- transferred viscous traction = momentum wall flux (≤ 1e-10);
- applied solid traction = transferred fluid traction (normal ≤ 1e-4 of |p|,
  shear ≤ 1e-8);
- the shear traction converges with observed order ≥ 1.5 and a Richardson
  ratio close to a positive real (|arg| ≤ 30°).

Everything else is reported, not asserted. Observed orders are used for
extrapolation only if p > 0.3 and |arg ratio| ≤ 30°; otherwise p = 2 is used
and the series is flagged.

## Results

OpenFOAM v2512, GCC, PETSc 3.24, this branch merged with development
`e123ecf75` (which includes [#671][pr671]), Robin–Neumann, μ = 5 Pa s.
`./AllverifyTraction --cores 20` ran end-to-end with exit status 0 in
1 h 57 min, and every check passed.

<!-- markdownlint-disable MD013 -->
| Quantity | Value | Basis |
|---|---|---|
| Transferred viscous traction vs momentum wall flux | ≤ 3.3e-14 (all runs) | measured |
| Applied solid vs transferred fluid traction | shear ≤ 7e-12 of \|σ_rx\|; normal 1.5e-6 of \|p\| in every run, unchanged (4e-11) with tolerances × 0.01 | measured |
| Normal viscous stress / exact 2μ∂u_r/∂r | 0.55 − 0.35i, 0.78 − 0.19i, 0.89 − 0.10i on m1/m2/m4 | measured |
| | → 1.01 + 0.00i (observed order 0.96) | extrapolated |
| Shear traction σ_rx | observed order 1.99; limit within 5e-5 of exact | measured order; extrapolated limit |
| Pressure | observed order 1.98; limit within 6e-5 of \|p\| | measured order; extrapolated limit |
| Total normal traction error | 1.2e-3 of \|p\| at m4/400; observed order 1.99; limit within 1.1e-4 of \|p\| | measured; extrapolated limit |
| Effect of coupling tolerances × 0.01 | ≤ 8e-7 of \|p\| on every traction component | measured |
<!-- markdownlint-enable MD013 -->

The extrapolated limits combine space (observed order), time (q = 2, from
m4/200 and m4/400) and the amplitude slope (a → 0).

The normal viscous stress converges to the exact value at first order. Before
[#671][pr671] the wall velocity conditions set ∂u_n/∂n = 0, so the transferred normal
viscous stress converged to about one third of the exact value (0.20, 0.27,
0.30 on m1/m2/m4); see [#670][i670]. With the continuity-based gradient,
`tr(grad(U)_b)` is zero to round-off by construction.

`tSigma` differs from the applied traction by 31 %, 16 % and 8 % of |p| on
m1/m2/m4 (first order). It is reconstructed from the cell gradient, whereas
the applied traction enters the solid residual directly; read interface loads
from the applied traction.

## Attenuation error budget

The attenuation (Im k) error at m4, 200 steps/period is −2.05e-3:

<!-- markdownlint-disable MD013 -->
| Part | Value | Basis |
|---|---|---|
| temporal | −0.67e-3 | m4/200 vs m4/400 with q = 2 |
| spatial | unresolved, \|·\| ≤ 0.18e-3 | m1/m2/m4 series non-monotone (−2.10, −2.23, −2.05e-3) |
| amplitude-dependent | −1.10e-3 | linear fit over a = 0.5–4 at m2/200 (exponent 1.00–1.01); the m1 slope differs by 0.5 % |
| residual | −0.22e-3 ± 0.20e-3 | remainder |
<!-- markdownlint-enable MD013 -->

The amplitude-dependent part is linear and odd in the amplitude, and does not
decay with time:

- a = −1 (m2/200): the odd part (E(1) − E(−1))/2 is 1.00 times the fitted
  slope (0.99–1.04 over the six quantities) and the even part equals the
  fit's intercept to within 0.04 times the slope. A periodic steady state of
  the symmetric problem would give E(−1) = E(1).
- 10 periods (analysing periods 9–10): the slope is 1.002 times that over
  periods 5–6 (0.993–1.004), so it is not a decaying start-up transient.

Its origin is not identified. A persistent effect of the initial state, such
as the static offset of the initial solid displacement carried by the fluid
mesh ([#650][i650]), would fit, but these runs do not test it.

On m2 the attenuation error is −4.28e-3, −2.23e-3 and −1.73e-3 at 100, 200
and 400 steps/period (observed order 2.05).

## Limitations

- One geometry, one solid model (linear geometry) and one platform. The
  Robin–Neumann coupling was studied in full; IQN-ILS only at m2/100.
- The normal viscous stress is only about 5.5e-4 of the pressure in this tube.
- The amplitude slope is assumed independent of Δt (not checked) and the
  error parts are assumed additive.

[i583]: https://github.com/solids4foam/solids4foam/issues/583
[i650]: https://github.com/solids4foam/solids4foam/issues/650
[i670]: https://github.com/solids4foam/solids4foam/issues/670
[pr671]: https://github.com/solids4foam/solids4foam/pull/671
