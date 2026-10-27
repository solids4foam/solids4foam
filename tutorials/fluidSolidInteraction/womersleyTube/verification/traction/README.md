# womersleyTube: interface-traction diagnostic (optional)

This opt-in diagnostic measures each part of the fluid–solid interface traction
in the womersleyTube FSI tutorial and compares it with the exact linear
solution. It also decomposes the tutorial's attenuation error.

- It is separate from the standard verification (`../Allverify`), which it does
  not change.
- It does not change the solver or the tutorial defaults.
- Its function object is included only in the study copies it creates.

## Requirements

- OpenFOAM.com with solids4foam built with PETSc, which needs MUMPS and HYPRE
  (the tutorial's requirement). Tested with OpenFOAM v2512, GCC 13.3, PETSc
  3.19.6 and solids4foam development `1576edcc4`.
- python3 with numpy. matplotlib is optional: without it, the plots are
  skipped.
- The `coded` function object is compiled at run time (the usual OpenFOAM
  `dynamicCode` requirements).

## Usage

```bash
./AllverifyTraction --quick            # one m1/50 run + consistency checks (~1 min)
./AllverifyTraction --cores 6          # full study (~2-3 CPU-hours)
./AllverifyTraction --analysis-only    # re-analyse existing runs in work/
```

Runs are made in `work/`, and results are written to `postProcessing/`. Both
are git-ignored.

Outputs:

- `interfaceTraction_cases.csv`, `traction_cases_all.csv`: per run;
- `traction_convergence.{json,md}`;
- `findings_summary.{json,txt}`;
- `traction_checks.txt`;
- plots, all generated from the JSON/CSV files.

## What is measured

`system/interfaceTractionDiagnostic` is a `coded` function object (region
`fluid`) that records, every time step over the last two periods:

<!-- markdownlint-disable MD013 -->
| Field | Definition | Purpose |
|---|---|---|
| `p_b` | boundary pressure on `interface` | pressure traction, as used by `patchPressureForce` and the momentum equation |
| `tTr` | `n & -devReff_b` from the registered turbulence model | exactly the viscous traction `pimpleFluid::patchViscousForce` transfers (before ρ) |
| `tLap`, `tExp` | `gradientInternalCoeffs*U_P + gradientBoundaryCoeffs`, `n & dev2(T(grad(U)_b))` | the momentum equation's own viscous wall flux |
| `trGradU` | `tr(grad(U)_b)` | `-2μ tr` estimates the continuum `2μ ∂u_n/∂n` from the fluid fields |
| `tApplied` | `traction - n pressure` of the solid `solidTraction` patch, read from the patch's own written entries | the traction applied to the solid |
| `tSigma` | `n & sigma_b` (solid) | the reconstructed boundary stress, for reference only |
<!-- markdownlint-enable MD013 -->

`scripts/womersley_interface_traction.py`:

1. takes the first temporal Fourier coefficient of each face value;
2. projects it on the exact travelling wave over the central half of the tube,
   `â = Σ a(x) e^{ikx} Δx / Σ Δx`;
3. compares the result with the exact `p`, `2μ ∂u_r/∂r` and `σ_rx` at
   `r = R_i` from `../scripts/womersley_exact.py`.

Integrated (modal) quantities are used because the finite-volume wall flux is
conservative: pointwise first order, second order in its integral (solids4foam
issue #583).

## Study

<!-- markdownlint-disable MD013 -->
| Copy | Runs | Purpose |
|---|---|---|
| `mu5` | standard study (m1/m2/m4 at 200 steps/period; m2 at 50/100/200/400; IQN-ILS and Robin) plus m4/400 | spatial and temporal error; non-intrusiveness |
| `amp05`, `amp2`, `amp4` | m2/200 (`amp4` also m1/200) | finite-amplitude part; mesh dependence of its slope |
| `tol` | m2/200, coupling tolerances × 0.01 | coupling-iteration error |
<!-- markdownlint-enable MD013 -->

The amplitude variants regenerate `constant/womersleyProperties` and the
reference JSON with the copy's own exact script. At 4× amplitude the standard
driver's check that the profile sample points stay within 5 % of the cell
spacing fails, because the wall moves further. `amplitude_series.py` relaxes
that check, and only that check, to 50 %. The axial pressure samples, from
which the wave number is obtained, are unaffected.

## Pass/fail checks

These cover consistency only:

- transferred viscous traction = momentum wall flux (≤ 1e-10);
- applied solid traction = transferred fluid traction (normal ≤ 1e-4 of |p|;
  shear ≤ 1e-8; a constant 1.5e-6 normal difference is observed);
- the shear traction converges with observed order ≥ 1.5.

Everything else is reported, not asserted.

## Results (OpenFOAM v2512, development 1576edcc4, Robin–Neumann, μ = 5 Pa s)

The table separates directly measured values from those obtained by
extrapolation.

<!-- markdownlint-disable MD013 -->
| Claim | Value | Basis |
|---|---|---|
| Transferred viscous traction = momentum wall flux | ≤ 5e-14 (all runs) | measured |
| Applied solid traction = transferred fluid traction | ≤ 5e-12 shear; 1.5e-6 of \|p\| normal, the same in every run (all meshes, time steps, amplitudes, both couplings, and with coupling tolerances × 0.01): not a coupling residual, origin not identified, negligible for every QoI | measured |
| Transferred normal viscous / exact 2μ∂_r u_r | 0.198, 0.270, 0.302 on m1/m2/m4 (200 steps/period) | measured |
| | → 0.332 (observed spatial order 1.06; 0.313 with p = 2) | extrapolated |
| Fields' own −2μ tr(∇U_b) / exact | 0.59, 0.81, 0.91 | measured |
| | → 0.996 | extrapolated |
| Shear traction σ_rx | observed order 2.00; limit within 4e-5 of exact | measured order; extrapolated limit |
| Total normal traction (pressure + viscous) error | 1.0e-3 of \|p\| at m4/400 | measured |
| | −1.2e-4 − 0.4e-4i of \|p\| | extrapolated (space p = 2.00, time q = 2 (observed 1.99), amplitude slope) |
| Uncompensated deficit implied by a one-third viscous part | 0.36e-4 + 3.68e-4i of \|p\| | computed |
| Residual along the deficit direction | −0.13 of the deficit (the pressure over-compensates by about 13 % of the deficit); 1.1e-4 of \|p\| across it, origin unidentified | extrapolated |
<!-- markdownlint-enable MD013 -->

**Interpretation.**

- The wall condition sets the normal component of the wall-normal velocity
  gradient to zero. With `devTwoSymm`, the transferred normal viscous stress is
  therefore (2/3)μ∂u_n/∂n instead of 2μ∂u_n/∂n.
- The momentum equation uses the same boundary flux, so this is a consistent,
  conservative split.
- The wall pressure absorbs the difference: in the extrapolated limit the total
  normal traction agrees with the exact solution along the deficit direction to
  about 13 % of the deficit.
- This is a statement about how the normal traction is split between pressure
  and viscous parts. **No error in the transferred total traction is
  demonstrated.**

## Attenuation error budget (m4, 200 steps/period: −1.961e-3)

<!-- markdownlint-disable MD013 -->
| Part | Value | Basis |
|---|---|---|
| temporal error at 200 steps/period | −0.67e-3 | measured m4/200 vs m4/400 with q = 2 |
| spatial error | −0.07e-3 | m1/m2/m4 series non-monotone (−2.08, −2.18, −1.96e-3): bounded by ±0.22e-3 |
| finite-amplitude part | −1.09e-3 | linear fit over a = 0.5–4 at m2/200 (observed exponent 1.00, misfit 4e-6); slope on m1 within 2 % |
| residual | −0.13e-3 ± 0.24e-3 | remainder |
<!-- markdownlint-enable MD013 -->

The finite-amplitude part is an **estimated model-form contribution**: the
difference between the nonlinear problem solved at the tutorial amplitude and
the linear reference. Assumptions:

- the slope measured at m2/200 holds at m4 and in the Δt → 0 limit (checked
  for m1 vs m2; not checked in Δt);
- the error parts are additive.

Its physical mechanism is not identified. The mean wall displacement is O(a²)
and only 2.7e-4 of the first harmonic.

On m2 (the standard study's time series) the attenuation error is −4.24, −2.18,
−1.69e-3 at 100/200/400 steps/period (time order 2.07). It is mesh-independent
at fixed 200 steps/period, but not time-step-independent.

## Limitations

- The normal viscous stress is only about 5.5e-4 of the pressure in this tube,
  so the compensation is tested against a small signal.
- A μ = 20 Pa s variant (deficit 4.6× larger) **did not converge**: Robin
  m4/200, Robin m2/50 and IQN-ILS m4/200 diverged in the first time step. Its
  m1/m2 results (two mesh levels) are not conclusive about compensation.
- The spatial part of the attenuation error is not resolved (non-monotone
  series).
- `tSigma` differs from the applied traction by 31 %, 16 % and 8 % of \|p\| on
  m1/m2/m4, i.e. it converges at first order (ratio ≈ 2 per refinement). It is
  a reconstruction of the boundary stress from the cell gradient: the applied
  traction enters the solid residual directly as a face force. Use the applied
  traction, not boundary `sigma`, for interface loads. No defect is implied.
- Only one geometry, one solid model (linear geometry) and the Robin–Neumann
  coupling were studied in full. IQN-ILS was compared at m2/100 only.
