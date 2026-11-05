# womersleyTube: interface-traction diagnostic (optional)

This opt-in diagnostic measures each part of the fluid–solid interface traction
in the womersleyTube FSI tutorial and compares it with the exact linear
solution. It also decomposes the tutorial's attenuation error.

- It is separate from the standard verification (`../Allverify`), which it does
  not change.
- It does not change the solver or the tutorial defaults.
- Its function object (`scripts/interfaceTractionDiagnostic`) is not part of
  the tutorial's `system/`; it is installed and included only in the study
  copies it creates.

## Requirements

- OpenFOAM.com with solids4foam built with PETSc, which needs MUMPS and HYPRE
  (the tutorial's requirement). The end-to-end run below used OpenFOAM v2512, GCC 11.4, PETSc
  3.24 and solids4foam development `1576edcc4`; earlier development used GCC 13.3 and
  PETSc 3.19.6.
- python3 with numpy. matplotlib is optional: without it, the plots are
  skipped.
- The `coded` function object is compiled at run time (the usual OpenFOAM
  `dynamicCode` requirements). It is serial only: it stops with a fatal error
  in a decomposed run, because every rank would write the same files with
  processor-local face indices.

## Usage

```bash
./AllverifyTraction --quick            # one m1/50 run + consistency checks (~1 min)
./AllverifyTraction --cores 6          # full study (~4 CPU-hours, about 1.8 h on 6 cores)
./AllverifyTraction --analysis-only    # re-analyse existing runs in work/
./AllverifyTraction --force ...        # discard stale study copies (and their runs)
```

Each study copy in `work/` is built in a staging directory (`<name>.partial`),
every copied file and edit is validated, and `diagnostic_settings.json` is
written last with a hash of what the solver reads, a hash of the analysis files
(`verification/scripts/*.py`, reference JSON), the diagnostic and
`prepare_study.py`, and the requested settings. The staging directory replaces
the old copy by renaming only: the old copy is renamed aside and deleted after
success, so an interrupted preparation leaves a complete copy (an absent one,
for an instant between the two renames, is restored from the renamed one at the
next call). An existing copy is kept if the record is complete and matches. If
only the analysis files changed it is rebuilt and its runs are moved into the
new copy (they are reused only if their solver-run key, below, still matches).
A copy that differs in what the solver reads, or lacks a complete record, is
prepared again if it holds no runs; if it holds runs the script stops (exit
status 3) and keeps them unless `--force` is given. `--force` only re-prepares
such a copy (discarding its runs); it does not rerun valid runs elsewhere.

Runs are made in `work/`, and results are written to `postProcessing/`. Both
are git-ignored by the repository's root `.gitignore` and removed by
`../Allclean` (and so `./Allclean all`).

Completed runs in `work/` are reused, not repeated, when their stored
solver-run key matches: case, mesh, physical parameters, time schedule,
content hashes of the study copy's tutorial files, a digest of the case inputs
as the driver configures them, and content hashes (not size or time stamps) of
the `solids4Foam` executable and `libsolids4Foam*` libraries; and the solver
log shows a completed run and the diagnostic output is non-empty. A run that
does not match is run again. Editing a README, the analysis code or the exact
values of the reference JSON does not change this key. The analysis key (content
hashes of the verification driver, the reference JSON and these scripts) is
written with the solver-run keys to `postProcessing/analysis_provenance.json`.
Runs are made by `run_single.py` without analysing them, so one failed case does
not prevent the remaining runs.

Before any analysis (also with `--analysis-only`) `scripts/check_runs.py`
applies the same validation to every required run (stored settings, tutorial
files, solver build, log, diagnostic output), and stops the script with a list
of the affected runs and exit status 1 if one is missing or not valid; no
partial analysis is made. `--skip-build-check` omits only the build
comparison, for analysing stored runs where the solver environment is not
available.

The full study was run end-to-end as the single command
`./AllverifyTraction --cores 6` (OpenFOAM v2512, PETSc 3.24, solids4foam
`1576edcc4`): exit status 0, 14,300 CPU-seconds (3.97 CPU-hours), 1 h 47 min
wall clock on 6 cores. A second `./AllverifyTraction --analysis-only` started
no solver run. The tutorial files were those of this branch before the commit
that updated this text.

Outputs:

- `interfaceTraction_cases.csv`, `traction_cases_all.csv`: per run;
- `traction_convergence.{json,md}`;
- `findings_summary.{json,txt}`;
- `amplitude_series_m1.json`: amplitude slope on m1 (a = 1 and 4);
- `traction_checks.txt`;
- plots, all generated from the JSON/CSV files.

## What is measured

`scripts/interfaceTractionDiagnostic` is a `coded` function object (region
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
   `â = Σ a(x) e^{ikx} w / Σ w`, with the face areas `magSf` of the solid patch
   as weights `w` and its undeformed face positions `x` (linear geometry: the
   solid mesh does not move). Solid faces are matched one-to-one to fluid
   faces on the full face-centre coordinates. The moving fluid mesh is not
   used: it oscillates axially by ~7e-4 m, and its window-mean position carries
   the static offset of the initial solid displacement (the fluid mesh moves by
   D(t) − D(0)); `--fluid-mean-positions` uses it to quantify the effect;
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
| `amp05`, `amp2`, `amp4` | m2/200 (`amp4` also m1/200) | amplitude-dependent part; mesh dependence of its slope |
| `ampm1` (optional) | m2/200, a = −1 | odd and even part of the amplitude dependence |
| `p10`, `p10a4` (optional) | m2/200, 10 periods, a = 1 and 4 | does the amplitude slope decay with time |
| `tol` | m2/200, coupling tolerances × 0.01 | coupling-iteration error |
<!-- markdownlint-enable MD013 -->

The optional copies are not run by `AllverifyTraction`; `ampm1`, `p10` and
`p10a4` are prepared into `work/` with `scripts/prepare_study.py --ensure
TUTORIAL work/ampm1/womersleyTube --amplitude=-1`, `... work/p10/womersleyTube
--periods 10` and `... work/p10a4/womersleyTube --periods 10 --amplitude 4`,
run with `scripts/run_single.py work/<copy>/womersleyTube robin:2:200` (about
13, 21 and 21 CPU-min), and analysed with `scripts/amplitude_series.py work
OUT.json` (adds the odd/even parts when `ampm1` is present) and
`scripts/periods_slope.py work OUT.json`. `--periods N` keeps the analysis window at the
last two periods.

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
- the shear traction converges with observed order ≥ 1.5 and a ratio
  (E1−E2)/(E2−E4) close to a positive real (|arg| ≤ 30°; otherwise the series
  is non-asymptotic and the check fails).

Everything else is reported, not asserted. Observed orders of the complex
errors are used for extrapolation only if p > 0.3 and |arg ratio| ≤ 30°
(`scripts/richardson.py`); otherwise the p = 2 limit is used and the series is
flagged.

## Results (OpenFOAM v2512, development 1576edcc4, Robin–Neumann, μ = 5 Pa s)

**Status.** The tables below were checked against the output of the end-to-end
run described in Usage, which used the current analysis code (projection on the
solid's undeformed positions with face-area weights; Richardson order accepted
only if p > 0.3 and |arg ratio| <= 30 deg). Every number in the tables was
reproduced to the digits given; the shear-order check passes (order 2.00,
arg(ratio) -2.3 deg). The numbers are in `postProcessing/findings_summary.txt`,
`traction_convergence.md`, `amplitude_series_m1.json` and `traction_checks.txt`;
to reproduce them, run:

```bash
./AllverifyTraction --cores 6           # or --analysis-only if the runs are in work/
```

The table separates directly measured values from those obtained by
extrapolation.

<!-- markdownlint-disable MD013 -->
| Claim | Value | Basis |
|---|---|---|
| Transferred viscous traction = momentum wall flux | ≤ 5e-14 (all runs) | measured |
| Applied solid traction = transferred fluid traction | ≤ 5e-12 shear; 1.5e-6 of \|p\| normal, the same in every run (all meshes, time steps, amplitudes, both couplings, and with coupling tolerances × 0.01): not a coupling residual, origin not identified, negligible for every QoI | measured |
| Transferred normal viscous / exact 2μ∂_r u_r | 0.198, 0.270, 0.303 on m1/m2/m4 (200 steps/period) | measured |
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
- The wall pressure can therefore absorb the difference. This pressure
  compensation is physically plausible but **not directly established** by
  these calculations:
  - at the computed resolutions the total normal traction error (1.0e-3 of
    \|p\| at m4/400) is larger than the deficit (3.7e-4 of \|p\|), so the
    compensation cannot be observed there;
  - the near-agreement applies only to the extrapolated limit (space, time and
    amplitude). There the residual is nonzero: about 1.2e-4 of \|p\|, roughly
    one third of the deficit. Along the deficit direction it is −0.13 of the
    deficit; 1.1e-4 of \|p\| lies across it, origin unidentified.
- This is a statement about how the normal traction is split between pressure
  and viscous parts. **No error in the transferred total traction is
  demonstrated.**

## Attenuation error budget (m4, 200 steps/period: −1.961e-3)

Checked against the end-to-end run, as for the results above.

<!-- markdownlint-disable MD013 -->
| Part | Value | Basis |
|---|---|---|
| temporal error at 200 steps/period | −0.67e-3 | measured m4/200 vs m4/400 with q = 2 |
| spatial error | −0.07e-3 | m1/m2/m4 series non-monotone (−2.08, −2.18, −1.96e-3): bounded by ±0.22e-3 |
| amplitude-dependent part | −1.09e-3 | linear fit over a = 0.5–4 at m2/200 (observed exponent 1.00, misfit 4e-6); the m1 slope (a = 1 and 4) differs by 2.0 % |
| residual | −0.13e-3 ± 0.24e-3 | remainder |
<!-- markdownlint-enable MD013 -->

The amplitude-dependent part is an **amplitude-dependent contribution, linear
in the tested amplitude range, odd in the sign of the amplitude and not
decaying over 10 periods; origin unidentified** (see the two tests below).
Assumptions:

- the slope measured at m2/200 holds at m4 and in the Δt → 0 limit (checked
  for m1 vs m2; not checked in Δt);
- the error parts are additive.

Why its origin is open:

- The forcing is a pure harmonic. For a smooth periodic response, reversing
  the sign of the amplitude is equivalent to a half-period time shift, so the
  leading relative nonlinear correction to the first harmonic is expected to
  be O(a²).
- The observed, approximately linear dependence on positive amplitude is not
  explained by this simple smooth periodic-response argument. The mean wall
  displacement is O(a²) and only 2.7e-4 of the first harmonic.
- The two tests below exclude a decaying start-up mode as the cause of the
  amplitude-linear part, and show that it is odd in the sign of the amplitude.
  They do not identify the mechanism.

**Sign of the amplitude (a = -1, m2/200, 6 periods).** For a smooth periodic
steady state the response to −a is the response to a shifted by half a period,
so E(−1) = E(1) would follow. Observed (attenuation error, fit of a = 0.5–4:
E0 = −1.087e-3, slope s = −1.0905e-3 per Pa): the odd part (E(1) − E(−1))/2 =
−1.0854e-3 = 1.00 s and the even part (E(1) + E(−1))/2 = −1.0911e-3 =
E0 + 0.00 s. All six QoIs behave the same (odd part = 0.99–1.03 s; even part
equal to E0 within 0.03 s). The error is therefore linear **and odd** in a: the
half-period symmetry of a periodic steady state is broken, so the error is not
a property of the periodic steady state of the symmetric problem.

**Duration (10 periods, analysis periods 9–10, m2/200, a = 1 and 4).** If the
slope came from a free mode excited at start-up and decaying, it would fall with
the number of periods. The two-point slope (E(4) − E(1))/3 of the attenuation
error is −1.0912e-3 per Pa over periods 5–6 and −1.0951e-3 per Pa over periods
9–10 (ratio 1.004; the other five QoIs 0.993–1.003), and E(1) is −2.1765e-3 and
−2.1819e-3. The amplitude-linear part does not decay over four more periods
(a decay on a time scale much longer than 10 periods is not excluded).

**What remains.** An error that is odd in a and persistent is compatible with
an effect of the initial state that does not decay, such as the static offset of
the initial solid displacement carried by the moving fluid mesh (the fluid mesh
moves by D(t) − D(0), closed issue #650, linear in the amplitude), but these
runs do not show that this is the cause, nor its size. A run with the offset
removed, or a smaller-amplitude series, would test it. The free-mode
explanation proposed earlier is not supported (slope ratio 1.004).

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
- The full driver was run once end-to-end (see Usage), on one machine; an
  independent rerun has not been made.
