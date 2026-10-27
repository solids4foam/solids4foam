# womersleyTube: interface-traction diagnostic (optional)

This opt-in diagnostic measures every part of the fluid–solid interface traction in the womersleyTube FSI tutorial and compares it with the exact linear solution. It is separate from the standard verification (`../Allverify`). It changes neither the solver nor the tutorial defaults. Its function object is included only in the study copies it creates.

```bash
./AllverifyTraction --quick            # one m1/50 run and the consistency checks (about 1 min)
./AllverifyTraction --cores 6          # full study (about 2–3 CPU-hours)
./AllverifyTraction --analysis-only    # re-analyse existing runs in work/
```

Outputs go to `postProcessing/`: CSV/JSON tables, convergence plots, `findings_summary.json` and `traction_checks.txt`.

## What is measured

The function object `system/interfaceTractionDiagnostic` (a `coded` function object, region `fluid`, compiled in each case's `dynamicCode`) writes, every time step over the last two periods:

| Quantity | Definition | Why |
|---|---|---|
| `p_b` | boundary pressure on `interface` | pressure traction, as used by `patchPressureForce` and the momentum equation |
| `tTr` | `n & -devReff_b` from the registered turbulence model | exactly the viscous traction `pimpleFluid::patchViscousForce` transfers (before the factor ρ) |
| `tLap`, `tExp` | `gradientInternalCoeffs*U_P + gradientBoundaryCoeffs` and `n & dev2(T(grad(U)_b))` | the momentum equation's own viscous wall flux, which the fluid actually loses |
| `trGradU` | `tr(grad(U)_b)` | `−2 μ tr` estimates the continuum normal viscous stress `2 μ ∂u_n/∂n` from the fields |
| `tApplied` | `traction − n pressure` of the solid `solidTraction` patch, read from the patch's own written entries | the traction the solid condition applies |
| `tSigma` | `n & sigma_b` (solid) | reference only. With the PETSc SNES solid, the boundary value of `sigma` is **not** the enforced traction (about 8–30 % different in these runs) |

`scripts/womersley_interface_traction.py` takes the first temporal Fourier coefficient of each face value. It then projects it on the exact travelling wave over the central half of the tube, `â = Σ a(x) e^{ikx} Δx / Σ Δx`. This is an integrated (modal) quantity, consistent with the conservative finite-volume wall flux (see solids4foam issue #583: the wall traction is first order pointwise but second order in its integral). The script compares the result with the exact `p`, `2μ ∂u_r/∂r` and `σ_rx` at `r = R_i` from `../scripts/womersley_exact.py`.

## Study

| Copy | Runs | Purpose |
|---|---|---|
| `mu5` | the standard study (m1/m2/m4 at 200 steps/period; m2 at 50/100/200/400; IQN-ILS vs Robin) plus m4 at 400 | spatial and temporal error; non-intrusiveness (results identical to the standard study) |
| `amp05`, `amp2`, `amp4` | m2/200 (`amp4` also m1/200) | finite-amplitude part of the error, and whether its slope depends on the mesh |
| `tol` | m2/200 with coupling tolerances × 0.01 | coupling-iteration error |

Each amplitude variant regenerates `constant/womersleyProperties` and the reference JSON with the copy's own exact-solution script, so the case and the reference stay consistent.

## Checks (pass/fail)

These cover consistency only:

- transferred viscous traction = momentum wall flux (≤ 1e-10);
- applied solid traction = transferred fluid traction (normal ≤ 1e-4 of |p|, coupling-tolerance level; shear ≤ 1e-8);
- the shear traction converges with observed order ≥ 1.5.

The normal-viscous ratio, the Richardson limits and the error budget are reported, not asserted.

## Findings (OpenFOAM v2512, development 1576edcc4)

- **Transfer chain:** consistent. The transferred traction equals the momentum flux to about 1e-14. The applied traction equals the transferred traction to the coupling tolerance (1.5e-6 normal; ≤ 5e-12 shear; 3e-11 with tolerances × 0.01).
- **Normal viscous part:** the transferred normal viscous stress converges to **one third** of `2μ ∂u_r/∂r` (space–time limit 0.332 − 0.002i). The fields' own `∂u_n/∂n` converges to the exact value (0.996). The one third follows from the wall condition, which sets the normal component of the wall-normal velocity gradient to zero, and `devTwoSymm`: (2/3)μ ∂u_n/∂n is transferred instead of 2μ ∂u_n/∂n.
- **Total normal traction:** the wall pressure compensates this split. At zero amplitude the total normal traction (pressure + viscous) converges to the exact value within about 1.2e-4 of |p|, against an uncompensated deficit of 3.7e-4 of |p|. The wave QoIs show no effect of the split. The shear traction converges at second order to the exact value.
- **Attenuation error:** the −1.96e-3 attenuation error at m4/200 decomposes as:
  - −0.67e-3 temporal error at 200 steps/period;
  - −0.07e-3 spatial error;
  - −1.09e-3 a finite-amplitude part, linear in the amplitude (observed exponent 1.00) and mesh-independent (m1 and m2 slopes within 2 %);
  - −0.13e-3 ± 0.15e-3 residual.

## Limitations

- **Single tube:** the normal viscous stress is about 5.5e-4 of the pressure here, so the compensation is resolved only to about 30 % of the deficit. A μ = 20 variant (deficit 4.6 times larger) diverged at start-up on the m4 mesh with both couplings, so that check has two mesh levels only.
- **Mesh series:** the spatial part of the attenuation error comes from a non-monotone mesh series.
- **Amplitude:** the amplitude slope is measured at 200 steps/period and assumed independent of the time step.
- **Mechanism:** the finite-amplitude mechanism is not identified. The mean wall displacement is O(a²) and too small to explain it.
