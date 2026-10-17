# womersleyTube verification summary

Errors are relative to the exact linear solution unless stated.

| Case | profile | flow_amp | flow_phase | wallMid_amp | wallMid_phase | speed | attenuation | wallQuarter_amp | wallThreeQuarter_amp | axialMid_amp | periodicity | Iter./step | Time (s) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| iqnils_m2_n100 | 4.98e-03 | 5.48e-04 | 1.09e-03 | -4.47e-04 | -2.40e-03 | -3.72e-05 | -4.41e-03 | -7.44e-04 | 5.62e-04 | 2.57e-03 | 1.50e-05 | 14.56 | 524 |
| monolithic_m1_n200 | 1.76e-02 | 6.08e-03 | 6.08e-03 | -1.04e-02 | -2.89e-03 | -2.35e-03 | -2.80e-03 | -7.03e-03 | -6.03e-03 | 1.17e-02 | 3.01e-05 | 2.00 | 98 |
| monolithic_m2_n50 | 6.22e-03 | -1.27e-03 | 1.39e-04 | 4.94e-03 | -6.59e-03 | 1.34e-03 | -1.35e-02 | 1.40e-03 | 5.64e-03 | 3.58e-03 | 6.38e-05 | 2.00 | 77 |
| monolithic_m2_n100 | 4.70e-03 | 4.73e-04 | 1.07e-03 | -9.96e-04 | -2.82e-03 | -1.36e-04 | -4.64e-03 | -1.24e-03 | 1.82e-04 | 2.94e-03 | 1.21e-05 | 2.00 | 149 |
| monolithic_m2_n200 | 4.69e-03 | 9.01e-04 | 1.32e-03 | -2.52e-03 | -1.96e-03 | -5.19e-04 | -2.58e-03 | -1.95e-03 | -1.18e-03 | 2.77e-03 | 8.46e-06 | 2.00 | 304 |
| monolithic_m2_n400 | 4.68e-03 | 1.01e-03 | 1.39e-03 | -2.90e-03 | -1.76e-03 | -6.16e-04 | -2.09e-03 | -2.13e-03 | -1.52e-03 | 2.72e-03 | 1.08e-05 | 2.00 | 595 |
| monolithic_m4_n200 | 9.97e-04 | -3.42e-04 | 6.30e-05 | -5.30e-04 | -1.62e-03 | -5.96e-05 | -2.31e-03 | -6.29e-04 | 1.05e-05 | 5.37e-04 | 6.91e-06 | 2.00 | 1056 |

## Checks

- PASS: iqnils_m2_n100: last two periods differ by 1.50e-05 (within 0.0005)
- PASS: monolithic_m1_n200: last two periods differ by 3.01e-05 (within 0.0005)
- PASS: monolithic_m2_n50: last two periods differ by 6.38e-05 (within 0.0005)
- PASS: monolithic_m2_n100: last two periods differ by 1.21e-05 (within 0.0005)
- PASS: monolithic_m2_n200: last two periods differ by 8.46e-06 (within 0.0005)
- PASS: monolithic_m2_n400: last two periods differ by 1.08e-05 (within 0.0005)
- PASS: monolithic_m4_n200: last two periods differ by 6.91e-06 (within 0.0005)
- PASS: finest mesh (monolithic_m4_n200): velocity profile at L/2, max |error| / max |u| 9.97e-04 within 0.003
- PASS: finest mesh (monolithic_m4_n200): flow rate amplitude at L/2 -3.42e-04 within 0.003
- PASS: finest mesh (monolithic_m4_n200): flow rate phase at L/2 (rad) 6.30e-05 within 0.003
- PASS: finest mesh (monolithic_m4_n200): wall radial displacement amplitude at L/2 -5.30e-04 within 0.003
- PASS: finest mesh (monolithic_m4_n200): wall radial displacement phase at L/2 (rad) -1.62e-03 within 0.005
- PASS: finest mesh (monolithic_m4_n200): wave speed from the pressure along the tube -5.96e-05 within 0.003
- PASS: finest mesh (monolithic_m4_n200): attenuation rate, Im(k) -2.31e-03 within 0.01
- PASS: mesh: observed order 2.23 of the profile error at least 1.5
- PASS: mesh: observed order 2.06 of the flow_amp error at least 1.5
- PASS: mesh: observed order 1.98 of the wallMid_amp error at least 1.5
- PASS: mesh: observed order 2.00 of the speed error at least 1.5
- PASS: time step: observed order 2.01 of the flow_amp (finest three steps) at least 1.7
- PASS: time step: observed order 1.94 of the flow_phase (finest three steps) at least 1.7
- PASS: time step: observed order 1.99 of the wallMid_amp (finest three steps) at least 1.7
- PASS: time step: observed order 2.07 of the wallMid_phase (finest three steps) at least 1.7
- PASS: time step: observed order 1.98 of the speed (finest three steps) at least 1.7
- PASS: time step: observed order 2.05 of the attenuation (finest three steps) at least 1.7
  - time error at 200 steps per period: wall amplitude 3.84e-04, wave speed 9.71e-05
- PASS: IQN-ILS (iqnils_m2_n100): velocity profile at L/2, max |error| / max |u| 4.98e-03 within 0.006
- PASS: IQN-ILS (iqnils_m2_n100): wall radial displacement amplitude at L/2 -4.47e-04 within 0.003
- PASS: IQN-ILS (iqnils_m2_n100): wave speed from the pressure along the tube -3.72e-05 within 0.003
- PASS: quasi-monolithic (monolithic_m2_n100): velocity profile at L/2, max |error| / max |u| 4.70e-03 within 0.006
- PASS: quasi-monolithic (monolithic_m2_n100): wall radial displacement amplitude at L/2 -9.96e-04 within 0.003
- PASS: quasi-monolithic (monolithic_m2_n100): wave speed from the pressure along the tube -1.36e-04 within 0.003
- PASS: quasi-monolithic and IQN-ILS agree on the velocity profile at L/2, max |difference| / max |u| to 7.21e-04 (within 0.001)
- PASS: quasi-monolithic and IQN-ILS agree on the flow rate amplitude at L/2 to 7.49e-05 (within 0.001)
- PASS: quasi-monolithic and IQN-ILS agree on the flow rate phase at L/2 (rad) to 2.54e-05 (within 0.001)
- PASS: quasi-monolithic and IQN-ILS agree on the wall radial displacement amplitude at L/2 to 5.50e-04 (within 0.001)
- PASS: quasi-monolithic and IQN-ILS agree on the wall radial displacement phase at L/2 (rad) to 4.17e-04 (within 0.001)
- PASS: quasi-monolithic and IQN-ILS agree on the wave speed from the pressure along the tube to 9.90e-05 (within 0.001)
- PASS: quasi-monolithic and IQN-ILS agree on the attenuation rate, Im(k) to 2.25e-04 (within 0.001)
  - mean coupling iterations per step: 14.56 (IQN-ILS), 2.00 (quasi-monolithic); run time 524 s and 149 s

## Observed orders

Unsigned errors (profile): from the errors of the two finest runs. Signed errors: from the differences between the three runs, which cancel the error shared by the series, such as the time error of the mesh study. Time-step series of more than three runs give one order per three successive steps.

| Series | Quantity | Errors | Order |
|---|---|---|---:|
| mesh | profile | 1.76e-02, 4.69e-03, 9.97e-04 | 2.23 |
| mesh | flow_amp | 6.08e-03, 9.01e-04, -3.42e-04 | 2.06 |
| mesh | flow_phase | 6.08e-03, 1.32e-03, 6.30e-05 | 1.92 |
| mesh | wallMid_amp | -1.04e-02, -2.52e-03, -5.30e-04 | 1.98 |
| mesh | wallMid_phase | -2.89e-03, -1.96e-03, -1.62e-03 | 1.46 |
| mesh | speed | -2.35e-03, -5.19e-04, -5.96e-05 | 2.00 |
| mesh | attenuation | -2.80e-03, -2.58e-03, -2.31e-03 | -0.30 |
| time step | flow_amp | -1.27e-03, 4.73e-04, 9.01e-04, 1.01e-03 | 2.03, 2.01 |
| time step | flow_phase | 1.39e-04, 1.07e-03, 1.32e-03, 1.39e-03 | 1.86, 1.94 |
| time step | wallMid_amp | 4.94e-03, -9.96e-04, -2.52e-03, -2.90e-03 | 1.96, 1.99 |
| time step | wallMid_phase | -6.59e-03, -2.82e-03, -1.96e-03, -1.76e-03 | 2.14, 2.07 |
| time step | speed | 1.34e-03, -1.36e-04, -5.19e-04, -6.16e-04 | 1.94, 1.98 |
| time step | attenuation | -1.35e-02, -4.64e-03, -2.58e-03, -2.09e-03 | 2.11, 2.05 |

PASSED
