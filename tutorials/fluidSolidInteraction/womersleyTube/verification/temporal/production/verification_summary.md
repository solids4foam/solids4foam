# womersleyTube verification summary

Errors are relative to the exact linear solution unless stated.

| Case | profile | flow_amp | flow_phase | wallMid_amp | wallMid_phase | speed | attenuation | wallQuarter_amp | wallThreeQuarter_amp | axialMid_amp | periodicity | Iter./step | Time (s) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| iqnils_m2_n100 | 4.98e-03 | 5.48e-04 | 1.09e-03 | -4.47e-04 | -2.40e-03 | -3.72e-05 | -4.41e-03 | -7.44e-04 | 5.62e-04 | 2.57e-03 | 1.50e-05 | 14.56 | 805 |
| robin_m1_n200 | 1.79e-02 | 6.10e-03 | 6.05e-03 | -1.03e-02 | -2.37e-03 | -2.33e-03 | -2.08e-03 | -6.84e-03 | -6.09e-03 | 1.12e-02 | 2.78e-05 | 12.76 | 339 |
| robin_m2_n50 | 6.18e-03 | -1.23e-03 | 1.43e-04 | 5.28e-03 | -6.12e-03 | 1.41e-03 | -1.31e-02 | 1.75e-03 | 5.81e-03 | 3.17e-03 | 6.61e-05 | 10.25 | 226 |
| robin_m2_n100 | 4.98e-03 | 5.14e-04 | 1.07e-03 | -6.63e-04 | -2.34e-03 | -6.62e-05 | -4.24e-03 | -9.00e-04 | 3.46e-04 | 2.53e-03 | 1.42e-05 | 9.25 | 381 |
| robin_m2_n200 | 4.96e-03 | 9.39e-04 | 1.33e-03 | -2.19e-03 | -1.48e-03 | -4.50e-04 | -2.18e-03 | -1.61e-03 | -1.02e-03 | 2.36e-03 | 6.11e-06 | 10.96 | 771 |
| robin_m2_n400 | 4.95e-03 | 1.04e-03 | 1.40e-03 | -2.58e-03 | -1.28e-03 | -5.48e-04 | -1.69e-03 | -1.80e-03 | -1.37e-03 | 2.31e-03 | 7.88e-06 | 9.46 | 1127 |
| robin_m4_n200 | 1.20e-03 | -3.09e-04 | 8.23e-05 | -1.76e-04 | -1.15e-03 | 1.40e-05 | -1.96e-03 | -2.80e-04 | 1.96e-04 | 1.45e-04 | 7.35e-06 | 8.38 | 2239 |

## Checks

- PASS: iqnils_m2_n100: last two periods differ by 1.50e-05 (within 0.0005)
- PASS: robin_m1_n200: last two periods differ by 2.78e-05 (within 0.0005)
- PASS: robin_m2_n50: last two periods differ by 6.61e-05 (within 0.0005)
- PASS: robin_m2_n100: last two periods differ by 1.42e-05 (within 0.0005)
- PASS: robin_m2_n200: last two periods differ by 6.11e-06 (within 0.0005)
- PASS: robin_m2_n400: last two periods differ by 7.88e-06 (within 0.0005)
- PASS: robin_m4_n200: last two periods differ by 7.35e-06 (within 0.0005)
- PASS: finest mesh (robin_m4_n200): velocity profile at L/2, max |error| / max |u| 1.20e-03 within 0.003
- PASS: finest mesh (robin_m4_n200): flow rate amplitude at L/2 -3.09e-04 within 0.003
- PASS: finest mesh (robin_m4_n200): flow rate phase at L/2 (rad) 8.23e-05 within 0.003
- PASS: finest mesh (robin_m4_n200): wall radial displacement amplitude at L/2 -1.76e-04 within 0.003
- PASS: finest mesh (robin_m4_n200): wall radial displacement phase at L/2 (rad) -1.15e-03 within 0.005
- PASS: finest mesh (robin_m4_n200): wave speed from the pressure along the tube 1.40e-05 within 0.003
- PASS: finest mesh (robin_m4_n200): attenuation rate, Im(k) -1.96e-03 within 0.01
- PASS: mesh: observed order 2.04 of the profile error at least 1.5
- PASS: mesh: observed order 2.05 of the flow_amp error at least 1.5
- PASS: mesh: observed order 2.01 of the wallMid_amp error at least 1.5
- PASS: mesh: observed order 2.02 of the speed error at least 1.5
- PASS: time step: observed order 2.05 of the flow_amp (finest three steps) at least 1.7
- PASS: time step: observed order 1.91 of the flow_phase (finest three steps) at least 1.7
- PASS: time step: observed order 1.97 of the wallMid_amp (finest three steps) at least 1.7
- PASS: time step: observed order 2.06 of the wallMid_phase (finest three steps) at least 1.7
- PASS: time step: observed order 1.97 of the speed (finest three steps) at least 1.7
- PASS: time step: observed order 2.07 of the attenuation (finest three steps) at least 1.7
  - time error at 200 steps per period: wall amplitude 3.90e-04, wave speed 9.83e-05
- PASS: IQN-ILS (iqnils_m2_n100): velocity profile at L/2, max |error| / max |u| 4.98e-03 within 0.006
- PASS: IQN-ILS (iqnils_m2_n100): wall radial displacement amplitude at L/2 -4.47e-04 within 0.003
- PASS: IQN-ILS (iqnils_m2_n100): wave speed from the pressure along the tube -3.72e-05 within 0.003
- PASS: Robin (robin_m2_n100): velocity profile at L/2, max |error| / max |u| 4.98e-03 within 0.006
- PASS: Robin (robin_m2_n100): wall radial displacement amplitude at L/2 -6.63e-04 within 0.003
- PASS: Robin (robin_m2_n100): wave speed from the pressure along the tube -6.62e-05 within 0.003
- PASS: Robin and IQN-ILS agree on the velocity profile at L/2, max |difference| / max |u| to 5.27e-05 (within 0.001)
- PASS: Robin and IQN-ILS agree on the flow rate amplitude at L/2 to 3.36e-05 (within 0.001)
- PASS: Robin and IQN-ILS agree on the flow rate phase at L/2 (rad) to 1.82e-05 (within 0.001)
- PASS: Robin and IQN-ILS agree on the wall radial displacement amplitude at L/2 to 2.16e-04 (within 0.001)
- PASS: Robin and IQN-ILS agree on the wall radial displacement phase at L/2 (rad) to 5.83e-05 (within 0.001)
- PASS: Robin and IQN-ILS agree on the wave speed from the pressure along the tube to 2.90e-05 (within 0.001)
- PASS: Robin and IQN-ILS agree on the attenuation rate, Im(k) to 1.77e-04 (within 0.001)
  - mean coupling iterations per step: 14.56 (IQN-ILS), 9.25 (Robin); run time 805 s and 381 s

## Observed orders

Unsigned errors (profile): from the errors of the two finest runs. Signed errors: from the differences between the three runs, which cancel the error shared by the series, such as the time error of the mesh study. Time-step series of more than three runs give one order per three successive steps.

| Series | Quantity | Errors | Order |
|---|---|---|---:|
| mesh | profile | 1.79e-02, 4.96e-03, 1.20e-03 | 2.04 |
| mesh | flow_amp | 6.10e-03, 9.39e-04, -3.09e-04 | 2.05 |
| mesh | flow_phase | 6.05e-03, 1.33e-03, 8.23e-05 | 1.92 |
| mesh | wallMid_amp | -1.03e-02, -2.19e-03, -1.76e-04 | 2.01 |
| mesh | wallMid_phase | -2.37e-03, -1.48e-03, -1.15e-03 | 1.40 |
| mesh | speed | -2.33e-03, -4.50e-04, 1.40e-05 | 2.02 |
| mesh | attenuation | -2.08e-03, -2.18e-03, -1.96e-03 | – |
| time step | flow_amp | -1.23e-03, 5.14e-04, 9.39e-04, 1.04e-03 | 2.04, 2.05 |
| time step | flow_phase | 1.43e-04, 1.07e-03, 1.33e-03, 1.40e-03 | 1.86, 1.91 |
| time step | wallMid_amp | 5.28e-03, -6.63e-04, -2.19e-03, -2.58e-03 | 1.96, 1.97 |
| time step | wallMid_phase | -6.12e-03, -2.34e-03, -1.48e-03, -1.28e-03 | 2.14, 2.06 |
| time step | speed | 1.41e-03, -6.62e-05, -4.50e-04, -5.48e-04 | 1.94, 1.97 |
| time step | attenuation | -1.31e-02, -4.24e-03, -2.18e-03, -1.69e-03 | 2.10, 2.07 |

PASSED
