# Failed 0.3 m/s production attempt (jobs 5323918-5323921, 2026-10-10)

Per-iteration coupling residuals (`max(r1, r2)`, strict test, tolerance 1e-6,
PIMPLE `nOuterCorrectors 1`) of the four aborted MeluXina runs (original
`N_t = 4f` family, dt = 0.00625 s). Each run reached `nOuterCorr 200` with r1
stuck at 1.6e-6 (L0, t = 1.78 s), 3.7e-6 (L1, 5.09 s), 2.8e-5 (L2, 3.47 s)
and 2.6e-6 (L3, 1.10 s) while r2 was ~2e-11. Full logs are archived on
MeluXina in `/project/home/p201287/philip/beam-richter-opus/archive/`.
