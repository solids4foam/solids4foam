# radialTorsionalRing verification

This opt-in study runs the tutorial on up to three meshes (M1, M2, M3) and
three drives (torsion, radial, combined) and compares the results with the
semi-analytical reference. Nothing here is run by `tutorials/Alltest` or
`tutorials/Alltest-regression`. It needs `python3` with numpy and scipy.

```bash
./Allverify                         # M1, all three drives
./Allverify torsion M1 M2 M3        # torsion mesh study
./Allverify combined M3 np=8        # M3 runs in parallel (default 4 ranks)
./Allverify analyse                 # analyse the existing cases in work/
./Allclean                          # remove work/
```

Each case is a copy of the tutorial in `work/<drive>_<mesh>`. Its dictionaries
and reference are written by `scripts/makeCase.py`, it is run with the
tutorial `Allrun`, and `scripts/analyse_case.py` then writes `analysis.json`
and `series.csv`.
`scripts/convergence.py` then writes the errors, observed and self-convergence
orders and the mixed terms (combined minus torsion) for all cases in `work/`
to `work/convergence.{txt,json,csv}`.

## Scripts

- `scripts/rt_reference.py`: the reduced finite-amplitude model and its solver.
- `scripts/make_reference.py`: writes `reference.json` at every solver time
  step (about 2 s).
- `scripts/makeCase.py`: writes every dictionary of the case for a mesh, drive
  and rank count. The stored tutorial dictionaries are its M1 combined output:

  ```bash
  python3 scripts/makeCase.py /tmp/chk --mesh M1 --config combined \
      --no-reference
  ```

  then `diff -r` the `0`, `constant` and `system` directories against the
  tutorial; only the last digits of the two drive tables may differ between
  numpy versions.
- `scripts/analyse_case.py`: quantities of interest, errors against the
  reference, and the angular-momentum and energy balances of one case.
- `scripts/convergence.py`: the mesh-convergence and mixed-term analysis.

Errors are $$E = \max_t |q - q_{ref}| / \max_t |q_{ref}|$$ over
$$t > 2\Delta t$$, and observed orders $$\log_2(E_1/E_2)$$ are reported only
when the errors decrease monotonically. The full study (numerical acceptance
driver, thresholds, time-step and strict-coupling checks, and the independent
checks of the reference) is kept with the verification paper.
