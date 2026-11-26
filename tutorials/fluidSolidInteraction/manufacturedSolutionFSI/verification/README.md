# Manufactured-solution FSI verification study

This opt-in study refines the mesh and the time step of the
`manufacturedSolutionFSI` case together and reports the observed orders of
accuracy. Each level has $$n$$ cells along the interface, $$n \times n/2$$
square cells in the fluid and $$n \times n/4$$ in the solid, and the time step
$$\Delta t = 0.32/n$$. The default levels are $$n$$ = 16, 32 and 64.

Source an OpenFOAM.com environment with a PETSc-enabled solids4foam build and
run:

```bash
cd tutorials/fluidSolidInteraction/manufacturedSolutionFSI/verification
./Allverify
```

The default runs the Aitken coupling and takes about six minutes, most of it
on the 64-cell level. Use `--couplings aitken,robin` to add the Robin-Neumann
coupling, `--levels` for other cell counts (multiples of 4), `--quick` for the
first two levels, and `--reuse` to keep completed cases:

```bash
./Allverify --couplings aitken,robin
./Allverify --quick
./Allverify --levels 16,32,64,128
```

Each level is a copy of the case under the ignored `verification/work/`
directory, with `system/caseParameters` set for the level; the copies share
the tutorial library in `../src`. The errors of the `mmsFsiErrors` function
object at $$t = 0.2$$ and 0.4 are written to
`verification/postProcessing/<coupling>_errors.csv`, and the error tables,
with the orders between successive levels and the net order between the
coarsest and finest levels, to
`verification/postProcessing/verification_summary.md`.

## Acceptance

The criteria are in
`reference/manufacturedSolutionFSI_verification_references.json`. For the L2
errors of `U`, `p`, `D`, `intD`, `intTrac` and `intP` at both times:

- every error is finite and positive;
- on the 16, 32 and 64-cell levels, every error agrees with the reference
  error to within 5%;
- in a full sweep, the net order is at least 0.5.

`--quick` checks the errors but not the orders. The reference errors are
those of the verification paper (OpenFOAM v2512, macOS); the same study on
the `development` branch reproduces them to within 0.5%.

## Expected orders

The partitioned `pimpleFluid` model is first order under this combined
refinement. With the Aitken coupling, the L2 errors and orders are:

| quantity | t | n16 | n32 | n64 | orders |
|---|---|---|---|---|---|
| `U` | 0.2 | 1.09e-2 | 4.56e-3 | 1.90e-3 | 1.26, 1.26 |
| `p` | 0.2 | 3.53e-2 | 1.66e-2 | 8.35e-3 | 1.09, 0.99 |
| `D` | 0.2 | 2.65e-3 | 1.47e-3 | 7.61e-4 | 0.85, 0.95 |
| `intTrac` | 0.2 | 3.33e-2 | 1.58e-2 | 9.28e-3 | 1.08, 0.76 |
| `U` | 0.4 | 6.85e-2 | 2.58e-2 | 1.05e-2 | 1.41, 1.30 |
| `p` | 0.4 | 6.41e-2 | 3.48e-2 | 1.82e-2 | 0.88, 0.94 |
| `D` | 0.4 | 3.72e-3 | 6.96e-4 | 1.59e-4 | 2.42, 2.13 |
| `intTrac` | 0.4 | 7.06e-2 | 2.14e-2 | 8.04e-3 | 1.72, 1.41 |

The Robin-Neumann coupling gives similar orders, with an interface traction
order of about 0.6 at $$t = 0.2$$.

The study is not run by `tutorials/Alltest` or `tutorials/Alltest-regression`.
`verification/Allclean` removes the generated `work/` and `postProcessing/`
directories, as does `./Allclean all` in the case directory.
