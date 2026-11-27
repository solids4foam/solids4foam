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

The default runs the Aitken coupling and takes about seven minutes, most of it
on the 64-cell level. Use `--couplings aitken,robin` to add the Robin-Neumann
coupling, `--levels` for other cell counts (multiples of 8), `--quick` for the
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
- in a full sweep, the net order is at least 1.3.

`--quick` checks the errors but not the orders. There are two sets of
reference errors (OpenFOAM v2512, macOS):

- `development` (the default): the case as committed, on the `development`
  branch;
- `withFixes`: with solids4foam PRs
  [#756](https://github.com/solids4foam/solids4foam/pull/756), which makes
  `tractionPatchSurfaceGradient` take effect, and
  [#776](https://github.com/solids4foam/solids4foam/pull/776), with
  `fluxConsistentPatches (outlet);` uncommented in
  `system/fluid/fvSolution`. Select it with `--references withFixes`.

## Expected orders

The L2 errors and orders on `development`:

| coupling | quantity | t | n16 | n32 | n64 | orders |
|---|---|---|---|---|---|---|
| Aitken | `U` | 0.2 | 1.01e-02 | 3.49e-03 | 1.17e-03 | 1.53, 1.58 |
| Aitken | `p` | 0.2 | 2.29e-02 | 6.97e-03 | 2.18e-03 | 1.72, 1.67 |
| Aitken | `D` | 0.2 | 1.28e-03 | 3.01e-04 | 6.60e-05 | 2.09, 2.19 |
| Aitken | `intTrac` | 0.2 | 2.92e-02 | 8.67e-03 | 2.59e-03 | 1.75, 1.75 |
| Aitken | `U` | 0.4 | 3.84e-02 | 1.29e-02 | 4.19e-03 | 1.57, 1.63 |
| Aitken | `p` | 0.4 | 6.74e-02 | 2.17e-02 | 7.55e-03 | 1.64, 1.52 |
| Aitken | `D` | 0.4 | 3.72e-03 | 7.56e-04 | 1.77e-04 | 2.30, 2.10 |
| Aitken | `intTrac` | 0.4 | 8.89e-02 | 2.88e-02 | 8.89e-03 | 1.62, 1.70 |
| Robin | `U` | 0.2 | 2.49e-02 | 9.36e-03 | 2.32e-03 | 1.41, 2.01 |
| Robin | `p` | 0.2 | 3.91e-02 | 1.45e-02 | 4.54e-03 | 1.43, 1.68 |
| Robin | `D` | 0.2 | 1.38e-03 | 3.50e-04 | 8.17e-05 | 1.98, 2.10 |
| Robin | `intTrac` | 0.2 | 2.61e-02 | 8.57e-03 | 2.49e-03 | 1.61, 1.78 |
| Robin | `U` | 0.4 | 4.40e-02 | 1.72e-02 | 5.52e-03 | 1.35, 1.64 |
| Robin | `p` | 0.4 | 8.66e-02 | 3.01e-02 | 1.02e-02 | 1.52, 1.56 |
| Robin | `D` | 0.4 | 3.81e-03 | 8.05e-04 | 2.05e-04 | 2.24, 1.97 |
| Robin | `intTrac` | 0.4 | 1.14e-01 | 3.42e-02 | 9.09e-03 | 1.74, 1.91 |

With #756 and #776, the orders between n = 16, 32, 64 and 128 at
$$t = 0.4$$, and the n128 errors, are:

| coupling | quantity | n128 | orders |
|---|---|---|---|
| Aitken | `U` | 7.24e-04 | 1.71, 1.83, 1.96 |
| Aitken | `p` | 1.12e-03 | 1.82, 2.05, 1.89 |
| Aitken | `D` | 3.83e-05 | 2.29, 2.20, 2.11 |
| Aitken | `intD` | 8.09e-05 | 1.90, 2.04, 2.02 |
| Aitken | `intTrac` | 1.34e-03 | 1.73, 2.14, 1.87 |
| Aitken | `intP` | 2.44e-03 | 1.69, 1.60, 1.28 |
| Robin | `U` | 1.09e-03 | 1.37, 1.73, 2.07 |
| Robin | `p` | 2.45e-03 | 1.61, 1.77, 1.66 |
| Robin | `D` | 3.91e-05 | 2.29, 2.19, 2.11 |
| Robin | `intD` | 7.09e-05 | 1.93, 2.15, 2.10 |
| Robin | `intTrac` | 1.53e-03 | 1.79, 2.30, 1.91 |
| Robin | `intP` | 2.32e-03 | 1.87, 2.07, 1.25 |

On a build with both changes, this case reproduces the n16 errors of these
tables to ten digits. Without #756 the solid displacement falls below second
order on finer meshes (order 1.6 between n = 64 and 128 with Aitken), and
without #776 a first-order pressure
error layer in the outlet cells limits the fluid pressure. The fluid
interface pressure `intP`, and the Robin `p`, stay below second order.

The study is not run by `tutorials/Alltest` or `tutorials/Alltest-regression`.
`verification/Allclean` removes the generated `work/` and `postProcessing/`
directories, as does `./Allclean all` in the case directory.
