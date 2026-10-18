# staticCylinderInChannel verification

Opt-in verification of the four forcing methods of `immersedBoundaryForce`:

- `cutLink`;
- `ghostCell`;
- `penalty`;
- legacy incremental, `method incremental`.

They are verified on the steady 2D-1 benchmark of Schäfer and Turek (1996),
on the four `MESH_LEVEL`s of the tutorial. It is not run by
`tutorials/Alltest` or `tutorials/Alltest-regression`.

| File | Content |
| --- | --- |
| `Allverify` | runs the verification (wrapper of the driver) |
| `Allclean` | removes `work/` and `postProcessing/` |
| `results.md` | meshes, results of all the studies, frozen settings |
| `reference/*.csv` | compact results of every reference case |
| `scripts/static_cylinder_verification.py` | the driver |
| `scripts/plot_verification.py` | summary figures from `reference/*.csv` |
| `scripts/plot_mesh.py` | close-up of a generated mesh |
| `images/` | previews of the summary figures and of the L2 mesh |

The cases are generated in `work/` and the results are written to
`postProcessing/`; neither is tracked, and `./Allclean all` in the tutorial
removes them.

## Summary figures

![Legacy incremental sensitivity](images/legacy_incremental_sensitivity.png)

The legacy incremental drag on L2 and L4 against (a) the outer correctors,
(b) the coupling coefficient and (c) the time step. Five outer correctors
are the historical setting; they do not make the drag independent of the
outer iterations.

![Four-method comparison](images/four_method_comparison.png)

(a) Drag and (b) CPU hours of the four methods on the identical
non-body-fitted L1-L4 meshes. `cutLink`, `ghostCell` and `penalty` use
`nOuterCorrectors 1`, and legacy incremental uses `nOuterCorrectors 5`. The
reference drag is a numerical benchmark value, not an analytical solution.

Hollow markers are limit cycles; the full captions are in `results.md`.
`python3 scripts/plot_verification.py` draws the same data from
`reference/*.csv`, in the Matplotlib style of the script.

## Frozen settings

All the methods use the same meshes, physical settings, time step
(`deltaT 0.001` s, maximum Courant number 0.42 on Level 4), end time (15 s)
and `nCorrectors 3`. They differ only in their method settings:

| Method | Settings |
| --- | --- |
| `cutLink`, `ghostCell`, `penalty` | `nOuterCorrectors 1` |
| legacy incremental | `nOuterCorrectors 5`, `couplingCoeff 0.2` |

The legacy incremental method also uses `occupancy vertexFraction`.
`nOuterCorrectors 5` is its historical setting, and its drag still depends
on the number of outer correctors (see `results.md`).

## Running

Source OpenFOAM and solids4foam, then run the verification locally:

```bash
./Allverify --dry-run        # list the cases and their estimated cost
./Allverify                  # cutLink on levels 1-4, then the check
./Allverify --methods all    # the four methods: the final comparison
./Allverify --help
```

`./Allverify` copies the tutorial for each case into
`work/<study>/<method>-L<level>`, with the frozen settings, and runs it:

- serial cases with the tutorial's `Allrun`;
- parallel cases with `solids4Foam::runParallel`.

It then analyses the cases into `postProcessing/<study>/cases.csv` and
`verification_summary.md`, and prints the result of the check:

- **PASS:** the time-averaged drag of `cutLink` on Level 4, with the frozen
  settings, lies within the Schäfer and Turek range `5.57 <= C_d <= 5.59`.
- **FAIL:** that drag lies outside the range (exit status 1).
- **SKIPPED:** that case was not run or did not complete.

The other methods are compared with their reference runs but not checked
against the range.

The default run (cutLink, about 5 CPU hours on the reference machine) takes
a few hours on a workstation. `--methods all` runs 16 cases, about 72 CPU
hours, of which 44 are for the legacy incremental method on Level 4.

Other useful options:

| Option | Effect |
| --- | --- |
| `--levels 1 2` | run only these mesh levels |
| `--cores N` | MPI ranks per case (default 1, 1, 2 and 8 on L1-L4) |
| `--reuse` | keep completed cases |
| `--force` | regenerate existing cases |

The sensitivity studies override one setting:

| Option | Overrides |
| --- | --- |
| `--n-outer N` | `nOuterCorrectors` |
| `--coupling-coeff C` | `couplingCoeff` of the legacy incremental method |
| `--deltaT DT` | the time step |
| `--end-time T` | the end time |

Overriding a setting adds a tag to the case names (e.g.
`incremental-L2-n3`), so the frozen cases are never replaced. Only the
frozen `cutLink` Level 4 case is checked. The commands of each study are
listed in `results.md`.

Unlike the other verification drivers, this one also has subcommands,
because the full study is long:

```bash
./Allverify generate [options]   # generate the cases only
./Allverify analyse [--study S]  # analyse and check the completed cases
./Allverify submit [options]     # optional: one SLURM job per case
```

`submit` is an optional addition to the usual convention of running
`./Allverify` locally or inside a batch job of your own. It submits each
generated case to SLURM, where the case runs its `Allrun.verify` script.
Without `--yes` it only prints a summary. Its options are:

| Option | Effect |
| --- | --- |
| `--env-script` | script sourced in each job (default `$S4F_ENV_SCRIPT`) |
| `--partition` | SLURM partition |
| `--max-cpus` | refuse to exceed this many queued CPUs |
| `--time`, `--mem` | job wall time and memory (default by level) |
| `--sbatch-args` | further site-specific `sbatch` options |

MPI binding follows the MPI environment, e.g.
`OMPI_MCA_hwloc_base_binding_policy` for Open MPI.

`python3 scripts/plot_mesh.py --preview` draws the mesh of a generated
`cutLink` case around the cylinder (default `work/final/cutLink-L2`), with
the cells it penalises, from its recorded fields. It runs no OpenFOAM
application.

## References

- Schäfer, M., Turek, S. (1996). Benchmark computations of laminar flow
  around a cylinder. In: Hirschel, E.H. (ed.) Flow Simulation with
  High-Performance Computers II. Notes on Numerical Fluid Mechanics, 52,
  547-566.
- John, V., Matthies, G. (2001). Higher-order finite element discretizations
  in a benchmark problem for incompressible flows. International Journal for
  Numerical Methods in Fluids, 37, 885-903.
  <https://doi.org/10.1002/fld.195>
- FeatFlow, DFG benchmark 2D-1 (Re = 20), reference values:
  <https://wwwold.mathematik.tu-dortmund.de/~featflow/en/benchmarks/cfdbenchmarking/flow/dfg_benchmark1_re20.html>
- Nabh, G. (1998). On high order methods for the stationary incompressible
  Navier-Stokes equations. PhD thesis, University of Heidelberg, Preprint
  42/98 (the reference values, as cited by John and Matthies, 2001).
