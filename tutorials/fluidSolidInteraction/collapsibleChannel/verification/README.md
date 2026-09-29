# collapsibleChannel verification study

This opt-in study compares the wall displacement of the `collapsibleChannel`
tutorial with solutions of the same problem computed with
[oomph-lib](https://oomph-lib.github.io/oomph-lib/), and compares the
second-order (linear) and high-order (cubic) solid discretisations. It is
separate from `regressionTest.sh`, and nothing here is run by
`tutorials/Alltest` or `tutorials/Alltest-regression`.

## Running

Source an OpenFOAM environment with solids4foam built with PETSc, and run:

```bash
cd tutorials/fluidSolidInteraction/collapsibleChannel/verification
./Allverify                    # the static, solid and mesh studies
./Allverify --study solid      # one study: static, solid, mesh or time
./Allverify --quick            # smoke test: the coarsest cases to t = 0.5 s
./Allverify --cores 6          # run six cases at a time
./Allverify --reuse            # resume a sweep without re-running cases
./Allverify --list             # list the cases and exit
```

Each case is a complete copy of the tutorial under the ignored
`verification/work/` directory. Results are written to the ignored
`verification/postProcessing/` directory: a CSV per study, a
`verification_summary.md`, and, when `gnuplot` is available, a plot of each
FSI study's wall-midpoint history against its reference.

Each case runs in serial, because the solid does not converge in parallel
(see the tutorial README); `--cores` sets how many cases run at the same
time. The full sweep takes several hours; the recorded results below were run
one case per core on a Slurm node.

## The reference solutions

The reference beam is computed with oomph-lib's monolithic Newton solver:
Crouzeix-Raviart (Q2P-1) Navier-Stokes elements on an algebraically updated
mesh, geometrically nonlinear Hermite Kirchhoff-Love beam elements, and BDF2
in time. `reference/oomph-lib/collapsible_channel_ref.cc` is the oomph-lib
demo driver `fsi_collapsible_channel.cc` with the parameters, a clamped
option and the pressure ramp exposed on the command line; the header of the
file lists every change. To regenerate the references, build oomph-lib
(its `oomph_build.py`; on macOS pass
`--ext-OOMPH_BUILD_OPENBLAS OFF --ext-OOMPH_USE_OPENBLAS_FROM $(brew --prefix openblas)`),
then:

```bash
cd reference/oomph-lib
cmake -G Ninja -B build -DOOMPH_INSTALL=<oomph-lib>/install
cmake --build build
mkdir -p run
./build/collapsible_channel_ref_CR --out run --refine 1 --dt 0.003125 \
    --sigma0 0 --h 0.05 --q 4e-11 --pext 4e-7 --clamped 1 --tramp 0.25 \
    --tmax 3.5
```

`--refine 1` is the resolution of the oomph-lib tutorial: 100 x 16 fluid
elements and 40 beam elements, 17,356 unknowns. `--steady 20` gives the
static beam deflection used by the `static` study.

`scripts/oomph_trace_to_csv.py` writes the reference CSV files:

- `collapsibleChannel_oomph_reference.csv`, the converged reference, is the
  Richardson extrapolation in time of the runs at
  $$\Delta t = 0.003125$$ and $$0.0015625\,\mathrm{s}$$;
- `collapsibleChannel_oomph_dt40.csv`, `_dt80.csv` and `_dt160.csv` are the
  runs at $$\Delta t = 0.025$$, $$0.0125$$ and $$0.00625\,\mathrm{s}$$, so that
  a solids4foam case can be compared with oomph-lib at its own time step.

Precision of the references, on the wall midpoint up to
$$t = 1.5\,\mathrm{s}$$, against a peak displacement of $$0.216\,\mathrm{m}$$:

| Check | Largest difference (m) | Relative to peak |
| --- | ---: | ---: |
| Refinement 1 vs 2 (200 x 32 elements), $$\Delta t = 0.025$$ | 5.7e-4 | 0.26% |
| Crouzeix-Raviart vs Taylor-Hood elements | 1.4e-7 | < 0.001% |
| Richardson correction of the converged reference | 1.3e-4 | 0.06% |
| $$\Delta t = 0.0031$$ vs $$0.0016$$ | 2.4e-4 | 0.11% |
| Static deflection, refinement 1 vs 3 | 1.8e-4 | 0.13% |

The observed temporal order of oomph-lib is 1.6 to 2.0, as expected of BDF2.
The references are therefore good to about 0.3% of the peak displacement,
dominated by the spatial error at the oomph-lib tutorial resolution.

## Studies and acceptance criteria

All FSI cases run to $$t = 1.5\,\mathrm{s}$$, which covers the collapse, the
first trough, the rebound and the second trough. Errors are the largest
difference in the wall-midpoint displacement history, relative to the peak
displacement of the reference, $$0.216\,\mathrm{m}$$.

- `static`: the wall alone under the full external pressure, raised in 20
  steps, against the static oomph-lib beam deflection
  $$-0.143628\,\mathrm{m}$$. The cubic solid must be within 0.5% on every
  mesh; the linear solid within 1% on its finest mesh, 640 x 8.
- `solid`: the tutorial fluid mesh and time step with solid meshes from
  80 x 8 to 640 x 8, for both solids. On this fluid mesh the comparison with
  oomph-lib carries a fluid discretisation error of about 4%, common to all
  cases, so the solid discretisation error is measured as the difference from
  the finest cubic solid: at most 1% for the cubic solid on every mesh, and
  for the linear solid on 640 x 8. Every cubic case must also be within 5% of
  oomph-lib at the same time step.
- `mesh`: fluid meshes refined 1, 2 and 3 times, with the 160 x 8 cubic solid
  and $$\Delta t = 0.025\,\mathrm{s}$$, against oomph-lib at the same time
  step. The error must decrease and the finest must be within 1%, i.e. three
  times the reference precision.
- `time` (not part of the default sweep, see below): the tutorial meshes at
  $$\Delta t = 0.025$$, $$0.0125$$ and $$0.00625\,\mathrm{s}$$ against the
  converged reference.

The tolerances follow from the results below and the precision of the
references: 1% is about three times the reference precision, and 5% leaves
room for the 4% fluid discretisation error of the tutorial mesh.

## Recorded results

Recorded with OpenFOAM v2412, PETSc 3.24 on Linux (one case per core on a
Slurm node) and PETSc 3.22 on macOS for the static study; `./Allverify`
passes.

### Static wall

| Solid mesh | Linear | Cubic |
| --- | ---: | ---: |
| 80 x 8 (tutorial) | -16.7% | +0.21% |
| 160 x 8 | -6.0% | +0.08% |
| 320 x 8 | -1.1% | +0.06% |
| 640 x 8 | -0.22% | +0.05% |
| 160 x 4 | +1.7% | -0.25% |

The cubic solid converges to 0.05% of the beam, which is the continuum-beam
difference at $$h/L = 0.005$$. The linear solid is too stiff on elongated
cells, converging at about second order in the in-plane spacing.

### Solid discretisation in the coupled problem

| Solid mesh | Solid | vs oomph-lib | vs finest cubic | Trough (m) | It. | s |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| 80 x 8 | linear | 15.6% | 17.5% | -0.18511 | 12.3 | 134 |
| 160 x 8 | linear | 5.8% | 5.6% | -0.20198 | 14.6 | 355 |
| 320 x 8 | linear | 4.1% | 1.1% | -0.21044 | 17.7 | 898 |
| 640 x 8 | linear | 4.2% | 0.26% | -0.21220 | 18.0 | 1793 |
| 80 x 8 | cubic | 4.1% | 0.57% | -0.21303 | 18.9 | 412 |
| 160 x 8 | cubic | 4.3% | 0.16% | -0.21280 | 18.5 | 806 |
| 320 x 8 | cubic | 4.3% | 0.10% | -0.21276 | 19.1 | 1647 |
| 640 x 8 | cubic | 4.4% | - | -0.21275 | 19.2 | 3251 |

"Trough" is the first trough of the wall midpoint, "It." the mean number of
FSI iterations per time step and "s" the run time in seconds. oomph-lib at
the same time step has its first trough at
$$-0.21293\,\mathrm{m}$$. On the tutorial mesh the cubic solid is as accurate
as the linear solid on a mesh eight times finer, at a quarter of its cost;
the high-order solid pays off for this bending wall, unlike the
tension-dominated Mok membrane.

![Solid study](reference/solid_study_wallMid.png)

### Fluid mesh

| Fluid level | Cells | vs oomph-lib | Trough (m) | It. | s |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 1 (tutorial) | 3,200 | 4.3% | -0.21280 | 18.5 | 806 |
| 2 | 12,800 | 1.4% | -0.21315 | 20.1 | 932 |
| 3 | 28,800 | 0.51% | -0.21330 | 19.6 | 1177 |

The error falls by 3.1 and 2.8 times, an observed order of 1.6 and 2.5, to
0.5% on the finest mesh, close to the 0.3% precision of the reference.

![Mesh study](reference/mesh_study_wallMid.png)

## Known limitations

- **Fluid level 4 (51,200 cells) does not run.** The coupling stalls at
  $$t \approx 0.7\,\mathrm{s}$$, during the rebound, with every coupling
  setting tried (secant reuse or none, `relMinSignificant` from 1e-3 to 1e-1,
  regularised QR, looser tolerances, up to 200 iterations); it is not in the
  mesh study.
- **The time-step study does not run to completion.** The added-mass coupling
  of the massless wall grows as $$1/\Delta t^2$$; at
  $$\Delta t = 0.0125$$ and $$0.00625\,\mathrm{s}$$ IQN-ILS stalls or the
  solid Newton solve fails during the collapse. It is kept, outside the
  default sweep, as `./Allverify --study time`. The oomph-lib references
  themselves converge at second order in time: the tutorial time step costs
  them 1.4e-2 m (6.3%) against the converged reference, which is why the
  solid and mesh studies compare at the same time step.
- The solid does not converge in parallel, and the cubic solid does not
  converge with sixteen or more cells through the wall thickness.
