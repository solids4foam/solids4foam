# channelLeaflet code-to-code comparison with oomph-lib

This opt-in study compares the leaflet-tip displacement of the
`channelLeaflet` tutorial with solutions of the same problem computed with
[oomph-lib](https://oomph-lib.github.io/oomph-lib/). oomph-lib models the
leaflet as a Kirchhoff-Love beam loaded on its midplane, and solids4foam as a
plane-strain continuum of thickness $$h = 0.05\,\mathrm{m}$$, so this is a
code-to-code comparison with a model difference, not a verification against a
converged solution of the same model. It is separate from `regressionTest.sh`,
and nothing here is run by `tutorials/Alltest` or
`tutorials/Alltest-regression`.

## Running

Source an OpenFOAM.com environment with solids4foam built with PETSc, and
run:

```bash
cd tutorials/fluidSolidInteraction/channelLeaflet/verification
./Allverify                    # the default studies: static, comparison, time
./Allverify --study time       # one study
./Allverify --study mesh       # an optional, exploratory study (see below)
./Allverify --quick            # smoke test: the tutorial case to t = 0.5 s
./Allverify --cores 4          # run four cases at a time
./Allverify --reuse            # resume without re-running completed cases
./Allverify --list             # list the cases and exit
```

Each case is a copy of the tutorial under the ignored `verification/work/`
directory; results are written to the ignored `verification/postProcessing/`
directory: a CSV per study, a `verification_summary.md` and, when `gnuplot`
is available, a plot of each FSI study against its reference. `--reuse`
accepts a case only if a fingerprint of the tutorial inputs, the case
settings, the set-up code and the solver build matches, and the run finished
(`End` in the log) with a complete, finite tip history. `./Allverify` exits
with 0 only when every check passes.

Each case runs in serial; `--cores` sets how many run at the same time. The
default studies take about 5 hours with one case per core, the time step
0.005 s being the longest case.

## The oomph-lib reference

`reference/oomph-lib/channel_with_leaflet_ref.cc` is the oomph-lib demo
driver `fsi_channel_with_leaflet.cc` with the physical and numerical
parameters on the command line, a start from rest with a ramped flux, and
leaflet inertia (`--lambdasq`, with a Newmark timestepper); the header of the
file lists every change. It solves the fully coupled problem with Newton's
method, Taylor-Hood elements on a Z2-adaptive algebraic mesh, BDF2 for the
fluid and Hermite beam elements. To regenerate the references, build oomph-lib
(its `oomph_build.py`; on macOS pass `--ext-OOMPH_BUILD_OPENBLAS OFF`
`--ext-OOMPH_USE_OPENBLAS_FROM $(brew --prefix openblas)`), then:

```bash
cd reference/oomph-lib
cmake -G Ninja -B build -DOOMPH_INSTALL=<oomph-lib>/install
cmake --build build
mkdir -p run
./build/channel_with_leaflet_ref --out run --lright 7 --nre 0 --tramp 1 \
    --lambdasq 2e-4 --dt 0.01 --tmax 10
```

`scripts/oomph_trace_to_csv.py --tmax 10` writes the reference CSV files:

- `channelLeaflet_oomph_reference.csv`: the Richardson extrapolation in time
  of the runs at $$\Delta t = 0.02$$ and $$0.01\,\mathrm{s}$$;
- `channelLeaflet_oomph_dt100.csv`: the run at $$\Delta t = 0.01\,\mathrm{s}$$,
  the tutorial time step;
- `channelLeaflet_oomph_h<h>.csv`: leaflets of thickness 0.05, 0.025 and
  0.0125 m at the same bending stiffness and mass per unit length
  ($$Q \propto h^3$$, $$\Lambda^2 \propto h^2$$), Richardson-extrapolated from
  $$\Delta t = 0.05$$ and $$0.025\,\mathrm{s}$$.

Precision of the reference, against the largest tip displacement of 0.43 m:

| Check | Largest difference |
| --- | ---: |
| Z2 tolerance 3e-4 (default 1e-3), 20 beam elements (default 5) | 0.2% |
| Downstream length 11 instead of 7, $$t \le 2.4\,\mathrm{s}$$ | 0.01% |
| Richardson correction of the converged reference | 0.05% |
| $$\Delta t = 0.01$$ against the converged reference, periodic state | 0.01% |
| Change of the periodic state from the fourth to the fifth period | 0.02% |
| Steady-ramp start against start from rest, from the fourth period | 0.04% |

The reference is good to about 0.2%, dominated by its spatial resolution.
The regularising density, $$\rho_s = \rho_f$$, moves the oomph-lib periodic
state by 0.4% from the massless leaflet, and 2.4% during the start-up; the
references are computed with it, so the comparison is like for like.

## Studies and acceptance criteria

Errors are the largest difference of the tip displacement vector from the
reference, relative to the largest tip displacement of the reference, 0.43 m.
All FSI cases run to $$t = 8\,\mathrm{s}$$, four periods.

Default studies:

- `static`: the leaflet alone, clamped, under a uniform traction of 0.01 Pa
  on its upstream face, against the Euler-Bernoulli cantilever deflection
  $$q L^4/(8 E_{eff} h^3/12)$$, for $$h = 0.05$$ m on solid meshes 4 x 40 to
  16 x 160 and for $$h = 0.025$$ and $$0.0125$$ m (4 x 40) at the same
  bending stiffness. The cubic solid must be within 1.5% on every mesh.
- `comparison`: the tutorial (cubic solid) and its linear variant against
  oomph-lib at the same time step. Over the last period, $$6 < t < 8$$, the
  cubic solid must be within 2% and the linear solid within 5%, and the last
  period must repeat the one before to 1%.
- `time`: the tutorial at $$\Delta t = 0.01$$ and $$0.005\,\mathrm{s}$$; over
  the periodic state, $$t \ge 4\,\mathrm{s}$$, the two must agree to 0.5%.

Optional studies, run only when named with `--study`; they have no acceptance
criteria beyond completing, and do not currently complete:

- `mesh`: the fluid mesh refined 1.5 and 2 times;
- `thickness`: leaflets of thickness 0.05, 0.025 and 0.0125 m at the same
  bending stiffness and mass per unit length, against oomph-lib at the same
  thickness.

## Recorded results

Recorded with OpenFOAM v2412 and PETSc 3.24 on Linux, one case per core on a
Slurm node; `./Allverify` passes.

### Static leaflet

| Solid mesh | $$h$$ (m) | Linear | Cubic |
| --- | ---: | ---: | ---: |
| 4 x 40 (tutorial) | 0.05 | -4.7% | +0.42% |
| 8 x 80 | 0.05 | -2.1% | +0.84% |
| 16 x 160 | 0.05 | -0.03% | +0.87% |
| 4 x 40 | 0.025 | | -0.40% |
| 4 x 40 | 0.0125 | | -1.2% |

The cubic solid converges to 0.9% above the beam at $$h = 0.05$$ m, the
shear and finite-thickness correction of a cantilever with $$h/L = 0.1$$.

### Comparison with oomph-lib

Tutorial mesh (16,160 fluid cells, 4 x 40 solid cells),
$$\Delta t = 0.01\,\mathrm{s}$$, against oomph-lib at the same time step:

| Solid | Last period | Whole run | Periodicity | FSI it. | Run time (h) |
| --- | ---: | ---: | ---: | ---: | ---: |
| cubic | 1.1% | 5.5% | 0.45% | 16 | 2.3 |
| linear | 3.6% | 3.6% | 0.51% | 15 | 2.1 |

"Whole run" includes the start-up transient, where the difference is
largest, during the first large deflection. The periodic-state difference of
the cubic solid, 1.1-1.3%, is five times the reference precision. It is the
model difference of this discretisation of the finite-thickness leaflet, not
a discretisation error that the study shows to converge: see below.

### Time step

Halving the time step to $$0.005\,\mathrm{s}$$ changes the periodic state by
0.3% (0.9% during the start-up). $$\Delta t = 0.02\,\mathrm{s}$$ and
$$0.0025\,\mathrm{s}$$ fail in the solid Newton solve during the first large
deflection, with every line search, trust-region Newton, a tighter Krylov
tolerance and without the solid predictor, so no order of convergence is
claimed.

## What does not work

These results were found while preparing the study and are the reason the
comparison is code to code.

- **The leaflet-thickness study does not converge to the beam model.** At the
  same bending stiffness and mass per unit length, oomph-lib's solution
  changes by 0.4% at $$h = 0.025$$ m and 0.6% at $$0.0125$$ m (only the
  extensional stiffness of the beam changes). The solids4foam solution changes
  by 3% when $$h$$ is halved, and its difference from oomph-lib grows, to
  2.5-2.8% at $$h = 0.025$$ m and about 3% at $$0.0125$$ m (which failed at
  $$t = 2.7\,\mathrm{s}$$). The difference is therefore not $$O(h)$$ on the
  tutorial mesh; it probably includes the fluid discretisation error near the
  thin leaflet, whose cells scale with $$h$$, but the refined meshes that
  would separate the two do not run.
- **Refined fluid meshes stall in the coupling.** With the fluid mesh refined
  1.5 and 2 times, the Aitken iterations stall at $$t \approx 0.65\,\mathrm{s}$$
  during the first large deflection (final residual 1.6e-4 and 2.6e-3 after
  200 iterations, tolerance 1e-4).
- **Only Aitken completes.** IQN-ILS diverges between $$t = 0.36$$ and
  $$0.54\,\mathrm{s}$$ with every setting tried (initial relaxation 0.005 and
  0.05, reuse of 0, 2 and 5 time steps, with and without filtering, column
  normalisation, residual-sum preconditioning, a tighter tolerance): the
  solid Newton solve or the fluid pressure solve fails after an interface
  update. Robin-Neumann does not converge in the first time step from rest
  with the secant, thickness-limited or constant Robin coefficient, but
  restarted from the Aitken solution at $$t = 1\,\mathrm{s}$$ it converges and
  agrees with Aitken to 0.1%.
- **Lighter leaflets fail.** $$\rho_s/\rho_f = 0.3$$ and $$0.1$$ fail in the
  solid Newton solve in the first period; $$\rho_s/\rho_f = 0.1$$ would move
  the oomph-lib periodic state by only 0.04% from the massless leaflet,
  against 0.4% for $$\rho_s/\rho_f = 1$$.
