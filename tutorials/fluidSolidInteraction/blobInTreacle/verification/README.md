# blobInTreacle verification study

This opt-in study checks the `blobInTreacle` tutorial against the published
solution of the elastic half-cylinder case of Liu, Jaiman and Gurugubelli
(2014) and its preprint (Liu, arXiv:1401.0082). It has two parts:

- a **time-step study**, which repeats the tutorial up to `t = 1 s` with
  successively halved time steps and checks the observed temporal order, the
  quantity the case was designed to test;
- a **mesh study**, which refines both meshes, runs each level into the steady
  state, and compares the top-point displacement and the interface shape with
  the published solution.

It is deliberately separate from `regressionTest.sh`: the regression test
checks that the tutorial result does not change, whereas this study checks
convergence towards the published solution. Nothing here is run by
`tutorials/Alltest` or `tutorials/Alltest-regression`.

## Running

Source an OpenFOAM environment with solids4foam built, and run:

```bash
cd tutorials/fluidSolidInteraction/blobInTreacle/verification
./Allverify
```

Useful options:

```bash
./Allverify --quick              # two time steps and two coarse meshes only
./Allverify --study time         # the time-step study only
./Allverify --study mesh         # the mesh study only
./Allverify --levels 0.5,1,2     # a subset of the levels, each double the last
./Allverify --cores 4            # MPI ranks for the levels finer than the tutorial
./Allverify --reuse              # resume a sweep without re-running cases
```

Each run is a complete copy of the tutorial under the ignored
`verification/work/` directory, so the tutorial itself and its regression test
are never modified. Results are written to the ignored
`verification/postProcessing/` directory: `time_convergence.csv`,
`mesh_convergence.csv`, the interface shapes in `interfaces/`,
`verification_summary.md` and, when `gnuplot` is available,
`blobInTreacle_interfaces.png`. The exit code is zero only if every check
passes.

## Reference data

Liu et al. (2014) report only the self-convergence of their errors in time,
not the solution itself. The preprint solves the same case, with the same
geometry, inflow ramp and material, and describes the solid as linear. Two of
its figures are vector graphics, from which the reference data were extracted
exactly rather than digitised:

- `reference/LiuInterface_t1.csv`: the deformed interface at `t = 1 s` from
  Fig. 5.3, computed on a P5/P4/P5 mesh with a first-order scheme and
  `dt = 0.025 s`. The preprint gives the error of this curve against its
  `dt = 5e-5 s` solution as up to `0.0025 m`, about 4% of the interface
  displacement.
- `reference/LiuInterface_t10.csv`: the deformed solid boundary in the steady
  state, `t = 10 s`, from Fig. 5.2, computed by a domain-decomposition solver
  on the coarse mesh of Fig. 5.1 (1 239 fluid and 370 solid triangles, P2/P1/P2
  elements for the combined-field solver). The combined-field solution in the
  same figure agrees with it to `1e-4 m`.
- The steady displacement of the top point, `(0.1297, -0.0117) m`, stored in
  `reference/blobInTreacle_verification_references.json`. It is the
  displacement of the solid mesh vertex at `(1.5043, -0.0006) m`, found by
  matching the triangles of the reference mesh in Fig. 5.1 with those of the
  deformed mesh in Fig. 5.2 one by one.

The coordinates were scaled with the axis ticks of Fig. 5.3 and with the fixed
base of the half cylinder in Figs. 5.1 and 5.2; the extraction error is well
below `1e-3 m`. There are no published histories, forces or tables for this
case, so the comparison is of these shapes and this displacement.

## Time-step study

The tutorial mesh is run to `t = 1 s`, during the inflow ramp, with
`deltaT = 0.25, 0.125, 0.0625` and `0.03125 s`. Both regions use the
second-order `backward` scheme. The order is measured from successive
differences, which needs no reference solution, for the top-point
displacement `u_x` and for the l2 norm of the change in the interface vertex
positions, as in Liu et al. (2014).

The displacement changes very little with the time step: the flow is highly
viscous, and the fluid damps the structure. Below `deltaT = 0.03125 s` the
changes, about `1e-5 m`, are at the level of the coupling and solver
tolerances and no longer measure the time error, which is why the study uses
large steps.

Acceptance criteria:

- the observed orders of `u_x` and of the interface position from the three
  coarsest pairs of steps must lie between `1.6` and `3.0`, for an expected
  order of 2;
- the last change in `u_x` must be below `1e-3` of `u_x`, so the time error of
  the tutorial step, `0.01 s`, is negligible.

## Mesh study

Each level scales the in-plane block divisions of both meshes, leaving the
single spanwise cell alone. Level `x1` is the tutorial mesh and level `x0.5` is
level 1 of the `solid-benchmarks` campaign.

| Level | Fluid cells | Solid cells | `deltaT` (s) |
| ---: | ---: | ---: | ---: |
| x0.5 | 120 | 30 | 0.01 |
| x1 | 480 | 120 | 0.01 |
| x2 | 1 920 | 480 | 0.01 |
| x4 | 7 680 | 1 920 | 0.005 |

Each level runs to `t = 3 s`. The inflow reaches its full value at `t = 2 s`,
and a level only counts as steady if `u_x` changes by less than `1e-3` of its
value over the last `0.5 s`. The finest level needs the smaller time step: at
`deltaT = 0.01 s` its IQN-ILS coupling stagnates.

Acceptance criteria and their basis:

- The change in the steady `u_x` between successive levels must decrease, with
  an observed order of at least `1.5` from the three finest levels, and the
  finest `u_x` must be within 1% of the Richardson extrapolation. The
  discretisation is second order; the observed order is 2.3.
- The finest steady `u_x` must be within 10% and `u_y` within 25% of the
  published top-point displacement. The published steady state comes from a
  coarse mesh, and its `u_y` is small, so these are comparisons
  rather than tight agreements. The converged solids4foam value is 7% below
  the published `u_x`.
- The largest distance of the finest interface from the published one,
  divided by the largest interface displacement, must be below 5% at
  `t = 1 s`, where the reference itself carries about 4% of time error, and
  below 10% in the steady state.

A `--quick` run only exercises the two coarsest meshes and the two largest
time steps, which are not expected to meet the accuracy criteria, so it checks
only that every run completes and reaches a steady state.

## Reference results

Recorded with OpenFOAM v2412 on an Apple M1 Ultra, with `--cores 4`: levels
x0.5 and x1 and the time-step study in serial, levels x2 and x4 on four MPI
ranks. The whole study took about 40 minutes.

Time-step study, tutorial mesh, `t = 1 s`:

| `deltaT` (s) | `u_x` (m) | Change (m) | Order | Interface order | FSI its |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 0.25 | 0.0558539 | – | – | – | 10.8/12 |
| 0.125 | 0.0563958 | 5.42e-4 | – | – | 9.2/11 |
| 0.0625 | 0.0565075 | 1.12e-4 | 2.28 | 2.27 | 10.3/11 |
| 0.03125 | 0.0565286 | 2.10e-5 | 2.41 | 2.26 | 10.8/12 |

Mesh study, top-point displacement in the steady state, and the largest
distance of the interface from the published one divided by the largest
interface displacement:

| Level | `u_x` (m) | `u_y` (m) | Change (m) | Shape, 1 s | Shape, steady |
| ---: | ---: | ---: | ---: | ---: | ---: |
| x0.5 | 0.108186 | -0.012647 | – | 5.5% | 12.6% |
| x1 | 0.117226 | -0.013731 | 9.04e-3 | 1.6% | 7.8% |
| x2 | 0.119711 | -0.014015 | 2.49e-3 | 1.6% | 6.5% |
| x4 | 0.120210 | -0.014027 | 4.99e-4 | 1.6% | 6.3% |

| Level | Cores | Clock time (s) | FSI iterations per step (mean/max) |
| ---: | ---: | ---: | ---: |
| x0.5 | 1 | 45 | 7.7/11 |
| x1 | 1 | 95 | 8.6/13 |
| x2 | 4 | 282 | 9.8/15 |
| x4 | 4 | 1876 | 15.7/53 |

The steady `u_x` converges at an observed order of 2.32 to an extrapolated
`0.12034 m`; the finest level is 0.1% from it and 7.3% below the published
`0.1297 m`, and its `u_y` is 20% larger in magnitude than the published
`-0.0117 m`. The interface at `t = 1 s` agrees with the published one to 1.6%
of the interface displacement on every level from x1, which is within the 4%
time error of the reference curve itself. In the steady state the published
half cylinder leans further downstream, by up to 6% of the displacement once
the mesh is converged.

The steady offset does not come from the time step or the mesh of the
solids4foam solution, both of which are converged well below it. The
published steady state was computed on the coarse mesh of the preprint,
whereas the `t = 1 s` curve, which solids4foam reproduces, was computed with
P5/P4/P5 elements on a refinement of it; the offset is therefore attributed
mainly to the reference, although this cannot be confirmed from the published
data.

Because the published data are limited to self-convergence norms and the
coarse-mesh figures of the preprint, an independent reference solution, e.g.
from COMSOL, LS-DYNA or ANSYS, would be valuable; one is being arranged.

## References

J. Liu, R.K. Jaiman, P.S. Gurugubelli, A stable second-order scheme for
fluid-structure interaction with strong added-mass effects, *Journal of
Computational Physics*, 270, 687-710, 2014.

J. Liu, Combined field formulation and a simple stable explicit interface
advancing scheme for fluid structure interaction, arXiv:1401.0082, 2014.
