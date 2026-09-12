# cavityFlexibleBottom verification study

This opt-in study migrates the mesh study from
`solid-benchmarks/papers/JFNK_quasi_monolithic_FSI/fluidSolidInteraction/cavityFlexibleBottom`
into the tutorial itself. It uniformly refines the tutorial fluid and solid
meshes, runs each level to a steady FSI response, and compares the steady
vertical displacement at `(4 -1 0.5)` and the steady vertical interface force
with the mesh study of Tuković et al. (2018).

It is deliberately separate from `regressionTest.sh`: the regression test
checks that the tutorial remains numerically stable, whereas this study checks
convergence towards published steady values. Nothing here is run by
`tutorials/Alltest` or `tutorials/Alltest-regression`.

## Running

```bash
source ~/bin/load-openfoam v2512
cd tutorials/fluidSolidInteraction/cavityFlexibleBottom/verification
./Allverify
```

Useful options:

```bash
./Allverify --quick              # levels 1 and 2 only, smoke test
./Allverify --variants iqnils    # IQN-ILS instead of Aitken coupling
./Allverify --end-time 400       # longer pseudo-transient
./Allverify --delta-t 0.1        # smaller time step, needed by level 4
./Allverify --reuse              # resume a sweep without re-running cases
```

Each level is a complete copy of the tutorial under the ignored
`verification/work/` directory, so the tutorial itself and its regression test
are never modified. Results are written to the ignored
`verification/postProcessing/` directory as `mesh_convergence.csv`,
`verification_summary.md`, and, when `gnuplot` is available,
`cavityFlexibleBottom_meshConvergence.pdf`.

The study is serial, like the campaign it is migrated from.

## Mesh levels

Each level doubles the in-plane block divisions of both the fluid and the solid
mesh, leaving the single spanwise cell alone. The resulting spacings match the
four published mesh-study points:

| Level | ΔX (m) | Fluid cells | Solid cells |
|---:|---:|---:|---:|
| 1 | 0.1 | 1 688 | 160 |
| 2 | 0.05 | 6 752 | 640 |
| 3 | 0.025 | 27 008 | 2 560 |
| 4 | 0.0125 | 108 032 | 10 240 |

Level 1 is the mesh shipped with the tutorial. Levels 1 to 3 form the default
sweep.

Level 4 is reachable with `--levels 1,2,3,4`, but **not at the tutorial time
step**: partitioned coupling on that mesh diverges at `deltaT = 0.5` with both
the Aitken and the IQN-ILS variant, because the added-mass effect grows as the
interface cells shrink. Run it with a smaller `--delta-t`, and expect a
correspondingly long serial run. The driver reports a diverged level as a
failure rather than reading its blown-up final value as a steady result.

## Steady state

The tutorial is a pseudo-transient route to a steady response, and the
published quantities are steady values. The driver runs each level to a fixed
end time, then checks that the monitored displacement and force have actually
plateaued: the spread of each quantity over the closing tenth of the run must
be within `1e-4` of its final value. A level that has not settled fails rather
than being compared with the published steady value.

The default end time is the tutorial's own `200`. Use `--end-time` when a finer
mesh needs longer.

## Acceptance criteria

- Every level must reach a steady response, as defined above.
- The change in the steady displacement between successive levels must decrease
  monotonically with a positive net order.
- On the finest mesh, the steady displacement must be within 10% and the steady
  interface force within 5% of the published value at the same spacing.

The displacement tolerance is looser than the force tolerance because the two
quantities converge very differently here: the interface force agrees with the
published value to a fraction of a percent on every mesh, whereas the
displacement carries a visible discretisation error on the coarser meshes and
approaches the published curve as the mesh is refined.

The published force values are quoted for an out-of-plane thickness of
`0.05 m`; the tutorial is one metre thick, so the driver compares against
twenty times the published values.

A `--quick` run only exercises the two coarsest meshes, which are not expected
to meet the accuracy tolerances, so it checks only that both levels reach a
steady response.

## Reference results

Recorded with OpenFOAM v2512 and the default `aitken` variant on an Apple M1
Ultra. The whole default sweep took about ten minutes in serial.

| Level | ΔX (m) | Steady u<sub>y</sub> (m) | vs. reference | Steady F<sub>y</sub> (N) | vs. reference | u<sub>y</sub> change (m) |
|---:|---:|---:|---:|---:|---:|---:|
| 1 | 0.1 | -0.20246 | 13.2% | -5.15312 | 0.13% | – |
| 2 | 0.05 | -0.22883 | 6.88% | -5.23894 | 0.19% | 2.64e-02 |
| 3 | 0.025 | -0.23563 | 5.59% | -5.25632 | 0.14% | 6.79e-03 |

The interface force agrees with the published mesh study to within 0.2% on
every level. The displacement converges cleanly at a net order of 1.96, but
towards a value about 5% away from the published curve rather than onto it.
That residual offset is consistent with the difference already noted in the
tutorial README between the solids4foam prediction and the published
displacement at the coarsest spacing, and it is why the displacement tolerance
is looser than the force tolerance.

## References

Ž. Tuković, A. Karač, P. Cardiff, H. Jasak, A. Ivanković, OpenFOAM finite
volume solver for fluid-solid interaction, *Transactions of FAMENA*, 42(3),
1–31, 2018.
