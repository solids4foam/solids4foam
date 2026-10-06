# Manufactured-solution verification studies

These opt-in studies measure the spatial and temporal convergence of the total
Lagrangian solid model with the neo-Hookean manufactured solution. They are
separate from the normal tutorial regression suite because they run many cases.

Source an OpenFOAM.com, OpenFOAM.org, or foam-extend environment, ensure the
tutorial library can be built, and run:

```bash
cd tutorials/solids/hyperelasticity/manufacturedSolution/verification
./Allverify
```

`./Allverify` runs both studies. Use `--study spatial` or `--study temporal`
to run one of them, `--variants` for a comma-separated subset, `--levels` for
custom levels, `--quick` for the first two levels, `--reuse` when resuming a
sweep, and `--keep-going` to continue after a failed case. For example:

```bash
./Allverify --study spatial --variants hex-segregated,hex-petscSnes
./Allverify --study spatial --variants distHex-segregated --levels 5,10,20
./Allverify --study temporal --variants Euler-segregated,backward-segregated
./Allverify --study temporal --levels 10,20,40,80,160,320
```

The cases run under the ignored `verification/work/` directory and the results
are written to the ignored `verification/postProcessing/` directory:
`spatial_convergence.csv`, `temporal_convergence.csv`, and
`verification_summary.md`.

## Spatial study

The case is run in the `steady` mode on 5, 10, 20, and 40 cells per coordinate
direction by default. The displacement and stress errors with respect to the
exact solution are taken at the end of the load ramp from the mean L2 and
L-infinity norms printed by the `neoHookeanManufacturedSolution` function
object. The default variants combine the segregated and PETSc SNES solution
procedures with the regular (`hex`) and distorted (`distHex`) hexahedral meshes.
The structured tetrahedral (`tet`) variants are available explicitly and require
Gmsh. The PETSc SNES variants require a PETSc-enabled solids4foam build.

The net order is measured from the coarsest and finest errors and the effective
cell spacing. Each variant has `minimum_net_order` entries in the reference JSON
for the displacement and the stress, applied to both the L2 and L-infinity
norms; a full sweep passes when all measured net orders meet them and the finest
errors are lower than the coarsest ones. The initial thresholds are 1.5 for the
displacement and 0.5 for the stress. Quick runs check only that the errors are
finite and positive.

## Temporal study

The case is run in the `dynamic` mode on a fixed 10x10x10 hexahedral mesh to the
end time 0.05 s with 10, 20, 40, 80, 160, 320, and 640 time steps by default,
for the first-order Euler, second-order backward (BDF2), and trapezoidal
Newmark-beta second time derivative schemes. Because the spatial error is fixed,
the errors with respect to the exact solution level off once the temporal error
falls below it. The temporal order is therefore measured from the differences
between the final displacement fields of successive time-step sizes: with `e_k =
||D(dt_k) - D(dt_{k+1})||`, the observed order between levels `k` and `k + 1` is
`log(e_k/e_{k+1})/log(dt_k/dt_{k+1})`. A full sweep passes when the order
measured from the finest three levels meets the `minimum_order` entry of the
variant, 0.8 for Euler and 1.5 for the second-order schemes. The errors with
respect to the exact solution are reported alongside. The final displacement
fields are written in ASCII with 12 significant digits for this comparison.

Each temporal variant sets the `timeFunction` in the dynamic parameters file of
the copied case. The default variants use `cosineSquared`,
`T = 0.5*(1 - cos(omegaT*t))^2`, whose displacement, velocity, and
acceleration are all zero at `t = 0`, consistent with the zero initial state of
the time schemes, so that each scheme shows its design order. The
`*-cosine` variants use the `cosine` time function of the dynamic mode,
`T = 1 - cos(omegaT*t)`, whose initial acceleration is nonzero. The three-level
Euler and backward schemes are self-starting and keep their orders, but
`NewmarkBeta` starts from a zero stored acceleration and its error then
converges at first order only; its `minimum_order` is 0.8 in those variants. To
recover second order with `cosine`, the exact initial acceleration would have to
be supplied in the `NewmarkA(D)` field at the start time.

The study is not run by `tutorials/Alltest`.
