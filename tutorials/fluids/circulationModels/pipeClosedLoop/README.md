# pipeClosedLoop: a CFD pipe in a closed circulation loop

A 3D pipe (4 cm x 30 cm, 24,675 cells) whose **both ends** are coupled to a
lumped closed circulation: the outlet drains into the arterial side (C1, R1,
L1, C2, R2, C3, R3), the heart (atrium, mitral valve, ventricle, aortic valve)
pumps, and the blood returns to the pipe inlet (C4).

```text
outlet -> C1 -R1-L1- C2 -R2- C3 -R3- atrium -mitral- ventricle
       -aortic(+R4)- C4 -> inlet
```

The CFD gives the flow through each patch; the circuit returns the pressure of
the node each patch is attached to. The circuit is one dictionary,
`constant/circulationProperties`.

## Run it

```bash
# blockMesh, circuit check, decomposePar, solids4Foam on the
# numberOfSubdomains (8) ranks of system/decomposeParDict; 1 runs serially
./Allrun
# figures open on screen and are saved to plots/   (--no-show: only save)
python3 plotter.py
# the same pipe solved by stock pimpleFoam
./Allrun pimpleFoam
# the 0D reference: the same loop with the pipe replaced by its R-L
# surrogate, solved as one circuit by circulationFoam (no mesh, seconds)
./Allrun 0D
# back to a clean case
./Allclean
```

The pipe is solved by `solids4Foam` with the `pimpleFluid` fluid model, or by
stock `pimpleFoam` with `./Allrun pimpleFoam`; the two agree to the PIMPLE
`residualControl` tolerance of 1e-6, which is also the run-to-run scatter of
one solver on 8 ranks (the Windkessel case, converged to 1e-10, agrees to all
digits written). Only the fluid solver is needed: the coupling
library (`libs` in `system/controlDict`) builds the circuit inside the solver
and integrates its ODEs at every exchange. `circulationFoam` is for
circuit-only cases; here it is used only for `circulationFoam -check` (what
the program understood, and that both patches exist in the mesh). `plotter.py`
needs matplotlib (`pip3 install matplotlib`).

Settings for this demo: one cardiac cycle (`endTime 0.8`), `deltaT 5e-3`,
explicit coupling of order 2. About a minute on 8 free cores (the three-cycle
run about 2.5 minutes); watch it with `tail -f log.solids4Foam`.

## Things to show

* `circulationFoam -check`: the circuit as understood (7 states, 7 nodes, 2
  coupling ports).
* `constant/circulationProperties`: the whole circuit; change `scheme` /
  `order` in `coupling`.
* `0/p`: `circulationPressure` on `inlet` and `outlet`; the only change from an
  ordinary pipe case.
* Results: `postProcessing/circulationCoupling/0/coupling.dat`, one row per
  time step (the flow and imposed pressure at each port, then every element's
  pressure, volume and flow).
* `plots/interface.png`, `heart.png`, `circuit.png`.

## Variations

* **three cycles (the periodic state)**: `endTime 2.4` in `system/controlDict`,
  then `python3 plotter.py` (last cycle)
* **implicit coupling (iterates the exchange within each step)**: `scheme implicit;`
* **higher-order coupling**: `order 3;` (limited to second order by the
  `backward` CFD time scheme)
* **a stiffer ventricle**: `ventricle { ... active { E1 3.0e9; ... } }`

## The 0D reference: `./Allrun 0D`

There is no closed form for this problem. `reference0D/` is the same circuit
with the 3D pipe replaced by its rigid, laminar Poiseuille surrogate, a
resistor `Rpipe` in series with an inductor `Lpipe`, so that the whole loop is
one `circulationNetwork` solved by `circulationFoam` with no CFD and no mesh:

```text
C1 -R1-L1- C2 -R2- C3 -R3- atrium -mitral- ventricle -aortic(+R4)- C4 -pipe- C1
```

* `Rpipe = 8 mu L / (pi r^4)`, with a viscosity `mu` of 0.004 Pa s, a pipe
  length `L` of 0.3 m and a radius `r` of 0.02 m;
* `Lpipe = rho L / A`, with a density `rho` of 1060 kg/m3 and the pipe
  cross-sectional area `A`.

The circuit is `reference0D/constant/circulationProperties`: the circuit of
the coupled case with the two coupling ports replaced by `Rpipe` and `Lpipe`.
`reference0D/system/controlDict` runs three cardiac cycles (`endTime 2.4`)
with `deltaT 1e-4`, which is only the **output interval**: the ODE solver
(`RKDP45`, tight tolerances) takes whatever internal steps it needs, many of
them small because of the hard switching of the valves.

`./Allrun 0D` runs `circulationFoam -check`, which describes the circuit and
checks the set-up, and then solves it. The results are in
`reference0D/postProcessing/circulation/0/circulation.dat`: a row every
`deltaT`, with the pressure, volume and flow of every element.

The solution has no partitioned coupling error, so it is the reference for
the coupled run, which differs from it by the modelling error of the
surrogate: the resolved pipe is not exactly its rigid laminar surrogate. At
the coarse time steps used in the tutorial, the coupling error dominates this
difference. It is also a test of the chamber and valve components and the
network assembly on a stiff, non-smooth, closed-loop circuit.

## Checks and known limits

Explicit order 1 is unstable at `deltaT 1e-3` and stops at `t = 1.997`; use
a smaller `deltaT`. Implicit coupling needs many outer correctors on this
loop: at `deltaT 5e-3` every step reaches the limit of 120, so the run
completes but the coupling is not converged. Use a smaller `deltaT` (about
1e-3 or below) for a converged implicit solution, or the explicit scheme.
With the `backward` CFD time scheme the coupled solution is limited to second
order for coupling orders 3 and 4 (the coupler prints a warning). A third- or
fourth-order BDF CFD time scheme is not A-stable, so it can amplify the
pipe-compliance oscillation of this loop (about 30 Hz): with explicit
coupling of order 2, BDF3 and BDF4 diverged at `deltaT 5e-3`.

The circuit's ODE solver is `RKDP45` with very tight tolerances (relTol 1e-13),
and its step cap is raised (`maxSteps 1e8`) because hard-valve switching needs
many small steps.

Every accepted step is logged in
`postProcessing/circulationCoupling/0/coupling.dat`.
