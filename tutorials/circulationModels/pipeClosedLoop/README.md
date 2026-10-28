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
# blockMesh, circuit check, decomposePar, pimpleFoam on 8 ranks
# (./Allrun 4 for 4)
./Allrun
# figures open on screen and are saved to plots/   (--no-show: only save)
python3 plotter.py
# back to a clean case
./Allclean
```

Only `pimpleFoam` is needed: the coupling library (`libs` in
`system/controlDict`) builds the circuit inside the solver and integrates its
ODEs at every exchange. `circulationFoam` is for circuit-only cases; here it is
used only for `circulationFoam -check` (what the program understood, and that
both patches exist in the mesh). `plotter.py` needs matplotlib
(`pip3 install matplotlib`).

Settings for this demo: one cardiac cycle (`endTime 0.8`), `deltaT 5e-3`,
explicit coupling of order 2. About a minute on 8 free cores (the three-cycle
run about 2.5 minutes); watch it with `tail -f log.pimpleFoam`.

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

## Checks and known limits

There is no closed form for this problem. `../pipeClosedLoopMonolithic` is the
same circuit with the pipe replaced by its rigid Poiseuille R-L surrogate,
solved as one system by `circulationFoam`. It is a reference for the coupled
run, which differs from it by the modelling error of the surrogate, and a
test of the chamber and valve components and the network assembly on a stiff,
non-smooth, closed-loop circuit.

Explicit order 1 is unstable at `deltaT 1e-3` and stops at `t = 1.997`; use
a smaller `deltaT`. Implicit coupling of orders 2, 3 and 4 runs through three
cycles at `deltaT 5e-3` with the `backward` CFD time scheme. That scheme limits
the coupled solution to second order (the coupler prints a warning for orders 3
and 4). A third- or fourth-order BDF CFD time scheme is not A-stable, so it can
add energy to the pipe-compliance oscillation of this loop (about 30 Hz) and
needs a smaller `deltaT`: BDF4 at `deltaT 5e-3` blew up, and was stable at
`3e-3`.

The circuit's ODE solver is `RKDP45` with very tight tolerances (relTol 1e-13),
and its step cap is raised (`maxSteps 1e8`) because hard-valve switching needs
many small steps.

Every accepted step is logged in
`postProcessing/circulationCoupling/0/coupling.dat`.
