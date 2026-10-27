# pipeClosedLoop: a CFD pipe in a closed circulation loop

A 3D pipe (4 cm x 30 cm, 24,675 cells, `pimpleFoam`) whose **both ends** are
coupled to a lumped circulation: the outlet drains into the arterial side, the
inlet is fed by the aortic valve, and the heart (atrium, ventricle, two
valves) closes the loop. The CFD delivers the flow through each patch; the
circuit returns the pressure of the node each patch is attached to. Three
cardiac cycles (0.8 s each), `t = 0 .. 2.4 s`.

The circuit is a dictionary, `constant/circulationProperties`:

```text
outlet -> C1 -R1-L1- C2 -R2- C3 -R3- atrium -mitral- ventricle
       -aortic(+R4)- C4 -> inlet
```

```c++
components
{
    // CFD patch
    outlet  { type externalPort; direction outlet; patch outlet; }
    inlet   { type externalPort; direction inlet;  patch inlet;  }
    C1      { type compliance; C 1.5e-10; initialPressure 12720; }
    ...
    atrium  { type mixedChamber; V0 1e-6;
              active { E1 9.8e8; }  passive { P0 22.66; K 4.84e5; }
              activation { type sineSquared; offset 0.544; duration 0.24; }
              initialVolume 1.408e-05; }
    mitral  { type orificeValve; Rhat 5.333e10; }
    // orifice in series with R4
    aortic  { type orificeValve; Rhat 5.333e9; R 1.2e7; }
}
connections
(
    (outlet C1 R1 L1 C2 R2 C3 R3 atrium mitral ventricle aortic C4 inlet)
);
// explicit | implicit, 1..4
coupling    { scheme implicit; order 2; rho 1060; }
```

and `0/p` is `circulationPressure` on both `inlet` and `outlet`. Nothing about
this circuit is in C++: the seven states, their equations, the valve switching
and the two-port exchange come from the library. Two general components were
added for it, `mixedChamber` (activation x active polynomial + exponential
passive filling; also the form used by ModularCirc and Mynard) and a series
resistance on `orificeValve`.

## Running

```bash
# -> runs/implicit/2/5e-3
./makeCase implicit 2 5e-3
# blockMesh, decomposePar, mpirun -np 8 pimpleFoam
./runLadder "implicit" "2" "5e-3" 1 8
./runLadder "explicit implicit" "2 3 4" "5e-3" 1 8
```

A coupled case at `dt = 5e-3` takes about 2.5 minutes on 8 cores. Orders 3 and
4 use the same `backward` CFD time scheme, which limits the coupled solution
to second order (the coupler prints a warning). Every order uses `RKDP45` for
the circuit. The step cap of the circuit's ODE solver is raised
(`maxSteps 1e8`) because hard-valve switching at the very tight tolerances
(relTol 1e-13) needs many small steps.

## Checks and known limits

There is no closed form for this problem. `../pipeClosedLoopMonolithic` is the
same circuit with the pipe replaced by its rigid Poiseuille R-L surrogate,
solved as one system by `circulationFoam`. It is a reference for the coupled
runs, which differ from it by the modelling error of the surrogate, and a
test of the chamber and valve components and the network assembly on a stiff,
non-smooth, closed-loop circuit. Comparisons with an independent RK4
integrator and with earlier simulations of the same problem are not part of
this directory yet.

Some combinations are unstable at the coarse time steps and stop: explicit
order 1 at `dt = 1e-3`, and implicit orders 3 and 4 at `dt = 5e-3`. They need a
smaller `dt`.

Every accepted step is logged in
`postProcessing/circulationCoupling/0/coupling.dat`.
