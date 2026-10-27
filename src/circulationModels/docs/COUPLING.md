# Swapping 0D for 3D, and the solids4foam merge

## The idea

A 3D domain is, to the rest of the circuit, a thing with terminals. Here the
swappable unit is the left ventricle with its two valves. The lumped version
and the CFD version expose the same terminals `in` (mitral side) and `out`
(aortic side), so the `connections` list never changes:

```c++
LV { type valvedChamber;                                   // 0D
     chamber { V0 5; P0 1; Emin 0.1; Emax 2.5; Tes 0.3; Ted 0.45; initialVolume 500; }
     inletValve { CQ 400; }  outletValve { CQ 350; } }

LV { type coupledDomain;  provider external;               // 3D
     inlet  { patch mitral; valve { CQ 400; reverseFraction 0.005; } }
     outlet { patch aorta;  valve { CQ 350; reverseFraction 0.005; } } }

connections ( (LV SAS SAT SAR SCP SVN RA TV RV PV PAS PAT PAR PCP PVN LA LV) );
```

`tutorials/healthyHumanLVModule` and `tutorials/healthyHumanCoupledReplay` differ
in nothing but that one entry (`diff` them).

## What crosses the interface

Each `coupledDomain` creates two **ports**, `<name>.inlet` and `<name>.outlet`,
one per coupled patch. The valves stay in the 0D network, as in your existing
LV case. At each port:

| CFD gives | network returns |
|---|---|
| the port flow `Q`, positive in the natural direction of the patch (inlet: into the domain, `Qmv`; outlet: out of the domain, `Qav`) | the pressure to impose on that patch |

with the interface law (the bounded bidirectional inverse valve law of
the original left-ventricle coupling, identical to `valveDPfromQ`):

```
outlet   Ppatch = P_node + dP(Q)        inlet   Ppatch = P_node - dP(Q)
dP(Q) =  (Q/CQ)^2          Q >= 0
      = -(|Q|/(fr CQ))^2   Q <  0       fr = reverseFraction, clamped to [1e-3, 1]
```

C++ interface of `circulationNetwork` (all SI), the same shape as
`lumpedCircuitODE::setForcing` plus `lumpedCircuitSolver`:

```c++
circulationNetwork net(dict, t0);
net.setPortForcing("LV.outlet", t0, t1, {c0, c1, c2, c3});   // Q = c0 + c1 a + c2 a^2 + c3 a^3, a=(t-t0)/(t1-t0)
net.setPortForcing("LV.inlet",  t0, t1, {...});
net.advance(t1);
scalar Pav = net.portPressure("LV.outlet");                  // Pa, at the state now
scalar Pmv = net.portPressure("LV.inlet");
```

`applications/testPortCoupling` is a complete working example. It is checked
against the table replay by `verification/Allverify`.

## Why this gives you the coupling algorithms for free

The eight algorithms of the original left-ventricle coupling (`leftVentricle/src/finalAlgorithms`) differ only in
*what forcing polynomial they pass and when they integrate*. They hold the
circuit by composition. Swapping `lumpedCircuitSolver` for `circulationNetwork`
changes the circuit, not the algorithms:

| existing | new |
|---|---|
| `lumpedCircuitODE::setForcing(t0,t1,Qav,Qmv)` | `setPortForcing("LV.outlet" / "LV.inlet", ...)`, same polynomial form, so the third and fourth order schemes still apply |
| `Pav = P1 + dPav(Qav)`, `Pmv = P9 - dPmv(Qmv)` | `portPressure("LV.outlet" / "LV.inlet")` |
| `lumpedCircuitState` checkpoint | `stateDict()` / `readState()` |
| `lumpedParameters` (13 hand-coded states) | `circulationProperties` (any circuit) |

The pipe case is the same with one port and no valve: a bare `externalPort`
`{ direction outlet; provider external; }` attached to a Windkessel.

## Order of work from here

1. **Reproduce the LV circuit exactly.** Write the old 13-state circuit
   (smooth tricuspid/pulmonary valves with `epsilon`, atria, SI parameters from
   `lumpedParameters`) in `circulationProperties` and compare with
   the original `lv0d.py` of the left-ventricle study. This is the gate for retiring
   `lumpedCircuitODE`.
2. **Point the eight algorithms at `circulationNetwork`.** Mechanical: the
   forcing code is already polynomial. Re-run `meshStudy` and
   `orderAccuracyTime`; the numbers should not move.
3. **Pressure-driven ports.** `externalPressurePort` (the network supplies
   Q, the CFD supplies P) for 3D inlets fed by a lumped upstream. Tried and
   removed: on the closed-loop pipe, at either end, the plain partitioned
   iteration diverges (interface gain about 1.5 C L/dt^2 + 0.75 L/L1 > 1,
   because both sides then differentiate: the pipe's pressure follows dQ/dt and
   the compliance's flow follows dP/dt). Needs an interface acceleration, or a
   domain that is not inertia-dominated; revisit with other geometries.
4. **solids4foam.** The library has no OpenFOAM-version dependence beyond
   `ODE` and `OpenFOAM`.
   * Done: the library is in `solids4foam/src/circulationModels`, built by `src/Allwmake`.
   * Add a `physicsModel` type `circulation` (0D only; `physicsModel` needs just a
     `Time`, so the stock `solids4Foam` solver can run a mesh-free case) that wraps
     `circulationNetwork` and implements `evolve()`.
   * For fluid and fluid-solid cases, the pressure boundary condition on each
     patch reads the `portPressure` and posts the patch flux through
     `setPortForcing`; the network lives in the `fluidModel`, not in the BC.
5. **Physiology.** Shi's dynamic leaflet valve (`Shi.pdf`, parameters in
   `ShiParam.jl`) and the passive exponential chamber of ModularCirc as further
   component types, each a single class.
