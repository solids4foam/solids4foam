# Swapping 0D for 3D

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

The CFD coupling tutorials in `tutorials/circulationModels` (`pipeWindkessel`,
`pipeClosedLoop`) use the simplest case of this, an `externalPort` attached to a
circuit, bound to a patch by name.

## What crosses the interface

Each `coupledDomain` creates two **ports**, `<name>.inlet` and `<name>.outlet`,
one per coupled patch. The valves stay in the 0D network. At each port:

| CFD gives | network returns |
|---|---|
| the port flow `Q`, positive in the natural direction of the patch (inlet: into the domain, `Qmv`; outlet: out of the domain, `Qav`) | the pressure to impose on that patch |

with the interface law (a bounded bidirectional inverse valve law):

```
outlet   Ppatch = P_node + dP(Q)        inlet   Ppatch = P_node - dP(Q)
dP(Q) =  (Q/CQ)^2          Q >= 0
      = -(|Q|/(fr CQ))^2   Q <  0       fr = reverseFraction, clamped to [1e-3, 1]
```

C++ interface of `circulationNetwork` (all SI):

```c++
circulationNetwork net(dict, t0);
net.setPortForcing("LV.outlet", t0, t1, {c0, c1, c2, c3});   // Q = c0 + c1 a + c2 a^2 + c3 a^3, a=(t-t0)/(t1-t0)
net.setPortForcing("LV.inlet",  t0, t1, {...});
net.advance(t1);
scalar Pav = net.portPressure("LV.outlet");                  // Pa, at the state now
scalar Pmv = net.portPressure("LV.inlet");
```

The flow over a step is passed as a polynomial of up to third order, so that
the coupling schemes of up to fourth order can be built on it. The state of the
network is saved and restored with `stateDict()` / `readState()`.

The pipe case is the same with one port and no valve: a bare `externalPort`
`{ direction outlet; provider external; }` attached to a Windkessel. The coupling
algorithms (explicit or implicit, orders 1 to 4) differ only in *what forcing
polynomial they pass and when they integrate*; they are implemented in
`circulationCoupling` and chosen in the `coupling` dictionary.

## Possible future work

1. **Pressure-driven ports.** `externalPressurePort` (the network supplies
   Q, the CFD supplies P) for 3D inlets fed by a lumped upstream. Tried and
   removed: on the closed-loop pipe, at either end, the plain partitioned
   iteration diverges (interface gain about 1.5 C L/dt^2 + 0.75 L/L1 > 1,
   because both sides then differentiate: the pipe's pressure follows dQ/dt and
   the compliance's flow follows dP/dt). Needs an interface acceleration, or a
   domain that is not inertia-dominated; revisit with other geometries.
2. **A `physicsModel` for solids4foam.** Add a `physicsModel` type `circulation`
   (0D only; `physicsModel` needs just a `Time`, so the stock `solids4Foam`
   solver can run a mesh-free case) that wraps `circulationNetwork` and
   implements `evolve()`. For fluid and fluid-solid cases, the pressure
   boundary condition on each patch would read the `portPressure` and post the
   patch flux through `setPortForcing`; the network would live in the
   `fluidModel`, not in the BC.
3. **Physiology.** Dynamic leaflet valves (Shi) and the passive exponential
   chamber of ModularCirc as further component types, each a single class.
4. **1D vessel models** alongside the 0D ones, in `src/circulationModels`.
