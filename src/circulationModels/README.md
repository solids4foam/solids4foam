# circulationModels

Models of the circulatory system for solids4foam, from a mesh-free network of
lumped (0D) components, such as chambers, valves, resistances and compliances,
to the coupling of that network with a 3D CFD domain. It is currently only
available for OpenFOAM.com (tested with v2412), where the libraries are built
with solids4foam. 1D vessel models are planned and will live here too.

## Contents

- `circulation0D` (library `libcirculation0D`): Lumped-parameter (0D)
  circulation network: components, waveforms, per-cycle metrics, presets, and
  the ODE solution of the network. See
  [circulation0D/README.md](circulation0D/README.md)
- `circulationCoupling` (library `libcirculationCoupling`): Coupling of a
  circulation network to a CFD domain: the `circulationPressure` boundary
  condition and the coupler (explicit/implicit, orders 1 to 4)
- `tools` (no library): `circulationPlot` (plots result columns, needs
  gnuplot)

The `circulationCoupling` library is not specific to the 0D models: the
coupler talks to the network through ports. The coupling is described in
[Swapping 0D for 3D](#swapping-0d-for-3d) below.

## Time scheme

Coupling orders 3 and 4 want a higher-order time scheme in the CFD, so that the
CFD time error does not dominate. OpenFOAM and solids4foam have no
first-derivative `ddtScheme` above second order (the solids4foam `bdfD2dt2`
scheme is a second-derivative scheme, `d2dt2Schemes`, for solid dynamics), so
the tutorials use `backward`, which limits the coupled solution to second
order. The coupler prints a warning when the CFD time scheme is of lower order
than the coupling. A `BDF` ddt scheme of order 1 to 6 is being proposed
separately.

The application is `applications/solvers/circulationFoam`, which solves a
circulation network with no mesh.

The tutorials are in `tutorials/fluids/circulationModels`.

## Building

The libraries are built by `src/Allwmake` and the applications by
`applications/Allwmake`, along with the rest of solids4foam. To build only this
module, after sourcing OpenFOAM:

```bash
(cd src/circulationModels && ./Allwmake)
(cd applications/solvers/circulationFoam && ./Allwmake)
```

The libraries are written to `$FOAM_MODULE_LIBBIN`, like the other solids4foam
libraries. For a CFD case, load the coupling library in `system/controlDict`:

```c++
libs ("libcirculationCoupling.so");
```

## Quick start

```bash
cd tutorials/fluids/circulationModels/healthyHuman
./Allrun    # 20 cardiac cycles of a 4-chamber closed loop, about 20 s
# results: postProcessing/circulation/0/circulation.dat
```

## Swapping 0D for 3D

### The idea

A 3D domain is, to the rest of the circuit, a thing with terminals. Here the
swappable unit is the left ventricle with its two valves. The lumped version
and the CFD version expose the same terminals `in` (mitral side) and `out`
(aortic side), so the `connections` list never changes:

```c++
LV { type valvedChamber;                                   // 0D
     chamber { V0 5; P0 1; Emin 0.1; Emax 2.5; Tes 0.3; Ted 0.45;
               initialVolume 500; }
     inletValve { CQ 400; }  outletValve { CQ 350; } }

LV { type coupledDomain;  provider external;               // 3D
     inlet  { patch mitral; valve { CQ 400; reverseFraction 0.005; } }
     outlet { patch aorta;  valve { CQ 350; reverseFraction 0.005; } } }

connections ( (LV SAS SAT SAR SCP SVN RA TV RV PV PAS PAT PAR PCP PVN LA LV) );
```

The CFD coupling tutorials in `tutorials/fluids/circulationModels`
(`pipeWindkessel`, solved by `pimpleFoam` or, with `./Allrun solids4Foam`, by
`solids4Foam`, and `pipeClosedLoop`) use the simplest case of this, an
`externalPort` attached to a circuit, bound to a patch by name.

### What crosses the interface

Each `coupledDomain` creates two **ports**, `<name>.inlet` and `<name>.outlet`,
one per coupled patch. The valves stay in the 0D network. At each port:

- CFD gives: the port flow `Q`, positive in the natural direction of the
  patch (inlet: into the domain, `Qmv`; outlet: out of the domain, `Qav`)
- network returns: the pressure to impose on that patch

with the interface law (a bounded bidirectional inverse valve law):

```text
outlet   Ppatch = P_node + dP(Q)        inlet   Ppatch = P_node - dP(Q)
dP(Q) =  (Q/CQ)^2          Q >= 0
      = -(|Q|/(fr CQ))^2   Q <  0
fr = reverseFraction, clamped to [1e-3, 1]
```

C++ interface of `circulationNetwork` (all SI):

```c++
circulationNetwork net(dict, t0);
// Q = c0 + c1 a + c2 a^2 + c3 a^3, a = (t - t0)/(t1 - t0)
net.setPortForcing("LV.outlet", t0, t1, {c0, c1, c2, c3});
net.setPortForcing("LV.inlet",  t0, t1, {...});
net.advance(t1);
// Pa, at the state now
scalar Pav = net.portPressure("LV.outlet");
scalar Pmv = net.portPressure("LV.inlet");
```

The flow over a step is passed as a polynomial of up to third order, so that
the coupling schemes of up to fourth order can be built on it. The state of the
network is saved and restored with `stateDict()` / `readState()`.

The pipe case is the same with one port and no valve: a bare `externalPort`
`{ direction outlet; provider external; }` attached to a Windkessel. The
coupling algorithms (explicit or implicit, orders 1 to 4) differ only in *what forcing
polynomial they pass and when they integrate*; they are implemented in
`circulationCoupling` and chosen in the `coupling` dictionary.

### Implicit coupling

Implicit coupling re-integrates the network from the start of the step at
every new CFD flux, so the outer correctors of PIMPLE form a fixed-point
iteration on the interface flow. There is no relaxation and no stopping test
on the interface: PIMPLE's `residualControl` only measures the CFD residuals.
The coupler logs the relative change of the port flows at every exchange,

```text
circulationCoupler: exchange k, flow residual |Q^k - Q^(k-1)|/|Q^k| = ...
```

(the L2 norm over the ports; the first exchange of a step is compared with
the flow at the start of the step), and the last value of each step is the
`residual` column of `coupling.dat`. Choose `nOuterCorrectors` so that it is
small.

The boundary condition is evaluated when the pressure equation is assembled,
and `phi` is corrected after that pressure solve. Each exchange therefore uses
the flux of the previous pressure correction, and the state accepted at `t1` is
integrated with the flux of the last correction but one: the flux left by the
final pressure solve is not passed to the network in that step. It is the
first iterate of the next step, so the residual of the first exchange of a
step measures this lag for the step before.

It follows that implicit coupling of order 1 with one outer corrector (and
one pressure corrector) is effectively lagged explicit: its single exchange
holds the flux of the end of the previous step constant over the step, as
explicit order 1 does. Implicit coupling needs several outer correctors to be
implicit.

### Possible future work

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

## Status

The step-by-step tutorial and the verification scripts and reference data for
these models are not included yet. The tutorials in
`tutorials/fluids/circulationModels` are self-contained.
