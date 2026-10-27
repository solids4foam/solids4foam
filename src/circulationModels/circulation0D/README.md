# circulation0D

A lumped-parameter (0D) model of the circulatory system for OpenFOAM v2412,
built from named components that you connect the way you would draw the
circuit: *a resistor connected to a capacitor connected to ...*. It is solved
with OpenFOAM's own ODE solvers (`RKDP45`, `Rosenbrock34`, `seulex`, ...), needs
no mesh, and is part of the solids4foam `circulationModels` module, see
[../README.md](../README.md). Any lumped piece of the
circuit, such as the left ventricle, can be swapped for a 3D CFD domain by
changing one entry; see [../docs/COUPLING.md](../docs/COUPLING.md).

```
LV -> AV -> SAS -> SAT -> SAR -> SCP -> SVN -> RA -> TV -> RV -> PV
   -> PAS -> PAT -> PAR -> PCP -> PVN -> LA -> MV -> LV            (healthyHuman)
```

## Quick start

```bash
# after building the module, see ../README.md
cd tutorials/circulationModels/healthyHuman
circulationFoam                     # 20 cardiac cycles in under a second
# results: postProcessing/circulation/0/circulation.dat
```

## Writing a circuit: `constant/circulationProperties`

```c++
units       mmHg_mL;        // SI (default) or mmHg_mL; time is always seconds
Tcycle      0.8333;         // cardiac period, the default for the chambers

components
{
    ground  { type ground; }
    Qin     { type flowSource;
              waveform { type sine; mean 100; amplitude 100; period 0.8333; } }
    rc      { type resistor;   R 0.05; }
    C       { type compliance; C 1.1; initialPressure 73; }
    R       { type resistor;   R 0.9; }
}

connections ( (ground Qin rc C R ground) );     // a chain: A.out -> B.in -> ...

solver { type RKDP45; absTol 1e-12; relTol 1e-10; }
probes ( Qin.out );             // optional: also report the pressure at a terminal
```

* **Connections.** `(A B C)` joins `A.out` to `B.in` and `B.out` to `C.in`. Use
  explicit terminals when you need them, `(Rc.out Rp.in)` and `(Rc.out C.in)`
  for a junction. Every component has terminals `in` and `out`; flow is
  positive from `in` to `out`. Terminals joined together share one pressure.
* **Units.** Write the case in the units of the paper you are copying from.
  Everything is converted to SI on reading, and the results are written in the
  same units as the case.
* **Errors** name the component and the line: unknown type (with the list of
  valid types), unknown terminal, dead end, non-positive parameter, a node with
  no pressure reference.

### Components

| `type` | what it is | parameters |
|---|---|---|
| `resistor` | `Q = dP/R` | `R` |
| `inductor` | `L dq/dt = dP`; flow is a state | `L`, `initialFlow` |
| `compliance` | storage node, `P = (V-V0)/C + P0` | `C`, `V0`, `P0`, `initialPressure` or `initialVolume` |
| `elastanceChamber` | heart chamber, `P = (V-V0)E(t) + P0`, half-cosine activation | `V0 P0 Emin Emax Tes Ted Eshift T initialVolume` |
| `mixedChamber` | `P = A(t)(E1 dV + E2 dV^2) + P0 (exp(K dV) - 1)`, `A` = sin^2 or half-cosine window | `V0`, `active {E1 E2}`, `passive {P0 K}`, `activation {...}`, `initialVolume` |
| `orificeValve` | `Q = CQ sqrt(dP)` forward, closed otherwise | `CQ` or `Rhat` (`dP = Rhat Q^2`), `R` (resistance in series), `epsilon` (smooth switch), `reverseFraction` (leak) |
| `flowSource` | prescribed flow in -> out | `waveform {...}` |
| `pressureSource`, `ground` | node at prescribed / fixed pressure | `waveform {...}` / `P` |
| `CR`, `CRL` | `in -> C -> R (-> L) -> out` vessel segment | `C R L` and the compliance options |
| `valvedChamber` | inlet valve + chamber + outlet valve | `chamber {} inletValve {} outletValve {}` |
| `coupledDomain` | external (3D) replacement for `valvedChamber` | `provider`, `inlet {} outlet {}` |

Waveforms: `constant`, `sine`, `halfSine`, `table` (whitespace file,
`column`, optional `periodic`).

Results for `CRL`/`CR` parts appear as `SAS.C.P`, `SAS.C.V`, `SAS.R.Q`,
`SAS.L.Q`; chambers as `LV.P`, `LV.V`; valves and resistors as `AV.Q`.

Adding a component is one class deriving from `circulationComponent` with
`TypeName` and `addToRunTimeSelectionTable`, with no other file touched.

### How it is solved

Connected terminals form *nodes*. A node holding a compliance or chamber has a
volume state. Every other node is solved for its pressure so that the flows
into it balance (Newton, one step when only resistors meet there). Inductors
add a flow state. That gives an explicit ODE `dy/dt = f(t,y)` handed to the
`ODESolver` of your choice. A group of nodes with no pressure reference whose
flow is carried by one inductor, such as the 4-element Windkessel
`source -> r -> L -> C||R`, is handled by pinning it with `L dQ/dt`.
Implicit solvers use a finite-difference Jacobian; with hard-diode valves
prefer `RKDP45`, or set `epsilon > 0` on the valves.

## Presets and your own modules

A sub-circuit can be a one-liner. These ship with the library, and any parameter is
named to set it:

| `type` | what | parameters |
|---|---|---|
| `windkessel2`, `windkessel3`, `windkessel4` | the Windkessels, ending on their own ground; terminal `in` | `R C` (+ `Rch`, `L`), `initialPressure` |
| `healthyHuman` | the closed 4-chamber healthy-human loop (`units mmHg_mL`; Healy 2024, after Korakianitis & Shi 2006) | every chamber, valve and vessel value, e.g. `Emax_lv 2.0;` |
| `pipeLoopHeart` | the heart and arteries of the closed-loop pipe study (SI); terminals `in`, `out` | the study's own names (`C1 R1 ... Ea Ev1 ...`) |

```c++
components
{
    outlet { type externalPort; direction outlet; patch outlet; }
    wk     { type windkessel3; R 10; Rch 1; C 0.01; }    // the whole outlet model
}
connections ( (outlet wk) );
```

Write your own in a `modules` dictionary, with `@name` where a parameter goes. A module
can use other modules, and a case's module replaces a preset of the same name:

```c++
modules
{
    myOutlet
    {
        required    ( R C );
        defaults    { initialPressure 0; }
        components  { wk { type windkessel2; R @R; C @C; initialPressure @initialPressure; } }
        terminals   { in wk.in; }
    }
}
```

The same circuit from a preset gives identical numbers to the hand-written one.

## Checking and drawing a set-up

```bash
circulationFoam -check      # elements, nodes, states, ports, solver; warnings; then exit
circulationFoam -diagram    # the circuit as an electrical circuit diagram: circuit.json, circuit.svg
```

`-check` works on any case with a `circulationProperties`, CFD cases included. It lists
every element and node, the ODE states with their initial values, the coupling ports
and, when the case has a mesh, whether each port's patch exists. It warns about parameter
values outside human-scale ranges, which is how a unit mix-up shows up (SI numbers in a
`units mmHg_mL;` case, or mmHg values in an SI one). `-diagram` draws the circuit as an electrical circuit diagram, in the electrical analogy of the model: a
resistor is a zigzag, an inductor a coil, a valve a diode, a flow source a circle with an arrow, a
compliance a capacitor to ground, a heart chamber a variable capacitor, a fixed pressure a ground
symbol, and a CFD port a blue box. Closed loops are drawn as a ring, chains left to right. Each
label carries the element's name and value in the case's units. The drawing is made by
`tools/circulationSchematic` (put `src/circulationModels/tools` on the PATH; it needs only Python 3; Graphviz,
if installed, is used for the node positions). Open the SVG with `open circuit.svg`.

## Dimensioned parameters

A parameter may carry OpenFOAM dimensions, which are checked against what it must be:
`R [1 -4 -1 0 0 0 0] 10;`. A dimensioned value is SI, whatever the case's `units`; a plain
number is in the case's units. A wrong dimension is rejected with OpenFOAM's own message.

## Per-cycle metrics

Name output columns in a `metrics` dictionary and get one row per cardiac cycle in
`metrics.dat` (and the latest pressure-volume loop in `<name>.pvLoop.dat`), in 0D and CFD
runs alike:

```c++
metrics
{
    LV  { type chamber;  volume heart.LV.V; pressure heart.LV.P; }   // EDV ESV SV EF CO Pmax Pmin
    SAT { type pressure; column heart.SAT.C.P; }                     // mean systolic diastolic pulse
    AV  { type flow;     column heart.AV.Q; }                        // mean max min net/fwd/back volume, regurgitant fraction
}
```

Cycles are `[k Tcycle, (k+1) Tcycle]`, and only complete ones are reported (a run that
starts or stops mid-cycle skips the partial one; a CFD run started at t = 0 reports from
its second cycle, because the initial state is not an accepted step). The healthy-loop
values agree with those of an independent reference implementation to 2e-6 relative.

## Restart

At every output time the 0D solver writes `<time>/circulationState`, and the CFD coupler
`<time>/circulationState-circulationProperties` (the state, the flow history and the ODE
step estimate, at full precision). A case that starts from a time holding one resumes the
circuit with the fields. A CFD run resumed this way continues to 1e-15 of the continuous
run, with identical coupling iterations (with time schemes that need only the old time, e.g.
Euler; `backward` and BDF also need their extra old-time fields, which is a CFD matter).
A resumed 0D run agrees on the states to 3e-4 % of range; it is not bit-exact because the
adaptive solver places each hard-valve switching instant slightly differently.

## When something goes wrong

* **ODE solver step cap**: the message says which interval failed, lists the usual causes
  (an unstable coupling, tolerances too tight for hard valves) and prints the state.
* **A floating point exception in the CFD solution** (the usual sign of an unstable coupling,
  e.g. explicit order 1 at a large step) is preceded by an explanation naming the scheme, the
  order, the time step and the last port pressure.
* **Divergence**: if a port pressure becomes non-finite or exceeds `coupling { maxPressure 1e8; }`
  (Pa) the run stops with advice.
* **A CFD time scheme of lower order than the coupling** (for example `backward` with order 3)
  prints a warning at start-up, since it limits the coupled solution.
* **Starting from rest**: `coupling { seed zero; }` states that every earlier flow is zero,
  so orders 3 and 4 run at full order from the first step instead of ramping up. (A flow
  that starts abruptly, like the sine inlet of the Windkessel test, needs its exact earlier
  values as a `seed {...}` list.)

## Coupling to CFD

A patch of a `pimpleFoam` case is coupled by naming it in the circuit and using one
boundary condition:

```c++
// constant/circulationProperties
outlet { type externalPort; direction outlet; patch outlet; }     // a component
coupling { scheme implicit; order 2; rho 1060; }                  // explicit|implicit, 1..4

// 0/p
outlet { type circulationPressure; value uniform 0; }
```

and `libs ("libcirculationCoupling.so");` in `controlDict`. See
`tutorials/circulationModels/pipeWindkessel`.

## Tutorials

| case | what |
|---|---|
| `healthyHuman` | 4-chamber closed loop, healthy-human parameters |
| `pipeClosedLoop` | a CFD pipe with **both** ends coupled to a closed circulation loop (atrium, ventricle, two valves, 7 states): same scheme/order choice; verified against the original runs and an independent RK4 reference |
| `pipeClosedLoopMonolithic` | the same circuit with the pipe replaced by its R-L surrogate, no CFD (the reference for the above) |
| `pipeWindkessel` | a CFD pipe (`pimpleFoam`) coupled to a Windkessel; `scheme explicit/implicit` and `order 1..4` chosen in the dictionary, verified against the exact solution and the original hand-written couplings |

## Verification

During development the models were checked against reference solutions: a
first-principles Windkessel, an independent implementation of the healthy-human
loop, and, for the coupled CFD cases, exact solutions or a monolithic
reference (see the READMEs in `tutorials/circulationModels`). The verification
scripts and the reference data are not included yet.

Last-cycle healthy loop: LV EDV/ESV 141.7/55.8 mL, SV 85.9 mL, EF 60.6 %, aortic
129/89 mmHg, MAP 109.5. Note that this healthy parameter set is itself hypertensive and
high-venous-pressure, so check parameters before trusting physiology.

## Layout

```
src/circulationModels/circulation0D/           this library (context, waveform, components, network)
src/circulationModels/circulationCoupling/     CFD coupling: circulationPressure BC and the coupler (scheme, order)
src/circulationModels/circulationDdtSchemes/   BDF3 / BDF4 time schemes
applications/solvers/circulationFoam, applications/test/testPortCoupling
tutorials/circulationModels/                   cases above
docs/COUPLING.md         swapping 0D for 3D (in ../docs)
```

## Known limits

* `epsilon = 0` valves are non-smooth; expect ~0.03 % of range noise at
  switching instants regardless of tolerance.
* A node with more than one compliance or chamber is rejected (give each its
  own node).
* The parameter-range warnings of `-check` assume a human-scale circulation; other
  scales (the unit-flow test pipe, say) give harmless warnings.
* Metrics are configured in the dictionary, not as an OpenFOAM function object.
* The Shi dynamic-leaflet valve, passive exponential chambers and
  pressure-driven ports are planned, not built.
