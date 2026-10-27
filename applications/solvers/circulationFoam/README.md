---
sort: 5
---

# `circulationFoam` solver

`circulationFoam` solves a lumped-parameter (0D) model of the circulatory
system: named components, such as heart chambers, valves, resistances and
compliances, joined by connections, integrated in time with an OpenFOAM ODE
solver. It needs no mesh. It is part of the
[`circulationModels`](../../../src/circulationModels/README.md) module and is
currently only built for OpenFOAM.com.

The circuit is a dictionary, `constant/circulationProperties`. Its syntax, the
available components and the presets are described in
[`circulation0D`](../../../src/circulationModels/circulation0D/README.md).

---

## Usage

```bash
circulationFoam              # solve the circuit
circulationFoam -check       # describe and check the set-up, then exit
circulationFoam -diagram     # draw the circuit, then exit
```

A case needs only two files, as there is no mesh:

* `constant/circulationProperties`: the circuit, i.e. the components, the
  connections, the ODE solver and, optionally, `metrics`.
* `system/controlDict`: `startTime`, `endTime`, `deltaT` and the write times.

`deltaT` is the **output interval**: the ODE solver takes whatever internal
steps its tolerances need inside it. See
[`tutorials/circulationModels/healthyHuman`](../../../tutorials/circulationModels/healthyHuman)
for a complete case.

---

## Options

* `-check`: nothing is solved. The circuit is described (elements, nodes,
  states, coupling ports and the ODE solver) and the set-up is checked,
  including the magnitude of the parameters. If the case has a mesh, it also
  checks that the patch of each coupled port exists. It works on any case
  with a `constant/circulationProperties`, CFD cases included.
* `-diagram`: writes `circuit.json` and, if `circulationSchematic` is on the
  `PATH`, draws the circuit as an electrical circuit diagram, `circuit.svg`.
  Then exits.

`circulationSchematic` is in
[`src/circulationModels/tools`](../../../src/circulationModels/tools); it
needs only Python 3.

---

## Output

Every `deltaT`, a row is appended to
`postProcessing/circulation/<startTime>/circulation.dat`, with the pressure,
volume and flow of every element, in the units the case is written in. The
first lines name the units and the columns:

```text
# circulationFoam, units <units>
# t <column> <column> ...
# s <unit> <unit> ...
```

If the circuit has a `metrics` dictionary, the per-cycle metrics are written
to the same directory. Columns can be plotted by name with `circulationPlot`
(in the same `tools` directory; it needs gnuplot).

---

## Restart

At each write time the state is written, at full precision, to
`<time>/circulationState`. If this file exists in the start time, the solver
reads it, so a long spin-up can be continued.

---

## CFD coupling

`circulationFoam` solves circuits on their own. A circuit coupled to a CFD
domain is solved inside the CFD solver (for example `pimpleFoam`), through the
`circulationCoupling` library. See the tutorials `pipeWindkessel` and
`pipeClosedLoop` in
[`tutorials/circulationModels`](../../../tutorials/circulationModels).
`circulationFoam -check` still works on those cases, to see what was
understood from the circuit and that the coupled patches are in the mesh.

---

## Building

The solver is built by `applications/Allwmake` along with the rest of
solids4foam. It requires `libcirculation0D`, built by `src/Allwmake`. On
other OpenFOAM flavours, the build is skipped. To build it alone, after
sourcing OpenFOAM and building `src/circulationModels`:

```bash
cd applications/solvers/circulationFoam
./Allwmake
```

---

## Source

[`applications/solvers/circulationFoam`](https://github.com/solids4foam/solids4foam/tree/development/applications/solvers/circulationFoam)
