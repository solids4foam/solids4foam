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
- `tools` (no library): `circulationSchematic` (draws a circuit as an SVG),
  `circulationPlot` (plots result columns, needs gnuplot) and `addHeader`
  (adds the solids4foam banner to case dictionaries)
- `docs` (no library): The description of the 0D/3D coupling,
  [docs/COUPLING.md](docs/COUPLING.md)

The `circulationCoupling` library is not specific to the 0D models: the
coupler talks to the network through ports.

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

The tutorials are in `tutorials/circulationModels`.

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
cd tutorials/circulationModels/healthyHuman
./Allrun    # 20 cardiac cycles of a 4-chamber closed loop in under a second
# results: postProcessing/circulation/0/circulation.dat
```

Add `src/circulationModels/tools` to your `PATH` so that
`circulationFoam -diagram` can call `circulationSchematic` to draw the circuit.

## Status

The step-by-step tutorial and the verification scripts and reference data for
these models are not included yet. The tutorials in
`tutorials/circulationModels` are self-contained.
