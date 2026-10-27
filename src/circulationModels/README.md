# circulationModels

Models of the circulatory system for solids4foam, from a mesh-free network of
lumped (0D) components, such as chambers, valves, resistances and compliances,
to the coupling of that network with a 3D CFD domain. It is currently only
available for OpenFOAM.com (tested with v2412), where the libraries are built
with solids4foam. 1D vessel models are planned and will live here too.

## Contents

| Directory | Library | Description |
|---|---|---|
| `circulation0D` | `libcirculation0D` | Lumped-parameter (0D) circulation network: components, waveforms, per-cycle metrics, presets, and the ODE solution of the network. See [circulation0D/README.md](circulation0D/README.md) |
| `circulationCoupling` | `libcirculationCoupling` | Coupling of a circulation network to a CFD domain: the `circulationPressure` boundary condition and the coupler (explicit/implicit, orders 1 to 4) |
| `circulationDdtSchemes` | `libcirculationDdtSchemes` | The `BDF3` and `BDF4` `ddtScheme`s, used by the CFD time discretisation at higher coupling orders. See [circulationDdtSchemes/README.md](circulationDdtSchemes/README.md) |
| `tools` | | `circulationSchematic` (draws a circuit as an SVG), `circulationPlot` (plots result columns, needs gnuplot) and `addHeader` (adds the solids4foam banner to case dictionaries) |
| `docs` | | A step-by-step tutorial and the description of the 0D/3D coupling |

The `circulationCoupling` and `circulationDdtSchemes` libraries are not specific
to the 0D models. The coupler talks to the network through ports, and the BDF
schemes are general time schemes.

The `BDF3` and `BDF4` schemes here are first-derivative `ddtScheme`s, for
`fvm::ddt`/`fvc::ddt`. They are unrelated to the solids4foam `bdfD2dt2` scheme,
which is a second-derivative scheme (`d2dt2Schemes`) for solid dynamics.

The applications are in:

* `applications/solvers/circulationFoam`: solves a circulation network, with no mesh.
* `applications/test/testPortCoupling`: tests the coupling API by driving a
  network from tabulated port flows.

The tutorials are in `tutorials/circulationModels`.

## Building

The libraries are built by `src/Allwmake` and the applications by
`applications/Allwmake`, along with the rest of solids4foam. To build only this
module, after sourcing OpenFOAM:

```bash
(cd src/circulationModels && ./Allwmake)
(cd applications/solvers/circulationFoam && ./Allwmake)
(cd applications/test/testPortCoupling && ./Allwmake)
```

The libraries are written to `$FOAM_MODULE_LIBBIN`, like the other solids4foam
libraries. For a CFD case, load the coupling library in `system/controlDict`,
and the BDF schemes if used:

```c++
libs ("libcirculationCoupling.so" "libcirculationDdtSchemes.so");
```

## Quick start

```bash
cd tutorials/circulationModels/healthyHuman
./Allrun    # 20 cardiac cycles of a 4-chamber closed loop in under a second
# results: postProcessing/circulation/0/circulation.dat
```

Add `src/circulationModels/tools` to your `PATH` so that `circulationFoam -diagram`
can call `circulationSchematic` to draw the circuit.
