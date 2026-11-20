# pipeWindkesselSolids4Foam: the Windkessel pipe solved by solids4Foam

The case of [`pipeWindkessel`](../pipeWindkessel/README.md), solved by
`solids4Foam` with the solids4foam `pimpleFluid` fluid model instead of stock
`pimpleFoam`. The pipe, the inflow `Q(t) = sin(2 pi t)` m3/s, the Windkessel
at the outlet and the coupling (explicit, order 2, `deltaT 0.05`) are the
same. It shows that the `circulationPressure` boundary condition works with
the solids4foam fluid models, and so with fluid-solid interaction cases.

## What differs from pipeWindkessel

```c++
// constant/physicsProperties
type    fluid;

// constant/fluidProperties
fluidModel pimpleFluid;
```

together with `constant/dynamicMeshDict` (`staticFvMesh`), `constant/g`, the
density in `constant/transportProperties` (the same 1060 kg/m3 as the
coupling `rho`), and `application solids4Foam;` in `system/controlDict`. The
circuit, `0/p` and `0/U` are those of `pipeWindkessel`.

## Running

```bash
# blockMesh, solids4Foam (under a minute), then the comparison with the exact
# solution by ../pipeWindkessel/analyse.py
./Allrun
# back to a clean case
./Allclean
```

## Results

The RMSE of the outlet pressure against the exact solution is 0.3385, as for
`pipeWindkessel`, and `postProcessing/circulationCoupling/0/coupling.dat`
agrees with that of the `pimpleFoam` case to the 12 digits written. See
[`pipeWindkessel`](../pipeWindkessel/README.md) for the exact solution and the
other coupling schemes and orders.
