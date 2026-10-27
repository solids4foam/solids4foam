# pipeClosedLoopSample: a ready-to-run demo

A 3D pipe (4 cm x 30 cm, 24,675 cells) whose **both ends** are coupled to a lumped closed
circulation: the outlet drains into the arterial side (C1, R1, L1, C2, R2, C3, R3), the heart
(atrium, mitral valve, ventricle, aortic valve) pumps, and the blood returns to the pipe inlet (C4).

```
outlet -> C1 -R1-L1- C2 -R2- C3 -R3- atrium -mitral- ventricle -aortic(+R4)- C4 -> inlet
```

The CFD gives the flow through each patch; the circuit returns the pressure of the node each
patch is attached to. The circuit is one dictionary, `constant/circulationProperties`.

## Run it

```bash
./Allrun            # blockMesh, circuit check, decomposePar, pimpleFoam on 8 ranks (./Allrun 4 for 4)
python3 plotter.py  # figures open on screen and are saved to plots/   (--no-show: only save)
./Allclean          # back to a clean case
```

Only `pimpleFoam` is needed: the coupling library (`libs` in `system/controlDict`) builds the
circuit inside the solver and integrates its ODEs at every exchange. `circulationFoam` is for
circuit-only cases; here it is used only for `circulationFoam -check` (what the program understood,
and that both patches exist in the mesh). `plotter.py` needs matplotlib (`pip3 install matplotlib`).

Settings for this demo: one cardiac cycle (`endTime 0.8`), `deltaT 5e-3`, implicit coupling of
order 2. About a minute on 8 free cores (the three-cycle run about 2.5 minutes); watch it with `tail -f log.pimpleFoam`.

## Things to show

* `circulationFoam -check`: the circuit as understood (7 states, 7 nodes, 2 coupling ports).
* `constant/circulationProperties`: the whole circuit; change `scheme` / `order` in `coupling`.
* `0/p`: `circulationPressure` on `inlet` and `outlet`; the only change from an ordinary pipe case.
* Results: `postProcessing/circulationCoupling/0/coupling.dat`, one row per time step (the flow
  and imposed pressure at each port, then every element's pressure, volume and flow).
* `plots/interface.png`, `heart.png`, `circuit.png`.

## Variations

| to try | change |
|---|---|
| three cycles (the periodic state) | `endTime 2.4` in `system/controlDict`, then `python3 plotter.py` (last cycle) |
| explicit coupling (cheaper, less stable) | `scheme explicit;` |
| higher-order coupling | `order 3;` plus `libs ("libcirculationCoupling.so" "libcirculationDdtSchemes.so");` and `default BDF3;` in `fvSchemes` |
| a stiffer ventricle | `ventricle { ... active { E1 3.0e9; ... } }` |

The wider set of schemes and orders is in `../README.md` and `../makeCase`, `../runLadder`.
