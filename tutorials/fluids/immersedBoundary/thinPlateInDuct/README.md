# thinPlateInDuct

**Experimental.** Check of `thinBodyBlocking` of `immersedBoundaryForce`
(`cutLink`); see "Experimental: thin bodies" in `src/immersedBoundary/README.md`.

A closed plate thinner than a cell spans a 2 m x 0.2 m duct (10 mm cells, one
cell in the empty direction, slip walls, open ends) and moves along it at
0.1 m/s from rest (impulsive start). An impermeable plate pushes the fluid at
its own speed, so the volume flux is exactly `w A` at every plane. Q is the
outlet flux; Q/(wA) is the time-average of the values of every time step
between 2 and 8 s, with a plate of 0.43 cells:

| `thinBodyBlocking` | Q/(wA), slip walls | Q/(wA), no-slip walls |
| --- | --- | --- |
| `off` (default) | 0.769 | 0.442 |
| `crossedFaces` | 1.0000 | 1.0000 |

`./Allrun` runs both cases (`./runCase RATIO WALLS MODE [NPROCS]` runs one,
for example `./runCase 0.5 noSlip crossedFaces 4` in four processes) and
checks the numbers (slip walls);
`RATIOS="0.5 1 3" WALLS="slip noSlip" ./Allrun` runs other thicknesses and
no-slip walls (`checkFlux.py CASE` prints the average).

`STATIC=1 PIN=0.01 END=6 ./runCase 0.43 slip crossedFaces` runs the plate at
rest under a pressure drop (outlet flux at the solver tolerance level with
the default `thinBodyDiffusivity zero`, see the README of the library).
The header of `runCase` lists the other options: `APERTURE=exact`,
`DIFFUSIVITY=keep`, `REFINE`, `YRANGES` (slots and plate edges), `STLFILE`,
`SHEAR` (non-orthogonal mesh), `X0`, `NPROCS`.

The test does not cover holes, gaps or free edges thinner than a cell,
oblique or curved thin shells, or a closed valve.
