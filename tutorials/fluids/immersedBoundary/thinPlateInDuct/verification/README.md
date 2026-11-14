# Verification cases of the thin-body option

Scripts that reproduce the numbers of the section "Experimental: thin bodies"
of `src/immersedBoundary/README.md` (OpenFOAM.com v2412; each script needs the
OpenFOAM environment sourced and solids4foam built, and writes its cases next
to itself). Everything here is experimental; the results are in the README,
with the failures. Serial, small meshes (seconds to an hour per case, the
`REFINE=4` openings take about an hour).

- `../runCase`, `../Allrun`, `../checkFlux.py`: moving and static plates;
  the options (APERTURE, DIFFUSIVITY, REFINE, SHEAR, NPROCS, ...) are listed in
  the header of `runCase`.
- `pocketCase`, analysed by `pockets.py` and `pockets.dat`: the pressure
  pocket and the topology change (plate A through the static plate B), with
  the policies and apertures.
- `gateCase`, analysed by `gate.py`: the sliding gate with a prescribed
  opening (`off`, `binary`, `exact`).
- `openings`, analysed by `openings.py`: slots and plate edges, static plate
  under a pressure drop.
- `shells`, `shellStl.py`, analysed by `shells.py`: the tilted and arc shells
  of 0.3 cells.
- `ale/runAle`, analysed by `ale/analyze.py`: the ALE test cases A-G (moving
  mesh, deforming body) with the switch off, binary zero and exact.

Not included: the 4-rank runs (use `NPROCS=4` in `../runCase`, with
`decomposeParDict` of the case), the cost measurements (the `search time` line
of the log and `/usr/bin/time -v` on `oscillatingCylinderInChannel` with
`MESH_LEVEL=3`), the non-orthogonal mesh (`SHEAR=0.1` in `../runCase`), and any
valve geometry.
