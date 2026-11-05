# thinPlateInDuct

Check of `thinBodyBlocking` of `immersedBoundaryForce` (cutLink): a closed
plate 0.43 cells thick spans a 2 m x 0.2 m duct (10 mm cells, one cell in the
empty direction, slip walls, open ends) and moves along it at 0.1 m/s. An
impermeable plate pushes the fluid at its own speed, so the volume flux is
exactly `w*A` at every plane. Mean over 2-6 s, Q/(wA):

| `thinBodyBlocking` | Q/(wA) |
|---|---|
| `off` (default) | 0.67 (oscillating with the position of the plate in the mesh) |
| `crossedFaces` | 1.0000 |

`./Allrun` runs both and checks the numbers. This test does not cover
holes, gaps or free edges thinner than a cell, oblique or curved thin shells
(which still transmit less than their swept volume during the start-up),
no-slip walls (a wall-layer deficit of about 1-5 % of the flux remains, also
for thick plates) or a closed valve that has to be gap-free.
