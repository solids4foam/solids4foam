---
sort: 2
---

# Cantilever bent into a circle by an end rotation: `beamCantilever`

---

Prepared by Seevani Bali

---

## Tutorial Aims

- Demonstrate a geometrically exact beam analysis with the `beamSolid` solid
  model, which runs the beamFoam beam formulation through `solids4Foam`.
- Verify large in-plane rotations against an analytical solution.

## Case Overview

A straight cantilever of length $$L = 10$$ m is clamped at its left end. The
rotation of its right end is prescribed about the $$-y$$ axis and increases
linearly in pseudo-time from zero to $$\psi = 2\pi$$ at $$t = 1$$, so the beam
bends into a full circle. Both time-derivative schemes are `steadyState`, and
pseudo-time only sets the load increments (100 increments of
$$\Delta t = 0.01$$).

| Property                 | Value                 |
| ------------------------ | --------------------: |
| Length $$L$$             | 10 m                  |
| Cross-section            | circle, radius 0.2 m  |
| Young's modulus $$E$$    | 7.9577472e4 Pa        |
| Shear modulus $$G$$      | 3.9788736e4 Pa        |
| Beam segments            | 50                    |

The beam is meshed with `createBeamMesh` from `constant/beamProperties`. The
boundary conditions are:

- `left`: `fixedValue` displacement `W` and rotation `Theta` (clamped);
- `right`: zero force (`forceBeamDisplacementNR`) and a prescribed rotation
  (`fixedBeamRotation`, read from `constant/timeVsTheta`).

## Running the Case with solids4foam

The beam model, its solver settings, the cross-section and the material are
read from `constant/beamProperties`, as in beamFoam. The case selects the beam
solid model with three additions:

- `system/controlDict` loads the library: `libs ("libS4FBeamSolid.so");`
- `constant/physicsProperties` selects `type solid;`
- `constant/solidProperties` selects `solidModel beamSolid;`

`constant/g` is required by the solid model base class and is zero here.

```bash
./Allrun
```

## Expected Results

An end rotation $$\psi$$ bends the beam into a circular arc, and the free end
moves by

$$
w_x = -\left(L - \frac{L \sin\psi}{\psi}\right), \qquad
w_z = \frac{L\,(1 - \cos\psi)}{\psi}.
$$

At $$\psi = 2\pi$$ the free end returns to the clamped end, $$w_x = -L$$. The
free-end history is written to
`postProcessing/0/beamDisplacements_right.dat`. Over the whole history it stays
within $$2.4 \times 10^{-4} L$$ of the analytical position, and the results are
identical to those of standalone beamFoam. `regressionTest.sh` checks both.

To view the deformed beam in ParaView, apply **Warp By Vector** with `pointW`.

## References

- Simo, J.C., & Vu-Quoc, L. (1986). A three-dimensional finite-strain rod
  model. Part II: Computational aspects. _Computer Methods in Applied Mechanics
  and Engineering_, 58, 79–116.
- Bali, S., Taran, A., Tuković, Ž., Pakrashi, V., & Cardiff, P. (2025).
  beamFoam: A cell-centred finite volume solver for nonlinear
  geometrically-exact beams in OpenFOAM. _OpenFOAM Journal_, 5, 180–210.
  <https://doi.org/10.51560/ofj.v5.170>
