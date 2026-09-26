---
sort: 7
---

# Immersed heart valve in a duct: `heartValveInDuct`

You can find the files for this tutorial under
[`tutorials/fluids/immersedBoundary/heartValveInDuct`](https://github.com/solids4foam/solids4foam/tree/master/tutorials/fluids/immersedBoundary/heartValveInDuct).

---

## Tutorial Aims

- Demonstrate the heart valve motions of the `immersedBoundaryForce` finite
  volume option: `valveSliceAxis`, the collapse of the leaflets towards the
  valve axis, and `valveMorph`, the interpolation between a closed and an
  open valve.

## Case Overview

A synthetic valve with three leaflets of thickness 1.5 mm sits in a square
duct of side 24 mm, from $$z = -10$$ mm to $$z = 40$$ mm, on an annular plate
at $$z = 0$$ with a circular orifice of radius $$R = 10$$ mm, which closes the
duct around the valve. The flow, of kinematic viscosity
$$10^{-5}$$ m$$^2$$/s, is driven by a constant kinematic pressure drop of
0.05 m$$^2$$/s$$^2$$ between the inlet and the outlet, and passes through the
valve. The surfaces in `constant/triSurface` (binary STL) are written by
`makeValveStl.py`: they are closed, with a non-zero thickness.

With `VALVE_MOTION=sliceAxis` (the default), the leaflets
(`valveTube.stl`) lie on the cylinder of radius $$R$$ from the annulus to
$$z = 16$$ mm, and collapse towards the axis of the valve with the
`valveSliceAxis` motion: each section normal to the axis shrinks towards the
axis by the fraction $$F = F_{max} \, g(\xi) \, G(t)$$, where $$\xi$$ is the
distance from the annulus, the space law $$g$$ rises from 0 at the annulus to
1 at 8 mm, and the time law $$G$$ rises from 0 (open) to 1 (closed) during
the closing window of the cycle of period 0.8 s, $$0.3 < t/T < 0.4$$, and
falls back to 0 during the opening window, $$0.8 < t/T < 0.95$$. The vertices
of the leaflets on the annulus (`annulus.stl`) are held fixed:

```c++
valve
{
    surface     "valveTube.stl";

    motion
    {
        type            valveSliceAxis;
        bottomCentre    (0 0 0);
        topCentre       (0 0 0.016);
        fixedEnd        bottom;
        Fmax            0.95;
        spaceLaw        smootherstep;
        nearRamp        0.008;
        activeLen       0.016;

        period          0.8;
        timeLaw         twoWindow;
        closeWindow     (0.3 0.4);
        openWindow      (0.8 0.95);

        annulusSurface   "annulus.stl";
        annulusTolerance 0.0003;
    }
}
```

With `VALVE_MOTION=morph`, the leaflets are flat and close the orifice
(`valveClosed.stl`), and open with the `valveMorph` motion towards
`valveOpen.stl`, the same leaflets rotated by 80 degrees about their hinges
on the annulus, with the same vertices; the valve opens during
$$0.05 < t/T < 0.15$$ and closes during $$0.35 < t/T < 0.45$$.

The velocity of the leaflets at a point is interpolated from the velocities
of the vertices of the nearest triangle of their surface. Neither motion
preserves the volume of the leaflets: the `pimpleFluid` fluid model does not
impose continuity in the cells entirely inside a body.

The motions follow the valve kinematics of the immersed boundary library
contributed to cardiacFoam by Sairam Pamulaparthi Venkata (see
`src/immersedBoundary/README.md`), here applied to synthetic valves.

## Running the Case

```bash
./Allrun
```

or `VALVE_MOTION=morph ./Allrun`, with `MESH_LEVEL=2 ./Allrun` etc. for the
finer meshes (0.5 mm cells for `MESH_LEVEL=1`, about 230 000 cells). The
flow rate through the outlet is written every time step to
`postProcessing/flowRate/0/surfaceFieldValue.dat`, and plotted in
`flowRate.pdf` if gnuplot is installed.

## Expected Results

The largest flow rate through the open valve, and the flow rate through the
closed valve, for `MESH_LEVEL=1`:

| Motion | Open (ml/s) | Closed (ml/s) |
| ------ | ----------- | ------------- |
| `sliceAxis`, first cycle | 61 | 17 |
| `sliceAxis`, second cycle | 76 | 17 |
| `morph` | 56 | 3 |

The flow through the open valve is still developing when the valve closes in
the first cycle. While the valve closes, the flow rate through the outlet is
briefly negative (about -13 ml/s). With `sliceAxis`, the leaflets collapse
into a cone that meets near the axis at about 6 mm from the annulus; beyond the
ramp of the space law, the sections shrink to 5% of their radius, and the
leaflets there, whose thickness shrinks with them, are thinner than the
cells, so the closed valve leaks 20-30% of the open flow rate on this
mesh. The flat closed leaflets of `morph` leak about 6%, through the
gaps between them.

The regression test runs `MESH_LEVEL=0` (1 mm cells) to $$t = 0.4$$ s, and
checks the flow rates through the open valve at $$t = 0.2$$ s and through
the closed valve at $$t = 0.4$$ s (57 and 22 ml/s).
