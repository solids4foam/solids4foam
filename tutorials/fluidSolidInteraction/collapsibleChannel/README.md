---
sort: 10
---

# Flow in a collapsible channel: `collapsibleChannel`

---

## Tutorial Aims

- Demonstrates a strongly coupled internal-flow fluid-solid interaction case
  with a massless, thin elastic wall.
- Verifies the partitioned IQN-ILS coupling against a converged monolithic
  solution computed with the open-source finite-element library
  [oomph-lib](https://oomph-lib.github.io/oomph-lib/).
- Compares the second-order and the high-order (cubic) solid discretisations
  on a thin wall that deforms in combined bending and stretching.

## Case Overview

Viscous fluid flows through a two-dimensional channel of width
$$a = 1\,\mathrm{m}$$. Part of the upper wall, $$5 < x < 15\,\mathrm{m}$$, is a
thin elastic wall of thickness $$h = 0.05\,\mathrm{m}$$, clamped at both ends
to rigid channel sections of length $$L_{up} = 5\,\mathrm{m}$$ upstream and
$$L_{down} = 10\,\mathrm{m}$$ downstream. The flow is driven by the inlet
pressure $$p_{in} = 12 \nu U (L_{up} + L + L_{down})/a^2 = 6\,\mathrm{m^2/s^2}$$
(kinematic) that sustains a Poiseuille flow of mean velocity
$$U = 1\,\mathrm{m/s}$$, with $$p = 0$$ at the outlet and parallel flow at both
ends. The fluid starts as that Poiseuille flow, the wall starts flat, and an
external pressure on the wall is ramped up as
$$p_{ext}(t) = 200\,(1 - \cos(\pi t/0.25))/2\,\mathrm{Pa}$$ over the first
$$0.25\,\mathrm{s}$$. The wall collapses into the channel, overshoots, and
performs a damped oscillation of period about $$0.8\,\mathrm{s}$$ about its
new equilibrium, pumping fluid in and out of the channel as it does so.

| Quantity | Value |
| --- | --- |
| Fluid density, kinematic viscosity | 1 kg/m^3, 0.02 m^2/s |
| Reynolds number $$Re = U a/\nu$$ (and $$Re\,St$$) | 50 |
| Wall thickness $$h$$, length $$L$$ | 0.05 m, 10 m |
| Wall Young's modulus $$E$$, Poisson's ratio $$\nu_s$$ | 455 MPa, 0.3 |
| Plane-strain wall modulus $$E/(1-\nu_s^2)$$ | 500 MPa |
| Wall density | massless (quasi-static wall) |
| External pressure $$p_{ext}$$ | 200 Pa, ramped over 0.25 s |
| Time step, end time | 0.025 s, 3.5 s |

The monitored quantity is the vertical displacement of the wall at 25%, 50%
and 75% of its length, $$x = 7.5$$, $$10$$ and $$12.5\,\mathrm{m}$$, written by
the `solidPointDisplacement` function objects `wallQuarter`, `wallMid` and
`wallThreeQuarter`.

### Relation to the published case

The configuration is that of the oomph-lib tutorial
[Flow in a 2D collapsible channel](https://oomph-lib.github.io/oomph-lib/doc/interaction/fsi_collapsible_channel/html/index.html),
itself a version of the problem studied by Jensen & Heil (2003), Heil (2004)
and Heil, Hazel & Boyle (2008): the same geometry, Reynolds number, flow
driving, boundary conditions, massless wall, initial condition and monitored
wall displacement. Three things are changed, each so that a finite-volume
continuum can represent the wall:

1. **No pre-stress.** In the published cases the wall is a Kirchhoff-Love
   beam with an axial pre-stress $$\sigma_0 = 10^3$$ on the scale of its own
   modulus, i.e. a pre-strain of a thousand. A continuum cannot carry that as
   strain. solids4foam can impose an initial stress `sigma0` on a hyperelastic
   law, but the high-order solid supports no law that carries it, and the
   standard solid needs the geometric stiffness in its implicit operator, so
   this tutorial uses $$\sigma_0 = 0$$. The wall then deforms in combined
   bending and stretching, rather than as the tensioned membrane of the
   published cases, whose bending stiffness is $$10^{-8}$$ of its tension.
2. **A thicker wall, $$h/a = 1/20$$,** the value of Heil, Hazel & Boyle
   (2008), and an external pressure scaled to keep $$p_{ext}/(\mu U/a) = 10^4$$,
   as in the oomph-lib tutorial. The collapse at equilibrium is then about 14% of
   the channel width, against 12% there.
3. **Clamped ends and a ramped external pressure.** A continuum wall fixed on
   its end faces is clamped, not pinned; the ramp replaces the impulsive
   start, which limits the time-step convergence of the second-order schemes.

The non-dimensional parameters, on the oomph-lib scales (lengths on $$a$$,
velocities on $$U$$, stresses on $$E/(1-\nu_s^2)$$), are $$Re = Re\,St = 50$$,
$$Q = \mu U/(a E_{eff}) = 4\times10^{-11}$$, $$H = h/a = 0.05$$,
$$\sigma_0 = 0$$, $$P_{ext} = 4\times10^{-7}$$ and
$$P_{up} = 12(L_{up} + L + L_{down}) = 300$$.

### From the beam to a continuum

The reference wall is oomph-lib's geometrically nonlinear Kirchhoff-Love beam
with an incrementally linear constitutive law: the axial second Piola-Kirchhoff
stress is $$E_{eff}\,\gamma$$, with $$\gamma$$ the Green strain of the
midplane. The solids4foam wall is a two-dimensional (plane-strain)
St Venant-Kirchhoff continuum of the same thickness, whose second
Piola-Kirchhoff stress is linear in the Green strain. In plane strain with a
traction-free top and bottom, its axial stiffness and bending stiffness per
unit width are $$E h/(1-\nu_s^2)$$ and $$E h^3/(12(1-\nu_s^2))$$, which equal the
beam's $$E_{eff} h$$ and $$E_{eff} h^3/12$$ when $$E = E_{eff}(1 - \nu_s^2)$$.
The fluid acts on the bottom face of the continuum, at $$y = 1$$, and on the
midplane of the beam, also at $$y = 1$$, so both see the same fluid domain.

The continuum converges to the beam only as $$h/L \to 0$$. Here
$$h/L = 0.005$$: the neglected transverse shear, the offset between the loaded
face and the midplane, and the clamping of the whole end face rather than the
midplane are all of order $$h/L$$ or smaller, below the tolerances used in the
verification.

## Solid discretisations

Both solid variants use the total-Lagrangian solid model solved with PETSc
SNES (Jacobian-free Newton-Krylov) and the least-squares reconstruction of the
high-order finite-volume method, preconditioned with its assembled Jacobian:

- `./Allrun` uses a cubic reconstruction, `polynomialOrder 3`
  (`constant/solid/solidProperties.highOrder`);
- `./Allrun linear` uses a linear reconstruction, `polynomialOrder 1`, i.e. a
  second-order discretisation (`constant/solid/solidProperties.linear`).

The high-order solid is the default because the wall bends. Solved alone
under the full external pressure, the wall-midpoint deflection of the clamped
oomph-lib beam is $$-0.14363\,\mathrm{m}$$, and the continuum gives:

| Solid mesh (length x thickness) | Linear | Cubic |
| --- | ---: | ---: |
| 80 x 8 (tutorial) | -16.7% | +0.21% |
| 160 x 8 | -6.0% | +0.08% |
| 320 x 8 | -1.1% | +0.06% |
| 640 x 8 | -0.22% | +0.05% |
| 160 x 4 | +1.7% | -0.25% |

A negative error means the wall is too stiff. The cubic solid is within 0.25%
on every mesh, and converges to a value 0.05% from the beam, which is the
continuum-beam difference. The linear solid stiffens as the cells become
elongated, and needs eight times the in-plane resolution to match the cubic
solid on the tutorial mesh; in the coupled problem its too-stiff wall
collapses 14% too little. The `static` study in `verification/`
reproduces this table.

The standard finite-volume solid is not practical on this wall. With its
segregated solution algorithm it needs about 34,000 outer iterations per solid
solve on the tutorial mesh. With PETSc SNES, preconditioned by its
compact-stencil Jacobian, it needs thousands of Krylov iterations per Newton
step, fails on refined meshes, and is 9% too stiff on the tutorial mesh.

The cubic solid does not converge through sixteen or more cells across the
wall thickness, where its linear solves diverge; the verification therefore
keeps eight cells across the thickness.

## Fluid-solid coupling

The coupling is partitioned (Dirichlet-Neumann) IQN-ILS
(`constant/fsiProperties`), with the fluid interface velocity from the mesh
motion (`newMovingWallVelocity`). Because the wall is massless, the fluid's
added mass is the only inertia in the system and the coupling is as strong as
it can be:

- `predictor yes` extrapolates the interface at the start of each time step;
- `predictSolid no`: solving the massless wall first with the old fluid load
  would jump it to the new external pressure, and passing that jump to the
  fluid unrelaxed gives an added-mass pressure about thirty times $$p_{ext}$$;
- `relaxationFactor 0.005` for the first two iterations of each time step,
  roughly the wall stiffness over the added mass divided by $$\Delta t^2$$;
  it must shrink with $$\Delta t^2$$ when the time step is refined;
- `relMinSignificant 1e-2` drops secant modes that are small relative to the
  newest one, which otherwise let round-off in the sub-solvers blow up the
  least-squares update;
- `couplingReuse 0`: re-using the secant modes of previous time steps halves
  the iteration count, but with the massless wall the old modes are nearly
  parallel to the new ones, and the coupling then diverged or stalled on some
  meshes and on one of the two platforms tested, whatever the
  `relMinSignificant` filter;
- `outerCorrTolerance 1e-4`, relative to the largest interface displacement,
  i.e. about $$2\times10^{-5}\,\mathrm{m}$$; much tighter tolerances reach
  the round-off floor of the high-order solid residual.

On the tutorial mesh the coupling converges in about 16 iterations per time
step (at most 54) with the cubic solid and 11 (at most 19) with the linear
one.

The added-mass coupling of the massless wall grows as $$1/\Delta t^2$$, so the
coupling becomes harder as the time step is refined; see the time-step study
in `verification/`.

## Running the case

```bash
./Allrun         # high-order (cubic) solid, about 11 minutes
./Allrun linear  # second-order (linear) solid, about 4 minutes
```

The case runs in serial: in parallel, the Krylov solves of the solid fail to
converge, even with an overlapping additive-Schwarz preconditioner.

## Expected results

Figure 1 compares the wall-midpoint displacement with the oomph-lib solution
of the same problem, both at the tutorial time step and converged in space
and time (see `verification/`). With the high-order solid the first trough,
$$-0.2129\,\mathrm{m}$$ at $$t = 0.47\,\mathrm{s}$$, is within 1.4% of the
converged $$-0.2159\,\mathrm{m}$$; the later oscillations are slightly
over-damped and lag by a few hundredths of a second, as the oomph-lib solution
at the same time step does. The linear solid, too stiff on this mesh,
collapses 14% too little and oscillates too fast.

![Wall-midpoint displacement](images/collapsibleChannel-wallMid.png)

### Figure 1: Vertical displacement of the wall midpoint, $$x = 10\,\mathrm{m}$$

## Verification

The `verification/` directory holds an opt-in study that compares the wall
displacement history with a converged oomph-lib solution under mesh and
time-step refinement, and compares the two solid discretisations; see
[`verification/README.md`](verification/README.md).

## Regression test

`regressionTest.sh` runs the tutorial to $$t = 1\,\mathrm{s}$$ and checks the
wall-midpoint displacement at the first trough against stored values.

## References

- O. E. Jensen and M. Heil, High-frequency self-excited oscillations in a
  collapsible-channel flow, *Journal of Fluid Mechanics*, 481, 235-268, 2003.
- M. Heil, An efficient solver for the fully coupled solution of
  large-displacement fluid-structure interaction problems, *Computer Methods
  in Applied Mechanics and Engineering*, 193, 1-23, 2004.
- M. Heil, A. L. Hazel and J. Boyle, Solvers for large-displacement
  fluid-structure interaction problems: segregated versus monolithic
  approaches, *Computational Mechanics*, 43, 91-101, 2008.
- oomph-lib, [Flow in a 2D collapsible channel](https://oomph-lib.github.io/oomph-lib/doc/interaction/fsi_collapsible_channel/html/index.html).
