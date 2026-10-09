# womersleyTube verification study

This opt-in study verifies pulsatile flow in a compliant tube, with the
partitioned fluid-solid coupling, against the exact linear solution for a
harmonic pressure wave travelling along an infinitely long elastic tube: the
velocity profile, the flow rate, the wall displacement and the complex wave
number, which gives the wave speed and the viscous attenuation. The reference
has no discretisation error of its own.

The driver runs copies of the tutorial under mesh and time-step refinement,
takes the Fourier coefficient at the forcing frequency of each sampled history
over the last period, and compares it with the complex amplitude of the exact
solution. It also compares the IQN-ILS and Robin-Neumann couplings.

It is deliberately separate from `regressionTest.sh`, and nothing here is run by
`tutorials/Alltest` or `tutorials/Alltest-regression`.

## Exact solution

Every field is the real part of a complex amplitude times
$$e^{i(\omega t - k x)}$$, where $$x$$ is the axial coordinate and the wave
number $$k$$ is complex: $$\omega / \mathrm{Re}(k)$$ is the wave speed and
$$\mathrm{Im}(k) < 0$$ the viscous attenuation.

### Womersley's thin-walled tube

In the long-wave limit the axial velocity in a tube of radius $$R$$ is

$$
u_x(r) = \frac{k P}{\rho \omega}
\left[1 - \frac{J_0(\Lambda r / R)}{J_0(\Lambda)}\right]
+ i \omega \xi \frac{J_0(\Lambda r / R)}{J_0(\Lambda)},
\qquad \Lambda = i^{3/2} \alpha, \qquad \alpha = R \sqrt{\omega / \nu_f},
$$

with $$P$$ the pressure amplitude and $$\xi$$ the axial wall displacement, and
the flow rate is
$$Q = \pi R^2 [k P / (\rho \omega) (1 - g) + i \omega \xi g]$$ with
$$g = F_{10} = 2 J_1(\Lambda) / (\Lambda J_0(\Lambda))$$. A thin elastic
membrane wall that is free to move axially gives Womersley's frequency
equation, in the form of Filonova et al. (2020), Eq. 41,

$$
(1 - g)(1 - \nu^2) v^2 - \left[2 + m (1 - g) + g (\tfrac{1}{2} - 2 \nu)\right] v
+ g + 2 m = 0,
\qquad v = \frac{E h}{(1 - \nu^2) \rho R c^2},
\qquad m = \frac{\rho_s h}{\rho R}.
$$

Of its two roots, the pressure (Young) wave is the one nearer the tethered
value; the other is the axial wave of the wall. For a longitudinally tethered
wall, $$\xi = 0$$, the equation reduces to
$$c = c_0 \sqrt{(1 - g) / (1 - \nu^2)}$$ for a massless wall, with the
Moens-Korteweg speed $$c_0 = \sqrt{E h / (2 \rho R)}$$.

### Exact continuum solution

solids4foam models the wall as an elastic continuum of finite thickness, and
the fluid by the full Navier-Stokes equations rather than their long-wave
limit, so the reference values are the exact linear solution of this
continuum problem, computed by `scripts/womersley_exact.py`:

- The fluid is the linearised incompressible Navier-Stokes solution,
  $$p = A I_0(k r)$$,
  $$u_x = k A I_0(k r) / (\rho \omega) + B I_0(s r)$$,
  $$u_r = i k A I_1(k r) / (\rho \omega) + i k B I_1(s r) / s$$,
  with $$s^2 = k^2 + i \omega / \nu_f$$.
- The wall, $$R_i < r < R_o$$, is a linear elastic cylinder whose radial
  displacement profiles are found by Chebyshev collocation of the Navier
  equations; sixteen intervals give $$k$$ to round-off.
- The fluid and the wall share the velocity and the traction at
  $$r = R_i$$, linearised about the undeformed position, and the outer
  surface is traction free.

The frequency equation is solved for $$k$$ by the secant method from
Womersley's value. The collocation solution is checked in two limits: as
$$h / R \to 0$$ at fixed $$E h$$ and $$\rho_s h$$, the difference between the
exact $$k$$ and Womersley's falls in proportion to $$h / R$$, for both the
free and the tethered wall, to the long-wave correction of order
$$(k R)^2$$ that remains in the limit; and for the tethered wall at long
wavelength (a stiffer wall, $$k R \approx 0.015$$) it agrees
with Womersley's equation with the plane-strain (Lamé) stiffness of the thick
cylinder to $$2 \times 10^{-5}$$. Womersley's quadratic was also checked
against an independent derivation from the membrane equations, which agrees
to round-off for a massless wall and to $$4 \times 10^{-4}$$ with the wall
mass, whose axial inertia the quadratic approximates.

### Assumptions

- Linear elasticity and small displacements: the radial wall displacement is
  $$5 \times 10^{-4} R$$, and the solid uses the linear-geometry model.
- Linearised flow: the peak velocity is $$1.5 \times 10^{-3}$$ of the wave
  speed, so the convective term is of that relative order. The fluid mesh
  moves with the wall, a geometric effect of the order of the wall
  displacement over the radius.

## The thickness effect

The tutorial wall has $$h / R = 0.1$$, and Womersley's thin-wall wave number
is 4.7% too small in its real part (the wave speed is 4.7% too high). The
table gives the ratio of the exact to Womersley's wave number as the wall is
thinned at fixed $$E h$$ and $$\rho_s h$$ (`reference/*.json`, `thickness`);
the difference falls in proportion to $$h / R$$, towards about 0.15% (a
long-wave correction of order $$(k R)^2$$) as $$h / R \to 0$$.

| h/R | Free wall | Tethered wall |
|---:|---:|---:|
| 0.2 | 1.0928 - 0.0068i | 1.0951 - 0.0003i |
| 0.1 | 1.0473 - 0.0040i | 1.0482 - 0.0005i |
| 0.05 | 1.0245 - 0.0024i | 1.0247 - 0.0006i |
| 0.02 | 1.0109 - 0.0014i | 1.0107 - 0.0006i |
| 0.01 | 1.0064 - 0.0010i | 1.0060 - 0.0006i |
| 0.005 | 1.0041 - 0.0009i | 1.0037 - 0.0006i |

The simulations are therefore compared with the exact continuum solution,
which has no thin-wall error.

## Boundary conditions and start

The exact solution is a single travelling wave of an infinite tube. The
tutorial's finite segment, $$0 < x < L$$, reproduces it exactly because the
same wave is imposed at both ends, so nothing reflects and no second wave is
generated:

- Fluid: the exact pressure (which varies slightly over the section, as
  $$I_0(k r)$$) and the exact normal gradient of the velocity,
  $$\partial \mathbf{u} / \partial x = -i k \mathbf{u}$$, at both ends.
- Solid: the exact radial and axial displacement across the wall at both
  ends. The outer surface is traction free, and the inner surface is the
  fluid-solid interface.

The fields are set to the exact solution at $$t = 0$$ (and the solid's old-time
displacements at $$-\Delta t$$, $$-2 \Delta t$$ and $$-3 \Delta t$$, as the
`backward` d2dt2 scheme applies the backward ddt twice), in the cells and on
the boundaries, so there is no start-up ramp. The boundary values of the
old-time displacements matter with Robin-Neumann coupling, which imposes the
solid's boundary acceleration on the fluid: with the placeholder zero values
an earlier version had, the first-period transient grew as the time-step was
reduced (to 3.4% of the amplitude at 400 steps per period); with exact values
it is 0.2-0.3% from 100 to 400 steps per period (1.6% at 50, a resolution
effect) (`womersley_temporal_investigation.md`, sections 6.2 and 15). The
transient decays slowly, so the driver runs six periods and analyses the last
two; `periodicity` is the change in the wall displacement coefficient from the
fifth to the sixth period.

At the tube ends the pressure is fixed and the velocity condition (mixed,
imposing the gradient) fixes the boundary value, so without further
treatment the pressure equation's end flux differs from the boundary
velocity flux by a term proportional to the time-step, which limits the
time-step study to first order. The tutorial uses the `pimpleFluid` entry
`fluxConsistentPatches (inlet outlet);`, which makes the end flux consistent
with the boundary velocity to O(Δt h): the observed time order is then two
over the study's steps and under joint refinement, while at a fixed mesh it
would tend to one at much smaller steps (beyond about 6400 steps per period on
the tutorial mesh; `womersley_temporal_investigation.md`, section 15).

### Why the wall is not tethered

Womersley's classical benchmark, and Filonova et al. (2020), use a
longitudinally tethered wall. A tethered outer surface needs a condition that
fixes the axial displacement and leaves the radial traction free. The only
such solids4foam condition, `displacementOrTraction` (with
`specifyNormalDirection -1`), is not enforced by the residual of the PETSc
SNES solid solver, which only enforces `solidTraction` boundaries directly: in
a test its SNES iterations diverged from the first time-step. The segregated
solid solver diverged with it too; see
[#511](https://github.com/solids4foam/solids4foam/issues/511). The wall is
therefore free, and the free
wall of the same parameters is solved exactly instead; the tethered exact and
thin-wall values are kept in the reference file for comparison (the tethered
wave is 0.9% faster and 13% more strongly attenuated).

## Parameters

| Quantity | Value |
|---|---|
| Inner radius $$R_i$$, thickness $$h$$, length $$L$$ | 1 m, 0.1 m, 15 m |
| Wall: $$E$$, $$\nu$$, $$\rho_s$$ | 20 kPa, 0.3, 1000 kg/m³ |
| Fluid: $$\rho$$, $$\mu$$ | 1000 kg/m³, 5 Pa s |
| Frequency, period | 0.02 Hz, 50 s |
| Pressure amplitude on the axis at $$x = 0$$ | 1 Pa |
| Womersley number $$\alpha$$ | 5.01 |
| Mass ratio $$\rho_s h / (\rho R)$$ | 0.1 |
| Moens-Korteweg speed $$c_0$$ | 1.000 m/s |
| Exact wave speed, wavelength | 0.8728 m/s, 43.64 m (43.6 R) |
| Amplitude ratio per wavelength | 0.405 |
| Phase lag and amplitude ratio over $$L$$ | 124°, 0.733 |
| Wall displacement amplitude at $$x = 0$$, radial and axial | 0.53 mm, 0.73 mm |

The problem is scaled to metre dimensions: the same dimensionless groups
($$\alpha$$, $$h / R$$, $$\lambda / R$$, the density and mass ratios and
$$\nu$$) describe, for example, a 1 cm tube at 2 Hz with $$\mu = 0.05$$ Pa s.
At the centimetre scale the displacements (about 5 μm) fall to the level of
the absolute tolerances of the PETSc solid solver and of the IQN-ILS
significance test, which then stall.

## Running

```bash
./Allverify                  # full study
./Allverify --quick          # coarsest mesh, 50 steps per period
./Allverify --study main     # mesh and time-step studies only
./Allverify --study coupling # IQN-ILS against Robin-Neumann only
./Allverify --cores 5        # run up to 5 serial cases at once
./Allverify --reuse          # reuse completed runs whose settings match
```

The driver works on copies under `verification/work/` and writes
`verification/postProcessing/verification_summary.md` and
`womersleyTube_results.csv`; both directories are ignored by git. A run is
reused only if its stored settings, including the tutorial file hashes and the
OpenFOAM version and solids4foam executable and library fingerprint, match,
and its solver log ends with `End`. Histories are rejected if they are
truncated, non-finite or non-increasing in time. The exit status is zero only
when every check passes.

`scripts/womersley_exact.py` (which needs numpy) writes the reference values,
`--json`, checks the stored ones, `--check`, and writes the exact solution
coefficients read by the tutorial, `--write-case <tutorial>`.

## Study

The mesh study refines every direction together by factors of 1, 2 and 4 from
32 axial cells, 10 radial fluid cells and 4 cells through the wall, at 200
time-steps per period; the tutorial mesh is the middle one. The time-step
study uses 50, 100, 200 and 400 steps per period on the tutorial mesh. Both use
the Robin-Neumann coupling, which is about twice as fast as IQN-ILS here. The
coupling comparison runs both couplings on the tutorial mesh at 100 steps per
period. The seven runs take 4 to 37 minutes each in serial; the whole study
takes about 40 minutes with `--cores 7`.

### Measured quantities

For each quantity, the first Fourier coefficient at the forcing frequency over
the last two periods (the fifth and sixth) is compared with the exact complex
amplitude:

- `profile`: the axial velocity of the cells just downstream of $$x = L/2$$,
  compared with the exact profile at their centres; the largest error over the
  radius, relative to the largest velocity. This covers every phase: the
  errors at phases 0, 90°, 180° and 270° are also in the CSV file.
- `flow_amp`, `flow_phase`: the flow rate through the section $$x = L/2$$.
- `wallMid_amp`, `wallMid_phase`: the radial displacement of the inner wall at
  $$x = L/2$$ (the CSV file also has $$L/4$$, $$3L/4$$ and the axial
  displacement).
- `speed`, `attenuation`: the wave speed and $$\mathrm{Im}(k)$$ from a
  least-squares fit of $$\log \hat{p}(x)$$ along the tube at $$r = R/2$$.

## Results

OpenFOAM-v2512 on Linux (x86_64), solids4foam with `fluxConsistentPatches`.
Errors are relative, and phases in radians.

| Case | profile | flow_amp | flow_phase | wallMid_amp |
|---|---:|---:|---:|---:|
| m1, n200 | 1.79e-2 | 6.10e-3 | 6.05e-3 | -1.03e-2 |
| m2, n200 | 4.96e-3 | 9.39e-4 | 1.33e-3 | -2.19e-3 |
| m4, n200 | 1.20e-3 | -3.09e-4 | 8.23e-5 | -1.76e-4 |
| m2, n50 | 6.18e-3 | -1.23e-3 | 1.43e-4 | 5.28e-3 |
| m2, n100 | 4.98e-3 | 5.14e-4 | 1.07e-3 | -6.63e-4 |
| m2, n400 | 4.95e-3 | 1.04e-3 | 1.40e-3 | -2.58e-3 |
| m2, n100, IQN-ILS | 4.98e-3 | 5.48e-4 | 1.09e-3 | -4.47e-4 |

| Case | wallMid_phase | speed | attenuation | Iter./step |
|---|---:|---:|---:|---:|
| m1, n200 | -2.37e-3 | -2.33e-3 | -2.08e-3 | 12.8 |
| m2, n200 | -1.48e-3 | -4.50e-4 | -2.18e-3 | 11.0 |
| m4, n200 | -1.15e-3 | 1.40e-5 | -1.96e-3 | 8.4 |
| m2, n50 | -6.12e-3 | 1.41e-3 | -1.31e-2 | 10.3 |
| m2, n100 | -2.34e-3 | -6.62e-5 | -4.24e-3 | 9.3 |
| m2, n400 | -1.28e-3 | -5.48e-4 | -1.69e-3 | 9.5 |
| m2, n100, IQN-ILS | -2.40e-3 | -3.72e-5 | -4.41e-3 | 14.6 |

`m` is the mesh factor and `n` the number of time-steps per period; the runs
use the Robin-Neumann coupling unless stated. The change in the wall
displacement coefficient from the fifth to the sixth period is below
$$7 \times 10^{-5}$$ in every run.

### Observed orders

| Quantity | Mesh | Time step (50/100/200, 100/200/400) |
|---|---:|---:|
| profile | 2.04 | – |
| flow_amp | 2.05 | 2.04, 2.05 |
| flow_phase | 1.92 | 1.86, 1.91 |
| wallMid_amp | 2.01 | 1.96, 1.97 |
| wallMid_phase | 1.40 | 2.14, 2.06 |
| speed | 2.02 | 1.94, 1.97 |
| attenuation | – | 2.10, 2.07 |

The profile order is from its errors on the two finest meshes. The other
orders are from the differences between three successive runs of a series,
which cancel the error the series shares, such as the time error of the mesh
study.

- The velocity profile, the flow rate, the wall displacement amplitude and the
  wave speed converge at second order in space, the wave speed to
  $$1.4 \times 10^{-5}$$ on the finest mesh.
- The wall displacement phase converges more slowly in space (order 1.4). The
  attenuation error is about $$-2 \times 10^{-3}$$ on all three meshes and
  no longer changes with the mesh: it is at the level of the linear theory
  (below), amplified because $$\mathrm{Im}(k)$$ is only 14% of
  $$\mathrm{Re}(k)$$.
- In time, every quantity converges at about second order over the study's
  steps. Without `fluxConsistentPatches` the orders were 1.3-1.8 from
  50/100/200 and fell to about one from 100/200/400, because of an O(Δt)
  inconsistency of the end flux. The remaining boundary term is O(Δt h), so
  the order stays two under joint refinement but, at a fixed mesh, would tend
  to one at much smaller steps. See `womersley_temporal_investigation.md`.
- IQN-ILS and Robin-Neumann agree to $$2.2 \times 10^{-4}$$ or better on every
  quantity. Robin-Neumann takes 9.3 coupling iterations per time-step against
  14.6 for IQN-ILS, and half the run time.

### Tolerances

| Check | Tolerance |
|---|---|
| Finest mesh (m4, n200): profile, flow rate, wall amplitude, wave speed | 3e-3 |
| Finest mesh: wall phase | 5e-3 rad |
| Finest mesh: attenuation | 1e-2 |
| Mesh order of the profile, flow rate, wall amplitude and wave speed | 1.5 |
| Time order of the signed quantities (finest three steps) | 1.7 |
| Coupling runs (m2, n100): profile | 6e-3 |
| Coupling runs: wall amplitude and wave speed | 3e-3 |
| Robin-Neumann against IQN-ILS, every quantity | 1e-3 |
| Change from the fifth to the sixth period | 5e-4 |
| `--quick` (m1, n50): profile, wall amplitude, wave speed | 3e-2, 2e-2, 1e-2 |

The linear theory neglects terms of the relative size of the wall displacement
over the radius ($$5 \times 10^{-4}$$) and of the velocity over the wave speed
($$1.5 \times 10^{-3}$$), so agreement much better than $$10^{-3}$$ cannot be
expected; the finest-mesh tolerances allow this and the measured error, with a
margin, and the looser attenuation tolerance reflects its larger sensitivity.
The order thresholds sit below the observed orders; the time-order threshold
(1.7, on the 100/200/400 triplet) fails without `fluxConsistentPatches`,
where those orders are 0.75-1.79. The wall phase and the attenuation are not
checked for their mesh order; every signed quantity is checked for its time
order.

The driver also rejects histories that are truncated, non-finite, not
increasing or not uniformly spaced in time, and fluid samples with an
unexpected number of points, and a run whose solver log lacks `End` or
reports a fatal error, whether freshly run or reused.

## Solid model

The standard (second-order) solid residual is used. The cubic high-order
residual was also tried: its moving-least-squares reconstruction stops with
"Empty direction should be vector::Z", as it requires an empty third
direction and OpenFOAM treats the wedge direction as a solution direction
([#512](https://github.com/solids4foam/solids4foam/issues/512)).
The standard solid is accurate here: the wall is loaded in hoop tension, not
bending, and eight cells through the thickness on the tutorial mesh give a
wall displacement within 0.1% of the exact value. With two cells through the
thickness (an earlier coarsest mesh) the coupled solution diverged.

On the coarsest meshes the SNES residual can stall at a round-off level, about
$$10^{-6}$$ of the cell forces, before the relative tolerance is met, so the
solid uses `stopOnPetscError false` and a relative step tolerance of
$$10^{-6}$$; the coupling tolerance then decides convergence.

## References

- J. R. Womersley, Oscillatory motion of a viscous liquid in a thin-walled
  elastic tube, I: The linear approximation for long waves, Philosophical
  Magazine 46 (1955) 199-221.
- J. R. Womersley, Oscillatory flow in arteries: the constrained elastic tube
  as a model of arterial flow and pulse transmission, Physics in Medicine and
  Biology 2 (1957) 178-187.
- V. Filonova, C. J. Arthurs, I. E. Vignon-Clementel, C. A. Figueroa,
  Verification of the coupled-momentum method with Womersley's Deformable Wall
  analytical solution, International Journal for Numerical Methods in
  Biomedical Engineering 36 (2020) e3266.
