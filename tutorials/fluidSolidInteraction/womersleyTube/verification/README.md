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
`backward` d2dt2 scheme applies the backward ddt twice), so there is no
start-up ramp. A small start-up transient remains, of about 0.5% in the wall
displacement over the first period, and it grows as the time-step is reduced.
With Robin-Neumann coupling it grows as the time-step is reduced because the
boundary values of the old-time displacements are zero rather than exact, and
the Robin condition imposes the solid's boundary acceleration on the fluid;
with exact boundary values, or with IQN-ILS, it does not grow
(`womersley_temporal_investigation.md`, section 6.2). The transient decays
slowly, so the driver runs six periods and analyses the last two;
`periodicity` is the change in the wall displacement coefficient from the
fifth to the sixth period.

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
study uses 50, 100 and 200 steps per period on the tutorial mesh. Both use the
Robin-Neumann coupling, which is about twice as fast as IQN-ILS here. The
coupling comparison runs both couplings on the tutorial mesh at 100 steps per
period. The runs take 2.5 to 25 minutes each in serial; the whole study takes
about 30 minutes with `--cores 5`.

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

OpenFOAM-v2412 on macOS (arm64). Errors are relative, and phases in radians.

| Case | profile | flow_amp | flow_phase | wallMid_amp |
|---|---:|---:|---:|---:|
| m1, n200 | 1.78e-2 | 5.62e-3 | 6.24e-3 | -9.73e-3 |
| m2, n200 | 4.91e-3 | 4.40e-4 | 1.31e-3 | -1.09e-3 |
| m4, n200 | 1.18e-3 | -7.15e-4 | 1.53e-5 | 8.57e-4 |
| m2, n50 | 7.05e-3 | -2.06e-3 | 1.23e-4 | 7.14e-3 |
| m2, n100 | 4.90e-3 | -1.69e-4 | 1.05e-3 | 8.40e-4 |
| m2, n100, IQN-ILS | 4.91e-3 | -1.41e-4 | 1.06e-3 | 1.12e-3 |

| Case | wallMid_phase | speed | attenuation | Iter./step |
|---|---:|---:|---:|---:|
| m1, n200 | -3.71e-3 | -1.80e-3 | -8.46e-3 | 12.3 |
| m2, n200 | -2.50e-3 | 4.83e-4 | -6.04e-3 | 10.3 |
| m4, n200 | -1.88e-3 | 8.93e-4 | -4.47e-3 | 7.9 |
| m2, n50 | -7.87e-3 | 2.97e-3 | -1.99e-2 | 10.4 |
| m2, n100 | -3.74e-3 | 1.21e-3 | -9.63e-3 | 8.7 |
| m2, n100, IQN-ILS | -3.64e-3 | 1.27e-3 | -9.46e-3 | 14.8 |

`m` is the mesh factor and `n` the number of time-steps per period; the runs
use the Robin-Neumann coupling unless stated. The change in the wall
displacement coefficient from the fifth to the sixth period is below
$$10^{-4}$$ in every run.

### Observed orders

| Quantity | Mesh | Time step |
|---|---:|---:|
| profile | 2.06 | – |
| flow_amp | 2.16 | 1.64 |
| flow_phase | 1.92 | 1.84 |
| wallMid_amp | 2.15 | 1.71 |
| wallMid_phase | 0.97 | 1.74 |
| speed | 2.48 | 1.29 |
| attenuation | 0.63 | 1.51 |

The profile order is from its errors on the two finest meshes. The other
orders are from the differences between the three runs of a series, which
cancel the error the series shares, such as the time error of the mesh study
(at 200 steps per period it is of the same size as the mesh error on the
finest mesh).

- The velocity profile, the flow rate, the wall displacement amplitude and the
  wave speed converge at second order in space.
- The wall displacement phase and the attenuation converge more slowly, at
  about first order, and the attenuation error is the largest, -0.45% on the
  finest mesh. $$\mathrm{Im}(k)$$ is only 14% of $$\mathrm{Re}(k)$$, so a given
  error in the fitted wave number is seven times larger relative to
  $$\mathrm{Im}(k)$$.
- In time, the orders are 1.3 to 1.8, below the nominal second order of the
  `backward` scheme, which is used in both regions and in the interface
  velocity condition; it is caused by an O(dt) mass-flux
  inconsistency at the tube ends, where the exact pressure and the exact
  normal velocity gradient are imposed (OpenFOAM treats the mixed velocity
  condition as fixing the value), not by the coupling; with 400 steps per
  period the orders fall to about one. The `pimpleFluid` option
  `fluxConsistentPatches (inlet outlet);` removes the O(dt) term (orders
  1.87-2.08 over the steps used); an O(dt h) boundary term remains, so at
  fixed mesh the order tends to one only at much smaller steps. See
  `womersley_temporal_investigation.md`.
- IQN-ILS and Robin-Neumann agree to $$3 \times 10^{-4}$$ or better on every
  quantity. Robin-Neumann takes 8.7 coupling iterations per time-step against
  14.8 for IQN-ILS, and half the run time.

### Tolerances

| Check | Tolerance |
|---|---|
| Finest mesh (m4, n200): profile, flow rate, wall amplitude, wave speed | 3e-3 |
| Finest mesh: wall phase | 5e-3 rad |
| Finest mesh: attenuation | 1e-2 |
| Mesh order of the profile, flow rate, wall amplitude and wave speed | 1.5 |
| Time order of the signed quantities | 1.2 |
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
The order thresholds sit below the observed orders. The wall phase and the
attenuation are not checked for their mesh order, which is about one; every
signed quantity is checked for its time order.

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
