# ringAddedMass verification study

This opt-in study verifies the partitioned fluid-solid coupling against an
exact solution: the free n = 2 (ovalling) vibration of an elastic ring in a
fluid-filled rigid annulus. The fluid adds inertia to the ring and lowers its
natural frequency by an amount that linear theory gives exactly, and the
strength of this added-mass effect is set by the fluid density. Unlike the
published numerical benchmarks, the reference has no discretisation error of
its own.

The driver runs dry (in vacuo) and wet (fluid-loaded) versions of the tutorial
under mesh and time-step refinement for three added-mass levels, measures the
frequencies, and compares the dry and wet frequencies and their ratio with the
exact values. It also compares the standard and high-order solid residuals and
the IQN-ILS and Robin-Neumann couplings.

It is deliberately separate from `regressionTest.sh`, and nothing here is run by
`tutorials/Alltest` or `tutorials/Alltest-regression`.

## Exact solution

### Thin ring

A thin circular ring of mean radius $$a$$, thickness $$h$$, density
$$\rho_s$$ and plane-strain bending rigidity $$D = E h^3 / (12 (1 - \nu^2))$$
has inextensional modes $$w = A \cos n\theta$$, $$v = -(A/n) \sin n\theta$$
(radial and tangential displacement). The strain energy of bending,
$$\tfrac{1}{2} \pi D (n^2 - 1)^2 A^2 / a^3$$ per unit length, and the kinetic
energy, $$\tfrac{1}{2} \pi \rho_s h a (1 + 1/n^2) \dot{A}^2$$, which includes
the tangential motion, give the classical

$$
\omega_{dry}^2 = \frac{D}{\rho_s h a^4} \frac{n^2 (n^2 - 1)^2}{n^2 + 1}.
$$

Rotary inertia and shear deformation are of relative order $$(h/a)^2$$.

### Added mass

For inviscid incompressible flow in the annulus between the ring's outer
surface, radius $$R$$, and a rigid wall at radius $$b$$, the potential
$$\phi = \dot{A} (C r^n + C' r^{-n}) \cos n\theta$$ with
$$\partial \phi / \partial r = \dot{w}$$ at $$r = R$$ and
$$\partial \phi / \partial r = 0$$ at $$r = b$$ gives the interface pressure
$$p = -\rho_f \partial \phi / \partial t = m_a \ddot{w}$$, with the added mass
per unit area of ring surface

$$
m_a = \frac{\rho_f R}{n} \frac{b^{2n} + R^{2n}}{b^{2n} - R^{2n}}.
$$

Only the radial motion couples: an inviscid fluid exerts no shear on the
ring's tangential motion. For $$n = 1$$, a rigid cylinder translating inside a
concentric cylinder, the added mass per unit length is
$$\pi R m_a = \rho_f \pi R^2 (b^2 + R^2) / (b^2 - R^2)$$, the entry for a
cylinder coaxial with an enclosing cylinder in Table II of Brennen (1982). As
$$b \to \infty$$, $$m_a$$ tends to the unconfined value $$\rho_f R / n$$.

The fluid kinetic energy, $$\tfrac{1}{2} \pi m_a R \dot{A}^2$$, is added to the
ring's, so the wet-to-dry frequency ratio of a thin ring is

$$
\frac{\omega_{wet}}{\omega_{dry}} = (1 + \mu)^{-1/2}, \qquad
\mu = \frac{n^2}{n^2 + 1} \frac{m_a R}{\rho_s h a},
$$

where $$\mu$$ is the ratio of the added to the structural modal mass. The
factor $$R/a$$ accounts for the fluid acting on the outer surface rather than
the mean radius; without it the ratio is off by $$O(h/a)$$ (2% here at
$$\mu = 10$$).

### Continuum ring and viscous fluid

solids4foam models the ring as a two-dimensional plane-strain continuum of
finite thickness, which only tends to a thin ring as $$h/a \to 0$$. The
reference values are therefore not the thin-ring ones but the exact linear
eigenvalues of the continuum: `scripts/ring_exact_frequencies.py` solves the
plane-strain Navier equations for $$u_r = U(r) \cos n\theta$$,
$$u_\theta = V(r) \sin n\theta$$ across the thickness by Chebyshev collocation,
with a traction-free inner surface and the exact added-mass traction on the
outer surface. Twenty-four collocation intervals give the frequencies to about
$$10^{-8}$$. With the fluid acting on the outer surface, the thin-ring
formula matches the continuum dry frequency to $$0.5 (h/a)^2$$ and the
frequency ratio to $$0.12 (h/a)^2$$ at $$\mu = 10$$, which checks the
collocation solution; the ratio cancels most of the thickness effect.

The same script also solves the linearised Navier-Stokes equations in the
annulus, with no slip on both walls, through the stream function
$$\Psi = A r^n + B r^{-n} + C I_n(\kappa r) + D K_n(\kappa r)$$,
$$\kappa = \sqrt{s / \nu_f}$$, where $$s$$ is the complex eigenvalue. The
viscous tractions then couple the tangential motion too, and $$s$$ gives the
viscous frequency shift and damping. This solution reproduces an independent
collocation of the fluid to $$10^{-5}$$.

### Assumptions

- Linear elasticity and small displacements: the ring vibrates with an
  amplitude of about $$5 \times 10^{-4} a$$, and the solid uses the linear
  geometry model.
- Incompressible fluid, linearised: the convective term is of relative order
  amplitude/gap, below $$10^{-3}$$.
- Inviscid fluid: the viscosity, $$\nu_f = 10^{-7}$$ m²/s, keeps the Stokes
  layer $$\sqrt{2 \nu_f / \omega}$$ below 0.12% of the gap. The exact viscous
  correction is a frequency shift of -0.004%, -0.026% and -0.073%, and a
  damping ratio of the same size, for the three levels. The meshes do not
  resolve this layer, so the computed flow is effectively inviscid; the
  inviscid values are the references.

## Parameter sets

The geometry and the solid are those of the tutorial: $$R_i = 0.95$$ m,
$$R = 1.05$$ m, $$b = 1.5$$ m, $$E = 1$$ MPa, $$\nu = 0.3$$,
$$\rho_s = 1000$$ kg/m³, so $$h/a = 0.1$$. The fluid density sets the
added-mass level.

| Level | ρf (kg/m³) | μ | ω (rad/s) | Ratio | Thin-ring ratio |
|---|---:|---:|---:|---:|---:|
| dry | – | – | 2.555042 | – | – |
| weak | 14 | 0.101 | 2.435615 | 0.953258 | 0.953136 |
| moderate | 140 | 1.008 | 1.804513 | 0.706256 | 0.705776 |
| strong | 1400 | 10.08 | 0.768654 | 0.300838 | 0.300482 |

The thin-ring dry frequency is 2.567763 rad/s, 0.5% above the continuum one.
All values are stored in `reference/ringAddedMass_verification_references.json`,
and `scripts/ring_exact_frequencies.py --check` (which needs numpy) recomputes
them.

## Running

Source an OpenFOAM environment with solids4foam built with PETSc, and run:

```bash
cd tutorials/fluidSolidInteraction/ringAddedMass/verification
./Allverify --cores 6
```

Useful options:

```bash
./Allverify --quick              # dry and moderate wet runs on the coarsest mesh
./Allverify --study main         # high-order solid mesh and time-step studies only
./Allverify --study solid        # standard against high-order solid
./Allverify --study coupling     # IQN-ILS against Robin-Neumann
./Allverify --study monolithic   # quasi-monolithic Newton against IQN-ILS
./Allverify --levels moderate    # one added-mass level
./Allverify --reuse              # reuse completed runs with matching settings
./Allverify --cores 32           # number of serial runs at the same time
```

Each run is a copy of the tutorial under the ignored `verification/work/`
directory, without any results the tutorial itself may contain; dry runs use the
solid region of the tutorial as a solid-only case. The driver writes a
`verification_settings.json` in each run; `--reuse` only reuses a run whose
settings, including the content of the tutorial's `0`, `constant` and `system`
directories and `Allrun`, the OpenFOAM version and the solids4foam executable
and library, match the requested ones, whose solver log ends with
`End`, and whose displacement histories are complete and finite. Results are
written to the ignored `verification/postProcessing/` directory as
`ringAddedMass_results.csv` and `verification_summary.md`. The driver exits
with 0 only if every check passes.

## Method

- **Runs.** Each run lasts four periods of the exact frequency, with a fixed
  number of time-steps per period. The time-step thus differs between the dry
  and wet runs, and between the levels, but the time-discretisation error,
  which depends on $$\omega \Delta t$$, is similar.
- **Meshes.** Refinement factors 1, 2 and 4 of the tutorial meshes: 4 × 48 to
  16 × 192 solid cells (radial × circumferential) and 12 × 48 to 48 × 192
  fluid cells. The interface meshes match.
- **Time-steps.** 100 steps per period for the mesh study; 25, 50 and 100
  steps per period on the refinement-factor-2 mesh for the time-step study.
- **Frequency.** The zero crossings of the ovalling displacement,
  $$(u_r(0) - u_r(90°))/2$$ on the outer surface, are found by linear
  interpolation, and the damped frequency is fitted to them; the zero
  crossings of a damped sinusoid are exactly half a damped period apart. The
  first crossing is skipped. The damping is fitted to the peaks, and reported
  as the amplitude loss per period. The damping ratios are below 0.6%, so the
  damped and undamped frequencies agree to $$2 \times 10^{-5}$$.
- **Ratio.** Each wet frequency is divided by the dry frequency on the same
  mesh, with the same solid model and number of steps per period, which
  cancels most of the structural and time-discretisation errors.

## Acceptance criteria

The tolerances are set by the measured convergence and the model
assumptions, and are stored in the reference JSON.

- **Frequency, finest mesh, 100 steps per period: 0.3%.** The second-order
  backward scheme underestimates the frequency by 0.12% at 100 steps per
  period (measured, and consistent with the observed second order), the high
  order solid is converged in space to 0.02%, and the viscous correction that
  the computation neglects is up to 0.07%.
- **Time-extrapolated frequency: 0.1%.** Richardson extrapolation of the 50 and
  100 steps-per-period frequencies with order 2 removes the time error; the
  remaining allowance covers the neglected viscous correction and the mesh
  error.
- **Mesh changes in the frequency: 0.1%.** The high-order solid is converged in
  space on the coarsest mesh, so the frequency changes between meshes are at
  the level of the remaining fluid and coupling errors, too small for a
  meaningful order; they must stay below 0.1%.
- **Wet/dry ratio, finest mesh: 0.1%**, and the ratio error must decrease
  under mesh refinement.
- **Ratio against thin-ring theory: 0.25%.** This is the exact thin-ring
  $$O((h/a)^2)$$ difference, up to 0.12% here, plus the 0.1% ratio
  tolerance.
- **Observed time order at least 1.5** (formally 2).
- **Standard solid: observed mesh order at least 1.5** and a finest-mesh ratio
  error within 0.1%.
- **Amplitude loss per period below 2%** at 100 steps per period on the
  time-study mesh (refinement factor 2), where it is at most 1.14%; the loss
  is larger on the coarsest mesh and halves with each refinement. A free
  vibration must not grow, so no run may gain more than 0.2% per period,
  which allows for the measurement noise of an undamped run.
- **Robin-Neumann and IQN-ILS frequencies agree to 0.1%.**
- **Monolithic study:** at 100 steps per period, the monolithic and IQN-ILS
  frequencies (standard solid) agree to 0.1%, the monolithic wet/dry ratio is
  within 0.1% of the exact value, and the monolithic amplitude loss is below
  2% per period.

A `--quick` run only checks that the coarse dry and moderate wet frequencies
are within 1% of the exact values.

`verification/Allclean` removes the generated `work/` and `postProcessing/`
directories, as does `./Allclean all` in the case directory; the plain case
`./Allclean` keeps them.

## Reference results

Recorded with OpenFOAM v2412 on xenosim with `./Allverify --cores 32`; all
checks pass. The 35 serial runs took about 30 minutes; the longest run, the
strong level on the finest mesh, about 25 minutes. Mesh is the refinement
factor, Steps the time-steps per period, Error the frequency error against the
exact continuum value, Ratio and Thin ring the errors of the wet/dry ratio
against the exact and the thin-ring values, Loss the amplitude loss per
period, and Iter. the mean coupling iterations per time-step.

### High-order solid, IQN-ILS

| Level | Mesh | Steps | ω (rad/s) | Error | Ratio | Thin ring | Loss |
|---|---:|---:|---:|---:|---:|---:|---:|
| dry | 1 | 100 | 2.55254 | -0.098% | – | – | 0.04% |
| dry | 2 | 100 | 2.55222 | -0.110% | – | – | 0.04% |
| dry | 4 | 100 | 2.55173 | -0.129% | – | – | 0.04% |
| dry | 2 | 25 | 2.50576 | -1.929% | – | – | 2.18% |
| dry | 2 | 50 | 2.54240 | -0.495% | – | – | 0.30% |
| weak | 1 | 100 | 2.43341 | -0.091% | +0.007% | +0.020% | 0.25% |
| weak | 2 | 100 | 2.43298 | -0.108% | +0.002% | +0.015% | 0.15% |
| weak | 4 | 100 | 2.43248 | -0.129% | +0.001% | +0.013% | 0.10% |
| weak | 2 | 25 | 2.38894 | -1.916% | +0.013% | +0.025% | 2.30% |
| weak | 2 | 50 | 2.42378 | -0.486% | +0.009% | +0.022% | 0.41% |
| moderate | 1 | 100 | 1.80344 | -0.059% | +0.038% | +0.106% | 1.18% |
| moderate | 2 | 100 | 1.80273 | -0.099% | +0.011% | +0.079% | 0.64% |
| moderate | 4 | 100 | 1.80224 | -0.126% | +0.003% | +0.071% | 0.35% |
| moderate | 2 | 25 | 1.77093 | -1.861% | +0.069% | +0.137% | 2.85% |
| moderate | 2 | 50 | 1.79645 | -0.447% | +0.048% | +0.116% | 0.92% |
| strong | 1 | 100 | 0.76843 | -0.029% | +0.068% | +0.187% | 2.11% |
| strong | 2 | 100 | 0.76796 | -0.090% | +0.021% | +0.139% | 1.14% |
| strong | 4 | 100 | 0.76770 | -0.124% | +0.006% | +0.124% | 0.64% |
| strong | 2 | 25 | 0.75477 | -1.806% | +0.125% | +0.244% | 3.41% |
| strong | 2 | 50 | 0.76552 | -0.408% | +0.087% | +0.206% | 1.43% |

- **Time step.** The observed orders of the frequency are 1.90, 1.92, 2.02
  and 2.14 (dry, weak, moderate, strong), and Richardson extrapolation brings
  every frequency to within +0.018% of the exact value. The -0.12% error at
  100 steps per period is therefore the backward scheme's phase error.
- **Mesh.** The high-order solid is converged in space on the coarsest mesh:
  the frequency changes by at most 0.06% between meshes, too little for a
  meaningful order. The wet/dry ratio error, which removes the common time
  error, converges at orders 1.77, 1.73 and 1.68, to +0.001%, +0.003% and
  +0.006% on the finest mesh.
- **Thin ring.** Against thin-ring theory the finest-mesh ratio differs by
  +0.013%, +0.071% and +0.124%, which is the exact $$O((h/a)^2)$$ thickness
  effect of the continuum.
- **Damping.** No physical damping is modelled (the viscous damping ratio is
  below 0.07%). The dry ring loses 0.04% of its amplitude per period at 100
  steps per period, the backward scheme's dissipation. The wet runs lose more,
  up to 2.1% on the coarsest mesh at the strong level, and this halves with
  each mesh refinement, so it is a spatial error of the fluid or the moving
  interface, not of the time scheme.

### Standard against high-order solid

| Level | Mesh | ω standard (rad/s) | Error | Ratio error |
|---|---:|---:|---:|---:|
| dry | 1 | 3.08290 | +20.66% | – |
| dry | 2 | 2.71258 | +6.17% | – |
| dry | 4 | 2.59528 | +1.58% | – |
| weak | 4 | 2.47400 | +1.58% | +0.001% |
| moderate | 4 | 1.83303 | +1.58% | +0.005% |
| strong | 4 | 0.78083 | +1.58% | +0.009% |

The standard second-order solid is much too stiff in bending on these meshes,
20.7% on the tutorial mesh, and converges at an observed order of 1.66; the
high-order solid is within 0.13% on every mesh. The wet/dry ratio, however,
agrees with the high-order one to within 0.003% on the finest mesh, because
the solid error cancels between the dry and the wet runs: the added-mass
coupling is verified independently of the solid's accuracy.

### IQN-ILS against Robin-Neumann

Refinement factor 2, 100 steps per period, high-order solid:

| Level | Coupling | ω (rad/s) | Error | Loss | Iter. | Time (s) |
|---|---|---:|---:|---:|---:|---:|
| weak | IQN-ILS | 2.43298 | -0.108% | 0.15% | 2.99 | 125 |
| weak | Robin | 2.43291 | -0.111% | 0.04% | 3.87 | 156 |
| moderate | IQN-ILS | 1.80273 | -0.099% | 0.64% | 3.95 | 185 |
| moderate | Robin | 1.80242 | -0.116% | 0.04% | 4.85 | 213 |
| strong | IQN-ILS | 0.76796 | -0.090% | 1.14% | 4.91 | 406 |
| strong | Robin | 0.76773 | -0.121% | 0.05% | 12.29 | 575 |

The two couplings agree to 0.03% in frequency. IQN-ILS needs fewer coupling
iterations at every level, and its advantage grows with the added mass
(4.9 against 12.3 per step at the strong level). The Robin-Neumann solution
is almost undamped, whereas the IQN-ILS one shows the mesh-dependent damping
described above.

### Quasi-monolithic Newton against IQN-ILS

`./Allverify --study monolithic` runs the quasi-monolithic Newton coupling
(`./Allrun monolithic`, OpenFOAM.com only) under time-step refinement on the
refinement-factor-2 mesh, together with the IQN-ILS coupling and the dry ring.
The monolithic velocity solid has the standard residual only, so every run of
this study uses the standard solid, and the comparison is through the wet/dry
ratio, which removes the solid error (see above). It is not part of
`--study all`. Recorded with OpenFOAM v2512 on macOS (Apple M-series), 10
serial runs at a time; the 21 runs took about 2 minutes, and the times below
are therefore only indicative:

| Level | Coupling | Steps | ω (rad/s) | Ratio | Loss | Iter. | Time (s) |
|---|---|---:|---:|---:|---:|---:|---:|
| weak | IQN-ILS | 25 | 2.53346 | +0.015% | 2.70% | 3.01 | 18 |
| weak | IQN-ILS | 50 | 2.57483 | +0.010% | 0.47% | 3.00 | 29 |
| weak | IQN-ILS | 100 | 2.58588 | +0.004% | 0.16% | 2.98 | 51 |
| weak | monolithic | 25 | 2.53306 | -0.000% | 2.57% | 2.00 | 23 |
| weak | monolithic | 50 | 2.57455 | -0.000% | 0.36% | 2.00 | 46 |
| weak | monolithic | 100 | 2.58578 | -0.000% | 0.045% | 2.00 | 86 |
| moderate | IQN-ILS | 25 | 1.87813 | +0.075% | 3.27% | 3.57 | 25 |
| moderate | IQN-ILS | 50 | 1.90853 | +0.056% | 0.98% | 3.92 | 44 |
| moderate | IQN-ILS | 100 | 1.91612 | +0.018% | 0.65% | 3.98 | 73 |
| moderate | monolithic | 25 | 1.87666 | -0.003% | 2.56% | 2.00 | 24 |
| moderate | monolithic | 50 | 1.90741 | -0.003% | 0.35% | 2.00 | 47 |
| moderate | monolithic | 100 | 1.91574 | -0.002% | 0.040% | 2.00 | 89 |
| strong | IQN-ILS | 25 | 0.80042 | +0.126% | 3.77% | 3.86 | 29 |
| strong | IQN-ILS | 50 | 0.81333 | +0.101% | 1.50% | 4.63 | 64 |
| strong | IQN-ILS | 100 | 0.81632 | +0.033% | 1.15% | 4.92 | 120 |
| strong | monolithic | 25 | 0.79936 | -0.007% | 2.57% | 2.00 | 26 |
| strong | monolithic | 50 | 0.81246 | -0.005% | 0.36% | 2.00 | 48 |
| strong | monolithic | 100 | 0.81602 | -0.004% | 0.049% | 2.00 | 91 |

Ratio is the error of the wet/dry ratio against the exact value, and Iter. the
mean coupling (IQN-ILS) or Newton (monolithic) iterations per time-step. The
monolithic coupling converges in two Newton iterations per step at every
level, as the problem is linear. Every check of the study passes.

- **Time step.** The monolithic frequency converges at order 1.88-1.89, and
  its wet/dry ratio is within 0.007% of the exact value at every level and
  time-step. Its amplitude loss, 2.57%, 0.36% and 0.045% per period at 25, 50
  and 100 steps per period at every level (order 2.8), is that of the dry
  ring, 2.57%, 0.36% and 0.046%, i.e. the dissipation of the backward scheme in
  the solid: the coupling adds none. The IQN-ILS coupling adds a loss that
  grows with the added mass, up to 1.15% per period at the strong level.
- **Mesh.** At 100 steps per period, the monolithic ratio error is -0.016%,
  -0.002% and -0.0001% (moderate) and -0.031%, -0.004% and -0.0004% (strong)
  on the refinement-factor-1, 2 and 4 meshes, an observed order of about 3,
  against +0.069%, +0.018% and +0.005% and +0.122%, +0.033% and +0.009% for
  IQN-ILS (recorded separately; not part of the study).

#### Defects found by this study

The study found two defects in the quasi-monolithic coupling, both now
corrected.

- **One-step lag of the interface force (time).** The first version lost
  2.0%, 10.6% and 18.3% of the amplitude per period at 100 steps per period
  (weak, moderate, strong), a loss that grew with the added mass, fell only at
  first order in the time-step and hardly depended on the mesh. The velocity
  solid averaged the stress at the old and the extrapolated future time, but
  the old-time half of the interface force was the traction from two
  time-steps before, so that the interface force lagged by one time-step: a
  damping force of $$m_a \omega^2 \Delta t$$ times the velocity, which gives a
  loss of $$\pi (m_a/M) \omega \Delta t$$ per period, 1.8%, 9.9% and 17.9%
  here. With the lag removed, the centred scheme dissipated nothing at any
  frequency, and a high-frequency mode grew in `blobInTreacle`. The velocity
  solid now uses the backward (BDF2) scheme for the displacement, consistent
  with `ddt(U)`, and evaluates the stress and the interface traction once, at
  the new time.
- **Zero-gradient interface pressure (space).** With the BDF2 solid, the
  amplitude loss matched IQN-ILS, but the wet/dry ratio error converged only
  at first order in the mesh, -0.194%, -0.117% and -0.063% (moderate) and
  -0.355%, -0.213% and -0.115% (strong), too large an added mass by about a
  quarter of a fluid cell. A `zeroGradient` pressure on the interface halves
  the normal pressure gradient that the least-squares gradient gives in the
  cells next to the moving wall, and makes the cell-centre pressure the
  interface pressure; the two first-order errors partly cancel. The
  monolithic variants now use a `fixedGradient` pressure on the interface,
  whose gradient the coupling sets from the normal momentum equation at the
  no-slip wall, $$\partial p/\partial n = n \cdot (\nu \nabla^2 U - a)$$, with
  $$a$$ the acceleration of the interface, and the interface velocity passed
  to the fluid is extrapolated from the solid cell centres.

## References

C. E. Brennen, A review of added mass and fluid inertial forces, Report CR
82.010, Naval Civil Engineering Laboratory, Port Hueneme, California, 1982.
Freely available from CaltechAUTHORS.
