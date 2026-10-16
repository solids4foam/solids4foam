# springSupportedBar

## Aims

- Verify the `solidSpringDashpot` boundary condition against exact
  solutions.
- Show the condition with linear geometry and with a total Lagrangian
  formulation at large strain.
- Verify the dashpot term with a spring-dashpot (Kelvin-Voigt) support.

## Case overview

A bar of length L = 10 mm along x has a spring support on its `springEnd`
patch (x = 0) and is loaded on its `loadedEnd` patch (x = L). The material has
E = 100 kPa and ν = 0.3.

### linearGeometry (default)

- `linearGeometryTotalDisplacement`, `linearElastic`.
- Cross-section 0.5 mm × 0.5 mm, 20 × 1 × 1 cells. The `sides` are symmetry
  planes, so the bar is in uniaxial strain, with constrained modulus
  M = E(1 − ν)/((1 + ν)(1 − 2ν)) = 134.6 kPa.
- `springEnd`: normal spring, k = 1e7 Pa/m.
- `loadedEnd`: traction t = 1 kPa along x.

The spring carries the full traction, so the exact displacements are:

$$
D_x(0) = \frac{t}{k} = 1.0 \times 10^{-4}\ \mathrm{m},\qquad
D_x(L) = \frac{t}{k} + \frac{t\,L}{M} = 1.742857 \times 10^{-4}\ \mathrm{m}
$$

### totalLagrangian

- `nonLinearGeometryTotalLagrangianTotalDisplacement`, `neoHookeanElastic`.
- Cross-section 3 mm × 3 mm (A₀ = 9 mm²), 20 × 3 × 3 cells, with
  traction-free `sides`. The bar stretches by about 20 % and contracts
  laterally, so the spring face area changes by about 13 %.
- `springEnd`: kNormal = 1e8 Pa/m, kTangential = 1e3 Pa/m. The coefficients
  act per unit reference area.
- `loadedEnd`: `solidForce`, 0.02 N per face along x, F = 0.18 N in total.

The force balance on the spring face gives the exact mean displacement:

$$
\bar{D}_x(0) = \frac{F}{k_N A_0} = 2.0 \times 10^{-4}\ \mathrm{m}
$$

### dashpot

- The linearGeometry case with a normal dashpot, c = 1e7 Pa·s/m, in parallel
  with the spring, so the support relaxes with time constant c/k = 1 s.
- The load is applied at time 0 and held for 10 time steps of Δt = 0.1 s.
  Inertia is neglected (`d2dt2` is `steadyState`), so the bar transmits the
  traction t to the support at every step, and the time step is still a
  physical time for the dashpot. The condition prints its warning about a
  dashpot in a `steadyState` run, which is expected here.

The dashpot velocity is (D − D.oldTime())/Δt, so the support follows the
backward-Euler recurrence $$k D_n + c (D_n - D_{n-1})/\Delta t = t$$ with
$$D_0 = 0$$, whose exact solution after n steps is

$$
D_x(0) = \frac{t}{k}\left[1 - r^n\right],\qquad
r = \frac{c/\Delta t}{k + c/\Delta t} = \frac{10}{11}
$$

giving $$D_x(0) = 6.144567 \times 10^{-5}\ \mathrm{m}$$ and
$$D_x(L) = D_x(0) + t L/M = 1.357314 \times 10^{-4}\ \mathrm{m}$$ at the end time
of 1 s.

## Running

```bash
./Allrun                    # linearGeometry
./Allrun totalLagrangian
./Allrun dashpot
```

`Allrun` copies the files of the chosen `caseOptions/` entry into the case,
and `Allclean` restores the `linearGeometry` files. Each case runs in about a
second. The patch displacements are written to
`postProcessing/0/solidDisplacementsspringEnd.dat` and
`postProcessing/0/solidDisplacementsloadedEnd.dat`; the average x displacement
is column 8.

## Regression test

`regressionTest.sh` runs the three cases and compares the mean x displacement
of `springEnd` (all cases) and `loadedEnd` (linearGeometry and dashpot) with
the exact values above, with a relative tolerance of 1e-6.
