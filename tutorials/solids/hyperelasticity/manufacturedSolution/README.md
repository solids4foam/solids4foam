---
sort: 17
---

# Method of Manufactured Solutions at Finite Strain: `manufacturedSolution`

Prepared by Ivan Batistić, with the manufactured solution derived by Pablo
Castrillo

---

## Tutorial Aims

- Demonstrate verification of the total Lagrangian solid model with the
  method of manufactured solutions (MMS) at finite strain, using the
  compressible neo-Hookean law.
- Measure the spatial order of accuracy of the displacement and stress with a
  quasi-static load ramp (`steady` mode).
- Measure the temporal order of accuracy of the `Euler`, `backward` (BDF2),
  and `NewmarkBeta` second time derivative schemes on a fixed mesh
  (`dynamic` mode).
- Compare the segregated, Jacobian-free Newton-Krylov (PETSc SNES), and
  high-order (cubic `movingLeastSquares` and `kExactLeastSquares`
  reconstruction) solution approaches.

---

## Case Overview

The domain is the unit cube $[0, 1]^3$ with Young's modulus 1 MPa, Poisson's
ratio 0.3, and density 1000 kg/m$^3$. The manufactured displacement is

$$
\boldsymbol{u}(\boldsymbol{X}, t) = \boldsymbol{a}\, S(\boldsymbol{X})\, T(t),
\qquad
S(\boldsymbol{X}) = \sin(\pi X)\sin(\pi Y)\sin(\pi Z),
$$

so that $\boldsymbol{u} = \boldsymbol{0}$ on the whole boundary. The
deformation gradient is $\boldsymbol{F} = \boldsymbol{I} + \boldsymbol{a}
\otimes \nabla_0 q$ with $q = S\,T$, and its Jacobian
$J = 1 + \boldsymbol{a} \cdot \nabla_0 q$ departs from unity by up to
approximately 0.34 at full load, so the strains are finite. The body force
per unit reference volume that balances the momentum equation,

$$
\boldsymbol{f}_0 = \rho_0 \boldsymbol{a}\, \ddot{q}
- \nabla_0 \cdot \boldsymbol{P}(\boldsymbol{u}),
$$

is evaluated in closed form from the rank-one structure of $\nabla_0
\boldsymbol{u}$, for the neo-Hookean strain energy with isochoric-volumetric
split used by the `neoHookeanElastic` law,

$$
\Psi = \frac{\mu}{2}\left(J^{-2/3} I_1 - 3\right)
+ \frac{\kappa}{4}\left(J^2 - 1 - 2\ln J\right).
$$

The derivation is given in
[`docs/neoHookean_isoVol_dynamic_MMS_2D_3D.tex`](docs/neoHookean_isoVol_dynamic_MMS_2D_3D.tex).
The tutorial-local library in `src/` implements the solution:

- `neoHookeanManufacturedSolution`: the exact displacement, Cauchy stress,
  and body force;
- `neoHookeanManufacturedSolutionSource`: an `fvOptions` source that adds
  the body force to the momentum equation on OpenFOAM.com;
- `neoHookeanManufacturedSolutionSolid`: a total Lagrangian solid model that
  adds the same source on OpenFOAM.org and foam-extend, where `fvOptions` is
  unavailable;
- `neoHookeanManufacturedSolution` boundary condition: the exact displacement
  at the patch faces;
- `neoHookeanManufacturedSolution` function object: the error norms of the
  cell displacement, point displacement, and Cauchy stress at every time
  step, and the analytical and difference fields at write times.

All parameters are set once in
`constant/neoHookeanManufacturedSolutionProperties`, which
`constant/mechanicalProperties` includes so that the solid model and the
manufactured body force use the same material.

The case runs in one of two modes, selected by `Allrun`, which differ only in
`controlDict`, `fvSchemes`, and the parameters file, each kept with a
`.steady` or `.dynamic` suffix:

- `steady`: $T = \min(t, 1)$, a quasi-static load ramp in four increments
  with the `steadyState` second time derivative scheme and
  $\boldsymbol{a} = (0.08, 0.06, 0.04)$ m. The errors at $t = 1$ measure the
  spatial discretisation error.
- `dynamic`: $T = 1 - \cos(\omega_t t)$ with $\omega_t = 20\pi$ rad/s and
  $\boldsymbol{a} = (0.04, 0.03, 0.02)$ m, so that the solid starts from
  rest, consistent with the solver, and the displacement reaches its peak
  $2\boldsymbol{a}S$ at the end time $t = 0.05$ s. With the mesh fixed, the
  errors measure the temporal discretisation error. The alternative
  `timeFunction cosineSquared`, $T = \tfrac{1}{2}(1 - \cos(\omega_t t))^2$,
  also has zero initial acceleration, which matters for the Newmark scheme
  (see below).

---

## Running the Case

The default run is the steady mode on a $10 \times 10 \times 10$ hexahedral
mesh with the segregated solution procedure, which takes seconds:

```bash
./Allrun
```

The run script accepts the mode, the solution approach, the mesh type, and,
for the dynamic mode, the time scheme:

```bash
./Allrun steady segregated hex
./Allrun steady petscSnes distHex
./Allrun steady petscSnes structTet
./Allrun steady segregated unstructTet
./Allrun steady highOrder-movingLeastSquares
./Allrun steady highOrder-kExactLeastSquares
./Allrun dynamic
./Allrun dynamic segregated hex Euler
./Allrun dynamic petscSnes hex NewmarkBeta
./Allrun dynamic highOrder-kExactLeastSquares hex backward
```

The `petscSnes` and high-order approaches require a PETSc-enabled solids4foam
build. The high-order approaches use a cubic displacement reconstruction with
face quadrature; the manufactured body force is then integrated over each
cell with the cell quadrature, the boundary condition evaluates the exact
displacement at the face quadrature points, and, because `kExactLeastSquares`
stores cell averages, its displacement error is measured against the cell
average of the exact solution. The `distHex` mesh perturbs the hexahedral
mesh with `perturbMeshPoints`, and the `structTet` and `unstructTet`
tetrahedral meshes require Gmsh; both read their cell size from
`gmsh/meshSpacing.geo`. The tutorial supports OpenFOAM.com, OpenFOAM.org, and
foam-extend.

The error norms are printed in the solver log after each time step, for
example, for the default steady run at the end of the ramp:

```text
Writing DDifference field
    Displacement error norms: mean L1, mean L2, LInf:
    Magnitude: 0.00118961959573 0.00133089091585 0.00356930685812
```

---

## Expected Results

The tables below come from `verification/Allverify` on OpenFOAM.com v2412.
The segregated and PETSc SNES approaches give the same errors to at least
five significant digits, so only the segregated results are shown for the
second-order discretisation.

### Spatial convergence (`steady` mode)

Mean L2 and L-infinity errors with respect to the exact solution at the end of
the load ramp, and the net order of accuracy between the coarsest and the
finest meshes. The structured tetrahedral meshes (`structTet`) have six cells
per hexahedron and the given numbers of cells per side:

| Mesh | Cells/side | D L2 [m] | D Linf [m] | sigma L2 [Pa] | sigma Linf [Pa] |
| --- | --- | --- | --- | --- | --- |
| hex | 5 | 5.62e-3 | 1.29e-2 | 4.54e4 | 1.01e5 |
| hex | 10 | 1.33e-3 | 3.57e-3 | 1.58e4 | 6.43e4 |
| hex | 20 | 3.15e-4 | 9.42e-4 | 5.41e3 | 3.39e4 |
| hex | 40 | 7.30e-5 | 2.39e-4 | 1.85e3 | 1.68e4 |
| hex | net order | 2.09 | 1.92 | 1.54 | 0.86 |
| distHex | 5 | 5.87e-3 | 1.25e-2 | 4.25e4 | 9.07e4 |
| distHex | 10 | 1.57e-3 | 4.76e-3 | 1.79e4 | 6.16e4 |
| distHex | 20 | 3.74e-4 | 1.30e-3 | 6.73e3 | 4.12e4 |
| distHex | 40 | 9.09e-5 | 3.19e-4 | 2.76e3 | 2.38e4 |
| distHex | net order | 2.00 | 1.76 | 1.31 | 0.64 |
| structTet | 5 | 8.20e-3 | 1.84e-2 | 4.77e4 | 1.14e5 |
| structTet | 10 | 1.65e-3 | 4.51e-3 | 1.53e4 | 5.83e4 |
| structTet | 20 | 3.89e-4 | 1.18e-3 | 6.30e3 | 2.73e4 |
| structTet | net order | 1.96 | 1.76 | 1.30 | 0.92 |

With the cubic high-order reconstructions, `movingLeastSquares` (MLS) and
`kExactLeastSquares` (kExact), on the regular hexahedral meshes:

| Method | Cells/side | D L2 [m] | D Linf [m] | sigma L2 [Pa] | sigma Linf [Pa] |
| --- | --- | --- | --- | --- | --- |
| MLS | 5 | 1.43e-3 | 3.18e-3 | 7.73e3 | 2.50e4 |
| MLS | 10 | 5.86e-5 | 1.62e-4 | 6.36e2 | 5.09e3 |
| MLS | 20 | 2.90e-6 | 7.02e-6 | 5.81e1 | 4.12e2 |
| MLS | 40 | 1.91e-7 | 5.26e-7 | 5.47e0 | 3.08e1 |
| MLS | net order | 4.29 | 4.19 | 3.49 | 3.22 |
| kExact | 5 | 1.22e-3 | 4.52e-3 | 5.81e3 | 1.78e4 |
| kExact | 10 | 6.53e-5 | 1.50e-4 | 8.69e2 | 3.59e3 |
| kExact | 20 | 4.42e-6 | 9.31e-6 | 8.17e1 | 2.97e2 |
| kExact | 40 | 3.04e-7 | 7.30e-7 | 7.55e0 | 3.31e1 |
| kExact | net order | 3.99 | 4.20 | 3.20 | 3.02 |

With the second-order discretisation the displacement converges at second
order and the cell-centred stress, which involves the gradient of the
displacement, at between first and second order in the L2 norm, as for the
small-strain `manufacturedSolution` tutorial. The cubic reconstructions give
fourth order in displacement and about third order in stress, and on the
$40^3$ mesh their errors are 300 to 400 times smaller.

### Temporal convergence (`dynamic` mode)

On the fixed $10^3$ mesh, the error with respect to the exact solution levels
off at the spatial error of about $1.5 \times 10^{-3}$ m, so the temporal
order is measured from the differences between the final displacement fields
of successive time-step sizes. The observed orders from the L2 norm of these
differences are, for the end time 0.05 s divided into the given numbers of
time steps:

| Scheme | Time function | 20-40-80 | 40-80-160 | 80-160-320 | 160-320-640 |
| --- | --- | --- | --- | --- | --- |
| Euler | cosineSquared | 0.62 | 0.77 | 0.88 | 0.94 |
| backward | cosineSquared | 0.54 | 1.57 | 1.84 | 1.93 |
| NewmarkBeta | cosineSquared | 1.90 | 1.97 | 1.99 | 2.00 |
| Euler | cosine | 3.59 | -0.69 | 0.22 | 0.70 |
| backward | cosine | 1.84 | 1.98 | 2.00 | 2.01 |
| NewmarkBeta | cosine | -0.45 | 0.52 | 0.83 | 0.93 |

With the `cosineSquared` time function, whose displacement, velocity, and
acceleration are all zero at $t = 0$, every scheme shows its design order:
first for `Euler` and second for `backward` and `NewmarkBeta`. With the
`cosine` time function of the dynamic mode, the initial acceleration is
nonzero, $\rho_0 \boldsymbol{a} S \omega_t^2$. The three-level `Euler` and
`backward` schemes are self-starting and keep their orders, although `Euler`
needs more than 320 steps to reach its asymptotic range, while `NewmarkBeta`
starts from a zero stored acceleration and the resulting velocity error makes
its displacement error first order. The exact initial acceleration can be
supplied to the Newmark scheme in a `NewmarkA(D)` field at the start time.

---

## Verification and Regression

The regression test, which `tutorials/Alltest-regression` also runs, runs the
default coarse case in both modes with the segregated and PETSc SNES
approaches and checks the final displacement and stress L2 errors:

```bash
./regressionTest.sh
```

The opt-in [`verification/`](verification/) directory runs the spatial and
temporal convergence sweeps summarised above. It is separate from the normal
tutorial regression suite because it runs many cases. See the verification
README for variants, commands, and acceptance criteria.
