---
sort: 6
---

# Fluid-solid coupling schemes: `fluidSolidInterfaces`

---

Fluid-solid coupling schemes are implemented in solids4foam via the `fluidSolidInterface`
base class and derived runtime-selectable fluid-solid interface coupling schemes.

The base class is implemented in
`fluidSolidInterface/fluidSolidInterface.{H,C}`. It manages the fluid and solid
models, interface mapping, residual evaluation, outer-correction controls, and
shared coupling data.

---

## Available Schemes

The following schemes are available through the `fluidSolidInterface` entry in
`constant/fsiProperties`.

### `weakCoupling`

Weak Dirichlet-Neumann coupling without outer FSI correctors in each time step.
This is the cheapest option, but it is also the least robust for strongly
coupled problems.

Implementation:
`weakCouplingInterface/weakCouplingInterface.{H,C}`

### `fixedRelaxation`

Strong Dirichlet-Neumann coupling with a fixed under-relaxation factor. This is
a simple strong-coupling scheme and is often used as a baseline or for mildly
coupled problems.

Implementation:
`fixedRelaxationCouplingInterface/fixedRelaxationCouplingInterface.{H,C}`

### `Aitken`

Strong Dirichlet-Neumann coupling with Aitken dynamic under-relaxation. This
improves on fixed relaxation by adapting the relaxation factor during the outer
iterations.

Implementation:
`AitkenCouplingInterface/AitkenCouplingInterface.{H,C}`

### `IQNILS`

Strong Dirichlet-Neumann coupling accelerated with the interface quasi-Newton
inverse least-squares (IQN-ILS) method. This scheme reuses secant information
from previous coupling iterations and time steps to improve convergence for
challenging FSI problems. The implementation includes filtering of repeated or
near-linearly-dependent secant modes before the least-squares solve to improve
numerical robustness.

Common `IQNILS` controls in `constant/fsiProperties` include:

- `relaxationFactor`: startup relaxation used before enough secant information
  is available.
- `couplingReuse`: number of previous time steps whose secant information is
  retained.
- `minSignificant`: absolute filtering tolerance for repeated or
  near-dependent modes.
- `relMinSignificant`: optional relative QR2-style filtering tolerance based on
  the fraction of new information in each older secant mode.
- `normalizeCouplingColumns`: optionally normalize each `V/W` secant column
  pair before the QR solve to improve the conditioning of the least-squares
  system.
- `residualSumPreconditioning`: optionally scale each `V` column by the inverse
  of the current residual magnitude before the QR solve, mimicking preCICE's
  residual-sum approach.
- `residualPreconditioningEpsilon`: regularization added to the residual
  magnitude when computing the above scaling weight.
- `reusePreviousStepModes`: optionally reuse the previous time step's cached
  secants only in the first IQN-ILS solve of the new time step. If the
  immediately previous time step converged in one iteration and generated no
  new secants, the implementation falls back to the latest non-empty cached
  time step instead, which is closer to preCICE's backup LS-system behavior.
- `maxReuseUpdateNormRatio`: optional safeguard that caps the norm of that
  first reused IQN-ILS correction relative to the current interface residual.
- `combinedCouplingSystem`: optionally assemble one aligned least-squares
  system across all coupled interfaces instead of solving one interface at a
  time. On single-interface cases this should reproduce the per-interface
  result; it is mainly a structural step toward a more preCICE-like assembly.
  This can also be combined with `preciceStyleCouplingQR`.
- `requireAllResidualMeasures`: optionally require all normalized residual
  measures in the shared FSI convergence check to satisfy the tolerance.
  Legacy behavior uses the minimum of the available residual measures, which
  can stop earlier when only one measure has decayed.
- `qrSolveTolerance`: relative cutoff used to ignore nearly singular
  directions during the triangular back-substitution.
- `reorthogonalizeCouplingColumns`: optionally apply a second
  Gram-Schmidt-style correction while assembling the QR system.
- `predictSolid`: optionally solve the solid once before the outer FSI loop.

Implementation:
`IQNILSCouplingInterface/IQNILSCouplingInterface.{H,C}`

### `oneWayCoupling`

Pseudo-coupling scheme for one-way FSI. The fluid solution is assumed to be
known already, and the fluid fields are read from a pre-run fluid case and
applied to the solid.

Implementation:
`oneWayCouplingInterface/oneWayCouplingInterface.{H,C}`

### `thermal`

Strong thermal coupling interface with fixed under-relaxation. This extends the
mechanical coupling infrastructure to temperature and heat-flux transfer and
can optionally include mechanical coupling as well.

Implementation:
`thermalCouplingInterface/thermalCouplingInterface.{H,C}`

---

## Immersed interfaces

A solid interface patch may drive an immersed body of the fluid (the
`immersedBoundaryForce` finite volume option of `src/immersedBoundary`,
OpenFOAM.com only) instead of a fluid patch: the fluid mesh is then fixed and
needs no motion solver or interface mapping. The fluid patch of the interface
is `none`, and the `immersedInterfaces` dictionary, keyed by the solid patch
name, gives the body, whose motion must be `fsiDriven`, and optionally the
`closurePatches`: other patches of the solid (e.g. a clamped root) such that
the union of the interface and closure patches is a closed surface with
outward normals. In a two-dimensional fluid mesh the patches of the empty
direction are not needed: the surface is extended through the mesh and
capped in that direction.

```c++
solidPatch      plate;
fluidPatch      none;

immersedInterfaces
{
    plate
    {
        body            flag;          // body of the fluid fvOptions
        closurePatches  (plateFix);    // optional
    }
}
```

The lists `solidPatches`/`fluidPatches` may mix immersed (`none`) and
body-fitted interfaces. In each coupling iteration the immersed surface is
moved to the solid interface displacement of the coupling scheme (relaxed or
accelerated as for a fluid patch), with the velocity of the displacement
increment over the time step, and the traction on the solid patch is the
surface traction of the immersed boundary averaged over the quadrature
points of each face. The log compares its total with the momentum exchange
between the fluid and the body. The immersed interfaces work with
`fixedRelaxation`, `Aitken` and `IQNILS`, and not with `weakCoupling`,
`oneWayCoupling`, `thermal` or the Robin (`elasticWallPressure`) interface
conditions. See the `fluids/immersedBoundary/immersedHronTurekFsi2` tutorial.

Implementation: `fluidSolidInterface/fsiImmersedInterface.{H,C}` builds the
surface and maps the tractions, and `fluidSolidInterface/fsiImmersedBoundary.H`
is the header-only interface implemented by `immersedBoundaryForce`.

---

## Notes

- The partitioned schemes share the common `fluidSolidInterface` machinery for
  transferring tractions and displacements across the interface.
- Different schemes have different stability and cost trade-offs. In general,
  `weakCoupling` is the cheapest, `fixedRelaxation` and `Aitken` are simple
  strong-coupling options, and `IQNILS` is the most sophisticated partitioned
  acceleration scheme in this directory.
- `fixedRelaxation`, `Aitken` and `IQNILS` start each time step from the
  solid predictor (`predictSolid`, default `yes`): the solid is solved once
  with the fluid force of the previous time step. With `predictor` (default
  `yes`), the first iteration moves the fluid interface with the full
  predicted solid displacement; later iterations are relaxed or accelerated as
  usual. The interface displacement is an increment within the time step, so
  with `predictor no` the first iteration moves the fluid interface by only
  `relaxationFactor` times the solid's step, and the interface stops, or with
  backward time differencing reverses, as seen by the fluid. On cases with a
  strong added-mass effect and a small `relaxationFactor`, the resulting
  first-iteration fluid force can be hundreds of times the converged one and
  make the solid diverge (#489).
- For `IQNILS`, increasing `couplingReuse` can improve convergence, but it is
  not guaranteed to make a case monotonically more robust. The best reuse level
  is case dependent and interacts with time-step size and the quality of the
  inner fluid and solid solves.
