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

## Notes

- The partitioned schemes share the common `fluidSolidInterface` machinery for
  transferring tractions and displacements across the interface.
- Different schemes have different stability and cost trade-offs. In general,
  `weakCoupling` is the cheapest, `fixedRelaxation` and `Aitken` are simple
  strong-coupling options, and `IQNILS` is the most sophisticated partitioned
  acceleration scheme in this directory.
- For `IQNILS`, increasing `couplingReuse` can improve convergence, but it is
  not guaranteed to make a case monotonically more robust. The best reuse level
  is case dependent and interacts with time-step size and the quality of the
  inner fluid and solid solves.
