# Splitting this branch into reviewable pull requests

A proposal, not a decision. Measured against `origin/development` at
`a029e4be`: 176 files, +38 800 / -2 500.

That is not one pull request. It is also not one piece of work: about a tenth
of it fixes things that were already wrong and would be worth merging if the
framework did not exist at all, and those should not wait behind seventeen
thousand lines of new subsystem.

The ordering below is by dependency, and every stage is meant to build on all
three forks and leave the sixteen-tutorial suite passing on its own.

## Stage A - fixes that do not depend on the framework

These touch the legacy code and stand alone. They should go first, separately,
because each changes behaviour and each deserves to be seen on its own.

| | what | size | results move? |
|---|---|---|---|
| A1 | `anisotropicBiotElastic` chose its reduced plane model when z **is** solved - inverted. Also refuses a mesh empty in x or y, and a 2-D case asking for plane strain from a plane-stress reduction. | +44 / -4, plus the `rodAndSeabed` reference | **yes** - see below |
| A2 | `normalDisplacement` and `normalDisplacementZeroShear` read `normalDisp`, not the documented keyword (issue #409) | +11 / -4 | no |
| A3 | The pore pressure field is called `p`. Rename to `porePressure`. | small, but touches tutorial `0/` and `boundaryData` | no |

A1 moves `rodAndSeabed`: `epsilonEq` 9.31e-5 to 1.78e-3, `sigmaEq` 80.1 kPa to
49.1 kPa. The old numbers were an artifact. The tutorial declares `Ez`,
`nuyz`, `nuzx`, `Gyz` and `Gzx` and was ignoring all five, and because
`poroMechanicalLaw` seeds the effective stress as `sigma + b*(p + p0)*I`, a
sub-law that wrote only xx and yy left an effective stress of `diag(0, 0, p)`.
With this case's initial pore pressure of 79.29 kPa that is a von Mises stress
of 79.29 kPa at zero strain, against the 80.1 kPa observed.

## Stage B - the framework, with nothing using it

Nothing selects any of this. The solid-model switch that reaches it does not
exist until stage C, so the merged result is byte-identical at runtime.

| | what | size |
|---|---|---|
| B1 | Integration-point topologies, state, state spec, kinematics, inputs, response, tangent request | ~2 000 |
| B2 | The manager: law selection, evaluation, tangents | ~5 100 |
| B3 | Constitutive state restart: `stateIO`, the decomposition identity in the file header, and the tests | ~1 550 |
| B4 | `Test-mechanicalConstitutiveLaw` | ~1 980 |
| B5 | `README.md`, `DESIGN-tangents.md`, `DESIGN-state-io.md` | ~5 000 |

B5 is documentation and could ride with B2, but it is what makes B2 reviewable,
so it should not land after it.

## Stage C - the laws

32 files, +7 800, and each law is independent of every other. Group them so a
reviewer can hold one physics question at a time:

- **C1** linear elastic, and the small-strain plastic laws
- **C2** the hyperelastic laws: neo-Hookean, Mooney-Rivlin, St Venant-Kirchhoff
- **C3** the composites: poro, thermo, electro
- **C4** the fibre laws: Guccione, and `setFibreField`

Every law arrives with the unit-test coverage that proves it reproduces the
legacy law it ports, so these are self-checking.

## Stage D - solid models opt in, one at a time

Each adds `useMechanicalConstitutiveLawManager`, defaulting off, and a
regression arm that runs the same tutorial both ways and requires them to
agree. This is where a reviewer can actually see the framework do something,
and where a mistake would show.

| | what | size |
|---|---|---|
| D1 | `solidModel` base plumbing and `linGeomTotalDispSolid` | ~800 |
| D2 | `nonLinGeomTotalLagTotalDispSolid`, `nonLinGeomUpdatedLagSolid` | ~700 |
| D3 | `poroLinGeomSolid`, `thermalLinGeomSolid` | ~280 |
| D4 | the two vertex-centred models | ~2 060 / -1 880 |

D4 is a rewrite rather than an addition and should go last of these.

## Stage E - what the framework then makes possible

| | what | depends on |
|---|---|---|
| E1 | The mixed displacement-pressure formulation on the declared isochoric/volumetric split | D2, C2, C4 |
| E2 | The `LandEtAl2015/problem3` tutorial | C4, E1 |

E1 is the one that changes an answer rather than reproducing one. On problem3,
with `dev()` projection `|D|` is 0.00113675 and with the declared split it is
0.00115475 - a 1.6% difference which is entirely the spherical part of the
active tension that a projection discards and the pressure never restores.
Two independent reviews agreed the split is the correct reading for this
formulation, on the grounds that the pressure residual ties `p` to `dU/dJ` and
to nothing else, so a term with no potential behind it belongs in the stress.

Nothing had ever run that combination, and getting it to run turned up four
things, which is worth knowing before this stage is reviewed:

  - `GuccioneElastic` ignored `tangentRequest::scalarDeviatoric` and returned
    an implicit stiffness about seventy times too large for a mixed
    formulation, so the linear solve never converged.
  - the response's fourth-order-tangent constructor left `volumetricPtr_`
    uninitialised, which nothing had dereferenced before.
  - the high-order Jacobian recovers `mu` from `impK` assuming the
    displacement form, and would have got a negative shear modulus.
  - the pressure equation hard-codes the volumetric energy it assumes, which
    `providesVolumetricSplit()` does not promise. It is now checked against
    what the law returns.

## What is not ready, and should not be in any of these

- The small-strain path has no isochoric/volumetric split, so
  `linGeomTotalDispSolid` still projects with `dev()`. Exact for isotropic
  linear elasticity, wrong for an anisotropic small-strain law under a solved
  pressure. No case exercises it.
- Point-collapse accumulators are not synchronised across processor
  boundaries (`syncTools::syncPointList`), and there is no shared-point
  precedence policy. No active tutorial reaches it.
- A law's `endTimeStep` hook may call `reduce()`, so boundary states cannot be
  visited there without deadlocking on differing patch counts per rank.
- Nine laws have no framework port and no tutorial to validate one against:
  `diffusionElastic`, `diffusionHyperElastic`, `orthotropicLinearElastic`,
  `StVenantKirchhoffOrthotropicElastic`, `GentElastic`, `isotropicFungElastic`,
  `YeohElastic`, `viscoNeoHookeanElastic`, and `HolzapfelGasserOgdenElastic`
  (which does have `ratCarotid`).

## Deprecating the legacy path

Not yet. The legacy `mechanicalLaw` hierarchy is still what every tutorial
runs by default, and the framework is opt-in everywhere. A deprecation needs
the switch defaulted on for at least one release with the regression arms
demanding agreement, which is what stage D builds. The nine unported laws
above are the actual blocker.
