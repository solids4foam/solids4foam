# Splitting this branch into pull requests

Measured against `origin/development`, not estimated:

    220 files changed, 42660 insertions(+), 662 deletions(-)     207 commits

## The two facts that make this reviewable

**82% of it is new files.** 152 of the 220 files are additions - 36 656 of the
42 660 lines - and a new file cannot change what an existing case does. The
part where behaviour can change at all is 64 modified files, +6 002 / -660.

**Exactly one tutorial's published numbers change.** `rodAndSeabed`, and only
because `anisotropicBiotElastic` was selecting its reduced plane model on
three-dimensional meshes. Every other bound in the suite is either untouched or
*added* alongside the existing one - `perforatedPlate` keeps its legacy
28-to-32 yielding cells and gains a separate 38-to-44 for the framework, which
counts boundary integration points as well.

So a reviewer's attention belongs almost entirely on two small things: the
standalone fixes in stage A, and the opt-in wiring in stages D and E. The bulk
is inert.

## Why the bulk is inert

The framework is opt-in everywhere. `useMechanicalConstitutiveLawManager`
defaults to false, no shipped tutorial sets it except in arms added by this
branch, and with it unset every code path is the one that ran before. Stages B
and C can therefore be merged on "it builds on three forks and the suite still
passes" without anyone having to reason about numerics - which is the whole
point of ordering them first.

## Stage A - fixes that stand alone

Independent of the framework, and each changes behaviour, so each wants its own
review rather than being buried.

| | what | size | results move? |
|---|---|---|---|
| A1 | `anisotropicBiotElastic` chose its reduced plane model when z **is** solved - inverted. Also refuses a mesh empty in x or y, and a 2-D case asking for plane strain from a plane-stress reduction. | +44 / -4 plus the `rodAndSeabed` bounds | **yes** |
| A2 | `normalDisplacement` and `normalDisplacementZeroShear` read `normalDisp`, not the documented keyword (issue #409) | +11 / -4 | no |
| A3 | The pore pressure field is called `p`; rename to `porePressure` | small, touches tutorial `0/` and `boundaryData` | no |
| A4 | The legacy `HolzapfelGasserOgdenElastic` gains `shearModulus()`, previously `notImplemented`, which is much of why it only ever ran under one solid model | +28 | no |
| A5 | `solidModel::newDeltaT` removed - no caller, no override | -9 | no |

A1 moves `rodAndSeabed`: `epsilonEq` 9.31e-5 to 1.78e-3, `sigmaEq` 80.1 kPa to
49.1 kPa. The old numbers were an artifact. The tutorial declares `Ez`,
`nuyz`, `nuzx`, `Gyz` and `Gzx` and was ignoring all five, and because
`poroMechanicalLaw` seeds the effective stress as `sigma + b*(p + p0)*I`, a
sub-law writing only xx and yy left `diag(0, 0, p)`. With this case's initial
pore pressure of 79.29 kPa that is a von Mises stress of 79.29 kPa at zero
strain, against the 80.1 kPa observed.

## Stage B - the framework, with nothing using it

| | what | size |
|---|---|---|
| B1a | Integration-point topologies | +1 486 |
| B1b | State, state spec, kinematics, inputs, response, tangent request, diagnostics | +2 779 |
| B2 | The manager: law selection, evaluation, tangents | +5 469 |
| B3 | Constitutive state restart: `stateIO` and the decomposition identity in the file header | +1 536 |
| B4 | `Test-mechanicalConstitutiveLaw` | +2 151 |
| B5 | `README.md`, `DESIGN-tangents.md`, `DESIGN-state-io.md` | +5 574 |

B5 is documentation, and it is what makes B2 reviewable, so it should not land
after it.

## Stage C - the laws

34 files, +8 539, each law independent of every other. Grouped so a reviewer
holds one physics question at a time:

- **C1** linear elastic, and the small-strain plastic laws
- **C2** the hyperelastic laws: neo-Hookean, Mooney-Rivlin, St Venant-Kirchhoff
- **C3** the composites: poro, thermo, electro
- **C4** the fibre laws: Guccione, Holzapfel-Gasser-Ogden, and `setFibreField`

Every law arrives with the unit-test coverage that shows it reproduces the
legacy law it ports, so these are largely self-checking. C4's two laws are the
exception and carry their own closed-form checks instead: Guccione because its
reformulation on the isochoric strain is deliberately *not* the legacy law away
from incompressibility, and HGO because its tutorial does not run.

## Stage D - solid models opt in, one at a time

Each adds `useMechanicalConstitutiveLawManager`, defaulting off, and a
regression arm running the same tutorial both ways. This is where a reviewer
can see the framework do something, and where a mistake would show.

| | what | size |
|---|---|---|
| D1 | `solidModel` base plumbing and `linGeomTotalDispSolid` | +1 081 / -67 |
| D2 | `nonLinGeomTotalLagTotalDispSolid`, `nonLinGeomUpdatedLagSolid` | +1 163 / -40 |
| D3 | `poroLinGeomSolid`, `thermalLinGeomSolid` | +364 / -13 |
| D4 | `vertexCentredLinGeomSolid` - Jacobian tangent only, not a stress port | +182 / -17 |

D4 is deliberately narrow and its header says so: the switch moves only the
tangent, the residual stress still comes from `dualMechanicalModel`, so the
converged answer is the legacy one and only the convergence path changes.

The one piece of D that is not a switch: a solid model on the framework
computes its own gradient and point interpolation rather than calling
`mechanical()`, because for more than one material those route through the
legacy per-material subMeshes that the framework exists to replace. On
`layeredPipe` that is worth 0.0192 against 0.0494 in radial stress error, on a
tolerance of 0.03. Multi-material with a gradient scheme that is not
material-aware is refused rather than warned about.

## Stage E - what the framework then makes possible

| | what | depends on |
|---|---|---|
| E1 | The mixed displacement-pressure formulation on the declared isochoric/volumetric split, at finite strain and small strain | D1, D2, C2, C4 |
| E2 | The `LandEtAl2015/problem3` tutorial | C4, E1 |
| E3 | `ratCarotid` on v2512 - the case does not run at all on its legacy path | C4, E1 |

E1 is the only stage that changes an answer rather than reproducing one, and it
has evidence in both directions:

  - `problem3`, where the split *must* change the answer: `dev()` projection
    gives 0.00113675 and the declared split 0.00115475, and the 1.6% between
    them is the spherical part of the active tension a projection discards.
  - `plateHole`, where it must change *nothing*: for isotropic linear
    elasticity a projection and the split are the same operation, and the two
    arms agree to every digit reported.

## What is deliberately left undone

- Eight laws have no framework port and no tutorial: `diffusionElastic`,
  `diffusionHyperElastic`, `orthotropicLinearElastic`,
  `StVenantKirchhoffOrthotropicElastic`, `GentElastic`, `isotropicFungElastic`,
  `YeohElastic`, `viscoNeoHookeanElastic`. Flagged rather than ported: a port
  nothing can check is how a wrong port gets merged looking right. A later PR
  can take any of them when a case exists to validate it.
- Point-centred stress collapse refuses on a decomposed mesh rather than
  guessing. Nothing calls it - section 24.1 explains why the obvious fix is
  wrong.
- The updated-Lagrangian mixed formulation is refused, not approximated.
- `coupledPressureDisplacementSolid` and `coupledUnsLinGeomLinearElasticSolid`
  are not ported. Both are foam-extend-only in practice; see below.

## Deprecating the legacy path

Not from this branch. Two things have to happen first, and neither is large.

`coupledPressureDisplacementSolid` is the substantial remaining model - five
tutorials, and foam-extend only. Its pressure equation already carries a
penalty term, so the shapes look compatible, but it solves a pressure
*increment* and accumulates `p = p.oldTime() + Dp` with no residual visibly
enforcing the framework's `p = -dU/dJ`. That increment relation needs deriving
before porting, not a coefficient adjusted.

`coupledUnsLinGeomLinearElasticSolid` is best left. It downcasts to
`linearElastic`, `FatalError`s on anything else, reads `mu()` and `lambda()`
off it to build a block matrix, and is single-material by construction - so
none of the framework's value applies to it, and its `evolve()` body is
foam-extend-only anyway.

Then the blocker is the eight unported laws, and that is a question rather than
a task: port them blind and mark them unvalidated, drop them, or write cases
for them. That is a decision about what the library supports.

## One practical note for whoever splits this

A single edit of mine normalised the line endings of
`vertexCentredNonLinGeomTotalLagSolid.C`, which is stored CRLF upstream. That
made a sixteen-line change read as 1 875 insertions against 1 860 deletions -
the largest diff on the branch, and entirely phantom. It is fixed, and it took
the branch's deletions from 2 521 to 662. Worth checking for again if any
further rebasing happens: `git diff -w --shortstat` against the same paths
will show it immediately.
