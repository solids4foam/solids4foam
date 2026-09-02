# Design: material stiffness and Jacobian tangents for solid models

Status: proposal. No implementation code has been written.
Applies to: branch `feature-mechanicalConstitutiveLaw` (PR #198), commit `5c4a03c7`.
Resolves the open question in `README.md`, section *Tangents and mixed formulations*
(lines 234-253).

**Section 8 records decisions taken after review and supersedes the earlier text
where it says so. Read it alongside sections 3 and 5.**

---

## 0. Summary of the decision

`impK` is doing two unrelated jobs, and the code already knows this even though
the name does not.

1. It is the **coefficient of a discretisation splitting**: the `gamma` in
   `fvm::laplacian(gamma, D) - fvc::laplacian(gamma, D)`, the scale of the
   stabilisation traction `impKf_*momentumStabilisation().faceVector()`, and the
   divisor in `tractionBoundarySnGrad`. This is a property of the *solid model's
   discretisation*, evaluated on cells and faces of the solution mesh.

2. It is a **surrogate Jacobian**: the coefficient handed to
   `momentumStabilisation().vectorJacobian(D, &impKf_)`, to `precondition()`,
   and back-derived into `mu`/`lambda` for `hofvm::laplacian*IntoPETScMatrix`.

The recommendation is therefore:

> **Keep (1) as a solid-model-owned scalar field pair (`impK_`, `impKf_`,
> `rImpK_`) whose values are obtained from the manager via
> `tangentRequest::scalar`, with today's lifecycle preserved exactly. Make (2)
> a per-solid-model dictionary choice `jacobianTangent`, which subsumes and
> replaces `approximateJacobian`, and which is orthogonal to (and does not
> subsume) `highOrderJacobian`.**

The user's proposed separation is **confirmed with one important correction**:
`impK` is *not* free to change in the solvers that add a stabilisation traction.
See §1.3. This is the binding constraint on "no change to numerical results".

---

## 1. The central analysis

### 1.1 Role 1 is real: the `fvm`/`fvc` pair cancels for any coefficient

`linGeomTotalDispSolid.C:229-231`

```c++
tmp<fvVectorMatrix> tRhsEqn
(
    fvm::laplacian(impKf_, D(), "laplacian(DD,D)")
);
tmpRef(tRhsEqn) -= fvc::laplacian(impKf_, D(), "laplacian(DD,D)");
```

Same field, same scheme name, same operator. At the fixed point
`D_new == D_old` the two are identical to round-off and cancel exactly. The
value of `impKf_` affects only the spectral radius of the segregated iteration.
It must, however, be *identical* on both sides: a refreshed `impKf_` between the
`fvm` and `fvc` calls would leave a spurious residual. The same pattern appears
in `unsLinGeomSolid.C:127-129`, `unsNonLinGeomTotalLagSolid.C:269-271`,
`unsNonLinGeomUpdatedLagSolid.C:274-276`, `thermalLinGeomSolid.C:320-322`,
`poroLinGeomSolid.C:302-304`, `nonLinGeomTotalLagTotalDispSolid.C:402-404`
and `nonLinGeomUpdatedLagSolid.C:363-365`.

### 1.2 Role 2 is real: the Jacobian is a pure convergence-rate object

`linGeomTotalDispSolid.C:1227-1245` assembles `approxJ` and inserts it into the
PETSc matrix; `linGeomTotalDispSolid.C:1284` builds an entirely separate
Laplacian for `precondition()`. Neither appears in `evaluateResidual`. In a SNES
(or any inexact-Newton) setting the converged solution is defined solely by
`F(x) = 0`; the Jacobian controls only the path taken to get there.

**This is the licence to fix the `highOrderJacobian` hack (§4) without violating
the "no change to results" constraint** — provided the nonlinear solve is
actually converged to tolerance. With a loose `snes_rtol` or a capped iteration
count, a different Jacobian does change the answer. Regression checks must
therefore compare converged fields *and* report iteration counts separately.

### 1.3 Correction: `impKf_` is **not** free in the stabilised cell-centred solvers

This is the finding that most affects the migration plan.

`linGeomTotalDispSolid.C:220` (segregated), `:533` (explicit acceleration) and
`:1020` (SNES residual) all do:

```c++
momentumStabilisation().updateVector(D, &gradD());
traction += impKf_*momentumStabilisation().faceVector();
```

The default momentum stabilisation is `diffStencilLaplacian` with
`scaleFactor 0.1` (`solidModel.C:1327-1329`), whose face vector is
(`diffStencilLaplacianStabTemplates.H:42-43`):

```c++
faceFlux = scaleFactor*(fvc::snGrad(field) - (n & fvc::interpolate(gradField)));
```

The compact `snGrad` and the interpolated full gradient do **not** agree on a
discrete mesh, so this term is non-zero at convergence. It enters the residual
through `fvc::div(mesh().magSf()*traction)`. Therefore `impKf_` scales a
solution-changing term: **changing `impKf_` changes the converged discrete
answer** in `linGeomTotalDispSolid`, `nonLinGeomTotalLagTotalDispSolid` and
`nonLinGeomUpdatedLagSolid`.

Two independent pieces of evidence that the codebase already knows this:

* `thermalLinGeomSolid.C:325` and `poroLinGeomSolid.C:307` use the **frozen cell
  field** `impK_` for the stabilisation, not `impKf_`, and consequently they are
  the only solvers that dare to refresh `impKf_` mid-solve
  (`thermalLinGeomSolid.C:359-365`):

  ```c++
  // Update impKf to improve convergence
  // Note: impK and rImpK are not updated as they are used for traction
  // boundaries
  if (iCorr % 10 == 0)
  {
      impKf_ = mechanical().impKf();
  }
  ```

* `linGeomTotalDispSolid.H:72` declares `const volScalarField impK_` and
  `:80` `const volScalarField rImpK_` — deliberately immutable, built once at
  `linGeomTotalDispSolid.C:571-578` and never refreshed, even for plasticity.

**Consequence for the design:** the *value* of the scalar tangent that feeds
`impK` must reproduce `mechanicalModel::impK()` bit-for-bit at the same point in
the lifecycle. Fortunately it does — see §1.6.

### 1.4 `tractionBoundarySnGrad` *is* free (contrary to first appearance)

`linGeomTotalDispSolid.C:1346-1355`:

```c++
(
    (traction - n*pressure)
  - (n & (pSigma - impK*pGradD))
)*rImpK
```

Rearranged, the fixed point satisfies

```text
impK*(snGrad(D) - (n & gradD_patch)) = traction_applied - (n & sigma_patch)
```

For a single material `mechanicalModel::grad` is `fvc::grad(D)`
(`mechanicalModel.C:665-668`), and every OpenFOAM gradient scheme ends with
`gaussGrad<Type>::correctBoundaryConditions`, which on every **non-coupled**
patch overwrites the normal component of the boundary gradient with `snGrad`
(verified in `$FOAM_SRC/finiteVolume/finiteVolume/gradSchemes/gaussGrad/gaussGrad.C`,
v2512):

```c++
        if (!vsf.boundaryField()[patchi].coupled())
        {
            ...
            gGradbf[patchi] += n *
            (
                vsf.boundaryField()[patchi].snGrad()
              - (n & gGradbf[patchi])
            );
        }
```

Traction boundary conditions are always on non-coupled patches, so
`n & gradD_patch == snGrad(D)_patch` identically, the left side above is exactly
zero, and the fixed point is `n & sigma_patch == traction_applied`
**independently of `impK`**. `impK` and `rImpK` here are pure relaxation
parameters.

The multi-material branch reaches the same place explicitly, by calling
`fv::gaussGrad<vector>(mesh()).correctBoundaryConditions(D, gradD)` at
`mechanicalModel.C:696-699`.

This is precisely what the `thermalLinGeomSolid.C:360` comment asserts, and it is
why freezing `impK_`/`rImpK_` while refreshing `impKf_` is legitimate *there*.

### 1.5 Consumers that straddle

<!-- markdownlint-disable MD013 -->

| Consumer | Straddles? | Why |
|---|---|---|
| `traction += impKf_*...faceVector()` | **Yes** | Written as a residual term; behaves as a stabilisation coefficient; changes converged answer. |
| `rAUf()` from `fvm::laplacian(impKf_,D)` (`linGeomTotalDispSolid.C:980-988`) | **Yes** | Built like a preconditioner (`approxMomJ.A()`), but feeds `pressureStabilisation().cellScalar(&rAUf(), true)` inside the **pressure residual** (`:994`). Solution-changing when `solvePressure()`. |
| `sqrt(impK/rho)` time step (`:741`) | **Yes** | Only used for `solutionAlgorithm::EXPLICIT`. Does not enter any equation, but it sets `deltaT`, so any change changes the answer. |
| `tractionBoundarySnGrad` | No | Cancels at the fixed point (§1.4). |
| `precondition()` (`:1284`) | No | Pure preconditioner. |
| `momentumStabilisation().vectorJacobian(D, &impKf_)` (`:1229`) | No | Jacobian only. |
| `fixedDofScale_` (`vertexCentredLinGeomSolid.C:816-826`) | No | Diagonal magnitude for `MatZeroRowsColumnsIS`; conditioning only. |

<!-- markdownlint-enable MD013 -->

### 1.6 The new framework already reproduces `impK` exactly

<!-- markdownlint-disable MD013 -->

| Law | Legacy `impK()` | New `tangentRequest::scalar` | Match |
|---|---|---|---|
| `linearElastic` (nu < 0.5) | `2*mu_ + lambda_` (`linearElastic.C:268`) | `2*mu + lambda` (`linearElasticMechanicalConstitutiveLaw.C:127`) | exact |
| `linearElastic` (nu == 0.5) | `2*mu_` (`linearElastic.C:249`) | n/a — new law rejects nu >= 0.5 - SMALL (`:71-78`) | see OQ-9 |
| `linearElasticMisesPlastic`, elastic point | `(4/3)mu + K` | `(4/3)mu + kappa` (`...MisesPlastic...C:310`) | exact |
| `linearElasticMisesPlastic`, plastic point | `scaleFactor*(4/3)mu + K` (`linearElasticMisesPlastic.C:739-777`) | `theta*(4/3)mu + kappa` (`...C:416`) | exact |
| `linearElastic` `mat66` | `linearElastic.C:295-309` | `linearElasticMechanicalConstitutiveLaw.C:191-205` | identical, same Voigt convention (`mu` on the shear diagonal) |

<!-- markdownlint-enable MD013 -->

`2.0*mechanical().shearModulus()` (the `solvePressure()` branch at
`linGeomTotalDispSolid.C:574`) equals `2*mu`, which is exactly
`tangentRequest::scalarDeviatoric` for `linearElastic`
(`linearElasticMechanicalConstitutiveLaw.C:130-132`).

---

## 2. Classification table of every `impK` consumer

Line numbers are on `feature-mechanicalConstitutiveLaw` @ `5c4a03c7`.
"Free" = the value may change without changing the converged solution.

### 2.1 `solidModels/linGeomTotalDispSolid/linGeomTotalDispSolid.C|.H`

<!-- markdownlint-disable MD013 -->

| Line | Expression | Role | Free? |
|---|---|---|---|
| `.H:72` | `const volScalarField impK_` | storage, immutable by design | — |
| `.H:75` | `surfaceScalarField impKf_` | storage | — |
| `.H:80` | `const volScalarField rImpK_` | storage | — |
| `.C:571-576` | `solvePressure() ? 2*shearModulus() : mechanical().impK()` | construction | — |
| `.C:577` | `impKf_(fvc::interpolate(impK_))` | construction | — |
| `.C:578` | `rImpK_(1.0/impK_)` | construction | — |
| `.C:220` | `traction += impKf_*momentumStabilisation().faceVector()` | **residual (stabilisation scale)** | **No** |
| `.C:229` | `fvm::laplacian(impKf_, D())` | residual-consistency | Yes (must match `:231`) |
| `.C:231` | `-= fvc::laplacian(impKf_, D())` | residual-consistency | Yes (must match `:229`) |
| `.C:533` | `traction += impKf_*...faceVector()` (explicit `A_`) | **residual** | **No** |
| `.C:741` | `sqrt(mechanical().impK()/rho())` | **explicit time-step estimate** | **No** (sets `deltaT`) |
| `.C:983` | `fvm::laplacian(impKf_, D)` -> `rAUf()` -> pressure residual `:994` | **residual (u-p only)** | **No** when `solvePressure()` |
| `.C:1020` | `traction += impKf_*...faceVector()` (SNES residual) | **residual** | **No** |
| `.C:1120` | `fvm::laplacian(impKf_, D)` -> `rAUf()` -> `approxPressureJ` | Jacobian | Yes |
| `.C:1173` | `tMu = (impK_ - K)*(3.0/4.0)` | Jacobian (high-order) | Yes |
| `.C:1176` | `tLambda = impK_ - 2.0*mu` | Jacobian (high-order) | Yes |
| `.C:1229` | `momentumStabilisation().vectorJacobian(D, &impKf_)` | Jacobian | Yes |
| `.C:1284` | `fvm::laplacian(impKf_, D, "preconditionD")` | preconditioner | Yes |
| `.C:1331,1334,1352,1353` | `tractionBoundarySnGrad` | boundary condition | Yes (§1.4) |

<!-- markdownlint-enable MD013 -->

### 2.2 Other solid models (same taxonomy)

<!-- markdownlint-disable MD013 -->

| File | Lines | Role | Free? |
|---|---|---|---|
| `nonLinGeomTotalLagTotalDispSolid.C` | 388, 1155 | residual stabilisation | **No** |
| | 402-404, 417-418, 1099, 1262 | residual-consistency pair | Yes |
| | 727-734 | construction | — |
| | 1315, 1318 | high-order Jacobian back-derivation | Yes |
| | 1371 | `vectorJacobian(D, &impKf_)` | Yes |
| | 1406, 1435 | `tractionBoundarySnGrad` | Yes |
| `nonLinGeomUpdatedLagSolid.C` | 351, 1094 | residual stabilisation | **No** |
| | 363-365, 1034, 1200 | residual-consistency pair | Yes |
| | 669-671 | construction | — |
| | 1252, 1255 | high-order Jacobian back-derivation | Yes |
| | 1308 | `vectorJacobian(DD, &impKf_)` | Yes |
| | 1344, 1410 | `tractionBoundarySnGrad` | Yes |
| `unsLinGeomSolid.C` | 78-80 | construction | — |
| | 127-129 | residual-consistency pair (**no stabilisation term**) | Yes |
| | 209, 230 | `tractionBoundarySnGrad` | Yes |
| `unsNonLinGeomTotalLagSolid.C` | 189-191, 269-271, 440, 464, 488 | as `unsLinGeomSolid` | Yes |
| `unsNonLinGeomUpdatedLagSolid.C` | 218-220, 274-276, 392, 423 | as `unsLinGeomSolid` | Yes |
| `thermalLinGeomSolid.C` | 252-254 construction; 320-322 pair; **325 stabilisation uses `impK_`**; 364 refresh of `impKf_` every 10 iters; 404, 425 BC | mixed, deliberately split | `impKf_` Yes; `impK_` **No** (`:325`) |
| `poroLinGeomSolid.C` | 195-197, 302-304, **307**, 346, 383, 404 | identical structure to `thermalLinGeomSolid` | same |
| `coupledPressureDisplacementSolid` | `.C:466-476` construction; `.C:1688` refresh in `updateTotalFields()`; `completeMomentumEqAssembly.H:22-57,150-297` block-coupled implicit coefficients; `addBlockCoupledBC.H:15-44`; `calcTraction.H:4-6`, `calcForce.H:56-58`; `.C:973-977` (`cellVector(&impKf_,true)` — **residual stabilisation**); `.C:1275, 1316, 1413` | mixed | `.C:977` **No**; rest Yes |
| `vertexCentredLinGeomSolid.C` | 517-537 `dualImpKf()`; **1324** `vfvm::laplacian(..., dualImpKf(), ...)` | **Jacobian only** | **Yes** |
| | 822 `average(mechanical().impK())` -> `fixedDofScale_` | conditioning | Yes |
| | 1037 `sqrt(impK/rho)` | explicit time step | **No** |
| | 1331-1332 `materialTangentFaceField(List<mat66>&)` | **Jacobian, mat66** | Yes |
| `vertexCentredNonLinGeomTotalLagSolid.C` | 760, 1082, 1380, 1684-1685 | same as above | same |

<!-- markdownlint-enable MD013 -->

### 2.3 Consumers outside `solidModels/`

All of these obtain `impK` by **registry lookup of a `volScalarField` literally
named `"impK"`**. This is a hard constraint on any redesign.

<!-- markdownlint-disable MD013 -->

| File | Line | Mechanism | Role | Free? |
|---|---|---|---|---|
| `.../normalContactModels/standardPenalty/standardPenalty.C` | 61, 67-68, 118 | `mesh_.lookupObject<volScalarField>("impK")` | penalty stiffness | **No** (penalty scale changes converged gap) |
| `.../pointStandardPenalty/pointStandardPenalty.C` | 64, 70-71, 121 | same | penalty stiffness | **No** |
| `.../standardPenaltyFriction/standardPenaltyFriction.C` | 66, 72-73, 121 | same | penalty stiffness | **No** |
| `fluidModels/fvPatchFields/elasticWallPressure/...C` | 79-87 | registry `"impK"`, `ap = sqrt(impK/rho)` | FSI Robin coefficient | **No** |
| `crackerFvMesh/.../modeICohesiveZoneModel.C` | 64-69, 95 | registry `"impK"` | cohesive penalty | **No** |
| `crackerFvMesh/.../fixedMixedModeCohesiveZoneModel.C`, `variableMixedModeCohesiveZoneModel.C` | (3 hits each) | registry `"impK"` | cohesive penalty | **No** |
| `materialModels/.../solidSubMeshes.C` | 1667 | registry `"impK"` | sub-mesh mapping | — |
| `numerics/mechanicalEnergies/mechanicalEnergies.C` | 171 | argument | **unused** (only in commented-out laplacian-smoothing energy at `:271`) | — |
| `materialModels/mechanicalModel/mechanicalModel.C` | 354-418 | `impK()`, `impKf()` | provider | — |
| `materialModels/mechanicalModel/dualMechanicalModel.C` | 334-358, 361-413 | `materialTangentFaceField()`, `impKf()` | provider (dual mesh) | — |
| `abaqusUMATs/abaqusUmatLinearElastic` | 4 hits | `impK()` | provider | — |

<!-- markdownlint-enable MD013 -->

**Constraint C1: a `volScalarField` registered under the name `"impK"` must
continue to exist on the solid mesh for the whole run.** Six independent
consumers depend on it, three of them solution-changing.

---

## 3. Recommended design

### 3.1 Ownership and lifecycle (question a)

**Decision: the solid model owns `impK_`, `impKf_` and `rImpK_` exactly as it
does today. The manager never retains a "current tangent".**

Reasons:

* Constraint C1 requires a registered `volScalarField` named `"impK"` on the
  solid mesh. That is a solid-model artefact, not a constitutive one.
* `impKf_` is a *face-interpolated* quantity whose interpolation scheme is a
  discretisation choice (`fvc::interpolate(impK_)` at
  `linGeomTotalDispSolid.C:577` vs the named scheme
  `"interpolate(impK)"` at `mechanicalModel.C:416`). The manager must not own
  an interpolation scheme.
* Different solid models want different refresh cadences and different
  *definitions* (`scalar` vs `scalarDeviatoric`, `linGeomTotalDispSolid.C:573`).
  A manager-owned singleton would have to guess.
* It keeps the one-way dependency the README asserts (lines 60-74).

**The manager gains exactly one new capability to support this: a
state-preserving, stress-free tangent query.** This is required because `impK`
is needed *at construction time*, before any stress update, and the only way to
obtain a tangent today is to also perform a stress update, which would write to
`sigma` and commit constitutive state.

```c++
        //- Evaluate only the requested scalar tangent at cell centres.
        //  The supplied kinematics are evaluated against a copy of each law's
        //  constitutive state, and the stress is written to manager-owned
        //  scratch storage. Neither the caller's stress field nor any
        //  constitutive history variable is modified.
        //  Intended for solid models that require an implicit stiffness
        //  coefficient at construction time, before any stress update has
        //  been performed, or at a chosen refresh cadence between updates.
        void updateScalarTangent
        (
            const volTensorField& gradD,
            const volTensorField& gradD0,
            const scalar dt,
            volScalarField& scalarTangent,
            const tangentRequest tangentReq = tangentRequest::scalar
        );
```

**Correction, from PR-3.** The copy described above is both heavier than
necessary and was not expressible: `mechanicalConstitutiveLawState` had copy
construction and assignment `= delete`d. What is implemented instead is a
*shadow*: a state that aliases the parent's old-time fields, which are
read-only through it, and owns its own current-time fields.

This works because a constitutive law is a pure function of the kinematics and
the old-time state. `evaluate()` reads only `*0` fields and writes only
current-time fields; the plastic law's `sigmaY` is an output, with the trial
value coming from `yieldStress(epsilonPEq0[i])`. History is therefore read-only
during evaluation, so copying it is waste. A shadow costs one set of
current-time fields and no copy of history, and the invariant becomes
structural: a shadow *cannot* write history because it does not own it, and
says so with a `FatalError` if asked.

The same mechanism serves the finite-difference tangent, which is the other
consumer that must not disturb state. See §8.9.

**Lifecycle, stated explicitly (this reproduces today's behaviour exactly):**

<!-- markdownlint-disable MD013 -->

| Field | Built | Refreshed | Invalidated by |
|---|---|---|---|
| `impK_` | solid model constructor, from `updateScalarTangent(gradD(), gradD0(), 0, impK_, req)` where `req = solvePressure() ? scalarDeviatoric : scalar` | **never** | topological mesh change only (`crackerFvMesh`) |
| `impKf_` | constructor, `fvc::interpolate(impK_)` | `linGeom*`, `nonLinGeom*`, `uns*`: **never**. `thermalLinGeomSolid`, `poroLinGeomSolid`: every 10th outer iteration. `coupledPressureDisplacementSolid`: once per `updateTotalFields()`. | as above |
| `rImpK_` | constructor, `1.0/impK_` | **never** | as above |

<!-- markdownlint-enable MD013 -->

The cadence is *not* generalised into a new dictionary entry in this design.
Each solid model keeps its own hard-coded cadence, because that cadence is
currently load-bearing for bit-for-bit reproduction (§1.3) and because
`linGeomTotalDispSolid`'s `impK_` is `const` for a reason. See OQ-3 if you want
it configurable.

Note that this preserves an existing wrinkle: for `linearElasticMisesPlastic`,
`impK()` depends on `DLambda_` (`linearElasticMisesPlastic.C:739-777`), which is
zero on a cold start but non-zero on restart. Freezing at construction therefore
makes `impK_` restart-dependent today. See OQ-1.

### 3.2 Selection mechanism (question b)

**Decision: one new entry, `jacobianTangent`, in the solid model coefficients
dictionary. It replaces `approximateJacobian` entirely. It is orthogonal to
`highOrderJacobian` and does not replace it.**

```text
linearGeometryTotalDisplacementCoeffs
{
    // Fidelity of the material tangent used to build the Jacobian.
    // scalar                      : scalar approximate tangent (default)
    // scalarDeviatoric            : deviatoric-only scalar tangent
    // fourthOrder                 : full 6x6 Voigt material tangent
    // fourthOrderFiniteDifference : as above, by finite differences
    jacobianTangent  scalar;
}
```

Why `approximateJacobian` is subsumed but `highOrderJacobian` is not:

* `approximateJacobian` (`vertexCentredLinGeomSolid.C:1312`,
  `vertexCentredNonLinGeomTotalLagSolid.C:1664`) is a two-valued switch between
  `vfvm::laplacian(..., dualImpKf(), ...)` and
  `vfvm::divSigma(..., materialTangent, ...)`. That is **exactly** the
  scalar-vs-`mat66` axis, expressed as a Boolean. It is a strict subset of
  `jacobianTangent`. Map: `yes -> scalar`, `no -> fourthOrder`.
* `highOrderJacobian` (`solidModel.C:1279-1280`, read from the
  `highOrderCoeffs` sub-dictionary) selects a **different discrete operator**
  (MLS/`hofvm` vs compact `fvm::laplacian`). It is only meaningful alongside
  `highOrderResidual` and only with `solutionAlgorithm PETScSNES`
  (`solidModel.H:595-600`). Folding an operator choice and a tangent-fidelity
  choice into one enum would produce a Cartesian product of names
  (`highOrderFourthOrder`, `compactScalar`, ...) and would put a
  discretisation switch in the material section of the dictionary. Keep them
  separate; they compose:

  <!-- markdownlint-disable MD013 -->

  | `highOrderJacobian` | `jacobianTangent` | Assembled Jacobian |
  |---|---|---|
  | no | `scalar` | `momentumStabilisation().vectorJacobian(D, &impKf_)` — today's default |
  | no | `fourthOrder` | **not supported** — no cell-centred operator consumes `mat66`; fatal error (see §3.4) |
  | yes | `scalar` | `hofvm::laplacian*IntoPETScMatrix` with back-derived `mu`/`lambda` — today's behaviour, retained for one release as deprecated |
  | yes | `fourthOrder` | `hofvm::divSigmaIntoPETScMatrix(C)` — the new correct path (§3.5) |

  <!-- markdownlint-enable MD013 -->

**Defaults are per solid model, not global**, because today's defaults differ:

```c++
// solidModel.H
            //- Requested fidelity of the material tangent used in the Jacobian.
            //  Read from the "jacobianTangent" entry of the solid model
            //  coefficients dictionary. Each solid model supplies its own
            //  default so that existing cases reproduce their current
            //  behaviour without dictionary changes.
            tangentRequest jacobianTangent
            (
                const tangentRequest deflt
            ) const;
```

* `linGeomTotalDispSolid` and every other cell-centred model pass
  `tangentRequest::scalar`.
* `vertexCentredLinGeomSolid` and `vertexCentredNonLinGeomTotalLagSolid` pass
  `tangentRequest::fourthOrder`, matching `approximateJacobian false`.

Backward compatibility (protects `tutorials/solids/linearElasticity/cantilever2d/constant/solidProperties.vertexCentred:28`):

```c++
// In solidModel::jacobianTangent(...)
if (solidModelDict().found("approximateJacobian"))
{
    if (solidModelDict().found("jacobianTangent"))
    {
        FatalIOErrorInFunction(solidModelDict())
            << "Both 'approximateJacobian' and 'jacobianTangent' are set. "
            << "'approximateJacobian' is deprecated: use 'jacobianTangent' "
            << "only." << exit(FatalIOError);
    }

    WarningInFunction
        << "'approximateJacobian' is deprecated. Use "
        << "'jacobianTangent scalar;' (was yes) or "
        << "'jacobianTangent fourthOrder;' (was no)." << endl;
    ...
}
```

### 3.3 Should anything ever be additively combined? (question c)

**Decision: coefficients are never summed. Operators may be summed in principle,
but are not summed today and should not start being summed in this work.**

The premise that there are "two Jacobians to combine" is a category error, and
the code proves it three ways:

1. **The stabilisation model is an operator factory, not a coefficient
   provider.** Every implementation of `vectorJacobian` is the same expression
   with the caller's `gamma`:
   `laplacianStab.C:172-201`, `diffStencilLaplacianStab.C:194-233`,
   `alphaStab.C:~228`, `JamesonSchmidtTurkelStab.C:~187`,
   `generalisedEvenOrderLaplacianStab.C:~202` all reduce to
   `scaleFactorJacobian()*fvm::laplacian(*gammaPtr, field, schemeName)`.
   `volStrainRateStab.C:237-256` is the only one that scales differently
   (`scaleFactorJacobian*tau/deltaT`), and even it takes `gamma` from the
   caller. There is no stabilisation-specific stiffness to add.

2. **`approxJ` is a surrogate for the whole residual, not for the stabilisation
   term.** `linGeomTotalDispSolid.C:1227-1231` builds
   `momentumStabilisation().vectorJacobian(D, &impKf_) - rho()*fvm::d2dt2(D)`
   and inserts it as the *entire* displacement block. The linearisation of
   `div(sigma)` is not present anywhere else. The name is misleading: it is a
   compact Laplacian obtained through the stabilisation model's factory.

3. **The framework already distinguishes summing residuals from selecting
   Jacobians.** `multipleStabilisationModels::updateVector`
   (`multipleStabilisationModels.C:127-167`) **sums** every sub-model's
   `faceVector()`. `multipleStabilisationModels::vectorJacobian`
   (`:182-191`) **picks one** by `jacobianModelName_`:

   ```c++
   const label idx = findIndex(modelNames_, jacobianModelName_);
   return models_[idx].vectorJacobian(field, gammaPtr, rebuild);
   ```

   So `multipleStabilisationModels` does not change the conclusion; it confirms
   it. Residual contributions are additive; the Jacobian is a single chosen
   operator.

The one place where addition *would* be mathematically correct: the residual is
`R(D) = div(sigma(D)) + s(D) - rho*d2dt2(D)`, so an exact Jacobian is
`d(div sigma)/dD + ds/dD - rho*d(d2dt2)/dD`. When `jacobianTangent fourthOrder`
assembles `d(div sigma)/dD` properly, `ds/dD` is genuinely missing. That
omission is `O(scaleFactor)` = `O(0.1)` of a Laplacian and is a legitimate
inexact-Newton choice. **Do not add a switch for it now** — it would invite the
double-counting it is meant to prevent (adding `ds/dD` *and* keeping `impKf_` as
the material coefficient is exactly the double count). Record the omission in
the header comment of the new Jacobian path and revisit only if SNES iteration
counts regress. See OQ-6.

What is *always wrong*: summing `impKf_` (a scalar stiffness) into a `mat66`
material tangent, or summing a scalar and a fourth-order tangent into one
coefficient. Neither is proposed.

### 3.4 API asymmetry and the topology rule (question e)

**Decision: `fourthOrder` is *not* made uniformly available. It is available
exactly on topologies whose integration points are flux-evaluation points and
whose integration-point-to-material map is a function.**

Encode the rule rather than documenting it:

```c++
// integrationPointTopology.H, new virtual
        //- Can this topology carry a fourth-order (mat66) material tangent?
        //  Two conditions must hold:
        //    1. the integration points are the locations at which a Jacobian
        //       operator evaluates fluxes (i.e. face-like, not cell-like), and
        //    2. each integration point belongs to exactly one mechanical
        //       constitutive law, so that no collapse is required.
        //  Condition 2 implies requiresUniqueIntegrationPointsPerMaterial()
        //  must be false. There is no meaningful collapse of two fourth-order
        //  tangents onto a shared integration point: the correct treatment of
        //  a material interface is a normal-direction traction/displacement
        //  continuity problem, not an average.
        virtual bool supportsFourthOrderTangent() const
        {
            return false;
        }
```

Applying the rule to the existing topologies:

<!-- markdownlint-disable MD013 -->

| Topology | `supportsFourthOrderTangent()` | Why |
|---|---|---|
| `cellCentredIntegrationPointTopology` | `false` | The cell-centred Jacobian is always operator-based (`fvm::laplacian` with a face-interpolated coefficient, or `hofvm` with per-face `gamma`). Fluxes are never evaluated at cell centres. **Your suspicion is correct: the cell-centred solvers genuinely never need a `mat66`.** |
| `faceCentredIntegrationPointTopology` | `false` **on the primary mesh** | `cellToIntegrationPointIDs_` appends each internal face to *both* owner and neighbour (`faceCentredIntegrationPointTopology.H:74-79`) and `requiresUniqueIntegrationPointsPerMaterial()` returns `true` (`:108-110`). A face on a material interface has two laws. Current `notImplemented` at `mechanicalConstitutiveLawManager.C:993-1000` becomes a structural invariant instead of a TODO. |
| `pointCentredIntegrationPointTopology` | `false` | Same sharing problem (`requiresUnique... == true`), and no operator consumes a point-based `mat66`. |
| `compactCellIntegrationPointTopology` | `false` | Cell-interior quadrature points; `hofvm` evaluates on *face* quadrature points. |
| **`dualFaceIntegrationPointTopology`** (new, §3.6) | `true` | Each dual face lies in exactly one primary cell, hence exactly one law. |
| `compactFaceIntegrationPointTopology` (future) | `true` when built | Face quadrature points; needed if `hofvm` moves from per-face to per-quadrature-point `gamma`. |

<!-- markdownlint-enable MD013 -->

The manager should assert this once, at the top of any update that is passed a
fourth-order request:

```c++
if (needsFourthOrderTangent(tangentReq) && !topo.supportsFourthOrderTangent())
{
    FatalErrorInFunction
        << "A fourth-order tangent was requested but topology "
        << topo.type() << " does not support one." << nl
        << "Fourth-order tangents require a topology whose integration "
        << "points each belong to exactly one mechanical constitutive law."
        << exit(FatalError);
}
```

**Keep `List<mat66>*` rather than a `GeometricField`.** `mat66` deliberately has
no `pTraits` and no `tmp<>` support (`mat66.H:30-33`: "a portable replacement
for `scalarSquareMatrix` ... avoids the need for `pTraits` and `tmp<>`"), so
there is no `surfaceMat66Field` to pass. The asymmetry in the signature is a
consequence of a deliberate portability choice and should stay. What *should*
change is that the flat-list primitive (§3.7) takes `UList<scalar>*` and
`UList<mat66>*` symmetrically, so the asymmetry lives only in the
convenience overloads.

### 3.5 The `highOrderJacobian` back-derivation (question d)

`linGeomTotalDispSolid.C:1170-1206` (and identically
`nonLinGeomTotalLagTotalDispSolid.C:1312-1318`,
`nonLinGeomUpdatedLagSolid.C:1249-1255`) does:

```c++
tmp<volScalarField> tK = mechanical().bulkModulus();
tmp<volScalarField> tMu = (impK_ - K)*(3.0/4.0);
tmp<volScalarField> tLambda = impK_ - 2.0*mu;
```

This inverts `impK = (4/3)*mu_eff + K` to recover the isotropic tangent pair
that `impK` and `bulkModulus()` encode between them.

**Correction, from a later reading.** An earlier version of this section said the
inversion is wrong for `linearElasticMisesPlastic` once a point yields. It is
not. That law returns `impK = scaleFactor*(4/3)*mu_ + K_` and
`bulkModulus() = K_`, so the back-derivation gives
`mu_back = scaleFactor*mu_`, exactly the softened shear modulus, and
`lambda_back = K_ - (2/3)*scaleFactor*mu_`, exactly the matching lambda for
deviatoric softening at unchanged bulk modulus - which is the right structure
for J2 plasticity. `linearElastic` and `neoHookeanElastic` define `impK` in the
same form, so all three are inverted exactly.

The real objection is structural rather than numerical. The solid model assumes
every law's `impK()` has that algebraic form and that `bulkModulus()` is the
matching `K`. `impK` is explicitly the "whatever converges best" coefficient, so
a law is entitled to return something else, and nothing declares, checks or
documents the assumption. A law that exercises that freedom gets a wrong pair
with no diagnostic.

Where an isotropic pair cannot represent the tangent at all - the orthotropic
and anisotropic hyperelastic laws - no way of obtaining `mu` and `lambda` helps.
`orthotropicLinearElastic` is protected only by an accident: it implements
`impK` but not `bulkModulus`, which is `notImplemented` in the base class, so
the back-derivation fatal-errors there rather than computing nonsense.

**Note also that the obvious repair is a trap.** `shearModulus()` exists, but
`linearElasticMisesPlastic::shearModulus()` returns the elastic `mu_` with no
`scaleFactor`, so substituting it would discard the softening the
back-derivation captures. It is also `notImplemented` in the base class and
overridden by only seven of about twenty-five laws. Adding a
`tangentShearModulus()` would re-encode the same implicit contract under a new
name and still could not express an anisotropic tangent.

**Decision: replace the three separate `hofvm` calls with one `mat66`-driven
call. No new mathematics is required — the kernel already exists.**

The three coefficient kernels
(`hofvm.C:123-172`: `laplacianCoeff`, `laplacianTransposeCoeff`,
`laplacianTraceCoeff`) together compute, for isotropic `C`,

```text
coeff_ij = w*( mu*(n.g)*delta_ij + mu*g_i*n_j + lambda*n_i*g_j )
```

which is exactly `Sf_m C_mikl g_k delta_lj` with `Sf = w*n`. That expression is
already implemented, and already used by the vertex-centred solver, as

```c++
    // numerics/multiplyCoeff/multiplyCoeff.H:54-60
    //     coeff_ij = Sf_m C_mikl g_k delta_lj
    void multiplyCoeff
    (
        tensor& coeff,
        const vector& Sf,
        const mat66& C,
        const vector& g
    );
```

So the change inside `hofvm.C` is a one-line substitution in the assembly loop
(`hofvm.C:236-262`):

```c++
    // was:
    //     const tensor coeff = calcCoeff
    //     (gammaFace, quadPointW, gradCoeffs[faceI][qpI][cI], faceNormal);
    tensor coeff;
    multiplyCoeff
    (
        coeff,
        quadPointW*faceNormal,              // Sf
        materialTangent[faceI],             // C
        gradCoeffs[faceI][qpI][cI]          // g
    );
```

New public signature:

```c++
    //- Add the coefficients of the linearisation of div(sigma) to the PETSc
    //  matrix using a full fourth-order material tangent.
    //  materialTangent must be indexed by mesh face index and sized
    //  mesh.nFaces(): internal faces first, then boundary faces in patch
    //  order, matching faceCentredIntegrationPointTopology.
    //  This replaces the combination of laplacianIntoPETScMatrix,
    //  laplacianTransposeIntoPETScMatrix and laplacianTraceIntoPETScMatrix,
    //  which are only valid for isotropic linear elasticity.
    void divSigmaIntoPETScMatrix
    (
        Mat jac,
        const foamPetscSnesHelper& petscSnesHelper,
        const movingLeastSquares& mls,
        const volVectorField& D,
        const List<mat66>& materialTangent,
        const label rowOffset = 0,
        const label colOffset = 0
    );
```

**Do the `hofvm` operators need `mat66` variants? Yes, one, and only one.**
The three scalar entry points can then be deleted once no caller remains.

**Cost.** `mat66` is 36 `scalar` = 288 B. One per mesh face (internal +
boundary), not per quadrature point, because `hofvm` currently uses a single
face value of `gamma` (`hofvm.C:215-218`). For a 1M-cell hex mesh
(~3M faces) that is ~865 MB, against ~24 MB for the scalar field. That is a
real cost and it is the reason `jacobianTangent` must default to `scalar`.
Two mitigations, in order of preference:

1. Do not materialise the list at all for laws whose tangent is spatially
   uniform. Add `virtual bool hasUniformFourthOrderTangent() const` to
   `mechanicalConstitutiveLaw` (true for `linearElastic`,
   `neoHookeanElastic` in the small-strain limit) and a manager query that
   returns a single `mat66`; `hofvm` then takes `const mat66&` in a second
   overload. This covers `plateHole` and `wobblyNewton` at zero memory cost.
2. Only allocate the `List<mat66>` when `jacobianTangent` is `fourthOrder*`,
   and free it between Jacobian evaluations if `snes_lag_jacobian` is large.

Note that once `hofvm` interpolates `gamma` to quadrature points ("If Gamma is
not constant, next step is to interpolate diffusivity to quad points using
`hofvc::interpolate`", `hofvm.C:215-216`), the storage multiplies by the number
of face quadrature points. Plan for the uniform-tangent fast path now.

**Result impact.** For `linearElastic` the back-derived `mu`/`lambda` and the
true `mat66` produce *identical* matrix coefficients, so `plateHole` and
`wobblyNewton` are bit-for-bit. For other laws the Jacobian changes, which by
§1.2 changes only the iteration count, not the converged solution. Report both.

### 3.6 The vertex-centred path (question f)

This was flagged as the hardest sub-problem. **It is not, and the reason is
worth stating plainly: there is no multi-material dual face.**

**Evidence.** `dualMeshToMeshMap.C:102-107`:

```c
// Set dualFaceToCell
// All internal dual faces should uniquely lie within one primary cell. For
// boundary faces, each dual face will uniquely lie on the boundary of one
// cell.
```

`dualFaceToCell_` is a `labelList` (`dualMeshToMeshMap.H:95`), i.e. a
single-valued function of dual face index, assigned once per dual face
(`dualMeshToMeshMap.C:169-186, 236`). Materials are defined per primary cell
(`mechanicalConstitutiveLawManager.C:557-571`). Composition of two functions is
a function: **dual face -> primary cell -> mechanical law is single-valued.**
`dualMechanicalModel::impKf()` already relies on this
(`dualMechanicalModel.C:396-397`: "Each face in the dualMesh lies in one cell
(and hence one material) in the primary mesh") and handles the multi-material
case correctly. The `notImplemented` at `dualMechanicalModel.C:356` is a
limitation of the sub-mesh plumbing, not of the geometry — as its own comment
admits ("This does work but the problem can be that the maps for the boundary
faces sometimes struggle to be defined").

**So: no collapse rule is required on the dual path, and removing sub-meshes
removes the reason the multi-material case was `notImplemented`.**

**Indexing contract.** `vfvmCellPoint.C:87-93`:

```c++
forAll(dualOwn, dualFaceI)
{
    const label cellID = dualFaceToCell[dualFaceI];
    const mat66& materialTangent = materialTangentField[dualFaceI];
```

The `List<mat66>` is indexed by **dual mesh face index** and sized
`dualMesh().nFaces()` (`vertexCentredNonLinGeomTotalLagSolid.C:1684`; note
`vertexCentredLinGeomSolid.C:1331` sizes it with the *primary* `mesh.nFaces()`
and gets away with it only because `materialTangentFaceField` resizes — fix
that while you are there).

**Replacement design.**

New topology, in `integrationPointTopology/`:

```c++
/*---------------------------------------------------------------------------*\
              Class dualFaceIntegrationPointTopology Declaration
\*---------------------------------------------------------------------------*/

//- Integration points located at the faces of a dual mesh.
//  Used by the vertex-centred solid models, where stresses and material
//  tangents are evaluated at dual-mesh faces while materials remain defined
//  per primary-mesh cell.
//  Each dual face lies within (or on the boundary of) exactly one primary
//  cell, so each integration point belongs to exactly one mechanical
//  constitutive law and no stress or tangent collapse is required.
class dualFaceIntegrationPointTopology
:
    public integrationPointTopology
{
    // Private Data

        //- Number of dual mesh faces
        const label nDualFaces_;

        //- Dual faces associated with each primary cell.
        //  This is the inverse of the dualFaceToCell map.
        labelListList cellToIntegrationPointIDs_;

public:

    TypeName("dualFaceIntegrationPointTopology");

        //- Construct from the primary mesh and the dual-face-to-cell map
        dualFaceIntegrationPointTopology
        (
            const fvMesh& mesh,
            const labelList& dualFaceToCell
        );

        virtual label nIntegrationPoints() const override
        {
            return nDualFaces_;
        }

        virtual const labelUList cellIntegrationPointIDs
        (
            const label cellI
        ) const override
        {
            return cellToIntegrationPointIDs_[cellI];
        }

        //- Dual-face boundary states are not required: boundary dual faces
        //  are already covered by the owning primary cell
        virtual bool boundaryAware() const override
        {
            return false;
        }

        //- Each dual face belongs to exactly one primary cell
        virtual bool requiresUniqueIntegrationPointsPerMaterial() const override
        {
            return false;
        }

        //- Supported: the map to materials is single-valued
        virtual bool supportsFourthOrderTangent() const override
        {
            return true;
        }
};
```

It cannot be built by the existing run-time selection table, which takes only
`(const fvMesh&)` (`integrationPointTopology.H:81-88`). Hence:

```c++
// mechanicalConstitutiveLawManager.H
        //- Register an externally constructed integration-point topology.
        //  Required for topologies that cannot be built from the mesh alone,
        //  e.g. dualFaceIntegrationPointTopology, which needs the
        //  dual-face-to-cell map.
        //  The manager takes ownership; the key must be unique.
        const integrationPointTopology& registerTopology
        (
            const word& key,
            autoPtr<integrationPointTopology> topoPtr
        ) const;
```

Everything downstream already works unchanged:
`mechanicalConstitutiveLawManager.C:137-269` builds
`lawIntegrationPointIDs_`, `states_` and (skipped here) `boundaryStates_`
purely from the virtual interface.

**The values live on the dual mesh, not on `mesh_`.** `dualGradDf_` and
`dualSigmaf_` are `surfaceTensorField`/`surfaceSymmTensorField` on
`dualMesh()` (`vertexCentredLinGeomSolid.H:150,153`), so every existing
`GeometricField` overload would trip `checkMeshConsistency`
(`mechanicalConstitutiveLawManager.H:191-207`). This forces the flat-list
primitive of §3.7 — which is the same addition that removes the API asymmetry
of question (e). One change, two problems solved.

**Geometric stiffness stays solid-model-owned.**
`vertexCentredNonLinGeomTotalLagSolid.C:1688-1691` builds a
`List<mat39> geometricStiffness` from `dualSf` and `dualGradDf_`. That is a
function of kinematics and the current stress, not of the constitutive law, and
it must not be pushed behind `tangentRequest`. See OQ-7.

### 3.7 The flat-list primitive

```c++
        //- Update the engineering (small-strain) stress at the integration
        //  points of the supplied topology.
        //  This is the primitive form of the small-strain update: every
        //  GeometricField and CompactListList overload is implemented in
        //  terms of it.
        //  The caller owns all storage and is responsible for its indexing,
        //  which must match the integration-point indices returned by the
        //  topology. The lists need not be associated with the manager's mesh,
        //  which allows evaluation on a dual or otherwise derived mesh.
        //  No stress or tangent collapse is performed: the topology must
        //  return requiresUniqueIntegrationPointsPerMaterial() == false.
        //  A fourth-order tangent request additionally requires
        //  topo.supportsFourthOrderTangent().
        void updateStressSmallStrain
        (
            const integrationPointTopology& topo,
            const UList<tensor>& gradD,
            const UList<tensor>& gradD0,
            const scalar dt,
            UList<symmTensor>& stress,
            UList<scalar>* scalarTangentPtr = nullptr,
            UList<mat66>* fourthOrderTangentPtr = nullptr,
            const tangentRequest tangentReq = tangentRequest::none
        );

        //- Finite-strain counterpart of the above
        void updateStressFiniteStrain
        (
            const integrationPointTopology& topo,
            const UList<tensor>& F,
            const UList<tensor>& F0,
            const UList<tensor>& Finv,
            const UList<tensor>& Finv0,
            const UList<scalar>& J,
            const UList<scalar>& J0,
            const scalar dt,
            UList<symmTensor>& stress,
            UList<scalar>* scalarTangentPtr = nullptr,
            UList<mat66>* fourthOrderTangentPtr = nullptr,
            const tangentRequest tangentReq = tangentRequest::none
        );
```

This is a pure refactor of what
`mechanicalConstitutiveLawManager.C:1002-1063` already does: build
`UIndirectList` views per law, wrap in kinematics and a response, call
`evaluate`. It introduces no temporaries, keeps the loop-over-integration-points
style, and is the natural kernel boundary for later GPU work. The existing
`GeometricField` overloads keep their signatures and add the collapse and
boundary handling on top.

### 3.8 `solvePressure` / mixed u-p (question g)

`linGeomTotalDispSolid.C:571-579` becomes:

```c++
    impK_
    (
        IOobject("impK", ...),   // Constraint C1: must keep this name
        mesh(),
        dimensionedScalar("0", dimPressure, 0.0)
    ),
    impKf_(...),
    rImpK_(...),
    rKappa_(1.0/mechanicalManager().kappa()),
```

with, in the constructor body:

```c++
    // Build the implicit stiffness. A deviatoric-only tangent is used in the
    // mixed displacement-pressure formulation, where the hydrostatic response
    // is carried by the pressure equation.
    mechanicalManager().updateScalarTangent
    (
        gradD(),
        gradD().oldTime(),
        0.0,
        impK_,
        solvePressure()
      ? tangentRequest::scalarDeviatoric
      : tangentRequest::scalar
    );
```

* `2.0*mechanical().shearModulus()` -> `tangentRequest::scalarDeviatoric`.
  Bit-for-bit for `linearElastic` (`2*mu` both sides, §1.6).
* `1.0/mechanical().bulkModulus()` -> `1.0/mechanicalManager().kappa()`
  (`mechanicalConstitutiveLawManager.H:352-356`, implemented at `.C:683-724`).
* `mechanical().shearModulus()` and `mechanical().bulkModulus()` then have no
  remaining consumer in this solid model, which is the point.

**Should the law define the hydrostatic response? Yes, and it must eventually —
but not in this change.** Today the solid model hard-codes the constitutive
split:

* `linGeomTotalDispSolid.C:963`: `sigma() = dev(sigma()) - p*I;`
* `linGeomTotalDispSolid.C:991-996`: `-p*rKappa_ + stab - tr(gradD())`

i.e. `p = -kappa*tr(eps)` with an additive volumetric/deviatoric split. That is
correct for linear elasticity and wrong for, e.g., a neo-Hookean law where
`p = kappa*ln(J)/J`. `tangentRequest::scalarDeviatoric` gets the *tangent* side
of the split right; the *stress* side stays hard-coded.

The design must not block the fix. The natural future addition, deliberately
**out of scope** here, is symmetric with the existing tangent request:

```c++
        //- Evaluate the hydrostatic (pressure) response and, optionally, its
        //  tangent. Only meaningful for laws with a well-defined
        //  volumetric-deviatoric split.
        virtual void evaluateHydrostatic
        (
            const smallStrainMechanicalConstitutiveLawKinematics& kin,
            mechanicalConstitutiveLawState& state,
            UIndirectList<scalar>& p,
            UIndirectList<scalar>* dpdVolStrainPtr = nullptr
        ) const;
```

Keeping `p` solid-model-owned now, and routing only the tangent through the
manager, leaves that door open without a second migration.

---

## 4. Rejected alternatives

**A. The manager owns and caches a "current tangent" field.** Rejected because
`impK` is not one quantity: `linGeomTotalDispSolid` wants `scalarDeviatoric`
when `solvePressure()` and `scalar` otherwise (`linGeomTotalDispSolid.C:573`),
`thermalLinGeomSolid` wants a frozen cell field *and* a periodically refreshed
face field with different values (`thermalLinGeomSolid.C:325` vs `:364`), and
the vertex-centred solvers want it on a different mesh entirely. A cached
singleton would have to be keyed by request type, mesh and refresh epoch, which
is a worse version of letting the caller own it. It also breaks the one-way
dependency the README establishes (lines 60-74) by forcing the manager to know
about interpolation schemes.

**B. One unified enum that also absorbs `highOrderJacobian`, e.g.
`jacobian compactScalar | compactFourthOrder | highOrderScalar |
highOrderFourthOrder`.** Rejected because it is a Cartesian product of two
genuinely orthogonal choices (which discrete operator; which material stiffness)
and would grow multiplicatively with every new operator. It would also move a
discretisation switch out of `highOrderCoeffs`, where it sits next to its
sibling `highOrderResidual` and its `PETScSNES` guard
(`solidModel.H:586-600`), into the material-facing part of the dictionary.

**C. Keep `approximateJacobian` alongside a new `jacobianTangent`.** Rejected:
they are the same axis with different spellings, and two spellings for one
concept is how `stabilisation/type` vs `stabilisation/momentum/type` became a
fatal-error migration (`solidModel.C:1346-1352`). Deprecate with a shim and
delete after one release.

**D. Define a collapse rule for `mat66` at multi-material faces (arithmetic or
harmonic mean of the two tangents).** Rejected as physically wrong.
Continuity at a bi-material interface is a normal-direction traction and
displacement matching problem; the correct "combined" stiffness is a
Schur complement of the two tangents restricted to the interface normal, not a
componentwise average of two 6x6 matrices. An arithmetic mean of tangents does
not even preserve positive-definiteness of the resulting operator in the
presence of strong contrast. Since the only real `mat66` consumer (the dual
path) is single-valued by construction (§3.6), pay nothing for a rule nobody
needs, and make the unsupported combination a structural invariant
(`supportsFourthOrderTangent()`) rather than a runtime `notImplemented`.

**E. Give `mat66` a `GeometricField` specialisation so all overloads look the
same.** Rejected: `mat66.H:30-33` is explicit that avoiding `pTraits` and
`tmp<>` is the point, precisely so the type works across foam-extend,
OpenFOAM.org and OpenFOAM.com without three sets of `#ifdef`s. Fixing a
cosmetic signature asymmetry by taking on that portability burden is a bad
trade. The flat-list primitive (§3.7) makes the asymmetry disappear where it
matters.

**F. Rewrite all solid models onto the manager in one change.** Rejected on the
evidence of §2: `impK` has 220 occurrences across 12 solid models with at least
four different lifecycles and three different stabilisation couplings. A single
change could not be regression-tested per solver. The staged plan below lets
one solid model at a time move while `mechanicalModel` remains the default.

**G. Make the `impKf_` refresh cadence a dictionary entry now.** Rejected for
this change: in `linGeomTotalDispSolid` the cadence is not a tuning parameter
but a correctness constraint (§1.3, `impKf_` scales a residual term), and
exposing it would invite users to change converged answers while believing they
were tuning convergence. Revisit once the stabilisation coefficient is
decoupled from `impKf_`. See OQ-3.

---

## 5. Staged implementation plan

Each stage is independently mergeable and independently revertible.
Every stage must build on OpenFOAM.com (v2512), OpenFOAM.org and foam-extend,
with PETSc code under `#ifdef USE_PETSC` and version-specific code under the
existing `OPENFOAM_COM` / `OPENFOAM_ORG` / `OPENFOAM_NOT_EXTEND` / `FOAMEXTEND`
guards, as the current files do.

### Stage 0 — framework only; nothing consumes it

* `integrationPointTopology::supportsFourthOrderTangent()` (default `false`).
* `mechanicalConstitutiveLawManager::registerTopology(word, autoPtr<...>)`.
* Flat-list `updateStressSmallStrain`/`updateStressFiniteStrain` primitives
  (§3.7); re-express the existing `GeometricField` and `CompactListList`
  overloads on top of them.
* `mechanicalConstitutiveLawManager::updateScalarTangent(...)` (§3.1).
* Guard fourth-order requests against `supportsFourthOrderTangent()`.
* Fix two defects found while reading (§6).

**Checkpoint 0.** `wmake` on all three forks. Run `plateHole`,
`wobblyNewton`, `perforatedPlate`. All three must be **byte-identical** to
`development`, since no solid model references any of the new code. Diff the
final time directories with `foamListTimes` + `diff -r`.

### Stage 1 — selection API in `solidModel`, unused

* `solidModel::jacobianTangent(const tangentRequest deflt) const`.
* Deprecation shim reading `approximateJacobian` (§3.2).

**Checkpoint 1.** Same three tutorials, plus
`tutorials/solids/linearElasticity/cantilever2d` with
`solidProperties.vertexCentred` (which sets `approximateJacobian no` at line
28). All identical; the shim must emit exactly one warning and select
`fourthOrder`.

### Stage 2 — `linGeomTotalDispSolid` sources `impK` from the manager

* Construct a `mechanicalConstitutiveLawManager` alongside the existing
  `mechanicalModel`, behind
  `solidModelDict().lookupOrDefault<Switch>
  ("useMechanicalConstitutiveLawManager", false)`.
* When on: fill `impK_` via `updateScalarTangent`, with
  `solvePressure() ? scalarDeviatoric : scalar` (§3.8); keep `impKf_` and
  `rImpK_` derived exactly as now; keep `impK_` registered under `"impK"`
  (Constraint C1).
* `rKappa_` from `mechanicalManager().kappa()`.
* Stress still comes from `mechanicalModel::correct()` in this stage.

**Checkpoint 2.** Run `plateHole` and `perforatedPlate` twice, switch off then
on. Off must be byte-identical to Stage 1. On must agree to solver tolerance,
and the `impK` field itself must agree to `1e-12` relative — write it out and
compare directly, because §1.3 means any `impK` difference propagates into the
converged answer. `perforatedPlate` is the discriminating case: it is the only
one of the three where `impK()` is state-dependent
(`linearElasticMisesPlastic.C:739-777`). Also run one restart to expose OQ-1.

### Stage 3 — `fourthOrder` on the high-order Jacobian path

* `hofvm::divSigmaIntoPETScMatrix(..., const List<mat66>&)` plus the
  uniform-tangent overload taking `const mat66&` (§3.5).
* `linGeomTotalDispSolid::formJacobian` selects on
  `jacobianTangent(tangentRequest::scalar)`; `scalar` keeps the back-derivation
  (deprecated), `fourthOrder` uses the new call.
* Same for `nonLinGeomTotalLagTotalDispSolid` and
  `nonLinGeomUpdatedLagSolid`.

**Checkpoint 3.** `wobblyNewton` (PETScSNES, `petscOptions.mf.hypre`) plus any
case setting `highOrderCoeffs/highOrderJacobian yes`. Converged fields must
match to solver tolerance in both settings and be *identical* for
`linearElastic` (§3.5). Log and compare SNES iteration counts separately —
they are expected to change for non-`linearElastic` laws and that is the point.

### Stage 4 — vertex-centred solvers off `dualMechanicalModel`

* `dualFaceIntegrationPointTopology` (§3.6), registered via
  `registerTopology`.
* `vertexCentredLinGeomSolid` and `vertexCentredNonLinGeomTotalLagSolid` call
  the flat-list update on dual faces for stress, scalar tangent
  (replacing `dualImpKf()`) and `mat66` (replacing
  `dualMechanicalPtr_().materialTangentFaceField(...)`).
* `mat39` geometric stiffness unchanged, still solid-model-owned.
* Delete `dualMechanicalModel` from these two solvers.

**Checkpoint 4.** `cantilever2d` with `solidProperties.vertexCentred`, run with
both `jacobianTangent scalar` and `jacobianTangent fourthOrder`. Add a
**new** two-material vertex-centred regression case — this stage is the first
time multi-material `mat66` is possible at all
(`dualMechanicalModel.C:356` is `notImplemented` today), so it needs coverage
that does not exist yet.

### Stage 5 — remaining cell-centred solid models, one PR each

Order by risk, easiest first, because the stabilisation coupling determines how
tightly `impK` is pinned:

1. `unsLinGeomSolid`, `unsNonLinGeomTotalLagSolid`,
   `unsNonLinGeomUpdatedLagSolid` — **no stabilisation term**, so `impKf_` is
   genuinely free (§2.2); lowest risk.
2. `thermalLinGeomSolid`, `poroLinGeomSolid` — already split frozen `impK_`
   from refreshed `impKf_`; preserve the `iCorr % 10` cadence verbatim.
3. `nonLinGeomTotalLagTotalDispSolid`, `nonLinGeomUpdatedLagSolid` — same
   stabilisation coupling as `linGeomTotalDispSolid`.
4. `coupledPressureDisplacementSolid` — most `impK` sites (48 across seven
   files), block-coupled assembly, and a refresh in `updateTotalFields()`.

**Checkpoint 5.** Per solver: its own tutorials, on/off comparison as in
Checkpoint 2.

### Stage 6 — retire the legacy interface

Only once no solid model uses it: remove `mechanicalModel::impK/impKf/
bulkModulus/shearModulus` and `mechanicalLaw::impK/impKf/materialTangent/
materialTangentField`. Migrate the six registry consumers
(`standardPenalty`, `pointStandardPenalty`, `standardPenaltyFriction`,
`elasticWallPressure`, the three cohesive zone models, `solidSubMeshes`) to an
explicit query, or keep `"impK"` registered for good and document it as public
API. See OQ-8.

---

## 6. Defects found while reading (fix in Stage 0)

**D1. `stressCollapseRule::average` returns the reciprocal of the tangent.**
`mechanicalConstitutiveLawManager.C:1076-1079` accumulates `sum(K)`, then
`:1207-1210` returns `w/sum(K)`. For a single contribution (`w == 1`) that is
`1/K`, with units of `1/Pa`. The `harmonic` branch accumulates `sum(1/K)` and
returns `w/sum(1/K)`, which is correct. The `average` branch should return
`sum(K)/w`. Latent today because no solid model calls this overload with a
scalar tangent.

**D2. `scalarDeviatoric` is self-inconsistent in the plastic law.**
`linearElasticMisesPlasticMechanicalConstitutiveLaw.C:302-306` returns `2*mu`
at an elastic point (with `(4/3)*mu` commented out immediately above), while
`:410-412` returns `theta*(4/3)*mu` at a plastic point. At `theta == 1` these
differ by a factor of `3/2`, so the tangent is discontinuous across first yield.
`2*mu` matches the legacy `2.0*mechanical().shearModulus()`
(`linGeomTotalDispSolid.C:574`), so the elastic branch is the one to keep and
the plastic branch should be `theta*2*mu`. Must be resolved before
`scalarDeviatoric` is used with plasticity (see OQ-4).

**D3. `vertexCentredLinGeomSolid.C:1331` sizes the `mat66` list with the
primary mesh's `nFaces()`** while its non-linear sibling
(`vertexCentredNonLinGeomTotalLagSolid.C:1684`) correctly uses
`dualMesh().nFaces()`. Harmless only because `materialTangentFaceField`
resizes internally (`dualMechanicalModel.C:342`). Removed by Stage 4 anyway,
but worth a one-line fix now so the two solvers read the same.

**D4. `finiteDifferenceFourthOrder` uses an absolute perturbation.**
`mechanicalConstitutiveLaw.C:56`: `const scalar h = 1e-8;  // configurable
later`, applied directly to `gradD` components regardless of their magnitude.
For strains below ~1e-6 this is not a small perturbation, and for strains above
~1e-2 it is close to round-off relative to the stress. It should be scaled, e.g.
`h = max(1e-8, 1e-6*max(mag(gradD), SMALL))`. Not blocking, but do not
advertise `fourthOrderFiniteDifference` as production-ready until it is fixed
(see OQ-5).

**D6. The two cell-centred topologies declare each other's type name.**
`cellCentredIntegrationPointTopology.H:68` declares
`TypeName("compactCellIntegrationPointTopology")` and
`compactCellIntegrationPointTopology.H:69` declares
`TypeName("cellCentredIntegrationPointTopology")`. Nothing misbehaves, because
the run-time selection entry and the manager's lookup both go through
`cellCentredIntegrationPointTopology::typeName` and so agree with each other.
But `topo.type()` names the wrong class in every diagnostic, including the
tangent guards of PR-1, and both classes register the same debug switch. Found
and fixed in PR-2.

**D5. `README.md` is stale.** Lines 139-144 list the supported tangent modes as
`none`, `scalar`, `scalarDeviatoric` only, and lines 283-286 list face- and
vertex-based topologies as "planned" although
`faceCentredIntegrationPointTopology` and
`pointCentredIntegrationPointTopology` both exist. Update alongside Stage 0 and
replace lines 234-253 with a pointer to this document.

---

## 7. Open questions

Each is phrased so it can be answered in one line.

**OQ-1.** On restart, `impK_` is built from a state-dependent
`linearElasticMisesPlastic::impK()` that reads `DLambda_` from disk, so a
restarted run has a different `impK_` from a continuous one. Reproduce that
restart-dependence, or always build `impK_` from the elastic tangent?

**OQ-2.** `mechanicalConstitutiveLawState` has no `IOobject` and is not
written. Should it become IO-aware before Stage 2, or is loss of constitutive
history across restart acceptable while the manager runs alongside
`mechanicalModel`?

**OQ-3.** Keep the `impKf_` refresh cadence hard-coded per solid model (my
recommendation, §3.1), or expose one `impKUpdateFrequency` entry defaulting to
each solver's current behaviour?

**OQ-4.** Is `tangentRequest::scalarDeviatoric` intended to be `2*mu` (matching
`2.0*mechanical().shearModulus()`) or `(4/3)*mu`? See D2.

**OQ-5.** Is `tangentRequest::fourthOrderFiniteDifference` intended for
production, or only for verifying analytical tangents? It costs six extra law
evaluations per integration point
(`mechanicalConstitutiveLaw.C:61-73`).

**OQ-6.** Should the `fourthOrder` Jacobian also include
`d(stabilisation)/dD`, or is omitting it an acceptable inexact-Newton choice
(my recommendation, §3.3)?

**OQ-7.** Should the `mat39` geometric stiffness stay solid-model-owned (my
recommendation, §3.6), or move behind the manager as a
`tangentRequest::fourthOrderWithGeometric`?

**OQ-8.** Should a `volScalarField` named `"impK"` remain permanently
registered as public API for contact, cohesive-zone and FSI models, or should
those six consumers be migrated to an explicit query in Stage 6?

**OQ-9.** `linearElastic::impK()` has a `nu == 0.5` branch returning `2*mu_`
(`linearElastic.C:249`), but `linearElasticMechanicalConstitutiveLaw` rejects
`nu >= 0.5 - SMALL` outright (`:71-78`). Is the incompressible-limit case meant
to be reachable in the new framework, and if so via `scalarDeviatoric` plus
`solvePressure`?

**OQ-10.** `momentumStabilisation().vectorJacobian(D, &impKf_)` is documented as
"an approximate Jacobian of the stabilisation term"
(`stabilisationModel.H:327-330`) but is used as the Jacobian of the *entire*
momentum residual (`linGeomTotalDispSolid.C:1227-1231`). Is the comment wrong,
or is the usage? This determines whether §3.3's "operators may be summed"
caveat is live.

---

## 8. Decisions taken after review

Four questions were put to the author and answered. This section records the
answers, what they supersede, and what follows from them.

### 8.1 `planeStress`: the manager injects it into each law's sub-dictionary

**Decision.** The manager copies the top-level `planeStress` entry of
`mechanicalProperties` into a copy of each law's sub-dictionary before calling
`mechanicalConstitutiveLaw::New(dict)`. Laws that care read
`dict.lookupOrDefault<Switch>("planeStress", false)`.

**Rationale and scope.** `planeStress` is a *constructor-time parameter
conversion*, not an evaluation-time branch. In `linearElastic.C` it appears only
at lines 77 and 117, both in the constructor, where it selects between
`lambda_ = nu*E/((1+nu)*(1-nu))` and `lambda_ = nu*E/((1+nu)*(1-2*nu))`
(and correspondingly for `K_`); `correct()` never sees it. The same holds for
`linearElasticMisesPlastic.C:117-129`, which additionally forbids the
combination with `solvePressureEqn`. Thirteen legacy laws implement it.

This choice requires **no RTS signature change** (`New(dict)` is already shipped
in this PR) and **no case dictionary change**. It also avoids reintroducing the
registry coupling that legacy uses: `mechanicalLaw::planeStress()`
(`mechanicalLaw.C:524-545`) reaches into the object registry for the
`mechanicalProperties` `IOdictionary`, which is exactly the hidden global
dependency the new `(dict)`-only constructor was designed to remove, and which
would make the manager unusable on a mesh with no registered
`mechanicalProperties`.

**Implementation notes.**

```c++
// mechanicalConstitutiveLawManager constructor, replacing the direct
// lawEntries[lawI].dict() use at mechanicalConstitutiveLawManager.C:542-549
const Switch planeStress(dict.lookupOrDefault<Switch>("planeStress", false));

dictionary lawDict(lawEntries[lawI].dict());

if (lawDict.found("planeStress"))
{
    FatalIOErrorInFunction(lawDict)
        << "'planeStress' is set by the mechanicalProperties dictionary and "
        << "must not be given inside the '" << lawName << "' sub-dictionary."
        << exit(FatalIOError);
}

lawDict.add("planeStress", planeStress);

laws_.set(lawI, mechanicalConstitutiveLaw::New(lawDict));
```

**Known limitation to document, not to fix.** Legacy plane stress for a *plastic*
law is only a 2D parameter reduction: a true plane-stress return map enforces
`sigma_zz == 0` during the return, which `linearElasticMisesPlastic` does not
do. The new law must reproduce the legacy behaviour exactly, but the README
should not describe it as plane-stress plasticity.

Closes OQ-9 partially: the incompressible-limit branch remains open.

### 8.2 `scalarDeviatoric` is `(4/3)*mu` — supersedes D2, §3.8 and OQ-4

**Decision.** `tangentRequest::scalarDeviatoric` returns `(4/3)*mu` for an
elastic point and `theta*(4/3)*mu` for a plastic point. The plastic branch
(`linearElasticMisesPlasticMechanicalConstitutiveLaw.C:410-412`) is already
correct; the **elastic branch at `:302-306` is the one to fix**, by restoring the
`(4.0/3.0)*mu` that is currently commented out at `:305` and deleting the
`2.0*mu` at `:306`. `linearElasticMechanicalConstitutiveLaw.C:130-132` must
change from `2.0*mu_.value()` to `(4.0/3.0)*mu_.value()`.

This is the principled scalar Laplacian surrogate: for
`dev(sigma) = 2*mu*dev(eps)`,
`div(dev(sigma)) = mu*lap(D) + (1/3)*mu*grad(div(D))`, whose scalar surrogate
coefficient is `mu + (1/3)*mu = (4/3)*mu`.

**Consequence, and it is not small.** `linGeomTotalDispSolid.C:574` currently
uses `2.0*mechanical().shearModulus()` = `2*mu` for `impK_` when
`solvePressure()`. The new value differs by a factor of `3/2`, and by §1.3
`impKf_` scales a solution-changing stabilisation traction, so **converged
mixed u-p results will change**. This is a deliberate correction, not a
regression, but it means:

* the mixed u-p path cannot be validated by an on/off comparison — the two
  systems are *meant* to disagree there;
* `linGeomTotalDispSolid.C:574` should be changed to `(4.0/3.0)*shearModulus()`
  in the **same** commit that fixes the constitutive laws, so that legacy and
  new agree and the on/off check stays meaningful. Do this in PR-1, before any
  solid model adopts the manager;
* **no tutorial exercises this particular path.** `solvePressure` is a
  `solidProperties` entry read at `solidModel.C:1163-1165` with default
  `false`, and no case in `tutorials/` sets it.

  Two things that look like coverage are not. The `solvePressureEqn yes`
  entries in five tutorials are the unrelated *law-level* pressure equation of
  the legacy mechanical laws. And `plateHole/regressionTest.sh` does run four
  `pressureDisplacement{Compressible,Incompressible} {coarse,medium}` cases —
  but those select `solidModel coupledPressureDisplacementSolid`
  (`caseOptions/pressureDisplacement/hex/common/constant/solidProperties:20`),
  are foam-extend only, and that solid model never calls `shearModulus()`.
  So they do not touch `linGeomTotalDispSolid.C:574`.

  The line being changed was therefore genuinely untested. **Resolved in PR-1**
  by adding a `petscSnesPressure` approach to the plateHole tutorial and its
  regression test. The path runs and converges; measured against the
  plate-with-hole analytical solution, `(4/3)*mu` improves on `2*mu` in every
  metric at the same SNES iteration count:

  | metric | `(4/3)*mu` | `2*mu` |
  |---|---|---|
  | `DDifference` LInf | 4.22837e-08 | 4.30140e-08 |
  | `pointDDifference` LInf | 5.03405e-08 | 5.68787e-08 |
  | stress component-0 LInf | 71259.9 | 82738.6 |

  This also confirms the expected behaviour of the stabilisation term: it is a
  consistency error that shrinks under refinement, so scaling it differently
  moves the discrete answer slightly without changing the continuum limit.

Supersedes the "2*mu both sides, bit-for-bit" claims in §1.6 and §3.8, and
resolves D2 and OQ-4.

### 8.3 Migration granularity: thin slice, tangent only

**Decision.** The first adopting solid model takes the **tangent only**. Stress
continues to come from the legacy path. The manager runs alongside
`mechanicalModel`, constructed from the identical `mechanicalProperties`
dictionary, behind a per-case switch defaulting to off.

This is possible with no dictionary change at all because the new laws register
the **same type names** as the legacy ones (`linearElastic`,
`linearElasticMisesPlastic`, `neoHookeanElastic`) in a **separate** run-time
selection table, and the manager reads the same `mechanical (...)` list.

Consequences:

* `mechanicalConstitutiveLawState` does **not** need to become IO-aware yet,
  because `updateScalarTangent` evaluates against a state *copy* and commits
  nothing. Defers OQ-2 to the stage where stress moves.
* Only three laws exist in the new framework. Opting in with any other law must
  fail at construction with a message listing the supported types, not silently.
* Two law objects per material are constructed per case. Negligible: laws hold
  only `dimensionedScalar` parameters.

Resolves OQ-2 for the near term.

### 8.4 First mover: `vertexCentredLinGeomSolid`

Supersedes the Stage 2/4 order in §5.

**Decision.** `vertexCentredLinGeomSolid` adopts first, not
`linGeomTotalDispSolid`.

The two answers compose unusually well here, because in this solver the
tangent *is* the whole of `impK`: both consumers are Jacobian-only, so a
tangent-only slice is the complete `impK` migration for this solver, at zero
bit-for-bit risk.

<!-- markdownlint-disable MD013 -->

| Site | Today | After |
|---|---|---|
| `vertexCentredLinGeomSolid.C:1324` (`jacobianTangent scalar`) | `dualImpKf()` = `dualMechanicalPtr_().impKf()` | manager scalar tangent on dual faces |
| `:1331-1332` (`jacobianTangent fourthOrder`, the default) | `dualMechanicalPtr_().materialTangentFaceField(List<mat66>&)` | manager `mat66` on dual faces |
| `:1310` | `dualMechanicalPtr_().correct(dualSigmaf_)` | **unchanged** — stress stays legacy in this slice |
| `:822`, `:1037` | `average(mechanical().impK())`, `sqrt(impK/rho)` on the **primary** mesh | **unchanged** |

<!-- markdownlint-enable MD013 -->

`dualMechanicalModel` is therefore *not* removed by this PR; it is removed when
stress moves, in the following one. That is the intended thin slice.

**Correction to §3.6.** All three `vfvm` assembly routines iterate internal dual
faces only — `forAll(dualOwn, dualFaceI)` with `dualOwn = dualMesh.owner()`, at
`vfvmCellPoint.C:87` (`divSigma`), `:324` (`divSigma` with geometric
stiffness) and `:690` (`laplacian`). No boundary dual face is ever read.
So `dualFaceIntegrationPointTopology::nIntegrationPoints()` should return
`dualMesh().nInternalFaces()`, not `nFaces()`. This also matches
`dualGradDf_.primitiveField().size()` exactly, so the flat-list call needs no
sub-ranging. The legacy `List<mat66>` sizing to `nFaces()`
(`dualMechanicalModel.C:342`) was simply wasteful; D3 becomes moot.

**Consequence for the framework PRs.** The thin slice needs the tangent-only
query on the **flat-list** primitive, not only on the `volScalarField`
convenience overload of §3.1, because the values live on the dual mesh:

```c++
        //- Evaluate only the requested tangent at the integration points of
        //  the supplied topology, leaving the caller's stress storage and all
        //  constitutive state unmodified.
        void updateTangent
        (
            const integrationPointTopology& topo,
            const UList<tensor>& gradD,
            const UList<tensor>& gradD0,
            const scalar dt,
            UList<scalar>* scalarTangentPtr,
            UList<mat66>* fourthOrderTangentPtr,
            const tangentRequest tangentReq
        );
```

The `volScalarField` form of §3.1 becomes a thin wrapper over this.

### 8.5 Revised increment order

Supersedes the stage numbering in §5. Content is otherwise as described there.

<!-- markdownlint-disable MD013 -->

| PR | Content |
|---|---|
| 1 | **Done.** Defect fixes D1, D4, D5 and the `(4/3)*mu` change of §8.2 including `linGeomTotalDispSolid.C:574`; the `petscSnesPressure` u-p regression case; `planeStress` injection (§8.1) and support in the three new laws; `supportsFourthOrderTangent()` + manager guard. D3 is left to PR-2, which deletes the line. plateHole, wobblyNewton and perforatedPlate all pass. |
| 2 | **Done.** Flat-list `updateStressSmallStrain`/`updateStressFiniteStrain`/`updateTangentSmallStrain`/`updateTangentFiniteStrain` primitives, with the two `CompactListList` overloads and the internal-field half of the two `volTensorField` overloads re-expressed on them; `registerTopology()` and `topologyFor()` made public; `dualFaceIntegrationPointTopology`; defect D6. Plus `Test-mechanicalConstitutiveLaw`, run by the `layeredPipe` regression test. See §8.6. |
| 3 | **Done.** `solidModel::jacobianTangent(deflt)` and the `approximateJacobian` deprecation shim, used by both vertex-centred solvers; the shadow state of §3.1; a working `fourthOrderFiniteDifference` tangent. The optional manager on `solidModel` moves to PR-4. See §8.8. |
| 4 | **Done.** `vertexCentredLinGeomSolid` tangent-only adoption (§8.4), plus the optional manager deferred from PR-3. See §8.9. |
| 5 | Stress for `vertexCentredLinGeomSolid`, removing its dependence on `dualMechanicalModel`. Its nonlinear sibling is disabled, see §8.11. |
| 6 | **Done for `linGeomTotalDispSolid`.** `hofvm::divSigmaIntoPETScMatrix(mat66)`; the back-derivation is retained as the default rather than deleted, so existing cases are unchanged. The two `nonLinGeom*` solvers still carry it. See §8.13. |
| 7…n | Cell-centred solvers and their high-order variants, the ones in common use and now the priority, risk-ordered: `uns*` → `thermal`/`poro` → `nonLinGeom*` → `linGeomTotalDispSolid` → `coupledPressureDisplacementSolid`. |
| final | Retire the legacy tangent interface and migrate or bless the six `"impK"` registry consumers (OQ-8). |

<!-- markdownlint-enable MD013 -->

Every increment from PR-3 onward carries a foam-extend follow-on step, per
§8.7: verify on OpenFOAM.com first, then make that increment's code
foam-extend-clean and add it to `files.foamextend`. PR-1 and PR-2 owe that
step too; §8.7 lists exactly what it costs.

### 8.6 Notes from PR-2

**The tangent-only update is named per formulation.** §8.4 proposed a single
`updateTangent`. It is implemented as `updateTangentSmallStrain` and
`updateTangentFiniteStrain`, matching `updateStressSmallStrain` /
`updateStressFiniteStrain`. The finite-strain form is not needed until PR-5,
but it is the same shape and adding it now avoids revisiting the API. Both are
one-line wrappers over the corresponding stress primitive with manager-owned
scratch stress storage.

**The flat-list primitives are strict where the GeometricField overloads were
lenient.** A tangent request with no matching storage is a `FatalError` rather
than a silently skipped tangent. Because the GeometricField overloads now
route through the primitives, they inherit that strictness. No consumer exists
yet, so this costs nothing and removes a class of silent-wrong-Jacobian bug.

**A new defect, D6: the cell-centred topology type names were swapped.**
`cellCentredIntegrationPointTopology` declared
`TypeName("compactCellIntegrationPointTopology")` and vice versa. Self-
consistent, because the run-time selection entry and the manager's lookup both
went through `cellCentredIntegrationPointTopology::typeName`, so nothing
misbehaved — but `topo.type()` reported the wrong class in every diagnostic,
including the tangent guards added in PR-1, and the two classes registered the
same debug switch. Fixed.

**`registerTopology` is idempotent by design.** Constitutive state is keyed on
the topology *object*, so replacing the topology held under a key would
silently start a second, empty set of history variables. Re-registering a key
returns the topology already held; registering a key with a different topology
type is an error.

**The framework had no runtime coverage at all.** Nothing in `src/` or
`applications/` constructs a `mechanicalConstitutiveLawManager`, so everything
added in PR-1 and PR-2 was compile-tested only. `Test-mechanicalConstitutiveLaw`
now covers the closed-form stress and scalar tangent, agreement between the
flat-list, `CompactListList` and `GeometricField` paths, the tangent-only
update, the dual-face topology addressing and its fourth-order tangent, and the
four misuse guards. It runs on `layeredPipe` because that is the only tutorial
with two materials, and it is checked by that tutorial's `regressionTest.sh`.
A deliberate mutation of the flat-list primitive was confirmed to fail it.

Two limits worth stating. First, the path-agreement checks compare the
overloads with each other, and the overloads now share the primitive, so a
fault inside the primitive moves them together; only the closed-form check sees
that. Second, the test requires every material to be `linearElastic`, so it
does not exercise the plastic law or the finite-strain path.

### 8.7 foam-extend support is part of the plan, not an afterthought

`src/solids4FoamModels/Make/files.foamextend` contains **no**
`mechanicalConstitutiveLaw` entries, while `files.openfoam` contains fourteen.
The framework has therefore never been compiled on foam-extend, and "it builds
on foam-extend" has been true only by omission.

**Decision.** foam-extend is a supported fork and must remain one. The endgame
is that the legacy `mechanicalModel` and `mechanicalLaw` classes are deprecated
and removed; if the new framework never builds on foam-extend, that endgame
silently removes mechanical modelling from foam-extend altogether. That is not
acceptable, so the framework must end up building and running there.

**Sequencing.** Portability is a *follow-on step within each increment*, taken
once that increment's design is agreed and its OpenFOAM.com behaviour is
verified — not before, and not deferred to one large sweep at the end. Doing it
before the design settles means redoing it; doing it only at the end means an
undirected sweep across code nobody remembers. So from PR-3 onward, each
increment has two parts:

1. design, implement and verify on OpenFOAM.com;
2. make that increment's code foam-extend-clean, and extend `files.foamextend`
   to cover it.

**Status: done.** The framework now builds *and runs* on foam-extend-4.1. Its
fourteen sources are in `files.foamextend`, the test application's `FOAMEXTEND`
guard is gone, and all 23 checks pass there, with the finite-difference tangent
matching the analytical one to 3.8e-10 - the same figure as on OpenFOAM.com.

It was forced sooner than the per-increment plan intended: once a compiled
solid model referenced the manager, the foam-extend link broke, so the choice
was to finish the portability work or guard the reference away. The survey
below had already costed it, and the estimate held.

The `CompactListList` blocker turned out to be avoidable rather than fatal.
Two internal uses of its const `operator[]` are gone -
`compactCellIntegrationPointTopology` now stores its addressing flat, and
`compactCellTopologyFor` uses `sizes()`. **The container remains the interface
for higher-order discretisations with quadrature points per cell or per face**;
the manager's `CompactListList` overloads are untouched.

**What broke, as surveyed before the work.** A survey was run for PR-2: the fourteen
framework sources were added to `files.foamextend` in a throwaway worktree and
built against foam-extend-4.1. Every source was compiled. The result is much
better than feared — **all five constitutive laws, all five integration-point
topologies, the law base class and its run-time selector compile clean.** Only
two files fail, in three classes of defect, all in our own code and all
mechanical:

<!-- markdownlint-disable MD013 -->

| Class | Sites | Fix |
|---|---|---|
| `autoPtr<T>::operator*` does not exist in foam-extend | 10, all in `mechanicalConstitutiveLawManager.C` | use `ptr()`, or add an `operator*` shim to `compatibilityFunctions.H` |
| `GeometricField::boundaryFieldRef()` does not exist | 3, in `mechanicalConstitutiveLawManager.C` | use the existing `Foam::boundaryFieldRef()` from `compatibilityFunctions.H` |
| `Field<T>(label, const zero&)` and `setSize(label, const zero&)` | 9, in `mechanicalConstitutiveLawState.C` | use `pTraits<T>::zero` |

<!-- markdownlint-enable MD013 -->

Plus one blocker that is **not ours**: foam-extend's own
`CompactListList<T>::operator[](const label) const` returns
`UList<T>(m_.begin(), ...)`, passing a `const T*` where a `T*` is wanted. It
fires whenever that operator is instantiated. Options, in preference order:
avoid const-indexing a `CompactListList` in anything foam-extend compiles and
go through `m()` and the offsets directly; or carry a patched copy; or, as a
last resort, relax the flag for that translation unit. Note the `CompactListList`
overloads are a convenience layer over the flat-list primitive of §3.7, so the
first option is cheap.

Two caveats on that survey. It fixed nothing, so a second class of error may sit
behind the first in the two failing files. And it says nothing about whether the
framework *runs* correctly on foam-extend, only that it compiles — the
`Test-mechanicalConstitutiveLaw` application of §8.6 is the instrument for that,
and its `FOAMEXTEND` guard should be removed as soon as the framework builds
there.

### 8.8 Notes from PR-3

**The shadow state, and why not a copy.** §3.1 is corrected above. The one
thing worth restating: the threat a tangent query poses is not to history, which
`evaluate()` never writes, but to the **current-time** fields, which
`endTimeStep()` reads and which `storeOldTime()` promotes to history at the next
time step. A tangent query therefore needs somewhere to put a law's *outputs*,
not a copy of its *inputs*.

**A latent contract violation this exposed.** The plastic law reached its
old-time fields through the non-const accessor, e.g.
`state.symmTensorField0("epsilonP")`, because `state` is a non-const reference
and overload resolution then picks the create-and-modify overload. The shadow's
guard rejected it immediately. Fixed by reading history through
`getSymmTensorField0()` / `getScalarField0()`, which say what is meant. The rule
is now enforced rather than assumed: **inside `evaluate()`, read `*0` fields
only through the `get*` accessors, and write only current-time fields.**

**`fourthOrderFiniteDifference` did not work at all.** It contained a
`FatalErrorInFunction` placeholder reading "In FD material tangent, make copy of
state", and separately computed each finite-difference column but never wrote it
into the `mat66`. This settles OQ-5: the mode was advertised in `README.md` but
could not run. It is now implemented.

**It is evaluated field-wise, not point-wise.** The original sketch looped six
Voigt components inside a loop over integration points, i.e. `6*N` calls into
the law with one-element views. It now perturbs every integration point at once
and calls the law six times in total. Same arithmetic, one sixth of the calls,
better vectorisation, and one shadow for the whole assembly instead of one per
point.

**The perturbation convention was wrong for shear.** `gradDPerturbed` added `h`
to both off-diagonals of `gradD`, which moves `eps_xy` by `h` and therefore the
engineering shear `gamma_xy = 2*eps_xy` by `2h`. Against an analytical `mat66`
where `C(XY,XY)` is `mu`, the finite difference came out a factor of two high on
every shear column. The perturbation is now defined on the Voigt strain vector:
a shear component moves each off-diagonal by `h/2`. Verified by comparing the
finite-difference and analytical tangents for `linearElastic`, which now agree
to 4e-10.

**The plastic law now supports a finite-difference tangent** and refuses an
analytical `fourthOrder` one with a message pointing at the alternative.
Previously it silently ignored a fourth-order request and left the caller's
`List<mat66>` holding uninitialised memory.

**A test that looked right and proved nothing.** The first version of the
state-preservation check evaluated a stress, took tangent queries at wildly
perturbed kinematics, evaluated the same stress again, and required the two to
agree. They always do - with or without the shadow - precisely because a law
recomputes current-time state from old-time state. Disabling the shadow did not
fail it. The check now straddles a time step, so that `storeOldTime()` commits
whatever the queries left behind; without the shadow it fails by 114%.

**Deferred from PR-3.** §8.5 also listed an optional manager owned by
`solidModel` behind `useMechanicalConstitutiveLawManager`. It is left to PR-4,
which is where the first consumer appears: it is inert on its own, and the
question of which dictionary object the manager is constructed from is better
answered next to a caller than in the abstract.

### 8.9 Notes from PR-4

`vertexCentredLinGeomSolid` now takes its Jacobian material tangent from the
framework, behind `useMechanicalConstitutiveLawManager` (default `no`). Both
arms move: the scalar tangent that fed `vfvm::laplacian` and the `mat66` field
that fed `vfvm::divSigma`. The residual stress at dual faces still comes from
`dualMechanicalModel`, so this is the thin slice §8.3 asked for.

**Verified byte-identical.** On `cantilever2d` with `solidProperties.vertexCentred`,
the manager and legacy paths produce identical output fields and the same SNES
iteration count. Confirmed the switch actually engaged, rather than the result
being identical because nothing happened: the new framework constructs a law
only in the manager run.

**The manager path no longer depends on `dualImpKf()`.** The first version
built its scalar field by copying it, which pulled in
`dualMechanicalModel::impKf()` and with it the `interpolate(impK)`
interpolation scheme - on a path whose whole purpose is to need neither. The
field is now constructed directly.

**Caveat worth carrying into PR-5: the manager's constitutive state is never
advanced in a tangent-only slice.** Nothing calls a stress update on the
manager, so its current-time fields stay at their initialised values and each
time step commits those to history. For `linearElastic` this is vacuous, since
the law has no state. For a history-dependent material the manager would
return the tangent at the *initial* state rather than the consistent tangent,
which is a Jacobian-quality issue rather than a correctness one, but it means
the tangent-only slice is only fully meaningful for elastic materials. It
resolves itself in PR-5, when stress moves and the state starts advancing.

**Not exercised end to end: the scalar Jacobian arm.** No tutorial runs
`vertexCentredLinGeomSolid` with a scalar tangent, and doing so needs two
dictionary entries no case provides - `compactImplicitStencil`, which has no
default, and `interpolate(impK)` in `system/dualMesh/fvSchemes`. Both
requirements pre-date this work. The scalar tangent *value* is covered by
`Test-mechanicalConstitutiveLaw` instead.

**From an independent review of PR-3.** A shadow constructed from a temporary
would dangle; the rvalue overload is now deleted. And the tangent-only updates
do participate in the once-per-time-step old-time rollover, so if one is the
first call of a new time step it is what commits the previous step's converged
state. That is bookkeeping owed to whoever evaluates first rather than
something the query computed - shadows never write current-time fields, so the
committed values are the same whichever call triggers it - but the contract now
says so instead of claiming the query modifies nothing at all.

### 8.10 Correction: boundary dual faces *are* read, by the residual

§8.4 states that "no boundary dual face is ever read" and concludes that
`dualFaceIntegrationPointTopology::nIntegrationPoints()` should be
`dualMesh().nInternalFaces()`. **That is true of the Jacobian and false of the
residual**, and the distinction was missed because only the Jacobian had moved.

The three `vfvm` assembly routines do iterate `dualMesh.owner()`, i.e. internal
dual faces only. But the residual path in
`vertexCentredLinGeomSolid::updatePointDivSigma` does this:

```c++
surfaceVectorField dualTraction(dualN & dualSigmaf_);
enforceTractionBoundaries(pointD, dualTraction, mesh(), ...);
const vectorField dualDivSigma = fvc::div(dualTraction*dualMesh().magSf());
```

`fvc::div` of a surface field sums fluxes over every face of each cell,
boundary faces included, so the **boundary field** of `dualSigmaf_` is read.
`enforceTractionBoundaries` overwrites it only on patches whose `pointD`
boundary condition is a `solidTractionPointPatchVectorField`; a
fixed-displacement patch keeps `dualN & dualSigmaf_`. `cantilever2d` has a
clamped end, so this is live in the one case that covers this solver.

**Consequence.** The tangent-only slice of PR-4 is unaffected: it only ever
needed internal faces. Moving *stress* onto the manager needs boundary dual
faces as well, which the present topology deliberately excludes. The stress
move is therefore blocked on a decision about how to represent them.

**Recommended.** One topology covering all dual faces, indexed
`[0, nInternalFaces)` for internal and `[nInternalFaces, nFaces)` for boundary,
with the solid model gathering `dualGradDf_` and scattering `dualSigmaf_`
across the internal field and the boundary patches. One topology means one set
of constitutive state, which is the property that matters: a history-dependent
law must not have its plastic strain split across two state objects. The
Jacobian then reads the first `nInternalFaces` entries of the same arrays, and
`vfvm::divSigma` is indifferent to a longer list.

Rejected: a second, boundary-only topology, because it gives each dual face
family its own `mechanicalConstitutiveLawState` and so splits history. Also
rejected as an endpoint: leaving boundary stress with `dualMechanicalModel`,
because then it can never be removed - though it would work as a stepping
stone.

### 8.11 `vertexCentredNonLinTotalLagGeometry` is disabled

Its compilation is commented out in `Make/files.openfoam` and
`Make/files.foamextend`, and the reason is recorded in the solver's own
`README.md`. No tutorial ever selected it, and an attempt to give it one failed:
three of its dictionary entries have no defaults and are set by no case, and
with those supplied the run dies in a floating point exception inside the PETSc
solve on the first Newton step.

That state pre-dates this work and is independent of it. Fixing it is separate
from the constitutive law migration, so the migration skips it rather than
carrying an untestable solver.

**Consequence for the plan.** PR-5's original content - the nonlinear
vertex-centred tangent, then stress for both vertex-centred models - reduces to
stress for `vertexCentredLinGeomSolid` alone. The priority after that moves to
the cell-centred solid models and their high-order variants, which are the ones
in common use.

### 8.12 Notes from the first cell-centred adoption

`linGeomTotalDispSolid` now sources `impK` - and `rKappa` - from the framework
behind `useMechanicalConstitutiveLawManager`, default `no`. `impKf_` and
`rImpK_` stay derived from `impK_` exactly as before, `impK_` is still
registered under `"impK"` for the contact and cohesive-zone models that look it
up by name, and the stress still comes from `mechanicalModel`.

**Byte-identical with the switch on**, on `plateHole` (linear elastic) and on
`perforatedPlate` (elastoplastic, twenty time steps), including the yielding
diagnostics. The switch was confirmed to engage rather than the agreement being
vacuous: the framework constructs a law only in the manager runs.

The elastoplastic case agreeing exactly is worth understanding rather than just
noting. `impK_` is frozen at construction, and on a cold start the constitutive
state is zero, so the legacy state-dependent `impK()` and the framework's
tangent query are evaluated at the same state. That is OQ-1 restated: the
agreement is exact from a cold start and would not be from a restart, for
either implementation, because both freeze a state-dependent quantity at
construction.

`plateHole` gains a `segregatedManager` approach so the path has continuous
coverage rather than a one-off manual check.

**A gap this found in `updateScalarTangent`.** The flat-list primitive fills
internal integration points only, so the returned `volScalarField` had a
zero boundary field, and the first thing this caller does is form `1/impK_` -
a floating point exception in the constructor. The manager now fills each
non-coupled patch from its patch-internal value, which is exact for a scalar
tangent because a boundary face belongs to its owner cell's material, then
syncs the coupled patches. Any cell-centred caller would have hit this.

### 8.13 Notes from the high-order Jacobian

`hofvm::divSigmaIntoPETScMatrix` assembles the high-order Jacobian from a full
`mat66` per mesh face, and `linGeomTotalDispSolid` uses it when
`jacobianTangent fourthOrder` is set. The `(impK_ - K)*3/4` back-derivation
stays as the default, so existing cases are untouched.

**The case against the back-derivation is structural, not numerical.** See the
correction in §3.5: it inverts exactly for every law currently defining `impK`
in the standard form, plastic points included. What is wrong with it is that it
depends on an undocumented contract about the form of `impK()`, and that no
isotropic pair can represent an anisotropic tangent however it is obtained. A
solid model does not need material constants; it needs a Jacobian operator, and
the material tangent is the general way to supply one.

**No new mathematics was needed, as §3.5 predicted.** The three isotropic
kernels sum to `w*(mu*(n.g)*I + mu*g_i*n_j + lambda*n_i*g_j)`, which is
`Sf_m C_mikl g_k delta_lj` with `Sf = w*n` - already implemented as
`multiplyCoeff` and already used by the vertex-centred solver. The assembler
was generalised to take either a scalar diffusivity with one of the kernels or
a material tangent, rather than being duplicated, so the internal, processor,
symmetry and generic patch handling stays in one place.

**Verified byte-identical**, with the same Newton iteration count and the same
convergence reason, on `plateHole` with a linear elastic material - which is
precisely the case where the back-derivation is valid, and therefore the only
case where the two *should* agree. `plateHole` gains a `highOrderFourthOrder`
approach so this stays covered.

**A guard was relaxed to make this possible.** The flat-list update refused any
topology whose integration points are shared between cells. The reason for that
rule is stress collapse at a material interface, which cannot arise with a
single law, so it now refuses only when there is more than one. The face-centred
topology is shared by construction, so without this the material tangent could
not be obtained per face at all.

The test application asserted the old blanket rule and duly failed on
`perforatedPlate`, which has one material. It now asserts the real contract:
rejected with several materials, allowed with one.

### 8.14 Boundary integration points in the flat-list update

`faceCentredIntegrationPointTopology::nIntegrationPoints()` returns
`mesh.nFaces()`, but its cell-to-integration-point map covers **internal faces
only**, and says so. Laws are matched to integration points through that map,
so a flat-list update left every boundary entry of the caller's storage exactly
as it was found. `linGeomTotalDispSolid::faceMaterialTangent()` sizes its
`List<mat66>` to `nFaces` and hands it to an assembler that does iterate the
boundary patches, so it was reading unwritten memory.

**Fixed** by giving the topology a `boundaryIntegrationPointIDs(patchI)` query
and having both flat-list primitives evaluate the boundary points of a
`boundaryAware` topology, using the per-law, per-patch boundary state the
manager already keeps. This is the general fix: any topology that declares
itself boundary-aware now has its boundary points covered, rather than the
caller working around the gap.

**Why this went unnoticed, which is the part worth remembering.** The
`highOrderFourthOrder` variant was byte-identical to `highOrder`, with the same
Newton iteration count, and that was reported as strong verification. It is
not. This is a *Jacobian*: a wrong one changes the path to the solution, not the
solution, and PETSc SNES converges to its tolerance either way. **A
byte-identical converged result is exactly what a wrong Jacobian looks like
when the residual is right.** No end-to-end solver comparison can validate a
Jacobian; only a direct comparison of the tangent itself can.

It was found by the finite-strain tangent test, which sized its comparison to
`nIntegrationPoints()` and got a relative error of exactly 1 on the unwritten
entries - a number too clean to be anything but "never assigned".

### 8.15 The boundary check earns its place, and how that was established

The finite-strain tangent test now compares **every** integration point of the
face-centred topology, splitting the report into an internal and a boundary
check rather than stopping at `nInternalFaces`. That restriction existed only
as a workaround for the defect of section 8.14, so removing it is what turns
the test into the regression guard for the fix.

Two things make the boundary check trustworthy rather than decorative:

1. The comparison buffer is **poisoned** with `GREAT` before the call. `mat66`
   is a POD, so an unpoisoned `List<mat66>` holds whatever memory held, and an
   integration point the manager never reaches would be compared against
   arbitrary values - a test that passes or fails by luck.

2. It was **mutation tested**. With `boundaryIntegrationPointIDs` reverted to
   returning an empty list, i.e. exactly the pre-fix behaviour, the internal
   check still passes and the boundary check fails with a relative error of
   1.7e12: the poison, untouched. With the fix in place both pass at 2.2e-4.
   The check fails when, and only when, the defect is present.

Coverage is on `blockPunch`, whose law is `neoHookeanElastic` and which
therefore has no small-strain evaluation at all. It is the only runtime
coverage of the finite-strain kinematics, the finite-difference fourth-order
tangent, and the boundary integration points. Its material is given as `mu` and
`K`, which is why the law was extended to accept that form alongside `E` and
`nu`.

**A build lesson that cost time twice.** The first mutation run reported a
pass. The mutation was genuinely in the library, but the test application had
not been relinked, because the library and the applications were built in
separate steps and only the library step was repeated. The library and the
applications must be rebuilt **together** before any result from the test
application means anything - the same trap as the vtable mismatch in section
8.14. A test result from a half-rebuilt tree is not evidence.

### 8.16 Empty patches: the boundary fix was still incomplete

Section 8.14 fixed the flat-list update to evaluate boundary integration
points. It was still wrong, and the first two-dimensional case to exercise it
said so immediately.

`lawBoundaryFaces_` is built from `mesh.boundary()`, the **fv** boundary. An
`empty` patch has no fvPatch faces at all, so that list is empty for it, and
the code skipped the patch. But an empty patch's **polyPatch** faces still
occupy slots in the face-centred topology's index space, because
`nIntegrationPoints()` is `mesh.nFaces()`. Those slots were left unwritten -
precisely the defect of section 8.14, surviving in the two-dimensional case.

`blockPunch` is three-dimensional and has no empty patch, so the regression
guard added in section 8.15 could not see it. `rotatingCylinder` is
two-dimensional, and its boundary check failed on the first run with a
relative error of 3714 - which is `GREAT/(lambda + 2 mu)`, i.e. the poison,
untouched. The poison is what made the diagnosis immediate rather than a
puzzle: an unpoisoned buffer would have shown an arbitrary number.

**Fixed** by addressing such patches through the polyPatch and taking the law
from the owner cell. Two consequences worth stating:

- These faces take no part in the finite volume discretisation - that is what
  `empty` means - so they are evaluated against a **scratch state**, not the
  per-patch state, which does not exist for them. They get a well-defined
  value without committing history for a face the discretisation does not
  have.
- That scratch state must be handed to `initialiseState`, exactly as a real
  state is. Skipping that made `perforatedPlate` fail with "Requested state
  field 'epsilonP' does not exist" - a history-dependent law failing to find
  its own history. A bare state of the right size is not a usable state.

**The lesson repeats.** Each round of this defect was found by a case with a
mesh feature the previous cases lacked, never by reasoning about the code.
Coverage of a topology contract means coverage of the mesh variety the
contract has to survive: 3D, 2D with empty patches, and - still untested here
- wedge, cyclic and processor patches.

### 8.17 The blocker is law coverage, not solver plumbing

Migrating `nonLinGeomTotalLagTotalDispSolid` and `nonLinGeomUpdatedLagSolid`
turned out to be straightforward: both hold their kinematics as volFields and
evaluate stress at cell centres only, so no face-based finite-strain
kinematics were needed after all. What stopped both from being *verifiable*
was the set of laws the framework implements.

Counting the tutorials by law:

| Solver | Tutorials | Laws they use |
|---|---|---|
| `linGeomTotalDispSolid` | 16 | `linearElastic`, `linearElasticMisesPlastic` |
| `nonLinGeomUpdatedLagSolid` | 9 | `neoHookeanElasticMisesPlastic`, `MooneyRivlinElastic` |
| `nonLinGeomTotalLagTotalDispSolid` | 8 | `MooneyRivlinElastic`, `StVenantKirchhoffElastic` |

The framework had `linearElastic`, `linearElasticMisesPlastic` and
`neoHookeanElastic`. That covers the first row completely and the other two
not at all: the only neo-Hookean cases using the finite-strain solvers are
FSI. Adding `StVenantKirchhoffElastic` bought coverage of one case,
`rotatingCylinder`, and that single case immediately found the empty-patch
defect of section 8.16.

**`MooneyRivlinElastic` does not fit the current law contract.** It calls
`updateSigmaHyd`, which, when `pressureSmoothingScaleFactor` is set, *solves a
Laplacian equation* for the hydrostatic stress over the whole field. A
`mechanicalConstitutiveLaw` is a pure function of integration-point kinematics
and old-time state; a field-level implicit solve is not expressible in that
interface.

**Resolved: the smoothing does not come across.** Pressure smoothing is the
solid model's business, not the law's - it stabilises the discretisation, it
does not describe the material - and the field-level solve is to be added at
the solid model level separately. The framework law therefore drops
`updateSigmaHyd` entirely and returns the unsmoothed hydrostatic stress
`p = K (J^2 - 1)/2`, which keeps `mechanicalConstitutiveLaw` a pure function
of integration-point kinematics and old-time state.

That decision is what made `MooneyRivlinElastic` implementable here. Note the
consequence: until the solid model provides the smoothing, a case that relied
on `pressureSmoothingScaleFactor` will converge differently under the
framework. The two migrations in this PR are unaffected, because they take
only `impK` from the framework and never call a law for a stress.

With `MooneyRivlinElastic` in place the coverage gap partly closes.
`longWall` exercises the total Lagrangian solver's framework path in CI.

The updated Lagrangian solver was verified **by hand, not in CI**, on
`cylinderCrush` under foam-extend, which is the only tutorial pairing that
solver with a law the framework implements. Both arms agree exactly:
force_y = -70145.7 N and disp_y = -0.0998977 m with the legacy impK and with
the framework impK, and the framework's own constitutive checks pass on that
mesh.

It is deliberately **not** wired into CI, for a reason worth recording. Run
locally the case completes all 30 time steps of its 30 s ramp. The regression
bands it is checked against, force_y in [-550, -545] N and disp_y in
[-0.0034, -0.0032] m, correspond instead to roughly t = 1 s: the displacement
ramps to about -0.0999 m at t = 30, and one thirtieth of that is -0.0033,
which is what CI reports (-0.00329221). So in CI the case appears to stop
after about one time step, and the expected values encode that truncated run
rather than the completed one.

Until that discrepancy is understood, this case is not a trustworthy vehicle
for a two-arm comparison, and doubling a long contact run in CI on top of an
unexplained early stop is not a good trade. The single-arm test is left
exactly as it was.

`neoHookeanElasticMisesPlastic`, which nine tutorials need, remains the other
gap. It is history-dependent finite-strain plasticity, a substantial law in
its own right, and is left for a later increment.

One further limitation worth recording: the framework's material constants are
uniform per law, whereas the legacy `MooneyRivlinElastic` holds `c10`, `c01`
and `c11` as volScalarFields that may be read per cell. Spatially varying
constants belong with the wider question of how the framework carries
per-integration-point material data, alongside initial stress and fibre
directions - see `DESIGN-state-io.md`.

### 8.18 impK under the mixed formulation: the two solvers disagreed

Review caught that the framework arm of `nonLinGeomTotalLagTotalDispSolid`
does not reproduce its legacy `impK` when `solvePressure()` is set. The legacy
code uses `2*mu`; the framework request `scalarDeviatoric` yields `(4/3)*mu`.

That is a real difference, and worth stating rather than quietly matching,
because the two migrated solvers did not agree with each other to begin with:
`linGeomTotalDispSolid` already used `(4/3)*mu` for exactly the same
formulation. So there was no single legacy convention to preserve.

`(4/3)*mu` is kept. It is the value the surrogate implies - the scalar
Laplacian standing in for `div(dev(sigma))` is
`mu*lap(D) + (1/3)*mu*grad(div(D))` - and it makes the two solvers agree.
`impK` is the coefficient of a Laplacian that is added implicitly and
subtracted explicitly, so it sets how the solution is reached and not what it
is: only the convergence path changes, and only when the framework is switched
on. Nothing in the default configuration moves.

The general point: "matches the legacy value" is not always a well-defined
target, and where the legacy values disagree among themselves it is better to
pick the defensible one and say so than to reproduce an arbitrary choice per
solver.

### 8.19 The stress comes from the framework, and what that cost

Every migration before this took only `impK`, which cannot change results, so
the framework's central path had never been driven by a solver. All three
cell-centred solid models now take their stress from it behind the same
switch: `linGeomTotalDispSolid`, `nonLinGeomTotalLagTotalDispSolid` and
`nonLinGeomUpdatedLagSolid`.

Three results, in increasing order of interest:

* `plateHole`, small strain, `linearElastic`: the legacy and framework arms
  agree to **every digit**.
* `rotatingCylinder`, finite strain, `StVenantKirchhoffElastic`: likewise
  exact, sigmaEq 3700.73 both ways.
* `longWall`, finite strain, `MooneyRivlinElastic`: agree to about **2e-6
  relative**, not exactly. This case sets `solvePressureEqn`, so the legacy law
  solves a Laplacian for its hydrostatic stress and the framework law does not.

That last one is the useful one. The difference is the pressure smoothing whose
removal §8.17 justified on the grounds that it stabilises the discretisation
rather than describing the material. A 2e-6 change in the converged answer is
what that claim predicts, and is the first evidence for it. One case is not a
proof - a nearly incompressible material could be far more sensitive - so the
claim stays provisional until the solid model owns the smoothing and the two
can be compared with it present.

**Identical numbers are not evidence on their own.** They are equally what a
switch that does nothing produces, which has already happened once in this
work. So the framework path was confirmed independently: the solver log shows
the manager constructed, and perturbing the framework stress by 0.1% moves
`plateHole`'s converged result from 2.22246e+06 to 2.22248e+06, which a dead
branch cannot do.

**What is not migrated.** The high-order quadrature stress stays on the legacy
path in both cell-centred solvers that have it. The framework's
`CompactListList` overload exists, but evaluating a law through it needs the
old-time gradient at the quadrature points, and neither solver keeps one:
`gradDQuad()` is rebuilt each call. A history-independent law would not notice;
a history-dependent one would silently read the wrong state. Closing that needs
an old-time quadrature gradient, which is a solver change, not a framework one.

### 8.20 Sharing a TypeName with a legacy law has a portability trap

Every framework law deliberately takes the same `TypeName` as the legacy law
it replaces, so that a case dictionary needs no change. That is the right
choice, and it has one consequence that is invisible on OpenFOAM.com.

Debug switches are registered **globally by name**, not per runtime-selection
table. Two classes sharing a name therefore register the same switch, and if
they declare *different* defaults, OpenFOAM.org fails at start-up with

```text
    Multiple defaults set for debug switch neoHookeanElasticMisesPlastic
```

before any case runs. Every tutorial in the CI job failed, including ones that
never touch the framework, because the error is thrown during static
initialisation. OpenFOAM.com tolerates the same mismatch silently, so the whole
thing is invisible locally.

It happened because `neoHookeanElasticMisesPlastic` was scaffolded from
`linearElasticMisesPlastic`, which legitimately uses a debug default of 1,
while its own legacy counterpart uses 0.

**Rule: a framework law's debug default must equal that of the legacy law
whose TypeName it shares.** All seven pairs were checked; only that one
differed. Worth re-checking whenever a law is added, because the symptom is a
total, unrelated-looking CI failure on one fork only.

### 8.21 The plastic law was wrong, and only a real solve found it

Running `perforatedPlate` through both paths - the first end-to-end comparison
of a *plastic* framework law against its legacy counterpart - showed the
framework giving about 5% more equivalent strain and 2.6% more stress.

Isolating the cause mattered, because two things had changed at once. Forcing
the legacy `impK` while keeping the framework stress reproduced the discrepancy
exactly, so it was not a Jacobian or convergence artefact: the stress itself
was wrong.

**The defect.** `linearElasticMisesPlastic` disagreed with itself about what
its plastic multiplier meant. It set

```text
    DEpsilonP = 1.5 dLambda n      so  |DEpsilonP| = 1.5 dLambda
```

whose equivalent plastic strain increment is
`sqrt(2/3) |DEpsilonP| = sqrt(3/2) dLambda`, but then accumulated
`epsilonPEq += sqrt(2/3) dLambda` - short by a factor of 3/2 - and the
closed-form hardening denominator carried the same error.

**Why every existing check missed it.** Under perfect plasticity the *stress*
is still right: the multiplier and the stress update are consistent with each
other, and only the accumulated equivalent strain is wrong. It goes wrong when
hardening is active, because the yield stress is then evaluated at the wrong
accumulated strain. The closed-form tangent check, the finite-difference
tangent check and the saturation check are all blind to that.

**The fix** is the standard convention in which `dLambda` *is* the equivalent
plastic strain increment:

```text
    q           = qTrial - 3 mu dLambda
    DEpsilonP   = sqrt(3/2) dLambda n
    DEpsilonPEq = dLambda
    dLambda     = (qTrial - sigmaY0) / (3 mu + Hp)
```

which is self-consistent. The framework then reproduces the legacy result
exactly: epsilonEq 0.00507147, sigmaEq 1.79161e+08, and 30 yielding
integration points against 30 yielding cells.

**The lesson.** Unit checks on a constitutive law verify the parts of it you
thought to check. A history-dependent law compared against a trusted
implementation over a real load path checks the parts you did not. This is the
argument for migrating the stress rather than stopping at `impK`: an `impK`
migration cannot change results, and therefore cannot detect this class of
error at all.

Two supporting gaps were fixed alongside it.
`mechanicalConstitutiveLawManager::endTimeStep()` existed but nothing called
it, so framework laws had no end-of-step hook and no diagnostic; it is now
called from all three migrated solid models. And the regression test's
extraction preferred the legacy yielding message, which in the framework arm is
still emitted by the idle legacy law and truthfully reports zero, masking the
framework's own count.

### 8.15 Open questions still outstanding

OQ-1 (restart `impK` state-dependence), OQ-3 (configurable `impKf_` cadence),
OQ-5 (`fourthOrderFiniteDifference` production status), OQ-6 (stabilisation term
in the `fourthOrder` Jacobian), OQ-7 (`mat39` ownership), OQ-8 (`"impK"` as
public API), and the incompressible-limit half of OQ-9.

OQ-2 is deferred to PR-5 by §8.3. OQ-4 is closed by §8.2.
