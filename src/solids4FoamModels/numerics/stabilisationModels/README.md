---
sort: 1
---

# Stabilisation Models

This directory contains runtime-selectable stabilisation models used by
`solids4foam` to add artificial dissipation for suppressing grid-scale
oscillations in scalar and vector fields.

## Design Overview

The common base class is
[`stabilisationModel`](https://github.com/solids4foam/solids4foam/blob/development/src/solids4FoamModels/numerics/stabilisationModels/stabilisationModel/stabilisationModel.H).

Each concrete stabilisation model is responsible for:

- reading its configuration from a dictionary,
- computing a face-based stabilisation contribution,
- optionally providing an approximate Jacobian for implicit use.

The classes are intentionally lightweight. The current interface is based on
an explicit update step followed by read access to cached results. This lets a
solid model control when the stabilisation is recalculated and re-use the
result more than once without repeating the work.

## Runtime Selection

Models are selected from dictionary input using the standard OpenFOAM runtime
selection mechanism:

```text
stabilisation
{
    momentum
    {
        type        diffStencilLaplacian;
        scaleFactor 0.1;
    }

    pressure
    {
        type        diffStencilLaplacian;
        scaleFactor 0.1;
    }
}
```

The currently available model types in this directory are:

- `alpha`
- `diffStencilLaplacian`
- `generalisedEvenOrderLaplacian`
- `gradDiv`
- `JamesonSchmidtTurkel`
- `laplacian`
- `multiple` — composite model combining two or more models (see below)
- `volStrainRate` — temporal stabilisation based on the rate of change of
  volumetric strain (see below)

`RhieChow` is implemented as a helper under
[`diffStencilLaplacianStab/RhieChowStab`](https://github.com/solids4foam/solids4foam/blob/development/src/solids4FoamModels/numerics/stabilisationModels/diffStencilLaplacianStab/RhieChowStab).

## Spectral Normalisation

The Laplacian-family models support optional reference-mode normalisation:

```text
spectralNormalisation        true;
referenceNyquistDirections   2;
```

`spectralNormalisation` defaults to `false`, preserving the legacy formulas
and dictionary behaviour. Its earlier spelling, `normalise`, is still read,
with a deprecation warning; a model dictionary that gives both is refused,
whatever their values. `referenceNyquistDirections`, written `rStar` below,
is required when normalisation is enabled and must be in the range
`1 .. min(3, mesh.nGeometricD())`: the number of mesh directions in which the
isotropic orthogonal reference mode is at the Nyquist frequency. Empty and
wedge directions, which have one cell, do not count.

The order is m = `laplacianPower`, counted from zero: m = 0 is the Laplacian,
m = 1 the bi-Laplacian and m = 2 the sixth-order operator, with m + 1
repeated Laplacian factors. The normalisation divides the operator by its
raw response to the reference mode, so that the normalised operator's symbol
at that mode equals `scaleFactor`:

<!-- markdownlint-disable MD013 -->

| Model | m | Spectral relation | Raw reference response | Normalisation |
| --- | ---: | --- | ---: | ---: |
| `laplacian` | 0 | Coupled repeated-Laplacian | `4*rStar` | `1/(4*rStar)` |
| `generalisedEvenOrderLaplacian` | `laplacianPower` | Coupled repeated-Laplacian | `(4*rStar)^(m+1)` | `(4*rStar)^(-(m+1))` |
| `JamesonSchmidtTurkel` | 1 | Coupled repeated-Laplacian | `(4*rStar)^2` | `(4*rStar)^(-2)` |
| `diffStencilLaplacian`, `RhieChow` | — | Directionally split m = 1 shape, `sum_q a_q^2` | `4*rStar` | `1/(4*rStar)` |

<!-- markdownlint-enable MD013 -->

The symbol is that of the operator applied to a pressure field (the face
flux, or the cell field built from it), on a uniform, orthogonal, periodic
Cartesian mesh, where the coupled formulas are exact. It is not the magnitude
of the complete stabilisation term in a solver residual, which is further
weighted by `rAUf`, the face areas, the residual scale and the cell volume
(see m = 0 below). `scaleFactor` is therefore the strength of the normalised
operator, not the absolute damping of the residual. On such a mesh, with a
centred cell gradient and linear interpolation, the
`diffStencilLaplacian`/`RhieChow` operator has a normalised spectral shape
equal to that of the directionally split m = 1 (bi-Laplacian) family; its raw
prefactor differs, and no higher-order directionally split family is
implemented. `stabilisationFourierCheck` and the
`tutorials/solids/linearElasticity/pressureStabilisationFourier` case verify
these symbols against the implemented operators.

### m = 0

For `laplacian`, and `generalisedEvenOrderLaplacian` with
`laplacianPower 0`, the face flux is computed with the orthogonal normal
gradient, `orthogonalSnGrad(p) = deltaCoeffs*(p_N - p_P)`. On an internal
face `1/deltaCoeffs = |d|`, the distance between the centroids of the two
cells sharing the face; on a boundary face of a non-coupled patch it is the
patch-normal distance from the owner-cell centroid to the face,
`n & (Cf - C_P)`, and the boundary value of `p` is that of its boundary
condition. The legacy flux is `scaleFactor*|d|^2*orthogonalSnGrad(p)`, with
the same `|d|`. A scalar coefficient cannot remove this mesh-dependent
factor, so with normalisation the m = 0 branch uses
`scaleFactor/(4*rStar)*orthogonalSnGrad(p)`. With normalisation disabled the
legacy expression is unchanged. A `scaleFactor` tuned with
`spectralNormalisation false` therefore does not carry over to
`spectralNormalisation true` (on a uniform mesh,
`s_normalised ≈ 4*rStar*h^2*s_legacy`).

Dropping `|d|^2` does not remove the mesh scaling of the stabilisation,
because the face flux is weighted by `rAUf`. The complete contribution to the
residual of the pressure equation of cell P, as assembled by the mixed solid
models (`linearGeometryTotalDisplacement` and
`nonLinearGeometryTotalLagrangianTotalDisplacement`), is

```text
R_P += s_r * sum_f rAUf_f * |Sf| * scaleFactor/(4*rStar) * deltaCoeffs_f*(p_N - p_P)
```

where `s_r` is the pressure-equation scale of `scaleMixedPETScFields` (1 when
scaling is off), and
`rAUf = -1/interpolate(A(laplacian(impKf, D) - rho*d2dt2(D)))` is the
reciprocal of the diagonal of the approximate momentum operator. On a uniform
orthogonal Cartesian mesh of spacing h in nD dimensions, with constant
implicit stiffness `impK` and no inertia or relaxation,
`|rAUf| = h^2/(2*nD*impK)`, and per unit volume the term is

```text
R_P/V = s_r * scaleFactor/(4*rStar) /(2*nD*impK) * sum_f (p_N - p_P)
```

So:

- for a smooth pressure, `sum_f (p_N - p_P) = h^2*laplacian(p) + O(h^4)`, and
  the term is `O(h^2*scaleFactor*p/impK)`, of the order of the discretisation
  error, so the formal order of accuracy is retained;
- for the reference grid-scale mode (Nyquist in `rStar` directions),
  `sum_f (p_N - p_P) = -4*rStar*p`, and the term is
  `-s_r*scaleFactor*p/(2*nD*impK)` per unit volume, independent of h: the
  damping of grid-scale modes does not vanish under refinement, which is the
  purpose of the normalisation. Its size relative to the other terms of the
  pressure equation is set by `scaleFactor`, `impK` (2μ with
  `mixedImplicitStiffness twoMu`) and `s_r`, not by `scaleFactor` alone.

The legacy form carries a further `|d|^2`: `O(h^4)` for smooth fields and
`O(h^2)` for grid-scale modes, so its damping of grid-scale modes vanishes
under refinement.

These estimates assume a uniform orthogonal mesh, a constant `impK`, and that
the momentum diagonal is dominated by the stiffness. Inertia adds
`rho/deltaT^2` to the diagonal and so reduces `|rAUf|` (and the
stabilisation) when `deltaT` is small compared with `h*sqrt(rho/impK)`;
under-relaxation of the approximate momentum equation and mesh
non-uniformity modify `rAUf` locally.

Supporting tests: the reference response `4*rStar` and its normalisation are
checked by `stabilisationFourierCheck` (the `pressureStabilisationFourier`
tutorial), and the manufactured-solution case with normalised m = 0
(`scaleFactor 1`, `referenceNyquistDirections 3`) converged on 5^3 to 40^3
meshes with net orders 1.50 (RMS) and 0.91 (maximum) for the pressure and
2.17 for the displacement, as with its unnormalised m = 1 configuration.

### Approximate Jacobian

Normalisation applies to the residual only. `scalarJacobian` is
`scaleFactorJacobian*laplacian(gamma, p)`, with the `laplacian` scheme
selected in `fvSchemes`, whether or not `spectralNormalisation` is on, so
`scaleFactorJacobian` remains an independent preconditioner control.

For m = 0 with normalisation, the assembled contribution equals the
derivative of the residual term only where the two use the same stencil and
coefficients:

- `scaleFactorJacobian = scaleFactor/(4*rStar)`;
- the Jacobian's face coefficient, `|Sf|` times the `laplacian` scheme's
  (orthogonal) delta coefficient, equals `|Sf|*deltaCoeffs` of
  `orthogonalSnGrad`: true on orthogonal meshes, where both distances are the
  cell-centroid separation, but not on non-orthogonal meshes, where a
  `corrected` scheme uses `1/(n & d)` and adds an explicit correction that the
  residual's `orthogonalSnGrad` does not have;
- the boundary treatment of `p` in the residual flux and in the Jacobian's
  boundary coefficients agree (as for a zero-gradient pressure);
- `gamma` (`rAUf`) does not depend on the unknowns, or its dependence is
  negligible: the Jacobian holds `rAUf` fixed. `rAUf` is built from the
  diagonal of the approximate momentum operator, which depends on the
  implicit stiffness, the mesh and the time step; it varies with the
  solution only where the implicit stiffness does.

Agreement to round-off (4e-16 relative, against a finite-difference Jacobian)
was established for one configuration: `linearGeometryTotalDisplacement`,
`linearElastic`, a 2-D orthogonal 6 x 5 graded Cartesian mesh,
`laplacian` with `scaleFactor 10`, `referenceNyquistDirections 2` and
`scaleFactorJacobian 1.25`, `Gauss linear corrected`, with the mixed system
scaling on. It is not claimed for non-orthogonal meshes or for
solution-dependent implicit stiffness. In every case the assembled matrix of
the mixed solid models remains an approximate preconditioning matrix, used
with a matrix-free Newton operator (`snes_mf_operator`). For m >= 1 the
Laplacian is a surrogate, and `scaleFactor/(4*rStar)` is a starting point
rather than an exact value.

## Usage Pattern

The intended solver-side usage is:

1. Construct the model through `stabilisationModel::New(...)`.
2. Call `updateScalar(...)` or `updateVector(...)` when the stabilisation
   needs to be refreshed.
3. Read the cached result using `faceScalar()` or `faceVector()`.
4. Optionally derive a cached cell field using `cellScalar(...)` or
   `cellVector(...)`.
5. If supported by the model, use `scalarJacobian(...)` or
   `vectorJacobian(...)` to obtain an approximate implicit contribution.

For example, in
[`linGeomTotalDispSolid.C`](https://github.com/solids4foam/solids4foam/blob/development/src/solids4FoamModels/solidModels/linGeomTotalDispSolid/linGeomTotalDispSolid.C),
the momentum stabilisation is updated before its face contribution is added to
the traction.

If `faceScalar()` or `faceVector()` is accessed before the corresponding
`update...()` call, the code intentionally aborts with a clear fatal error.
This is by design and is used to keep the interface simple.

The `cellScalar(...)` and `cellVector(...)` accessors follow the same pattern.
They build and cache a cell-centred divergence of the current face
stabilisation. If the corresponding face field has not been initialised yet,
they also abort with a clear fatal error.

The cell-field interface accepts:

```text
cellScalar(gammaPtr, rebuild)
cellVector(gammaPtr, rebuild)
```

The optional `gammaPtr` lets a solver include a face diffusivity inside the
divergence, for example `div(gamma*faceScalar*magSf)`. The optional `rebuild`
flag defaults to `false` and forces the cached cell field to be recalculated
from the current cached face field when set to `true`.

## Gradient Arguments

The `updateScalar(...)` and `updateVector(...)` functions optionally accept a
pointer to the gradient of the primary field.

This is intentionally a minimal interface. Some models need the gradient to
build the stabilisation efficiently or accurately, while others can reconstruct
what they need internally and therefore ignore the pointer.

Current expectation:

- gradient-based stencil/jump models may require `gradPtr`,
- Laplacian-like models may ignore `gradPtr`,
- callers should pass the available gradient when they already have it.

This avoids forcing more specialised interfaces for each stabilisation family.

## Jacobian Caching

Some models provide an approximate Jacobian through `scalarJacobian(...)` or
`vectorJacobian(...)`. These matrices are cached internally after first
construction.

The Jacobian interface also accepts an optional `rebuild` flag:

```text
scalarJacobian(field, gammaPtr, rebuild)
vectorJacobian(field, gammaPtr, rebuild)
```

The optional `scaleFactorJacobian` entry defaults to `1.0`, not to the
explicit stabilisation `scaleFactor`. This lets the approximate Jacobian act
as an independently tuned implicit contribution when needed.

The default is `false`, which preserves the existing behaviour and reuses the
cached matrix. If a solver needs the approximate Jacobian to be reconstructed
because the relevant coefficients or matrix structure have changed, it can pass
`true`.

For the current development work, the solid-model call sites are left on the
default behaviour, so no existing solver logic changes.

## Combining Multiple Models

The `multiple` type combines two or more stabilisation models at run time. Their
face contributions are summed and an optional outer `scaleFactor` (defaults to
`1.0`) is applied to the total. The Jacobian is delegated to a single
user-nominated sub-model.

```text
stabilisation
{
    momentum
    {
        type            multiple;
        // scaleFactor 1.0;  // optional outer scale; defaults to 1.0
        models          ( stab1  stab2 );
        jacobianModel   stab1;

        stab1
        {
            type        laplacian;
            scaleFactor 0.1;
        }

        stab2
        {
            type        gradDiv;
            scaleFactor 0.05;
        }
    }
}
```

Each entry in `models` must correspond to a named sub-dictionary in the same
dictionary. The `jacobianModel` entry names whichever sub-model provides the
Jacobian (it must be one of the entries in `models`).

The `highOrderResidual` option is allowed with `multiple` only if every
sub-model returns `true` from `supportsHighOrderResidual()`. Currently only the
`alpha` model satisfies this requirement.

## Volumetric Strain Rate Stabilisation

The `volStrainRate` type targets temporal oscillations in the volumetric
strain `tr(gradD)` rather than spatial oscillations.  The face stabilisation
vector is:

```text
faceVector_f =
    scaleFactor * C(mode)
  * [tr(gradD_f_future) - tr(gradD_f_old)] * n_f
```

Five modes are available, selected via the `mode` entry:

| mode | C(mode) | Vanishes as |
| --- | --- | --- |
| `physicalDamping` | `tau/deltaT` | never — permanent physical damping |
| `firstOrderTemporal` | `1` | O(deltaT) |
| `secondOrderTemporal` | `deltaT` | O(deltaT^2) |
| `spatioTemporal` | `h^2/m^2` | O(h^2 * deltaT) |
| `h2PhysicalDamping` | `h^2/(m^2*deltaT)` | O(h^2) as h→0; O(1) as deltaT→0 |

The gradient pointer (`gradPtr`) passed to `updateVector` must be non-null
and must carry an old-time level (i.e. `gradPtr->oldTime()` must be valid).
In `nonLinGeomTotalLagVelocitySolid` this is satisfied by passing `&gradD()`.

```text
stabilisation
{
    momentum
    {
        type            volStrainRate;
        scaleFactor     0.1;
        mode            firstOrderTemporal;
    }
}
```

For `physicalDamping` mode, an additional `tau` entry (in seconds) sets the
damping time scale:

```text
        type            volStrainRate;
        scaleFactor     0.1;
        mode            physicalDamping;
        tau             1e-4;
```

The `h2PhysicalDamping` mode combines the O(h²) spatial scaling of
`spatioTemporal` with the Δt-persistence of `physicalDamping`.  It is the
recommended mode when the instability is a *temporal* pressure pulsation
(whole-domain oscillation) that worsens as Δt decreases, and when spatial
consistency under mesh refinement is also required.  No `tau` parameter is
needed; the time scale is folded into `scaleFactor`.

```text
        type            volStrainRate;
        scaleFactor     0.1;
        mode            h2PhysicalDamping;
```

An approximate Laplacian Jacobian is provided via `vectorJacobian`.  For
modes whose coefficient depends on `deltaT` (`physicalDamping`,
`secondOrderTemporal`, and `h2PhysicalDamping`), the Jacobian should be
rebuilt each time step by calling `vectorJacobian(field, gammaPtr, true)`.

## Notes For Extension

When adding a new stabilisation model:

- derive from `stabilisationModel`,
- add `TypeName("...")` in the header,
- register with `addToRunTimeSelectionTable(...)` in the source,
- add the source to the relevant build lists,
- keep the interface consistent with the existing `update...()` and cached
  accessor pattern unless there is a strong reason to change the framework,
- if the new model is compatible with `highOrderResidual` calculation, override
  `supportsHighOrderResidual()` to return `true`.

The current framework is being validated first with
`linGeomTotalDispSolid`. Wider adoption across the other solid models can
follow once the design settles.
