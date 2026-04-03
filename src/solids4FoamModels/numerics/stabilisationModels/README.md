# Stabilisation Models

This directory contains runtime-selectable stabilisation models used by
`solids4foam` to add artificial dissipation for suppressing grid-scale
oscillations in scalar and vector fields.

## Design Overview

The common base class is
[`stabilisationModel`](./stabilisationModel/stabilisationModel.H).

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

`RhieChow` is implemented as a helper under
[`diffStencilLaplacianStab/RhieChowStab`](./diffStencilLaplacianStab/RhieChowStab).

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
[`linGeomTotalDispSolid.C`](../../solidModels/linGeomTotalDispSolid/linGeomTotalDispSolid.C),
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
