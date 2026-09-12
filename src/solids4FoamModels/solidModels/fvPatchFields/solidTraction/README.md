# `solidTraction` boundary condition

`solidTraction` applies a prescribed traction component and optional pressure
to a solid displacement boundary. It is used by the standard cell-centred
solid models and by the high-order total-displacement models.

## Case setup

Use it on a displacement field, normally `D`:

```foam
interface
{
    type            solidTraction;
    traction        uniform (0 0 0);
    pressure        uniform 0;
    value           uniform (0 0 0);
}
```

`traction` and `pressure` may instead be supplied using `tractionSeries` /
`pressureSeries` or `tractionField` / `pressureField`. Only one input form is
allowed for each quantity. `useUndeformedArea`, `nonOrthogonalCorrections`,
`secondOrder`, `extrapolateValue`, and `relaxationFactor` retain their existing
dictionary meanings.

The prescribed traction is the non-pressure component. The effective boundary
load is `traction - n*pressure`, where `n` is the patch face normal.

## Partitioned FSI and preCICE

Existing partitioned FSI interfaces, including `IQNILS`, transfer face-based
tractions. Calling `solidModel::setTraction(...)` stores the face values and
expands each value to all quadrature points on that face. Consequently, current
face-centre fluid solvers can use high-order solid residuals without a case
dictionary change; the traction is piecewise constant over each solid face.

The preCICE `solidForce` boundary condition follows the same path. It reads the
per-face force field, converts it to traction using the appropriate current or
reference area, and supplies constant quadrature traction. Existing preCICE
case setup therefore remains unchanged.

## Direct quadrature input for developers

Higher-order fluid couplers can set the traction component directly with:

```cpp
solid.setTractionQuadrature(tractionPatch, tractionQuadrature);
```

`tractionQuadrature` is a `CompactListList<vector>` in local patch-face order.
For every patch face, its number and order of values must exactly match the
MLS face-quadrature points returned by the solid model. Invalid layouts stop
with a fatal error.

Direct quadrature data is the canonical input for high-order assembly.
`solidTraction` also calculates the face traction as the quadrature-weighted
average, so face and quadrature representations remain consistent. A later
face-based update replaces the direct data with its piecewise-constant
equivalent. Pressure is currently face-based in both paths.

`evaluateQuadrature()` returns the traction component with the face pressure
already applied. Couplers that provide a pressure-inclusive total traction must
separate the pressure component before calling `setTractionQuadrature()`.

## Compatibility

Quadrature storage and direct quadrature input are available for OpenFOAM.com
and OpenFOAM.org builds. High-order MLS residuals are not supported by
foam-extend, where the existing face-based boundary-condition behaviour is
unchanged.
