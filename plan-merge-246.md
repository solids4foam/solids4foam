# Plan to Resolve PR 246 Merge Conflicts

PR: https://github.com/solids4foam/solids4foam/pull/246

## Goal

Resolve PR 246 against current `development` while keeping the stabilisation
changes aligned with the current base-class API.

## Plan

1. Resolve `stabilisationModel.H` by rejecting the PR's new
   `usesAdaptiveScalarGamma()` / `adaptiveScalarGamma()` API. Current
   `development` already routes problem-dependent pressure stabilisation
   through the existing `cellScalar(..., gammaPtr)` and
   `scalarJacobian(..., gammaPtr)` interfaces.

2. Resolve `linGeomTotalDispSolid` by using the `solidModel::rAUf()` field now
   present on `development`. Do not keep the PR's duplicate `rAUfPtr_`,
   `makeRAUf()`, or `rAUf()` declarations/definitions in
   `linGeomTotalDispSolid`.

3. Keep the current base sign convention:

   ```cpp
   rAUf() = -1.0/(fvc::interpolate(approxMomJ.A())*one);
   ```

   and keep the pressure stabilisation terms as:

   ```cpp
   + pressureStabilisation().cellScalar(&rAUf(), true)*one
   + one*pressureStabilisation().scalarJacobian(p, &rAUf())
   ```

   This matches current `development` and `nonLinGeomTotalLagTotalDispSolid`.

4. Keep the PR's useful `linGeomTotalDispSolid` change that recalculates the
   face traction after pressure is applied to `sigma()`.

5. Keep the PR's `JamesonSchmidtTurkel` and
   `generalisedEvenOrderLaplacian` `scalarJacobian` / `vectorJacobian`
   additions. These fit the current `stabilisationModel` Jacobian interface.

6. Resolve `GuccioneElastic.C` by keeping `development`'s OpenFOAM-version
   compatible `typeHeaderOk/headerOk` checks and
   `io.readOpt() = IOobject::MUST_READ`, while preserving the PR's intended
   fibre-field fallback/interpolation and tensor-index fixes.

7. Verify the result by building `libsolids4FoamModels` for v2512:

   ```sh
   source /usr/lib/openfoam/openfoam2512/etc/bashrc
   cd src/solids4FoamModels
   wmake libso
   ```

8. If v2512 builds, repeat focused compatibility builds for
   `../solids4foam.of9` and `../solids4foam.fe41` where practical.

9. For behavioural validation, rerun the Cook's membrane PETScSNES
   pressure-stabilisation checks for `laplacian`, `JST`, and
   `generalisedEvenOrderLaplacian` with `laplacianPower 0/1`.
