# vertexCentredNonLinTotalLagGeometry

**This solid model is currently disabled and does not compile into the
library.** See the commented-out entry in `src/solids4FoamModels/Make/files.*`.

## Why

No tutorial in the repository selects it, so nothing has ever exercised it. An
attempt to give it a regression case, on a hyperelastic cantilever beam with a
`neoHookeanElastic` material, did not succeed:

- `useGeometricStiffness` and `compactImplicitStencil` are looked up with no
  default, and no case sets either.
- The material tangent path needs `tangentEps` in the material sub-dictionary,
  which no tutorial in the repository sets. It is an absolute perturbation, and
  so has the scale-dependence problem that the equivalent in the mechanical
  constitutive law framework was fixed for.
- The scalar Jacobian arm additionally needs an `interpolate(impK)` entry in
  `system/dualMesh/fvSchemes`, which no case provides.

With all of those supplied, the run still fails with a floating point exception
inside the PETSc solve, on the first Newton step.

## What would re-enable it

A case that runs it, and whatever fix the floating point exception turns out to
need. This is independent of the mechanical constitutive law framework: the
model was already in this state before that work began, and the framework
neither caused nor depends on it.

Its linear-geometry sibling, `vertexCentredLinearGeometry`, is unaffected,
compiles, and is covered by the `cantilever2d` tutorial.
