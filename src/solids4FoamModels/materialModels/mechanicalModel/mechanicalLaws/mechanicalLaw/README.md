---
sort: 1
---

# mechanicalLaw: the base class

`mechanicalLaw` is the abstract base class for every constitutive law in
solids4foam. A law's job is narrow: given the current deformation, return the
stress, plus the handful of scalars the solid models need in order to build
and solve the momentum equation. Everything else — reading
`mechanicalProperties`, partitioning the mesh by material, mapping fields
between sub-meshes — is done by `mechanicalModel` and is invisible to the law.

Source:
[`mechanicalLaw.H`](https://github.com/solids4foam/solids4foam/blob/development/src/solids4FoamModels/materialModels/mechanicalModel/mechanicalLaws/mechanicalLaw/mechanicalLaw.H),
[`mechanicalLaw.C`](https://github.com/solids4foam/solids4foam/blob/development/src/solids4FoamModels/materialModels/mechanicalModel/mechanicalLaws/mechanicalLaw/mechanicalLaw.C)
and
[`newMechanicalLaw.C`](https://github.com/solids4foam/solids4foam/blob/development/src/solids4FoamModels/materialModels/mechanicalModel/mechanicalLaws/mechanicalLaw/newMechanicalLaw.C).

---

## User Guide

This is a developer-facing page. As a user you never select `mechanicalLaw`
itself; you select one of its derived types with the `type` keyword inside a
`mechanical` sub-dictionary. The only part of this class you interact with
directly is the small set of entries it reads for every law:

| Entry | Default | Description |
| --- | --- | --- |
| `rho` | none | Density, as a `dimensionedScalar` |
| `solvePressureEqn` | `no` | Smooth the hydrostatic stress with a Laplacian |
| `pressureSmoothingScaleFactor` | `100` | Scale factor for that smoothing |
| `regionName` | see below | Name of the base solid mesh region |

`rho` is read only when the law uses the base-class `rho()`/`rhoScalar()`,
which nearly all do. `solvePressureEqn` and `pressureSmoothingScaleFactor` are
read with `lookupOrAddDefault`, so their values are echoed into
`constant/mechanicalProperties.withDefaultValues`. `regionName` defaults to
`solid` if a `solid` region exists, otherwise `region0`; set it only in
unusual multi-region set-ups, such as a law attached to a mesh-motion mesh.

Enabling `solvePressureEqn` solves an additional Laplacian equation for the
hydrostatic stress. It can quell checkerboarding in the hydrostatic stress at
the cost of an extra linear solve per outer iteration; it is off by default.

---

## Developer Notes

### The interface a law must implement

Only two functions are pure virtual:

| Function | Contract |
| --- | --- |
| `correct(volSymmTensorField& sigma)` | Set `sigma` from the deformation |
| `impK()` | Return the implicit stiffness field |

Everything else has a default implementation, which is either a sensible
fallback or a `notImplemented` that fires only if a solid model actually
requests it.

**`correct(volSymmTensorField& sigma)`** is the constitutive update. It must
overwrite `sigma` on the mesh the law was constructed on — which, in a
multi-material case, is the material's sub-mesh, not the base mesh. Whether
`sigma` is the Cauchy stress or a work-conjugate measure is fixed by the solid
model the law is registered for: linear-geometry and total-Lagrangian laws in
solids4foam return the Cauchy stress, and it is the law's responsibility to
perform any push-forward internally. `correct()` is called every outer
iteration, so it must be idempotent in the sense that calling it twice with
the same deformation gives the same answer; any path-dependent state (plastic
strain, viscous internal variables) must be advanced in `updateTotalFields()`,
not here.

**`impK()`** returns the diffusivity used for the implicit Laplacian term in
the segregated momentum equation. It is not a physical property: it is a
numerical stiffness that only has to make the outer iterations converge, since
the discretisation is written as an implicit Laplacian plus an explicit
deferred correction that cancels it at convergence. In practice laws return
something of the order of `2*mu + lambda`, and it is legitimate to return a
larger value to damp the iterations at the cost of slower convergence. Two
things do depend on `impK()` beyond convergence rate, so it cannot be
arbitrary:

- the traction boundary conditions use it to convert a traction into a
  surface-normal gradient;
- the bi-material interface correction in `solidSubMeshes` uses the base-mesh
  `impK` field to weight the two sides of an interface, so a badly scaled
  `impK` degrades the interface treatment in multi-material cases.

`impKf()` is the surface-field equivalent. The base class implements it as
`fvc::interpolate(impK())`; override it only if the law can produce a genuinely
better face value.

**`rho()`** returns a `volScalarField` built from `rhoScalar()`, which reads
the `rho` entry. The field is created with `READ_IF_PRESENT`, so a `rho` field
in the time directory takes precedence over the dictionary value, and with
`zeroGradient` boundaries which are corrected before the field is returned.
Override `rho()` for a genuinely non-uniform density, or `rhoScalar()` alone
if the law computes a uniform density from other entries.

**`bulkModulus()`** and **`shearModulus()`** default to `notImplemented`. They
are not needed by every solid model, but the ones that do call them will
abort if the law has not overridden them, so implement them for any new law
intended for general use. `bulkModulus()` is what the compressibility-related
machinery uses; a nearly incompressible law must return the correct value here
or the mixed pressure-displacement solid models cannot be used with it.

### Optional interface

| Function | Default | Purpose |
| --- | --- | --- |
| `correct(surfaceSymmTensorField&)` | `notImplemented` | The `uns` models |
| `correct(pointSymmTensorField&, ...)` | `notImplemented` | Vertex-centred |
| `correct(CompactListList<...>)` | `notImplemented` | Face quadrature |
| `materialTangent()` | `notImplemented` | Newton-Raphson (PETSc SNES) |
| `materialTangentField(...)` | uniform | Per-face material tangent |
| `residual()` | `0` | Material convergence measure |
| `newDeltaT()` | `endTime` | Law-requested time-step limit |
| `updateTotalFields()` | no-op | End-of-time-step state update |
| `setRestart()` | no-op | Adjust field write options on restart |

A few notes on these. The surface-field `correct()` must be implemented for a
law to work with the `uns` family of solid models, and the point-field variant
for the vertex-centred ones; the face-quadrature variant is not available on
foam-extend at all. `materialTangentField()` defaults to filling
`mesh().nFaces()` entries with the single value from `materialTangent()`, so
overriding just `materialTangent()` is enough for a homogeneous law.
`residual()` returning a non-zero value makes the solid model keep iterating
until the material has converged as well as the momentum equation — this is
what plasticity laws use — and `mechanicalModel::residual()` takes the maximum
over all laws. `newDeltaT()` is likewise reduced to the minimum over all laws.

### Total- versus updated-Lagrangian variants

A law is registered in exactly one of the two run-time selection tables:
`linGeomMechLaw` or `nonLinGeomMechLaw`. Both tables have the same constructor
signature, and the `nonLinearGeometry::nonLinearType` enumerator passed to it
is what tells a nonlinear law which formulation it is being used in. It is
available to the law through the protected `nonLinGeom()` accessor and takes
one of `LINEAR_GEOMETRY`, `TOTAL_LAGRANGIAN` or `UPDATED_LAGRANGIAN`.

The base class provides both deformation gradients so that a single law class
can serve both formulations:

- `F()` and `Ff()`: the total deformation gradient, relative to the initial
  configuration;
- `relF()` and `relFf()`: the relative deformation gradient, relative to the
  configuration at the end of the previous time step.

The `updateF(...)` family of protected functions is the intended entry point.
It updates the appropriate gradient for the current formulation from the
registered `grad(D)`/`grad(DD)` fields, and it also implements the
*enforce-linear* escape hatch: if the solid model has tripped its
`enforceLinear` switch — typically because the outer iterations are diverging
— `updateF` overwrites `sigma` with the linear elastic response built from the
linearised `mu` and `K` passed in, and returns `true`. A law should check that
return value and skip its own constitutive update:

```text
if (updateF(sigma, mu_, K_))
{
    return;
}
```

This is why laws must call `mu(...)` and `K(...)` to register their linearised
moduli even when the constitutive model itself has no such constants: those
fields are the fallback.

Whether the law is used incrementally is a separate question from the
formulation, and is answered by the protected `incremental()` accessor, which
asks the solid model whether it solves for `DD` rather than `D`. The base
class uses this to decide which gradient field to look up, so a law that goes
through `updateF()` does not need to care.

### Other base-class services

`mechanicalLaw` also manages a number of demand-driven fields so that derived
laws do not each reimplement them: `mu`/`muf` and `K`/`Kf` (linearised shear
and bulk moduli), `epsilon`/`epsilonf` (small strain, via `updateEpsilon()`),
`sigma0`/`sigma0f` (initial or residual stress), and
`sigmaHyd`/`gradSigmaHyd` with the `updateSigmaHyd(...)` overloads that
implement the optional pressure smoothing equation. The accessors are
protected and each field is only allocated when first used.

`planeStress()` returns the `planeStress` switch from `mechanicalProperties`,
searching the law's own region first and then the `region0` or `solid`
sub-registry, so it works both for standalone solid cases and for the solid
region of a fluid-solid interaction case.

### Adding a new mechanical law

1. **Create the directory and class.** Place it under
   `src/solids4FoamModels/materialModels/mechanicalModel/mechanicalLaws/`,
   in `linearGeometryLaws/` or `nonLinearGeometryLaws/` according to which
   selection table it belongs to. Copying `linearElastic` or
   `neoHookeanElastic` is the fastest start.

2. **Derive from `mechanicalLaw`** and declare the run-time type name in the
   header:

   ```text
   TypeName("myElasticLaw");
   ```

   The string is what users write as `type` in `mechanicalProperties`; it does
   not have to match the class name, though by convention it does. The one
   deliberate exception in the toolbox is the thermal law `constantThermal`,
   whose type name is `constant`.

3. **Provide the standard constructor**, taking
   `(const word& name, const fvMesh& mesh, const dictionary& dict,
   const nonLinearGeometry::nonLinearType& nonLinGeom)`, and forward all four
   arguments to the base class. Read the law's own entries from `dict` in the
   initialiser list, using `lookup` for required entries and
   `lookupOrDefault`/`lookupOrAddDefault` for optional ones. Disallow copy
   construction and assignment, as the base class does.

4. **Implement `correct(volSymmTensorField&)` and `impK()`**, and override
   `bulkModulus()` and `shearModulus()` unless the law is only ever going to
   be used with solid models that do not need them.

5. **Register the law** in the `.C` file:

   ```text
   #include "addToRunTimeSelectionTable.H"

   namespace Foam
   {
       defineTypeNameAndDebug(myElasticLaw, 0);
       addToRunTimeSelectionTable
       (
           mechanicalLaw, myElasticLaw, linGeomMechLaw
       );
   }
   ```

   Use `nonLinGeomMechLaw` instead for a finite-strain law. Registering in the
   wrong table is the common mistake; the selector's error message will then
   tell users the law "can only be used with a linear geometry solid model".

6. **Add the source to the build.** There is a single `Make/files` per build
   flavour in `src/solids4FoamModels/Make/`: add one line to **both**
   `files.openfoam` and `files.foamextend`, next to the sibling laws, e.g.

   ```text
   $(linGeomLaws)/myElasticLaw/myElasticLaw.C
   ```

   Header files are not listed. If the law cannot compile on foam-extend,
   omit it from `files.foamextend` rather than guarding the whole file.

7. **Rebuild** with `./Allwmake` in the repository root. Because the laws are
   pulled in by the run-time selection table rather than by any explicit
   reference, a law that is missing from `Make/files` will compile in a unity
   build but never appear in the table — the symptom is "Unknown mechanicalLaw
   type" listing every law except yours.

8. **Add a tutorial and a `README.md`** in the law's directory, following the
   pattern of the existing law pages.
