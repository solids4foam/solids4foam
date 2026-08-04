---
sort: 9
---

# Under the hood: `solids4FoamModels`

This directory holds the sources of `libsolids4FoamModels`, the main
solids4foam library. It is built by the `Allwmake` script in this directory,
which is called in turn by `src/Allwmake`.

The `solids4Foam` application in `applications/solvers/solids4Foam` is a thin
front end: it creates a physics model and advances it in time. Everything else
lives in this library.

---

## Library overview

### Physics models

`physicsModel` is the virtual base class for the three families of solution
procedure, each of which is selected at run time from its own dictionary:

- **Solid**: base class `solidModel`, read from `constant/solidProperties`,
  in `solidModels`;
- **Fluid**: base class `fluidModel`, read from `constant/fluidProperties`,
  in `fluidModels`;
- **Fluid-solid**: base class `fluidSolidInterface`, read from
  `constant/fsiProperties`, in `fluidSolidInterfaces`.

Which family is used for a case is set by the `type` entry in
`constant/physicsProperties`, i.e. `solid`, `fluid` or `fluidSolidInteraction`.

A **solid model** solves the governing momentum equation in a solid domain.
The available models differ in whether the geometry is linear or nonlinear,
whether a total or updated Lagrangian formulation is used, whether the
discretisation is cell-centred or vertex-centred, and whether the solution
algorithm is segregated, coupled or explicit. They are described in
[solid models](https://www.solids4foam.com/documentation/solid-models/).

A **fluid model** solves the flow equations, using either a solver built on
the standard OpenFOAM ones or a solver specific to solids4foam; see
[fluid models](https://www.solids4foam.com/documentation/fluid-models/).

A **fluid-solid interface** owns one fluid model and one solid model and
implements the partitioned coupling algorithm between them, for example fixed
relaxation, Aitken or IQN-ILS. See
[fluid-solid interfaces](https://www.solids4foam.com/documentation/fluid-solid-interfaces/).

### Material behaviour

The constitutive behaviour of a solid is separate from the solid model that
uses it, so the two can be varied independently:

- `materialModels/mechanicalModel` reads `constant/mechanicalProperties` and
  creates the mechanical laws. A `mechanicalLaw` returns the stress for a
  given deformation; the laws are grouped into `linearGeometryLaws` and
  `nonLinearGeometryLaws`. Multiple materials are supported, each region
  taking its own law, with corrections at bi-material interfaces to keep the
  stress continuous without oscillations.
- `materialModels/thermalModel` does the same for thermal properties, reading
  `constant/thermalProperties` and creating `thermalLaw` objects.

### Supporting components

- `physicsModel`: the base class of the three physics model families;
- `solidModels/fvPatchFields` and `solidModels/pointPatchFields`: solid
  boundary conditions, including tractions, contact, symmetry and prescribed
  displacement and rotation;
- `numerics`: discretisation and linear algebra support, including gradient
  schemes, stabilisation models, interpolation, additional tensor types and
  PETSc helpers;
- `higherOrderHelpers`: higher-order finite volume support, including moving
  least squares, quadrature and the associated operators;
- `dynamicFvMesh`: meshes which change topology, such as `crackerFvMesh` for
  crack propagation;
- `functionObjects`: run-time post-processing, and the analytical solutions
  used by the tutorials; see
  [function objects](https://www.solids4foam.com/documentation/function-objects/).

### Adding a class

New classes follow the usual OpenFOAM run-time selection pattern: derive from
the relevant base class, add the `TypeName` and `addToRunTimeSelectionTable`
entries, and add the source to the canonical list described below. Please also
see [the contributing guide](https://www.solids4foam.com/documentation/under-the-hood/contributing.html).

---

## Source lists

The sources are listed in one of two canonical lists, which differ because the
library supports both OpenFOAM and foam-extend:

- `Make/files.openfoam` for OpenFOAM.com and OpenFOAM.org;
- `Make/files.foamextend` for foam-extend.

`Allwmake` selects the appropriate one and links it to `Make/files`, which is
the file `wmake` reads and which is not tracked by Git. A new class is added to
the canonical list for the versions it supports, exactly as in any other
OpenFOAM library.

---

## Unity builds

Compiling the library in the conventional way means each of its ~230 sources is
a separate translation unit, so the large set of OpenFOAM headers is re-parsed
once per source. That parsing, rather than the solids4foam code itself,
dominates the build time.

A "unity build" (also called a "jumbo build") instead groups the sources into a
small number of buckets, where each bucket is a generated source file that
`#include`s the sources assigned to it. The headers are then parsed once per
bucket. The resulting library is the same either way.

To enable it, set `S4F_UNITY_BUILD` before building:

```bash
export S4F_UNITY_BUILD=1
./Allwmake -j 16
```

As a guide, a from-scratch build of the library with `-j 16` on a 192-core
machine takes 119 s conventionally and 41 s as a unity build, a speed-up of
about 2.9.

Unsetting `S4F_UNITY_BUILD` returns to a conventional build; the generated
sources and their objects are discarded automatically.

### Number of buckets

By default, the number of buckets follows the number of compilation processes,
clamped to between 8 and 32. Fewer buckets means less repeated header parsing,
but larger translation units and so a higher peak memory use per compiler
process; more buckets gives diminishing returns in wall-clock time while the
total compilation cost keeps growing.

The default can be overridden with `S4F_UNITY_NBUCKETS`:

```bash
export S4F_UNITY_BUILD=1
export S4F_UNITY_NBUCKETS=24
./Allwmake -j 24
```

### What is generated

`Make/makeUnityFiles` reads the canonical source list, expands its `$(var)`
path variables, and writes:

- `unityBuild/unityNN.C`, the bucket sources;
- `Make/files`, listing the buckets in place of the individual sources.

Both are generated and neither is tracked by Git. The canonical list remains
the single point of maintenance: new classes are added there as before and are
picked up automatically.

### Sharing a translation unit

Sources which share a translation unit also share its file scope, so a macro or
a file-scope name defined in one source is visible to those which follow it in
the same bucket. Macros defined in a `.C` file should therefore be `#undef`ined
at the end of that file, as is done in `numerics/logExpVolFields/eig3/eig3.C`.

Where that is not practical, a source can be listed in `Make/unity-exclude`,
using the path exactly as it appears in the canonical list with any `$(var)`
references expanded. Listed sources are compiled on their own, as they are in a
conventional build.

### When not to use it

Unity builds pay off when building the library from scratch, for example in
continuous integration, in a container image, or after a `wclean`. They are a
poor fit for day-to-day development of the library itself: editing any source
rebuilds its whole bucket rather than the single file.
