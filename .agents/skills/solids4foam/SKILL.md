---
name: solids4foam
description: >
  Use this skill for any task involving the solids4foam repository: adding new
  classes, editing solvers, modifying build files, writing patches, following
  OpenFOAM runtime selection table patterns, or maintaining coding style.
  Trigger whenever the user mentions solids4foam, OpenFOAM, mechanicalLaw,
  solidModel, fluidModel, wmake, or requests C++ changes in this codebase.
---

# Codex Guidelines for solids4foam

This document defines how automated coding changes should be made in this repository.

## 1) Core Principles

- Make the smallest correct change.
- Preserve existing architecture and naming patterns.
- Prefer consistency with nearby code over introducing new style variants.
- Avoid broad refactors unless explicitly requested.
- Keep compatibility across supported OpenFOAM variants in mind.
- State assumptions before implementation when behavior, target OpenFOAM
  variant, or expected numerical result is unclear.
- If multiple interpretations are plausible, ask or present the tradeoff instead
  of silently choosing.
- Push back on requests that would cause unnecessary redesign, broad
  compatibility risk, or speculative behavior.

## 2) Coding Style Rules

### C++ style

- Follow existing OpenFOAM/solids4foam style in surrounding files.
- Use the same indentation, brace style, comment style, and naming conventions
  as local code.
- Keep lines and expressions readable; avoid clever/condensed code.
- Prefer explicit, local, maintainable changes over abstraction-heavy rewrites.
- Add explanatory comments only when behavior is non-obvious; do not add
  redundant explanatory comments. This restriction does not apply to mandatory
  OpenFOAM structural comments and `//-` API documentation comments, which must
  be preserved and reproduced.

### OpenFOAM structural comments

- Preserve and reproduce the standard OpenFOAM/solids4foam structural comment
  layout when creating or editing C++ classes in both `.H` and `.C` files.
- These comments are mandatory style, not optional explanatory comments.
- Do not remove existing OpenFOAM section separators, class declaration banners,
  namespace separators, end-of-file banners, or `//-` API documentation comments
  when editing a class unless the user explicitly asks for style cleanup.
- When creating a new `.H`/`.C` class pair, copy the structural comment pattern
  from the closest existing class in the same directory or model family.

Header files should include the applicable standard banners:

```cpp
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class myClass Declaration
\*---------------------------------------------------------------------------*/

// class declaration here

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
```

Source files should include the applicable standard banners:

```cpp
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// ************************************************************************* //
```

- Header declarations should use OpenFOAM `//-` documentation comments for
  private data, constructors, runtime type information, and public member
  functions, following nearby classes.
- Do not add empty section banners for sections that do not exist, but do add
  the normal banner whenever a section exists.

### File/header conventions

- Preserve existing license/header block format in C++ files.
- Keep include ordering consistent with nearby files.
- Do not change copyright headers unless explicitly asked.

### Scripts and docs

- Match existing shell script style in `Allwmake`, `Allrun`, `Alltest`, etc.
- Keep Markdown concise, practical, and repository-specific.
- All Markdown files must pass the Lint checker, e.g. `markdownlint README.md`

## 3) OpenFOAM Conventions to Follow

- Respect OpenFOAM runtime type conventions:
  - `TypeName("...")` in class headers.
  - Registration via `addToRunTimeSelectionTable(...)` in source files.
- Preserve dictionary-driven behavior and runtime configurability.
- Keep OpenFOAM flavor compatibility intact:
  - OpenFOAM.com
  - OpenFOAM.org
  - foam-extend
- When touching build lists, update both where relevant:
  - `Make/files.openfoam`
  - `Make/files.foamextend`
- Avoid introducing dependencies that break existing wmake workflows.

## 4) Runtime Selection Tables (How They Work)

solids4foam relies heavily on OpenFOAM runtime selection tables to instantiate
models from dictionaries at runtime.

Typical pattern:

1. Base class declares/defines a selection table.
2. Derived class provides `TypeName("...")`.
3. Derived class registers itself with `addToRunTimeSelectionTable(...)`.
4. User selects the type in case dictionaries.

Key base tables in this repository include:

- `physicsModel`
- `fluidModel`
- `solidModel`
- `fluidSolidInterface`
- `interfaceToInterfaceMapping`
- `mechanicalLaw` (linear/nonlinear geometry variants)
- `thermalLaw`
- contact and dynamic-mesh related tables

Rules when adding new runtime-selectable classes:

- Add `TypeName("...")` in header.
- Add `addToRunTimeSelectionTable(...)` in source.
- Add source file to relevant `Make/files*`.
- Ensure dictionary `type` string exactly matches `TypeName` value.

## 5) Minimise Changes

- Only modify files necessary for the requested task.
- Do not reformat unrelated code.
- Do not rename symbols/files unless required.
- Do not alter behavior outside requested scope.
- Prefer targeted edits over cleanup passes.
- Every changed line should trace directly to the user's request, a required
  build fix, or a test/documentation update needed to verify the change.
- If unrelated dead code or cleanup is noticed, mention it separately; do not
  remove it unless asked.
- Remove only unused code, imports, includes, or variables introduced by the
  current change.

Before finalizing, verify:

- Build impact is localized.
- No unrelated files changed.
- No compatibility regressions introduced by style-only edits.

## 6) Change Delivery Format

- Keep changes patch-oriented and reviewable.
- When asked for a patch, provide unified diff output suitable for application.
- When editing directly in the workspace, summarize changed files and
  verification performed.
- Keep patches focused and reviewable.
- If multiple concerns are required, split into logical commits/patches.

## 7) Practical Workflow for Codex

1. Define the concrete success criteria for the request.
2. Read nearby code and follow local patterns.
3. Implement the minimal patch.
4. Update runtime registration/build lists if adding a new class.
5. Verify with the narrowest relevant build/test/check.
6. Report what was verified and what was not.

## 8) What to Avoid

- Large-scale refactors without explicit request.
- API redesigns when a local fix is sufficient.
- New dictionary options, runtime switches, fallback behavior, or
  configurability unless requested or required by existing repository patterns.
- Introducing new style conventions inconsistent with repository norms.
- Silent behavioral changes not documented in the patch summary.

## 9) Reference Files (Commonly Relevant)

- Main solver: `applications/solvers/solids4Foam/solids4Foam.C`
- Core model base classes:
  - `src/solids4FoamModels/physicsModel/physicsModel.H`
  - `src/solids4FoamModels/fluidModels/fluidModel/fluidModel.H`
  - `src/solids4FoamModels/solidModels/solidModel/solidModel.H`
  - `src/solids4FoamModels/fluidSolidInterfaces/fluidSolidInterface/fluidSolidInterface.H`
- Build lists:
  - `src/solids4FoamModels/Make/files.openfoam`
  - `src/solids4FoamModels/Make/files.foamextend`
- Build scripts:
  - `Allwmake`, `src/Allwmake`, `applications/Allwmake`
