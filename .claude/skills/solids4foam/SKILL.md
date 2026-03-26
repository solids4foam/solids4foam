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

## 2) Coding Style Rules

### C++ style

- Follow existing OpenFOAM/solids4foam style in surrounding files.
- Use the same indentation, brace style, comment style, and naming conventions as local code.
- Keep lines and expressions readable; avoid clever/condensed code.
- Prefer explicit, local, maintainable changes over abstraction-heavy rewrites.
- Add comments only when behavior is non-obvious; do not add redundant comments.

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

Before finalizing, verify:

- Build impact is localized.
- No unrelated files changed.
- No compatibility regressions introduced by style-only edits.

## 6) Change Delivery Format (Mandatory)

All code changes must be returned as **git patches**.

- Provide changes in patch form (`git diff`/unified diff) suitable for application.
- Keep patches focused and reviewable.
- If multiple concerns are required, split into logical commits/patches.
- Do not provide only prose summaries when code edits were requested.

## 7) Practical Workflow for Codex

1. Read nearby code and follow local patterns.
2. Implement minimal patch.
3. Update runtime registration/build lists if adding a new class.
4. Run or describe relevant checks/tests.
5. Return changes as patch-oriented output.

## 8) What to Avoid

- Large-scale refactors without explicit request.
- API redesigns when a local fix is sufficient.
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
