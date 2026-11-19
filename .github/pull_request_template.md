# solids4foam pull request

## Problem, solution, and verification

<!-- Use short sentences. Explain terms that other developers may not know.
Link related issues. Keep this summary useful without the technical details. -->

**Problem:**

**Solution:**

**Verification:** State what you tested and the result. Say what is untested.

Closes #<!-- issue number, or delete this line -->

## Before review

<!-- Reviewers, human and AI, check these. Tick them or say why they do not
apply. Keep the PR in draft until they are done. -->

- [ ] User-facing changes have an entry under `[Unreleased]` in `CHANGELOG.md`.
- [ ] New or renamed files in `src/solids4FoamModels` are listed in
  `Make/files.openfoam` and/or `Make/files.foamextend`, for the variants they
  support.
- [ ] Added and changed lines follow the
  [OpenFOAM coding style](https://openfoam.org/dev/coding-style-guide/).

## Technical details

<!-- Keep only the details that apply. Use “Not applicable” where needed. -->

- **Implementation:** Explain the main changes and why you chose this approach.
- **Numerical considerations:** State assumptions and effects on accuracy,
  stability, conservation, and convergence. Link the method or reference.
- **Evidence:** Give case names, commands, OpenFOAM variants and versions,
  reference results, and error tolerances. Link logs or plots where useful.
- **Compatibility:** Describe changes to dictionaries, defaults, dependencies,
  existing cases, or solver behaviour.
- **Review focus:** Identify uncertainties and decisions that need human review.

See the [contribution guide][contributing] for review roles and expectations.

[contributing]: https://github.com/solids4foam/solids4foam/blob/development/CONTRIBUTING.md
