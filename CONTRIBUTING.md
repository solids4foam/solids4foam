---
sort: 2
---

# Contributing to solids4foam

Thank you for considering contributing to **solids4foam**!

Thanks for helping create a world where fluid and structure are solved using
your favourite CFD code.

By participating in this project, you agree to abide by our
[Code of Conduct](CODE_OF_CONDUCT.md).

## How Can I Contribute?

### 1. Reporting Bugs

If you encounter a bug, please open a
[GitHub issue](https://github.com/solids4foam/solids4foam/issues) and include
the following:

- A concise and descriptive title.

- Steps to reproduce the issue, along with explaining the expected behaviour.

- Relevant details, such as screenshots, logs, or error messages and an
  environment description in the following form:

  ```plaintext
  - OS version: Ubuntu 22.04.5 LTS
  - OpenFOAM version: openfoam2312
  - solids4foam version: 2.1
  - Docker image (if applicable)
  ```

### 2. Suggesting Features

We welcome ideas for new features! To suggest a feature:

- Open a [GitHub issue](https://github.com/solids4foam/solids4foam/issues) with
  the title: **Feature request: feature title**.
- Provide a detailed description of the proposed feature.
- Explain its potential benefits for the solids4foam toolbox.

### 3. Submitting Code Changes

#### Step 1: Fork and Clone

1. Fork the repository by clicking the "Fork" button on GitHub.

2. Clone your fork locally:

   ```bash
   git clone https://github.com/your-username/solids4foam.git
   cd solids4foam
   ```

#### Step 2: Make Your Changes

- Follow the
  [OpenFOAM Coding Style Guide](https://openfoam.org/dev/coding-style-guide/)).
- Ensure your code can be compiled and adheres to the existing code structure
  and organisation.
- Test your changes thoroughly with relevant test cases. Include or update
  tutorial cases if applicable.
- Document your code with clear comments.
- Before submitting, check the code using available linters or style checkers
  ([shellcheck](https://www.shellcheck.net/#),
  [markdownlint](https://dlaa.me/markdownlint/)).

#### Step 3: Commit Your Changes

Write clear and descriptive commit messages and commit your changes. Reference
any related issues using hashtags if applicable::

```bash
git add file_name.C
git commit -m "Fix for issue #51"
```

#### Step 4: Push and Open a Pull Request

Push your branch to your fork:

```bash
git push
```

Go to the original solids4foam repository on GitHub and open a pull request
(target the `development` branch for your pull request). Include:

- **A clear title**, such as `Fix: Issue #123` or
  `Feature: Add new functionality`.
- **A detailed description** of your changes.
- Any **related issue numbers** (e.g., `Closes #123`).
- **Instructions or example cases** to test your contribution, if applicable.

## Reviewing Pull Requests

Every ready PR into `development` goes through the steps below. Core
maintainers (currently Philip Cardiff and Ivan Batistić) approve and merge.
Contributing reviewers are welcome at every risk level.

### 1. Automated review and risk label

When a PR is marked ready, and on every later push, a Claude review is posted
as one comment and updated in place. It labels the PR `ai-risk: low`,
`ai-risk: medium` or `ai-risk: high`. A failed review is labelled
`ai-risk: high` and `ai-review: incomplete`; re-run it or review by hand.
Check that the report's commit SHA matches the PR head before relying on it.
Iterate in draft to avoid a review on every push.

Risk follows the consequences of a change, not its size. Changes to core
numerics, solvers, FSI coupling, CI or agent configuration are at least HIGH
by path. If you think the risk is too low, say so on the PR; the higher risk
applies.

### 2. Merge rules

<!-- markdownlint-disable MD013 -->

| Risk | Second AI review | Human check | Who merges |
| --- | --- | --- | --- |
| LOW | Not needed | A maintainer works through the report's human review checklist | One core maintainer |
| MEDIUM | Required, posted on the PR | Author addresses findings; maintainer checks the checklist | One core maintainer |
| HIGH | Required, posted on the PR | As MEDIUM, plus numerical evidence | Both core maintainers approve |

<!-- markdownlint-enable MD013 -->

At every level you may check things privately with your own AI assistant. Do
not post another full AI review unless it adds something specific.

### 3. Second independent AI review (MEDIUM and HIGH)

A core maintainer runs it once the PR is close to final, and again only after
substantive changes. The reviewer gets the same instructions and patches as
the Claude report, but not the report or any PR comment, so the two are
independent. From an up-to-date `development` checkout:

```bash
# Preferred: Codex (gpt-6-astra, xhigh), a different model family
node .github/scripts/secondReview.cjs <PR number> --post

# Substitute: Claude Fable 5.2 (high); needs ANTHROPIC_API_KEY and
# npm ci --prefix .github/scripts --ignore-scripts
node .github/scripts/secondReview.cjs <PR number> --model fable --post
```

Fable shares the first report's model family, so agreement is weaker
evidence. It is fine for MEDIUM. For HIGH, use it only if a core maintainer
judges the classification borderline and says so on the PR. Findings from
either review are addressed by the author, with reviewers' help, before merge.
A report whose SHA no longer matches the PR head is historical only.

### Contributing reviewers

You can check a numerical method, reproduce results, or explain a
compatibility concern. State what you inspected and what remains uncertain. On
LOW and MEDIUM PRs you may do the human-check pass; a core maintainer still
merges. HIGH PRs still need both core maintainers.

Contributors develop into contributing reviewers through sustained work and
paired reviews. Progression to core maintainer depends on review quality,
architectural understanding, and sound technical judgement. Core maintainers
discuss that responsibility with the contributor; review counts do not grant
it.

The AI reviews are advisory. They can be wrong or miss bugs, and they never
replace human judgement or numerical verification. See the
[review guide](.github/review/README.md) for the workflow's security, costs
and setup.

## Questions?

If you have any questions or need assistance, we encourage you to use the
appropriate discussion platform based on your query:

- **User Support**: For general usage questions, user support, or
  troubleshooting, please post on the **CFD Online Forum**:
  [https://www.cfd-online.com/Forums/openfoam-cc-toolkits-fluid-structure-interaction/](https://www.cfd-online.com/Forums/openfoam-cc-toolkits-fluid-structure-interaction/)

- **Developer Discussions**: For discussions about development, contributions,
  or feature requests, please open a discussion on our **GitHub Discussions**
  page:
  <https://github.com/solids4foam/solids4foam/discussions>](<https://github.com/solids4foam/solids4foam/discussions>)

Thank you for contributing to solids4foam!
