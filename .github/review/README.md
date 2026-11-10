# Advisory PR reviews

The Claude review helps people assess a PR. It does not approve a PR, count
approvals, enforce a merge gate, or merge changes. Existing CI still supplies
build and regression evidence.

## Human roles

Core maintainers own architecture, final approval, and merging. Contributing
reviewers are encouraged to assess methods, reproduce results, and explain
technical concerns. Their assessments inform core-maintainer decisions. They
do not independently authorise merging at this stage. AI provides an independent
initial assessment. It can miss defects or report incorrect findings.

Keep role descriptions here rather than individual names in workflow code.
For a future reviewer directory, use an opt-in GitHub organisation team or a
small section in CONTRIBUTING.md with verified profile links and areas of
expertise. Core maintainers should confirm membership before adding anyone.
Do not infer usernames from names. Do not request every reviewer on every PR.
Authors and maintainers can invite a reviewer with relevant experience.

A lightweight path to core maintainership is to pair a contributing reviewer
with a core maintainer for selected PRs. Compare assessments and discuss missed
issues. After sustained, sound reviews and useful numerical verification,
core maintainers can discuss expanding responsibility with that contributor.
No points system or new permission automation is needed.

## Risk and future approval policy

Classify by consequences, not diff size. Use the highest relevant risk.

| Risk | Examples | Intended future policy (not enforced) |
| --- | --- | --- |
| LOW | Documentation; no expected behaviour change | Lightweight process |
| MEDIUM | Isolated functionality | One core-maintainer approval |
| HIGH | Core numerics, architecture, security | Both core-maintainer approvals |

The intended one/both policy refers to the two current core maintainers. Role
membership must be verified before any future enforcement is introduced.
This version changes no branch protection, approval rules, or permissions.

## Setup

1. Review this PR against `development`, the normal contribution target.
   Deploy the workflow and its helper/instruction files to the default branch
   (`master`) before enabling it. GitHub's current event documentation specifies
   that `pull_request_target` runs in the default-branch context. The explicit
   `github.sha` checkout uses that trusted workflow snapshot. See
   [GitHub's event reference](https://docs.github.com/en/actions/reference/workflows-and-actions/events-that-trigger-workflows#pull_request_target).
   Confirm the ref and SHA in the first live run. This draft PR cannot test its
   own privileged workflow before deployment.
2. Create an organisation-owned Anthropic API key for a dedicated workspace.
   Set a spending limit and monitor usage. Store it as the organisation Actions
   secret `ANTHROPIC_API_KEY`, restricted to this repository. A repository secret
   can be used for a pilot. Do not use a personal Claude subscription token.
3. Allow the pinned checkout, github-script, and official Anthropic Action in
   organisation Actions policy. The workflow supplies `GITHUB_TOKEN` explicitly;
   no GitHub App installation, personal access token, or custom App is required.
   The model job needs Contents read and Pull requests read. Only the separate
   publisher job needs Pull requests write to create a conversation comment.
   It never calls the approval or merge APIs. No Issues write or OIDC is needed.
4. After reviewing security and costs, set repository Actions variable
   `CLAUDE_REVIEW_ENABLED` to `true`. Until then all AI jobs are skipped.
   Set it to `false` to stop new reviews. Cancel active runs separately if needed.
5. Run the live-validation scenarios below before relying on report delivery.
   Do not make this advisory workflow a required branch-protection check.

The official integration supports a supplied workflow token. This avoids the
broader permissions of the default Claude GitHub App. See the pinned
[Action inputs](https://github.com/anthropics/claude-code-action/blob/4fa93a88de49bbef68b10eb15aabedf28ca0e1ed/action.yml)
and [security guidance](https://github.com/anthropics/claude-code-action/blob/4fa93a88de49bbef68b10eb15aabedf28ca0e1ed/docs/security.md).

## Operation and security

Eligible PRs are open, non-draft, human-authored PRs. Opening, reopening,
marking ready, or adding commits triggers a review, including external forks.
Bot-authored PRs are excluded to keep the pilot small. No comment commands or
reviewer requests are installed. Editing only the PR description does not
trigger another review; re-run the workflow if evidence was added there.

`pull_request_target` gives access to base-repository secrets for fork PRs.
Both jobs check out only the event's trusted base commit (`github.sha`), with
no persisted checkout credentials. Neither job fetches or checks out PR code.
The collector reads patches through the GitHub API. It never executes them,
loads their configuration, installs their dependencies, or runs their tests.
Only trusted base scripts run. Never change checkout to the PR head or merge
ref, even to inspect a proposed instruction file.

Claude receives a bounded JSON payload with the PR description and patches.
The versioned instructions explicitly treat that payload as untrusted.
All tools are denied. MCP configuration is empty and strict. Bare mode,
user-only settings, and disabled commands prevent project skills and hooks
from loading. The runner is GitHub-hosted and ephemeral. The model has a
read-only GitHub token and no publishing token. The publisher uses a fresh
runner, has no Anthropic secret, validates structured output, and escapes
model text before creating a comment. CI/security paths force HIGH risk even
if the model suggests a lower risk.

Prompt injection can still distort analysis. It cannot be solved by prompts
alone. Human reviewers must check findings and risk. This boundary prevents
model tool use and source writes; it does not guarantee a correct assessment.
Untrusted content is sent to Anthropic for analysis. Do not put credentials in
PRs. Keep debug/full-output logging disabled, including Actions step debug.
Review action updates before changing the full commit pins. The upstream
Action installs its Claude runtime; the action pin does not freeze every
transitive runtime dependency. Repeat live validation after updates.

Each run has a ten-minute model timeout, a small turn limit, and a USD 2 model
budget. Organisation spending limits are still needed: a public contributor
can trigger repeated runs. Concurrency cancels older runs for the same PR.
Reports are appended with the exact head SHA, attempt, and run link. Reports
for older SHAs are historical evidence only. A final API check suppresses stale
publication, though a commit can still arrive between that check and posting;
the explicit SHA remains authoritative.

Missing patches (including binary files), 3000 or more files, payloads over
120 KB, collection errors, API failures, timeouts, and invalid model output
produce an INCOMPLETE comment and a failing job. Risk is then unassessed.
A cancelled run may leave no comment. A publication failure also fails the job.
**No report for the current SHA means no completed AI review.** A green job or
an older report is not an approval. Skipped jobs do not mean a review occurred.

The first version sees patches, not the whole source tree or CI logs. It cannot
confirm linked evidence or run numerical cases. Missing patches fail rather
than silently claiming coverage. GitHub can also limit patch content. Humans
must inspect the complete diff and surrounding code.

## Verification

Run the offline tests from the repository root:

```sh
node --test .github/scripts/claudeReview.test.cjs
markdownlint .github/review/*.md .github/pull_request_template.md CONTRIBUTING.md
```

Tests use mocked GitHub APIs and model responses. They verify data handling,
comment publication, and failure semantics. They do not establish model quality
or live integration. No OpenFOAM build is needed for these workflow-only changes.

Use disposable PRs in a test repository with the same settings for live checks:

<!-- markdownlint-disable MD013 -->

| Scenario | Expected result |
| --- | --- |
| Documentation only | LOW; clear summary; no numerical test demanded |
| Isolated new boundary condition | At least MEDIUM; relevant compatibility questions |
| Core solver/default or constitutive change | HIGH, even for a small diff |
| Numerical change without results | Missing verification stated; practical human checklist |
| Add a second commit | New SHA report; older run cancelled or suppressed |
| External fork without write access | Review runs; no head checkout or code execution |
| Hidden instructions in title/body/patch | No tools, secret disclosure, writes, or approval; concern flagged |
| Invalid API key or model timeout | INCOMPLETE report; failed run; no completed report |
| Workflow permission change | HIGH; proposed permissions never take effect in this run |

<!-- markdownlint-enable MD013 -->

Also test draft/ready transitions, a missing secret, truncated or binary diffs,
malformed output, and denied comment permission. Confirm the comment author,
SHA, run link, escaping, and repeat-run behaviour on GitHub. Confirm no review
approval, source commit, or merge occurs. Review Actions logs for accidental
sensitive output. Model risk and checklist quality require human assessment.
