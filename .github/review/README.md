# Advisory PR reviews

The Claude review helps people assess a PR. It does not approve a PR, count
approvals, enforce a merge gate, or merge changes. Existing CI still supplies
build and regression evidence.

## Human roles

- **Core maintainers:** currently Philip Cardiff and Ivan Batistić. They own
  architectural decisions, final integration, approval, and merging.
- **Contributing reviewers:** experienced contributors who volunteer to review.
  Their substantive reviews are valued and should be recognised in discussions.
  They assess methods, reproduce results, and explain concerns. Their reviews
  inform but do not replace core approval.
- **Automated AI reviewers:** provide preliminary technical assessments and
  identify decisions that need human judgement. They can be wrong or miss bugs.

Names here describe current roles; they are not permission configuration. Do not
infer GitHub usernames. If a reviewer directory is useful later, add opt-in,
verified profile links and areas of expertise to CONTRIBUTING.md, or use an
organisation team. No permissions or automatic reviewer requests are changed.
Invite relevant expertise rather than requesting every reviewer on every PR.

Progression is lightweight and based on demonstrated work:

1. **Contributor:** build experience through sustained fixes, features,
   examples, and reproducible numerical verification.
2. **Contributing reviewer:** volunteer for reviews in familiar areas. Pair with
   a core maintainer; compare assessments and discuss missed issues. Recognise
   review quality, clear reasoning, and constructive follow-through.
3. **Core maintainer:** after sustained contributions and strong reviews, the
   core maintainers discuss progression with the contributor. Consider broad
   architectural understanding, sound numerical and technical judgement, and
   responsible integration decisions. Promotion is a human decision, not an
   automatic result of review counts. No new administrative system is needed.

## One standard AI review

Use this workflow as the standard AI review mechanism. Human reviews should add
independent technical judgement, numerical assessment, verification, and
interpretation of AI findings. Private AI assistance is welcome, but do not
routinely paste another full AI review into the PR. Own and check any findings
you bring into a human review.

For another AI assessment, ask a maintainer to re-run the structured workflow,
or agree a specific unresolved concern to investigate. State its scope and
report only new evidence or unresolved findings. No new comment-command system
is introduced. Broader multi-model reviews belong in the coordinated extension
below, with one consolidated report rather than competing bot comments.

## Risk and future approval policy

Classify by consequences, not diff size. Use the highest relevant risk.

| Risk | Examples | Intended future policy (not enforced) |
| --- | --- | --- |
| LOW | Documentation; no expected behaviour change | Lightweight process |
| MEDIUM | Isolated functionality | One core-maintainer approval |
| HIGH | Core numerics, architecture, security | Both core-maintainer approvals |

The intended one/both policy refers to the two current core maintainers. Role
membership must be verified before any future enforcement is introduced. This
version changes no branch protection, approval rules, or permissions.

The current pilot uses one strong reviewer at every risk level. See the [model
policy](model-policy.md) for the planned risk-dependent reviewers, verified API
options, independence rules, and consolidation safeguards.

## Setup

1. Review this PR against `development`, the normal contribution target. Deploy
   the workflow and its helper/instruction files to the default branch
   (`master`) before enabling it. GitHub's current event documentation specifies
   that `pull_request_target` runs in the default-branch context. The explicit
   `github.sha` checkout uses that trusted workflow snapshot. See [GitHub's
   event
   reference](https://docs.github.com/en/actions/reference/workflows-and-actions/events-that-trigger-workflows#pull_request_target).
   Confirm the ref and SHA in the first live run. This draft PR cannot test its
   own privileged workflow before deployment.
2. Create an organisation-owned Anthropic API key for a dedicated workspace. Set
   a spending limit and monitor usage. Store it as the organisation Actions
   secret `ANTHROPIC_API_KEY`, restricted to this repository. A repository
   secret can be used for a pilot. Do not use a personal Claude subscription
   token.
3. An organisation owner must review the Actions event policy before enabling
   this privileged workflow. GitHub's default policy blocks
   `pull_request_target` in affected public repositories, with enforcement from
   2 December 2026. If the organisation accepts this design, explicitly allow
   this event at the narrowest available scope. See [GitHub's policy
   guidance](https://docs.github.com/en/actions/reference/security/securely-using-pull_request_target).
   This PR does not change the policy or repository settings.
4. Create the `ai-review-approved` label before using the first-time-author
   path. A maintainer or triager applies it after checking the PR's relevance
   and cost. It authorises reviews on subsequent commits to that PR, not merge
   approval.
5. Allow the pinned checkout, github-script, upload-artifact, and Anthropic
   Action in organisation Actions policy. The workflow supplies `GITHUB_TOKEN`
   explicitly; no GitHub App installation, personal access token, or custom App
   is required. The model job needs Contents read and Pull requests read. Only
   the separate publisher job needs Pull requests write to create a conversation
   comment. It never calls the approval or merge APIs. The publisher also needs
   Actions read to detect explicit timeouts. No Issues write or OIDC is needed.
6. After reviewing security and costs, set repository Actions variable
   `CLAUDE_REVIEW_ENABLED` to `true`. Until then AI execution is skipped;
   finalisation still uploads a skipped status record. Set it to `false` to stop
   new reviews. Cancel active runs separately if needed.
7. Run the live-validation scenarios below before relying on report delivery. Do
   not make this advisory workflow a required branch-protection check.

The official integration supports a supplied workflow token. This avoids the
broader permissions of the default Claude GitHub App. See the pinned [Action
inputs](https://github.com/anthropics/claude-code-action/blob/4fa93a88de49bbef68b10eb15aabedf28ca0e1ed/action.yml)
and [security
guidance](https://github.com/anthropics/claude-code-action/blob/4fa93a88de49bbef68b10eb15aabedf28ca0e1ed/docs/security.md).

## Operation and security

Eligible PRs target `development` or `master`, are open, non-draft, and have a
human author. Opening, reopening, marking ready, or adding commits triggers a
review for OWNER, MEMBER, COLLABORATOR, or CONTRIBUTOR author associations.
Other authors, including first-time fork contributors, need the triage-applied
`ai-review-approved` label. Adding that label triggers a review. Unrelated label
events do not cancel an active review or invoke the model. Known contributors
and labelled PRs can still trigger repeated runs; spending limits remain needed.
These associations are cost controls, not proof that PR content is trusted. No
comment commands or reviewer requests are installed. Editing only the PR body
does not trigger review; re-run the workflow if evidence was added there.

`pull_request_target` gives access to base-repository secrets for fork PRs. Both
jobs check out only the event's trusted base commit (`github.sha`), with no
persisted checkout credentials. Neither job fetches or checks out PR code. The
collector reads patches through the GitHub API. It never executes them, loads
their configuration, installs their dependencies, or runs their tests. Only
trusted base scripts run. Never change checkout to the PR head or merge ref,
even to inspect a proposed instruction file.

Claude receives the review instructions, the repository guide, and a bounded
JSON payload with the PR description and patches. The guide is the solids4foam
agent skill (`.agents/skills/solids4foam/SKILL.md`, the same file every coding
agent uses) plus the latest released version from `CHANGELOG.md`. Both are read
from the trusted base checkout on every run, never copied, so a PR that edits
the skill cannot change its own review. The versioned instructions explicitly
treat the PR payload as untrusted. Only the
schema-return tool `StructuredOutput` is enabled. Shell, filesystem, network,
and GitHub tools are unavailable. MCP is empty, strict, and denied. Bare mode,
user-only settings, and disabled commands prevent project skills and hooks from
loading. The runner is GitHub-hosted and ephemeral. The model has a read-only
GitHub token and no publishing token. The publisher uses a fresh runner, has no
Anthropic secret, validates structured output, and escapes model text before
creating a comment. Path floors cover all changed and renamed paths, including
omitted files: core numerical/model/solver directories and CI, tests, build
scripts, and agent configuration force HIGH. Other `src/` changes are at least
MEDIUM. Non-document changes default to at least MEDIUM. These floors are
conservative heuristics, not an exhaustive semantic classifier.

Prompt injection can still distort analysis. It cannot be solved by prompts
alone. Human reviewers must check findings and risk. This boundary prevents
external tool use and source writes; it does not guarantee a correct assessment.
Untrusted content is sent to Anthropic for analysis. Do not put credentials in
PRs. Keep debug/full-output logging disabled, including Actions step debug.
Review action updates before changing the full commit pins. The upstream Action
installs its Claude runtime; the action pin does not freeze every transitive
runtime dependency. Repeat live validation after updates.

Each run has a ten-minute model timeout, a small turn limit, and a USD 2 model
budget. Opus 5.5 and explicit high effort are selected in the workflow.
Organisation spending limits are still needed: a public contributor can trigger
repeated runs. Concurrency cancels older runs for the same PR. Reports include
the requested model/effort, exact head SHA, attempt, and run link. One sticky
bot comment is updated across commits, including on failure. Its SHA must match
the current PR head; an older SHA is historical evidence only. A final API check
suppresses stale publication, though a commit can still arrive between that
check and posting; the explicit SHA remains authoritative.

The collector includes whole patches within a 90 KB JSON-escaped prompt budget,
including instructions, the repository guide (about 18 KB with the
instructions) and metadata. This leaves margin for the Action's
serialized environment inputs. Source files precede application code, other
text, and tutorial mesh/field/CSV data. Missing patches (binary, empty files,
renames) and budget-excluded files are listed as not reviewed in the prompt and
report. Their paths still contribute to the risk floor.

A successful assessment can have COMPLETE or PARTIAL **patch** coverage. Neither
means the surrounding source, linked plots, CI logs, or numerical cases were
independently inspected. Partial coverage alone does not fail the job. No
reviewable patches, an incomplete file listing (3000 or more files), oversized
metadata, collection/API errors, timeouts, or invalid model output produce an
INCOMPLETE comment with a safe reason code and a failing job. Failure signalling
is retained so it cannot be confused with success. Risk then defaults to HIGH. A
cancelled run attempts to publish INCOMPLETE if its SHA is still current, but
infrastructure/forced cancellation can prevent any comment. Publication failure
also fails the job. **No current-SHA report means no completed AI review.** A
green job is not an approval. Humans must inspect the complete diff and context.

## Attempt status

Each run attempt uploads a small machine-readable status artifact, including
skipped and failed reviews. Completion requires successful report publication.
See [status records](status-records.md) for fields, conservative error reasons,
optional numeric usage, retention, and retrieval limits. No fallback is enabled.

## Verification

Run the offline tests from the repository root:

```sh
node --test .github/scripts/claudeReview*.test.cjs
markdownlint .github/review/*.md .github/pull_request_template.md
CONTRIBUTING.md
```

Tests use mocked GitHub APIs and model responses. They verify data handling,
comment publication, and failure semantics. They do not establish model quality
or live integration. No OpenFOAM build is needed for these workflow-only
changes.

Use disposable PRs in a test repository with the same settings for live checks:

<!-- markdownlint-disable MD013 -->

| Scenario | Expected result |
| --- | --- |
| Documentation only | LOW; clear summary; no numerical test demanded |
| Isolated new boundary condition | At least MEDIUM; relevant compatibility questions |
| Core solver/default or constitutive change | HIGH, even for a small diff |
| Numerical change without results | Missing verification stated; practical human checklist |
| Add a second commit | Same sticky report; current SHA; stale results suppressed |
| First-time external fork | Requires approval label; no head checkout or PR code execution |
| Hidden instructions in title/body/patch | Only StructuredOutput; no disclosure, writes, or approval |
| Invalid API key or model timeout | INCOMPLETE report; failed run; no completed report |
| Workflow permission change | HIGH; proposed permissions never take effect in this run |

<!-- markdownlint-enable MD013 -->

Also test draft/ready transitions, label gating, a missing secret, binary files,
malformed output, and denied comment permission. Confirm the comment author,
SHA, model/effort, run link, escaping, and cross-SHA sticky comment updates on
GitHub. Confirm no review approval, source commit, or merge occurs. Review
Actions logs for accidental sensitive output. Model risk and checklist quality
require human assessment.

## Independent review before wider team discussion

Ask an independent reviewer to inspect the trusted-ref checkout, token scopes,
tool allowlist, injection handling, and whether the pinned Action correctly
passes model/effort and structured-output options. Verify failure reports, HIGH
fallback, comment update ownership checks, and fork/cancellation behaviour in a
test repo. Evaluate sample numerical PRs for risk and finding quality. Confirm
Opus API access, costs, and data handling with the organisation account owner.
Review the human progression language and the planned two-reviewer coverage with
the core maintainers before enabling the pilot or asking the wider team to rely
on it.
