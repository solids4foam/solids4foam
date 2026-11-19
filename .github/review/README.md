# Advisory PR reviews

The Claude review helps people assess a PR. It does not approve a PR, count
approvals, enforce a merge gate, or merge changes. Existing CI still supplies
build and regression evidence. The human review and merge process, by risk
level, is in [CONTRIBUTING.md](../../CONTRIBUTING.md#reviewing-pull-requests);
this guide covers the workflow itself.

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
is introduced. The only planned additional AI review is the second
independent review for MEDIUM/HIGH risk, described below.

## Risk and approval policy

Classify by consequences, not diff size. Use the highest relevant risk.

<!-- markdownlint-disable MD013 -->

| Risk | Examples | Agreed policy (not enforced by GitHub) |
| --- | --- | --- |
| LOW | Documentation; no expected behaviour change | One core maintainer |
| MEDIUM | Isolated functionality | Second AI review; one core maintainer |
| HIGH | Core numerics, architecture, security | Second AI review; both core maintainers |

<!-- markdownlint-enable MD013 -->

The one/both policy refers to the two current core maintainers. It is agreed
practice, not branch protection: no approval rules or permissions change.
Each report states the rule for its risk and links the process.

The publisher labels each PR with its risk: `ai-risk: low`, `ai-risk: medium`
or `ai-risk: high`, plus `ai-review: incomplete` when no successful review
exists (which is also labelled HIGH). It replaces only these four labels and
updates them on every published report; a suppressed (stale) run leaves them
unchanged. A label failure is a warning, not a failed run: the comment is
authoritative. Labels added with the workflow token do not trigger new runs.

The workflow runs one strong reviewer at every risk level. See the [model
policy](model-policy.md) for the planned risk-dependent reviewers, verified API
options, independence rules, and consolidation safeguards.

## Second independent AI review (interim, manual)

MEDIUM and HIGH PRs need a second AI review, preferably from a different
model family. Until that is automated, the Claude report for such a PR (and
any failed run) asks for it, gives both commands below, and adds a checklist
item. A core maintainer then runs, from an up-to-date `development` checkout:

```bash
node .github/scripts/secondReview.cjs <PR number> --post                # Codex
node .github/scripts/secondReview.cjs <PR number> --model fable --post  # Fable
```

The script builds exactly the prompt CI would build: the instructions, skill
and `CHANGELOG.md` are fetched from the PR's base branch on GitHub, not from
the local tree, and the PR is untrusted data. Codex (`gpt-6-astra`, `xhigh`
effort) reviews it without seeing the Claude report or any PR comment, so
agreement between the two is independent. The result is validated and escaped
like the CI report and posted as one separate comment per maintainer; a rerun
updates it. Without `--post`, it only writes the report locally.

Because the PR text is untrusted, Codex runs without tools: shell, code
execution, apps, plugins, browser, image, agent and web tools are disabled,
inside a read-only sandbox in an empty temporary directory, with no user
config, MCP servers or rules. Recheck the disabled feature list in the script
when Codex is upgraded. It needs `gh` (authenticated) and a Codex login with
access to `gpt-6-astra`; it uses that maintainer's Codex allowance. Codex
may warn that "Code Mode is unavailable because code-mode host is disabled"
and "will fail closed". This is expected: the script disables that feature on
purpose, and failing closed means no code can run. Do not enable it.

`--model fable` is the substitute when Codex is unavailable: Claude Fable 5.2
at high effort, through the same `analyse` function as CI (no tools,
structured output, validation, one retry, server-side fallback). It needs
`ANTHROPIC_API_KEY` set to a key with Fable access, billed to that key's
workspace, and the pinned SDK (`npm ci --prefix .github/scripts
--ignore-scripts`). Its comment says it shares the Claude report's model
family, so agreement is weaker evidence. It suits MEDIUM risk; for HIGH, use it
only when a core maintainer judges the classification borderline. Fable 5.2
requires 30-day data retention on the API account; confirm that is acceptable
before first use.

## Setup

1. Merge this PR into `development` as usual. GitHub always runs
   `pull_request_target` workflows from the default branch (`master`), even for
   PRs into `development`; see [GitHub's event
   reference](https://docs.github.com/en/actions/reference/workflows-and-actions/events-that-trigger-workflows#pull_request_target).
   So also copy **only** `.github/workflows/claudeReview.yml` to `master`, for
   example in a small maintainer PR. The workflow then checks out the tip of
   `development` for the scripts, instructions and skill, so PRs are reviewed
   with `development`'s current rules and `master` needs no other review file.
   Only PRs into `development` are reviewed: PRs into `master` are release
   merges of reviewed work, and `master` does not carry the review scripts.
   Later changes to the YAML itself (triggers, permissions, model flags) only
   take effect once copied to `master`. Confirm the checked-out commit in the
   first live run. This PR cannot test its own privileged workflow before
   deployment.
2. Create an organisation-owned Anthropic API key for a dedicated workspace. Set
   a spending limit and monitor usage. Store it as the organisation Actions
   secret `ANTHROPIC_API_KEY`, restricted to this repository. A repository
   secret can be used for a pilot. Prefer a service-account key over one linked
   to a person. Do not use a personal Claude subscription token. The workflow
   calls the Messages API directly, so it is billed as API usage; organisation
   API credit that covers "API" applies to it.
3. An organisation owner must review the Actions event policy before enabling
   this privileged workflow. GitHub's default policy blocks
   `pull_request_target` in affected public repositories, with enforcement from
   2 December 2026. If the organisation accepts this design, explicitly allow
   this event at the narrowest available scope. See [GitHub's policy
   guidance](https://docs.github.com/en/actions/reference/security/securely-using-pull_request_target).
   This PR does not change the policy or repository settings.
4. Create the `ai-risk: low`, `ai-risk: medium`, `ai-risk: high` and
   `ai-review: incomplete` labels (the publisher would otherwise create them
   uncoloured). Create the `ai-review-approved` label before using the
   first-time-author path. A maintainer or triager applies it after checking
   the PR's relevance and cost. It authorises reviews on subsequent commits to
   that PR, not merge approval.
5. Allow the pinned checkout, github-script, and upload-artifact actions in
   organisation Actions policy. The workflow supplies `GITHUB_TOKEN`
   explicitly; no GitHub App installation, personal access token, or custom App
   is required. The analysis job needs Contents read and Pull requests read. Only
   the separate publisher job needs Pull requests write to create a conversation
   comment and set the risk labels. It never calls the approval or merge APIs.
   The publisher also needs Actions read to detect explicit timeouts. No Issues
   write or OIDC is needed.
6. After reviewing security and costs, set repository Actions variable
   `CLAUDE_REVIEW_ENABLED` to `true`. Until then AI execution is skipped;
   finalisation still uploads a skipped status record. Set it to `false` to stop
   new reviews. Cancel active runs separately if needed.
7. Run the live-validation scenarios below before relying on report delivery. Do
   not make this advisory workflow a required branch-protection check.

The analysis calls the Claude Messages API directly through the official
Anthropic TypeScript SDK, pinned to an exact version in
`.github/scripts/package-lock.json` and installed with `npm ci --ignore-scripts`.
No Claude GitHub App, agent action, personal access token, or custom App is
used. The API key is a service-account key restricted to one workspace; only the
analysis step receives it.

## Operation and security

Eligible PRs target `development`, are open, non-draft, and have a
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
jobs check out only the tip of the PR's base branch, which the branch filter
limits to `development`. It is maintainer-controlled, so code merged there is
trusted as much as code on `master`. Checkout credentials are not
persisted. Neither job fetches or checks out PR code. The
collector reads patches through the GitHub API. It never executes them, loads
their configuration, installs their dependencies, or runs their tests. Only
trusted base scripts run. Never change checkout to the PR head or merge ref,
even to inspect a proposed instruction file.

The analysis sends one Messages API request. The system prompt holds the review
instructions and the repository guide: the solids4foam agent skill
(`.agents/skills/solids4foam/SKILL.md`, the same file every coding agent uses)
plus the latest released version from `CHANGELOG.md`. Both are read from the
trusted base checkout on every run, never copied, so a PR that edits the skill
cannot change its own review. The user message holds only a bounded JSON payload
with the PR description and patches, which the instructions treat as untrusted.
The request defines no tools, so the model has no shell, filesystem, network, or
GitHub access; it can only return JSON in the report schema (structured output).
The API cannot enforce length and count limits, so the schema states them and
the workflow validates them; a report that fails validation gets one retry that
names the failed check. The analysis step outputs only a validated report, token
counts and fixed reason codes, never prompt, model, or API error text, and SDK
logging is off. The runner is GitHub-hosted and ephemeral. The analysis job has
a read-only GitHub token and no publishing token. The publisher uses a fresh
runner, has no Anthropic secret, validates the report again, and escapes model
text before creating a comment. Path floors cover all changed and renamed paths,
including omitted files: core numerical/model/solver directories and CI, tests,
build scripts, and agent configuration force HIGH. Other `src/` changes are at
least MEDIUM. Non-document changes default to at least MEDIUM. These floors are
conservative heuristics, not an exhaustive semantic classifier.

If the model declines a request for safety reasons, the API's server-side
fallback re-runs it on a fallback model in the same call. The report then says
which model served it.

Prompt injection can still distort analysis. It cannot be solved by prompts
alone. Human reviewers must check findings and risk. This boundary prevents
external tool use and source writes; it does not guarantee a correct assessment.
Untrusted content is sent to Anthropic for analysis. Do not put credentials in
PRs. Review SDK updates (`package.json` and the lockfile) and action pin updates
before merging them, and repeat live validation after updates.

Each run has an 8.5-minute analysis deadline (inside a ten-minute step timeout)
and at most two Messages API requests: the second only to fix an invalid report
or after one transient API error (408, 409, 429, 5xx or overloaded), after
waiting as the server's Retry-After asks. SDK retries are off, so every request
is counted. Each request allows up to 40,000 output tokens (USD 0.80 at most
at Opus 5.5 prices of USD 4/20 per million input/output tokens) and a full
prompt budget of about 250,000 input tokens (about USD 1). A small PR costs a
few cents; a 335 KB patch such as #719 costs about USD 1. A safety fallback
can re-run a request on the fallback model inside the same request, which adds
its tokens. Set a monthly spending limit on the API workspace below the
available credit.
Opus 5.5 and explicit high effort are selected in the workflow. Organisation
spending limits are still needed: a public contributor can trigger repeated
runs. Concurrency cancels older runs for the same PR. Reports include the
requested model/effort, exact head SHA, attempt, and run link. One sticky bot
comment is updated across commits, including on failure. Its SHA must match the
current PR head; an older SHA is historical evidence only. A final API check
suppresses stale publication, though a commit can still arrive between that
check and posting; the explicit SHA remains authoritative.

The collector includes whole patches within a 600 KB JSON-escaped prompt budget,
including instructions, the repository guide (about 18 KB with the instructions)
and metadata. The collector writes the prompt to a file in the runner's
temporary directory for the analysis step; it is never passed through
environment variables, which Actions prints in the step log. Source files
precede application code, other text, and tutorial mesh/field/CSV data. Missing
patches (binary, empty files, renames) and budget-excluded files are listed as
not reviewed in the prompt and report. Their paths still contribute to the risk
floor.

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
npm ci --prefix .github/scripts --ignore-scripts
node --test .github/scripts/*.test.cjs
markdownlint .github/review/*.md .github/pull_request_template.md CONTRIBUTING.md
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
| Hidden instructions in title/body/patch | No tools; JSON report only; no disclosure, writes, or approval |
| Invalid API key, API error or timeout | INCOMPLETE report naming the reason code and HTTP status; failed run |
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
absence of tools, injection handling, and whether the analysis step passes
model, effort and the report schema correctly. Verify failure reports, HIGH
fallback, comment update ownership checks, and fork/cancellation behaviour in a
test repo. Evaluate sample numerical PRs for risk and finding quality. Confirm
Opus API access, costs, and data handling with the organisation account owner.
Review the human progression language and the planned two-reviewer coverage with
the core maintainers before enabling the pilot or asking the wider team to rely
on it.
