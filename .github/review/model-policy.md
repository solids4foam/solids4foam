# AI review policy and next extension

## Current pilot

One Anthropic reviewer runs for every eligible PR, irrespective of risk. The
workflow selects `claude-opus-5-5` with explicit `high` effort. Reports show the
requested API model and effort. This meets the intended single-reviewer model
selection; numerical review quality still needs evaluation on real PRs.

The MEDIUM/HIGH two-reviewer target below is **not automated**. As an interim
step, MEDIUM/HIGH and failed reports ask a core maintainer to run
`.github/scripts/secondReview.cjs`. It sends the same prompt to Codex
(`gpt-6-astra`, `xhigh`) without the first report and posts a separate,
validated comment. Two separate comments are not yet one consolidated report;
the [README](README.md#second-independent-ai-review-interim-manual) gives
details. These manual pairs will show how often the reviewers agree and differ,
which informs the automated extension below.

Classification, execution, presentation, and merge decisions are separate:

1. The collector obtains one bounded PR snapshot and deterministic numerical
   and security path floors. No PR code is executed.
2. The official Claude Action performs one assessment and proposes a risk.
   There is no separate classifier API call and no risk-based execution branch.
3. `classifyRisk` applies the conservative floor. Missing or invalid risk/floor
   data means HIGH. An invalid report or failed execution remains INCOMPLETE;
   HIGH fallback is a routing precaution, not a successful classification.
4. Findings have kind, severity, category, location, claim, evidence, and a
   suggested check. `formatAssessment` presents every finding without AI
   synthesis. The publisher validates and escapes data and updates one sticky
   report per PR, including on failure. Each update identifies its SHA.
5. Human reviewers assess the evidence. Core maintainers decide integration
   and merging. AI risk or completion never supplies merge authority.

Keeping a fixed strong reviewer avoids a failed classifier selecting a weaker
review path. Before enabling routing, move classification ahead of independent
review execution and use HIGH for every failure or unknown result. Reviewers
and humans may raise risk; a lower model suggestion must not override the floor.

## Target policy (planned)

| Risk | Independent reviewers |
| --- | --- |
| LOW | One strong reviewer: Opus 5.5 or Sol 6.1 class or better |
| MEDIUM | Two reviewers from different families/classes |
| HIGH | Two; at least one Fable 5.2 class or better at high effort |

Select explicit API IDs and effort after checking organisation access and
project-specific quality. Do not infer equivalence from a model's name alone.
An Anthropic/OpenAI pair offers useful family diversity. Two calls to the same
model, or one review that reads the other first, do not meet this policy.
For HIGH, both reviewers must examine numerical correctness, verification,
architecture, and regression risks. Cost caps must not trigger silent downgrades.

The next extension should:

1. Send the same immutable SHA, evidence, and review instructions to two isolated
   jobs. Neither receives the other's conclusions or PR review comments.
2. Validate each result and retain the original structured findings, model and
   effort metadata, input SHAs, and completion state in run artifacts. Link them
   from the report. Exclude credentials and internal reasoning traces.
3. Consolidate duplicate findings into one concise report. Preserve provenance,
   unique concerns, disagreements, and unresolved significant findings. A
   synthesis model must not silently dismiss such findings. Any proposed
   dismissal must remain visible with its reason for human resolution.
4. Mark partial coverage explicitly if either reviewer fails. Never present one
   surviving review as satisfying the two-reviewer target. Test preservation of
   significant findings and disagreements before enabling synthesis.

A natural first automated step is a deterministic publisher that lists both
reviewers' findings in one comment, grouped by file and line and tagged by
source, with no synthesis model. Duplicates then sit side by side and nothing
can be dropped.

This PR defers the automated second API integration, ahead-of-review classification,
original-result artifacts, and multi-model consolidation. They require another
credential/billing setup, schema and retention decisions, and failure tests.
These are planned extensions, not hidden optional paths in the pilot.

## Verified API options

Checked against official documentation on 10 November 2026. Availability in
public docs does not confirm access for the solids4foam organisation. Before
activation, test the exact model and effort with its organisation API credential.
Do not silently substitute an unavailable model.

- Anthropic documents `claude-opus-5-5` on the Claude API. This pilot passes
  `--effort high` through the official Action to Claude Code. The Messages API
  equivalent is `output_config.effort: "high"`. See the
  [model page](https://platform.claude.com/docs/en/models/opus-5-5/overview),
  [effort reference](https://platform.claude.com/docs/en/build-with-claude/effort),
  and [CLI flags](https://code.claude.com/docs/en/cli-reference).
- OpenAI documents `gpt-6.1-sol` with `reasoning.effort` values including `high`
  through the Responses API. It needs an OpenAI integration; the Anthropic
  Action does not run OpenAI models. See the
  [model page](https://developers.openai.com/api/docs/models/gpt-6.1-sol).
- Anthropic documents `claude-fable-5-2`, with high effort and adaptive thinking.
  It is a future HIGH-risk candidate, not enabled here. Its documented 30-day
  retention requirement needs review before adoption. See the
  [model page](https://platform.claude.com/docs/en/models/fable-5-2/overview).

The pilot only needs the existing organisation-owned `ANTHROPIC_API_KEY` and
Anthropic API billing with access to Opus 5.5. A future OpenAI reviewer needs a
separately managed organisation/project API key (for example `OPENAI_API_KEY`),
billing, spending limits, and model access. Fable may use the Anthropic account
but its access and terms need confirmation. Interactive Claude/Codex access is
not proof of API entitlement. Some plans offer API credits through explicit
account linking; verify eligibility and limits rather than assuming coverage.
See [Anthropic credits](https://platform.claude.com/docs/en/about-claude/api-credits-for-subscribers)
and [OpenAI credits](https://developers.openai.com/api/docs/guides/chatgpt-api-credits).
