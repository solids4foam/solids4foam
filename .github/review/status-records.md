# Review attempt status artifacts

Each workflow run attempt writes `review-status.json` under `RUNNER_TEMP` and
uploads only that file as `review-status-<run_id>-<run_attempt>`. It is not a PR
comment or a committed file. No fallback reviewer, scheduling, or retry across
accounts is implemented.

| Field | Meaning |
| --- | --- |
| `schema_version` | Integer format version; currently `1` |
| `repository`, `pr` | Repository name and numeric PR number |
| `head_sha` | Exact PR head from the triggering event, never the workflow SHA |
| `review_type` | `primary` in this pilot |
| `reviewer` | `claude-api` in this pilot |
| `status` | `completed`, `failed`, or `skipped` |
| `failure_reason` | Null only for completed; otherwise a machine-readable code |
| `run_id`, `run_attempt` | GitHub run ID string and numeric attempt number |
| `coverage` | Patch coverage: `complete`, `partial`, `incomplete`, or null |
| `usage` | Available numeric SDK token/cost totals, or null if unavailable |

`completed` requires successful model execution, a valid report, and a successful
GitHub comment create/update. A successful model invocation or INCOMPLETE comment
alone cannot produce it. `failed` includes model, validation, and publication
failures. `skipped` means disabled/ineligible (`not_eligible_or_disabled`) or a PR
that changed or became ineligible before publication (`pr_changed_or_ineligible`).
A completed record may have partial patch coverage. It is historical evidence
for its stated SHA, never an approval.

Failure reasons use evidence, not guesses:

- `rate_limited`: the SDK's last assistant envelope has `error: rate_limit`.
- `authentication_error`: that envelope has `error: authentication_failed`.
- `timeout`: GitHub explicitly reports `timed_out` for the review job or a step.
  An ambiguous cancellation or generic step failure does not establish timeout.
- `review_failed`: a known failed review/publication without a more specific
  reliable signal. This also covers invalid reports and generic billing errors.
- `unknown`: no reliable cause is available, including ambiguous cancellation.
- `api_quota_exceeded`: reserved for an explicit quota-exhaustion signal. The
  pinned Action does not guarantee one, so this version does not emit this code.
  Neither a 429/rate-limit code nor a generic billing error proves spent credits.

Only recognised structured SDK codes and non-negative finite token/cost numbers
are extracted from the local execution file. Message text, error descriptions,
authentication responses, prompts, and credentials are never copied into status
outputs or artifacts. Extraction failure yields unknown telemetry. Full/debug
logging must remain disabled. See the official
[SDK error definitions](https://code.claude.com/docs/en/agent-sdk/typescript#sdkassistantmessage).
Usage is best-effort SDK accounting, not a billing ledger or API request count.
It can be absent on early failure. Re-running only failed jobs can reuse earlier
model outputs; do not sum repeated usage as if every record were a new API call.

## Finalisation and retention limits

The existing publisher job uses `always()`, as do status generation and upload.
It runs even when review is skipped or fails. Publication has its own timeout so
there is time to finalise. A small `actions: read` permission on the publisher
allows one lookup of this run attempt's job conclusions; lookup failure does
not prevent writing the record. No job logs are downloaded or parsed.

Upload requests 30-day retention, subject to repository/organisation policy.
Artifacts can expire or be deleted. Missing files or upload errors fail the
upload step. No finalisation can guarantee an artifact after runner loss,
checkout failure, forced cancellation, GitHub outage, or a whole-job timeout.
**A missing artifact is unknown, never completed.** A record may say completed
while its upload fails: publication and artifact availability are separate facts.
Live tests must confirm upload after API/publication failures and skipped runs.

A future agent can list runs for the trusted `claudeReview.yml` workflow, list
that run's artifacts, and download the newest available attempt artifact. Check
its repository, PR, schema version, and `head_sha` against the current PR. Do not
filter workflow runs by PR head SHA: `pull_request_target` runs use a trusted
base-repository SHA. Ignore artifacts from other workflows and expired artifacts.
Use the highest run/attempt for that PR; do not fall back silently to an older
completed record if the newest attempt is missing or incomplete. Check current
PR state again before any action. Artifact access may require an authenticated
GitHub token with Actions read. See GitHub's
[artifact API](https://docs.github.com/en/rest/actions/artifacts).

Offline tests cover completed/failed/skipped records, publication failure,
conservative error mapping, timeout evidence, missing telemetry, safe usage
extraction, and writing the file when metadata lookup fails:

```sh
node --test .github/scripts/claudeReview*.test.cjs
```
