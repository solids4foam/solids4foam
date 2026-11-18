# Review attempt status artifacts

Each workflow run attempt writes `review-status.json` under `RUNNER_TEMP` and
uploads only that file as `review-status-<run_id>-<run_attempt>`. It is not a PR
comment or a committed file. No fallback reviewer, scheduling, or retry across
accounts is implemented.

| Field | Meaning |
| --- | --- |
| `schema_version` | Integer format version; currently `2` |
| `repository`, `pr` | Repository name and numeric PR number |
| `head_sha` | Exact PR head from the triggering event, never the workflow SHA |
| `review_type` | `primary` in this pilot |
| `reviewer` | `claude-api` in this pilot |
| `status` | `completed`, `failed`, or `skipped` |
| `failure_reason` | Null only for completed; otherwise a machine-readable code |
| `run_id`, `run_attempt` | GitHub run ID string and numeric attempt number |
| `coverage` | Patch coverage: `complete`, `partial`, `incomplete`, or null |
| `api_error_status` | HTTP status of a failed API request (100-599), or null |
| `served_model` | Model that served the report (a fallback, if used), or null |
| `model_requests` | Number of Messages API requests made (at most 2), or null |
| `usage` | Summed token counts from the API responses, or null |

`completed` requires successful model execution, a valid report, and a successful
GitHub comment create/update. A successful model invocation or INCOMPLETE comment
alone cannot produce it. `failed` includes model, validation, and publication
failures. `skipped` means disabled/ineligible (`not_eligible_or_disabled`) or a PR
that changed or became ineligible before publication (`pr_changed_or_ineligible`).
A completed record may have partial patch coverage. It is historical evidence
for its stated SHA, never an approval.

Failure reasons use evidence, not guesses. The analysis step maps typed SDK
errors and response stop reasons to fixed codes; it never copies error text:

- `authentication_error` (HTTP 401), `permission_error` (403), `rate_limited`
  (429), and `api_error` (any other HTTP error; see `api_error_status`). A
  billing problem such as no remaining credit is `api_error` with status 400.
- `connection_error` and `timeout`: the SDK could not reach the API, or the
  8.5-minute analysis deadline passed. A streamed error has no HTTP status; its
  structured error type is mapped instead and `api_error_status` stays null.
- `refusal`: the model, and any fallback model, declined the request.
- `max_tokens`: the report was cut off at the output limit.
- `invalid_report`: the report failed validation twice.
- `no_api_key`: the secret was not available to the analysis step.
- `no_input`: the collector's prompt file was missing or malformed.
- `timeout` is also used when GitHub explicitly reports `timed_out` for the
  review job or a step. An ambiguous cancellation does not establish timeout.
- `review_failed`: a known failed review or publication without a more specific
  reliable signal.
- `unknown`: no reliable cause is available, including ambiguous cancellation.

Only those codes, HTTP status numbers, validated model identifiers and
non-negative integer token counts are copied. Message text, error descriptions,
prompts, and credentials are never copied into status outputs or artifacts.
Usage is the sum over every sampling iteration of the run's API requests
(including a declined fallback attempt), not a billing ledger, and can be
absent on early failure. Re-running only failed jobs repeats the analysis; do
not sum records across attempts as if they were one review.

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
npm ci --prefix .github/scripts --ignore-scripts
node --test .github/scripts/claudeReview*.test.cjs
```
