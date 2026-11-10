// Read only structured SDK envelopes. Never copy messages, errors, or credentials.
const fs = require('node:fs');

function telemetry(messages) {
  if (!Array.isArray(messages)) return {failure_reason: 'unknown', usage: null};
  const result = messages.filter(item => item?.type === 'result').at(-1);
  const error = messages.filter(item => item?.type === 'assistant').at(-1)?.error;
  const reasons = {rate_limit: 'rate_limited', authentication_failed: 'authentication_error'};
  // billing_error covers more than quota exhaustion; do not guess from its text.
  const failure_reason = Object.hasOwn(reasons, error) ? reasons[error] : 'unknown';
  const usage = {};
  for (const key of ['input_tokens', 'output_tokens', 'cache_read_input_tokens',
    'cache_creation_input_tokens']) {
    const value = result?.usage?.[key];
    if (typeof value === 'number' && Number.isFinite(value) && value >= 0) usage[key] = value;
  }
  const cost = result?.total_cost_usd;
  if (typeof cost === 'number' && Number.isFinite(cost) && cost >= 0) usage.cost_usd = cost;
  return {failure_reason, usage: Object.keys(usage).length ? usage : null};
}

function capture(core, executionFile) {
  let summary = telemetry(null);
  try {
    summary = telemetry(JSON.parse(fs.readFileSync(executionFile, 'utf8')));
  } catch {
    // Missing or malformed telemetry must not prevent workflow finalisation.
  }
  core.setOutput('telemetry', JSON.stringify(summary));
}

function statusRecord(context, env, timedOut = false) {
  let summary = {};
  try { summary = JSON.parse(env.REVIEW_TELEMETRY); } catch { /* No reliable signal. */ }
  let publication = {};
  try { publication = JSON.parse(env.PUBLICATION_STATUS); } catch { /* Not published. */ }
  let status = 'failed';
  let reason = 'unknown';
  if (env.REVIEW_RESULT === 'skipped') {
    status = 'skipped';
    reason = 'not_eligible_or_disabled';
  } else if (publication?.status === 'skipped') {
    status = 'skipped';
    reason = 'pr_changed_or_ineligible';
  } else if (env.REVIEW_RESULT === 'success' && env.CLAUDE_CONCLUSION === 'success' &&
      env.PUBLICATION_OUTCOME === 'success' && publication?.status === 'completed' &&
      ['complete', 'partial'].includes(publication.coverage)) {
    status = 'completed';
    reason = null;
  } else if (timedOut) {
    reason = 'timeout';
  } else if (env.REVIEW_RESULT !== 'success' &&
      ['rate_limited', 'authentication_error'].includes(summary?.failure_reason)) {
    reason = summary.failure_reason;
  } else if (env.REVIEW_RESULT === 'failure' || env.PUBLICATION_OUTCOME === 'failure' ||
      publication?.status === 'failed') {
    reason = 'review_failed';
  }
  return {
    schema_version: 1,
    repository: `${context.repo.owner}/${context.repo.repo}`,
    pr: context.payload.pull_request.number,
    head_sha: context.payload.pull_request.head.sha,
    review_type: 'primary', reviewer: 'claude-api', status, failure_reason: reason,
    run_id: String(context.runId), run_attempt: Number(env.GITHUB_RUN_ATTEMPT),
    coverage: status === 'completed' && ['complete', 'partial'].includes(publication?.coverage)
      ? publication.coverage : status === 'skipped' ? null : 'incomplete',
    // Re-sanitise across the job boundary, rather than copying arbitrary JSON.
    usage: telemetry([{type: 'result', usage: summary?.usage,
      total_cost_usd: summary?.usage?.cost_usd}]).usage
  };
}

async function finalise({github, context, core}, env = process.env) {
  let timedOut = false;
  try {
    const jobs = await github.paginate(github.rest.actions.listJobsForWorkflowRunAttempt,
      {...context.repo, run_id: context.runId, attempt_number: Number(env.GITHUB_RUN_ATTEMPT),
        per_page: 100});
    timedOut = jobs.some(job => job.name === 'review' &&
      (job.conclusion === 'timed_out' || job.steps?.some(step => step.conclusion === 'timed_out')));
  } catch {
    core.notice('Job timing metadata unavailable; no timeout reason inferred.');
  }
  const record = statusRecord(context, env, timedOut);
  fs.writeFileSync(`${env.RUNNER_TEMP}/review-status.json`, JSON.stringify(record, null, 2) + '\n');
}

module.exports = {telemetry, capture, statusRecord, finalise};
