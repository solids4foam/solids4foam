// Only fixed reason codes and numbers cross into the record; never messages,
// errors, model text or credentials.
const fs = require('node:fs');
const {analysisReasons, httpStatus} = require('./claudeReview.cjs');

const usageKeys = ['input_tokens', 'output_tokens', 'cache_read_input_tokens',
  'cache_creation_input_tokens'];
function sanitiseUsage(usage) {
  const out = {};
  for (const key of usageKeys) {
    const value = usage?.[key];
    if (Number.isSafeInteger(value) && value >= 0) out[key] = value;
  }
  return Object.keys(out).length ? out : null;
}
// Parse the analysis step's telemetry output (see analyse() in claudeReview.cjs).
function telemetry(raw) {
  let summary = {};
  try { summary = JSON.parse(raw) || {}; } catch { /* No reliable signal. */ }
  return {
    failure_reason: analysisReasons.includes(summary.failure_reason) ? summary.failure_reason : null,
    api_error_status: httpStatus(summary.api_error_status),
    served_model: /^[a-z0-9.-]{1,64}$/.test(summary.served_model || '') ? summary.served_model : null,
    model_requests: Number.isSafeInteger(summary.attempts) && summary.attempts >= 0 ? summary.attempts : null,
    usage: sanitiseUsage(summary.usage)
  };
}

function statusRecord(context, env, timedOut = false) {
  const summary = telemetry(env.REVIEW_TELEMETRY);
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
  } else if (env.REVIEW_RESULT !== 'success' && summary.failure_reason) {
    reason = summary.failure_reason;
  } else if (env.REVIEW_RESULT === 'failure' || env.PUBLICATION_OUTCOME === 'failure' ||
      publication?.status === 'failed') {
    reason = 'review_failed';
  }
  return {
    schema_version: 2,
    repository: `${context.repo.owner}/${context.repo.repo}`,
    pr: context.payload.pull_request.number,
    head_sha: context.payload.pull_request.head.sha,
    review_type: 'primary', reviewer: 'claude-api', status, failure_reason: reason,
    run_id: String(context.runId), run_attempt: Number(env.GITHUB_RUN_ATTEMPT),
    coverage: status === 'completed' && ['complete', 'partial'].includes(publication?.coverage)
      ? publication.coverage : status === 'skipped' ? null : 'incomplete',
    // Re-sanitised across the job boundary, rather than copying arbitrary JSON.
    api_error_status: summary.api_error_status, served_model: summary.served_model,
    model_requests: summary.model_requests, usage: summary.usage
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

module.exports = {telemetry, statusRecord, finalise};
