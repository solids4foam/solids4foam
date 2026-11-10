const {test} = require('node:test');
const assert = require('node:assert/strict');
const fs = require('node:fs');
const os = require('node:os');
const path = require('node:path');
const {telemetry, capture, statusRecord, finalise} = require('./claudeReviewStatus.cjs');
const {publish} = require('./claudeReview.cjs');

const context = {repo: {owner: 'example', repo: 'test'}, runId: 123,
  payload: {pull_request: {number: 42, head: {sha: 'a'.repeat(40)}}}};
const env = extra => ({GITHUB_RUN_ATTEMPT: '2', REVIEW_RESULT: 'success',
  CLAUDE_CONCLUSION: 'success', PUBLICATION_OUTCOME: 'success',
  PUBLICATION_STATUS: '{"status":"completed"}', ...extra});

function check(record) {
  assert.equal(record.schema_version, 1);
  assert.equal(record.pr, 42);
  assert.equal(record.head_sha, 'a'.repeat(40));
  assert.equal(record.reviewer, 'claude-api');
  assert.equal(record.review_type, 'primary');
  assert.ok(['completed', 'failed', 'skipped'].includes(record.status));
  assert.equal(record.failure_reason === null, record.status === 'completed');
  assert.deepEqual(JSON.parse(JSON.stringify(record)), record);
}

test('completed requires both model success and confirmed publication', () => {
  const completed = statusRecord(context, env());
  check(completed);
  assert.equal(completed.status, 'completed');
  assert.equal(completed.failure_reason, null);
  for (const extra of [{PUBLICATION_OUTCOME: 'failure'}, {PUBLICATION_STATUS: ''},
    {PUBLICATION_STATUS: '{"status":"failed"}'}, {CLAUDE_CONCLUSION: 'failure'},
    {REVIEW_RESULT: 'failure'}, {PUBLICATION_STATUS: 'null'}]) {
    const record = statusRecord(context, env(extra));
    check(record);
    assert.equal(record.status, 'failed');
  }
});

test('disabled and stale PRs are skipped, not completed', () => {
  for (const extra of [{REVIEW_RESULT: 'skipped'}, {PUBLICATION_STATUS: '{"status":"skipped"}'}]) {
    const record = statusRecord(context, env(extra));
    check(record);
    assert.equal(record.status, 'skipped');
  }
});

test('only structured rate-limit and authentication codes get specific reasons', () => {
  for (const [error, reason] of [['rate_limit', 'rate_limited'],
    ['authentication_failed', 'authentication_error'], ['billing_error', 'unknown'],
    ['credit_balance_exhausted', 'unknown'], ['__proto__', 'unknown']]) {
    const summary = telemetry([{type: 'assistant', error}]);
    assert.equal(summary.failure_reason, reason);
    const record = statusRecord(context, env({REVIEW_RESULT: 'failure',
      REVIEW_TELEMETRY: JSON.stringify(summary)}));
    check(record);
    assert.equal(record.failure_reason, reason === 'unknown' ? 'review_failed' : reason);
  }
  assert.equal(telemetry([{type: 'assistant', message: {content: 'rate_limit'}}])
    .failure_reason, 'unknown');
  assert.equal(telemetry([{type: 'assistant', error: 'rate_limit'},
    {type: 'assistant', message: {content: 'Recovered'}}]).failure_reason, 'unknown');
  assert.equal(statusRecord(context, env({PUBLICATION_OUTCOME: 'failure',
    REVIEW_TELEMETRY: '{"failure_reason":"rate_limited"}'})).failure_reason, 'review_failed');
});

test('timeout needs explicit evidence; missing telemetry remains unknown or generic', () => {
  const timeout = statusRecord(context, env({REVIEW_RESULT: 'failure'}), true);
  check(timeout);
  assert.equal(timeout.failure_reason, 'timeout');
  assert.equal(statusRecord(context, env({REVIEW_RESULT: 'cancelled',
    PUBLICATION_OUTCOME: 'skipped', PUBLICATION_STATUS: ''})).failure_reason, 'unknown');
  assert.equal(statusRecord(context, env({REVIEW_RESULT: 'failure',
    REVIEW_TELEMETRY: 'not json'})).failure_reason, 'review_failed');
});

test('usage is numeric-only; raw responses, prompts and auth material are excluded', () => {
  const raw = [{type: 'assistant', message: {content: 'SECRET'}},
    {type: 'result', total_cost_usd: 0.2, errors: ['SECRET'], usage: {
      input_tokens: 123, output_tokens: 45, cache_read_input_tokens: 0,
      cache_creation_input_tokens: -1, authentication: 'SECRET'}}];
  const summary = telemetry(raw);
  assert.deepEqual(summary.usage, {input_tokens: 123, output_tokens: 45,
    cache_read_input_tokens: 0, cost_usd: 0.2});
  const record = statusRecord(context, env({REVIEW_TELEMETRY: JSON.stringify(summary)}));
  check(record);
  assert.deepEqual(record.usage, summary.usage);
  assert.ok(!JSON.stringify(record).includes('SECRET'));
  assert.equal(telemetry([{type: 'result', usage: {input_tokens: 'SECRET'},
    total_cost_usd: Infinity}]).usage, null);
});

test('missing execution file is safe and finalisation writes after lookup failure', async () => {
  const tmp = fs.mkdtempSync(path.join(os.tmpdir(), 's4f-review-status-'));
  try {
    const outputs = {};
    const core = {setOutput: (k, v) => {outputs[k] = v;}, notice: () => {}};
    capture(core, path.join(tmp, 'missing'));
    assert.deepEqual(JSON.parse(outputs.telemetry), {failure_reason: 'unknown', usage: null});
    const github = {rest: {actions: {listJobsForWorkflowRunAttempt: {}}},
      paginate: async () => {throw new Error('Sensitive API response');}};
    await finalise({github, context, core}, env({RUNNER_TEMP: tmp, REVIEW_RESULT: 'failure'}));
    const record = JSON.parse(fs.readFileSync(path.join(tmp, 'review-status.json')));
    check(record);
    assert.equal(record.status, 'failed');
    github.paginate = async () => [{name: 'review', conclusion: 'timed_out'}];
    await finalise({github, context, core}, env({RUNNER_TEMP: tmp, REVIEW_RESULT: 'failure'}));
    assert.equal(JSON.parse(fs.readFileSync(path.join(tmp, 'review-status.json'))).failure_reason,
      'timeout');
    github.paginate = async () => [];
    await finalise({github, context, core}, env({RUNNER_TEMP: tmp}));
    const completed = JSON.parse(fs.readFileSync(path.join(tmp, 'review-status.json')));
    check(completed);
    assert.equal(completed.status, 'completed');
  } finally {
    fs.rmSync(tmp, {recursive: true, force: true});
  }
});

test('publisher confirms completed only after successful API write', async () => {
  const data = {risk: 'LOW', summary: 'Summary', reason: 'Reason', findings: [],
    verification: 'Evidence', human_review: ['Check evidence']};
  const github = {rest: {pulls: {get: async () => ({data: {...context.payload.pull_request,
    state: 'open', draft: false}})}, issues: {listComments: {},
    createComment: async () => {throw new Error('Denied');}}}, paginate: async () => []};
  const core = {setFailed: () => {}, notice: () => {}};
  const params = {...env(), REVIEW_JSON: JSON.stringify(data), SECURITY_CHANGE: 'false'};
  await assert.rejects(publish({github, context, core}, params), /Denied/);
  github.rest.issues.createComment = async () => {};
  assert.deepEqual(await publish({github, context, core}, params), {status: 'completed'});
  params.REVIEW_JSON = 'invalid';
  assert.deepEqual(await publish({github, context, core}, params), {status: 'failed'});
});

test('workflow finalises skipped/failed attempts and uploads only the status file', () => {
  const yaml = fs.readFileSync('.github/workflows/claudeReview.yml', 'utf8');
  assert.match(yaml, /publish:\n    needs: review\n    if: \$\{\{ always\(\) \}\}/);
  assert.match(yaml, /Write machine-readable attempt status\n        if: \$\{\{ always\(\) \}\}/);
  assert.match(yaml, /Upload attempt status\n        if: \$\{\{ always\(\) \}\}/);
  assert.match(yaml, /path: \$\{\{ runner.temp \}\}\/review-status.json/);
  assert.match(yaml, /if-no-files-found: error/);
  assert.match(yaml, /retention-days: 30/);
  assert.ok(!yaml.includes('continue-on-error'));
});
