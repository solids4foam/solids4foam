const {test} = require('node:test');
const assert = require('node:assert/strict');
const fs = require('node:fs');
const os = require('node:os');
const path = require('node:path');
const {telemetry, statusRecord, finalise} = require('./claudeReviewStatus.cjs');
const {publish} = require('./claudeReview.cjs');

const context = {repo: {owner: 'example', repo: 'test'}, runId: 123,
  payload: {pull_request: {number: 42, head: {sha: 'a'.repeat(40)}}}};
const env = extra => ({GITHUB_RUN_ATTEMPT: '2', REVIEW_RESULT: 'success',
  CLAUDE_CONCLUSION: 'success', PUBLICATION_OUTCOME: 'success',
  PUBLICATION_STATUS: '{"status":"completed","coverage":"complete"}', ...extra});

function check(record) {
  assert.equal(record.schema_version, 2);
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
  assert.equal(completed.coverage, 'complete');
  assert.equal(statusRecord(context, env({PUBLICATION_STATUS:
    '{"status":"completed","coverage":"partial"}'})).coverage, 'partial');
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

test('only fixed analysis reason codes and HTTP status numbers reach the record', () => {
  for (const [reason, expected] of [['rate_limited', 'rate_limited'],
    ['authentication_error', 'authentication_error'], ['api_error', 'api_error'],
    ['refusal', 'refusal'], ['Credit balance is too low', 'review_failed'], ['__proto__', 'review_failed']]) {
    const record = statusRecord(context, env({REVIEW_RESULT: 'failure',
      REVIEW_TELEMETRY: JSON.stringify({failure_reason: reason, api_error_status: 400})}));
    check(record);
    assert.equal(record.failure_reason, expected);
    assert.equal(record.api_error_status, 400);
  }
  for (const bad of ['400', 4.5, -1, 99, 600, null]) {
    assert.equal(telemetry(JSON.stringify({api_error_status: bad})).api_error_status, null);
  }
  // A failed publication is not blamed on the analysis.
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

test('usage is numeric-only; model ids are validated; nothing else is copied', () => {
  const raw = JSON.stringify({failure_reason: null, served_model: 'claude-opus-5-5', attempts: 2,
    secret: 'SECRET', usage: {input_tokens: 123, output_tokens: 45, cache_read_input_tokens: 0,
      cache_creation_input_tokens: -1, authentication: 'SECRET'}});
  const summary = telemetry(raw);
  assert.deepEqual(summary.usage, {input_tokens: 123, output_tokens: 45, cache_read_input_tokens: 0});
  const record = statusRecord(context, env({REVIEW_TELEMETRY: raw}));
  check(record);
  assert.deepEqual(record.usage, summary.usage);
  assert.equal(record.served_model, 'claude-opus-5-5');
  assert.equal(record.model_requests, 2);
  assert.ok(!JSON.stringify(record).includes('SECRET'));
  assert.equal(telemetry(JSON.stringify({served_model: '<b>x</b>', usage: {input_tokens: 'SECRET'}}))
    .served_model, null);
  assert.equal(telemetry('{"usage":{"input_tokens":"SECRET"}}').usage, null);
});

test('finalisation writes a record after a job lookup failure', async () => {
  const tmp = fs.mkdtempSync(path.join(os.tmpdir(), 's4f-review-status-'));
  try {
    const outputs = {};
    const core = {setOutput: (k, v) => {outputs[k] = v;}, notice: () => {}};
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
  const params = {...env(), REVIEW_JSON: JSON.stringify(data), RISK_FLOOR: 'LOW',
    REVIEW_COVERAGE: JSON.stringify({status: 'complete', omitted: []})};
  await assert.rejects(publish({github, context, core}, params), /Denied/);
  github.rest.issues.createComment = async () => {};
  assert.deepEqual(await publish({github, context, core}, params), {status: 'completed', coverage: 'complete'});
  params.REVIEW_JSON = 'invalid';
  assert.deepEqual(await publish({github, context, core}, params), {status: 'failed', coverage: 'incomplete'});
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
