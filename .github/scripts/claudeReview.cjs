// Trusted base-branch code only. PR material is data, never executable code.
const fs = require('node:fs');

const schema = {
  type: 'object', additionalProperties: false,
  required: ['risk', 'summary', 'reason', 'findings', 'verification', 'human_review'],
  properties: {
    risk: {type: 'string', enum: ['LOW', 'MEDIUM', 'HIGH']},
    ...Object.fromEntries(['summary', 'reason', 'verification'].map(key =>
      [key, {type: 'string', minLength: 1, maxLength: 4000}])),
    ...Object.fromEntries(['findings', 'human_review'].map(key =>
      [key, {type: 'array', maxItems: key === 'findings' ? 15 : 5,
        items: {type: 'string', minLength: 1, maxLength: 2000}}]))
  }
};

function validate(raw) {
  const data = JSON.parse(raw);
  if (!data || typeof data !== 'object' || Array.isArray(data) ||
      Object.keys(data).sort().join() !== schema.required.slice().sort().join()) {
    throw new Error('Invalid report fields');
  }
  if (!schema.properties.risk.enum.includes(data.risk)) throw new Error('Invalid risk');
  for (const key of ['summary', 'reason', 'verification']) {
    if (typeof data[key] !== 'string' || !data[key].trim() || data[key].length > 4000) {
      throw new Error('Invalid report text');
    }
  }
  for (const key of ['findings', 'human_review']) {
    if (!Array.isArray(data[key]) || data[key].length > schema.properties[key].maxItems ||
        data[key].some(item => typeof item !== 'string' || !item.trim() || item.length > 2000)) {
      throw new Error('Invalid report list');
    }
  }
  if (!data.human_review.length) throw new Error('Missing human review checklist');
  return data;
}

// Do not render model-generated HTML, links, images, or mentions as active content.
function plain(value) {
  return value.replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;')
    .replace(/@/g, '@\u200b').replace(/[\\`*_{}[\]()#+.!|~-]/g, '\\$&')
    .replace(/[\r\n]+/g, ' ');
}

async function prepare({github, context, core}) {
  const pr = context.payload.pull_request;
  const args = {...context.repo, pull_number: pr.number};
  const {data: current} = await github.rest.pulls.get(args);
  if (current.head.sha !== pr.head.sha || current.state !== 'open' || current.draft) {
    throw new Error('PR changed or is no longer eligible');
  }
  const files = await github.paginate(github.rest.pulls.listFiles, {...args, per_page: 100});
  // GitHub caps this endpoint at 3000 files and may omit binary/large patches.
  if (!files.length || files.length >= 3000 || files.some(file => !file.patch)) {
    throw new Error('Incomplete diff; manual review required');
  }
  const payload = JSON.stringify({title: current.title, body: current.body,
    files: files.map(({filename, previous_filename, status, patch}) =>
      ({filename, previous_filename, status, patch}))});
  if (Buffer.byteLength(payload) > 120000) throw new Error('Diff exceeds review limit');
  // Guard against a synchronize event during pagination.
  const {data: after} = await github.rest.pulls.get(args);
  if (after.head.sha !== pr.head.sha || after.base.sha !== current.base.sha) {
    throw new Error('Diff changed during collection');
  }
  const instructions = fs.readFileSync('.github/review/instructions.md', 'utf8');
  core.setOutput('prompt', `${instructions}\n\nReviewed head: ${pr.head.sha}\n` +
    `Base: ${current.base.sha}\n\nUNTRUSTED PR DATA (JSON):\n${payload}`);
  core.setOutput('schema', JSON.stringify(schema));
  // Security/CI changes must never be downgraded by the model.
  core.setOutput('security_change', files.some(file =>
    [file.filename, file.previous_filename || ''].some(path =>
      path.startsWith('.github/') || /(^|\/)(SECURITY\.md|CODEOWNERS)$/.test(path))));
}

async function publish({github, context, core}, env = process.env) {
  const pr = context.payload.pull_request;
  const args = {...context.repo, pull_number: pr.number};
  const {data: current} = await github.rest.pulls.get(args);
  if (current.head.sha !== pr.head.sha || current.state !== 'open' || current.draft) {
    core.notice('No report published: the PR changed or is no longer eligible.');
    return;
  }
  const run = `${context.serverUrl}/${context.repo.owner}/${context.repo.repo}` +
    `/actions/runs/${context.runId}`;
  const heading = `## Claude advisory review\n\nReviewed commit: \`${pr.head.sha}\`\n` +
    `\n[Workflow run](${run}) · Attempt ${env.GITHUB_RUN_ATTEMPT || '1'}\n\n`;
  let body;
  let failed = false;
  try {
    if (env.REVIEW_RESULT !== 'success' || env.CLAUDE_CONCLUSION !== 'success') {
      throw new Error('Review did not complete');
    }
    const data = validate(env.REVIEW_JSON);
    const risk = env.SECURITY_CHANGE === 'true' ? 'HIGH' : data.risk;
    body = heading + `**Status: analysis completed. Advisory only; no approval.**\n\n` +
      `**Risk: ${risk}**\n\n${plain(data.reason)}\n\n` +
      (env.SECURITY_CHANGE === 'true' ? 'CI or security configuration changed; risk is at least HIGH.\n\n' : '') +
      `### Summary\n\n${plain(data.summary)}\n\n### Findings\n\n` +
      (data.findings.map(item => `- ${plain(item)}`).join('\n') ||
        'No actionable issue identified in the supplied patches.') +
      `\n\n### Verification evidence and limits\n\n${plain(data.verification)}\n\n` +
      `### Human review checklist\n\n` +
      data.human_review.map(item => `- [ ] ${plain(item)}`).join('\n') +
      '\n\nContributing reviewers are welcome. Core maintainers decide approval and merging.';
  } catch {
    failed = true;
    body = heading + '**Status: INCOMPLETE — no successful AI review.**\n\n' +
      'Risk is unassessed. Human review is required. Check the workflow run for ' +
      'missing credentials, API errors, timeout, incomplete diff, or invalid output. ' +
      'A previous report applies only to its stated commit.';
  }
  // A conversation comment, never a GitHub approval or request-changes review.
  await github.rest.issues.createComment({...context.repo, issue_number: pr.number, body});
  if (failed) core.setFailed('AI review incomplete; failure report published.');
}

module.exports = {schema, validate, plain, prepare, publish};
