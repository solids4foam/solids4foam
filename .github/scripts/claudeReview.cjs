// Trusted base-branch code only. PR material is data, never executable code.
const fs = require('node:fs');
const risks = ['LOW', 'MEDIUM', 'HIGH'];
// Size limits, in characters, for every model-written field. They are chosen so
// that ANY report the schema accepts also passes validate() and can be
// published (see maxReportBytes and fitComment); a schema-valid report is never
// discarded for its size.
const limits = {summary: 2000, reason: 1000, verification: 2000, file: 300,
  claim: 800, evidence: 1200, suggested_check: 600, checklist: 800,
  findings: 10, checklistItems: 5};
const findingProperties = {
  kind: {type: 'string', enum: ['defect', 'risk', 'question']},
  severity: {type: 'string', enum: ['blocking', 'important', 'minor']},
  category: {type: 'string', enum: ['numerics', 'architecture', 'compatibility',
    'verification', 'security', 'build', 'style']},
  file: {type: ['string', 'null'], maxLength: limits.file},
  line: {type: ['integer', 'null'], minimum: 1},
  ...Object.fromEntries(['claim', 'evidence', 'suggested_check'].map(key =>
    [key, {type: 'string', minLength: 1, maxLength: limits[key]}]))
};
const schema = {
  type: 'object', additionalProperties: false,
  required: ['risk', 'summary', 'reason', 'findings', 'verification', 'human_review'],
  properties: {
    risk: {type: 'string', enum: risks},
    ...Object.fromEntries(['summary', 'reason', 'verification'].map(key =>
      [key, {type: 'string', minLength: 1, maxLength: limits[key]}])),
    findings: {type: 'array', maxItems: limits.findings, items: {type: 'object',
      additionalProperties: false, required: Object.keys(findingProperties),
      properties: findingProperties}},
    human_review: {type: 'array', minItems: 1, maxItems: limits.checklistItems,
      items: {type: 'string', minLength: 1, maxLength: limits.checklist}}
  }
};
// Most model-written characters a schema-valid report can hold (38,000).
const maxReportChars = limits.summary + limits.reason + limits.verification +
  limits.findings*(limits.file + limits.claim + limits.evidence + limits.suggested_check) +
  limits.checklistItems*limits.checklist;
// The report reaches the publisher as an environment variable, and Linux
// rejects a single variable over 128 KiB. Up to 3 bytes per character (any
// non-ASCII text such as Greek symbols) plus keys, enums and whitespace keeps
// every schema-valid report below that limit.
const maxReportBytes = 3*maxReportChars + 6000;
function exactKeys(data, keys) {
  return data && typeof data === 'object' && !Array.isArray(data) &&
    Object.keys(data).sort().join() === keys.slice().sort().join();
}
function validText(value, max) {
  return typeof value === 'string' && value.trim() && value.length <= max;
}
function validate(raw) {
  if (typeof raw !== 'string' || Buffer.byteLength(raw) > maxReportBytes) throw new Error('Invalid report size');
  const data = JSON.parse(raw);
  if (!exactKeys(data, schema.required) || !risks.includes(data.risk)) throw new Error('Invalid report');
  for (const key of ['summary', 'reason', 'verification']) {
    if (!validText(data[key], limits[key])) throw new Error('Invalid report text');
  }
  if (!Array.isArray(data.findings) || data.findings.length > limits.findings) throw new Error('Invalid findings');
  for (const finding of data.findings) {
    if (!exactKeys(finding, Object.keys(findingProperties))) throw new Error('Invalid finding');
    for (const key of ['kind', 'severity', 'category']) {
      if (!findingProperties[key].enum.includes(finding[key])) throw new Error('Invalid finding type');
    }
    for (const key of ['claim', 'evidence', 'suggested_check']) {
      if (!validText(finding[key], limits[key])) throw new Error('Invalid finding text');
    }
    if (finding.file !== null && !validText(finding.file, limits.file)) throw new Error('Invalid file');
    if (finding.line !== null && (!Number.isSafeInteger(finding.line) || finding.line < 1)) {
      throw new Error('Invalid line');
    }
  }
  if (!Array.isArray(data.human_review) || !data.human_review.length ||
      data.human_review.length > limits.checklistItems ||
      data.human_review.some(item => !validText(item, limits.checklist))) throw new Error('Invalid checklist');
  return data;
}

// Render model-controlled text as inert content inside fixed report formatting.
// `max` shortens the raw text first (never an escape sequence); see fitComment.
function plain(value, max = Infinity) {
  if (value.length > max) value = Array.from(value).slice(0, max).join('') + '…';
  return value.replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;')
    .replace(/@/g, '@\u200b').replace(/[\\`*_{}[\]()#+.!|~-]/g, '\\$&')
    .replace(/[\r\n]+/g, ' ');
}
function classifyRisk(proposed, floor) {
  if (!risks.includes(proposed) || !risks.includes(floor)) return 'HIGH';
  return risks[Math.max(risks.indexOf(proposed), risks.indexOf(floor))];
}
function pathRisk(path) {
  if (/^(\.github\/|\.claude\/|\.agents\/|\.docker\/|tests\/|applications\/scripts\/)/.test(path) ||
      /(^|\/)(CLAUDE(\.local)?\.md|AGENTS\.md|SECURITY\.md|CODEOWNERS|\.gitmodules|\.mcp\.json|\.claude\.json|Allwmake|Allwclean)$/.test(path)) {
    return 'HIGH';
  }
  if (/^src\/(blockCoupledSolids4FoamTools\/|immersedBoundary\/|RBFMeshMotionSolver\/|solids4FoamModels\/(numerics|mechanicalConstitutiveLaw|materialModels|physicsModel|fluidSolidInterfaces|interfaceToInterfaceMapping|Make)[^/]*\/)/.test(path) ||
      (/^src\/solids4FoamModels\/(solidModels|fluidModels)\//.test(path) &&
       !/\/(fvPatchFields|pointPatchFields)\//.test(path)) ||
      /^applications\/solvers\//.test(path)) return 'HIGH';
  if (path.startsWith('src/')) return 'MEDIUM';
  return /\.(md|rst|txt|png|jpg|jpeg|svg|pdf)$/i.test(path) ? 'LOW' : 'MEDIUM';
}
function riskFloor(files) {
  return files.reduce((floor, file) => [file.filename, file.previous_filename].filter(Boolean)
    .reduce((risk, path) => classifyRisk(pathRisk(path), risk), floor), 'LOW');
}
function formatAssessment(data, risk, floor, max = Infinity) {
  const text = value => plain(value, max);
  const findings = data.findings.map(f => {
    const location = f.file === null ? 'PR-wide' : text(f.file) + (f.line === null ? '' : `:${f.line}`);
    return `- **${f.severity} ${f.kind} / ${f.category}** (${location})\n` +
      `  - Claim: ${text(f.claim)}\n  - Evidence: ${text(f.evidence)}\n` +
      `  - Suggested check: ${text(f.suggested_check)}`;
  }).join('\n');
  return `**Status: analysis completed. Advisory only; no approval.**\n\n` +
    (max === Infinity ? '' : '_Some report text was shortened to fit the GitHub comment ' +
      'size limit; no finding was removed._\n\n') +
    `**Risk: ${risk}** (path floor: ${risks.includes(floor) ? floor : 'HIGH, unavailable'})\n\n` +
    `${text(data.reason)}\n\n### Summary\n\n${text(data.summary)}\n\n### Findings\n\n` +
    (findings || 'No actionable issue identified in the supplied evidence.') +
    `\n\n### Verification evidence and limits\n\n${text(data.verification)}\n\n` +
    `### Human review checklist\n\n` + data.human_review.map(item => `- [ ] ${text(item)}`).join('\n') +
    '\n\nContributing reviewer: what did you independently check? ' +
    'Use this standard AI report; core maintainers decide approval and merging.';
}

// MEDIUM/HIGH need a second, independent AI review (model-policy.md). Until it
// is automated, the report asks a core maintainer to run it locally.
function secondReviewNotice(risk, number) {
  if (risk === 'LOW') return '';
  return '\n\n### Second independent AI review\n\n' +
    `Required for ${risk === 'MEDIUM' ? 'MEDIUM' : 'HIGH'} risk, but not yet automated. ` +
    'A core maintainer runs it locally from an up-to-date `development` checkout with ' +
    `\`node .github/scripts/secondReview.cjs ${number} --post\`. It gives Codex ` +
    '(gpt-6-astra) the same instructions, repository guide and patches, without this ' +
    'report or any PR comment, and posts its result as a separate comment. This report ' +
    'alone does not meet the two-reviewer policy.\n\n' +
    '- [ ] Second independent AI review posted (core maintainer)';
}
function formatCoverage(coverage, max = Infinity, listed = Infinity) {
  const omitted = coverage.omitted.slice(0, listed);
  const more = coverage.omitted.length - omitted.length;
  return `**Patch coverage: ${coverage.status.toUpperCase()}** ` +
    '(supplied hunks only; surrounding source and linked evidence are not inspected).\n\n' +
    (omitted.length ? 'Not reviewed:\n' + omitted.map(item =>
      `- ${plain(item.path, max)} (${item.reason === 'no_patch' ? 'no patch available' : 'prompt budget'})`)
      .join('\n') + (more ? `\n- and ${more} more files` : '') + '\n\n' : '');
}

// GitHub rejects comments over 65,536 characters, and escaping can expand
// model text up to five times. Rather than discard a valid report, shorten
// each text field (and the not-reviewed list) step by step until it fits. The
// last step fits any schema-valid report: 48 fields of 121 characters, escaped,
// plus 20 listed paths, is under 45,000 characters.
const commentLimit = 60000;
const fitSteps = [[Infinity, Infinity], [1000, 100], [500, 50], [250, 30], [120, 20]];
function fitComment(render) {
  for (const [max, listed] of fitSteps) {
    const body = render(max, listed);
    if (body.length <= commentLimit) return body;
  }
  throw new Error('Report too large');
}

// The agent skill shared by all solids4foam agents, read from the trusted base
// branch on every run (never copied), so the reviewer knows the same rules.
// `read` defaults to the trusted checkout; the local second review passes a
// reader that fetches the same files from the PR's base branch on GitHub.
const skillPath = '.agents/skills/solids4foam/SKILL.md';
const readCheckout = path => fs.readFileSync(path, 'utf8');
function repositoryGuide(read = readCheckout) {
  const skill = read(skillPath).replace(/^---\n[\s\S]*?\n---\n+/, '');
  const release = read('CHANGELOG.md').match(/^## \[(v[^\]]+)\]/m);
  return `## Repository guide (trusted, from \`${skillPath}\` on the base branch)\n\n` +
    (release ? `Latest released version in CHANGELOG.md: ${release[1]}\n\n` : '') +
    skill.trim();
}

// Budget the JSON-escaped *whole prompt*, leaving room for ALL_INPUTS and schema.
function buildPrompt(pr, files, read = readCheckout) {
  const instructions = read('.github/review/instructions.md');
  const guide = repositoryGuide(read);
  const selected = [], omitted = [];
  const priority = file => file.filename.startsWith('src/') ? 0 :
    file.filename.startsWith('applications/') ? 1 :
    /(^|\/)(polyMesh|0)\/|\.csv$/.test(file.filename) ? 3 : 2;
  const prompt = () => `${instructions}\n\n${guide}\n\n` +
    `Reviewed head: ${pr.head.sha}\nBase: ${pr.base.sha}\n` +
    `\nUNTRUSTED PR DATA (JSON):\n${JSON.stringify({title: pr.title, body: pr.body,
      files: selected, not_reviewed: omitted})}`;
  // Reserve manifest space before selecting any patches; never silently truncate a patch.
  // The limit covers instructions and guide (~18 KB escaped) plus the file manifest.
  for (const file of files) omitted.push({path: file.filename,
    reason: file.patch ? 'budget' : 'no_patch'});
  if (Buffer.byteLength(JSON.stringify(prompt())) > 45000) throw new Error('metadata_too_large');
  for (const file of files.slice().sort((a, b) => priority(a) - priority(b))) {
    if (!file.patch) continue;
    const item = {filename: file.filename, previous_filename: file.previous_filename,
      status: file.status, patch: file.patch};
    const index = omitted.findIndex(entry => entry.path === file.filename);
    const [entry] = omitted.splice(index, 1);
    selected.push(item);
    if (Buffer.byteLength(JSON.stringify(prompt())) > 90000) {
      selected.pop();
      omitted.push(entry);
    }
  }
  if (!selected.length) throw new Error('nothing_reviewable');
  return {prompt: prompt(), coverage: {status: omitted.length ? 'partial' : 'complete', omitted}};
}
async function prepare({github, context, core}) {
  try {
    const pr = context.payload.pull_request;
    const args = {...context.repo, pull_number: pr.number};
    const {data: current} = await github.rest.pulls.get(args);
    if (current.head.sha !== pr.head.sha || current.state !== 'open' || current.draft) {
      throw new Error('pr_changed_or_ineligible');
    }
    const files = await github.paginate(github.rest.pulls.listFiles, {...args, per_page: 100});
    if (!files.length || files.length >= 3000) throw new Error('file_list_incomplete');
    const result = buildPrompt(current, files);
    const {data: after} = await github.rest.pulls.get(args);
    if (after.head.sha !== pr.head.sha || after.base.sha !== current.base.sha) {
      throw new Error('pr_changed_during_collection');
    }
    core.setOutput('prompt', result.prompt);
    core.setOutput('coverage', JSON.stringify(result.coverage));
    core.setOutput('base_sha', current.base.sha);
    core.setOutput('schema', JSON.stringify(schema));
    core.setOutput('risk_floor', riskFloor(files));
  } catch (error) {
    const reasons = ['metadata_too_large', 'nothing_reviewable', 'pr_changed_or_ineligible',
      'file_list_incomplete', 'pr_changed_during_collection'];
    core.setOutput('collection_reason', reasons.includes(error.message) ? error.message : 'collection_failed');
    throw new Error('PR evidence collection failed; see the machine-readable collection reason.');
  }
}
async function publish({github, context, core}, env = process.env) {
  const pr = context.payload.pull_request;
  const args = {...context.repo, pull_number: pr.number};
  const {data: current} = await github.rest.pulls.get(args);
  if (current.head.sha !== pr.head.sha || current.state !== 'open' || current.draft ||
      (env.REVIEW_BASE_SHA && current.base.sha !== env.REVIEW_BASE_SHA)) {
    core.notice('No report published: the PR changed or is no longer eligible.');
    return {status: 'skipped'};
  }
  const run = `${context.serverUrl}/${context.repo.owner}/${context.repo.repo}/actions/runs/${context.runId}`;
  const heading = `## Claude advisory review\n\nReviewed commit: \`${pr.head.sha}\`\n\n` +
    '**Applies only if this SHA matches the current PR head.**\n\n' +
    `[Workflow run](${run}) · Attempt ${env.GITHUB_RUN_ATTEMPT || '1'}\n\n` +
    `Requested API model: ${plain(env.CLAUDE_REVIEW_MODEL || 'unspecified')}; ` +
    `effort: ${plain(env.CLAUDE_REVIEW_EFFORT || 'unspecified')}.\n\n`;
  let body, coverage;
  let failed = false;
  try {
    if (env.REVIEW_RESULT !== 'success' || env.CLAUDE_CONCLUSION !== 'success') throw new Error('execution_failed');
    const data = validate(env.REVIEW_JSON);
    coverage = JSON.parse(env.REVIEW_COVERAGE);
    if (!['complete', 'partial'].includes(coverage.status) || !Array.isArray(coverage.omitted) ||
        Buffer.byteLength(env.REVIEW_COVERAGE) > 30000) throw new Error('Invalid coverage');
    const risk = classifyRisk(data.risk, env.RISK_FLOOR);
    body = fitComment((max, listed) => heading + formatCoverage(coverage, max, listed) +
      formatAssessment(data, risk, env.RISK_FLOOR, max) + secondReviewNotice(risk, pr.number));
  } catch {
    failed = true;
    const known = ['metadata_too_large', 'nothing_reviewable', 'file_list_incomplete',
      'pr_changed_or_ineligible', 'pr_changed_during_collection', 'collection_failed'];
    const reason = known.includes(env.COLLECTION_REASON) ? env.COLLECTION_REASON : 'execution_or_report_failed';
    body = heading + '**Status: INCOMPLETE — no successful AI review.**\n\n' +
      '**Risk: HIGH (conservative fallback; classification unassessed).**\n\n' +
      `Reason: ${reason}. Human review is required. See the workflow run for details.` +
      secondReviewNotice('HIGH', pr.number);
  }
  // One sticky conversation comment, never an approval or merge API call.
  const marker = '<!-- solids4foam-ai-review -->';
  const comments = await github.paginate(github.rest.issues.listComments,
    {...context.repo, issue_number: pr.number, per_page: 100});
  const previous = comments.find(comment => comment.user?.login === 'github-actions[bot]' &&
    comment.body?.startsWith(marker));
  body = `${marker}\n${body}`;
  if (previous) await github.rest.issues.updateComment({...context.repo, comment_id: previous.id, body});
  else await github.rest.issues.createComment({...context.repo, issue_number: pr.number, body});
  if (failed) core.setFailed('AI review incomplete; failure report published.');
  return {status: failed ? 'failed' : 'completed', coverage: failed ? 'incomplete' : coverage.status};
}
module.exports = {schema, limits, maxReportChars, maxReportBytes, commentLimit, fitComment, validate, classifyRisk, pathRisk, riskFloor, buildPrompt,
  formatAssessment, formatCoverage, secondReviewNotice, plain, prepare, publish};
