#!/usr/bin/env node
// Interim, manual second AI review for MEDIUM/HIGH PRs (see
// .github/review/model-policy.md). A core maintainer runs it from an
// up-to-date development checkout:
//
//   node .github/scripts/secondReview.cjs <PR number> [--model fable] [--post]
//
// With --post it also raises the PR's ai-risk label if this review rates the
// risk higher than the label shows; it never lowers it.
//
// The default reviewer is Codex (gpt-6-astra, xhigh), a different model family
// from the CI report. --model fable uses Claude Fable 5.2 (high) through the
// Messages API instead: the same no-tool request, validation and fallback as
// CI, with the maintainer's ANTHROPIC_API_KEY. It shares the CI report's model
// family, so it is a substitute, not the preferred second reviewer. If Fable
// declines and the API's fallback serves another model (possibly the CI
// reviewer itself), the run fails rather than posting a non-independent review.
//
// --trusted-ref=<branch> reads the instructions and skill from another
// upstream branch instead of the PR's base; use it only for testing.
//
// The prompt is built exactly as in CI: the review instructions, the skill and
// CHANGELOG.md are fetched from the PR's base branch on GitHub (never from the
// local tree), and the PR title, body and patches are untrusted data. Codex
// never sees the Claude report or any PR comment. It runs in an empty
// temporary directory without tools (see runCodex), user config, MCP servers
// or rules. Needs `gh` (authenticated), and `codex` or, for --model fable, the
// pinned SDK (npm ci --prefix .github/scripts --ignore-scripts).
const {execFileSync} = require('node:child_process');
const fs = require('node:fs');
const os = require('node:os');
const path = require('node:path');
const review = require('./claudeReview.cjs');

const repo = 'solids4foam/solids4foam';
const model = 'gpt-6-astra';
const effort = 'xhigh';
const fable = {model: 'claude-fable-5-2', effort: 'high'};
const marker = review.secondMarker;

function run(cmd, args, input) {
  return execFileSync(cmd, args, {encoding: 'utf8', maxBuffer: 64 << 20, input,
    stdio: ['pipe', 'pipe', 'inherit']});
}
// Retry transient GitHub errors (for example HTTP 502/503) a few times.
function gh(args) {
  for (let attempt = 1; ; attempt++) {
    try {
      return run('gh', args);
    } catch (error) {
      if (attempt === 3) throw error;
      execFileSync('sleep', [String(5*attempt)]);
    }
  }
}
const githubApi = {
  json: endpoint => JSON.parse(gh(['api', endpoint])),
  list: endpoint => gh(['api', '--paginate', endpoint, '--jq', '.[] | tojson'])
    .split('\n').filter(Boolean).map(line => JSON.parse(line)),
  raw: endpoint => gh(['api', '-H', 'Accept: application/vnd.github.raw', endpoint])
};

// Same evidence and budget as prepare(); trusted files from the base branch.
// `trustedRef` overrides that branch only for testing before a merge, and must
// be a branch of the upstream repository that maintainers control.
function collect(number, api = githubApi, trustedRef = null) {
  const pr = api.json(`repos/${repo}/pulls/${number}`);
  if (pr.state !== 'open' || pr.draft || pr.base.ref !== 'development') {
    throw new Error('PR is not open, ready for review and against development');
  }
  const files = api.list(`repos/${repo}/pulls/${number}/files?per_page=100`);
  if (!files.length || files.length >= 3000) throw new Error('file_list_incomplete');
  const ref = encodeURIComponent(trustedRef || pr.base.ref);
  const read = file => api.raw(`repos/${repo}/contents/${file}?ref=${ref}`);
  const {prompt, system, data, coverage} = review.buildPrompt(pr, files, read);
  if (changed(pr, api.json(`repos/${repo}/pulls/${number}`))) {
    throw new Error('PR changed during collection; rerun');
  }
  return {pr, prompt, system, data, coverage, floor: review.riskFloor(files)};
}

// Like CI, a report only applies to the PR as collected: same head and base
// commits, same target branch, still open and not a draft.
function changed(pr, now) {
  return now.state !== 'open' || now.draft || now.head.sha !== pr.head.sha ||
    now.base.ref !== pr.base.ref || now.base.sha !== pr.base.sha;
}

// The PR text is untrusted, so Codex gets no way to read local files or reach
// the network: shell, code execution, apps, plugins, browser, image, agent and
// web tools are disabled (checked with codex-cli 0.167), on top of the
// read-only sandbox. Recheck this list when Codex is upgraded.
const disabledFeatures = ['shell_tool', 'unified_exec', 'code_mode_host', 'apps', 'plugins',
  'remote_plugin', 'browser_use', 'browser_use_external', 'computer_use', 'in_app_browser',
  'tool_suggest', 'skill_search', 'goals', 'image_generation', 'multi_agent', 'sleep_tool',
  'view_image'];
function runCodex(prompt, dir) {
  const schemaFile = path.join(dir, 'schema.json');
  const reportFile = path.join(dir, 'report.json');
  fs.writeFileSync(schemaFile, JSON.stringify(review.schema));
  run('codex', ['exec', '-m', model, '-c', `model_reasoning_effort="${effort}"`,
    '-c', 'web_search="disabled"', '-s', 'read-only', '-C', dir,
    '--skip-git-repo-check', '--ephemeral', '--ignore-user-config', '--ignore-rules',
    ...disabledFeatures.flatMap(feature => ['--disable', feature]),
    '--output-schema', schemaFile, '-o', reportFile, '-'], prompt);
  return fs.readFileSync(reportFile, 'utf8');
}

// Fable through the CI analysis function: system prompt and untrusted data as
// separate messages, no tools, structured output, server-side fallback.
async function runFable({system, data}, dir, env = process.env, analyse = review.analyse) {
  if (!env.ANTHROPIC_API_KEY) throw new Error('--model fable needs ANTHROPIC_API_KEY');
  const runEnv = {ANTHROPIC_API_KEY: env.ANTHROPIC_API_KEY, RUNNER_TEMP: dir,
    CLAUDE_REVIEW_MODEL: fable.model, CLAUDE_REVIEW_EFFORT: fable.effort};
  fs.writeFileSync(review.inputFile(runEnv), JSON.stringify({system, data}), {mode: 0o600});
  const outputs = {};
  let failure = 'analysis failed';
  await analyse({core: {setOutput: (key, value) => { outputs[key] = value; },
    setFailed: message => { failure = message; }}}, runEnv);
  if (outputs.conclusion !== 'success') throw new Error(failure);
  const served = JSON.parse(outputs.telemetry || '{}').served_model;
  if (served && served !== fable.model) {
    throw new Error(`${fable.model} declined and the fallback model ${served} answered. That is ` +
      'not an independent second review, so nothing was posted. Use Codex instead.');
  }
  return outputs.report;
}

// Validated and escaped exactly like the CI report.
function formatSecond({pr, coverage, floor}, raw, reviewer = null) {
  const data = review.validate(raw.trim());
  const risk = review.classifyRisk(data.risk, floor);
  return review.fitComment((max, listed) =>
    `${marker}\n## Second independent AI review (manual, interim)\n\n` +
    `Reviewed commit: \`${pr.head.sha}\`\n\n` +
    '**Applies only if this SHA matches the current PR head.**\n\n' +
    (reviewer ?
      `Model: ${reviewer.model} (Anthropic Messages API); effort: ${reviewer.effort}.` +
      ' Run locally by a core maintainer with the same instructions, repository guide and ' +
      'patches as the Claude report. It did not see that report or any PR comment, but it ' +
      'shares its model family, so overlapping findings are weaker evidence than agreement ' +
      'across families.\n\n' :
      `Model: ${model} (Codex CLI); reasoning effort: ${effort}. Run locally by a core ` +
      'maintainer with the same instructions, repository guide and patches as the Claude ' +
      'report. It did not see that report or any PR comment, so overlapping findings are ' +
      'independent agreement.\n\n') +
    review.formatCoverage(coverage, max, listed) + review.formatAssessment(data, risk, floor, max));
}

// One comment per maintainer: rerunning updates it instead of adding another.
function post(number, body, dir) {
  const me = githubApi.json('user').login;
  const comments = githubApi.list(`repos/${repo}/issues/${number}/comments?per_page=100`);
  const previous = comments.find(c => c.user?.login === me && c.body?.startsWith(marker));
  const file = path.join(dir, 'body.json');
  fs.writeFileSync(file, JSON.stringify({body}));
  const [method, endpoint] = previous ?
    ['PATCH', `repos/${repo}/issues/comments/${previous.id}`] :
    ['POST', `repos/${repo}/issues/${number}/comments`];
  // No retry: a repeated POST could add a duplicate comment.
  return JSON.parse(run('gh', ['api', '-X', method, endpoint, '--input', file])).html_url;
}

// The label change that raises the PR's ai-risk label to `risk`, or null if
// the label already shows that risk or higher. Other labels are untouched.
function labelRaise(names, risk) {
  const present = review.risks.filter(r => names.includes(review.riskLabels[r]));
  const shown = present.length ? present[present.length - 1] : null;
  if (shown && review.classifyRisk(shown, risk) === shown) return null;
  return {remove: present.map(r => review.riskLabels[r]), add: review.riskLabels[risk]};
}
function raiseLabel(number, risk, dir) {
  const names = githubApi.list(`repos/${repo}/issues/${number}/labels?per_page=100`).map(l => l.name);
  const change = labelRaise(names, risk);
  if (!change) return null;
  for (const name of change.remove) {
    run('gh', ['api', '-X', 'DELETE', `repos/${repo}/issues/${number}/labels/${encodeURIComponent(name)}`]);
  }
  const file = path.join(dir, 'labels.json');
  fs.writeFileSync(file, JSON.stringify({labels: [change.add]}));
  run('gh', ['api', '-X', 'POST', `repos/${repo}/issues/${number}/labels`, '--input', file]);
  return change.add;
}

// Options: --post, --model <astra|fable> or --model=<...>, --trusted-ref=<branch>.
// Anything else is an error, so a mistyped option never silently changes the run.
function parseArgs(argv) {
  const options = {model: 'astra', post: false, trustedRef: null, number: null};
  for (let i = 0; i < argv.length; i++) {
    const arg = argv[i];
    if (arg === '--post') options.post = true;
    else if (arg === '--model') options.model = argv[++i];
    else if (arg.startsWith('--model=')) options.model = arg.slice(8);
    else if (arg.startsWith('--trusted-ref=')) options.trustedRef = arg.slice(14) || null;
    else if (/^\d+$/.test(arg) && options.number === null) options.number = Number(arg);
    else return null;
  }
  if (!Number.isSafeInteger(options.number) || options.number < 1 ||
      !['astra', 'fable'].includes(options.model)) return null;
  return options;
}

async function main(argv) {
  const options = parseArgs(argv);
  if (!options) {
    console.error('Usage: node .github/scripts/secondReview.cjs <PR number> ' +
      '[--model astra|fable] [--post] [--trusted-ref=<branch>]');
    process.exit(2);
  }
  const {number, trustedRef} = options;
  const choice = options.model;
  const dir = fs.mkdtempSync(path.join(os.tmpdir(), 's4f-second-review-'));
  console.error(`Collecting PR #${number} from ${repo} (trusted files from its base branch)...`);
  if (trustedRef) console.error(`Testing: trusted files from ${trustedRef}, not the PR's base branch.`);
  const collected = collect(number, githubApi, trustedRef);
  let body, raw;
  if (choice === 'fable') {
    console.error(`Running ${fable.model} (${fable.effort}); this can take several minutes...`);
    raw = await runFable(collected, dir);
    body = formatSecond(collected, raw, fable);
  } else {
    console.error(`Running ${model} (${effort}); this can take several minutes...`);
    raw = runCodex(collected.prompt, dir);
    body = formatSecond(collected, raw);
  }
  fs.writeFileSync(path.join(dir, 'comment.md'), body);
  if (changed(collected.pr, githubApi.json(`repos/${repo}/pulls/${number}`))) {
    console.error('PR changed during the review (head, base, state or draft); not posting. ' +
      `Rerun the review. Report: ${dir}/comment.md`);
    process.exit(1);
  }
  if (options.post) {
    console.log(`Posted: ${post(number, body, dir)}`);
    const risk = review.classifyRisk(review.validate(raw.trim()).risk, collected.floor);
    const raised = raiseLabel(number, risk, dir);
    if (raised) console.log(`Raised the PR label to "${raised}".`);
  } else {
    console.log(`Report: ${dir}/comment.md\nTo publish it: gh pr comment ${number} ` +
      `-R ${repo} --body-file ${dir}/comment.md`);
  }
}

if (require.main === module) {
  main(process.argv.slice(2)).catch(error => { console.error(error.message); process.exit(1); });
}
module.exports = {collect, changed, formatSecond, runFable, labelRaise, parseArgs, marker, model,
  effort, fable};
