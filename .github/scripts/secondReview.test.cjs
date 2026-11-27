const {test} = require('node:test');
const assert = require('node:assert/strict');
const fs = require('node:fs');
const os = require('node:os');
const path = require('node:path');
const {collect, changed, formatSecond, runFable, labelRaise, parseArgs, marker, fable} =
  require('./secondReview.cjs');
const {buildPrompt, inputFile} = require('./claudeReview.cjs');

const pr = {number: 7, state: 'open', draft: false, head: {sha: 'c'.repeat(40)},
  base: {ref: 'development', sha: 'd'.repeat(40)}, title: 'T', body: 'B'};
const files = [{filename: 'src/solids4FoamModels/numerics/x.C', status: 'modified',
  patch: '@@ -1 +1 @@\n-a\n+b'}];
function fakeApi(prData = pr, contents = {}) {
  const raws = [], fetched = [];
  return {raws, fetched,
    file: endpoint => {
      fetched.push(endpoint);
      const name = decodeURIComponent(endpoint.replace(/^repos\/[^/]+\/[^/]+\/contents\//, '')
        .replace(/\?.*$/, ''));
      if (!(name in contents)) throw new Error('gh: Not Found (HTTP 404)');
      const bytes = Buffer.from(contents[name]);
      return {type: 'file', encoding: 'base64', size: bytes.length, content: bytes.toString('base64')};
    },
    json: () => prData,
    list: () => files,
    raw: endpoint => {
      raws.push(endpoint);
      if (endpoint.includes('CHANGELOG.md')) return '## [v9.9] - date\n';
      if (endpoint.includes('SKILL.md')) return '---\nname: s\n---\n\nBASE SKILL\n';
      return 'BASE INSTRUCTIONS';
    }};
}

test('second review builds the CI prompt from the base branch on GitHub, not the local tree', async () => {
  const api = fakeApi();
  const c = await collect(7, api);
  assert.ok(api.raws.length === 3 && api.raws.every(e => e.endsWith('?ref=development')));
  assert.match(c.prompt, /^BASE INSTRUCTIONS\n\n## Repository guide/);
  assert.match(c.prompt, /Latest released version in CHANGELOG.md: v9\.9/);
  assert.match(c.prompt, /BASE SKILL/);
  assert.equal(c.prompt, buildPrompt(pr, files, endpoint => api.raw(`x/${endpoint}?ref=development`)).prompt);
  assert.equal(c.floor, 'MEDIUM');
  assert.ok(!c.prompt.includes('Claude advisory review'));
});

test('second review adds full files from the PR head SHA, with encoded paths', async () => {
  const api = fakeApi(pr, {'src/solids4FoamModels/numerics/x.C': 'int whole() { return 1; }\n'});
  const c = await collect(7, api);
  assert.deepEqual(api.fetched, [`repos/solids4foam/solids4foam/contents/src/solids4FoamModels/numerics/x.C?ref=${pr.head.sha}`]);
  assert.ok(api.raws.length === 3);
  assert.deepEqual(c.coverage.full_files, ['src/solids4FoamModels/numerics/x.C']);
  assert.ok(JSON.parse(c.data.split('UNTRUSTED PR DATA (JSON):\n')[1]).full_files[0].content
    .includes('int whole()'));
  // An untrusted path cannot add query parameters or leave the contents endpoint.
  const odd = fakeApi();
  files.push({filename: 'src/a b/x?ref=main#y.C', status: 'modified', patch: '@@ -1 +1 @@\n-a\n+b'});
  try {
    await collect(7, odd);
  } finally {
    files.pop();
  }
  assert.ok(odd.fetched.some(e => e.includes('src/a%20b/x%3Fref%3Dmain%23y.C?ref=')));
});

test('second review refuses PRs that are closed, drafts or not into development', async () => {
  await assert.rejects(collect(7, fakeApi({...pr, state: 'closed'})), /not open/);
  await assert.rejects(collect(7, fakeApi({...pr, draft: true})), /not open/);
  await assert.rejects(collect(7, fakeApi({...pr, base: {ref: 'master', sha: pr.base.sha}})), /not open/);
  await assert.rejects(collect(7, fakeApi({...pr, base: {ref: 'feature-x'}})), /not open/);
});

test('second review is not posted if the PR changed in any way CI also checks', async () => {
  assert.equal(changed(pr, structuredClone(pr)), false);
  for (const now of [{...pr, state: 'closed'}, {...pr, draft: true},
    {...pr, head: {sha: 'e'.repeat(40)}}, {...pr, base: {...pr.base, ref: 'master'}},
    {...pr, base: {...pr.base, sha: 'f'.repeat(40)}}]) {
    assert.equal(changed(pr, now), true);
  }
  // During collection too: a retarget between the two reads aborts.
  let reads = 0;
  const api = {...fakeApi(), json: () => (reads++ ? {...pr, base: {...pr.base, ref: 'master'}} : pr)};
  await assert.rejects(collect(7, api), /changed during collection/);
});

test('second review report is validated, escaped and applies the path floor', () => {
  const raw = JSON.stringify({risk: 'LOW', summary: 'Ping @someone <b>x</b>', reason: 'R',
    findings: [], verification: 'V', human_review: ['Check']});
  const body = formatSecond({pr, coverage: {status: 'complete', omitted: []}, floor: 'HIGH'}, raw);
  assert.ok(body.startsWith(`${marker}\n## Second independent AI review`));
  assert.match(body, /gpt-6-astra/);
  assert.match(body, /\*\*Risk: HIGH\*\*/);
  assert.ok(!body.includes('@someone') && !body.includes('<b>'));
  assert.ok(!body.includes('Required for HIGH risk'));
  assert.throws(() => formatSecond({pr, coverage: {status: 'complete', omitted: []}, floor: 'LOW'},
    '{"risk":"LOW"}'));
});

const report = JSON.stringify({risk: 'LOW', summary: 'S', reason: 'R', findings: [],
  verification: 'V', human_review: ['Check']});

test('Fable second review runs the CI analysis with separate system and untrusted data', async () => {
  const dir = fs.mkdtempSync(path.join(os.tmpdir(), 's4f-fable-'));
  let seen;
  const analyse = async ({core}, env) => {
    seen = {env, input: JSON.parse(fs.readFileSync(inputFile(env), 'utf8'))};
    core.setOutput('report', report);
    core.setOutput('conclusion', 'success');
    core.setOutput('telemetry', JSON.stringify({served_model: 'claude-fable-5-2'}));
  };
  const c = await collect(7, fakeApi());
  const result = await runFable(c, dir, {ANTHROPIC_API_KEY: 'k', OTHER: 'x'}, analyse);
  assert.equal(result, report);
  assert.equal(seen.env.CLAUDE_REVIEW_MODEL, 'claude-fable-5-2');
  assert.equal(seen.env.CLAUDE_REVIEW_EFFORT, 'high');
  assert.equal(seen.env.OTHER, undefined);
  assert.equal(seen.input.system + '\n\n' + seen.input.data, c.prompt);
});

test('Fable second review needs a key; analysis failure and any fallback model fail', async () => {
  const dir = fs.mkdtempSync(path.join(os.tmpdir(), 's4f-fable-'));
  const c = await collect(7, fakeApi());
  await assert.rejects(runFable(c, dir, {}, async () => {}), /ANTHROPIC_API_KEY/);
  await assert.rejects(runFable(c, dir, {ANTHROPIC_API_KEY: 'k'}, async ({core}) => {
    core.setOutput('conclusion', 'failure');
    core.setFailed('AI analysis failed: timeout');
  }), /timeout/);
  await assert.rejects(runFable(c, dir, {ANTHROPIC_API_KEY: 'k'}, async ({core}) => {
    core.setOutput('report', report);
    core.setOutput('conclusion', 'success');
    core.setOutput('telemetry', JSON.stringify({served_model: 'claude-opus-5-5'}));
  }), /fallback model claude-opus-5-5 answered.*not an independent second review/);
});

test('Fable second review comment names the model and the weaker independence', () => {
  const body = formatSecond({pr, coverage: {status: 'complete', omitted: []}, floor: 'HIGH'},
    report, fable);
  assert.match(body, /Model: claude-fable-5-2 \(Anthropic Messages API\); effort: high/);
  assert.match(body, /shares its model family/);
  assert.ok(!body.includes('gpt-6-astra'));
});

test('options are parsed strictly; a mistyped option never changes the reviewer silently', () => {
  assert.deepEqual(parseArgs(['726']), {model: 'astra', post: false, trustedRef: null, number: 726});
  assert.equal(parseArgs(['726', '--model', 'fable', '--post']).model, 'fable');
  assert.equal(parseArgs(['--model=fable', '726']).model, 'fable');
  assert.equal(parseArgs(['726', '--trusted-ref=test-branch']).trustedRef, 'test-branch');
  for (const bad of [[], ['726', '--modle', 'fable'], ['726', '--model'], ['726', '--model=opus'],
    ['726', '727'], ['726', '--Post'], ['abc']]) {
    assert.equal(parseArgs(bad), null, bad.join(' '));
  }
});

test('a second review only ever raises the ai-risk label', () => {
  assert.deepEqual(labelRaise(['bug'], 'MEDIUM'), {remove: [], add: 'ai-risk: medium'});
  assert.deepEqual(labelRaise(['ai-risk: medium', 'bug'], 'HIGH'),
    {remove: ['ai-risk: medium'], add: 'ai-risk: high'});
  assert.equal(labelRaise(['ai-risk: high'], 'MEDIUM'), null);
  assert.equal(labelRaise(['ai-risk: medium'], 'MEDIUM'), null);
  assert.deepEqual(labelRaise(['ai-risk: low', 'ai-review: incomplete'], 'HIGH'),
    {remove: ['ai-risk: low'], add: 'ai-risk: high'});
});
