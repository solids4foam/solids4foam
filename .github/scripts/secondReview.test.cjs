const {test} = require('node:test');
const assert = require('node:assert/strict');
const {collect, formatSecond, marker} = require('./secondReview.cjs');
const {buildPrompt} = require('./claudeReview.cjs');

const pr = {number: 7, state: 'open', head: {sha: 'c'.repeat(40)},
  base: {ref: 'development', sha: 'd'.repeat(40)}, title: 'T', body: 'B'};
const files = [{filename: 'src/solids4FoamModels/numerics/x.C', status: 'modified',
  patch: '@@ -1 +1 @@\n-a\n+b'}];
function fakeApi(prData = pr) {
  const raws = [];
  return {raws,
    json: () => prData,
    list: () => files,
    raw: endpoint => {
      raws.push(endpoint);
      if (endpoint.includes('CHANGELOG.md')) return '## [v9.9] - date\n';
      if (endpoint.includes('SKILL.md')) return '---\nname: s\n---\n\nBASE SKILL\n';
      return 'BASE INSTRUCTIONS';
    }};
}

test('second review builds the CI prompt from the base branch on GitHub, not the local tree', () => {
  const api = fakeApi();
  const c = collect(7, api);
  assert.ok(api.raws.length === 3 && api.raws.every(e => e.endsWith('?ref=development')));
  assert.match(c.prompt, /^BASE INSTRUCTIONS\n\n## Repository guide/);
  assert.match(c.prompt, /Latest released version in CHANGELOG.md: v9\.9/);
  assert.match(c.prompt, /BASE SKILL/);
  assert.equal(c.prompt, buildPrompt(pr, files, endpoint => api.raw(`x/${endpoint}?ref=development`)).prompt);
  assert.equal(c.floor, 'HIGH');
  assert.ok(!c.prompt.includes('Claude advisory review'));
});

test('second review refuses PRs that are closed or target other branches', () => {
  assert.throws(() => collect(7, fakeApi({...pr, state: 'closed'})), /not open/);
  assert.throws(() => collect(7, fakeApi({...pr, base: {ref: 'feature-x'}})), /not open/);
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
