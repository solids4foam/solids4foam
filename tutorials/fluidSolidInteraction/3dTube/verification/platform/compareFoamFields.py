"""Compare two ASCII OpenFOAM fields: internal field and every patch value."""
import re, sys, math
def parse_list(txt, start):
    m = re.compile(r'(nonuniform\s+List<(\w+)>\s*(\d+)\s*\(|uniform\s+(\([^)]*\)|[-+0-9.eE]+))').search(txt, start)
    if m.group(4) is not None:
        v = m.group(4).strip('()').split(); return ('uniform', [tuple(map(float, v))]), m.end()
    n = int(m.group(3)); kind = m.group(2); i = m.end(); out = []
    if kind == 'scalar':
        body = txt[i:txt.index(')', i)]; out = [(float(x),) for x in body.split()[:n]]; i = txt.index(')', i)
    else:
        for _ in range(n):
            a = txt.index('(', i); b = txt.index(')', a); out.append(tuple(map(float, txt[a+1:b].split()))); i = b + 1
    return ('nonuniform', out), i
def read(path):
    t = open(path).read(); res = {}
    j = t.index('internalField'); res['internal'], _ = parse_list(t, j)
    bf = t.index('boundaryField')
    for m in re.finditer(r'\n    ([\w-]+)\s*\n    \{(.*?)\n    \}', t[bf:], re.S):
        body = m.group(2)
        k = body.find('value')
        if k >= 0:
            res[m.group(1)], _ = parse_list(body, k)
    return res
def cmp(a, b):
    (ka, va), (kb, vb) = a, b
    if len(va) != len(vb):
        if len(va) == 1: va = va * len(vb)
        elif len(vb) == 1: vb = vb * len(va)
    ncomp = len(va[0]); out = []
    for c in range(ncomp):
        d = [abs(x[c]-y[c]) for x, y in zip(va, vb)]
        mx = max(d); imax = d.index(mx); scale = max(max(abs(x[c]) for x in va), 1e-300)
        l2 = math.sqrt(sum(x*x for x in d)); l2ref = math.sqrt(sum(x[c]**2 for x in va)) or 1e-300
        out.append((c, mx, mx/scale, l2/l2ref, imax, va[imax][c], vb[imax][c]))
    return out
if __name__ == '__main__':
    A, B = read(sys.argv[1]), read(sys.argv[2])
    for key in A:
        if key not in B: continue
        for c, mx, rmx, rl2, im, x, y in cmp(A[key], B[key]):
            flag = '  <==' if rmx > 1e-9 else ''
            print(f"{key:12s} comp{c} maxAbs {mx:.3e} max/scale {rmx:.3e} relL2 {rl2:.3e} at {im} ({x:.6e} vs {y:.6e}){flag}")
