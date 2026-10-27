#!/usr/bin/env python3
"""Compare the coupled pipe-Windkessel runs with the exact solution and with
the stored results of the original hand-written Windkessel couplings (data
copied to verification/reference/pipeWindkessel).

    analyse.py [runs-dir]            (default: runs/)

For every runs/<scheme>/<order>/<dt>/ that has postProcessing results:
  * RMSE of the applied pressure Po against the closed form
        Pc(t) = A R/(1+(wRC)^2) [sin wt - wRC cos wt + wRC exp(-t/RC)]
        Po(t) = Pc(t) + Rch A sin wt          (A = 1, w = 2 pi, R = 10,
                                               Rch = 1, C = 0.01, Pc(0) = 0)
    over the logged steps 0 < t <= 1, exactly as the original study defined it;
  * the same RMSE from the original run, and the largest difference between
    the two runs' Pc and Po at every step;
  * the number of coupling updates per step, which must agree for implicit
    schemes (it counts outer correctors).
Also prints the order of accuracy between successive dt, new and old.

Standard library only.
"""
import math, os, sys, glob

OLD = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..",
                   "verification", "reference", "pipeWindkessel")
ORD = {1: "First", 2: "Second", 3: "Third", 4: "Fourth"}
R, RCH, C, A, W = 10.0, 1.0, 0.01, 1.0, 2 * math.pi


def exact(t):
    k = W * R * C
    pc = A * R / (1 + k * k) * (math.sin(W * t) - k * math.cos(W * t)
                                + k * math.exp(-t / (R * C)))
    return pc, pc + RCH * A * math.sin(W * t)


def read_new(path):
    names, rows = None, []
    for l in open(path):
        if l.startswith('#'):
            p = l[1:].split()
            if p and p[0] == 'time':
                names = p
            continue
        if l.strip():
            rows.append([float(x) for x in l.split()])
    ix = {n: i for i, n in enumerate(names)}
    return [(r[0], r[ix['C.P']], r[ix['outlet.Ppatch']], int(r[ix['updates']]))
            for r in rows]


def read_old(path):
    """Accepted values of every step, from the original integration log (the
    state log omits the final step): the last update of each time step."""
    last = {}
    for l in open(path.replace("windkesselPressure.dat", "windkesselIntegrationLog.dat")):
        if l.startswith('#') or not l.strip():
            continue
        p = l.split()
        # time timeIndex outerCorrector stepUpdate totalUpdates Qo0 Qo1 Pc Po
        last[float(p[0])] = (float(p[0]), float(p[7]), float(p[8]), int(p[3]))
    return [last[t] for t in sorted(last) if t > 0]


def rmse(rows):
    e = [(po - exact(t)[1]) ** 2 for t, pc, po, n in rows if 0 < t <= 1 + 1e-9]
    return math.sqrt(sum(e) / len(e)), len(e)


def old_path(scheme, order, dt):
    name = scheme.capitalize() + ORD[order] + "Order"
    p = os.path.join(OLD, name, dt, "windkesselPressure.dat")
    return p if os.path.exists(p) else None


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--tol")]
    tol = None
    for a in sys.argv[1:]:
        if a.startswith("--tol="):
            tol = float(a.split("=")[1])
    runs = args[0] if args else "runs"
    bad = 0
    found = {}
    for f in glob.glob(os.path.join(runs, "*", "*", "*", "postProcessing",
                                    "circulationCoupling", "*", "coupling.dat")):
        parts = f.split(os.sep)
        i = parts.index("postProcessing")
        scheme, order, dt = parts[i - 3], int(parts[i - 2]), parts[i - 1]
        found.setdefault((scheme, order), {})[dt] = f

    print("%-9s %-5s %-8s %5s | %-12s %-12s | %-10s %-10s | %s" % (
        "scheme", "order", "dt", "steps", "RMSE new", "RMSE old",
        "max|dPc|", "max|dPo|", "updates new/old"))
    summary = {}
    for (scheme, order), d in sorted(found.items()):
        for dt in sorted(d, key=float, reverse=True):
            new = read_new(d[dt])
            e, n = rmse(new)
            op = old_path(scheme, order, dt)
            if op:
                old = read_old(op)
                eo, _ = rmse(old)
                dn = {round(t, 9): (pc, po, u) for t, pc, po, u in old}
                dpc = dpo = 0.0
                same = 0
                cnt = 0
                for t, pc, po, u in new:
                    o = dn.get(round(t, 9))
                    if o is None:
                        continue
                    cnt += 1
                    dpc = max(dpc, abs(pc - o[0])); dpo = max(dpo, abs(po - o[1]))
                    same += (u == o[2])
                tail = "%d/%d rows equal" % (same, cnt)
                if tol is not None and (dpc > tol or dpo > tol or same != cnt):
                    bad += 1
                print("%-9s %-5d %-8s %5d | %-12.6g %-12.6g | %-10.3g %-10.3g | %s" % (
                    scheme, order, dt, n, e, eo, dpc, dpo, tail))
                summary.setdefault((scheme, order), []).append((float(dt), e, eo))
            else:
                print("%-9s %-5d %-8s %5d | %-12.6g %-12s | (no stored run)" % (
                    scheme, order, dt, n, e, "-"))
                summary.setdefault((scheme, order), []).append((float(dt), e, None))

    print("\norder of accuracy between successive dt (new | old):")
    for (scheme, order), pts in sorted(summary.items()):
        pts.sort(reverse=True)
        s = []
        for (d1, e1, o1), (d2, e2, o2) in zip(pts, pts[1:]):
            pn = math.log(e1 / e2) / math.log(d1 / d2)
            po = (math.log(o1 / o2) / math.log(d1 / d2)) if o1 and o2 else float('nan')
            s.append("%g->%g: %.3f|%.3f" % (d1, d2, pn, po))
        print("  %-9s order %d  expected %d   %s" % (scheme, order, order, "   ".join(s)))


    if tol is not None:
        print("%s: %d case(s) differ from the stored runs by more than %g Pa"
              % ("FAIL" if bad else "PASS", bad, tol))
        sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
