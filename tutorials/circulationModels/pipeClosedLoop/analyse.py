#!/usr/bin/env python3
"""Compare coupled closed-loop pipe runs with the stored results of the
original hand-written couplings (verification/reference/pipeClosedLoop) and
with the monolithic reference (pipe replaced by its R-L surrogate).

    analyse.py [runs-dir] [--tol=PCT]       (default runs/)

For every runs/<scheme>/<order>/<dt>/ with results it prints, per case, the
largest difference of each circuit state from the original run over the whole
simulation, as a percentage of the state's range, and the interface-flow and
interface-pressure differences from the 17-digit function-object records.
With --tol, exits non-zero if the worst state difference exceeds it (percent
of range).

Standard library only.
"""
import glob, math, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
OLD = os.path.join(HERE, "..", "..", "verification", "reference", "pipeClosedLoop")
ORD = {1: "First", 2: "Second", 3: "Third", 4: "Fourth"}

# original state-log column -> column of coupling.dat
STATES = [("Po", "C1.P"), ("Q1", "L1.Q"), ("P1", "C2.P"), ("P2", "C3.P"),
          ("Pa", "atrium.P"), ("Pv", "ventricle.P"), ("Va", "atrium.V"),
          ("Vv", "ventricle.V"), ("Pi", "C4.P"), ("Q2", "R2.Q"), ("Q3", "R3.Q"),
          ("Q4", "mitral.Q"), ("Q5", "aortic.Q")]


def table(path, key):
    names, rows = None, []
    for l in open(path):
        if l.startswith('#'):
            p = l[1:].split()
            if p and p[0] == key:
                names = p
            continue
        if l.strip():
            rows.append([float(x) for x in l.split()])
    return names, {round(r[0], 9): r for r in rows}


def monolithic():
    """Reference with the pipe replaced by its R-L surrogate: time -> (Po, Pi)."""
    path = os.path.join(HERE, "..", "..", "verification", "reference",
                        "closedLoop_monolithicRK4_dt6p25e-6.csv")
    out = {}
    if os.path.exists(path):
        for l in list(open(path))[1:]:
            p = l.split(",")
            out[round(float(p[0]), 9)] = (float(p[1]), float(p[7]))
    return out


def fo(path):
    d = {}
    for l in open(path):
        if l.startswith('#') or not l.strip():
            continue
        p = l.split()
        d[round(float(p[0]), 9)] = float(p[1])
    return d


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    tol = None
    for a in sys.argv[1:]:
        if a.startswith("--tol="):
            tol = float(a[6:])
    runs = args[0] if args else os.path.join(HERE, "runs")
    worst_all = 0.0
    mono = monolithic()
    summary = []
    print("%-9s %-5s %-9s %5s | %-11s %-11s | %-12s %-12s %-12s | %s" % (
        "scheme", "order", "dt", "steps", "worst state", "(which)",
        "max|dPi| Pa", "max|dPo| Pa", "max|dQ| m3/s", "updates equal"))
    for f in sorted(glob.glob(os.path.join(runs, "*", "*", "*", "postProcessing",
                                           "circulationCoupling", "*", "coupling.dat"))):
        parts = f.split(os.sep)
        i = parts.index("postProcessing")
        scheme, order, dt = parts[i - 3], int(parts[i - 2]), parts[i - 1]
        case = os.path.join(OLD, scheme + ORD[order] + "Order", "dt_" + dt)
        if not os.path.exists(os.path.join(case, "stateLog.dat")):
            print("%-9s %-5d %-9s  (no stored run)" % (scheme, order, dt))
            continue
        nn, new = table(f, "time")
        on, old = table(os.path.join(case, "stateLog.dat"), "time")
        ni = {n: k for k, n in enumerate(nn)}
        oi = {n: k for k, n in enumerate(on)}
        common = [t for t in sorted(old) if t in new]

        worst, wname = 0.0, ""
        for oname, nname in STATES:
            vals = [old[t][oi[oname]] for t in common]
            rng = max(vals) - min(vals)
            err = max(abs(new[t][ni[nname]] - old[t][oi[oname]]) for t in common)
            pct = 100 * err / rng if rng > 0 else 0.0
            if pct > worst:
                worst, wname = pct, oname
        same = sum(1 for t in common if int(new[t][ni["updates"]]) == int(old[t][oi["couplingSolveCount"]]))

        # interface records at full precision
        dpi = dpo = dq = float("nan")
        try:
            bp = fo(os.path.join(case, "fo_boundaryPressure.dat"))
            op = fo(os.path.join(case, "fo_outletPressure.dat"))
            of = fo(os.path.join(case, "fo_outletFlow.dat"))
            dpi = max(abs(new[t][ni["C4.P"]] - bp[t]*1060.0) for t in common if t in bp)
            dpo = max(abs(new[t][ni["C1.P"]] - op[t]*1060.0) for t in common if t in op)
            # The explicit schemes sample the flux at the start of the step; the
            # records hold the flux of each completed step, so shift by one.
            h = min(common[1] - common[0], 1.0) if len(common) > 1 else 0.0
            shift = h if scheme == "explicit" else 0.0
            dq = max(abs(new[t][ni["outlet.Qsampled"]] - of[round(t - shift, 9)])
                     for t in common if round(t - shift, 9) in of)
        except (OSError, KeyError):
            pass
        print("%-9s %-5d %-9s %5d | %-11.3g %-11s | %-12.3g %-12.3g %-12.3g | %d/%d" % (
            scheme, order, dt, len(common), worst, "(" + wname + " %)", dpi, dpo, dq, same, len(common)))
        worst_all = max(worst_all, worst)

        # distance from the monolithic reference over the final cycle, new and old
        if mono and any(t >= 1.6 - 1e-9 for t in common):
            fin = [t for t in common if t >= 1.6 - 1e-9 and t in mono]
            m_new = max(max(abs(new[t][ni["C1.P"]] - mono[t][0]), abs(new[t][ni["C4.P"]] - mono[t][1])) for t in fin)
            m_old = max(max(abs(old[t][oi["Po"]] - mono[t][0]), abs(old[t][oi["Pi"]] - mono[t][1])) for t in fin)
            summary.append((scheme, order, dt, m_new, m_old))

    if summary:
        print("\nmax |Po|,|Pi| distance from the monolithic R-L reference over the final cycle (Pa),"
              "\nthis code | original run:")
        for sch, o, dt, a, b in summary:
            print("  %-9s order %d  dt %-7s  %8.3f | %8.3f" % (sch, o, dt, a, b))

    if tol is not None:
        ok = worst_all <= tol
        print("%s: worst state difference %.3g %% of range (tolerance %g %%)"
              % ("PASS" if ok else "FAIL", worst_all, tol))
        sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
