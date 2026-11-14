#!/usr/bin/env python3
"""Analyse the ALE flux tests. usage: analyze.py <build-label> [<build-label> ...]
Reads runs/<label>/<case>/postProcessing and log.solids4Foam. Prints a markdown table."""
import os, re, sys, math
import numpy as np

R = os.path.dirname(os.path.abspath(__file__))
A = 0.2 * 0.01          # duct cross section
TSTART = 1.0            # window start for the steady metrics


def dat(path):
    rows = []
    for ln in open(path):
        if ln.startswith("#") or not ln.strip():
            continue
        v = [float(x) for x in re.sub(r"[()]", " ", ln).split()]
        rows.append(v)
    return np.array(rows)


def pp(d, name, f):
    p = os.path.join(d, "postProcessing", name, "0", f)
    return dat(p) if os.path.exists(p) else None


def contmetrics(d):
    sl, gl, cum = [], [], None
    for ln in open(os.path.join(d, "log.solids4Foam")):
        m = re.search(r"continuity errors : sum local = ([-+0-9.e]+), global = ([-+0-9.e]+), cumulative = ([-+0-9.e]+)", ln)
        if m:
            sl.append(float(m.group(1))); gl.append(float(m.group(2))); cum = float(m.group(3))
    if not sl:
        return None, None
    return max(sl), cum


# case -> (w body velocity (0 if static), U expected lab frame, mesh outlet velocity function)
def spec(name):
    base = name.split("_dt2")[0].split("_h2")[0]
    if base.startswith("A"):
        return dict(w=0.0, kind="A")
    if base.startswith("B") or base.startswith("F"):
        return dict(w=0.1, kind="fixedmesh")
    if base.startswith("C"):
        return dict(w=0.0, kind="C", vm=0.1)
    if base.startswith("D"):
        return dict(w=0.1, kind="D", vm=0.1)
    if base.startswith("G"):
        return dict(w=0.05, kind="G", vm=0.1)
    if base.startswith("E"):
        return dict(w=0.1, kind="E")


def dtof(d):
    c = open(os.path.join(d, "system", "controlDict")).read()
    return float(re.search(r"deltaT\s+([0-9.e-]+);", c).group(1))


def analyse(d, name):
    s = spec(name)
    w = s["w"]
    dt = dtof(d)
    out = dict(case=name, dt=dt)
    umin, umax = pp(d, "uMin", "volFieldValue.dat"), pp(d, "uMax", "volFieldValue.dat")
    t = umax[:, 0]
    sel = t >= TSTART - 1e-9
    # lab-frame velocity deviation from the exact uniform w (x) and 0 (y)
    ex = np.maximum(np.abs(umax[sel, 1] - w), np.abs(umin[sel, 1] - w))
    ey = np.maximum(np.abs(umax[sel, 2]), np.abs(umin[sel, 2]))
    refv = max(abs(w), 0.1 if s["kind"] in ("C", "A") else abs(w))
    out["maxdUx/v"] = ex.max() / refv
    out["maxUy/v"] = ey.max() / refv
    if s["kind"] in ("fixedmesh", "A"):
        fl = [pp(d, n, "surfaceFieldValue.dat") for n in ("fz_up", "fz_a", "fz_b")]
        fo = pp(d, "fluxOut", "surfaceFieldValue.dat")
        fi = pp(d, "fluxIn", "surfaceFieldValue.dat")
        ref = 0.1 * A if s["kind"] == "fixedmesh" else 0.1 * A
        Q = [f[f[:, 0] >= TSTART - 1e-9, 1] for f in fl + [fo]]
        out["Qmean/wA"] = np.mean([q.mean() for q in Q]) / (w * A) if w else float("nan")
        out["Qripple/wA"] = max(q.max() - q.min() for q in Q) / (w * A) if w else float("nan")
        out["maxdQ/wA"] = max(np.abs(q - w * A).max() for q in Q) / (0.1 * A)
        # continuity between planes: in - out
        n = min(len(fi), len(fo))
        out["|Qin+Qout|/wA"] = np.abs(fi[:n, 1] + fo[:n, 1])[fi[:n, 0] >= TSTART].max() / (0.1 * A)
    else:
        ai, ao = pp(d, "absIn", "surfaceFieldValue.dat"), pp(d, "absOut", "surfaceFieldValue.dat")
        ri, ro = pp(d, "fluxIn", "surfaceFieldValue.dat"), pp(d, "fluxOut", "surfaceFieldValue.dat")
        sel = ao[:, 0] >= TSTART - 1e-9
        scale = 0.1 * A
        # absolute flux out of outlet = +wA, into inlet = -wA
        out["maxdQout/scale"] = np.abs(ao[sel, 1] - w * A).max() / scale
        out["maxdQin/scale"] = np.abs(ai[sel, 1] + w * A).max() / scale
        # relative boundary fluxes: sum = -dV/dt exactly
        if s["kind"] == "E":
            amp = 0.1
            tt = ro[:, 0]
            tk = np.arange(0, 4.2, 0.05)
            dd = lambda x: np.interp(x, tk, amp * np.sin(np.pi * tk / 4.0))
            dVdt = A * (dd(tt) - dd(tt - dt)) / dt
        else:
            dVdt = np.zeros(len(ro))
        gcl = (ri[:, 1] + ro[:, 1]) + dVdt
        out["maxGCL/scale"] = np.abs(gcl[ro[:, 0] >= TSTART - 1e-9]).max() / scale
    sl, cum = contmetrics(d)
    out["contSumLocalMax"] = sl
    out["contCum"] = cum
    # impermeability proxy: flux the plate must displace is w*A
    return out


if __name__ == "__main__":
    labels = sys.argv[1:]
    for lab in labels:
        print("\n### build:", lab)
        keys = []
        res = []
        for c in sorted(os.listdir(os.path.join(R, "runs", lab))):
            d = os.path.join(R, "runs", lab, c)
            try:
                r = analyse(d, c)
            except Exception as e:
                r = dict(case=c, error=str(e)[:60])
            res.append(r)
            for k in r:
                if k not in keys:
                    keys.append(k)
        print("| " + " | ".join(keys) + " |")
        print("|" + "---|" * len(keys))
        for r in res:
            print("| " + " | ".join(("%.3g" % r[k]) if isinstance(r.get(k), float) else str(r.get(k, "")) for k in keys) + " |")
