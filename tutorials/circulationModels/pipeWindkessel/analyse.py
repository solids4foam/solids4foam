#!/usr/bin/env python3
"""Compare the coupled pipe-Windkessel run with the exact solution.

    python3 analyse.py [case-dir]            (default: this directory)

Reads postProcessing/circulationCoupling/*/coupling.dat and prints the RMSE of the
applied outlet pressure Po over the logged steps 0 < t <= 1 against the closed form

    Pc(t) = A R/(1+(wRC)^2) [sin wt - wRC cos wt + wRC exp(-t/RC)]
    Po(t) = Pc(t) + Rch A sin wt          (A = 1, w = 2 pi, R = 10, Rch = 1, C = 0.01, Pc(0) = 0)

Standard library only.
"""
import glob, math, os, sys

R, RCH, C, A, W = 10.0, 1.0, 0.01, 1.0, 2 * math.pi


def exact(t):
    k = W * R * C
    pc = A * R / (1 + k * k) * (math.sin(W * t) - k * math.cos(W * t)
                                + k * math.exp(-t / (R * C)))
    return pc, pc + RCH * A * math.sin(W * t)


def read(path):
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
    return [(r[0], r[ix['outlet.Ppatch']]) for r in rows]


case = sys.argv[1] if len(sys.argv) > 1 else os.path.dirname(os.path.abspath(__file__))
files = glob.glob(os.path.join(case, "postProcessing", "circulationCoupling", "*", "coupling.dat"))
if not files:
    sys.exit("no postProcessing/circulationCoupling/*/coupling.dat in " + case)
rows = read(files[0])
err = [(po - exact(t)[1]) ** 2 for t, po in rows if 0 < t <= 1 + 1e-9]
print("%d steps, RMSE of the outlet pressure Po against the exact solution: %.4g" % (
    len(err), math.sqrt(sum(err) / len(err))))
