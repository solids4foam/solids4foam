#!/usr/bin/env python3
"""Plot the coupled pipe closed loop after ./Allrun.

    python3 plotter.py                 # figures open on screen and are saved to plots/
    python3 plotter.py --no-show       # only save the files
    python3 plotter.py --cycle 3       # cycle 3 of a longer run (default: the last complete cycle,
                                       #  or everything, if the run is shorter than one cycle)

Writes, into ./plots/:

    interface.png   what the pipe and the circuit exchange: the flow through the outlet
                    and inlet patches and the pressure the circuit imposes there
    heart.png       the ventricle: pressure-volume loop, and pressures either side of the
                    aortic and mitral valves with the valve flows
    circuit.png     pressure in every compliance of the loop

Reads postProcessing/circulationCoupling/0/coupling.dat (SI: Pa, m3, m3/s) and shows mmHg,
mL and mL/s. Needs matplotlib:  pip3 install matplotlib
"""
import argparse
import os
import re
import sys

try:
    import matplotlib
except ImportError:
    sys.exit("plotter.py needs matplotlib:  pip3 install matplotlib")

HERE = os.path.dirname(os.path.abspath(__file__))
PA_TO_MMHG = 1.0/133.322387415
M3_TO_ML = 1.0e6


def read_coupling(path):
    names, rows = None, []
    with open(path) as f:
        for line in f:
            if line.startswith('#'):
                words = line[1:].split()
                if words and words[0] == 'time':
                    names = words
                continue
            if line.strip():
                rows.append([float(x) for x in line.split()])
    if names is None or not rows:
        sys.exit("%s is not a circulationCoupler results file (is the run finished?)" % path)
    return {n: [r[i] for r in rows] for i, n in enumerate(names)}


def read_period():
    try:
        text = open(os.path.join(HERE, 'constant', 'circulationProperties')).read()
        m = re.search(r'^\s*Tcycle\s+([0-9.eE+-]+)\s*;', text, re.M)
        return float(m.group(1)) if m else 0.8
    except OSError:
        return 0.8


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('results', nargs='?', default=os.path.join(
        HERE, 'postProcessing', 'circulationCoupling', '0', 'coupling.dat'))
    ap.add_argument('--cycle', type=int, help='cycle number to plot')
    ap.add_argument('--no-show', action='store_true', help='do not open windows, only save the files')
    args = ap.parse_args()

    if args.no_show:
        matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    if not os.path.exists(args.results):
        sys.exit("%s not found: run ./Allrun first" % args.results)

    d = read_coupling(args.results)
    T = read_period()
    t = d['time']
    complete = int((t[-1] + 1e-9)//T)
    cycle = args.cycle or complete
    if cycle >= 1:
        t0, t1 = (cycle - 1)*T, cycle*T
        title = 'cycle %d' % cycle
    else:
        t0, t1 = 0.0, t[-1]
        title = 't = 0 to %.3g s' % t[-1]
    idx = [i for i, x in enumerate(t) if t0 - 1e-9 <= x <= t1 + 1e-9]
    if len(idx) < 2:
        sys.exit("no data for %s (the run covers t = %g to %g s)" % (title, t[0], t[-1]))

    def get(name, scale=1.0):
        return [d[name][i]*scale for i in idx]

    P = lambda n: get(n, PA_TO_MMHG)
    Q = lambda n: get(n, M3_TO_ML)
    V = lambda n: get(n, M3_TO_ML)
    tt = get('time')
    out = os.path.join(HERE, 'plots')
    os.makedirs(out, exist_ok=True)
    written = []

    # ---- 1. The interface: flows out of / into the pipe, and the pressures imposed ---
    fig, axs = plt.subplots(2, 1, figsize=(8, 7), sharex=True)
    axs[0].plot(tt, Q('outlet.Q'), label='outlet (flow leaves the pipe)')
    axs[0].plot(tt, Q('inlet.Q'), '--', label='inlet (flow enters the pipe)')
    axs[0].set_ylabel('flow through the patch (mL/s)')
    axs[0].legend(fontsize=8)
    axs[1].plot(tt, P('outlet.Ppatch'), label='outlet (node C1)')
    axs[1].plot(tt, P('inlet.Ppatch'), label='inlet (node C4)')
    axs[1].set_ylabel('pressure the circuit imposes (mmHg)')
    axs[1].set_xlabel('time (s)')
    axs[1].legend(fontsize=8)
    for ax in axs:
        ax.grid(alpha=0.3)
    fig.suptitle('Pipe-circuit interface, ' + title)
    fig.tight_layout()
    written.append(os.path.join(out, 'interface.png'))
    fig.savefig(written[-1], dpi=150)

    # ---- 2. The heart -----------------------------------------------------------------
    fig = plt.figure(figsize=(11, 7))
    ax = fig.add_subplot(2, 2, 1)
    ax.plot(V('ventricle.V'), P('ventricle.P'), lw=1.8)
    ax.set_xlabel('ventricle volume (mL)')
    ax.set_ylabel('ventricle pressure (mmHg)')
    ax.set_title('Ventricle pressure-volume loop')
    ax.grid(alpha=0.3)

    ax = fig.add_subplot(2, 2, 2)
    ax.plot(tt, P('atrium.P'), label='atrium')
    ax.plot(tt, P('ventricle.P'), label='ventricle')
    ax.plot(tt, P('C4.P'), label='aortic side (C4)')
    ax.set_ylabel('pressure (mmHg)')
    ax.set_xlabel('time (s)')
    ax.set_title('Chamber pressures')
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)

    for k, (up, down, q, name) in enumerate([
            ('ventricle.P', 'C4.P', 'aortic.Q', 'Aortic valve: ventricle | C4'),
            ('atrium.P', 'ventricle.P', 'mitral.Q', 'Mitral valve: atrium | ventricle')]):
        ax = fig.add_subplot(2, 2, 3 + k)
        ax.plot(tt, P(up), label=up.split('.')[0] + ' (upstream)')
        ax.plot(tt, P(down), label=down.split('.')[0] + ' (downstream)')
        ax.set_ylabel('pressure (mmHg)')
        ax.set_xlabel('time (s)')
        ax.set_title(name)
        ax.grid(alpha=0.3)
        ax2 = ax.twinx()
        ax2.plot(tt, Q(q), color='0.45', lw=1.0, ls='--', label=q)
        ax2.set_ylabel('valve flow (mL/s)')
        h1, l1 = ax.get_legend_handles_labels()
        h2, l2 = ax2.get_legend_handles_labels()
        ax.legend(h1 + h2, l1 + l2, fontsize=8, loc='upper right')
    fig.suptitle('Heart, ' + title)
    fig.tight_layout()
    written.append(os.path.join(out, 'heart.png'))
    fig.savefig(written[-1], dpi=150)

    # ---- 3. Pressure in every compliance -----------------------------------------------
    fig, ax = plt.subplots(figsize=(8, 5))
    for n, label in [('C1.P', 'C1 (pipe outlet)'), ('C2.P', 'C2'), ('C3.P', 'C3'), ('C4.P', 'C4 (pipe inlet)')]:
        ax.plot(tt, P(n), label=label)
    ax.set_xlabel('time (s)')
    ax.set_ylabel('pressure (mmHg)')
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    fig.suptitle('Compliance pressures, ' + title)
    fig.tight_layout()
    written.append(os.path.join(out, 'circuit.png'))
    fig.savefig(written[-1], dpi=150)

    vv = V('ventricle.V')
    print("%s: ventricle EDV %.1f mL, ESV %.1f mL, SV %.1f mL, peak ventricle pressure %.0f mmHg"
          % (title, max(vv), min(vv), max(vv) - min(vv), max(P('ventricle.P'))))
    for w in written:
        print('wrote', os.path.relpath(w, HERE))

    if not args.no_show:
        plt.show()


if __name__ == '__main__':
    main()
