#!/usr/bin/env python3
"""Plot the healthy-human loop after running circulationFoam.

    circulationFoam
    python3 plotter.py                  # the last cycle: figures open on screen and are saved to plots/
    python3 plotter.py --cycle 12       # cycle 12 instead (the one ending at t = 12 T)
    python3 plotter.py --no-show        # only save the files (no windows, e.g. on a remote machine)
    python3 plotter.py other.dat        # another results file of the same circuit

The windows stay open until you close them; the files are written first.

Writes, into ./plots/:

    pv_loops.png        pressure-volume loops of the four chambers
    valve_pressures.png the pressure each side of the four valves, with the valve flow
    wiggers.png         left heart: pressures, left-ventricular volume, aortic flow

The results file is postProcessing/circulation/0/circulation.dat; the cycle is
one period Tcycle (read from constant/circulationProperties, 1 s if absent).
Needs matplotlib:  pip3 install matplotlib
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

NAMES = {'SAS.C.P': 'aorta (SAS)', 'PAS.C.P': 'pulmonary artery (PAS)'}


def read_results(path):
    """Columns by name (lists of floats) and their units, from a circulationFoam file."""
    names, units, rows = None, None, []
    with open(path) as f:
        for line in f:
            if line.startswith('#'):
                words = line[1:].split()
                if words and words[0] == 't' and names is None:
                    names = words
                elif words and words[0] == 's' and names is not None and units is None:
                    units = words
                continue
            if line.strip():
                rows.append([float(x) for x in line.split()])
    if names is None or not rows:
        sys.exit("%s is not a circulationFoam results file" % path)
    cols = {n: [r[i] for r in rows] for i, n in enumerate(names)}
    return cols, dict(zip(names, units or [''] * len(names)))


def read_period(case):
    """Tcycle from constant/circulationProperties (seconds), default 1."""
    try:
        text = open(os.path.join(case, 'constant', 'circulationProperties')).read()
        m = re.search(r'^\s*Tcycle\s+([0-9.eE+-]+)\s*;', text, re.M)
        return float(m.group(1)) if m else 1.0
    except OSError:
        return 1.0


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('results', nargs='?',
                    default=os.path.join(HERE, 'postProcessing', 'circulation', '0', 'circulation.dat'))
    ap.add_argument('--cycle', type=int, help='cycle number to plot (default: the last one)')
    ap.add_argument('--out', default=os.path.join(HERE, 'plots'), help='output directory')
    ap.add_argument('--no-show', action='store_true', help='do not open windows, only save the files')
    args = ap.parse_args()

    if args.no_show:
        matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    if not os.path.exists(args.results):
        sys.exit("%s not found: run circulationFoam first" % args.results)

    d, u = read_results(args.results)
    T = read_period(HERE)
    t = d['t']
    last = int((t[-1] + 1e-9) // T)
    cycle = args.cycle or last
    t0, t1 = (cycle - 1)*T, cycle*T
    idx = [i for i, x in enumerate(t) if t0 - 1e-9 <= x <= t1 + 1e-9]
    if len(idx) < 2:
        sys.exit("no data for cycle %d (the run covers t = %g to %g s)" % (cycle, t[0], t[-1]))

    def get(name):
        return [d[name][i] for i in idx]

    tt = [x - t0 for x in get('t')]
    unit = lambda c: u.get(c, '')
    os.makedirs(args.out, exist_ok=True)
    written = []

    # ---- 1. Pressure-volume loops ------------------------------------------------
    fig, axs = plt.subplots(2, 2, figsize=(9, 7))
    for ax, ch, title in zip(axs.flat, ['LV', 'RV', 'LA', 'RA'],
                             ['Left ventricle', 'Right ventricle', 'Left atrium', 'Right atrium']):
        ax.plot(get(ch + '.V'), get(ch + '.P'), lw=1.8)
        ax.set_title(title)
        ax.set_xlabel('volume (%s)' % unit(ch + '.V'))
        ax.set_ylabel('pressure (%s)' % unit(ch + '.P'))
        ax.grid(alpha=0.3)
    fig.suptitle('Pressure-volume loops, cycle %d' % cycle)
    fig.tight_layout()
    written.append(os.path.join(args.out, 'pv_loops.png'))
    fig.savefig(written[-1], dpi=150)

    # ---- 2. Pressures either side of each valve --------------------------------
    # (upstream pressure, downstream pressure, valve flow, title)
    valves = [
        ('LV.P', 'SAS.C.P', 'AV.Q', 'Aortic valve: LV | aorta'),
        ('LA.P', 'LV.P', 'MV.Q', 'Mitral valve: LA | LV'),
        ('RV.P', 'PAS.C.P', 'PV.Q', 'Pulmonary valve: RV | pulmonary artery'),
        ('RA.P', 'RV.P', 'TV.Q', 'Tricuspid valve: RA | RV'),
    ]
    fig, axs = plt.subplots(2, 2, figsize=(11, 7), sharex=True)
    for ax, (up, down, q, title) in zip(axs.flat, valves):
        ax.plot(tt, get(up), label=NAMES.get(up, up) + ' (upstream)', lw=1.6)
        ax.plot(tt, get(down), label=NAMES.get(down, down) + ' (downstream)', lw=1.6)
        ax.set_ylabel('pressure (%s)' % unit(up))
        ax.set_title(title)
        ax.grid(alpha=0.3)
        ax2 = ax.twinx()
        ax2.plot(tt, get(q), color='0.45', lw=1.0, ls='--', label=q)
        ax2.set_ylabel('valve flow (%s)' % unit(q))
        h1, l1 = ax.get_legend_handles_labels()
        h2, l2 = ax2.get_legend_handles_labels()
        ax.legend(h1 + h2, l1 + l2, loc='upper right', fontsize=8)
    for ax in axs[1]:
        ax.set_xlabel('time in cycle (s)')
    fig.suptitle('Pressure each side of the valves, cycle %d' % cycle)
    fig.tight_layout()
    written.append(os.path.join(args.out, 'valve_pressures.png'))
    fig.savefig(written[-1], dpi=150)

    # ---- 3. Wiggers-style diagram of the left heart ---------------------------
    fig, axs = plt.subplots(3, 1, figsize=(8, 8), sharex=True)
    axs[0].plot(tt, get('LV.P'), label='LV')
    axs[0].plot(tt, get('SAS.C.P'), label='aorta')
    axs[0].plot(tt, get('LA.P'), label='LA')
    axs[0].set_ylabel('pressure (%s)' % unit('LV.P'))
    axs[0].legend(loc='upper right', fontsize=8)
    axs[1].plot(tt, get('LV.V'), color='tab:green')
    axs[1].set_ylabel('LV volume (%s)' % unit('LV.V'))
    axs[2].plot(tt, get('AV.Q'), label='aortic valve', color='tab:red')
    axs[2].plot(tt, get('MV.Q'), label='mitral valve', color='tab:purple')
    axs[2].set_ylabel('flow (%s)' % unit('AV.Q'))
    axs[2].set_xlabel('time in cycle (s)')
    axs[2].legend(loc='upper right', fontsize=8)
    for ax in axs:
        ax.grid(alpha=0.3)
    fig.suptitle('Left heart, cycle %d' % cycle)
    fig.tight_layout()
    written.append(os.path.join(args.out, 'wiggers.png'))
    fig.savefig(written[-1], dpi=150)

    # ---- numbers worth a glance ---------------------------------------------------
    lv_v, lv_p = get('LV.V'), get('LV.P')
    edv, esv = max(lv_v), min(lv_v)
    print("cycle %d (t = %g to %g s): LV EDV %.1f %s, ESV %.1f %s, SV %.1f %s, EF %.0f %%, peak LV pressure %.0f %s"
          % (cycle, t0, t1, edv, unit('LV.V'), esv, unit('LV.V'), edv - esv, unit('LV.V'),
             100*(edv - esv)/edv, max(lv_p), unit('LV.P')))
    for w in written:
        print('wrote', os.path.relpath(w))

    if not args.no_show:
        plt.show()


if __name__ == '__main__':
    main()
