#!/usr/bin/env python3
"""Generator of reference/results/controls/velocity-regions.csv: the velocity
L2 error against the exact initial-value solution (fraction of U0) in the
interface row, within one Stokes thickness, in the bulk (s > delta) and in all
fluid rows, and where the largest error is (in cell heights), for the cutLink
momentum-exchange run and the other methods at the tutorial step and at
1.25e-4 s, over the eight snapshots of the fourth period.

    velocity_regions.py --work-dir D [D ...] [--out FILE]

The written fields of the completed cases of the work directories are the
input (the runs 'matrix' and 'dt-fine'); the default output is
verification/work/export/velocity-regions.csv. Standard library only.
"""

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import oscillating_wall_stokes_layer_verification as v  # noqa: E402
import export_results as ex  # noqa: E402

HEADER = ("method,level,dt,rows_all,rows_near,rows_bulk,L2_all,"
          "L2_interface_row,L2_near,L2_bulk,s_of_largest_error_over_h")
STEPS = (v.TUTORIAL_DT, 0.000125)
METHODS = (("cutLink", "ME", "cutLink"), ("penalty", "", "penalty"),
           ("incremental", "", "incremental"), ("ghostCell", "", "ghostCell"))


def region_rows(p: dict, level: int, case: Path, label: str, dt: float,
                modes: list) -> str:
    H = p["height"] - p["yw"]
    h, s = v.cell_heights(p, level)
    sets = {"all": [j for j, x in enumerate(s) if x > 0],
            "iface": [j for j, x in enumerate(s) if 0 < x < h],
            "near": [j for j, x in enumerate(s) if 0 < x <= p["delta"]],
            "bulk": [j for j, x in enumerate(s) if x > p["delta"]]}
    acc = {k: [] for k in sets}
    worst = (0.0, 0.0)
    for tk in v.snapshot_times():
        ux = v.row_profile(p, level, v.read_vector_field(case / f"{tk:g}" / "U"))
        err = [abs(ux[j] - v.exact_u_ivp(p, s[j], tk, H, modes))
               if s[j] > 0 else 0.0 for j in range(len(s))]
        for k, rows in sets.items():
            acc[k].append(v.rms([err[j] for j in rows]) / p["U0"])
        for j in sets["all"]:
            if err[j] > worst[0]:
                worst = (err[j], s[j])
    f = {k: v.rms(acc[k]) for k in sets}
    return (f"{label},{level},{dt:g},{len(sets['all'])},{len(sets['near'])},"
            f"{len(sets['bulk'])},{f['all']:.4e},{f['iface']:.4e},"
            f"{f['near']:.4e},{f['bulk']:.4e},{worst[1] / h:.2f}")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--work-dir", nargs="+", required=True)
    ap.add_argument("--out", default="")
    a = ap.parse_args(argv)
    p = v.parameters()
    modes = v.cached_modes(p)
    index = {}
    for name, e, case in ex.completed(a.work_dir):
        if ex.is_baseline(e):
            index[(e["method"], e["estimator"], e["level"],
                   e["dt"] or v.TUTORIAL_DT)] = case
    out = [HEADER]
    for method, est, label in METHODS:
        for dt in STEPS:
            for level in v.LEVELS:
                case = index.get((method, est, level, dt))
                if case is not None:
                    out.append(region_rows(p, level, case, label, dt, modes))
    path = Path(a.out) if a.out else v.WORK_DIR / "export/velocity-regions.csv"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(out) + "\n")
    print("\n".join(out))
    print(f"written {path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
