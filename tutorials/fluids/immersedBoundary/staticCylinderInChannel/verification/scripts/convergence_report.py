"""
Mesh convergence assessment, summary and plots of the staticCylinderInChannel
verification study, from the rows of runs.csv (see
static_cylinder_verification.py).

An observed order of convergence is only given for a triplet of consecutive
mesh levels (refinement ratio r = 2) whose three cases are steady and whose
differences converge monotonically, 0 < R = d_fine/d_coarse < 1, where
d_k = Cd_k - Cd_(k+1), and whose differences are at least 100 times the
variation of Cd over the final window of these cases. A method-level order,
Richardson extrapolation and grid convergence index are only given if both
triplets (1, 2, 3) and (2, 3, 4) satisfy these conditions, their orders
agree to within 0.3, and the finer order is between 0.5 and 3.
"""

import csv
import math
from pathlib import Path

R_REFINE = 2.0
P_AGREEMENT = 0.3
P_RANGE = (0.5, 3.0)
NOISE_FACTOR = 100.0
GCI_FS = 1.25

# Categorical palette (light mode, slots 1-4), with a marker shape for each
# method as secondary encoding
STYLE = {
    "cutLink": ("#2a78d6", "o"),
    "ghostCell": ("#eb6834", "s"),
    "penalty": ("#1baf7a", "^"),
    "incremental": ("#eda100", "D"),
}
INK = "#2b2b29"


def display(method):
    """Name of the formulation of a method in the figures and tables"""
    return "legacy incremental" if method == "incremental" else method

INK_MUTED = "#6b6a63"
GRID = "#e4e3de"
BAND = "#d9d8d2"


def _triplet(rows3):
    """Convergence assessment of three consecutive levels"""
    levels = [r["level"] for r in rows3]
    out = {
        "levels": "-".join(str(lv) for lv in levels),
        "Cd_coarse": rows3[0]["Cd_mean"], "Cd_medium": rows3[1]["Cd_mean"],
        "Cd_fine": rows3[2]["Cd_mean"],
        "d_coarse": math.nan, "d_fine": math.nan, "R": math.nan,
        "type": "", "p_obs": math.nan, "noise_ok": False, "valid": False,
        "reason": "",
    }
    if not all(r["steady"] for r in rows3):
        out["type"] = "not assessed"
        out["reason"] = "not all steady: " + ", ".join(
            f"L{r['level']}" for r in rows3 if not r["steady"])
        return out
    c1, c2, c3 = (r["Cd_mean"] for r in rows3)
    d1, d2 = c1 - c2, c2 - c3
    out["d_coarse"], out["d_fine"] = d1, d2
    noise = max(r["Cd_max"] - r["Cd_min"] for r in rows3)
    out["noise_ok"] = min(abs(d1), abs(d2)) >= NOISE_FACTOR*noise
    if d1 == 0.0:
        out["type"] = "converged (zero difference)"
        out["reason"] = "zero coarse difference"
        return out
    R = d2/d1
    out["R"] = R
    if R < 0:
        out["type"] = "oscillatory"
        out["reason"] = "differences change sign"
    elif R >= 1:
        out["type"] = "divergent"
        out["reason"] = "differences do not decrease"
    else:
        out["type"] = "monotone"
        out["p_obs"] = math.log(1.0/R)/math.log(R_REFINE)
        if not out["noise_ok"]:
            out["reason"] = "differences not >= 100x the window variation"
        else:
            out["valid"] = True
    return out


def assess(rows, methods, cd_ref):
    """Rows of convergence.csv: two triplets and a method-level verdict for
    each method"""
    result = []
    for method in methods:
        seq = sorted((r for r in rows if r["method"] == method),
                     key=lambda r: r["level"])
        levels = [r["level"] for r in seq]
        triplets = []
        for start in (1, 2):
            want = [start, start + 1, start + 2]
            if all(lv in levels for lv in want):
                triplets.append(_triplet([seq[levels.index(lv)]
                                          for lv in want]))
        verdict = {
            "accepted": False, "p": math.nan, "Cd_RE": math.nan,
            "GCI_fine": math.nan, "Cd_RE_minus_ref": math.nan, "reason": "",
        }
        if len(triplets) < 2:
            verdict["reason"] = "levels 1-4 not all available"
        else:
            t1, t2 = triplets
            if not (t1["valid"] and t2["valid"]):
                bad = [t for t in (t1, t2) if not t["valid"]]
                verdict["reason"] = "; ".join(
                    f"levels {t['levels']}: {t['type']}"
                    + (f" ({t['reason']})" if t["reason"] else "")
                    for t in bad)
            elif abs(t1["p_obs"] - t2["p_obs"]) > P_AGREEMENT:
                verdict["reason"] = (
                    f"orders of the triplets differ: {t1['p_obs']:.2f} and "
                    f"{t2['p_obs']:.2f} (not in the asymptotic range)")
            elif not P_RANGE[0] <= t2["p_obs"] <= P_RANGE[1]:
                verdict["reason"] = (f"order {t2['p_obs']:.2f} outside "
                                     f"{P_RANGE}")
            else:
                p = t2["p_obs"]
                fine = t2["Cd_fine"]
                cd_re = fine - t2["d_fine"]/(R_REFINE**p - 1.0)
                verdict.update({
                    "accepted": True, "p": p, "Cd_RE": cd_re,
                    "GCI_fine": GCI_FS*abs(t2["d_fine"]/fine)
                    / (R_REFINE**p - 1.0),
                    "Cd_RE_minus_ref": cd_re - cd_ref,
                    "reason": "monotone, consistent orders",
                })
        for t in triplets:
            row = {"method": method, "triplet": t["levels"]}
            row.update({k: v for k, v in t.items() if k != "levels"})
            row.update({f"method_{k}": v for k, v in verdict.items()})
            result.append(row)
        if not triplets:
            row = {"method": method, "triplet": ""}
            row.update({f"method_{k}": v for k, v in verdict.items()})
            result.append(row)
    return result


def _fmt(v, spec):
    if isinstance(v, float) and math.isnan(v):
        return "-"
    return format(v, spec)


def _axes_style(ax):
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(INK_MUTED)
    ax.tick_params(colors=INK_MUTED, labelcolor=INK)
    ax.grid(True, which="major", color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)


def _plot(rows, out_dir, cd_ref, cd_bounds, methods):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update({
        "font.size": 10, "axes.labelcolor": INK, "text.color": INK,
        "axes.titlesize": 11, "legend.frameon": False,
    })
    files = []

    def series(method, key_x, key_y, transform=lambda r, v: v):
        seq = sorted((r for r in rows if r["method"] == method),
                     key=lambda r: r["level"])
        return ([r[key_x] for r in seq],
                [transform(r, r[key_y]) for r in seq],
                [r["steady"] for r in seq], seq)

    def draw(ax, xs, ys, steady, method, labels, side="right"):
        colour, marker = STYLE[method]
        ax.plot(xs, ys, color=colour, linewidth=1.6, zorder=3,
                label=display(method),
                marker=marker, markersize=7, markeredgewidth=1.4)
        for x, y, s in zip(xs, ys, steady):
            if not s:
                ax.plot(x, y, marker=marker, markersize=7, linestyle="none",
                        color=colour, markerfacecolor="white",
                        markeredgewidth=1.4, zorder=4)
        labels.append((xs[-1], ys[-1], display(method), side))

    def place_labels(ax, fig, labels, gap_px=12):
        """Direct labels at the end of each series, moved apart vertically
        so that they do not overlap"""
        fig.canvas.draw()
        pts = []
        for x, y, text, side in labels:
            px, py = ax.transData.transform((x, y))
            pts.append([py, px, text, side])
        pts.sort()
        for i in range(1, len(pts)):
            if pts[i][0] - pts[i - 1][0] < gap_px:
                pts[i][0] = pts[i - 1][0] + gap_px
        inv = ax.transData.inverted()
        for py, px, text, side in pts:
            dx = 8 if side == "right" else -8
            xd, yd = inv.transform((px + dx, py))
            ax.text(xd, yd, text, va="center", fontsize=9, color=INK,
                    ha="left" if side == "right" else "right")

    note = "Hollow markers: not steady by the criterion of the study"

    # 1. Cd versus cells across D: all the methods, and the methods close to
    # the reference on a zoomed scale
    fig, (ax, az) = plt.subplots(1, 2, figsize=(10.5, 4.2))
    zoom = [m for m in methods
            if all(abs(r["Cd_mean"] - cd_ref) < 0.5 for r in rows
                   if r["method"] == m)]
    for axis, shown, title in (
        (ax, methods, "All the methods"),
        (az, zoom, "Close to the reference (zoomed)"),
    ):
        axis.axhspan(*cd_bounds, color=BAND, alpha=0.6, linewidth=0,
                     label="Schäfer and Turek (1996) bounds")
        axis.axhline(cd_ref, color=INK_MUTED, linestyle="--", linewidth=1.0,
                     label=f"reference {cd_ref:.5f}")
        labels = []
        for m in shown:
            xs, ys, st, _ = series(m, "cells_per_D", "Cd_mean")
            if xs:
                draw(axis, xs, ys, st, m, labels)
        axis.set_xscale("log", base=2)
        axis.set_xticks([10, 20, 40, 80])
        axis.set_xticklabels(["10", "20", "40", "80"])
        axis.set_xlim(8, 170)
        axis.set_xlabel("Cells across the cylinder diameter, D/h")
        axis.set_title(title, loc="left")
        _axes_style(axis)
        axis._s4f_labels = labels
    ax.set_ylabel(r"$C_d$ (mean over the final window)")
    ax.legend(loc="upper left", fontsize=8)
    fig.suptitle("Drag coefficient versus mesh resolution", x=0.01,
                 ha="left", fontsize=11)
    fig.text(0.01, 0.01, note, fontsize=8, color=INK_MUTED)
    fig.tight_layout(rect=(0, 0.03, 1, 0.97))
    for axis in (ax, az):
        place_labels(axis, fig, axis._s4f_labels)
    path = out_dir / "Cd_vs_cells_per_D.png"
    fig.savefig(path, dpi=200)
    plt.close(fig)
    files.append(path)

    # 2. |Cd - Cd_ref| versus h/D, with slope guides
    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    allx = []
    labels = []
    for m in methods:
        xs, ys, st, _ = series(m, "h_over_D", "Cd_mean",
                               lambda r, v: abs(v - cd_ref))
        if xs:
            draw(ax, xs, ys, st, m, labels, side="left")
            allx += xs
    ax.set_xscale("log")
    ax.set_yscale("log")
    if allx:
        x0, x1 = min(allx), max(allx)
        y_top = ax.get_ylim()[1]
        for order, ls in ((1, ":"), (2, "-.")):
            gx = [x0, x1]
            gy = [y_top*0.5*(x/x1)**order for x in gx]
            ax.plot(gx, gy, color=INK_MUTED, linestyle=ls, linewidth=1.0,
                    label=f"slope {order}")
        ax.set_xlim(x0/2.2, x1*1.15)
    band = max(abs(cd_bounds[0] - cd_ref), abs(cd_bounds[1] - cd_ref))
    ax.axhline(band, color=INK_MUTED, linestyle="--", linewidth=1.0)
    ax.annotate("half-width of the bounds", (ax.get_xlim()[1], band),
                xytext=(-4, 3), textcoords="offset points", fontsize=8,
                color=INK_MUTED, ha="right")
    ax.set_xlabel("h/D")
    ax.set_ylabel(r"$|C_d - C_{d,ref}|$")
    ax.set_title(r"Drag error versus cell size", loc="left")
    _axes_style(ax)
    ax.grid(True, which="minor", color=GRID, linewidth=0.4)
    ax.legend(loc="lower right", fontsize=8)
    fig.text(0.01, 0.01, note + f"; Cd_ref = {cd_ref}", fontsize=8,
             color=INK_MUTED)
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    place_labels(ax, fig, labels)
    path = out_dir / "Cd_error_vs_h_over_D.png"
    fig.savefig(path, dpi=200)
    plt.close(fig)
    files.append(path)

    # 3. |Cd - Cd_ref| versus CPU hours
    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    labels = []
    for m in methods:
        xs, ys, st, seq = series(m, "cpu_hours", "Cd_mean",
                                 lambda r, v: abs(v - cd_ref))
        if xs:
            draw(ax, xs, ys, st, m, labels)
            for x, y, r in zip(xs, ys, seq):
                ax.annotate(f"L{r['level']}", (x, y), xytext=(-4, 7),
                            textcoords="offset points", fontsize=7,
                            color=INK_MUTED, ha="right")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.axhline(band, color=INK_MUTED, linestyle="--", linewidth=1.0)
    ax.set_xlabel("CPU hours (ranks x wall time)")
    ax.set_ylabel(r"$|C_d - C_{d,ref}|$")
    ax.set_title("Drag error versus computational cost", loc="left")
    _axes_style(ax)
    ax.set_xlim(right=ax.get_xlim()[1]*4)
    ax.legend(loc="lower left", fontsize=8)
    fig.text(0.01, 0.01, note + "; dashed: half-width of the bounds",
             fontsize=8, color=INK_MUTED)
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    place_labels(ax, fig, labels)
    path = out_dir / "Cd_error_vs_cpu_hours.png"
    fig.savefig(path, dpi=200)
    plt.close(fig)
    files.append(path)
    return files


def report(rows, out_dir: Path, cd_ref, cd_bounds, cl_ref, cl_bounds,
           methods, levels):
    conv = assess(rows, methods, cd_ref)
    path = out_dir / "convergence.csv"
    keys = []
    for r in conv:
        keys += [k for k in r if k not in keys]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=keys)
        writer.writeheader()
        writer.writerows(conv)
    print(f"Wrote {path}")

    plots = _plot(rows, out_dir, cd_ref, cd_bounds, methods)
    for p in plots:
        print(f"Wrote {p}")

    lines = [
        "# staticCylinderInChannel: mesh convergence of the four methods",
        "",
        f"Schäfer and Turek (1996) 2D-1, Re = 20. Reference "
        f"`Cd_ref = {cd_ref}`; bounds `{cd_bounds[0]} <= Cd <= "
        f"{cd_bounds[1]}`, `{cl_bounds[0]} <= Cl <= {cl_bounds[1]}`.",
        f"Fixed time step {rows[0]['deltaT']} s; statistics are "
        "time-weighted over the final window "
        f"[{rows[0]['window_t0']:g}, {rows[0]['end_time']:g}] s or later.",
        "",
        "| Method | Level | Cells/D | End time (s) | Cd mean | Cd abs. error"
        " | Cd rel. error | Cl mean | Steady | Ranks | Wall time (s) |"
        " CPU hours |",
        "| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |"
        " --- |",
    ]
    for m in methods:
        for r in sorted((r for r in rows if r["method"] == m),
                        key=lambda r: r["level"]):
            lines.append(
                f"| {display(m)} | {r['level']} | {r['cells_per_D']} | "
                f"{r['end_time']:g} | {r['Cd_mean']:.5f} | "
                f"{r['Cd_abs_err']:+.5f} | {100*r['Cd_rel_err']:+.3f}% | "
                f"{r['Cl_mean']:.5f} | {'yes' if r['steady'] else 'NO'} | "
                f"{r['ranks']} | {_fmt(r['wall_s'], '.0f')} | "
                f"{_fmt(r['cpu_hours'], '.3f')} |")
    lines += [
        "",
        "## Convergence",
        "",
        "| Method | Levels | d_coarse | d_fine | R | Type | p (triplet) |"
        " Notes |",
        "| --- | --- | --- | --- | --- | --- | --- | --- |",
    ]
    for c in conv:
        if not c.get("triplet"):
            continue
        lines.append(
            f"| {display(c['method'])} | {c['triplet']} | "
            f"{_fmt(c['d_coarse'], '+.5f')}"
            f" | {_fmt(c['d_fine'], '+.5f')} | {_fmt(c['R'], '.3f')} | "
            f"{c['type']} | "
            f"{_fmt(c['p_obs'], '.2f') if c['valid'] else '-'} | "
            f"{c['reason']} |")
    lines += ["", "| Method | Observed order | Richardson Cd | GCI (fine) |"
              " Cd_RE - Cd_ref | Verdict |",
              "| --- | --- | --- | --- | --- | --- |"]
    seen = set()
    for c in conv:
        if c["method"] in seen:
            continue
        seen.add(c["method"])
        ok = c["method_accepted"]
        lines.append(
            f"| {display(c['method'])} | "
            f"{_fmt(c['method_p'], '.2f') if ok else '-'} | "
            f"{_fmt(c['method_Cd_RE'], '.5f') if ok else '-'} | "
            f"{_fmt(100*c['method_GCI_fine'], '.3f') + '%' if ok else '-'} | "
            f"{_fmt(c['method_Cd_RE_minus_ref'], '+.5f') if ok else '-'} | "
            f"{'accepted' if ok else 'not reported'}: {c['method_reason']} |")
    lines += [
        "",
        "An order is reported for a triplet only if its three cases are "
        "steady, its differences decrease monotonically (0 < R < 1) and "
        "they are at least 100 times the variation of Cd over the final "
        "window. A method-level order, Richardson extrapolation and GCI "
        "(Fs = 1.25) are reported only if both triplets satisfy this and "
        "their orders agree to within 0.3.",
        "",
        "Plots: " + ", ".join(f"`{p.name}`" for p in plots),
        "",
    ]
    path = out_dir / "summary.md"
    path.write_text("\n".join(lines))
    print(f"Wrote {path}")
