#!/usr/bin/env python3
"""
Write the compact results of the completed static cylinder studies to
verification/reference/ for version control, from the cases in work/ (see
results.md). One row per case, with its settings, the statistics of the
forces over the final 2 s, its status (steady, limit cycle or drifting) and
its cost. The cases of method "incremental" (the solver identifier) are
labelled "common-settings incremental baseline" with nOuterCorrectors 1 (the
earlier studies), and "legacy incremental" otherwise (the final campaign:
nOuterCorrectors 3-7, couplingCoeff and time step studies, and the frozen
configuration couplingCoeff 0.2, nOuterCorrectors 5).

    mesh_study.csv                    common settings, nOuterCorrectors 1
                                      for all four formulations (legacy
                                      incremental with couplingCoeff 0.8)
    incremental_sensitivity.csv       legacy incremental, nOuterCorrectors 1,
                                      couplingCoeff and time step studies
    outer_corrector_study.csv         nOuterCorrectors 1-7 (final-A)
    historical_parity.csv             historical pimpleHFDIBFoam cases with
                                      the current port (final-B)
    legacy_incremental_coupling.csv   couplingCoeff 0.05-0.8, nOuter 5 (D)
    legacy_incremental_timestep.csv   deltaT, couplingCoeff 0.2, nOuter 5 (E)
    final_method_comparison.csv       frozen settings, levels 1-4 (F)
    mpi_scaling.csv                   rank calibrations (calibration, C)

Run after the cases have completed; work/ is not needed afterwards.
"""

import csv
import glob
import json
import math
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import static_cylinder_verification as scv  # noqa: E402

REFERENCE = scv.VERIFICATION_DIR / "reference"
SOURCE = ("solids4foam f33b244a2, OpenFOAM v2412; Cd and Cl are time-weighted "
          "means over the final 2 s; wall_s and cpu_hours (ranks x wall time "
          "of the job script) are machine-dependent (2 x AMD EPYC 9684X)")

# Saved drag coefficients of the historical openHFDIB-DEM pimpleHFDIBFoam
# cases (mean over 19 < t < 20 s, Cd_comparison.dat of the manuscript-era
# Immersed_Cyl_Static study)
HISTORICAL_CD = {
    "ImmersedCyl-Static-Re20-Mesh20000": 5.72124,
    "ImmersedCyl-Static-Re20-Mesh80000": 5.70618,
}

LEGACY = "legacy incremental"
BASELINE = "common-settings incremental baseline"

COLUMNS = [
    ("study", None), ("case", None), ("method", None), ("formulation", None),
    ("level", None),
    ("cells_per_D", None), ("n_cells", None), ("ranks", None),
    ("couplingCoeff", None), ("occupancy", None), ("deltaT", None),
    ("nOuterCorrectors", None), ("nCorrectors", None), ("end_time", "g"),
    ("max_co", ".3f"), ("max_co_developed", ".3f"),
    ("Cd_mean", ".6f"), ("Cd_range", ".3e"), ("Cd_rel_range", ".2e"),
    ("Cd_rel_drift_10s", ".2e"), ("Cl_mean", ".6f"), ("Cl_range", ".2e"),
    ("Cd_range_ratio", ".3f"), ("Cl_range_ratio", ".3f"),
    ("Cl_step_alternation", ".2f"), ("status", None),
    ("Cd_surface_traction_last", ".5f"), ("p_iters_per_step", ".1f"),
    ("n_steps", None), ("wall_s", ".0f"), ("cpu_hours", ".3f"),
]


def formulation(method, n_outer):
    if method != "incremental":
        return method
    return BASELINE if n_outer == 1 else LEGACY


def rows_of(study):
    if not (scv.WORK_DIR / study).exists():
        return []
    rows = scv.case_rows(study, 2.0, 1e-5, 1e-5)
    for r in rows:
        r["nOuterCorrectors"] = r["nOuterCorrectors"] or 1
        r["formulation"] = formulation(r["method"], r["nOuterCorrectors"])
        r["Cd_range"] = r["Cd_max"] - r["Cd_min"]
        meta = json.loads((scv.WORK_DIR / study / r["case"]
                           / "case.json").read_text())
        r["occupancy"] = meta.get("body", {}).get("occupancy", "")
    return [r for r in rows if r["complete"]]


def compact(row, columns):
    out = {}
    for name, spec in columns:
        value = row.get(name, "")
        if isinstance(value, float) and math.isnan(value):
            value = ""
        elif spec and isinstance(value, float):
            value = format(value, spec)
        out[name] = value
    return out


def write(name, rows, columns=COLUMNS, note=SOURCE):
    path = REFERENCE / name
    with path.open("w", newline="") as handle:
        handle.write(f"# staticCylinderInChannel verification: {note}\n")
        writer = csv.DictWriter(handle, fieldnames=[c for c, _ in columns],
                                lineterminator="\n")
        writer.writeheader()
        writer.writerows(compact(r, columns) for r in rows)
    print(f"Wrote {path} ({len(rows)} rows)")


def find(rows, **kw):
    out = []
    for r in rows:
        ok = True
        for key, value in kw.items():
            have = r.get(key)
            if isinstance(value, float):
                ok = ok and have not in ("", None) and \
                    abs(float(have) - value) < 1e-9
            else:
                ok = ok and have == value
        if ok:
            out.append(r)
    return out


def one(rows, **kw):
    found = find(rows, **kw)
    if len(found) != 1:
        sys.exit(f"ERROR: {len(found)} cases match {kw}")
    return found[0]


def historical_parity():
    rows = []
    for case in sorted((scv.WORK_DIR / "final-B").glob("*")):
        if not (case / "case.json").exists():
            continue
        meta = json.loads((case / "case.json").read_text())
        files = scv.force_files(case)
        t, fx, fy, _ = scv.read_forces(files)
        scale = 1.0/(0.5*1000.0*1.0**2*1.0*0.1)
        cd = [scale*f for f in fx]
        cl = [scale*f for f in fy]
        # The historical value is the mean over the final 1 s
        scd = scv.window_stats(t, cd, t[-1] - 1.0)
        scl = scv.window_stats(t, cl, t[-1] - 1.0)
        old = Path(meta["historical_case"]).name
        wall = sum(int(w) for w in re.findall(
            r"wall_s=(\d+)", (case / "wallclock.txt").read_text()))
        stats = scv.log_stats(sorted(case.glob("log.solids4Foam*")))
        rows.append({
            "historical_case": old,
            "n_cells": re.search(r"nCells:\s*(\d+)",
                                 (case / "log.blockMesh").read_text()).group(1),
            "inlet": meta["inlet"], "formulation": LEGACY,
            "couplingCoeff": meta["couplingCoeff"],
            "deltaT": meta["deltaT"], "end_time": t[-1],
            "nOuterCorrectors": 5, "nCorrectors": 3,
            "nNonOrthogonalCorrectors": 2,
            "Cd_historical": HISTORICAL_CD[old],
            "Cd_current": scd[0],
            "Cd_current_minus_historical": scd[0] - HISTORICAL_CD[old],
            "Cl_current": scl[0], "max_co": stats["max_co"],
            "ranks": meta["ranks"], "wall_s": wall,
            "cpu_hours": wall*meta["ranks"]/3600.0,
        })
    columns = [(k, None) for k in rows[0]]
    columns = [(k, ".6f" if k.startswith(("Cd", "Cl")) else
                (".3f" if k in ("max_co", "cpu_hours") else None))
               for k, _ in columns]
    columns = [(k, ".2e" if k == "Cd_current_minus_historical" else s)
               for k, s in columns]
    write("historical_parity.csv", rows, columns,
          "historical pimpleHFDIBFoam cases (scaled by 10: D = 1, rho = "
          "1000, Ubar = 1) with the current legacy incremental port and the "
          "historical numerical settings; Cd = Fx/(0.5 rho Ubar^2 D Lz), mean "
          "over the final 1 s; Cd_historical as saved by pimpleHFDIBFoam")


def mpi_scaling():
    rows = []
    for name, path, case in (
        ("calibration", "work/calibration/scaling/scaling.csv",
         "cutLink L4, nOuter 1, from t = 3 s"),
        ("final-C", "work/final-C/scaling/scaling.csv",
         "legacy incremental L4, c 0.2, nOuter 5, from t = 15 s"),
    ):
        f = scv.VERIFICATION_DIR / path
        if not f.exists():
            continue
        for r in csv.DictReader(f.open()):
            ranks = int(r["ranks"])
            wall = float(r["solver_wall_s"])
            run = f.parent / f"ranks-{ranks}"
            t, fx, fy, _ = scv.read_forces(scv.force_files(run))
            skip = int(r["nsteps"])//4
            cd = [scv.COEFF_SCALE*x for x in fx[skip:]]
            cl = [scv.COEFF_SCALE*x for x in fy[skip:]]
            rows.append({
                "study": name, "case": case, "ranks": ranks,
                "time_steps": r["nsteps"], "wall_s": wall,
                "wall_s_per_step": wall/int(r["nsteps"]),
                "cpu_hours": ranks*wall/3600.0,
                "Cd_mean_last_75pct": sum(cd)/len(cd),
                "Cd_range_last_75pct": max(cd) - min(cd),
                "Cl_mean_last_75pct": sum(cl)/len(cl),
                "Cl_range_last_75pct": max(cl) - min(cl),
            })
    columns = [("study", None), ("case", None), ("ranks", None),
               ("time_steps", None), ("wall_s", ".1f"),
               ("wall_s_per_step", ".4f"), ("cpu_hours", ".4f"),
               ("Cd_mean_last_75pct", ".7f"), ("Cd_range_last_75pct", ".2e"),
               ("Cl_mean_last_75pct", ".7f"), ("Cl_range_last_75pct", ".2e")]
    write("mpi_scaling.csv", rows, columns,
          "the same time steps from the same developed state with each rank "
          "count (scripts/mpi_scaling.sh), mpirun --bind-to none")


def main():
    REFERENCE.mkdir(exist_ok=True)
    mesh = rows_of("mesh")
    inc1 = rows_of("incA") + rows_of("incB") + rows_of("incC")
    final_a = rows_of("final-A")
    final_d = rows_of("final-D")
    final_e = rows_of("final-E")
    final_f = rows_of("final-F")

    write("mesh_study.csv", mesh)
    write("incremental_sensitivity.csv", inc1)

    # Outer correctors: nOuter 1 from the common-settings studies
    outer = []
    for method in ("cutLink", "ghostCell", "penalty"):
        for level in (2, 4):
            outer += find(mesh, method=method, level=level)
    for level in (2, 4):
        outer += find(inc1, level=level, couplingCoeff=0.2, deltaT=0.001)
    outer += final_a
    outer.sort(key=lambda r: (r["method"], r["level"],
                              r["nOuterCorrectors"]))
    write("outer_corrector_study.csv", outer)

    if (scv.WORK_DIR / "final-B").exists():
        historical_parity()

    # Legacy incremental with nOuter 5: couplingCoeff and time step
    coupling = [r for r in final_d]
    coupling += find(final_a, method="incremental", couplingCoeff=0.2,
                     nOuterCorrectors=5)
    coupling.sort(key=lambda r: (r["level"], float(r["couplingCoeff"])))
    write("legacy_incremental_coupling.csv", coupling)

    timestep = [r for r in final_e]
    timestep += find(final_a, method="incremental", couplingCoeff=0.2,
                     nOuterCorrectors=5)
    timestep.sort(key=lambda r: (r["level"], -float(r["deltaT"])))
    write("legacy_incremental_timestep.csv", timestep)

    # Final comparison: frozen settings
    comparison = []
    for method in ("cutLink", "ghostCell", "penalty"):
        comparison += sorted(find(mesh, method=method),
                             key=lambda r: r["level"])
    legacy = find(final_f, couplingCoeff=0.2, nOuterCorrectors=5,
                  deltaT=0.001)
    legacy += find(final_a, method="incremental", couplingCoeff=0.2,
                   nOuterCorrectors=5)
    comparison += sorted(legacy, key=lambda r: r["level"])
    for r in comparison:
        expected = 5 if r["method"] == "incremental" else 1
        if r["nOuterCorrectors"] != expected or abs(r["deltaT"] - 0.001) > 1e-12:
            sys.exit(f"ERROR: {r['case']} does not have the frozen settings")
    write("final_method_comparison.csv", comparison)

    mpi_scaling()


if __name__ == "__main__":
    main()
