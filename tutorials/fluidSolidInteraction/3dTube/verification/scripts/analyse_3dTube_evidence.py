#!/usr/bin/env python3
"""Collect the 3dTube verification evidence into compact CSV and JSON files.

Run after Allverify has evaluated the studies below; it reads the CSV files
and point-A histories that Allverify writes to verification/postProcessing/
and writes verification/results/3dTube_evidence.{json,csv}:

- the mesh study (levels 1-3, time step halved with the mesh): values,
  successive changes, the observed order where it is defined and reliable,
  a Richardson estimate where that order is justified, and an uncertainty
  band otherwise;
- the separation of the spatial and temporal error: the backward time-step
  study on the level-1 mesh, and levels 1 and 2 at a fixed time step;
- the implicit-Euler comparison: Euler at dt = 1e-4 s on levels 1 and 2 and
  the Euler time-step study on level 1, and the fraction of the published
  finite element/finite volume difference in u_r,max(A) that it reproduces;
- the Robin-Neumann/IQN-ILS comparison, used as the iterative-error floor.

Every order is classified rather than forced: an order is reported only for a
monotone sequence whose finest change is above the iterative-error and
extraction floor and whose value lies in [0.5, 4] (nominal 2).
"""

from __future__ import annotations

import csv
import json
import math
import sys
from pathlib import Path

VERIFICATION = Path(__file__).resolve().parents[1]
POST = VERIFICATION / "postProcessing"
RESULTS = VERIFICATION / "results"
REFERENCE = json.loads(
    (VERIFICATION / "reference/3dTube_verification_references.json").read_text()
)

RATIO = 2.0
NOMINAL_ORDER = 2.0
ORDER_RANGE = (0.5, 4.0)
# Richardson extrapolation only within one of the nominal order
RICHARDSON_RANGE = (1.0, 3.0)

# Quantity, unit scale, unit, extraction floor (relative). The extraction
# floor is the resolution of the extraction itself: the parabolic peak fit
# and the linear interpolation of the arrival are accurate to well below
# 0.01% of the value for the displacements; the times of the flat radial peak
# move by tens of microseconds for a 0.01% change in the curve, so their
# floor is set at 0.5%.
QUANTITIES = [
    ("ur_max_A_m", 1e3, "mm", 1e-4),
    ("t_ur_max_A_s", 1e3, "ms", 5e-3),
    ("uz_min_A_m", 1e3, "mm", 1e-4),
    ("t_uz_min_A_s", 1e3, "ms", 1e-3),
    ("ur_min_late_A_m", 1e3, "mm", 1e-4),
    ("t_arrival_A_s", 1e3, "ms", 1e-4),
    ("wave_speed_pressure_m_s", 1.0, "m/s", 1e-4),
    ("wave_speed_wall_m_s", 1.0, "m/s", 1e-4),
]


def read_csv(name: str) -> list[dict]:
    path = POST / name
    if not path.is_file():
        return []
    with path.open() as handle:
        rows = list(csv.DictReader(handle))
    for row in rows:
        for key, value in row.items():
            try:
                row[key] = float(value)
            except (TypeError, ValueError):
                pass
    return rows


def read_history(name: str) -> dict[float, tuple[float, float, float]]:
    path = POST / "histories" / f"{name}.csv"
    history = {}
    with path.open() as handle:
        for row in csv.DictReader(handle):
            history[round(float(row["time_s"]), 9)] = (
                float(row["radial_displacement_m"]),
                float(row["axial_displacement_m"]),
                float(row["axis_pressure_Pa"]),
            )
    return history


def history_difference(coarse: str, fine: str) -> dict[str, float]:
    """Largest difference over the common time levels, as a fraction of the
    largest |value| of the finer history (radial, axial, pressure)."""
    a, b = read_history(coarse), read_history(fine)
    common = sorted(set(a) & set(b))
    result = {}
    for index, name in enumerate(("radial", "axial", "pressure")):
        scale = max(abs(value[index]) for value in b.values())
        result[f"{name}_history_max_difference"] = max(
            abs(a[t][index] - b[t][index]) for t in common
        ) / scale
    result["common_time_levels"] = len(common)
    return result


def classify(values: list[float], floor: float) -> dict:
    """Successive changes, observed order and Richardson estimate of a
    three-member sequence with ratio two, with an explicit verdict."""
    f1, f2, f3 = values
    e21, e32 = f2 - f1, f3 - f2
    scale = abs(f3)
    out = {
        "change_1_2": e21, "change_2_3": e32,
        "relative_change_1_2": abs(e21) / scale if scale else math.nan,
        "relative_change_2_3": abs(e32) / scale if scale else math.nan,
        "observed_order": None, "richardson_estimate": None,
        "richardson_relative_error_fine": None,
    }
    resolved = abs(e32) > floor * scale and abs(e21) > floor * scale
    if not resolved:
        verdict = ("order undefined: a change is below the iterative/"
                   "extraction floor")
    elif e21 * e32 < 0.0:
        verdict = "order undefined: non-monotone (the changes change sign)"
    else:
        order = math.log(abs(e21 / e32), RATIO)
        out["observed_order"] = order
        if abs(e32) >= abs(e21):
            verdict = "order undefined: the changes do not decrease"
        elif not ORDER_RANGE[0] <= order <= ORDER_RANGE[1]:
            verdict = (f"order unreliable: {order:.2f} is outside "
                       f"[{ORDER_RANGE[0]}, {ORDER_RANGE[1]}]")
        elif not RICHARDSON_RANGE[0] <= order <= RICHARDSON_RANGE[1]:
            verdict = (f"monotone; formal order {order:.2f} reported, far from "
                       f"the nominal {NOMINAL_ORDER:g}, so no Richardson "
                       "estimate")
        else:
            verdict = "monotone; order and Richardson estimate reported"
            extrapolated = f3 + e32 / (RATIO**order - 1.0)
            out["richardson_estimate"] = extrapolated
            out["richardson_relative_error_fine"] = (
                abs(f3 - extrapolated) / abs(extrapolated)
            )
    out["verdict"] = verdict
    # Roache's GCI with the safety factor 3 and the nominal order, the
    # recommended band when the observed order is not reliable; for a
    # non-monotone sequence the larger change is used.
    change = abs(e32) if e21 * e32 >= 0.0 else max(abs(e21), abs(e32))
    out["uncertainty_band_Fs3_p2"] = 3.0 * change / (RATIO**NOMINAL_ORDER - 1.0)
    out["relative_uncertainty_band_Fs3_p2"] = (
        out["uncertainty_band_Fs3_p2"] / scale if scale else math.nan
    )
    return out


def rel(value: float, reference: float) -> float:
    return (value - reference) / abs(reference)


def main() -> int:
    evidence: dict = {
        "description": "3dTube verification evidence (point A, z = 2.5 cm, "
                       "inner wall); generated by analyse_3dTube_evidence.py",
        "definitions": {
            "ur_max_A_m": "peak radial displacement at A, parabolic fit "
                          "through the largest sample and its neighbours",
            "uz_min_A_m": "incident trough of the axial displacement at A, "
                          "minimum over t < 6.5 ms",
            "ur_min_late_A_m": "trough of the radial displacement at A after "
                               "the outlet reflection, t > 14 ms",
            "t_arrival_A_s": "first time the radial displacement at A reaches "
                             "half its first peak",
            "wave_speed_pressure_m_s": "least-squares slope of the half-inlet-"
                                       "pressure arrival on the axis, "
                                       "z = 1-4 cm",
            "wave_speed_wall_m_s": "least-squares slope of the wall half-first-"
                                   "peak arrival, z = 1-4 cm",
        },
    }

    # Iterative-error floor from the coupling study
    coupling = {row["coupling"]: row for row in read_csv("coupling_backward.csv")}
    floors = {}
    if {"robin", "iqnils"} <= set(coupling):
        r, q = coupling["robin"], coupling["iqnils"]
        evidence["coupling"] = {
            name: {key: r_or_q[key] for key in (
                "ur_max_A_m", "uz_min_A_m", "t_arrival_A_s",
                "wave_speed_pressure_m_s", "mean_fsi_iterations",
                "maximum_fsi_iterations", "total_fsi_iterations",
                "clock_time_s", "cores")}
            for name, r_or_q in (("robin", r), ("iqnils", q))
        }
        evidence["coupling"]["radial_history_difference_vs_iqnils"] = r[
            "radial_A_history_difference_vs_iqnils"]
        for quantity, *_ in QUANTITIES:
            if quantity in r and quantity in q:
                floors[quantity] = abs(r[quantity] - q[quantity]) / abs(q[quantity])
                evidence["coupling"].setdefault("relative_difference", {})[
                    quantity] = floors[quantity]

    # Mesh study
    mesh = read_csv("mesh_robin_backward.csv")
    if len(mesh) < 3:
        print("ERROR: the mesh study CSV does not have three levels",
              file=sys.stderr)
        return 1
    levels = []
    for row in mesh:
        levels.append({key: row[key] for key in (
            "level", "fluid_cells", "solid_cells", "delta_t_s", "time_scheme",
            "coupling", "cores", "clock_time_s", "mean_fsi_iterations",
            "maximum_fsi_iterations", "total_fsi_iterations",
            "robin_stalled_steps", "unconverged_steps",
            "maximum_final_displacement_residual",
            "maximum_final_pressure_residual", "maximum_final_flux_residual",
            *[q for q, *_ in QUANTITIES])})
    evidence["mesh_study"] = {"levels": levels, "quantities": {}}
    for quantity, *_rest in QUANTITIES:
        extraction_floor = _rest[2]
        floor = max(extraction_floor, floors.get(quantity, 0.0))
        result = classify([row[quantity] for row in mesh[:3]], floor)
        result["floor_relative"] = floor
        evidence["mesh_study"]["quantities"][quantity] = result
    # The axis pressure probes of the default runs lie on cell faces, so
    # their pressure-front speed carries a sampling error of up to one axial
    # cell; the probe-shifted runs (identical otherwise) give c_p.
    shifted = read_csv("mesh_robin_backward_pz7.8125e-05.csv")
    if len(shifted) >= 3:
        quantity = "wave_speed_pressure_m_s"
        result = classify([row[quantity] for row in shifted[:3]],
                          max(1e-4, floors.get(quantity, 0.0)))
        result["values"] = [row[quantity] for row in shifted[:3]]
        result["floor_relative"] = max(1e-4, floors.get(quantity, 0.0))
        result["note"] = ("probes shifted 7.8125e-5 m axially, inside one "
                          "cell on every level; the default-probe values "
                          "carry a one-cell sampling ambiguity")
        evidence["mesh_study"]["quantities"][
            "wave_speed_pressure_shifted_probes_m_s"] = result
        evidence["mesh_study"]["shifted_probe_runs_identical"] = {
            q: max(abs(a[q] - b[q]) / abs(b[q]) for a, b in zip(shifted, mesh))
            for q in ("ur_max_A_m", "uz_min_A_m", "t_arrival_A_s")
        }
    evidence["mesh_study"]["history_differences"] = {
        "level_1_to_2": history_difference("robin_backward_mesh1",
                                           "robin_backward_mesh2"),
        "level_2_to_3": history_difference("robin_backward_mesh2",
                                           "robin_backward_mesh3"),
    }

    # Iterative error: every level again with tight tolerances
    tight = {int(row["level"]): row for row in read_csv("mesh_robin_backward_tight.csv")}
    if tight:
        evidence["iterative_error_tight_tolerances"] = {
            f"level{level}": {
                **{q: {"default": mesh[level - 1][q], "tight": row[q],
                       "relative_difference": rel(mesh[level - 1][q], row[q])}
                   for q, *_ in QUANTITIES},
                "mean_fsi_iterations_tight": row["mean_fsi_iterations"],
                "clock_time_s_tight": row["clock_time_s"],
            }
            for level, row in sorted(tight.items()) if level <= len(mesh)
        }

    # Spatial versus temporal error
    separation = {}
    timestep = read_csv("timestep_robin_backward.csv")
    if timestep:
        separation["timestep_level1_mesh_backward"] = [
            {key: row[key] for key in ("delta_t_s", "mean_fsi_iterations",
                                       "clock_time_s",
                                       *[q for q, *_ in QUANTITIES])}
            for row in timestep
        ]
    fixed = read_csv("mesh_robin_backward_dt2.5e-05.csv")
    if len(fixed) >= 2 and len(mesh) >= 2:
        f1, f2 = fixed[0], fixed[1]
        separation["fixed_dt_2.5e-5"] = {
            quantity: {
                "level1": f1[quantity], "level2": f2[quantity],
                "spatial_change_level1_to_2_relative": rel(f2[quantity], f1[quantity]),
                "temporal_change_on_level2_dt_2.5e-5_to_1.25e-5_relative":
                    rel(mesh[1][quantity], f2[quantity]),
                "combined_change_level1_to_2_relative":
                    rel(mesh[1][quantity], mesh[0][quantity]),
            }
            for quantity, *_ in QUANTITIES
        }
    if len(fixed) >= 3:
        separation["fixed_dt_2.5e-5_three_levels"] = {
            quantity: {
                "values": [row[quantity] for row in fixed[:3]],
                **classify([row[quantity] for row in fixed[:3]],
                           max(extraction, floors.get(quantity, 0.0))),
            }
            for quantity, _, _, extraction in QUANTITIES
        }
        # Time-step error on the level-3 mesh: dt 2.5e-5 against the
        # level-3 time step 6.25e-6 (and 1.25e-5 when available)
        level3 = {"2.5e-05": fixed[2], "6.25e-06": mesh[2]}
        halved = read_csv("mesh_robin_backward_dt1.25e-05.csv")
        if halved:
            level3["1.25e-05"] = halved[-1]
        separation["level3_time_step"] = {
            quantity: {dt: row[quantity] for dt, row in level3.items()}
            | {"relative_change_2.5e-5_to_6.25e-6":
               rel(mesh[2][quantity], fixed[2][quantity])}
            for quantity, *_ in QUANTITIES
        }
    evidence["space_time_separation"] = separation

    # Cross-platform replicate (MeluXina, OpenFOAM v2412 EasyBuild)
    platform = RESULTS / "3dTube_meluxina_runs.csv"
    if platform.is_file():
        lines = [line for line in platform.read_text().splitlines()
                 if not line.startswith("#")]
        runs = {row["run"]: row for row in csv.DictReader(lines)}
        pairs = {  # MeluXina run: (XenoSim row, label)
            "m_fx1": (fixed[0] if fixed else None, "level 1, dt 2.5e-5"),
            "m_fx2": (fixed[1] if len(fixed) > 1 else None, "level 2, dt 2.5e-5"),
            "m_fx3": (fixed[2] if len(fixed) > 2 else None, "level 3, dt 2.5e-5"),
            "m_l3": (mesh[2], "level 3, dt 6.25e-6"),
        }
        names = {"ur_max_mm": "ur_max_A_m", "uz_min_mm": "uz_min_A_m",
                 "t_arr_ms": "t_arrival_A_s", "ur_min_late_mm": "ur_min_late_A_m"}
        comparison = {}
        for run, (row, label) in pairs.items():
            if row is None or run not in runs:
                continue
            comparison[label] = {
                key: {"meluxina": float(runs[run][key]),
                      "xenosim": row[quantity] * 1e3,
                      "relative_difference":
                          rel(float(runs[run][key]), row[quantity] * 1e3)}
                for key, quantity in names.items()
            }
        evidence["platform_comparison"] = {
            "note": ("XenoSim: OpenFOAM v2512 (Ubuntu package), PETSc 3.24 "
                     "development; MeluXina: OpenFOAM v2412 EasyBuild "
                     "foss-2024a, PETSc 3.22. Each platform is invariant to "
                     "MPI ranks (1-8), solid preconditioner (hypre/LU) and "
                     "tight tolerances; XenoSim v2412 equals XenoSim v2512."),
            "pairs": comparison,
            "meluxina_fixed_dt_uz_min": classify(
                [float(runs[r]["uz_min_mm"]) for r in ("m_fx1", "m_fx2", "m_fx3")],
                1e-4),
            "meluxina_fixed_dt_ur_max": classify(
                [float(runs[r]["ur_max_mm"]) for r in ("m_fx1", "m_fx2", "m_fx3")],
                1e-4),
        }

    # Implicit Euler
    euler = {}
    euler_mesh = read_csv("mesh_robin_Euler_dt0.0001.csv")
    euler_dt = read_csv("timestep_robin_Euler.csv")
    published = REFERENCE["published"]
    fe_peaks = {
        "lozovskiy2019": published["lozovskiy2019"]["values"]["ur_max_A_m"]["value"],
        "eken2016": published["eken2016"]["values"]["ur_max_A_m"]["value"],
    }
    tukovic_peak = published["tukovic2018"]["values"]["ur_max_A_m"]["dt_2.5e-5_s"]
    reading = 3e-6  # raster-reading accuracy of the published FE histories
    bdf2_fine = mesh[2]["ur_max_A_m"]
    if euler_mesh:
        euler["euler_dt1e-4_by_level"] = [
            {key: row[key] for key in ("level", "fluid_cells", "delta_t_s",
                                       *[q for q, *_ in QUANTITIES])}
            for row in euler_mesh
        ]
    if euler_dt:
        euler["euler_timestep_level1_mesh"] = [
            {key: row[key] for key in ("delta_t_s",
                                       *[q for q, *_ in QUANTITIES])}
            for row in euler_dt
        ]
    if euler_mesh:
        fractions = {}
        for row in euler_mesh:
            damping = bdf2_fine - row["ur_max_A_m"]
            for source, peak in fe_peaks.items():
                for basis_name, basis in (("bdf2_level3", bdf2_fine),
                                          ("tukovic2018", tukovic_peak)):
                    gap = basis - peak
                    fractions[f"level{int(row['level'])}_{source}_vs_{basis_name}"] = {
                        "published_gap_m": gap,
                        "euler_damping_m": damping,
                        "fraction_reproduced": damping / gap,
                        "fraction_range_reading_accuracy": [
                            damping / (gap + reading), damping / (gap - reading)
                        ],
                        "euler_vs_published_relative": rel(row["ur_max_A_m"], peak),
                    }
        euler["fraction_of_published_gap_reproduced"] = fractions
        euler["bdf2_level3_ur_max_A_m"] = bdf2_fine
        euler["published_fe_ur_max_A_m"] = fe_peaks
        euler["tukovic2018_ur_max_A_m"] = tukovic_peak
    evidence["implicit_euler"] = euler

    # Literature comparison of the finest BDF2 level
    finest = mesh[2]
    evidence["literature_comparison_level3_bdf2"] = {
        "tukovic2018": {
            "ur_max_A_m": rel(finest["ur_max_A_m"], tukovic_peak),
            "t_ur_max_A_s": rel(finest["t_ur_max_A_s"],
                                published["tukovic2018"]["values"]["t_ur_max_A_s"]["value"]),
        },
        **{
            source: {
                quantity: rel(finest[quantity], published[source]["values"][quantity]["value"])
                for quantity in ("ur_max_A_m", "uz_min_A_m", "ur_min_late_A_m")
            }
            for source in ("lozovskiy2019", "eken2016")
        },
        "thick_wall_inertia_wave_speed": rel(
            finest["wave_speed_pressure_m_s"],
            finest.get("thick_wall_inertia_wave_speed_m_s", 4.81)),
    }

    RESULTS.mkdir(exist_ok=True)
    (RESULTS / "3dTube_evidence.json").write_text(
        json.dumps(evidence, indent=2, default=float) + "\n")
    with (RESULTS / "3dTube_mesh_study.csv").open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow([
            "quantity", "unit", "level1", "level2", "level3",
            "rel_change_1_2", "rel_change_2_3", "observed_order",
            "richardson_estimate", "uncertainty_band_Fs3_p2", "floor", "verdict"])
        for quantity, scale, unit, _ in QUANTITIES:
            q = evidence["mesh_study"]["quantities"][quantity]
            values = [row[quantity] * scale for row in mesh[:3]]
            writer.writerow([
                quantity, unit, *[f"{v:.6g}" for v in values],
                f"{q['relative_change_1_2']:.3e}", f"{q['relative_change_2_3']:.3e}",
                "" if q["observed_order"] is None else f"{q['observed_order']:.3f}",
                "" if q["richardson_estimate"] is None
                else f"{q['richardson_estimate'] * scale:.6g}",
                f"{q['uncertainty_band_Fs3_p2'] * scale:.3g}",
                f"{q['floor_relative']:.1e}",
                q["verdict"] + ("; default probes on cell faces, tie-affected"
                                if quantity == "wave_speed_pressure_m_s" else "")])
        shifted_cp = evidence["mesh_study"]["quantities"].get(
            "wave_speed_pressure_shifted_probes_m_s")
        if shifted_cp:
            q = shifted_cp
            writer.writerow([
                "wave_speed_pressure_shifted_probes_m_s", "m/s",
                *[f"{v:.6g}" for v in q["values"]],
                f"{q['relative_change_1_2']:.3e}", f"{q['relative_change_2_3']:.3e}",
                "" if q["observed_order"] is None else f"{q['observed_order']:.3f}",
                "" if q["richardson_estimate"] is None
                else f"{q['richardson_estimate']:.6g}",
                f"{q['uncertainty_band_Fs3_p2']:.3g}",
                f"{q['floor_relative']:.1e}",
                q["verdict"] + "; preferred c_p (probes off the cell faces)"])
    with (RESULTS / "3dTube_runs.csv").open("w", newline="") as handle:
        writer = csv.writer(handle)
        columns = ["study", "level", "time_scheme", "coupling", "fluid_cells",
                   "solid_cells", "delta_t_s", "cores", "clock_time_s",
                   "mean_fsi_iterations", "maximum_fsi_iterations",
                   *[q for q, *_ in QUANTITIES]]
        writer.writerow(["run", *columns])
        for name in ("mesh_robin_backward.csv",
                     "mesh_robin_backward_dt2.5e-05.csv",
                     "timestep_robin_backward.csv",
                     "timestep_robin_Euler.csv",
                     "mesh_robin_Euler_dt0.0001.csv",
                     "literature_robin_Euler.csv",
                     "coupling_backward.csv"):
            for row in read_csv(name):
                writer.writerow([name.removesuffix(".csv"),
                                 *[row.get(c, "") for c in columns]])
    print(json.dumps(evidence["mesh_study"]["quantities"], indent=1, default=float))
    print(f"Wrote {RESULTS}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
