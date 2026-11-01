#!/usr/bin/env python3
"""Extract the small CSV tables behind Figures 2 and 3 from the diagnostics (rtDiagFluid.dat, rtDiagSolid.dat, reference.json) of the runs.
Read-only; no CFD. Uses the tutorial's own verification/scripts/analyse_case.py loader (same as the verification driver).
Usage: make_figure_data.py DIAG_DIR TUTORIAL_VERIFICATION_SCRIPTS_DIR [OUT_DIR]
DIAG_DIR/<case>/ must contain rtDiagFluid.dat, rtDiagSolid.dat and reference.json (for parallel runs the per-rank files recombined into one set).
Formulations (by dictionaries): A = tutorial (leastSquares p and U, stock newMovingWallVelocity); B = A with grad(p) Gauss linear only."""
import csv, os, sys
import numpy as np
RUNS = [  # formulation, member, mesh, case, nominal dt
 ("A", "torsion", "M1", "F2_M1T_nNonOrth2", 0.0025), ("A", "torsion", "M2", "F2_M2T_nNonOrth2", 0.0025), ("A", "torsion", "M3", "F1_M3T_nNonOrth2_clean", 0.0025),
 ("A", "combined", "M1", "F3_M1C_nNonOrth2", 0.0025), ("A", "combined", "M2", "F3_M2C_nNonOrth2", 0.0025), ("A", "combined", "M3", "F3_M3C_nNonOrth2", 0.0025),
 ("B", "combined", "M1", "C1_M1C_Gauss_dt2p5e-3", 0.0025), ("B", "combined", "M2", "C2_M2C_Gauss_dt2p5e-3", 0.0025), ("B", "combined", "M3", "C3_M3C_Gauss_dt2p5e-3", 0.0025)]
TSTEP = [("P0", "P0_M2C_A_np2_dt2p5e-3", 0.0025, "control (2 ranks)"), ("P1", "P1_M2C_A_np2_dt1p25e-3", 0.00125, "dt/2"), ("P2", "P2_M2C_A_np2_dt0p625e-3", 0.000625, "dt/4"),
         ("P3", "P3_M2C_A_np2_strict_tol1e-7", 0.0025, "strict coupling, tol 1e-7")]
QOIS = ["alpha_B", "T_B", "sig_rt", "T_A", "H_f", "H_s", "w_wall", "P_wall", "KE_s", "u_B", "p_b", "sig_nv", "sig_rr"]
NOISE = {("torsion", "alpha_B"): 6.9e-07, ("torsion", "T_B"): 2.4e-4, ("combined", "alpha_B"): 5.7e-06, ("combined", "T_B"): 1.3e-3}   # of max|ref|, verification/thresholds.json numerical_noise
def main(diag, scripts, out="."):
    sys.path.insert(0, scripts); sys.dont_write_bytecode = True
    import analyse_case as AC
    def load_fine(case, dt):
        """Errors are always evaluated on the common 2.5 ms grid: finer runs are subsampled (every 2nd / 4th row; the diagnostics are
        time-level values, not averages) and compared with reference.json of the 2.5 ms control, which is exact at those times."""
        d = os.path.join(diag, case)
        if dt < 0.0025 - 1e-12:
            import tempfile, shutil
            sub = tempfile.mkdtemp(prefix="sub_")
            for fn in ("rtDiagFluid.dat", "rtDiagSolid.dat"):
                x = np.loadtxt(os.path.join(d, fn), ndmin=2); k = np.round(x[:, 0] / 0.0025, 6); m = np.abs(k - np.round(k)) < 1e-6
                np.savetxt(os.path.join(sub, fn), x[m], fmt="%.16e")
            shutil.copy(os.path.join(diag, "F3_M2C_nNonOrth2", "reference.json"), sub)
            d = sub
        t, num, ref, _ = AC.load(d); return t, num, ref
    rows = []
    for f, mem, mesh, case, dt in RUNS:
        t, num, ref = load_fine(case, dt)
        for q in QOIS:
            if q not in num or q not in ref: continue
            pk = float(np.max(np.abs(ref[q])))
            if not pk > 0: continue
            e = num[q][:len(t)] - ref[q][:len(t)]
            rows.append(dict(formulation=f, member=mem, mesh=mesh, case=case, dt=dt, qoi=q, peak=pk, E_max_over_peak=float(np.max(np.abs(e)) / pk), E_L2t_over_peak=float(np.sqrt(np.mean(e ** 2)) / pk), noise_over_peak=NOISE.get((mem, q), ""), ))
    with open(os.path.join(out, "data", "fig2_convergence.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, list(rows[0]), lineterminator="\n"); w.writeheader(); w.writerows(rows)
    rows = []
    base = {}
    for tag, case, dt, label in TSTEP:
        if not os.path.isdir(os.path.join(diag, case)): continue
        t, num, ref = load_fine(case, dt); n = len(t)
        t0, num0, ref0 = load_fine(TSTEP[0][1], TSTEP[0][2])          # P0 on the same window as this run
        for q in QOIS:
            if q not in num or q not in ref: continue
            pk = float(np.max(np.abs(ref0[q]))); e = num[q][:n] - ref[q][:n]; e0 = num0[q][:n] - ref0[q][:n]
            rows.append(dict(run=tag, case=case, dt=dt, label=label, window_end_s=round(float(t[n - 1]), 4), qoi=q, peak=pk, E_max_over_peak=float(np.max(np.abs(e)) / pk),
                             E_P0_same_window=float(np.max(np.abs(e0)) / pk), E_L2t_over_peak=float(np.sqrt(np.mean(e ** 2)) / pk),
                             max_abs_diff_to_P0_over_peak=float(np.max(np.abs(num[q][:n] - num0[q][:n])) / pk)))
    if rows:
        with open(os.path.join(out, "data", "fig2_timestep_strict.csv"), "w", newline="") as fh:
            w = csv.DictWriter(fh, list(rows[0]), lineterminator="\n"); w.writeheader(); w.writerows(rows)
    rows = []
    for f, case in (("A", "F3_M2C_nNonOrth2"), ("B", "C2_M2C_Gauss_dt2p5e-3")):
        t, num, ref = load_fine(case, 0.0025)
        for i in range(len(t)):
            rows.append(dict(formulation=f, case=case, t=round(float(t[i]), 4), ref_p_b=ref["p_b"][i], cfd_p_b=num["p_b"][i], ref_sig_nv=ref["sig_nv"][i], cfd_sig_nv=num["sig_nv"][i], ref_sig_rr=ref["sig_rr"][i], cfd_sig_rr=num["sig_rr"][i]))
    with open(os.path.join(out, "data", "fig3_traction_timeseries.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, list(rows[0]), lineterminator="\n"); w.writeheader()
        for r in rows: w.writerow({k: (("%.6g" % v) if isinstance(v, float) else v) for k, v in r.items()})
    print("written CSVs to", os.path.join(out, "data"))
if __name__ == "__main__":
    if len(sys.argv) < 3: print(__doc__); sys.exit(2)
    main(sys.argv[1], sys.argv[2], sys.argv[3] if len(sys.argv) > 3 else ".")
