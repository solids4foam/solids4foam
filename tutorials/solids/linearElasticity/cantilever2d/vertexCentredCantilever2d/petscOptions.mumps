# Direct solver
-sub_pc_type lu
-pc_factor_mat_solver_package mumps
-ksp_type preonly

# Monitor residuals
-ksp_monitor_short