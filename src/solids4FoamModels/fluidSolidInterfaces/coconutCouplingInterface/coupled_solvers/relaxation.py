from coconut.coupling_components.coupled_solvers.coupled_solver import CoupledSolver


def create(parameters):
    return CoupledSolverRelaxation(parameters)


class CoupledSolverRelaxation(CoupledSolver):
    def __init__(self, parameters):
        super().__init__(parameters)

        self.omega = self.settings['omega']

    def solve_solution_step(self):
        # initial value
        self.x = self.predictor.predict(self.x)
        # first coupling iteration
        self.y = self.solver_wrappers[0].solve_solution_step(self.x.copy()).copy()
        xt = self.solver_wrappers[1].solve_solution_step(self.y.copy()).copy()
        r = xt - self.x
        self.finalize_iteration(r)
        # coupling iteration loop
        while not self.convergence_criterion.is_satisfied():
            self.x += self.omega * r
            self.y = self.solver_wrappers[0].solve_solution_step(self.x.copy()).copy()
            xt = self.solver_wrappers[1].solve_solution_step(self.y.copy()).copy()
            r = xt - self.x
            self.finalize_iteration(r)

    def check_restart_data(self, restart_data, coupled_solver_settings=None):
        super().check_restart_data(restart_data, ['omega'])
