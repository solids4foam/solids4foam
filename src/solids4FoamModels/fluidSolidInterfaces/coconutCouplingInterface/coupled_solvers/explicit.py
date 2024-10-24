from coconut.coupling_components.coupled_solvers.coupled_solver import CoupledSolver
from coconut.tools import print_info


def create(parameters):
    return CoupledSolverExplicit(parameters)


class CoupledSolverExplicit(CoupledSolver):
    def __init__(self, parameters):
        super().__init__(parameters)
        print_info(f'CoupledSolverExplicit is chosen: the convergence criterion is not used', layout='info')

    def solve_solution_step(self):
        # initial value
        self.x = self.predictor.predict(self.x)
        # first coupling iteration
        self.y = self.solver_wrappers[0].solve_solution_step(self.x.copy()).copy()
        xt = self.solver_wrappers[1].solve_solution_step(self.y.copy()).copy()
        r = xt - self.x
        self.finalize_iteration(r)
        self.x = xt
