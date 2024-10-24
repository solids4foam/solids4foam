from coconut.tools import print_info
from coconut.coupling_components.coupled_solvers.coupled_solver import CoupledSolver


def create(parameters):
    return CoupledSolverOneWay(parameters)


class CoupledSolverOneWay(CoupledSolver):
    def __init__(self, parameters):
        super().__init__(parameters)
        print_info('CoupledSolverOneWay is chosen: convergence criterion and predictor are ignored', layout='info')

    def solve_solution_step(self):
        self.y = self.solver_wrappers[0].solve_solution_step(self.x.copy()).copy()
        xt = self.solver_wrappers[1].solve_solution_step(self.y.copy())
        r = xt - self.x
        self.finalize_iteration(r)
