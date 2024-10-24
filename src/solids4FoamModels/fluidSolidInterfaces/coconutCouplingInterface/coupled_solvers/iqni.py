from coconut import tools
from coconut.coupling_components.coupled_solvers.coupled_solver import CoupledSolver


def create(parameters):
    return CoupledSolverIQNI(parameters)


class CoupledSolverIQNI(CoupledSolver):
    def __init__(self, parameters):
        super().__init__(parameters)

        tools.pass_on_parameters(self.settings, self.settings['model']['settings'],
                                 ('timestep_start', 'delta_t', 'save_restart'))

        self.model = tools.create_instance(self.parameters['settings']['model'])
        self.omega = self.settings['omega']

        self.restart_model = False  # indicates if model has to be restarted

    def initialize(self):
        super().initialize(print_components=False)

        self.model.size_in = self.model.size_out = self.x.size
        self.model.out = self.x.copy()
        self.model.initialize()
        if self.restart_model:  # restart with the same model type
            self.model.restart(self.restart_data['model'])
        self.components += [self.model]

        if self.solver_level == 0:
            self.print_components_info(' ')

    def solve_solution_step(self):
        # initial value
        self.x = self.predictor.predict(self.x)
        # first coupling iteration
        self.y = self.solver_wrappers[0].solve_solution_step(self.x.copy()).copy()
        xt = self.solver_wrappers[1].solve_solution_step(self.y.copy()).copy()
        r = xt - self.x
        self.model.add(r.copy(), xt)
        self.finalize_iteration(r)
        # coupling iteration loop
        while not self.convergence_criterion.is_satisfied():
            if not self.model.is_ready():
                dx = self.omega * r
            else:
                dr = -1 * r
                dx = self.model.predict(dr.copy()) - dr
            self.x += dx
            self.y = self.solver_wrappers[0].solve_solution_step(self.x.copy()).copy()
            xt = self.solver_wrappers[1].solve_solution_step(self.y.copy()).copy()
            r = xt - self.x
            self.model.add(r.copy(), xt)
            self.finalize_iteration(r)

    def check_restart_data(self, restart_data, coupled_solver_settings=None):
        continue_check = super().check_restart_data(restart_data, ['omega'])
        if continue_check:
            old_model_type = restart_data['parameters']['settings']['model']['type']
            new_model_type = self.parameters['settings']['model']['type']
            if new_model_type != old_model_type:
                tools.print_info(f'Model type changed from "{old_model_type}" to "{new_model_type}"', layout='blue')
            else:
                self.restart_model = True
                self.model.check_restart_data(restart_data['parameters']['settings']['model'])

    def add_restart_data(self, restart_data):
        restart_data.update({'model': self.model.save_restart_data()})
