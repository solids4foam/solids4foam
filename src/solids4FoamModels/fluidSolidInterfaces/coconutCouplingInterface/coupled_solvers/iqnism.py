from coconut import tools
from coconut.coupling_components.coupled_solvers.coupled_solver import CoupledSolver

import time


def create(parameters):
    return CoupledSolverIQNISM(parameters)


class CoupledSolverIQNISM(CoupledSolver):
    def __init__(self, parameters):
        super().__init__(parameters)

        # add timestep_start, delta_t, save_restart and case_name to surrogate settings
        if 'settings' not in self.settings['surrogate']:
            self.settings['surrogate']['settings'] = {}
        tools.pass_on_parameters(self.settings, self.settings['surrogate']['settings'], ['timestep_start', 'delta_t',
                                                                                         'save_restart', 'case_name'])
        # add timestep_start, delta_t, save_restart to model settings
        tools.pass_on_parameters(self.settings, self.settings['model']['settings'],
                                 ('timestep_start', 'delta_t', 'save_restart'))

        self.model = tools.create_instance(self.settings['model'])
        self.settings['surrogate']['type'] = self.settings['surrogate'].get('type',
                                                                            'coupled_solvers.models.dummy_model')
        self.surrogate = tools.create_instance(self.settings['surrogate'], 'models.dummy_model')
        self.omega = self.settings.get('omega', 1)  # relaxation factor
        self.surrogate_modes = self.settings.get('surrogate_modes')  # number of surrogates modes to use
        self.surrogate_synchronize = self.settings.get('surrogate_synchronize', True)  # synchronize surrogate

        if (self.surrogate_modes == 0 or self.surrogate.dummy) and 'omega' not in self.settings:
            raise ValueError('A relaxation factor (omega) is required when no surrogate Jacobian is used')

        self.restart_model_surrogate = [False, False]  # indicates if model and/or surrogate has to be restarted

    def initialize(self):
        super().initialize(print_components=False)

        self.surrogate.solver_level = self.solver_level + 1

        # initialize mapper for surrogate if required
        if self.surrogate.mapped:
            interface_input_from = self.solver_wrappers[0].get_interface_input()
            interface_output_to = self.solver_wrappers[1].get_interface_output()  # same interface a input_from
            self.surrogate.initialize(interface_input_from, interface_output_to)
        else:
            self.surrogate.initialize()

        self.model.size_in = self.model.size_out = self.x.size
        self.model.out = self.x.copy()
        self.model.initialize()

        restart_model, restart_surrogate = self.restart_model_surrogate
        if restart_model:  # restart with the same model type
            self.model.restart(self.restart_data['model'])
        if restart_surrogate:  # restart with the same surrogate type
            self.surrogate.restart(self.restart_data['surrogate'])

        self.components += [self.model, self.surrogate]

        # set initial surrogate value in surrogate predictor
        if self.surrogate.provides_get_solution and not self.restart:
            if hasattr(self.predictor, 'update_surrogate') and not self.surrogate.dummy:
                self.predictor.update_surrogate(self.surrogate.get_interface_output())

        self.start_run_time = time.time()  # reset start of calculation
        self.init_time = self.start_run_time - self.start_init_time  # reset duration of initialization

        if self.solver_level == 0:
            self.print_components_info(' ')

    def solve_solution_step(self):
        # solve surrogate
        if self.surrogate.provides_get_solution:
            x_surrogate = self.surrogate.get_solution()
            if hasattr(self.predictor, 'update_surrogate') and not self.surrogate.dummy:
                self.predictor.update_surrogate(x_surrogate)
        # initial value
        self.x = self.predictor.predict(self.x)
        # first coupling iteration
        self.y = self.solver_wrappers[0].solve_solution_step(self.x.copy()).copy()
        xt = self.solver_wrappers[1].solve_solution_step(self.y.copy()).copy()
        r = xt - self.x
        self.model.add(r.copy(), xt.copy())
        self.surrogate.add(r.copy(), xt)  # only used when derivative info of surrogate is updated every iteration
        self.finalize_iteration(r)
        # coupling iteration loop
        while not self.convergence_criterion.is_satisfied():
            dr = -1 * r
            if not self.model.is_ready():
                if not self.surrogate.is_ready():
                    dx = -self.omega * dr
                else:
                    dx = self.surrogate.predict(dr.copy(), modes=self.surrogate_modes) - dr
                    # relax other modes
                    dx -= (self.omega - 1.0) * self.surrogate.filter_q(dr.copy(), modes=self.surrogate_modes)
            else:
                if not self.surrogate.is_ready():
                    dx = self.model.predict(dr.copy()) - dr
                else:
                    dx = self.model.predict(dr.copy()) + self.surrogate.predict(self.model.filter_q(dr.copy()),
                                                                                modes=self.surrogate_modes) - dr
            self.x += dx
            self.y = self.solver_wrappers[0].solve_solution_step(self.x.copy()).copy()
            xt = self.solver_wrappers[1].solve_solution_step(self.y.copy()).copy()
            r = xt - self.x
            self.model.add(r.copy(), xt.copy())
            self.surrogate.add(r.copy(), xt)  # only used when derivative information of surrogate is function of x
            self.finalize_iteration(r)
        # synchronize
        if self.surrogate_synchronize and self.surrogate.provides_set_solution:
            self.surrogate.set_solution(self.x.copy())

    def check_restart_data(self, restart_data, coupled_solver_settings=None):
        continue_check = super().check_restart_data(restart_data, ['omega', 'surrogate_modes', 'surrogate_synchronize'])
        if continue_check:
            for i, model_name in enumerate(['model', 'surrogate']):
                old_model_type = restart_data['parameters']['settings'][model_name]['type']
                new_model_type = self.parameters['settings'][model_name]['type']
                if new_model_type != old_model_type:
                    tools.print_info(f'Model type changed from "{old_model_type}" to "{new_model_type}"', layout='blue')
                else:
                    self.restart_model_surrogate[i] = True
            if self.restart_model_surrogate[0]:
                self.model.check_restart_data(restart_data['parameters']['settings']['model'])

    def add_restart_data(self, restart_data):
        restart_data.update({'model': self.model.save_restart_data(), 'surrogate': self.surrogate.save_restart_data()})

    def add_time_allocation(self, time_allocation):
        if hasattr(self.surrogate, 'get_time_allocation'):
            time_allocation_surrogate = self.surrogate.get_time_allocation()
            for time_type in ('init_time', 'run_time', 'save_time'):
                time_allocation_sub = time_allocation[time_type]
                tasur_sub = time_allocation_surrogate[time_type]
                time_allocation_sub.update({'surrogate': tasur_sub})
                time_allocation_sub['coupling'] = time_allocation_sub['coupling'] \
                    - tasur_sub if isinstance(tasur_sub, float) else tasur_sub['total']

    def add_time_distribution(self, ta, pre, reference=None):
        out = ''
        if 'surrogate' not in ta:
            return out
        if reference is None:
            reference = ta['total']
        ta_surrogate = ta['surrogate']
        if isinstance(ta_surrogate, float):
            out += f'{pre}├─{self.surrogate.__class__.__name__}: {ta_surrogate:.0f}s ' \
                  f'({ta_surrogate / reference * 100:0.1f}%)\n'
        else:
            out += f'{pre}├─{self.surrogate.__class__.__name__}: {ta_surrogate["total"]:.0f}s ' \
                   f'({ta_surrogate["total"] / reference * 100:0.1f}%)\n'
            if self.surrogate.mapped:
                ta_surrogate = ta_surrogate['surrogate']
                out += f'{pre}│ └─{self.surrogate.surrogate.__class__.__name__}: {ta_surrogate["total"]:.0f}s' \
                       f' ({ta_surrogate["total"] / reference * 100:0.1f}%)\n'
            out += self.surrogate.print_time_distribution(
                ta_surrogate, pre + f'│ {"  " if self.surrogate.mapped else ""}', reference=reference)
        return out
