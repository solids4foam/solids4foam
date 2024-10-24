from coconut.coupling_components.component import Component
from coconut import tools


def create(parameters):
    return ModelSurrogate(parameters)


class ModelSurrogate(Component):
    provides_get_solution = provides_set_solution = True

    @tools.time_initialize
    def __init__(self, parameters):
        super().__init__()

        self.parameters = parameters
        self.settings = parameters['settings']

        # change default case name to surrogate
        if 'case_name' not in self.parameters['coupled_solver']['settings']:
            self.parameters['coupled_solver']['settings']['case_name'] = f'{self.settings["case_name"]}_surrogate'
        # add timestep_start, delta_t and save_restart to surrogate settings
        tools.pass_on_parameters(self.settings, self.parameters['coupled_solver']['settings'],
                                 ('timestep_start', 'delta_t', 'save_restart'))

        self.coupled_solver = tools.create_instance(self.parameters['coupled_solver'])
        self.solver_level = None

        # time
        # noinspection PyUnresolvedReferences
        self.init_time = self.init_time  # created by decorator time_initialize
        self.run_time = 0.0
        self.save_time = 0.0

    @tools.time_initialize
    def initialize(self):
        super().initialize()

        self.coupled_solver.solver_level = self.solver_level
        self.coupled_solver.initialize()

    def initialize_solution_step(self):
        super().initialize_solution_step()

        self.coupled_solver.initialize_solution_step()

    def finalize_solution_step(self):
        super().finalize_solution_step()

        self.coupled_solver.finalize_solution_step()

    @tools.time_save
    def output_solution_step(self):
        super().output_solution_step()

        self.coupled_solver.output_solution_step()

    def finalize(self):
        super().finalize()

        self.coupled_solver.finalize()

    @tools.time_solve_solution_step
    def get_solution(self):
        # surrogate information
        pre = ' │' * (self.solver_level - 1)
        out = f'{pre} ┌{(78 - len(pre)) * "─"}\n' \
              f'{pre} │\tSurrogate\n' \
              f'{pre} ├{(78 - len(pre)) * "─"}\n' \
              f'{pre} │{"Iteration":<16}{"Norm residual":<28}'
        tools.print_info(out)
        self.coupled_solver.solve_solution_step()
        out = f'{pre} └{(78 - len(pre)) * "─"}\n' \
              f'{pre}{"Iteration":<16}{"Norm residual":<28}'
        tools.print_info(out)
        return self.coupled_solver.x.copy()

    @tools.time_solve_solution_step
    def set_solution(self, x):
        y = self.coupled_solver.solver_wrappers[0].solve_solution_step(x)
        xt = self.coupled_solver.solver_wrappers[1].solve_solution_step(y)
        pre = ' │' * (self.solver_level - 1)
        out = f'{pre} ┌{(78 - len(pre)) * "─"}\n' \
              f'{pre} │{"Synchronization":<16}{"Norm residual":<28}\n' \
              f'{pre} │{1:<16d}{(x - xt).norm():<28.17e}\n' \
              f'{pre} └{(78 - len(pre)) * "─"}\n' \
              f'{pre}'
        tools.print_info(out)

    def predict(self, dr, **kwargs):
        return self.coupled_solver.model.predict(dr, **kwargs)

    def add(self, r, xt):
        pass

    def is_ready(self):
        return self.coupled_solver.model.is_ready()

    def filter_q(self, dr, **kwargs):
        return self.coupled_solver.model.filter_q(dr, **kwargs)

    def get_interface_input(self):
        return self.coupled_solver.x.copy()

    def get_interface_output(self):
        return self.coupled_solver.x.copy()

    def restart(self, restart_data):
        pass

    def check_restart_data(self, restart_data):
        pass

    def save_restart_data(self):
        pass

    def get_time_allocation(self):
        time_allocation = {}
        for time_type in ('init_time', 'run_time', 'save_time'):
            total_time = self.__getattribute__(time_type)
            coupled_solver_time = self.coupled_solver.get_time_allocation()[time_type]
            solution_time = coupled_solver_time.pop('total') - coupled_solver_time.pop(
                'coupling')  # time for solver, surrogates etc
            other_time = total_time - solution_time
            time_allocation[time_type] = {'total': total_time}
            time_allocation[time_type].update(coupled_solver_time)
            time_allocation[time_type].update({'coupling': other_time})
        return time_allocation

    def print_components_info(self, pre):
        tools.print_info(pre, 'The component ', self.__class__.__name__, ' with the following CoupledSolver:')
        pre = tools.update_pre(pre)
        self.coupled_solver.print_components_info(pre + '└─')

    def print_time_distribution(self, ta, pre, reference=None, last=True):
        return self.coupled_solver.print_time_distribution(ta, pre, reference=reference, last=last)
