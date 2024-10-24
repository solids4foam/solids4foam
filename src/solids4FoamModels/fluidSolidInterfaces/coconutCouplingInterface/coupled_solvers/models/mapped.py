from coconut import tools
from coconut.coupling_components.component import Component


def create(parameters):
    return ModelMapped(parameters)


class ModelMapped(Component):
    mapped = True

    @tools.time_initialize
    def __init__(self, parameters):
        super().__init__()

        # read parameters
        self.parameters = parameters
        self.settings = parameters['settings']

        # create surrogate
        if 'settings' not in self.settings['surrogate']:
            self.settings['surrogate']['settings'] = {}
        tools.pass_on_parameters(self.settings, self.settings['surrogate']['settings'],
                                 ['timestep_start', 'delta_t', 'save_restart', 'case_name'])
        self.surrogate = tools.create_instance(self.settings['surrogate'])
        self.provides_get_solution = self.surrogate.provides_get_solution
        self.provides_set_solution = self.surrogate.provides_set_solution

        # create mappers
        self.mapper_interface_input = tools.create_instance(self.settings['mapper_interface_input'])
        self.mapper_interface_output = tools.create_instance(self.settings['mapper_interface_output'])

        self.interface_input_from = None
        self.interface_input_to = None
        self.interface_output_to = None
        self.out = None  # interface of output
        self.solver_level = None

        # time
        # noinspection PyUnresolvedReferences
        self.init_time = self.init_time  # created by decorator time_initialize
        self.run_time = 0.0
        self.save_time = 0.0

    @tools.time_initialize
    def initialize(self, interface_input_from, interface_output_to):
        super().initialize()

        self.surrogate.solver_level = self.solver_level
        self.surrogate.initialize()

        # initialize input mapper
        self.interface_input_from = interface_input_from.copy()
        self.interface_input_to = self.surrogate.get_interface_input()
        self.mapper_interface_input.initialize(self.interface_input_from, self.interface_input_to)

        # initialize output mapper
        self.interface_output_to = interface_output_to.copy()
        interface_output_from = self.surrogate.get_interface_output()
        self.mapper_interface_output.initialize(interface_output_from, self.interface_output_to)

    def initialize_solution_step(self):
        super().initialize_solution_step()

        self.surrogate.initialize_solution_step()

    def predict(self, interface_input_from, **kwargs):
        self.interface_input_from = interface_input_from.copy()
        self.mapper_interface_input(self.interface_input_from, self.interface_input_to)
        interface_output_from = self.surrogate.predict(self.interface_input_to, **kwargs)
        self.mapper_interface_output(interface_output_from, self.interface_output_to)
        return self.interface_output_to.copy()

    @tools.time_solve_solution_step
    def get_solution(self):
        interface_output_from = self.surrogate.get_solution()
        self.mapper_interface_output(interface_output_from, self.interface_output_to)
        return self.interface_output_to.copy()

    @tools.time_solve_solution_step
    def set_solution(self, interface_input_from):
        self.interface_input_from = interface_input_from.copy()
        self.mapper_interface_input(self.interface_input_from, self.interface_input_to)
        self.surrogate.set_solution(self.interface_input_to)

    def add(self, interface_input_from_1, interface_input_from_2):
        self.interface_input_from = interface_input_from_1.copy()
        self.mapper_interface_input(self.interface_input_from, self.interface_input_to)
        interface_input_to_1 = self.interface_input_to.copy()
        self.interface_input_from = interface_input_from_2.copy()
        self.mapper_interface_input(self.interface_input_from, self.interface_input_to)
        interface_input_to_2 = self.interface_input_to.copy()
        self.surrogate.add(interface_input_to_1, interface_input_to_2)

    def is_ready(self):
        return self.surrogate.is_ready()

    def filter_q(self, interface_input_from, **kwargs):
        self.interface_input_from = interface_input_from.copy()
        self.mapper_interface_input(self.interface_input_from, self.interface_input_to)
        interface_output_from = self.surrogate.filter_q(self.interface_input_to, **kwargs)
        self.mapper_interface_output(interface_output_from, self.interface_output_to)
        return self.interface_output_to.copy()

    def finalize_solution_step(self):
        super().finalize_solution_step()

        self.surrogate.finalize_solution_step()

    @tools.time_save
    def output_solution_step(self):
        super().output_solution_step()

        self.surrogate.output_solution_step()
        self.mapper_interface_input.output_solution_step()
        self.mapper_interface_output.output_solution_step()

    def finalize(self):
        super().finalize()

        self.surrogate.finalize()
        self.mapper_interface_input.finalize()
        self.mapper_interface_output.finalize()

    def get_interface_input(self):
        # does not contain most recent data
        return self.interface_input_from.copy()

    def get_interface_output(self):
        interface_output_from = self.surrogate.get_interface_output()
        self.mapper_interface_output(interface_output_from, self.interface_output_to)
        return self.interface_output_to.copy()

    def restart(self, restart_data):
        self.surrogate.restart(restart_data)

    def check_restart_data(self, restart_data):
        self.surrogate.check_restart_data(restart_data)

    def save_restart_data(self):
        self.surrogate.save_restart_data()

    def get_time_allocation(self):
        time_allocation = {}
        for time_type in ('init_time', 'run_time', 'save_time'):
            total_time = self.__getattribute__(time_type)
            surrogate_time = self.surrogate.get_time_allocation()[time_type]
            mapper_time = total_time - (surrogate_time['total'] if isinstance(surrogate_time, dict) else surrogate_time)
            time_allocation[time_type] = {'total': total_time, 'mapper': mapper_time, 'surrogate': surrogate_time}
        return time_allocation

    def print_time_distribution(self, ta, pre, reference=None, last=True):
        return self.surrogate.print_time_distribution(ta, pre, reference=reference, last=last)

    def print_components_info(self, pre):
        tools.print_info(pre, 'The component ', self.__class__.__name__, ' maps the following surrogate model:')
        pre = tools.update_pre(pre)
        self.surrogate.print_components_info(pre + '├─')
        tools.print_info(pre, '├─', 'Input mapper:')
        self.mapper_interface_input.print_components_info(pre + '│ └─')
        tools.print_info(pre, '└─', 'Output mapper:')
        self.mapper_interface_output.print_components_info(pre + '  └─')
