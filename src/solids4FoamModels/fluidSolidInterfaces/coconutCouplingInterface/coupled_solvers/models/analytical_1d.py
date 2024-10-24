from coconut import tools
from coconut.coupling_components.component import Component

import numpy as np


def create(parameters):
    return ModelAnalytical1D(parameters)


class ModelAnalytical1D(Component):
    provides_get_solution = provides_set_solution = False

    def __init__(self, parameters):
        super().__init__()

        # reading
        self.parameters = parameters
        self.settings = parameters['settings']

        self.update_every_iteration = self.settings.get('update_every_iteration', False)  # update Jacobian every iter

        # add timestep_start and delta_t to solver_models settings
        for solver_model_parameters in self.settings['solver_models']:
            tools.pass_on_parameters(self.settings, solver_model_parameters['settings'],
                                     ['timestep_start', 'delta_t', 'save_restart'])

        # initialization
        self.solver_models = []
        self.solver_models.append(tools.create_instance(self.settings['solver_models'][0]))  # flow model
        self.solver_models.append(tools.create_instance(self.settings['solver_models'][1]))  # structure model

        self.interface_input = self.solver_models[0].get_interface_input()  # residual of displacements
        self.interface_output = self.solver_models[1].get_interface_output()  # displacements
        if self.interface_input.size != self.interface_output.size:
            raise Exception('Input and output of surrogate model have to be the same size')

        self.x = None
        self.m = None
        self.out = None
        self.jit = None  # dxt/dr

    def initialize(self):
        super().initialize()

        for solver_model in self.solver_models:
            solver_model.initialize()

        self.x = self.solver_models[1].get_interface_output()
        self.m = int(self.x.get_interface_data().shape[0] / 3)
        self.out = self.x.copy()
        self.jit = self.get_inverse_jacobian()

    def finalize(self):
        super().finalize()

        for solver_model in self.solver_models:
            solver_model.finalize()

    def initialize_solution_step(self):
        super().initialize_solution_step()

        for solver_model in self.solver_models:
            solver_model.initialize_solution_step()

        if not self.update_every_iteration:
            self.jit = self.get_inverse_jacobian()

    def finalize_solution_step(self):
        super().finalize_solution_step()

        # update with correct final values
        self.solver_models[0].solve_solution_step(self.x)

        for solver_model in self.solver_models:
            solver_model.finalize_solution_step()

    def output_solution_step(self):
        super().output_solution_step()

        for solver_model in self.solver_models:
            solver_model.output_solution_step()

    def get_interface_input(self):
        return self.interface_input.copy()

    def get_interface_output(self):
        return self.interface_output.copy()

    def print_components_info(self, pre):
        tools.print_info(pre, 'The surrogate model ', self.__class__.__name__, ' has the following solver models:')
        tools.print_components_info(pre, self.solver_models)

    def predict(self, dr_in, modes=None):
        dr = dr_in.get_interface_data()[1::3]  # only y-displacement
        if modes == 0:
            tools.print_info("Jacobian of model Analytical1D disabled: zero returned", layout='warning')
            return dr_in * 0
        elif modes is not None:
            tools.print_info(f'Mode limiting not possible for Analytical1D model', layout='warning')

        # approximation for the inverse of the Jacobian from a surrogate model
        dxt = self.jit @ dr
        dxt_out = self.out.copy()
        dxt_full = np.zeros(3 * self.m)
        dxt_full[1::3] = dxt
        dxt_out.set_interface_data(dxt_full)
        return dxt_out

    def add(self, r_in, xt_in):
        self.x = xt_in - r_in
        if self.update_every_iteration:
            self.jit = self.get_inverse_jacobian()

    # noinspection PyMethodMayBeStatic
    def is_ready(self):
        return True

    def get_flow_jacobian(self):
        # Jacobian of flow model dp/dr

        # calculate yt
        self.solver_models[0].solve_solution_step(self.x)

        # calculation of Jacobian
        return self.solver_models[0].get_surrogate_jacobian()

    def get_structure_jacobian(self):
        # Jacobian of structure model

        # Calculation of Jacobian
        return self.solver_models[1].get_surrogate_jacobian()

    def get_inverse_jacobian(self):
        # calculation of Jacobian
        jf = self.get_flow_jacobian()  # Jacobian of flow model
        js = self.get_structure_jacobian()  # Jacobian of structure model

        # dr/dx
        j = js @ jf - np.identity(self.m)

        # dx/dr
        ji = np.linalg.inv(j)

        # dxt/dr
        jit = ji + np.identity(self.m)
        return jit

    def filter_q(self, dr_in, modes=None):
        dr = dr_in.get_interface_data().reshape(-1, 1)[1::3]
        dr_out = self.out.copy()
        if modes == 0:
            pass  # return copy of dr_in
        elif modes is None:
            qt, *_ = np.linalg.qr(self.jit.T)
            q = qt[:, :np.linalg.matrix_rank(self.jit)]
            dr = dr - q @ (q.T @ dr)
        else:
            tools.print_info(f'Mode limiting not possible for Analytical1D model', layout='warning')
        r_full = np.zeros(3 * self.m)
        r_full[1::3] = dr.flatten()
        dr_out.set_interface_data(r_full)
        return dr_out

    def restart(self, restart_data):
        pass

    def check_restart_data(self, restart_data):
        pass

    def save_restart_data(self):
        pass

    def get_time_allocation(self):
        time_allocation = {}
        for time_type in ('init_time', 'run_time', 'save_time'):
            time_allocation_sub = time_allocation[time_type] = {'total': 0.0}
            for i, solver_model in enumerate(self.solver_models):
                time_allocation_sub[f'solver_wrapper_{i}'] = solver_model.get_time_allocation()[time_type]
            time_allocation_sub['total'] = sum([s if isinstance(s, float) else s['total']
                                                for s in time_allocation_sub.values()])
        return time_allocation

    def print_time_distribution(self, ta, pre, reference=None, last=True):
        out = ''
        if reference is None:
            reference = ta['total']
        for i, solver in enumerate(self.solver_models):
            ta_solver = ta[f'solver_wrapper_{i}']
            conn_char = '└' if i == len(self.solver_models) - 1 and last else '├'
            if isinstance(ta_solver, float):
                out += f'{pre}{conn_char}─{solver.__class__.__name__}: {ta_solver:.0f}s ' \
                       f'({ta_solver / reference * 100:0.1f}%)\n'
            else:
                out += f'{pre}├─{solver.__class__.__name__}: {ta_solver["total"]:.0f}s ' \
                       f'({ta_solver["total"] / reference * 100:0.1f}%)\n'
                if solver.mapped:
                    out += f'{pre}{conn_char}  └─{solver.solver_wrapper.__class__.__name__}: ' \
                           f'{ta_solver["solver_wrapper"]:.0f}s ' \
                           f'({ta_solver["solver_wrapper"] / reference * 100:0.1f}%)\n'
        return out
