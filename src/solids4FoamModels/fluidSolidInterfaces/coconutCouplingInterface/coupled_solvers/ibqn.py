from coconut import tools
from coconut.coupling_components.coupled_solvers.coupled_solver import CoupledSolver

from scipy.sparse.linalg import gmres, LinearOperator


def create(parameters):
    return CoupledSolverIBQN(parameters)


class CoupledSolverIBQN(CoupledSolver):
    def __init__(self, parameters):
        super().__init__(parameters)

        for model in ('model_f', 'model_s'):
            tools.pass_on_parameters(self.settings, self.settings[model]['settings'],
                                     ('timestep_start', 'delta_t', 'save_restart'))

        self.model_f = tools.create_instance(self.parameters['settings']['model_f'])
        self.model_s = tools.create_instance(self.parameters['settings']['model_s'])
        self.models = [self.model_f, self.model_s]
        self.omega = self.settings['omega']
        self.atol = self.settings['absolute_tolerance_gmres']
        self.rtol = self.settings['relative_tolerance_gmres']

        self.xtemp = self.ytemp = None
        self.dxtemp = self.dytemp = None
        self.u = self.w = None
        self.ready = None

        self.restart_models = [False, False]  # indicates if model_f and/or model_s have to be restarted

    def initialize(self):
        super().initialize(print_components=False)

        self.dxtemp = self.x.copy()
        self.dytemp = self.y.copy()
        self.u = self.x.size
        self.w = self.y.size
        self.ready = False
        self.model_f.size_in = self.model_s.size_out = self.u
        self.model_f.size_out = self.model_s.size_in = self.w
        self.model_f.out = self.y.copy()
        self.model_s.out = self.x.copy()
        for i, model_name in enumerate(['model_f', 'model_s']):
            model = self.models[i]
            model.initialize()
            if self.restart_models[i]:  # restart with the same model type
                model.restart(self.restart_data[model_name])
            self.components += [model]

        if self.solver_level == 0:
            self.print_components_info(' ')

    def lop_f(self, dx):
        self.dxtemp.set_interface_data(dx.flatten())
        return self.model_f.predict(self.dxtemp).get_interface_data()

    def lop_s(self, dy):
        self.dytemp.set_interface_data(dy.flatten())
        return self.model_s.predict(self.dytemp).get_interface_data()

    # noinspection PyMethodMayBeStatic
    def identity_matvec(self, v):
        return v

    def callback(self, rk):
        pass

    def solve_solution_step(self):
        iu = LinearOperator((self.u, self.u), self.identity_matvec)
        iw = LinearOperator((self.w, self.w), self.identity_matvec)
        dx = self.x.copy()
        dy = self.y.copy()
        # initial value
        self.x = self.predictor.predict(self.x)
        # first coupling iteration
        self.y = yt = self.solver_wrappers[0].solve_solution_step(self.x.copy()).copy()
        self.model_f.add(self.x.copy(), yt.copy())
        xt = self.solver_wrappers[1].solve_solution_step(self.y.copy()).copy()
        r = xt - self.x
        self.model_s.add(self.y.copy(), xt)
        self.finalize_iteration(r)
        # coupling iteration loop
        while not self.convergence_criterion.is_satisfied():
            if not self.model_s.is_ready() or not self.model_f.is_ready:
                dx = self.omega * r
            else:
                mf = LinearOperator((self.w, self.u), self.lop_f)
                ms = LinearOperator((self.u, self.w), self.lop_s)
                a = iu - ms @ mf
                b = (xt - self.x).get_interface_data() + ms @ (yt - self.y).get_interface_data()
                dx_sol, exitcode = gmres(a, b, tol=self.rtol, atol=self.atol, maxiter=20, callback=self.callback)
                if exitcode != 0:
                    RuntimeError('GMRES failed')
                dx.set_interface_data(dx_sol)
            self.x += dx
            yt = self.solver_wrappers[0].solve_solution_step(self.x.copy()).copy()
            self.model_f.add(self.x.copy(), yt.copy())
            if not self.model_s.is_ready() or not self.model_f.is_ready:
                dy = yt - self.y
            else:
                a = iw - mf @ ms
                b = (yt - self.y).get_interface_data() + mf @ (xt - self.x).get_interface_data()
                dy_sol, exitcode = gmres(a, b, tol=self.rtol, atol=self.atol, maxiter=20, callback=self.callback)
                if exitcode != 0:
                    RuntimeError('GMRES failed')
                dy.set_interface_data(dy_sol)
            self.y += dy
            xt = self.solver_wrappers[1].solve_solution_step(self.y.copy()).copy()
            r = xt - self.x
            self.model_s.add(self.y.copy(), xt)
            self.finalize_iteration(r)

    def check_restart_data(self, restart_data, coupled_solver_settings=None):
        continue_check = super().check_restart_data(restart_data, ['omega', 'absolute_tolerance_gmres',
                                                                   'relative_tolerance_gmres'])
        if continue_check:
            for i, model_name in enumerate(['model_f', 'model_s']):
                old_model_type = restart_data['parameters']['settings'][model_name]['type']
                new_model_type = self.parameters['settings'][model_name]['type']
                if new_model_type != old_model_type:
                    tools.print_info(f'Model type changed from "{old_model_type}" to "{new_model_type}"', layout='blue')
                else:
                    self.restart_models[i] = True
                    self.models[i].check_restart_data(restart_data['parameters']['settings'][model_name])

    def add_restart_data(self, restart_data):
        restart_data.update({'model_f': self.model_f.save_restart_data(), 'model_s': self.model_s.save_restart_data()})
