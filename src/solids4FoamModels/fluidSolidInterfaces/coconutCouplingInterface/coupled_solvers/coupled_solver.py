from coconut import tools
from coconut.tools import create_instance
from coconut.coupling_components.component import Component

import numpy as np
import time
import pickle
import os
from datetime import datetime
import socket


class CoupledSolver(Component):
    def __init__(self, parameters):
        super().__init__()

        self.parameters = parameters
        self.settings = parameters['settings']
        self.start_init_time = time.time()  # start of initialization

        # read parameters
        self.case_name = self.settings.get('case_name', 'case')  # case name
        self.settings['case_name'] = self.case_name  # make sure a case name is present
        self.timestep_start_global = self.settings['timestep_start']  # time step for global calculation (restart)
        self.timestep_start_current = self.settings['timestep_start']  # time step start for this calculation (restart)
        self.restart = self.timestep_start_current != 0  # true if restart
        self.save_restart = self.settings.get('save_restart', -1)  # time step interval to save restart data
        self.settings['save_restart'] = self.save_restart  # in order to pass on default value
        self.save_results = self.settings.get('save_results', 0)  # time step interval to save results
        self.anonymous = self.settings.get('anonymous', False)  # disables saving 'info' in the pickle file
        self.time_step = self.timestep_start_current  # time step
        self.delta_t = self.settings['delta_t']  # time step size

        self.predictor = create_instance(self.parameters['predictor'], 'predictors.dummy_predictor')
        self.convergence_criterion = create_instance(self.parameters['convergence_criterion'],
                                                     'convergence_criteria.dummy_convergence_criterion')
        self.solver_wrappers = []
        self.index_mapped = None
        self.index_other = None
        for index in range(2):
            parameters = self.parameters['solver_wrappers'][index]
            # add timestep_start, delta_t and save_restart to solver_wrapper settings
            tools.pass_on_parameters(self.settings, parameters['settings'], ['timestep_start', 'delta_t',
                                                                             'save_restart'])
            self.solver_wrappers.append(create_instance(parameters))
            # determine index of mapped solver if present
            if self.solver_wrappers[-1].mapped:
                self.index_mapped = index
            else:
                self.index_other = index
        if self.index_other is None:
            raise ValueError('Not both solvers may be mapped solvers.')

        self.components = [self.predictor, self.convergence_criterion, self.solver_wrappers[0], self.solver_wrappers[1]]

        self.x = None  # input interface of solver 0
        self.y = None  # input interface of solver 1
        self.iteration = None  # iteration
        self.solver_level = 0  # 0 is main solver (time step is printed)
        self.init_time = None
        self.start_run_time = None
        self.run_time = None
        self.run_time_previous = 0
        self.save_time = 0
        self.time_allocation = {'previous_calculations': []}
        self.iterations = []

        # restart
        if self.restart:
            self.restart_case = self.settings.get('restart_case', self.case_name)  # case to restart from
            self.restart_data = self.load_restart_data()
            self.restart_predictor = None  # indicates if predictor has to be restarted

        # save results variables
        if self.save_results:
            self.complete_solution_x = None
            self.complete_solution_y = None
            self.residual = []
            self.info = None

        # debug
        self.debug = self.settings.get('debug', False)  # save results each iteration including residual interfaces
        if self.debug:
            self.complete_solution_r = None

    def initialize(self, print_components=True):
        super().initialize()

        self.convergence_criterion.initialize(self.solver_wrappers)

        # initialize mappers if required
        if self.index_mapped is not None:
            self.solver_wrappers[self.index_other].initialize()
            interface_input_from = self.solver_wrappers[self.index_other].get_interface_output()
            interface_output_to = self.solver_wrappers[self.index_other].get_interface_input()
            self.solver_wrappers[self.index_mapped].initialize(interface_input_from, interface_output_to)
        else:
            self.solver_wrappers[0].initialize()
            self.solver_wrappers[1].initialize()

        self.x = self.solver_wrappers[1].get_interface_output().copy()
        self.y = self.solver_wrappers[0].get_interface_output().copy()
        self.predictor.initialize(self.x)

        if self.solver_level == 0:
            title = '\n╔' + 78 * '═' + f'╗\n║{self.case_name.upper():^78}║\n╚' + 78 * '═' + '╝\n'
            tools.print_info(title)

        if print_components and self.solver_level == 0:
            self.print_components_info(' ')

        # restart
        if self.restart:
            if not (self.x.has_same_model_parts(self.restart_data['interface_x']) and
                    self.y.has_same_model_parts(self.restart_data['interface_y'])):
                raise ValueError('Restart not possible because model parts changed')
            self.check_restart_data(self.restart_data)
            if self.restart_predictor:
                self.predictor.restart(self.restart_data['predictor'])

        # update save results
        if self.save_results:
            results_data = None
            if self.restart:
                results_data = self.load_results_data()
            if results_data is None:  # no results file to append to
                self.info = f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")} : ' \
                            f'start calculation of time step {self.timestep_start_current} on {socket.gethostname()}\n'
                if self.debug:
                    self.complete_solution_x = np.empty((self.x.get_interface_data().shape[0], 0))
                    self.complete_solution_y = np.empty((self.y.get_interface_data().shape[0], 0))
                    self.complete_solution_r = np.empty((self.x.get_interface_data().shape[0], 0))
                else:
                    self.complete_solution_x = self.x.get_interface_data().reshape(-1, 1)
                    self.complete_solution_y = self.y.get_interface_data().reshape(-1, 1)

        self.start_run_time = time.time()  # start of calculation
        self.init_time = self.start_run_time - self.start_init_time  # duration of initialization

    def initialize_solution_step(self):
        super().initialize_solution_step()

        # update time step and iteration
        self.time_step += 1
        self.iteration = 0

        # print time step
        if not self.solver_level:
            self.print_header()

        for component in self.components:
            component.initialize_solution_step()

        # update save results
        if self.save_results:
            self.residual.append([])

    def solve_solution_step(self):
        pass

    def finalize_iteration(self, r):
        self.iteration += 1  # increment iteration
        self.convergence_criterion.update(r.copy())  # update convergence criterion
        self.print_iteration_info(r)  # print iteration information
        self.output_iteration(r)

    @tools.time_save
    def output_iteration(self, r):
        # update save results
        if self.save_results:
            self.residual[self.time_step - self.timestep_start_global - 1].append(r.norm())
            if self.debug:
                self.complete_solution_x = np.hstack((self.complete_solution_x,
                                                      self.x.get_interface_data().reshape(-1, 1)))
                self.complete_solution_y = np.hstack((self.complete_solution_y,
                                                      self.y.get_interface_data().reshape(-1, 1)))
                self.complete_solution_r = np.hstack((self.complete_solution_r, r.get_interface_data().reshape(-1, 1)))
                self.output_solution_step()

    def finalize_solution_step(self):
        super().finalize_solution_step()

        self.predictor.update(self.x)
        for component in self.components:
            component.finalize_solution_step()

    @tools.time_save
    def output_solution_step(self):
        super().output_solution_step()

        for component in self.components:
            component.output_solution_step()

        # save data for restart
        if self.save_restart != 0 and self.time_step % self.save_restart == 0:
            restart_data = self.save_restart_data()
            with open(self.case_name + f'_restart_ts{self.time_step}.pickle', 'wb') as file:
                pickle.dump(restart_data, file)
            if self.save_restart < 0 and self.time_step + self.save_restart > self.timestep_start_current:
                try:
                    os.remove(self.case_name + f'_restart_ts{self.time_step + self.save_restart}.pickle')
                except OSError:
                    pass

        # update save results
        self.iterations.append(self.iteration)
        if self.save_results:
            if not self.debug:
                self.complete_solution_x = np.hstack((self.complete_solution_x,
                                                      self.x.get_interface_data().reshape(-1, 1)))
                self.complete_solution_y = np.hstack((self.complete_solution_y,
                                                      self.y.get_interface_data().reshape(-1, 1)))

        # output save results
        if self.save_results != 0 and (self.time_step % self.save_results == 0 or
                                       (self.save_restart != 0 and self.time_step % self.save_restart == 0)):
            self.time_allocation.update(self.get_time_allocation())
            output = {'solution_x': self.complete_solution_x, 'solution_y': self.complete_solution_y,
                      'interface_x': self.x, 'interface_y': self.y, 'iterations': self.iterations,
                      'residual': self.residual, 'run_time': self.run_time + self.run_time_previous,
                      'time_allocation': self.time_allocation, 'delta_t': self.delta_t,
                      'timestep_start': self.timestep_start_global, 'case_name': self.case_name}
            if not self.anonymous:
                output['info'] = self.info
            if self.debug:
                output.update({'solution_r': self.complete_solution_r})
            with open(self.case_name + '_results.pickle', 'wb') as file:
                pickle.dump(output, file)

    def finalize(self):
        super().finalize()

        time_allocation = self.get_time_allocation()  # the save time for the coupling will have changed

        for component in self.components:
            component.finalize()

        # print summary
        if self.solver_level == 0:
            self.print_summary(time_allocation)

    def load_restart_data(self):
        restart_file_name = self.restart_case + f'_restart_ts{self.timestep_start_current}.pickle'
        if restart_file_name not in os.listdir(os.getcwd()):
            raise FileNotFoundError(f'Not able to perform restart because {restart_file_name} '
                                    f'not found in {os.getcwd()}')
        else:
            with open(restart_file_name, 'rb') as restart_file:
                restart_data = pickle.load(restart_file)
        return restart_data

    def check_restart_data(self, restart_data, coupled_solver_settings=None):
        # predictor type
        new_value = self.parameters['predictor']['type']
        old_value = restart_data['parameters']['predictor']['type']
        if new_value != old_value:
            tools.print_info(f'Predictor type changed from "{old_value}" to "{new_value}"', layout='blue')
        extrapolators = ('predictors.constant', 'predictors.linear', 'predictors.quadratic', 'predictors.cubic',
                         'predictors.legacy')
        if new_value in extrapolators and old_value in extrapolators:
            self.restart_predictor = True
            self.predictor.check_restart_data(restart_data['parameters']['predictor'])  # check settings
        # coupled solver type
        new_value = self.parameters['type']
        old_value = restart_data['parameters']['type']
        if new_value != old_value:
            tools.print_info(f'CoupledSolver type changed from "{old_value}" to "{new_value}"', layout='blue')
            continue_check = False
        else:
            continue_check = True
            if coupled_solver_settings is not None:
                coupled_solver_type = new_value
                for key in coupled_solver_settings:
                    new_value = self.parameters['settings'].get(key)  # None if not present
                    old_value = restart_data['parameters']['settings'].get(key)  # None if not present
                    if new_value != old_value:
                        tools.print_info(f'"{coupled_solver_type}" parameter "{key}" changed from {old_value} to '
                                         f'{new_value}', layout='blue')
        # delta_t
        if self.delta_t != restart_data['delta_t']:
            raise ValueError(f"Time step size has changed upon restart:\n\told: {restart_data['delta_t']}s"
                             f"\n\tnew: {self.delta_t}s")
        return continue_check

    @tools.time_save
    def save_restart_data(self):
        restart_data = {'predictor': self.predictor.save_restart_data(), 'interface_x': self.x, 'interface_y': self.y,
                        'parameters': {key: self.parameters[key] for key in ('type', 'settings', 'predictor')},
                        'delta_t': self.delta_t, 'time_step': self.time_step}
        self.add_restart_data(restart_data)
        return restart_data

    def add_restart_data(self, restart_data):
        pass

    def load_results_data(self):
        results_file_name = f'{self.case_name}_results.pickle'
        try:
            with open(results_file_name, 'rb') as results_file:
                results_data = pickle.load(results_file)
        except FileNotFoundError:
            tools.print_info(f'Not able to append results to {results_file_name} because file not found\n'
                             f' Saving results to new file: {results_file_name}', layout='warning')
            return
        if self.debug != ('solution_r' in results_data.keys()):
            raise ValueError(f'Value of debug attribute in {self.__class__.__name__} can not be changed upon restart')
        self.timestep_start_global = results_data['timestep_start']
        self.complete_solution_x = results_data['solution_x'][:, :self.timestep_start_current
                                                              - self.timestep_start_global + 1]
        self.complete_solution_y = results_data['solution_y'][:, :self.timestep_start_current
                                                              - self.timestep_start_global + 1]
        self.x = results_data['interface_x']
        self.y = results_data['interface_y']
        self.iterations = results_data['iterations'][:self.timestep_start_current - self.timestep_start_global]
        self.run_time_previous = results_data['run_time']
        self.residual = results_data['residual'][:self.timestep_start_current - self.timestep_start_global]
        self.time_allocation['previous_calculations'] = results_data['time_allocation'].pop('previous_calculations')
        self.time_allocation['previous_calculations'].append(results_data['time_allocation'])
        self.info = results_data.get('info', '') + '' + f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")} :' \
            f' restart calculation from time step {self.timestep_start_current} on {socket.gethostname()}\n'
        if self.debug:
            tools.print_info(f'Restart in debug mode may not append results to pickle file correctly', layout='warning')
            self.complete_solution_r = results_data['solution_r']
        return results_data

    def print_summary(self, time_allocation):
        pre = '║' + ' │' * self.solver_level

        out = '╔' + 79 * '═' + '\n║\tSummary\n╠' + 79 * '═' + '\n'
        out += f'{pre}Total calculation time{" (after restart)" if self.restart else ""}: ' \
               f'{time_allocation["total"]:.3f}s\n'

        # initialization time
        out += f'{pre}Initialization time: {time_allocation["init_time"]["total"]:0.3f}s\n'
        out += f'{pre}Distribution of initialization time:\n'
        out += self.print_time_distribution(time_allocation['init_time'], pre + '  ')

        # run time
        out += f'{pre}Run time{" (after restart)" if self.restart else ""}: ' \
               f'{time_allocation["run_time"]["total"]:0.3f}s\n'
        out += f'{pre}Distribution of run time:\n'
        out += self.print_time_distribution(time_allocation['run_time'], pre + '  ', last=False)

        # save time (part of the run time)
        reference = time_allocation['run_time']['total']
        ta_save = time_allocation['run_time']['save_time']
        out += f'{pre}  └─Save time: {ta_save["total"]:.0f}s ({ta_save["total"] / reference * 100:0.1f}%)\n'
        out += self.print_time_distribution(ta_save, pre + '    ', reference=reference)

        out += f'{pre}Average number of iterations per time step' \
               f'{" (including before restart)" if self.restart else ""}: {np.array(self.iterations).mean():0.2f}'
        out += '\n╚' + 79 * '═'
        tools.print_info(out)

    def get_time_allocation(self):
        self.run_time = time.time() - self.start_run_time  # duration of calculation
        time_allocation = {'total': self.run_time + self.init_time}
        for time_type in ('init_time', 'run_time', 'save_time'):
            time_allocation_sub = time_allocation[time_type] = {}
            time_allocation_sub['total'] = self.__getattribute__(time_type)
            for i, solver_wrapper in enumerate(self.solver_wrappers):
                time_allocation_sub[f'solver_wrapper_{i}'] = solver_wrapper.get_time_allocation()[time_type]
            time_allocation_sub['coupling'] = self.__getattribute__(time_type) - sum(
                [s.__getattribute__(time_type) for s in self.solver_wrappers])
        self.add_time_allocation(time_allocation)
        if self.solver_level == 0:
            time_allocation_save = time_allocation['run_time']['save_time'] = time_allocation.pop('save_time')
            time_allocation['run_time']['coupling'] = time_allocation['run_time']['coupling'] \
                - time_allocation_save['total']
        return time_allocation

    def add_time_allocation(self, time_allocation):
        pass

    def print_time_distribution(self, ta, pre, reference=None, last=True):
        out = ''
        if reference is None:
            reference = ta['total']
        for i, solver in enumerate(self.solver_wrappers):
            ta_solver = ta[f'solver_wrapper_{i}']
            if isinstance(ta_solver, float):
                out += f'{pre}├─{solver.__class__.__name__}: {ta_solver:.0f}s ({ta_solver / reference * 100:0.1f}%)\n'
            else:
                out += f'{pre}├─{solver.__class__.__name__}: {ta_solver["total"]:.0f}s ' \
                       f'({ta_solver["total"] / reference * 100:0.1f}%)\n'
                if solver.mapped:
                    out += f'{pre}│ └─{solver.solver_wrapper.__class__.__name__}: ' \
                           f'{ta_solver["solver_wrapper"]:.0f}s ' \
                           f'({ta_solver["solver_wrapper"] / reference * 100:0.1f}%)\n'
        out += self.add_time_distribution(ta, pre, reference)
        out += f'{pre}{"└" if last else "├"}─Coupling: {ta["coupling"]:.0f}s ' \
               f'({ta["coupling"] / reference * 100:0.1f}%)\n'
        return out

    @staticmethod
    def add_time_distribution(*_):
        return ''

    def print_header(self):
        header = (80 * '═' + f'\n\tTime step {self.time_step}\n' +
                  80 * '═' + f'\n{"Iteration":<16}{"Norm residual":<28}')
        tools.print_info(header, flush=True)

    def print_iteration_info(self, r):
        info = f'{self.iteration:<16d}{r.norm():<28.17e}'
        tools.print_info(' │' * self.solver_level, info, flush=True)

    def print_components_info(self, pre):
        tools.print_info(pre, 'The coupled solver ', self.__class__.__name__, ' has the following components:')
        tools.print_components_info(pre, self.components)

        # restart
        if self.restart:
            tools.print_info(80 * '═' + f'\n\tRestart from time step {self.timestep_start_current}\n' + 80 * '═')
