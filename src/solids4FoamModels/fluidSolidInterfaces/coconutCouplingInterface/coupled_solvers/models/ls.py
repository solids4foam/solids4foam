from coconut.coupling_components.component import Component
from coconut import tools

import numpy as np
from scipy.linalg import solve_triangular


def create(parameters):
    return ModelLS(parameters)


class ModelLS(Component):
    def __init__(self, parameters):
        super().__init__()

        self.settings = parameters['settings']
        self.min_significant = self.settings['min_significant']
        self.q = self.settings['q']

        self.size_in = None  # set by coupled solver
        self.size_out = None  # set by coupled solver
        self.out = None  # interface of output, set by coupled solver
        self.added = False
        self.rref = None
        self.xtref = None
        self.vcurr = None
        self.wcurr = None
        self.vprev = None
        self.wprev = None

    def initialize(self):
        super().initialize()

        self.vprev = [np.empty((self.size_in, 0))]
        self.wprev = [np.empty((self.size_out, 0))]

    def filter(self):
        v = np.hstack((self.vcurr, np.hstack(self.vprev)))
        if not v.shape[1]:
            raise RuntimeError('No information to filter')
        # remove columns resulting in small diagonal elements in R
        singular = True
        while singular and v.shape[1]:
            rr = np.linalg.qr(v, mode='r')
            diag = np.diagonal(rr)
            m = min(abs(diag))
            if m < self.min_significant:
                i = np.argmin(abs(diag))
                tools.print_info(f'Removing column {i}: {m} < min_significant', layout='info')
                if i < self.vcurr.shape[1]:
                    self.vcurr = np.delete(self.vcurr, i, 1)
                    self.wcurr = np.delete(self.wcurr, i, 1)
                else:
                    num_columns = self.vcurr.shape[1]
                    j = -1
                    while i >= num_columns:
                        j += 1
                        num_columns += self.vprev[j].shape[1]
                    num_columns -= self.vprev[j].shape[1]
                    self.vprev[j] = np.delete(self.vprev[j], i - num_columns, 1)
                    self.wprev[j] = np.delete(self.wprev[j], i - num_columns, 1)
                v = np.hstack((self.vcurr, np.hstack(self.vprev)))
            else:
                singular = False
        # remove columns if number of columns exceeds number of rows
        while v.shape[0] < v.shape[1]:
            if self.vcurr.shape[0] < self.vcurr.shape[1]:
                self.vcurr = np.delete(self.vcurr, -1, 1)
                self.wcurr = np.delete(self.wcurr, -1, 1)
            else:
                i = -1
                while self.vprev[i].shape[1] == 0:
                    i -= 1
                self.vprev[i] = np.delete(self.vprev[i], -1, 1)
                self.wprev[i] = np.delete(self.wprev[i], -1, 1)
                tools.print_info(f'Removing column {i}: too many columns', layout='warning')
            v = np.hstack((self.vcurr, np.hstack(self.vprev)))

    def predict(self, dr_in, modes=None):
        dr = dr_in.get_interface_data().reshape(-1, 1)
        self.filter()
        v = np.hstack((limit(self.vcurr, modes), np.hstack([limit(v, modes) for v in self.vprev])))
        w = np.hstack((limit(self.wcurr, modes), np.hstack([limit(w, modes) for w in self.wprev])))
        if not v.shape[1]:
            tools.print_info('Least-squares model has no information to predict: zero returned', layout='warning')
            return dr_in * 0
        # approximation for the inverse of the Jacobian from a least-squares model
        qq, rr = np.linalg.qr(v, mode='reduced')
        b = qq.T @ dr
        c = solve_triangular(rr, b)
        dxt = w @ c
        dxt_out = self.out.copy()
        dxt_out.set_interface_data(dxt.flatten())
        return dxt_out

    def add(self, r_in, xt_in):
        r = r_in.get_interface_data().reshape(-1, 1)
        xt = xt_in.get_interface_data().reshape(-1, 1)
        if self.added:
            dr = r - self.rref
            dxt = xt - self.xtref
            # Update V and W matrices
            self.vcurr = np.hstack((dr, self.vcurr))
            self.wcurr = np.hstack((dxt, self.wcurr))
        else:
            self.added = True
        self.rref = r
        self.xtref = xt

    def is_ready(self):
        v = np.hstack((self.vcurr, np.hstack(self.vprev)))
        return v.shape[1]

    def initialize_solution_step(self):
        super().initialize_solution_step()

        self.rref = None
        self.xtref = None
        self.vcurr = np.empty((self.size_in, 0))
        self.wcurr = np.empty((self.size_out, 0))
        self.added = False

        while len(self.vprev) > max(self.q, 1):  # for q = 0, there is one item: an empty array
            self.vprev.pop()
            self.wprev.pop()

    def finalize_solution_step(self):
        super().finalize_solution_step()

        if self.q > 0:
            self.vprev = [self.vcurr] + self.vprev
            self.wprev = [self.wcurr] + self.wprev

    def filter_q(self, dr_in, modes=None):
        dr = dr_in.get_interface_data().reshape(-1, 1)
        dr_out = dr_in.copy()
        v = np.hstack((limit(self.vcurr, modes), np.hstack([limit(v, modes) for v in self.vprev])))
        if v.shape[1]:
            qq, *_ = np.linalg.qr(v, mode='reduced')
            dr = dr - qq @ (qq.T @ dr)
            dr_out.set_interface_data(dr.flatten())
        return dr_out

    def restart(self, restart_data):
        if self.q != 0:
            self.vprev = restart_data['vprev']
            self.wprev = restart_data['wprev']

    def check_restart_data(self, restart_data):
        model_type = restart_data['type']
        for key in ['min_significant', 'q']:
            new_value = self.settings.get(key)
            old_value = restart_data['settings'].get(key)
            if new_value != old_value:
                tools.print_info(f'"{model_type}" parameter "{key}" changed from {old_value} to {new_value}',
                                 layout='blue')

    def save_restart_data(self):
        return {'settings': self.settings, 'vprev': self.vprev, 'wprev': self.wprev}


def limit(matrix, modes):
    if modes is None:
        return matrix
    elif modes == 0:
        return matrix[:, :0]
    return matrix[:, -min(matrix.shape[1], modes):]  # return oldest
    # return matrix[:, :min(matrix.shape[1], modes)]  # return most recent
