"""Define the structure of a lineal finite element solver
"""

import numpy as np
import logging
from scipy.sparse.linalg import spsolve
from .Solver import Solver
from tqdm import tqdm


class Parabolic(Solver):
    """Lineal Finite Element Solver.
    """

    def __init__(self, FEMObject: 'Core'):
        """Lineal Finite Element Solver

        Args:
            FEMObject (Core): Finite Element Problem
        """
        Solver.__init__(self, FEMObject)
        self.type = 'lineal-transient'
        self.solutions = []

    def run(self, t0, tf, steps, dt=None):
        if not dt:
            self.system.dt = (tf-t0)/steps
            self.system.t = t0
        else:
            self.system.dt = dt
            self.system.t = t0
        self.system.apply_initial_condition(t0, dt)
        for _ in tqdm(range(steps)):
            logging.info(
                f'Solving for time {self.system.t}. Using dt = {self.system.dt}')
            logging.info('Creating element matrices...')
            self.system.restartMatrix()
            self.system.elementMatrices()
            self.system.ensembling()
            self.system.boundaryConditions()
            logging.info('Solving equation system...')
            self.solutions.append(np.linalg.solve(
                self.system.K, self.system.S))
            self.system.t += self.system.dt
            self.solutions_info.append(
                {'solver-type': self.type, "time": self.system.t, "dt": self.system.dt})
            self.setSolution(elements=True)
        logging.info('Done!')
