"""Define the structure of a non lineal finite element solver
"""

import numpy as np
import logging
from tqdm import tqdm


class NonLinealSolver():
    """General class for non lineal solvers
    """

    def __init__(self):
        """General class for non lineal solvers
        """
        self.type = 'non-lineal'


class Newton(NonLinealSolver):
    """Creates a Newton Raphson iterative solver

        Args:
            FEMObject (Core): Finite Element Model. The model have to calculate tangent matrix T in the self.elementMatrices() method.
        """

    def __init__(self, FEMObject: 'Core', tol: float = 10**(-10), n: int = 50):
        """Creates a Newton Raphson iterative solver

        Args:
            FEMObject (Core): Finite Element Model. The model have to calculate tangent matrix T in the self.elementMatrices() method.
            tol (float, optional): Tolerance for the maximum absolute value for the delta vector. Defaults to 10**(-10).
            n (int, optional): Maximum number of iterations per step. Defaults to 50.
        """
        self.system = FEMObject
        self.maxiter = n
        self.tol = tol
        NonLinealSolver.__init__(self)
        self.type = 'non-lineal-newton'

    def run(self, **kargs):
        """Solves the equation system using newtons method
        """
        self.solve(**kargs)

    def solve(self, path: str = '', **kargs):
        """Solves the equation system using newtons method
        """

        logging.info('Starting newton iterations.')
        logging.info(f'tol: {self.tol}, maxiter: {self.maxiter}')
        self.system.U = np.zeros(self.system.U.shape)+1.0
        for i in self.system.cbe:
            self.system.U[int(i[0])] = i[1]

        for e in self.system.elements:
            e.restartMatrix()
            e.setUe(self.system.U)

        for i in tqdm(range(self.maxiter), disable=False):
            logging.debug(
                f'----------------- Newton iteration {i} -------------------')
            self.system.restartMatrix()
            logging.debug(f'Matrix at 0')
            self.system.elementMatrices()
            logging.debug(f'Calculating element matrix')
            self.system.ensembling()
            logging.debug(f'Matrices enssembling')
            self.system.borderConditions()
            logging.debug(f'Border conditions')
            R = self.system.K@self.system.U - self.system.S
            logging.debug(f'Residual')
            try:
                du = -np.linalg.solve(self.system.T, R)
            except Exception as e:
                logging.error(e)
                raise e

            logging.debug(f'delta u')
            self.system.U += du
            for e in self.system.elements:
                e.restartMatrix()
                e.setUe(self.system.U)
            logging.debug(f'Updated elements')
            err = np.max(np.abs(du))
            logging.info(
                f'----------------- Iteration error {err} -------------------')
            if err < self.tol:
                break
        logging.info('Done!')
