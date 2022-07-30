"""Define the structure of a non lineal finite element solver"""

import numpy as np
import copy
import logging
from tqdm import tqdm
from .Solver import Solver


class NonLinealSolver(Solver):
    """General class for non lineal solvers
    Args:
        tol (float): Tolerance for the maximum absolute value for the delta vector
        n (int): Maximum number of iterations per step
    """

    def __init__(self, FEMObject: "Core", tol: float, n: int) -> None:
        """General class for non lineal solvers
        Args:
            FEMObject (Core): FEM Object
            tol (float): Tolerance for the maximum absolute value for the delta vector
            n (int): Maximum number of iterations per step
        """
        Solver.__init__(self, FEMObject)
        self.maxiter = n
        self.tol = tol
        self.type = 'non-lineal'

    def run(self, **kargs) -> None:
        """Solves the equation system using newtons method
        """
        self.solve(**kargs)


class Newton(NonLinealSolver):
    """Creates a Newton Raphson iterative solver

        Args:
            FEMObject (Core): Finite Element Model. The model have to calculate tangent matrix T in the self.elementMatrices() method.
        """

    def __init__(self, FEMObject: 'Core', tol: float = 10**(-10), n: int = 50) -> None:
        """Creates a Newton Raphson iterative solver

        Args:
            FEMObject (Core): Finite Element Model. The model have to calculate tangent matrix T in the self.elementMatrices() method.
            tol (float, optional): Tolerance for the maximum absolute value for the delta vector. Defaults to 10**(-10).
            n (int, optional): Maximum number of iterations per step. Defaults to 50.
        """
        NonLinealSolver.__init__(self, FEMObject, tol, n)
        self.type = 'non-lineal-newton'

    def solve(self, path: str = '', **kargs) -> None:
        """Solves the equation system using newtons method
        """

        logging.info('Starting newton iterations.')
        logging.info(f'tol: {self.tol}, maxiter: {self.maxiter}')
        self.system.U = np.zeros(self.system.U.shape)+1.0
        # self.setSolution(0)
        for i in self.system.cbe:
            self.system.U[int(i[0])] = i[1]

        for e in self.system.elements:
            e.restartMatrix()
            e.setUe(self.system.U)
        warn = 'Max number of iterations. Not convergence achived!'
        for i in tqdm(range(self.maxiter), unit="Newton iteration", disable=False):
            logging.debug(
                f'----------------- Newton iteration {i} -------------------')
            self.system.restartMatrix()
            logging.debug('Matrix at 0')
            self.system.elementMatrices()
            logging.debug('Calculating element matrix')
            self.system.ensembling()
            logging.debug('Matrices enssembling')
            self.system.borderConditions()
            logging.debug('Border conditions')
            R = self.system.K@self.system.U - self.system.S
            logging.debug('Residual')
            try:
                du = -np.linalg.solve(self.system.T, R)
            except Exception as e:
                logging.error(e)
                raise e

            logging.debug('delta u')
            self.system.U += du
            for e in self.system.elements:
                e.restartMatrix()
                e.setUe(self.system.U)
            logging.debug('Updated elements')
            err = np.max(np.abs(du))
            logging.info(
                f'----------------- Iteration error {err} -------------------')
            if err < self.tol:
                warn = 'No warnings'
                break
        self.solutions = [self.system.U]
        self.solutions_info = [
            {'solver-type': self.type, 'last-it-error': err, 'n-it': i, 'warnings': warn}]
        logging.info('Done!')


class DirectIteration(NonLinealSolver):
    """docstring for DirectIteration
    """

    def __init__(self, FEMObject: 'Core', tol: float = 10**(-10), n: int = 50) -> None:
        """Creates a Direct Iteration iterative solver

        Args:
            FEMObject (Core): Finite Element Model. The model have to calculate tangent matrix T in the self.elementMatrices() method.
            tol (float, optional): Tolerance for the maximum absolute value for the delta vector. Defaults to 10**(-10).
            n (int, optional): Maximum number of iterations per step. Defaults to 50.
        """
        NonLinealSolver.__init__(self, FEMObject, tol, n)
        self.type = 'non-lineal-direct'

    def solve(self, path: str = '', guess=None, _guess=False, **kargs) -> None:
        """Solves the equation system using newtons method
        """

        logging.info('Starting iterations.')
        logging.info(f'tol: {self.tol}, maxiter: {self.maxiter}')
        if _guess:
            self.system.U = guess
        else:
            self.system.U = np.zeros(self.system.U.shape)

        for i in self.system.cbe:
            self.system.U[int(i[0])] = i[1]

        for e in self.system.elements:
            e.restartMatrix()
            e.setUe(self.system.U)
        warn = 'Max number of iterations. Not convergence achived!'
        for i in tqdm(range(self.maxiter), unit="Iteration", disable=False):
            logging.debug(
                f'----------------- Iteration {i+1} -------------------')
            self.system.restartMatrix()
            logging.debug('Matrix at 0')
            self.system.elementMatrices()
            logging.debug('Calculating element matrix')
            self.system.ensembling()
            logging.debug('Matrices enssembling')
            self.system.borderConditions()
            logging.debug('Border conditions')
            uim11 = self.system.U.copy()
            try:
                self.system.U = np.linalg.solve(self.system.K, self.system.S)
            except Exception as e:
                logging.error(e)
                raise e
            logging.debug('Equation system solved')

            R = (self.system.U - uim11)
            logging.debug('Residual')

            for e in self.system.elements:
                e.restartMatrix()
                e.setUe(self.system.U)
            logging.debug('Updated elements')
            err = np.max(np.abs(R))
            logging.info(
                f'----------------- Iteration error {err} -------------------')
            if err < self.tol:
                warn = 'No warnings'
                break
        self.solutions_info = [
            {'solver-type': self.type, 'last-it-error': err, 'n-it': i, 'warnings': warn}]
        self.solutions = [self.system.U]
        logging.info('Done!')


class LoadControl(DirectIteration):
    """General class for non lineal solvers
    Args:
        tol (float): Tolerance for the maximum absolute value for the delta vector
        n (int): Maximum number of iterations per step
    """

    def __init__(self, FEMObject, tol: float = 10**(-10), n: int = 500, nls=10) -> None:
        """General class for non lineal solvers
        Args:
            tol (float): Tolerance for the maximum absolute value for the delta vector
            n (int): Maximum number of iterations per step
            nls (int): Number of load steps
        """
        DirectIteration.__init__(self, FEMObject, tol=tol, n=n)
        self.nls = nls
        self.type += '-load-control'

    def run(self, **kargs) -> None:
        """Solves the equation system using newtons method
        """
        guess = None

        solutioms = []
        solutioms_info = []
        for i in tqdm(range(self.nls), unit="Load Step", disable=False):
            logging.info(f'================LOAD STEP {i+1}===================')
            # FIXME WTF IS THIS. Esto solamente funciona para la clase de EB no lineal
            self.system.fx = lambda x: self.system.fx0(x)/self.nls*(i+1)
            self.system.fy = lambda x: self.system.fy0(x)/self.nls*(i+1)
            guess = self.system.U
            self.solve(guess=guess, _guess=(i >= 1), **kargs)
            solutioms.append(self.system.U)
            solutioms_info.append(copy.deepcopy(self.solutions_info[-1]))
        self.solutions_info = solutioms_info
        self.solutions = solutioms
        self.setSolution()
