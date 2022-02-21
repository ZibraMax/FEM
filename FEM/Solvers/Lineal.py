"""Define the structure of a lineal finite element solver
"""

import numpy as np
import logging
from scipy.sparse.linalg import spsolve


class Lineal():
    """Lineal Finite Element Solver.
    """

    def __init__(self, FEMObject: 'Core'):
        """Lineal Finite Element Solver

        Args:
            FEMObject (Core): [Finite Element Problem
        """
        self.system = FEMObject
        self.type = 'lineal'

    def run(self, path: str = '', **kargs):
        """Solves the equation system using numpy's solve function

        Args:
            path (str, optional): Path where the solution is stored. Defaults to ''.
        """
        logging.info('Solving equation system...')
        self.system.U = np.linalg.solve(self.system.K, self.system.S)
        if not path == '':
            np.savetxt(path, self.system.U, delimiter=',')
        for e in self.system.elements:
            e.setUe(self.system.U)
        logging.info('Done!')


class LinealSparse():
    """Lineal Finite Element Solver using sparse matrix
    """

    def __init__(self, FEMObject: 'Core'):
        """Lineal Finite Element Solver

        Args:
            FEMObject (Core): [Finite Element Problem
        """
        self.system = FEMObject
        self.type = 'lineal-sparse'

    def run(self, path: str = '', **kargs):
        """Solves the equation system using numpy's solve function

        Args:
            path (str, optional): Path where the solution is stored. Defaults to ''.
        """
        logging.info('Converting to csr format')
        self.system.K = self.system.K.tocsr()
        logging.info('Solving...')
        self.system.U = spsolve(self.system.K, self.system.S)
        for e in self.system.elements:
            e.setUe(self.system.U)
        logging.info('Solved!')
