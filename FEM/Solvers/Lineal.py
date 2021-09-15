"""Define the structure of a lineal finite element solver
"""

import numpy as np
import logging


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
