"""Define the structure of a lineal finite element solver
"""

from scipy.linalg import eigh
from scipy.sparse.linalg import eigsh
import numpy as np
import logging
from scipy.sparse.linalg import spsolve
from .Solver import Solver


class Lineal(Solver):
    """Lineal Finite Element Solver.
    """

    def __init__(self, FEMObject: 'Core'):
        """Lineal Finite Element Solver

        Args:
            FEMObject (Core): Finite Element Problem
        """
        Solver.__init__(self, FEMObject)
        self.type = 'lineal'
        self.solutions = []

    def run(self, path: str = '', **kargs):
        """Solves the equation system using numpy's solve function

        Args:
            path (str, optional): Path where the solution is stored. Defaults to ''.
        """
        # TODO This should delete all matrices becaus the add feature
        logging.info('Creating element matrices...')
        self.system.elementMatrices()
        logging.info('Done!')
        self.system.ensembling()
        self.system.borderConditions()
        logging.info('Solving equation system...')
        self.solutions = [np.linalg.solve(self.system.K, self.system.S)]
        self.solutions_info = [{'solver-type': self.type}]
        self.setSolution()
        if not path == '':
            np.savetxt(path, self.system.U, delimiter=',')
        for e in self.system.elements:
            e.setUe(self.system.U)
        logging.info('Done!')


class LinealSparse(Lineal):
    """Lineal Finite Element Solver using sparse matrix
    """

    def __init__(self, FEMObject: 'Core'):
        """Lineal Finite Element Solver

        Args:
            FEMObject (Core): Finite Element Problem
        """
        Lineal.__init__(self, FEMObject)
        self.type = 'lineal-sparse'

    def run(self, path: str = '', **kargs):
        """Solves the equation system using scipy's spsolve function

        Args:
            path (str, optional): Path where the solution is stored. Defaults to ''.
        """
        logging.info('Creating element matrices...')
        self.system.elementMatrices()
        logging.info('Done!')
        self.system.ensembling()
        self.system.borderConditions()
        logging.info('Converting to csr format')
        self.system.K = self.system.K.tocsr()
        logging.info('Solving...')
        self.solutions = [spsolve(self.system.K, self.system.S)]
        self.solutions_info = [{'solver-type': self.type}]
        self.setSolution()
        if path:
            np.savetxt(path, self.system.U, delimiter=',')
        for e in self.system.elements:
            e.setUe(self.system.U)
        logging.info('Solved!')


class LinealEigen(Lineal):
    """Eigen value solver

    Args:
        FEMObject (Core): FEM problem
    """

    def __init__(self, FEMObject: 'Core'):
        """Eigen value solver

        Args:
            FEMObject (Core): FEM problem
        """
        Lineal.__init__(self, FEMObject)
        self.type = 'lineal-sparse-eigen'

    def run(self, path: str = '', k=20, **kargs):
        """Solves the smallest k eigenvalues using scipy's eigen value solver

        Args:
            path (str, optional): Path where the solution is stored. Defaults to ''.
            k (int, optional): Number of eigenvalues to calculate. Defaults to 20.
        """
        logging.info('Creating element matrices...')
        self.system.elementMatrices()
        logging.info('Done!')
        self.system.ensembling()
        self.system.condensedSystem()
        logging.info('Converting to csr format')
        K = self.system.K.tocsr()
        logging.info('Solving...')
        # eigv, eigvec = largest_eigsh(
        #     self.system.K, k, self.system.M, which='SM')
        # N = self.system.K.shape[0]
        # eigv, eigvec = eigh(
        #     self.system.K.todense(), self.system.M.todense(), eigvals=(N-k, N-1))
        eigv, eigvec = eigsh(
            K, k, self.system.M, which='SM')
        idx = eigv.argsort()
        eigv = eigv[idx]
        eigvec = eigvec[:, idx]
        self.system.eigv = eigv
        self.system.eigvec = eigvec
        if path:
            np.savetxt(path.replace('.', '_eigv.'),
                       self.system.eigv, delimiter=',', fmt='%s')
            np.savetxt(path.replace('.', '_eigvec.'),
                       self.system.eigvec, delimiter=',', fmt='%s')
        eeevalues = []
        for eigenvalue in eigvec.T:
            eeevalues.append(eigenvalue)

        self.solutions = np.array(eeevalues)
        self.solutions_info = [
            {'solver-type': self.type, 'eigv': ei} for ei in self.system.eigv]

        self.setSolution(0)
        logging.info('Solved!')
