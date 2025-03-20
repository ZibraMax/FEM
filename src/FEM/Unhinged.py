"""Origami Bar and Hinge Model
"""

from .Core import Core, Geometry, logging
from .Solvers import LinealSparse
from typing import Callable, Tuple
import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm
from scipy import sparse  # Creation of sparse matrices


class Truss(Core):
    """Truss analysis"""

    def __init__(self, geometry: Geometry, E: list, A: list, **kargs):
        """Truss analysis"""
        solver = LinealSparse  # Solver for the truss analysis
        Core.__init__(self, geometry, solver, sparse=True, **kargs)
        if isinstance(E, float) or isinstance(E, int):
            E = [E]*len(geometry.elements)
        if isinstance(A, float) or isinstance(A, int):
            A = [A]*len(geometry.elements)
        self.E = E
        self.A = A
        self.properties["E"] = E
        self.properties["A"] = A
        self.name = 'Truss analysis'
        self.K = sparse.lil_matrix((self.ngdl, self.ngdl))

    def addLoadNode(self, node: int, load: list):
        """Add a node load to the truss analysis"""
        self.cbn += [[node*3, load[0]],
                     [node*3+1, load[1]], [node*3+2, load[2]]]

    def elementMatrices(self):
        for ee, e in enumerate(tqdm(self.elements, unit='Element')):
            e.gdlm = e.gdl.T.flatten()
            Ke = np.zeros([6, 6])
            # Bar element matrix
            L = np.linalg.norm(e.coords[1] - e.coords[0])
            r = 2*e.jacs[0][0]/L

            K = np.array([[1, -1],
                          [-1, 1]]) * self.E[ee] * self.A[ee] / L
            o = [0.0, 0.0, 0.0]
            T = np.array([[*r, *o], [*o, *r]])
            Ke = T.T @ K @ T
            self.K[np.ix_(e.gdlm, e.gdlm)] += Ke

    def ensembling(self) -> None:
        """Creation of the system sparse matrix. Force vector is ensembled in integration method
        """
        logging.info('Ensembling equation system...')
        logging.info('Done!')


class BarAndHingeLineal(Core):
    def __init__(self, geometry: Geometry, E: list, A: list, k: Callable, **kargs):
        """Truss analysis"""
        solver = LinealSparse
        Core.__init__(self, geometry, solver, sparse=True, **kargs)
        if isinstance(E, float) or isinstance(E, int):
            E = [E]*len(geometry.elements)
        if isinstance(A, float) or isinstance(A, int):
            A = [A]*len(geometry.elements)
        self.E = E
        for e in geometry.elements:
            if e.__class__.__name__ == 'OriHinge':
                e.set_kf(k)
        self.A = A
        self.properties["E"] = E
        self.properties["A"] = A
        self.name = 'Truss analysis'
        self.K = sparse.lil_matrix((self.ngdl, self.ngdl))

    def addLoadNode(self, node: int, load: list):
        """Add a node load to the truss analysis"""
        self.cbn += [[node*3, load[0]],
                     [node*3+1, load[1]], [node*3+2, load[2]]]

    def barElementMatrix(self, e, ee):
        L = np.linalg.norm(e.coords[1] - e.coords[0])
        r = 2*e.jacs[0][0]/L
        K = np.array([[1, -1],
                      [-1, 1]]) * self.E[ee] * self.A[ee] / L
        o = [0.0, 0.0, 0.0]
        T = np.array([[*r, *o], [*o, *r]])
        Ke = T.T @ K @ T
        return Ke

    def elementMatrices(self):
        for ee, e in enumerate(tqdm(self.elements, unit='Element')):
            e.gdlm = e.gdl.T.flatten()  # Change in numbering...
            if e.__class__.__name__ == 'OriHinge':
                Ke = e.elementMatrix()
            else:
                Ke = self.barElementMatrix(e, ee)
            self.K[np.ix_(e.gdlm, e.gdlm)] += Ke

    def ensembling(self) -> None:
        """Creation of the system sparse matrix. Force vector is ensembled in integration method
        """
        logging.info('Ensembling equation system...')
        logging.info('Done!')
