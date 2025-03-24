"""Origami Bar and Hinge Model
"""

from .Core import Core, Geometry, logging
from .Solvers import LinearSparse, DirectIteration, Newton, MGDCM, LoadControl
from typing import Callable, Tuple
import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm
from scipy import sparse  # Creation of sparse matrices


class TrussLinear(Core):
    """Truss analysis"""

    def __init__(self, geometry: Geometry, E: list, A: list, **kargs):
        """Truss analysis"""
        solver = LinearSparse  # Solver for the truss analysis
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


class TrussNonLinear(Core):

    def __init__(self, geometry: Geometry, E: list, A: list, **kargs):
        solver = MGDCM
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
        self.cbn += [[node*3, load[0]],
                     [node*3+1, load[1]], [node*3+2, load[2]]]

    def elementMatrices(self):
        self.F_int = np.zeros([self.ngdl, 1])
        for ee, e in enumerate(tqdm(self.elements, unit='Element')):
            e.gdlm = e.gdl.T.flatten()
            L = np.linalg.norm(e.coords[1] - e.coords[0])
            C0 = self.E[ee]
            A = self.A[ee]
            o = [0, 0, 0.0]

            Ue = e.Ue.T.flatten().reshape([6, 1])
            # Undeformed shape
            B1 = np.array([[-1, 0, 0, 1.0, 0, 0]])/L
            B1P = 1/L * \
                np.array([[*-(e.coords[1] - e.coords[0])/L,
                         *(e.coords[1] - e.coords[0])/L]])
            B2 = np.array([[1, 0, 0, -1, 0, 0],
                           [0, 1, 0, 0, -1, 0],
                           [0, 0, 1, 0, 0, -1],
                           [-1, 0, 0, 1, 0, 0],
                           [0, -1, 0, 0, 1, 0],
                           [0, 0, -1, 0, 0, 1]])/(L**2)
            B2Ue = B2@Ue

            Ex = (B1P@Ue + 1/2 * Ue.T@B2Ue)[0, 0]
            S11 = C0*Ex
            Tbar = S11*A*L*(B1P.T+B2Ue)
            C = C0
            # Deformed shape

            Ke = np.zeros([6, 6])

            KE = C*A*L*B1P.T@B1P
            K1 = C*A*L*((B2Ue)@B1P + B1P.T@(B2Ue.T))
            K2 = C*A*L*B2Ue@B2Ue.T
            KG = S11*A*L*B2

            Ke = KE + K1 + K2 + KG
            self.K[np.ix_(e.gdlm, e.gdlm)] += Ke
            self.F_int[np.ix_(e.gdlm)] += Tbar

    def ensembling(self) -> None:
        logging.info('Ensembling equation system...')
        logging.info('Done!')


class BarAndHingeLinear(Core):
    def __init__(self, geometry: Geometry, E: list, A: list, k: Callable, **kargs):
        """Truss analysis"""
        solver = LinearSparse
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
                Ke, te = e.elementMatrix()
            else:
                Ke = self.barElementMatrix(e, ee)
            self.K[np.ix_(e.gdlm, e.gdlm)] += Ke

    def ensembling(self) -> None:
        """Creation of the system sparse matrix. Force vector is ensembled in integration method
        """
        logging.info('Ensembling equation system...')
        logging.info('Done!')


class BarAndHingeNonLinear(Core):

    def __init__(self, geometry: Geometry, E: list, A: list, **kargs):
        """Truss analysis"""
        solver = MGDCM
        Core.__init__(self, geometry, solver, sparse=True, **kargs)
        if isinstance(E, float) or isinstance(E, int):
            E = [E]*len(geometry.elements)
        if isinstance(A, float) or isinstance(A, int):
            A = [A]*len(geometry.elements)
        self.E = E
        self.A = A
        self.properties["E"] = E
        self.properties["A"] = A
        self.name = 'Bar and hinge analysis for origami. Non linear'
        self.K = sparse.lil_matrix((self.ngdl, self.ngdl))

    def addLoadNode(self, node: int, load: list):
        """Add a node load to the truss analysis"""
        self.cbn += [[node*3, load[0]],
                     [node*3+1, load[1]],
                     [node*3+2, load[2]]]

    def Ogden(self, Ex, C0):
        # Ogden hyperelastic constitutive model for bar elements
        alfa = [3, 1]  # Linear

        pstr = np.real(np.sqrt(2 * Ex + 1))
        C0 = np.where(pstr < 1, 1 * C0, C0)

        Ct = C0 / (alfa[0] - alfa[1]) * ((alfa[0] - 2) * pstr **
                                         (alfa[0] - 4) - (alfa[1] - 2) * pstr**(alfa[1] - 4))
        Sx = C0 / (alfa[0] - alfa[1]) * \
            (pstr**(alfa[0] - 2) - pstr**(alfa[1] - 2))

        Wb = C0 / (alfa[0] - alfa[1]) * ((pstr**alfa[0] - 1) /
                                         alfa[0] - (pstr**alfa[1] - 1) / alfa[1])

        return Sx, Ct, Wb

    def barElementMatrix(self, e, ee):
        L = np.linalg.norm(e.coords[1] - e.coords[0])
        C0 = self.E[ee]
        A = self.A[ee]
        o = [0, 0, 0.0]

        Ue = e.Ue.T.flatten().reshape([6, 1])
        # Undeformed shape
        B1 = np.array([[-1, 0, 0, 1.0, 0, 0]])/L
        B1P = 1/L * \
            np.array([[*-(e.coords[1] - e.coords[0])/L,
                       *(e.coords[1] - e.coords[0])/L]])
        B2 = np.array([[1, 0, 0, -1, 0, 0],
                       [0, 1, 0, 0, -1, 0],
                       [0, 0, 1, 0, 0, -1],
                       [-1, 0, 0, 1, 0, 0],
                       [0, -1, 0, 0, 1, 0],
                       [0, 0, -1, 0, 0, 1]])/(L**2)
        B2Ue = B2@Ue

        Ex = (B1P@Ue + 1/2 * Ue.T@B2Ue)[0, 0]
        # S11 = C0*Ex
        S11, C0, Wb = self.Ogden(Ex, C0)

        Tbar = S11*A*L*(B1P.T+B2Ue)
        C = C0
        # Deformed shape

        Ke = np.zeros([6, 6])

        KE = C*A*L*B1P.T@B1P  # This one works
        K1 = C*A*L*((B2Ue)@B1P + B1P.T@(B2Ue.T))
        K2 = C*A*L*B2Ue@B2Ue.T  # This one wordks!
        KG = S11*A*L*B2

        Ke = KE + K1 + K2 + KG
        return Ke, Tbar

    def elementMatrices(self):
        self.F_int = np.zeros([self.ngdl, 1])
        self.F = np.zeros([self.ngdl, 1])
        for ee, e in enumerate(tqdm(self.elements, unit='Element')):
            e.gdlm = e.gdl.T.flatten()  # Change in numbering...
            if e.__class__.__name__ == 'OriHinge':
                Ke, Te = e.elementMatrixNonLineal()
                if e.internal_force is not None:
                    self.F[np.ix_(e.gdlm)] += e.get_internal_force()
            else:
                Ke, Te = self.barElementMatrix(e, ee)
            self.K[np.ix_(e.gdlm, e.gdlm)] += Ke
            self.F_int[np.ix_(e.gdlm)] += Te

    def ensembling(self) -> None:
        """Creation of the system sparse matrix. Force vector is ensembled in integration method
        """
        logging.info('Ensembling equation system...')
        logging.info('Done!')
