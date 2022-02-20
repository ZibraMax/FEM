"""3D Isotropic elasticity"""

from typing import Callable, Tuple
import numpy as np
from tqdm import tqdm
from scipy import sparse
from scipy.sparse.linalg import spsolve
from .Elements.E2D import Quadrilateral
from .Core import Core, Geometry, logging


class Elasticity(Core):
    """Creates a 3D Elasticity problem

    Args:
        geometry (Geometry): Input 3D geometry
        E (Tuple[float, list]): Young Moduli
        v (Tuple[float, list]): Poisson coeficient
        rho (Tuple[float, list]): Density.
        fx (Callable, optional): Force in x direction. Defaults to lambdax:0.
        fy (Callable, optional): Force in y direction. Defaults to lambdax:0.
        fz (Callable, optional): Force in z direction. Defaults to lambdax:0.
    """

    def __init__(self, geometry: Geometry, E: Tuple[float, list], v: Tuple[float, list], rho: Tuple[float, list], fx: Callable = lambda x: 0, fy: Callable = lambda x: 0, fz: Callable = lambda x: 0, **kargs) -> None:
        """Creates a 3D Elasticity problem

        Args:
            geometry (Geometry): Input 3D geometry
            E (Tuple[float, list]): Young Moduli
            v (Tuple[float, list]): Poisson coeficient
            rho (Tuple[float, list]): Density.
            fx (Callable, optional): Force in x direction. Defaults to lambdax:0.
            fy (Callable, optional): Force in y direction. Defaults to lambdax:0.
            fz (Callable, optional): Force in z direction. Defaults to lambdax:0.
        """
        if type(E) == float or type(E) == int:
            E = [E]*len(geometry.elements)
        if type(v) == float or type(v) == int:
            v = [v]*len(geometry.elements)
        if type(rho) == float or type(rho) == int:
            rho = [rho]*len(geometry.elements)
        self.E = E
        self.v = v
        self.rho = rho
        self.fx = fx
        self.fy = fy
        self.fz = fz
        if not geometry.nvn == 3:
            print(
                'Border conditions lost, please usea a geometry with 3 variables per node (nvn=3)\nRegenerating Geoemtry...')
            geometry.nvn = 3
            geometry.cbe = []
            geometry.cbn = []
            geometry.initialize()
        Core.__init__(self, geometry, sparse=True, **kargs)

        self.I = []
        self.Im = []
        self.J = []
        self.Jm = []
        self.V = []
        self.Vm = []

    def elementMatrices(self) -> None:
        """Calculate the element matrices usign Reddy's (2005) finite element model
        """

        for ee, e in enumerate(tqdm(self.elements, unit='Element')):
            m = len(e.gdl.T)

            _x, _p = e.T(e.Z.T)
            jac, dpz = e.J(e.Z.T)
            detjac = np.linalg.det(jac)
            _j = np.linalg.inv(jac)
            dpx = _j @ dpz

            E = self.E[ee]
            v = self.v[ee]

            C = E/((1.0+v)*(1.0-2.0*v))*np.array([
                [1.0-v, v, v, 0.0, 0.0, 0.0],
                [v, 1.0-v, v, 0.0, 0.0, 0.0],
                [v, v, 1.0-v, 0.0, 0.0, 0.0],
                [0.0, 0.0, 0.0, (1.0-2.0*v)/2.0, 0.0, 0.0],
                [0.0, 0.0, 0.0, 0.0, (1.0-2.0*v)/2.0, 0.0],
                [0.0, 0.0, 0.0, 0.0, 0.0, (1.0-2.0*v)/2.0]])

            o = [0.0]*m
            Ke = np.zeros([3*m, 3*m])
            Me = np.zeros([3*m, 3*m])
            Fe = np.zeros([3*m, 1])

            for k in range(len(e.Z)):  # Iterate over gauss points on domain
                B = np.array([
                    [*dpx[k, 0, :], *o, *o],
                    [*o, *dpx[k, 1, :], *o],
                    [*o, *o, *dpx[k, 2, :]],
                    [*dpx[k, 2, :], *o, *dpx[k, 0, :]],
                    [*o, *dpx[k, 2, :], *dpx[k, 1, :]],
                    [*dpx[k, 1, :], *dpx[k, 0, :], *o]])
                P = np.array([
                    [*_p[k], *o, *o],
                    [*o, *_p[k], *o],
                    [*o, *o, *_p[k]]])

                Ke += (B.T@C@B)*detjac[k]*e.W[k]
                Me += self.rho[ee]*(P.T@P)*detjac[k]*e.W[k]

                _fx = self.fx(_x[k])
                _fy = self.fy(_x[k])
                _fz = self.fz(_x[k])

                F = np.array([[_fx], [_fy], [_fz]])
                Fe += (P.T @ F)*detjac[k]*e.W[k]

            self.F[np.ix_(e.gdlm)] += Fe

            for gdl in e.gdlm:
                self.I += [gdl]*(3*m)
                self.J += e.gdlm
            self.V += Ke.flatten().tolist()
            self.Vm += Me.flatten().tolist()

    def ensembling(self) -> None:
        """Creation of the system sparse matrix. Force vector is ensembled in integration method
        """
        logging.info('Ensembling equation system...')
        self.K = sparse.coo_matrix(
            (self.V, (self.I, self.J)), shape=(self.ngdl, self.ngdl)).tolil()
        self.M = sparse.coo_matrix(
            (self.Vm, (self.I, self.J)), shape=(self.ngdl, self.ngdl)).tocsr()
        logging.info('Done!')

    def solveES(self, **kargs) -> None:
        """Solve the sparse system
        """
        logging.info('Converting to csr format')
        self.K = self.K.tocsr()
        logging.info('Solving...')
        self.U = spsolve(self.K, self.S)
        for e in self.elements:
            e.setUe(self.U)
        logging.info('Solved!')

    def postProcess(self, mult: float = 1000, gs=None, levels=1000, **kargs) -> None:
        """Generate the stress surfaces and displacement fields for the geometry

        Args:
                mult (int, optional): Factor for displacements. Defaults to 1000.
                gs (list, optional): List of 4 gridSpec matplotlib objects. Defaults to None.
        """

    def profile(self, region: list[float], n: float = 10) -> None:
        """Creates a profile in a given region coordinates

        Args:
            region (list[float]): List of region coordinates (square region 2D)
            n (float, optional): Number of points. Defaults to 10.

        Returns:
            np.ndarray: Coordinates, displacements and second variable solution
        """
        coords = np.array(region)
        gdl = np.array([[-1]*len(coords)])
        e = Quadrilateral(coords, gdl, n)
        _x, _ = e.T(e.domain.T)
        valuesU = []
        valuesDU = []
        for x in tqdm(_x, unit='point'):
            for e in self.elements:
                if e.isInside([x])[0]:
                    np.array([[1.3, 2.5, 3.5], [1.5, 2.6, 8.5]])
                    z = e.inverseMapping(x.reshape([3, 1]))
                    _, u, du = e.giveSolutionPoint(z, True)
                    valuesU += [u]
                    valuesDU += [du]

        return _x, np.array(valuesU), np.array(valuesDU)


class NonLocalElasticity(Elasticity):
    """Creates a 3D Elasticity problem

    Args:
        geometry (Geometry): Input 3D geometry
        E (Tuple[float, list]): Young Moduli
        v (Tuple[float, list]): Poisson coeficient
        rho (Tuple[float, list]): Density.
        l (float): Internal lenght
        z1 (float): Z1 factor
        Lr (float): Influence distance
        af (Callable): Atenuation function
        fx (Callable, optional): Force in x direction. Defaults to lambdax:0.
        fy (Callable, optional): Force in y direction. Defaults to lambdax:0.
        fz (Callable, optional): Force in z direction. Defaults to lambdax:0.
    """

    def __init__(self, geometry: Geometry, E: Tuple[float, list], v: Tuple[float, list], rho: Tuple[float, list], l: float, z1: float, Lr: float, af: Callable, fx: Callable = lambda x: 0, fy: Callable = lambda x: 0, fz: Callable = lambda x: 0, **kargs) -> None:
        """Creates a 3D Elasticity problem

        Args:
            geometry (Geometry): Input 3D geometry
            E (Tuple[float, list]): Young Moduli
            v (Tuple[float, list]): Poisson coeficient
            rho (Tuple[float, list]): Density.
            l (float): Internal lenght
            z1 (float): Z1 factor
            Lr (float): Influence distance
            af (Callable): Atenuation function
            fx (Callable, optional): Force in x direction. Defaults to lambdax:0.
            fy (Callable, optional): Force in y direction. Defaults to lambdax:0.
            fz (Callable, optional): Force in z direction. Defaults to lambdax:0.
        """
        Elasticity.__init__(self, geometry, E, v, rho, fx, fy, fz, **kargs)
        self.l = l
        self.z1 = z1
        self.z2 = 1.0-z1

        self.af = af
        self.Lr = Lr
        nonlocals = self.geometry.detectNonLocal(Lr)
        for e, dno in zip(self.elements, nonlocals):
            e.enl = dno

    def elementMatrices(self) -> None:
        """Calculate the element matrices usign Reddy's (2005) finite element model
        """

        for ee, e in enumerate(tqdm(self.elements, unit='Element')):
            m = len(e.gdl.T)

            # Gauss points in global coordinates and Shape functions evaluated in gauss points
            _x, _p = e.T(e.Z.T)
            # Jacobian evaluated in gauss points and shape functions derivatives in natural coordinates
            jac, dpz = e.J(e.Z.T)
            detjac = np.linalg.det(jac)
            _j = np.linalg.inv(jac)  # Jacobian inverse
            dpx = _j @ dpz  # Shape function derivatives in global coordinates

            E = self.E[ee]
            v = self.v[ee]

            C = E/((1.0+v)*(1.0-2.0*v))*np.array([
                [1.0-v, v, v, 0.0, 0.0, 0.0],
                [v, 1.0-v, v, 0.0, 0.0, 0.0],
                [v, v, 1.0-v, 0.0, 0.0, 0.0],
                [0.0, 0.0, 0.0, (1.0-2.0*v)/2.0, 0.0, 0.0],
                [0.0, 0.0, 0.0, 0.0, (1.0-2.0*v)/2.0, 0.0],
                [0.0, 0.0, 0.0, 0.0, 0.0, (1.0-2.0*v)/2.0]])

            o = [0.0]*m
            Ke = np.zeros([3*m, 3*m])
            Me = np.zeros([3*m, 3*m])
            Fe = np.zeros([3*m, 1])

            for k in range(len(e.Z)):  # Iterate over gauss points on domain
                B = np.array([
                    [*dpx[k, 0, :], *o, *o],
                    [*o, *dpx[k, 1, :], *o],
                    [*o, *o, *dpx[k, 2, :]],
                    [*dpx[k, 2, :], *o, *dpx[k, 0, :]],
                    [*o, *dpx[k, 2, :], *dpx[k, 1, :]],
                    [*dpx[k, 1, :], *dpx[k, 0, :], *o]])
                P = np.array([
                    [*_p[k], *o, *o],
                    [*o, *_p[k], *o],
                    [*o, *o, *_p[k]]])

                Ke += (B.T@C@B)*detjac[k]*e.W[k]
                Me += self.rho[ee]*(P.T@P)*detjac[k]*e.W[k]

                _fx = self.fx(_x[k])
                _fy = self.fy(_x[k])
                _fz = self.fz(_x[k])

                F = np.array([[_fx], [_fy], [_fz]])
                Fe += (P.T @ F)*detjac[k]*e.W[k]

            self.F[np.ix_(e.gdlm)] += Fe

            for gdl in e.gdlm:
                self.I += [gdl]*(3*m)
                self.J += e.gdlm
                self.Im += [gdl]*(3*m)
                self.Jm += e.gdlm

            self.V += (Ke*self.z1).flatten().tolist()
            self.Vm += Me.flatten().tolist()

            e.knls = []
            for inl in tqdm(e.enl, unit=' Nolocal'):
                enl = self.elements[inl]
                mnl = len(enl.gdl.T)
                Knl = np.zeros([3*m, 3*mnl])
                _xnl, _ = enl.T(enl.Z.T)
                jacnl, dpznl = enl.J(enl.Z.T)
                detjacnl = np.linalg.det(jacnl)
                _jnl = np.linalg.inv(jacnl)
                dpxnl = _jnl @ dpznl

                for k in range(len(e.Z)):
                    B = np.array([
                        [*dpx[k, 0, :], *o, *o],
                        [*o, *dpx[k, 1, :], *o],
                        [*o, *o, *dpx[k, 2, :]],
                        [*dpx[k, 2, :], *o, *dpx[k, 0, :]],
                        [*o, *dpx[k, 2, :], *dpx[k, 1, :]],
                        [*dpx[k, 1, :], *dpx[k, 0, :], *o]])

                    for knl in range(len(enl.Z)):
                        ro = np.linalg.norm(_x[k]-_xnl[knl])/self.l
                        azn = self.af(ro)
                        Bnl = np.array([
                            [*dpxnl[knl, 0, :], *o, *o],
                            [*o, *dpxnl[knl, 1, :], *o],
                            [*o, *o, *dpxnl[knl, 2, :]],
                            [*dpxnl[knl, 2, :], *o, *dpxnl[knl, 0, :]],
                            [*o, *dpxnl[knl, 2, :], *dpxnl[knl, 1, :]],
                            [*dpxnl[knl, 1, :], *dpxnl[knl, 0, :], *o]])

                        Knl += azn*(Bnl.T@C@B)*detjac[k] * \
                            e.W[k]*detjacnl[knl]*enl.W[knl]
                for gdl in e.gdlm:
                    self.I += [gdl]*(3*m)
                    self.J += enl.gdlm
                self.V += (Knl*self.z2).flatten().tolist()

    def ensembling(self) -> None:
        """Creation of the system sparse matrix. Force vector is ensembled in integration method
        """
        logging.info('Ensembling equation system...')
        self.K = sparse.coo_matrix(
            (self.V, (self.I, self.J)), shape=(self.ngdl, self.ngdl)).tolil()
        self.M = sparse.coo_matrix(
            (self.Vm, (self.Im, self.Jm)), shape=(self.ngdl, self.ngdl)).tocsr()
        logging.info('Done!')
