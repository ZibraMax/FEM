"""3D Isotropic elasticity"""

from typing import Callable, Tuple
import numpy as np
from tqdm import tqdm
from scipy import sparse
from FEM.Solvers.Lineal import LinealSparse
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
        if isinstance(E, float) or isinstance(E, int):
            E = [E]*len(geometry.elements)
        if isinstance(v, float) or isinstance(v, int):
            v = [v]*len(geometry.elements)
        if isinstance(rho, float) or isinstance(rho, int):
            rho = [rho]*len(geometry.elements)
        self.E = E
        self.v = v
        self.C = []
        for E, v in zip(self.E, self.v):
            c = E/((1.0+v)*(1.0-2.0*v))*np.array([
                [1.0-v, v, v, 0.0, 0.0, 0.0],
                [v, 1.0-v, v, 0.0, 0.0, 0.0],
                [v, v, 1.0-v, 0.0, 0.0, 0.0],
                [0.0, 0.0, 0.0, (1.0-2.0*v)/2.0, 0.0, 0.0],
                [0.0, 0.0, 0.0, 0.0, (1.0-2.0*v)/2.0, 0.0],
                [0.0, 0.0, 0.0, 0.0, 0.0, (1.0-2.0*v)/2.0]])
            self.C.append(c)
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

        self.K = sparse.lil_matrix((self.ngdl, self.ngdl))
        self.M = sparse.lil_matrix((self.ngdl, self.ngdl))
        self.name = 'Isotropic Elasticity sparse'
        self.properties['E'] = self.E
        self.properties['v'] = self.v
        self.properties['fx'] = None
        self.properties['fy'] = None
        self.properties['fz'] = None
        self.properties['rho'] = self.rho

    def elementMatrices(self) -> None:
        """Calculate the element matrices usign Reddy's (2005) finite element model
        """

        for ee, e in enumerate(tqdm(self.elements, unit='Element')):
            m = len(e.gdl.T)

            _x, _p = e._x, e._p
            # jac, dpz = e.jacs, e.dpz
            detjac = e.detjac
            # _j = np.linalg.inv(jac)
            dpx = e.dpx

            C = self.C[ee]
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
            self.K[np.ix_(e.gdlm, e.gdlm)] += Ke
            self.M[np.ix_(e.gdlm, e.gdlm)] += Me

    def ensembling(self) -> None:
        """Creation of the system sparse matrix. Force vector is ensembled in integration method
        """
        logging.info('Ensembling equation system...')
        self.M = self.M.tocsr()
        logging.info('Done!')

    def postProcess(self, **kargs) -> None:
        """Calculate stress and strain for each element gauss point and vertice. The results are stored in each element sigmas and epsilons properties.
        """
        for i, e in enumerate(self.elements):
            _, _, du = e.giveSolution(True)
            exx = du[:, 0, 0]
            eyy = du[:, 1, 1]
            ezz = du[:, 2, 2]
            exy = du[:, 0, 1]+du[:, 1, 0]
            eyz = du[:, 1, 2]+du[:, 2, 1]
            exz = du[:, 0, 2]+du[:, 2, 0]
            epsilons = np.array([exx, eyy, ezz, exz, eyz, exy])
            C = self.C[i]
            e.sigmas = (C @ epsilons).T
            e.epsilons = epsilons.T

    # TODO Region 2D is a more robust class for this job
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
        e = Quadrilateral(coords, gdl, n, fast=True)
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

        self.properties['l'] = self.l
        self.properties['z1'] = self.z1
        self.properties['z2'] = self.z2
        self.properties['Lr'] = self.Lr
        self.properties['af'] = None

        nonlocals = self.geometry.detectNonLocal(Lr)
        for e, dno in zip(self.elements, nonlocals):
            e.enl = dno
        self.name = 'Non Local Elasticity sparse-lil'
        self.KL = sparse.lil_matrix((self.ngdl, self.ngdl))
        self.KNL = sparse.lil_matrix((self.ngdl, self.ngdl))
        self.M = sparse.lil_matrix((self.ngdl, self.ngdl))

    def elementMatrices(self) -> None:
        """Calculate the element matrices usign Reddy's (2005) finite element model
        """

        for ee, e in enumerate(tqdm(self.elements, unit='Element')):
            m = len(e.gdl.T)

            # Gauss points in global coordinates and Shape functions evaluated in gauss points
            _x, _p = e._x, e._p
            # Jacobian evaluated in gauss points and shape functions derivatives in natural coordinates
            # jac, dpz = e.J(e.Z.T)
            detjac = e.detjac
            # _j = np.linalg.inv(jac)  # Jacobian inverse
            dpx = e.dpx  # Shape function derivatives in global coordinates

            C = self.C[ee]
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
            self.KL[np.ix_(e.gdlm, e.gdlm)] += Ke
            self.M[np.ix_(e.gdlm, e.gdlm)] += Me
            # print('Ensembling equation system...')
            # for e in tqdm(self.elements, unit='Element'):
            #     self.K[np.ix_(e.gdlm, e.gdlm)] += e.Ke*self.z1
            #     for i, eee in enumerate(e.enl):
            #         enl = self.elements[eee]
            #         self.K[np.ix_(e.gdlm, enl.gdlm)] += e.knls[i]*self.z2
            #     self.F[np.ix_(e.gdlm)] += e.Fe
            #     self.Q[np.ix_(e.gdlm)] += e.Qe
            # print('Done!')

            e.knls = []
            for inl in tqdm(e.enl, unit=' Nolocal'):
                enl = self.elements[inl]
                mnl = len(enl.gdl.T)
                Knl = np.zeros([3*m, 3*mnl])
                _xnl = enl._x
                detjacnl = enl.detjac
                dpxnl = enl.dpx

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
                self.KNL[np.ix_(e.gdlm, enl.gdlm)] += Knl.T

    def ensembling(self) -> None:
        """Creation of the system sparse matrix. Force vector is ensembled in integration method
        """
        logging.info('Ensembling equation system...')
        self.K = self.KL*self.z1 + self.KNL*self.z2
        self.M = self.M.tocsr()
        logging.info('Done!')


class NonLocalElasticityFromTensor(NonLocalElasticity):

    def __init__(self, geometry: Geometry, C: np.ndarray, rho: Tuple[float, list], l: float, z1: float, Lr: float, af: Callable, fx: Callable = lambda x: 0, fy: Callable = lambda x: 0, fz: Callable = lambda x: 0, **kargs) -> None:
        NonLocalElasticity.__init__(
            self, geometry, 0.0, 0.0, rho, l, z1, Lr, af, fx, fy, fz, **kargs)
        self.properties['C'] = C.tolist()
        self.C = [C]*len(geometry.elements)


class ElasticityFromTensor(Elasticity):

    def __init__(self, geometry: Geometry, C: np.ndarray, rho: Tuple[float, list], fx: Callable = lambda x: 0, fy: Callable = lambda x: 0, fz: Callable = lambda x: 0, **kargs) -> None:
        Elasticity.__init__(self, geometry, 0.0, 0.0,
                            rho, fx, fy, fz, **kargs)
        self.properties['C'] = C.tolist()
        self.C = [C]*len(geometry.elements)


class NonLocalElasticityLegacy(Elasticity):
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
        self.name = 'Non Local Elasticity sparse'
        self.zs = []

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

            C = self.C[ee]
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
            Ke_flat = (Ke).flatten().tolist()
            self.zs += [1]*len(Ke_flat)
            self.V += Ke_flat
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
                Knl_flat = (Knl).flatten().tolist()
                self.V += (Knl).flatten().tolist()
                self.zs += [0]*len(Knl_flat)
        self.V_0 = np.array(self.V)
        self.zs = np.array(self.zs)

    def ensembling(self) -> None:
        """Creation of the system sparse matrix. Force vector is ensembled in integration method
        """
        logging.info('Ensembling equation system...')
        self.V = np.zeros(self.V_0.shape)
        self.V[self.zs == 1] = self.V_0[self.zs == 1]*self.z1
        self.V[self.zs == 0] = self.V_0[self.zs == 0]*self.z2

        self.K = sparse.coo_matrix(
            (self.V, (self.I, self.J)), shape=(self.ngdl, self.ngdl)).tolil()
        self.M = sparse.coo_matrix(
            (self.Vm, (self.Im, self.Jm)), shape=(self.ngdl, self.ngdl)).tocsr()
        logging.info('Done!')
