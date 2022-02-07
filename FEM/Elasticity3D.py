
from typing import Callable, Tuple

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import gridspec
from tqdm import tqdm
from scipy import sparse
from scipy.sparse.linalg import spsolve

from .Core import Core, Geometry, logging


class Elasticity(Core):

    def __init__(self, geometry: Geometry, E: Tuple[float, list], v: Tuple[float, list], rho: Tuple[float, list], fx: Callable = lambda x: 0, fy: Callable = lambda x: 0, fz: Callable = lambda x: 0, **kargs) -> None:
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
        self.J = []
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
        # X = []
        # Y = []
        # U1 = []
        # U2 = []
        # U3 = []
        # fig = plt.figure()
        # if not gs:
        #     gss = gridspec.GridSpec(3, 3)
        #     gs = [gss[0, 0], gss[0, 1], gss[0, 2], gss[1:, :]]
        # ax1 = fig.add_subplot(gs[0])
        # ax2 = fig.add_subplot(gs[1])
        # ax3 = fig.add_subplot(gs[2])
        # ax5 = fig.add_subplot(gs[3])
        # ee = -1
        # for e in tqdm(self.elements, unit='Element'):
        #     ee += 1
        #     _x, _u, du = e.giveSolution(True)
        #     X += _x.T[0].tolist()
        #     Y += _x.T[1].tolist()
        #     U1 += (self.C11[ee]*du[:, 0, 0]+self.C12[ee]*du[:, 1, 1]).tolist()
        #     U2 += (self.C12[ee]*du[:, 0, 0]+self.C11[ee]*du[:, 1, 1]).tolist()
        #     U3 += (self.C66[ee]*(du[:, 0, 1]+du[:, 1, 0])).tolist()
        #     coordsNuevas = e._coordsg + e._Ueg * mult
        #     ax5.plot(*e._coordsg.T, '--', color='gray', alpha=0.7)
        #     ax5.plot(*coordsNuevas.T, '-', color='black')
        # ax5.legend(['Original Shape', 'Deformed Shape (x'+format(mult)+')'])
        # ax5.set_aspect('equal')
        # ax1.set_aspect('equal')
        # ax3.set_aspect('equal')
        # ax2.set_aspect('equal')
        # cmap = 'rainbow'
        # def fmt(x): return format(x, '.3f')

        # surf = ax1.tricontourf(X, Y, U1, cmap=cmap, levels=levels, **kargs)
        # plt.colorbar(surf, ax=ax1)
        # ax1.set_title(r'$\sigma_{xx}$')

        # surf = ax2.tricontourf(X, Y, U2, cmap=cmap, levels=levels, **kargs)
        # plt.colorbar(surf, ax=ax2)
        # ax2.set_title(r'$\sigma_{yy}$')

        # surf = ax3.tricontourf(X, Y, U3, cmap=cmap, levels=levels, **kargs)
        # plt.colorbar(surf, ax=ax3)
        # ax3.set_title(r'$\sigma_{xy}$')

        # mask = self.geometry.mask
        # if self.geometry.holes:
        #     for hole in self.geometry.holes:
        #         Xs = np.array(hole['vertices'])[:, 0]
        #         Ys = np.array(hole['vertices'])[:, 1]
        #         ax1.fill(Xs, Ys, color='white', zorder=30)
        #         ax2.fill(Xs, Ys, color='white', zorder=30)
        #         ax3.fill(Xs, Ys, color='white', zorder=30)
        # if not mask == None:
        #     mask = np.array(mask)
        #     cornersnt = np.array(mask[::-1])

        #     xmin = np.min(cornersnt[:, 0])
        #     xmax = np.max(cornersnt[:, 0])

        #     ymin = np.min(cornersnt[:, 1])
        #     ymax = np.max(cornersnt[:, 1])

        #     Xs = [xmin, xmax, xmax, xmin]+cornersnt[:, 0].tolist()
        #     Ys = [ymin, ymin, ymax, ymax]+cornersnt[:, 1].tolist()
        #     ax1.fill(Xs, Ys, color='white', zorder=30)
        #     ax2.fill(Xs, Ys, color='white', zorder=30)
        #     ax3.fill(Xs, Ys, color='white', zorder=30)

    def giveStressPoint(self, X: np.ndarray) -> Tuple[tuple, None]:
        """Calculates the stress in a given set of points.

        Args:
            X (np.ndarray): Points to calculate the Stress. 2D Matrix. with 2 rows. First row is an array of 1 column with X coordinate. Second row is an array of 1 column with Y coordinate

        Returns:
            tuple or None: Tuple of stress (:math:`\sigma_x,\sigma_y,\sigma_{xy}`) if X,Y exists in domain.
        """
        # for ee, e in enumerate(self.elements):
        #     if e.isInside(X.T[0]):
        #         z = e.inverseMapping(np.array([X.T[0]]).T)
        #         _, _, du = e.giveSolutionPoint(z, True)
        #         # TODO Arreglar calculo de esfuerzos para PlaneStrain
        #         sx = (self.C11[ee]*du[:, 0, 0] +
        #               self.C12[ee]*du[:, 1, 1]).tolist()
        #         sy = (self.C12[ee]*du[:, 0, 0] +
        #               self.C11[ee]*du[:, 1, 1]).tolist()
        #         sxy = (self.C66[ee]*(du[:, 0, 1]+du[:, 1, 0])).tolist()
        # return sx, sy, sxy
        pass

    def profile(self, region: int, vn: list[int] = None, n: float = 100) -> None:
        """Generate a profile between selected points

        Args:
            region (int): Geometry region 
            n (int, optional): NUmber of samples for graph. Defaults to 100.

        """
        pass
        # _x = np.linspace(p0[0], p1[0], n)
        # _y = np.linspace(p0[1], p1[1], n)
        # X = np.array([_x, _y])
        # U = []
        # U1 = []
        # U2 = []
        # U3 = []
        # _X = []
        # def dist(X): return np.sqrt((p0[0]-X[0])**2+(p0[1]-X[1])**2)
        # for i in range(n):
        #     for ee, e in enumerate(self.elements):
        #         if e.isInside(X.T[i]):
        #             z = e.inverseMapping(np.array([X.T[i]]).T)
        #             _, u, du = e.giveSolutionPoint(z, True)
        #             # TODO Arreglar calculo de esfuerzos para PlaneStrain
        #             U += [u.tolist()]
        #             U1 += (self.C11[ee]*du[:, 0, 0] +
        #                    self.C12[ee]*du[:, 1, 1]).tolist()
        #             U2 += (self.C12[ee]*du[:, 0, 0] +
        #                    self.C11[ee]*du[:, 1, 1]).tolist()
        #             U3 += (self.C66[ee]*(du[:, 0, 1]+du[:, 1, 0])).tolist()
        #             _X.append(dist(X.T[i]))
        #             break
        # fig = plt.figure()
        # ax = fig.add_subplot(1, 3, 1)
        # ax.plot(_X, U1, color='black')
        # ax.grid()
        # ax.set_xlabel('d')
        # ax.set_ylabel(r'$\sigma_{xx}$')
        # ax = fig.add_subplot(1, 3, 2)
        # ax.plot(_X, U2, color='black')
        # ax.grid()
        # ax.set_xlabel('d')
        # ax.set_ylabel(r'$\sigma_{yy}$')
        # ax = fig.add_subplot(1, 3, 3)
        # ax.plot(_X, U3, color='black')
        # ax.grid()
        # ax.set_xlabel('d')
        # ax.set_ylabel(r'$\sigma_{xy}$')
        # return _X, U


class NonLocalElasticity(Elasticity):

    def __init__(self, geometry: Geometry, E: Tuple[float, list], v: Tuple[float, list], rho: Tuple[float, list], l: float, z1: float, Lr: float, af: Callable, fx: Callable = lambda x: 0, fy: Callable = lambda x: 0, fz: Callable = lambda x: 0, **kargs) -> None:
        Elasticity.__init__(self, geometry, E, v, rho, fx, fy, fz, **kargs)
        self.l = l
        self.z1 = z1
        self.af = af
        self.Lr = Lr

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
            self.V += (Ke*self.z1).flatten().tolist()
            self.Vm += Me.flatten().tolist()

            e.knls = []
            for inl in tqdm(e.enl, unit=' Nolocal'):
                enl = self.elements[inl]
                mnl = len(enl.gdl.T)
                Knl = np.zeros([2*m, 2*mnl])
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
                        azn = self.af(self.l0, ro)
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
                    self.I += [gdl]*(m+mnl)
                    self.J += enl.gdlm
                self.V += (Knl*self.z2).flatten().tolist()

                # TODO detección de elementos No Locales
