"""2D Elasticity
"""


from typing import Callable, Tuple

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import gridspec
from tqdm import tqdm
from scipy import sparse

from .Solvers import LinealSparse


from .Core import Core, Geometry, logging


class PlaneStressOrthotropic(Core):

    """Creates a plane stress problem with orthotropic formulation

        Args:
            geometry (Geometry): Input geometry
            E1 (Tuple[float, list]): Young moduli in direction 1 (x)
            E2 (Tuple[float, list]): Young moduli in direction 2 (y)
            G12 (Tuple[float, list]): Shear moduli
            v12 (Tuple[float, list]): Poisson moduli
            t (Tuple[float, list]): Thickness
            rho (Tuple[float, list], optional): Density. If not given, mass matrix will not be calculated. Defaults to None.
            fx (Callable, optional): Force in x direction. Defaults to lambdax:0.
            fy (Callable, optional): Force in y direction. Defaults to lambdax:0.
        """

    def __init__(self, geometry: Geometry, E1: Tuple[float, list], E2: Tuple[float, list], G12: Tuple[float, list], v12: Tuple[float, list], t: Tuple[float, list], rho: Tuple[float, list] = None, fx: Callable = lambda x: 0, fy: Callable = lambda x: 0, **kargs) -> None:
        """Creates a plane stress problem with orthotropic formulation

        Args:
            geometry (Geometry): Input geometry
            E1 (Tuple[float, list]): Young moduli in direction 1 (x)
            E2 (Tuple[float, list]): Young moduli in direction 2 (y)
            G12 (Tuple[float, list]): Shear moduli
            v12 (Tuple[float, list]): Poisson moduli
            t (Tuple[float, list]): Thickness
            rho (Tuple[float, list], optional): Density. If not given, mass matrix will not be calculated. Defaults to None.
            fx (Callable, optional): Force in x direction. Defaults to lambdax:0.
            fy (Callable, optional): Force in y direction. Defaults to lambdax:0.
        """
        if isinstance(t, float) or isinstance(t, int):
            t = [t]*len(geometry.elements)
        if isinstance(E1, float) or isinstance(E1, int):
            E1 = [E1]*len(geometry.elements)
        if isinstance(E2, float) or isinstance(E2, int):
            E2 = [E2]*len(geometry.elements)
        if isinstance(G12, float) or isinstance(G12, int):
            G12 = [G12]*len(geometry.elements)
        if isinstance(v12, float) or isinstance(v12, int):
            v12 = [v12]*len(geometry.elements)
        self.calculateMass = False
        if rho:
            if isinstance(rho, int) or isinstance(rho, float):
                self.rho = [rho]*len(geometry.elements)
                self.calculateMass = True
        self.t = t
        self.E1 = E1
        self.E2 = E2
        self.G12 = G12
        self.v12 = v12
        self.v21 = []
        self.C11 = []
        self.C22 = []
        self.C12 = []
        self.C66 = []
        self.fx = fx
        self.fy = fy
        for i in range(len(self.E1)):
            v21 = self.v12[i]*self.E2[i]/self.E1[i]
            self.v21.append(v21)
            C11 = self.E1[i] / (1 - self.v12[i]*self.v21[i])
            C22 = self.E2[i] / (1 - self.v12[i]*self.v21[i])
            C12 = self.v21[i] * C11
            C66 = G12[i]
            self.C11.append(C11)
            self.C22.append(C22)
            self.C12.append(C12)
            self.C66.append(C66)
        if not geometry.nvn == 2:
            logging.warning(
                'Border conditions lost, please usea a geometry with 2 variables per node (nvn=2)\nRegenerating Geoemtry...')
            geometry.nvn = 2
            geometry.cbe = []
            geometry.cbn = []
            geometry.initialize()
        Core.__init__(self, geometry, **kargs)
        self.name = 'Plane Stress Orthotropic'

    def elementMatrices(self) -> None:
        """Calculate the element matrices usign Reddy's (2005) finite element model
        """

        for ee, e in enumerate(tqdm(self.elements, unit='Element')):
            m = len(e.gdl.T)
            Fux = np.zeros([m, 1])
            Fvx = np.zeros([m, 1])
            # Gauss points in global coordinates and Shape functions evaluated in gauss points
            _x, _p = e.T(e.Z.T)
            # Jacobian evaluated in gauss points and shape functions derivatives in natural coordinates
            jac, dpz = e.J(e.Z.T)
            detjac = np.linalg.det(jac)
            _j = np.linalg.inv(jac)  # Jacobian inverse
            dpx = _j @ dpz  # Shape function derivatives in global coordinates
            c11 = self.C11[ee]
            c12 = self.C12[ee]
            c22 = self.C22[ee]
            c66 = self.C66[ee]
            C = np.array([
                [c11, c12, 0.0],
                [c12, c22, 0.0],
                [0.0, 0.0, c66]])

            Fe = np.zeros([2*m, 1])
            Ke = np.zeros([2*m, 2*m])
            if self.calculateMass:
                Me = np.zeros([2*m, 2*m])

            o = [0.0]*m
            for k in range(len(e.Z)):  # Iterate over gauss points on domain
                B = np.array([
                    [*dpx[k, 0, :], *o],
                    [*o, *dpx[k, 1, :]],
                    [*dpx[k, 1, :], *dpx[k, 0, :]]])

                P = np.array([
                    [*_p[k], *o],
                    [*o, *_p[k]]])
                Ke += self.t[ee]*(B.T@C@B)*detjac[k]*e.W[k]
                if self.calculateMass:
                    Me += self.rho[ee]*self.t[ee]*(P.T@P)*detjac[k]*e.W[k]
                Fe += self.t[ee]*(P.T@np.array([[self.fx(_x[k])],
                                                [self.fy(_x[k])]]))*detjac[k]*e.W[k]

            if e.intBorders:  # Cambiar esto a la notación matricial
                for j in range(len(e.borders)):
                    border = e.borders[j]
                    if (len(border.properties['load_x']) + len(border.properties['load_y'])):
                        _x, _p = e.T(e.Tj[j](border.Z.T))
                        _s = border.TS(border.Z.T)
                        detjac = border.coords[-1, 0]*0.5
                        for i in range(m):
                            for fx in border.properties['load_x']:
                                for k in range(len(border.Z)):
                                    Fux[i, 0] += fx(_s[k, 0])*_p[k, i] * \
                                        detjac*border.W[k]
                            for fy in border.properties['load_y']:
                                for k in range(len(border.Z)):
                                    Fvx[i, 0] += fy(_s[k, 0])*_p[k, i] * \
                                        detjac*border.W[k]
            subm = np.linspace(0, 2*m-1, 2*m).reshape([2, m]).astype(int)
            e.Fe = Fe
            e.Ke = Ke
            if self.calculateMass:
                e.Me = Me
            e.Fe[np.ix_(subm[0])] += Fux
            e.Fe[np.ix_(subm[1])] += Fvx

    def postProcess(self, mult: float = 1000, gs=None, levels=1000, **kargs) -> None:
        """Generate the stress surfaces and displacement fields for the geometry

        Args:
                mult (int, optional): Factor for displacements. Defaults to 1000.
                gs (list, optional): List of 4 gridSpec matplotlib objects. Defaults to None.
        """
        X = []
        Y = []
        U1 = []
        U2 = []
        U3 = []
        fig = plt.figure()
        if not gs:
            gss = gridspec.GridSpec(3, 3)
            gs = [gss[0, 0], gss[0, 1], gss[0, 2], gss[1:, :]]
        ax1 = fig.add_subplot(gs[0])
        ax2 = fig.add_subplot(gs[1])
        ax3 = fig.add_subplot(gs[2])
        ax5 = fig.add_subplot(gs[3])
        ee = -1
        for e in tqdm(self.elements, unit='Element'):
            ee += 1
            _x, _u, du = e.giveSolution(True)
            X += _x.T[0].tolist()
            Y += _x.T[1].tolist()
            U1 += (self.C11[ee]*du[:, 0, 0]+self.C12[ee]*du[:, 1, 1]).tolist()
            U2 += (self.C12[ee]*du[:, 0, 0]+self.C22[ee]*du[:, 1, 1]).tolist()
            U3 += (self.C66[ee]*(du[:, 0, 1]+du[:, 1, 0])).tolist()
            coordsNuevas = e._coordsg + e._Ueg * mult
            ax5.plot(*e._coordsg.T, '--', color='gray', alpha=0.7)
            ax5.plot(*coordsNuevas.T, '-', color='black')
        ax5.legend(['Original Shape', 'Deformed Shape (x'+format(mult)+')'])
        ax5.set_aspect('equal')
        ax1.set_aspect('equal')
        ax3.set_aspect('equal')
        ax2.set_aspect('equal')
        cmap = 'rainbow'

        surf = ax1.tricontourf(X, Y, U1, cmap=cmap, levels=levels, **kargs)
        plt.colorbar(surf, ax=ax1)
        ax1.set_title(r'$\sigma_{xx}$')

        surf = ax2.tricontourf(X, Y, U2, cmap=cmap, levels=levels, **kargs)
        plt.colorbar(surf, ax=ax2)
        ax2.set_title(r'$\sigma_{yy}$')

        surf = ax3.tricontourf(X, Y, U3, cmap=cmap, levels=levels, **kargs)
        plt.colorbar(surf, ax=ax3)
        ax3.set_title(r'$\sigma_{xy}$')

        mask = self.geometry.mask
        if self.geometry.holes:
            for hole in self.geometry.holes:
                Xs = np.array(hole['vertices'])[:, 0]
                Ys = np.array(hole['vertices'])[:, 1]
                ax1.fill(Xs, Ys, color='white', zorder=30)
                ax2.fill(Xs, Ys, color='white', zorder=30)
                ax3.fill(Xs, Ys, color='white', zorder=30)
        if mask:
            mask = np.array(mask)
            cornersnt = np.array(mask[::-1])

            xmin = np.min(cornersnt[:, 0])
            xmax = np.max(cornersnt[:, 0])

            ymin = np.min(cornersnt[:, 1])
            ymax = np.max(cornersnt[:, 1])

            Xs = [xmin, xmax, xmax, xmin]+cornersnt[:, 0].tolist()
            Ys = [ymin, ymin, ymax, ymax]+cornersnt[:, 1].tolist()
            ax1.fill(Xs, Ys, color='white', zorder=30)
            ax2.fill(Xs, Ys, color='white', zorder=30)
            ax3.fill(Xs, Ys, color='white', zorder=30)

    def giveStressPoint(self, X: np.ndarray) -> Tuple[tuple, None]:
        """Calculates the stress in a given set of points.

        Args:
            X (np.ndarray): Points to calculate the Stress. 2D Matrix. with 2 rows. First row is an array of 1 column with X coordinate. Second row is an array of 1 column with Y coordinate

        Returns:
            tuple or None: Tuple of stress (:math:`\sigma_x,\sigma_y,\sigma_{xy}`) if X,Y exists in domain.
        """
        for ee, e in enumerate(self.elements):
            if e.isInside(X.T[0]):
                z = e.inverseMapping(np.array([X.T[0]]).T)
                _, _, du = e.giveSolutionPoint(z, True)
                sx = (self.C11[ee]*du[:, 0, 0] +
                      self.C12[ee]*du[:, 1, 1]).tolist()
                sy = (self.C12[ee]*du[:, 0, 0] +
                      self.C22[ee]*du[:, 1, 1]).tolist()
                sxy = (self.C66[ee]*(du[:, 0, 1]+du[:, 1, 0])).tolist()
                return sx, sy, sxy

    def profile(self, p0: list, p1: list, n: float = 100) -> None:
        """
        Generate a profile between selected points

        Args:
            p0 (list): start point of the profile [x0,y0]
            p1 (list): end point of the profile [xf,yf]
            n (int, optional): NUmber of samples for graph. Defaults to 100.

        """
        _x = np.linspace(p0[0], p1[0], n)
        _y = np.linspace(p0[1], p1[1], n)
        X = np.array([_x, _y])
        U = []
        U1 = []
        U2 = []
        U3 = []
        _X = []
        def dist(X): return np.sqrt((p0[0]-X[0])**2+(p0[1]-X[1])**2)
        for i in range(n):
            for ee, e in enumerate(self.elements):
                if e.isInside(X.T[i]):
                    z = e.inverseMapping(np.array([X.T[i]]).T)
                    _, u, du = e.giveSolutionPoint(z, True)
                    U += [u.tolist()]
                    U1 += (self.C11[ee]*du[:, 0, 0] +
                           self.C12[ee]*du[:, 1, 1]).tolist()
                    U2 += (self.C12[ee]*du[:, 0, 0] +
                           self.C22[ee]*du[:, 1, 1]).tolist()
                    U3 += (self.C66[ee]*(du[:, 0, 1]+du[:, 1, 0])).tolist()
                    _X.append(dist(X.T[i]))
                    break
        fig = plt.figure()
        ax = fig.add_subplot(1, 3, 1)
        ax.plot(_X, U1, color='black')
        ax.grid()
        ax.set_xlabel('d')
        ax.set_ylabel(r'$\sigma_{xx}$')
        ax = fig.add_subplot(1, 3, 2)
        ax.plot(_X, U2, color='black')
        ax.grid()
        ax.set_xlabel('d')
        ax.set_ylabel(r'$\sigma_{yy}$')
        ax = fig.add_subplot(1, 3, 3)
        ax.plot(_X, U3, color='black')
        ax.grid()
        ax.set_xlabel('d')
        ax.set_ylabel(r'$\sigma_{xy}$')
        return _X, U


class PlaneStressOrthotropicSparse(PlaneStressOrthotropic):
    """Creates a plane stress problem with orthotropic formulation and sparce matrices

    Args:
        geometry (Geometry): Input geometry
        E1 (Tuple[float, list]): Young moduli in direction 1 (x)
        E2 (Tuple[float, list]): Young moduli in direction 2 (y)
        G12 (Tuple[float, list]): Shear moduli
        v12 (Tuple[float, list]): Poisson moduli
        t (Tuple[float, list]): Thickness
        rho (Tuple[float, list], optional): Density. If not given, mass matrix will not be calculated. Defaults to None.
        fx (Callable, optional): Force in x direction. Defaults to lambdax:0.
        fy (Callable, optional): Force in y direction. Defaults to lambdax:0.
    """

    def __init__(self, geometry: Geometry, E1: Tuple[float, list], E2: Tuple[float, list], G12: Tuple[float, list], v12: Tuple[float, list], t: Tuple[float, list], rho: Tuple[float, list] = None, fx: Callable = lambda x: 0, fy: Callable = lambda x: 0, **kargs) -> None:
        """Creates a plane stress problem with orthotropic formulation and sparce matrices

        Args:
            geometry (Geometry): Input geometry
            E1 (Tuple[float, list]): Young moduli in direction 1 (x)
            E2 (Tuple[float, list]): Young moduli in direction 2 (y)
            G12 (Tuple[float, list]): Shear moduli
            v12 (Tuple[float, list]): Poisson moduli
            t (Tuple[float, list]): Thickness
            rho (Tuple[float, list], optional): Density. If not given, mass matrix will not be calculated. Defaults to None.
            fx (Callable, optional): Force in x direction. Defaults to lambdax:0.
            fy (Callable, optional): Force in y direction. Defaults to lambdax:0.
        """
        if not geometry.nvn == 2:
            logging.warning(
                'Border conditions lost, please usea a geometry with 2 variables per node (nvn=2)\nRegenerating Geoemtry...')
            geometry.nvn = 2
            geometry.cbe = []
            geometry.cbn = []
            geometry.initialize()
        if not geometry.fast:
            logging.warning("Use fast elements")
            geometry.fast = True
            geometry.initialize()
        PlaneStressOrthotropic.__init__(
            self, geometry, E1, E2, G12, v12, t, rho, fx, fy, sparse=True, **kargs)
        self.K = sparse.lil_matrix((self.ngdl, self.ngdl))
        if self.calculateMass:
            self.M = sparse.lil_matrix((self.ngdl, self.ngdl))
        self.name = 'Plane Stress Orthotropic sparse'

    def elementMatrices(self) -> None:
        """Calculate the element matrices usign Reddy's (2005) finite element model
        """

        for ee, e in enumerate(tqdm(self.elements, unit='Element')):
            m = len(e.gdl.T)
            Fux = np.zeros([m, 1])
            Fvx = np.zeros([m, 1])
            # Gauss points in global coordinates and Shape functions evaluated in gauss points
            _x, _p = e.T(e.Z.T)
            # Jacobian evaluated in gauss points and shape functions derivatives in natural coordinates
            jac, dpz = e.J(e.Z.T)
            detjac = np.linalg.det(jac)
            _j = np.linalg.inv(jac)  # Jacobian inverse
            dpx = _j @ dpz  # Shape function derivatives in global coordinates
            c11 = self.C11[ee]
            c12 = self.C12[ee]
            c22 = self.C22[ee]
            c66 = self.C66[ee]
            C = np.array([
                [c11, c12, 0.0],
                [c12, c22, 0.0],
                [0.0, 0.0, c66]])

            Fe = np.zeros([2*m, 1])
            Ke = np.zeros([2*m, 2*m])
            if self.calculateMass:
                Me = np.zeros([2*m, 2*m])

            o = [0.0]*m
            for k in range(len(e.Z)):  # Iterate over gauss points on domain
                B = np.array([
                    [*dpx[k, 0, :], *o],
                    [*o, *dpx[k, 1, :]],
                    [*dpx[k, 1, :], *dpx[k, 0, :]]])

                P = np.array([
                    [*_p[k], *o],
                    [*o, *_p[k]]])
                Ke += self.t[ee]*(B.T@C@B)*detjac[k]*e.W[k]
                if self.calculateMass:
                    Me += self.rho[ee]*self.t[ee]*(P.T@P)*detjac[k]*e.W[k]
                Fe += self.t[ee]*(P.T@np.array([[self.fx(_x[k])],
                                                [self.fy(_x[k])]]))*detjac[k]*e.W[k]

            if e.intBorders:  # Cambiar esto a la notación matricial
                for j in range(len(e.borders)):
                    border = e.borders[j]
                    if (len(border.properties['load_x']) + len(border.properties['load_y'])):
                        _x, _p = e.T(e.Tj[j](border.Z.T))
                        _s = border.TS(border.Z.T)
                        detjac = border.coords[-1, 0]*0.5
                        for i in range(m):
                            for fx in border.properties['load_x']:
                                for k in range(len(border.Z)):
                                    Fux[i, 0] += fx(_s[k, 0])*_p[k, i] * \
                                        detjac*border.W[k]
                            for fy in border.properties['load_y']:
                                for k in range(len(border.Z)):
                                    Fvx[i, 0] += fy(_s[k, 0])*_p[k, i] * \
                                        detjac*border.W[k]
            subm = np.linspace(0, 2*m-1, 2*m).reshape([2, m]).astype(int)
            Fe[np.ix_(subm[0])] += Fux
            Fe[np.ix_(subm[1])] += Fvx
            self.F[np.ix_(e.gdlm)] += Fe

            self.K[np.ix_(e.gdlm, e.gdlm)] += Ke

            if self.calculateMass:
                self.M[np.ix_(e.gdlm, e.gdlm)] += Me

    def ensembling(self) -> None:
        """Creation of the system sparse matrix. Force vector is ensembled in integration method"""
        logging.info('Ensembling equation system...')
        if self.calculateMass:
            self.M = self.M.tocsr()
        logging.info('Done!')


class PlaneStress(PlaneStressOrthotropic):

    """Create a Plain Stress problem

    Args:
            geometry (Geometry): 2D 2 variables per node geometry
            E (int or float or list): Young Moduli. If number, all element will have the same young moduli. If list, each position will be the element young moduli, so len(E) == len(self.elements)
            v (int or float or list): Poisson ratio. If number, all element will have the same Poisson ratio. If list, each position will be the element Poisson ratio, so len(v) == len(self.elements)
            t (int or float or list): Element thickness. If number, all element will have the same thickness. If list, each position will be the element thickness, so len(t) == len(self.elements)
            rho (int or float or list, optional): Density. If not given, mass matrix will not be calculated. Defaults to None.
            fx (function, optional): Function fx, if fx is constant you can use fx = lambda x: [value]. Defaults to lambda x:0.
            fy (function, optional): Function fy, if fy is constant you can use fy = lambda x: [value]. Defaults to lambda x:0.
    """

    def __init__(self, geometry: Geometry, E: Tuple[float, list], v: Tuple[float, list], t: Tuple[float, list], rho: Tuple[float, list] = None, fx: Callable = lambda x: 0, fy: Callable = lambda x: 0, **kargs) -> None:
        """Create a Plain Stress problem

        Args:
                geometry (Geometry): 2D 2 variables per node geometry
                E (int or float or list): Young Moduli. If number, all element will have the same young moduli. If list, each position will be the element young moduli, so len(E) == len(self.elements)
                v (int or float or list): Poisson ratio. If number, all element will have the same Poisson ratio. If list, each position will be the element Poisson ratio, so len(v) == len(self.elements)
                t (int or float or list): Element thickness. If number, all element will have the same thickness. If list, each position will be the element thickness, so len(t) == len(self.elements)
                rho (int or float or list, optional): Density. If not given, mass matrix will not be calculated. Defaults to None.
                fx (function, optional): Function fx, if fx is constant you can use fx = lambda x: [value]. Defaults to lambda x:0.
                fy (function, optional): Function fy, if fy is constant you can use fy = lambda x: [value]. Defaults to lambda x:0.
        """
        G = E/2.0/(1.0+v)
        PlaneStressOrthotropic.__init__(
            self, geometry, E, E, G, v, t, rho, fx, fy, **kargs)
        self.name = 'Plane Stress Isotropic'


class PlaneStressSparse(PlaneStressOrthotropicSparse):

    """Create a Plain Stress problem using sparse matrices

    Args:
            geometry (Geometry): 2D 2 variables per node geometry
            E (int or float or list): Young Moduli. If number, all element will have the same young moduli. If list, each position will be the element young moduli, so len(E) == len(self.elements)
            v (int or float or list): Poisson ratio. If number, all element will have the same Poisson ratio. If list, each position will be the element Poisson ratio, so len(v) == len(self.elements)
            t (int or float or list): Element thickness. If number, all element will have the same thickness. If list, each position will be the element thickness, so len(t) == len(self.elements)
            rho (int or float or list, optional): Density. If not given, mass matrix will not be calculated. Defaults to None.
            fx (function, optional): Function fx, if fx is constant you can use fx = lambda x: [value]. Defaults to lambda x:0.
            fy (function, optional): Function fy, if fy is constant you can use fy = lambda x: [value]. Defaults to lambda x:0.
    """

    def __init__(self, geometry: Geometry, E: Tuple[float, list], v: Tuple[float, list], t: Tuple[float, list], rho: Tuple[float, list] = None, fx: Callable = lambda x: 0, fy: Callable = lambda x: 0, **kargs) -> None:
        """Create a Plain Stress problem using sparse matrices

        Args:
                geometry (Geometry): 2D 2 variables per node geometry
                E (int or float or list): Young Moduli. If number, all element will have the same young moduli. If list, each position will be the element young moduli, so len(E) == len(self.elements)
                v (int or float or list): Poisson ratio. If number, all element will have the same Poisson ratio. If list, each position will be the element Poisson ratio, so len(v) == len(self.elements)
                t (int or float or list): Element thickness. If number, all element will have the same thickness. If list, each position will be the element thickness, so len(t) == len(self.elements)
                rho (int or float or list, optional): Density. If not given, mass matrix will not be calculated. Defaults to None.
                fx (function, optional): Function fx, if fx is constant you can use fx = lambda x: [value]. Defaults to lambda x:0.
                fy (function, optional): Function fy, if fy is constant you can use fy = lambda x: [value]. Defaults to lambda x:0.
        """
        G = E/2.0/(1.0+v)
        PlaneStressOrthotropicSparse.__init__(
            self, geometry, E, E, G, v, t, rho, fx, fy, **kargs)
        self.name = 'Plane Stress Isotropic sparse'


class PlaneStrain(PlaneStress):
    """Create a Plain Strain problem

    Args:
            geometry (Geometry): 2D 2 variables per node geometry
            E (int or float or list): Young Moduli. If number, all element will have the same young moduli. If list, each position will be the element young moduli, so len(E) == len(self.elements)
            v (int or float or list): Poisson ratio. If number, all element will have the same Poisson ratio. If list, each position will be the element Poisson ratio, so len(v) == len(self.elements)
            rho (int or float or list, optional): Density. If not given, mass matrix will not be calculated. Defaults to None.
            fx (function, optional): Function fx, if fx is constant you can use fx = lambda x: [value]. Defaults to lambda x:0.
            fy (function, optional): Function fy, if fy is constant you can use fy = lambda x: [value]. Defaults to lambda x:0.
    """

    def __init__(self, geometry: Geometry, E: Tuple[float, list], v: Tuple[float, list], rho: Tuple[float, list] = None, fx: Callable = lambda x: 0, fy: Callable = lambda x: 0, **kargs) -> None:
        """Create a Plain Strain problem

        Args:
                geometry (Geometry): 2D 2 variables per node geometry
                E (int or float or list): Young Moduli. If number, all element will have the same young moduli. If list, each position will be the element young moduli, so len(E) == len(self.elements)
                v (int or float or list): Poisson ratio. If number, all element will have the same Poisson ratio. If list, each position will be the element Poisson ratio, so len(v) == len(self.elements)
                rho (int or float or list, optional): Density. If not given, mass matrix will not be calculated. Defaults to None.
                fx (function, optional): Function fx, if fx is constant you can use fx = lambda x: [value]. Defaults to lambda x:0.
                fy (function, optional): Function fy, if fy is constant you can use fy = lambda x: [value]. Defaults to lambda x:0.
        """

        PlaneStress.__init__(
            self, geometry, E, v, 1, rho, fx, fy, **kargs)
        self.C11 = []
        self.C22 = []
        self.C12 = []
        self.C66 = []
        for i in range(len(self.geometry.elements)):
            C11 = self.E1[i]*(1-self.v12[i])/(1+self.v12[i])/(1-2*self.v12[i])
            C12 = self.E1[i]*(self.v12[i])/(1+self.v12[i])/(1-2*self.v12[i])
            C66 = self.E1[i] / 2 / (1 + self.v12[i])
            self.C11.append(C11)
            self.C22.append(C11)
            self.C12.append(C12)
            self.C66.append(C66)
        self.name = 'Plane Strain Isotropic'


class PlaneStrainSparse(PlaneStressSparse):
    """Create a Plain Strain problem using sparse matrix

    Args:
            geometry (Geometry): 2D 2 variables per node geometry with fast elements
            E (int or float or list): Young Moduli. If number, all element will have the same young moduli. If list, each position will be the element young moduli, so len(E) == len(self.elements)
            v (int or float or list): Poisson ratio. If number, all element will have the same Poisson ratio. If list, each position will be the element Poisson ratio, so len(v) == len(self.elements)
            rho (int or float or list, optional): Density. If not given, mass matrix will not be calculated. Defaults to None.
            fx (function, optional): Function fx, if fx is constant you can use fx = lambda x: [value]. Defaults to lambda x:0.
            fy (function, optional): Function fy, if fy is constant you can use fy = lambda x: [value]. Defaults to lambda x:0.
    """

    def __init__(self, geometry: Geometry, E: Tuple[float, list], v: Tuple[float, list], rho: Tuple[float, list] = None, fx: Callable = lambda x: 0, fy: Callable = lambda x: 0, **kargs) -> None:
        """Create a Plain Strain problem using sparse matrix

        Args:
                geometry (Geometry): 2D 2 variables per node geometry with fast elements
                E (int or float or list): Young Moduli. If number, all element will have the same young moduli. If list, each position will be the element young moduli, so len(E) == len(self.elements)
                v (int or float or list): Poisson ratio. If number, all element will have the same Poisson ratio. If list, each position will be the element Poisson ratio, so len(v) == len(self.elements)
                rho (int or float or list, optional): Density. If not given, mass matrix will not be calculated. Defaults to None.
                fx (function, optional): Function fx, if fx is constant you can use fx = lambda x: [value]. Defaults to lambda x:0.
                fy (function, optional): Function fy, if fy is constant you can use fy = lambda x: [value]. Defaults to lambda x:0.
        """

        PlaneStressSparse.__init__(
            self, geometry, E, v, 1, rho, fx, fy, **kargs)
        self.C11 = []
        self.C22 = []
        self.C12 = []
        self.C66 = []
        for i in range(len(self.geometry.elements)):
            C11 = self.E1[i]*(1-self.v12[i])/(1+self.v12[i])/(1-2*self.v12[i])
            C12 = self.E1[i]*(self.v12[i])/(1+self.v12[i])/(1-2*self.v12[i])
            C66 = self.E1[i] / 2 / (1 + self.v12[i])
            self.C11.append(C11)
            self.C22.append(C11)
            self.C12.append(C12)
            self.C66.append(C66)
        self.name = 'Plane Strain Isotropic sparse'


class PlaneStressNonLocalSparse(PlaneStressSparse):
    """Create a Plain Stress nonlocal problem using sparse matrices and the Pisano 2006 formulation.

        Args:
                geometry (Geometry): 2D 2 variables per node geometry
                E (int or float or list): Young Moduli. If number, all element will have the same young moduli. If list, each position will be the element young moduli, so len(E) == len(self.elements)
                v (int or float or list): Poisson ratio. If number, all element will have the same Poisson ratio. If list, each position will be the element Poisson ratio, so len(v) == len(self.elements)
                t (int or float or list): Element thickness. If number, all element will have the same thickness. If list, each position will be the element thickness, so len(t) == len(self.elements)
                l (float): Internal lenght
                z1 (float): z1 factor
                Lr (float): Influence distance Lr
                af (Callable): Atenuation function
                rho (int or float or list, optional): Density. If not given, mass matrix will not be calculated. Defaults to None.
                fx (function, optional): Function fx, if fx is constant you can use fx = lambda x: [value]. Defaults to lambda x:0.
                fy (function, optional): Function fy, if fy is constant you can use fy = lambda x: [value]. Defaults to lambda x:0.
        """

    def __init__(self, geometry: Geometry, E: Tuple[float, list], v: Tuple[float, list], t: Tuple[float, list], l: float, z1: float, Lr: float, af: Callable, rho: Tuple[float, list] = None, fx: Callable = lambda x: 0, fy: Callable = lambda x: 0, **kargs) -> None:
        """Create a Plain Stress nonlocal problem using sparse matrices and the Pisano 2006 formulation.

        Args:
                geometry (Geometry): 2D 2 variables per node geometry
                E (int or float or list): Young Moduli. If number, all element will have the same young moduli. If list, each position will be the element young moduli, so len(E) == len(self.elements)
                v (int or float or list): Poisson ratio. If number, all element will have the same Poisson ratio. If list, each position will be the element Poisson ratio, so len(v) == len(self.elements)
                t (int or float or list): Element thickness. If number, all element will have the same thickness. If list, each position will be the element thickness, so len(t) == len(self.elements)
                l (float): Internal lenght
                z1 (float): z1 factor
                Lr (float): Influence distance Lr
                af (Callable): Atenuation function
                rho (int or float or list, optional): Density. If not given, mass matrix will not be calculated. Defaults to None.
                fx (function, optional): Function fx, if fx is constant you can use fx = lambda x: [value]. Defaults to lambda x:0.
                fy (function, optional): Function fy, if fy is constant you can use fy = lambda x: [value]. Defaults to lambda x:0.
        """

        self.l = l
        self.Lr = Lr
        self.af = af
        self.z1 = z1
        self.l0 = 0.5/np.pi/l/l/t
        self.z2 = 1.0-self.z1

        PlaneStressSparse.__init__(
            self, geometry, E, v, t, rho, fx, fy, **kargs)
        nonlocals = self.geometry.detectNonLocal(Lr)
        for e, dno in zip(self.elements, nonlocals):
            e.enl = dno
        self.name = 'Plane Stress Isotropic non local sparse'

        self.KL = sparse.lil_matrix((self.ngdl, self.ngdl))
        self.KNL = sparse.lil_matrix((self.ngdl, self.ngdl))
        if self.calculateMass:
            self.M = sparse.lil_matrix((self.ngdl, self.ngdl))

    def elementMatrices(self) -> None:
        """Calculate the elements matrices
        """
        for ee in tqdm(range(len(self.elements)), unit='Local'):
            self.elementMatrix(ee)

    def elementMatrix(self, ee: 'Element') -> None:
        """Calculates a single element local and nonlocal matrices

        Args:
            ee (Element): Element to be calculated
        """
        e = self.elements[ee]
        m = len(e.gdl.T)

        Fux = np.zeros([m, 1])
        Fvx = np.zeros([m, 1])

        _x, _p = e._x, e._p
        detjac = e.detjac
        dpx = e.dpx

        c11 = self.C11[ee]
        c12 = self.C12[ee]
        c22 = self.C22[ee]
        c66 = self.C66[ee]
        C = np.array([
            [c11, c12, 0.0],
            [c12, c22, 0.0],
            [0.0, 0.0, c66]])

        Fe = np.zeros([2*m, 1])
        Ke = np.zeros([2*m, 2*m])
        if self.calculateMass:
            Me = np.zeros([2*m, 2*m])

        o = [0.0]*m
        for k in range(len(e.Z)):  # Iterate over gauss points on domain
            B = np.array([
                [*dpx[k, 0, :], *o],
                [*o, *dpx[k, 1, :]],
                [*dpx[k, 1, :], *dpx[k, 0, :]]])

            P = np.array([
                [*_p[k], *o],
                [*o, *_p[k]]])
            Ke += self.t[ee]*(B.T@C@B)*detjac[k]*e.W[k]
            if self.calculateMass:
                Me += self.rho[ee]*self.t[ee]*(P.T@P)*detjac[k]*e.W[k]
            Fe += self.t[ee]*(P.T@np.array([[self.fx(_x[k])],
                                            [self.fy(_x[k])]]))*detjac[k]*e.W[k]

        if e.intBorders:  # Cambiar esto a la notación matricial
            for j in range(len(e.borders)):
                border = e.borders[j]
                if (len(border.properties['load_x']) + len(border.properties['load_y'])):
                    _x, _p = e.T(e.Tj[j](border.Z.T))
                    _s = border.TS(border.Z.T)
                    detjac = border.coords[-1, 0]*0.5
                    for i in range(m):
                        for fx in border.properties['load_x']:
                            for k in range(len(border.Z)):
                                Fux[i, 0] += fx(_s[k, 0])*_p[k, i] * \
                                    detjac*border.W[k]
                        for fy in border.properties['load_y']:
                            for k in range(len(border.Z)):
                                Fvx[i, 0] += fy(_s[k, 0])*_p[k, i] * \
                                    detjac*border.W[k]
        subm = np.linspace(0, 2*m-1, 2*m).reshape([2, m]).astype(int)
        Fe[np.ix_(subm[0])] += Fux
        Fe[np.ix_(subm[1])] += Fvx
        self.F[np.ix_(e.gdlm)] += Fe

        self.KL[np.ix_(e.gdlm, e.gdlm)] += Ke

        if self.calculateMass:
            self.M[np.ix_(e.gdlm, e.gdlm)] += Me

        # e.knls = []
        for inl in tqdm(e.enl, unit=' Nolocal'):
            enl = self.elements[inl]
            mnl = len(enl.gdl.T)
            onl = [0.0]*mnl
            Knl = np.zeros([2*m, 2*mnl])
            _xnl = enl._x
            detjacnl = enl.detjac
            dpxnl = enl.dpx

            for k in range(len(e.Z)):
                B = np.array([
                    [*dpx[k, 0, :], *o],
                    [*o, *dpx[k, 1, :]],
                    [*dpx[k, 1, :], *dpx[k, 0, :]]])

                for knl in range(len(enl.Z)):
                    ro = np.linalg.norm(_x[k]-_xnl[knl])/self.l
                    azn = self.af(self.l0, ro)
                    Bnl = np.array([
                        [*dpxnl[knl, 0, :], *onl],
                        [*onl, *dpxnl[knl, 1, :]],
                        [*dpxnl[knl, 1, :], *dpxnl[knl, 0, :]]])

                    Knl += self.t[ee]*self.t[inl]*azn*(Bnl.T@C@B)*detjac[k] * \
                        e.W[k]*detjacnl[knl]*enl.W[knl]
            # e.knls.append(Knl)
            self.KNL[np.ix_(e.gdlm, enl.gdlm)] += Knl

    def profile(self, p0: list, p1: list, n: float = 100) -> None:
        """Generate a profile between selected points

        Args:
            p0 (list): start point of the profile [x0,y0]
            p1 (list): end point of the profile [xf,yf]
            n (int, optional): NUmber of samples for graph. Defaults to 100.

        """

        _x = np.linspace(p0[0], p1[0], n)
        _y = np.linspace(p0[1], p1[1], n)
        X = np.array([_x, _y])
        U = []
        U1 = []
        U2 = []
        U3 = []
        _X = []
        def dist(X): return np.sqrt((p0[0]-X[0])**2+(p0[1]-X[1])**2)
        for i in range(n):
            for _, e in enumerate(self.elements):
                if e.isInside(X.T[i]):
                    z = e.inverseMapping(np.array([X.T[i]]).T)
                    _, u, du = e.giveSolutionPoint(z, True)
                    U += [u.tolist()]
                    U1 += (du[:, 0, 0]).tolist()
                    U2 += (du[:, 1, 1]).tolist()
                    U3 += (1/2*(du[:, 0, 1]+du[:, 1, 0])).tolist()
                    _X.append(dist(X.T[i]))
                    break
        fig = plt.figure()
        ax = fig.add_subplot(1, 3, 1)
        ax.plot(_X, U1, color='black')
        ax.grid()
        ax.set_xlabel('d')
        ax.set_ylabel(r'$\varepsilon_{xx}$')
        ax = fig.add_subplot(1, 3, 2)
        ax.plot(_X, U2, color='black')
        ax.grid()
        ax.set_xlabel('d')
        ax.set_ylabel(r'$\varepsilon_{yy}$')
        ax = fig.add_subplot(1, 3, 3)
        ax.plot(_X, U3, color='black')
        ax.grid()
        ax.set_xlabel('d')
        ax.set_ylabel(r'$\varepsilon_{xy}$')
        return _X, U1, U2, U3, U

    def ensembling(self) -> None:
        """Creation of the system sparse matrix. Force vector is ensembled in integration method
        """
        logging.info('Ensembling equation system...')
        self.K = self.KL*self.z1 + self.KNL*self.z2
        if self.calculateMass:
            self.M = self.M.tocsr()
        logging.info('Done!')
