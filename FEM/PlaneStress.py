"""2D Elasticity: Plane Stress
"""


from .Core import Core, Geometry
from tqdm import tqdm
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
from typing import Tuple, Callable


class PlaneStress(Core):
    """Create a Plain Stress problem

    Args:
            geometry (Geometry): 2D 2 variables per node geometry
            E (int or float or list): Young Moduli. If number, all element will have the same young moduli. If list, each position will be the element young moduli, so len(E) == len(self.elements)
            v (int or float or list): Poisson ratio. If number, all element will have the same Poisson ratio. If list, each position will be the element Poisson ratio, so len(v) == len(self.elements)
            t (int or float or list): Element thickness. If number, all element will have the same thickness. If list, each position will be the element thickness, so len(t) == len(self.elements)
            fx (function, optional): Function fx, if fx is constant you can use fx = lambda x: [value]. Defaults to lambda x:0.
            fy (function, optional): Function fy, if fy is constant you can use fy = lambda x: [value]. Defaults to lambda x:0.
    """

    def __init__(self, geometry: Geometry, E: Tuple[float, list], v: Tuple[float, list], t: Tuple[float, list], fx: Callable = lambda x: 0, fy: Callable = lambda x: 0) -> None:
        """Create a Plain Stress problem

        Args:
                geometry (Geometry): 2D 2 variables per node geometry
                E (int or float or list): Young Moduli. If number, all element will have the same young moduli. If list, each position will be the element young moduli, so len(E) == len(self.elements)
                v (int or float or list): Poisson ratio. If number, all element will have the same Poisson ratio. If list, each position will be the element Poisson ratio, so len(v) == len(self.elements)
                t (int or float or list): Element thickness. If number, all element will have the same thickness. If list, each position will be the element thickness, so len(t) == len(self.elements)
                fx (function, optional): Function fx, if fx is constant you can use fx = lambda x: [value]. Defaults to lambda x:0.
                fy (function, optional): Function fy, if fy is constant you can use fy = lambda x: [value]. Defaults to lambda x:0.
        """
        if type(t) == float or type(t) == int:
            t = [t]*len(geometry.elements)
        if type(E) == float or type(E) == int:
            E = [E]*len(geometry.elements)
        if type(v) == float or type(v) == int:
            v = [v]*len(geometry.elements)
        self.t = t
        self.E = E
        self.v = v
        self.C11 = []
        self.C12 = []
        self.C66 = []
        self.fx = fx
        self.fy = fy
        for i in range(len(self.E)):
            C11 = self.E[i] / (1 - self.v[i] ** 2)
            C12 = self.v[i] * C11
            C66 = self.E[i] / 2 / (1 + self.v[i])
            self.C11.append(C11)
            self.C12.append(C12)
            self.C66.append(C66)
        if geometry.nvn == 1:
            print(
                'Border conditions lost, please usea a geometry with 2 variables per node (nvn=2)\nRegenerating Geoemtry...')
            geometry.nvn = 2
            geometry.cbe = []
            geometry.cbn = []
            geometry.initialize()
        Core.__init__(self, geometry)

    def elementMatrices(self) -> None:
        """Calculate the element matrices usign Reddy's (2005) finite element model
        """

        for ee, e in enumerate(tqdm(self.elements, unit='Element')):
            m = len(e.gdl.T)
            Kuu = np.zeros([m, m])
            Kuv = np.zeros([m, m])
            Kvu = np.zeros([m, m])
            Kvv = np.zeros([m, m])
            Fu = np.zeros([m, 1])
            Fv = np.zeros([m, 1])
            Fux = np.zeros([m, 1])
            Fvx = np.zeros([m, 1])
            # Gauss points in global coordinates and Shape functions evaluated in gauss points
            _x, _p = e.T(e.Z.T)
            # Jacobian evaluated in gauss points and shape functions derivatives in natural coordinates
            jac, dpz = e.J(e.Z.T)
            detjac = np.linalg.det(jac)
            _j = np.linalg.inv(jac)  # Jacobian inverse
            dpx = _j @ dpz  # Shape function derivatives in global coordinates
            for i in range(m):  # self part must be vectorized
                for j in range(m):
                    for k in range(len(e.Z)):  # Iterate over gauss points on domain
                        Kuu[i, j] += (self.C11[ee]*dpx[k, 0, i]*dpx[k, 0, j] +
                                      self.C66[ee]*dpx[k, 1, i]*dpx[k, 1, j])*detjac[k]*e.W[k]
                        Kuv[i, j] += (self.C12[ee]*dpx[k, 0, i]*dpx[k, 1, j] +
                                      self.C66[ee]*dpx[k, 1, i]*dpx[k, 0, j])*detjac[k]*e.W[k]
                        Kvu[i, j] += (self.C12[ee]*dpx[k, 1, i]*dpx[k, 0, j] +
                                      self.C66[ee]*dpx[k, 0, i]*dpx[k, 1, j])*detjac[k]*e.W[k]
                        Kvv[i, j] += (self.C11[ee]*dpx[k, 1, i]*dpx[k, 1, j] +
                                      self.C66[ee]*dpx[k, 0, i]*dpx[k, 0, j])*detjac[k]*e.W[k]
                for k in range(len(e.Z)):  # Iterate over gauss points on domain

                    Fu[i][0] += _p[k, i] * self.fx(_x[k])*detjac[k]*e.W[k]
                    Fv[i][0] += _p[k, i] * self.fy(_x[k])*detjac[k]*e.W[k]

            if e.intBorders:
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
            e.Fe[np.ix_(subm[0])] += self.t[ee]*(Fu)+Fux
            e.Fe[np.ix_(subm[1])] += self.t[ee]*(Fv)+Fvx
            e.Ke[np.ix_(subm[0], subm[0])] += self.t[ee]*(Kuu)
            e.Ke[np.ix_(subm[0], subm[1])] += self.t[ee]*(Kuv)
            e.Ke[np.ix_(subm[1], subm[0])] += self.t[ee]*(Kvu)
            e.Ke[np.ix_(subm[1], subm[1])] += self.t[ee]*(Kvv)

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
            U2 += (self.C12[ee]*du[:, 0, 0]+self.C11[ee]*du[:, 1, 1]).tolist()
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
        def fmt(x): return format(x, '.3f')

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
        if not mask == None:
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
                # TODO Arreglar calculo de esfuerzos para PlaneStrain
                sx = (self.C11[ee]*du[:, 0, 0] +
                      self.C12[ee]*du[:, 1, 1]).tolist()
                sy = (self.C12[ee]*du[:, 0, 0] +
                      self.C11[ee]*du[:, 1, 1]).tolist()
                sxy = (self.C66[ee]*(du[:, 0, 1]+du[:, 1, 0])).tolist()
                return sx, sy, sxy

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
        U1 = []
        U2 = []
        U3 = []
        _X = []
        def dist(X): return np.sqrt((p0[0]-X[0])**2+(p0[1]-X[1])**2)
        for i in range(n):
            for ee, e in enumerate(self.elements):
                if e.isInside(X.T[i]):
                    z = e.inverseMapping(np.array([X.T[i]]).T)
                    _, _, du = e.giveSolutionPoint(z, True)
                    # TODO Arreglar calculo de esfuerzos para PlaneStrain
                    U1 += (self.C11[ee]*du[:, 0, 0] +
                           self.C12[ee]*du[:, 1, 1]).tolist()
                    U2 += (self.C12[ee]*du[:, 0, 0] +
                           self.C11[ee]*du[:, 1, 1]).tolist()
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
