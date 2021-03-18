"""Plane Stress Class
"""


from .Core import *
from tqdm import tqdm
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec


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

    def __init__(self, geometry, E, v, t, fx=lambda x: 0, fy=lambda x: 0):
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
            C12 = self.v[i] * self.E[i] / (1 - self.v[i] ** 2)
            C66 = self.E[i] / 2 / (1 + self.v[i])
            self.C11.append(C11)
            self.C12.append(C12)
            self.C66.append(C66)
        if geometry.nvn == 1:
            print(
                'Border conditions lost, please usea a geometry with 2 variables per node (nvn=2)')
            geometry.nvn = 2
            geometry.cbe = []
            geometry.cbn = []
            geometry.initialize()
        Core.__init__(self, geometry)

    def elementMatrices(self):
        """Calculate the element matrices usign Reddy's (2005) finite element model
        """
# EKuu = lambda z, n: (C11 * dfdx(z, n, i) * dfdx(z, n, j) + C66 * dfdy(z, n, i) * dfdy(z, n, j)) * np.linalg.det(J(z, n))
# EKuv = lambda z, n: (C12 * dfdx(z, n, i) * dfdy(z, n, j) + C66 * dfdy(z, n, i) * dfdx(z, n, j)) * np.linalg.det(J(z, n))
# EKvu = lambda z, n: (C12 * dfdy(z, n, i) * dfdx(z, n, j) + C66 * dfdx(z, n, i) * dfdy(z, n, j)) * np.linalg.det(J(z, n))
# EKvv = lambda z, n: (C11 * dfdy(z, n, i) * dfdy(z, n, j) + C66 * dfdx(z, n, i) * dfdx(z, n,j)) * np.linalg.det(J(z, n))
        ee = 0
        for e in tqdm(self.elements, unit='Element'):
            m = len(e.gdl.T)
            Kuu = np.zeros([m, m])
            Kuv = np.zeros([m, m])
            Kvu = np.zeros([m, m])
            Kvv = np.zeros([m, m])
            Fu = np.zeros([m, 1])
            Fv = np.zeros([m, 1])
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
                        Kuu[i, j] += self.t[ee]*(self.C11[ee]*dpx[k, 0, i]*dpx[k, 0, j] +
                                                 self.C66[ee]*dpx[k, 1, i]*dpx[k, 1, j])*detjac[k]*e.W[k]
                        Kuv[i, j] += self.t[ee]*(self.C12[ee]*dpx[k, 0, i]*dpx[k, 1, j] +
                                                 self.C66[ee]*dpx[k, 1, i]*dpx[k, 0, j])*detjac[k]*e.W[k]
                        Kvu[i, j] += self.t[ee]*(self.C12[ee]*dpx[k, 1, i]*dpx[k, 0, j] +
                                                 self.C66[ee]*dpx[k, 0, i]*dpx[k, 1, j])*detjac[k]*e.W[k]
                        Kvv[i, j] += self.t[ee]*(self.C11[ee]*dpx[k, 1, i]*dpx[k, 1, j] +
                                                 self.C66[ee]*dpx[k, 0, i]*dpx[k, 0, j])*detjac[k]*e.W[k]
                for k in range(len(e.Z)):  # Iterate over gauss points on domain
                    Fu[i][0] += self.t[ee]*_p[k, i] * \
                        self.fx(_x[k])*detjac[k]*e.W[k]
                    Fv[i][0] += self.t[ee]*_p[k, i] * \
                        self.fy(_x[k])*detjac[k]*e.W[k]
            subm = np.linspace(0, 2*m-1, 2*m).reshape([2, m]).astype(int)
            e.Fe[np.ix_(subm[0])] += Fu
            e.Fe[np.ix_(subm[1])] += Fv
            e.Ke[np.ix_(subm[0], subm[0])] += Kuu
            e.Ke[np.ix_(subm[0], subm[1])] += Kuv
            e.Ke[np.ix_(subm[1], subm[0])] += Kvu
            e.Ke[np.ix_(subm[1], subm[1])] += Kvv
            ee += 1
            # e.Fe[:,0] = 2*self.G*self._phi*detjac@_p
            # e.Ke = (np.transpose(dpx,axes=[0,2,1]) @ dpx).T @ detjac

    def postProcess(self, mult=1000):
        """Generate the stress surfaces and displacement fields for the geometry

        Args:
                mult (int, optional): Factor for displacements. Defaults to 1000.
        """
        X = []
        Y = []
        U1 = []
        U2 = []
        U3 = []
        fig = plt.figure()

        gs = gridspec.GridSpec(3, 3)

        ax1 = fig.add_subplot(gs[0, 0])
        ax2 = fig.add_subplot(gs[0, 1])
        ax3 = fig.add_subplot(gs[0, 2])
        ax5 = fig.add_subplot(gs[1:, :])
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
        surf = ax1.tricontourf(X, Y, U1, cmap='magma', levels=20)
        plt.colorbar(surf, ax=ax1)
        ax1.set_title(r'$\sigma_{xx}$')

        surf = ax2.tricontourf(X, Y, U2, cmap='magma', levels=20)
        plt.colorbar(surf, ax=ax2)
        ax2.set_title(r'$\sigma_{yy}$')

        surf = ax3.tricontourf(X, Y, U3, cmap='magma', levels=20)
        plt.colorbar(surf, ax=ax3)
        ax3.set_title(r'$\sigma_{xy}$')
        mask = self.geometry.mask
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

    def profile(self, p0, p1, n=100):
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
