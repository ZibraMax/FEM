"""Solve the torsion Poisson equation for a 2D domain.
"""


from .Core import *
from tqdm import tqdm
import numpy as np
import matplotlib.pyplot as plt


class Torsion2D(Core):
    """Create a torsional finite element problem

    Args:
            geometry (Geometry): 2D 1 variable per node geometry
            G (float): Shear moduli of elements
            phi (float): Rotation angle in radians
    """

    def __init__(self, geometry: Geometry, G: float, phi: float) -> None:
        """Create a torsional finite element problem

        Args:
                geometry (Geometry): 2D 1 variable per node geometry
                G (float): Shear moduli of elements
                phi (float): Rotation angle in radians
        """

        self.G = G
        self._phi = phi
        geometry.cbeAllBorders(0)
        Core.__init__(self, geometry)

    def elementMatrices(self) -> None:
        """Calculate the element matrices usign Reddy's (2005) finite element model
        """

        for e in tqdm(self.elements, unit='Element'):
            # Gauss points in global coordinates and Shape functions evaluated in gauss points
            _x, _p = e.T(e.Z.T)
            # Jacobian evaluated in gauss points and shape functions derivatives in natural coordinates
            jac, dpz = e.J(e.Z.T)
            detjac = np.linalg.det(jac)*e.W
            _j = np.linalg.inv(jac)  # Jacobian inverse
            dpx = _j @ dpz  # Shape function derivatives in global coordinates
            # for k in range(len(_x)): #Iterate over gauss points on domain
            # 	for i in range(e.n): #self part must be vectorized
            # 		for j in range(e.n):
            # 			e.Ke[i,j] += (dpx[k][0][i]*dpx[k][0][j] + dpx[k][1][i]*dpx[k][1][j])*detjac[k]*e.W[k]
            # 		e.Fe[i][0] += 2*self.G*self._phi*_p[k][i]*detjac[k]*e.W[k]
            e.Fe[:, 0] = 2*self.G*self._phi*detjac@_p
            e.Ke = (np.transpose(dpx, axes=[0, 2, 1]) @ dpx).T @ detjac

    def postProcess(self) -> None:
        """Create graphs for stress function and derivatives.
        """

        X = []
        Y = []
        U1 = []
        U2 = []
        U3 = []
        U4 = []
        fig = plt.figure()
        ax1 = fig.add_subplot(2, 2, 1, projection='3d')
        ax2 = fig.add_subplot(2, 2, 2)
        ax3 = fig.add_subplot(2, 2, 3)
        ax4 = fig.add_subplot(2, 2, 4)
        for e in tqdm(self.elements, unit='Element'):
            _x, _u, du = e.giveSolution(True)
            X += _x.T[0].tolist()
            Y += _x.T[1].tolist()
            U2 += du[:, 0, 0].tolist()
            U3 += du[:, 0, 1].tolist()
            U1 += _u[0].tolist()
            U4 += np.sqrt(du[:, 0, 0]**2 + du[:, 0, 1]**2).tolist()
        surf = ax1.plot_trisurf(X, Y, U1, cmap='magma')
        fig.colorbar(surf, ax=ax1)

        surf = ax2.tricontourf(X, Y, U2, cmap='magma')

        fig.colorbar(surf, ax=ax2)

        surf = ax3.tricontourf(X, Y, U3, cmap='magma')

        fig.colorbar(surf, ax=ax3)

        surf = ax4.tricontourf(X, Y, U4, cmap='magma')

        fig.colorbar(surf, ax=ax4)
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
            ax2.fill(Xs, Ys, color='white', zorder=30)
            ax3.fill(Xs, Ys, color='white', zorder=30)
            ax4.fill(Xs, Ys, color='white', zorder=30)
