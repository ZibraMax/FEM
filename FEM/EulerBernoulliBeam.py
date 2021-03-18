from .Core import *
from tqdm import tqdm
import numpy as np
import matplotlib.pyplot as plt


class EulerBernoulliBeam(Core):
    def __init__(self, geometry, EI, cf):
        if geometry.nvn == 1:
            print(
                'Border conditions lost, please usea a geometry with 2 variables per node (nvn=2)')
        Core.__init__(self, geometry)

    def elementMatrices(self):
        for e in tqdm(self.elements, unit='Element'):
            # Gauss points in global coordinates and Shape functions evaluated in gauss points
            _x, _p = e.T(e.Z.T)
            # Jacobian evaluated in gauss points and shape functions derivatives in natural coordinates
            jac, dpz = e.J(e.Z.T)
            detjac = np.linalg.det(jac)
            _j = np.linalg.inv(jac)  # Jacobian inverse
            dpx = _j @ dpz  # Shape function derivatives in global coordinates
            for i in range(e.n):  # self part must be vectorized
                for j in range(e.n):
                    for k in range(len(e.Z)):  # Iterate over gauss points on domain
                        e.Ke[i, j] += (self.a(_x[k])*dpx[k][0][i]*dpx[k][0]
                                       [j] + self.c(_x[k])*_p[k][i]*_p[k][j])*detjac[k]*e.W[k]
                for k in range(len(e.Z)):  # Iterate over gauss points on domain
                    e.Fe[i][0] += self.f(_x[k])*_p[k][i]*detjac[k]*e.W[k]
            # e.Fe[:,0] = 2*self.G*self._phi*detjac@_p
            # e.Ke = (np.transpose(dpx,axes=[0,2,1]) @ dpx).T @ detjac

    def postProcess(self):
        X = []
        U1 = []
        U2 = []
        fig = plt.figure()
        ax1 = fig.add_subplot(1, 2, 1)
        ax2 = fig.add_subplot(1, 2, 2)
        for e in self.elements:
            _x, _u, du = e.giveSolution(True)
            X += _x.T[0].tolist()
            U1 += _u[0].tolist()
            U2 += (du[:, 0, 0]).tolist()
        surf = ax1.plot(X, U1)
        surf = ax2.plot(X, U2)
        ax1.grid()
        ax2.grid()
        ax1.set_title(r'$U(x)$')
        ax2.set_title(r'$\frac{dU}{dx}$')
