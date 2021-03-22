"""Euler Bernoulli Beam implementation [WIP]
"""


from .Core import *
from tqdm import tqdm
import numpy as np
import matplotlib.pyplot as plt


class EulerBernoulliBeam(Core):
    def __init__(self, geometry: Geometry, EI: float, cf: float = 0) -> None:
        """Creates a Euler Bernoulli beam problem

        Args:
            geometry (Geometry): 1D 2 variables per node problem geometry. Geometry must have Euler Bernoulli elements.
            EI (float): Young's moduli multiplied by second moment of area (inertia).
            cf (float, optional): Soil coeficient. Defaults to 0.
        """
        if geometry.nvn == 1:
            print(
                'Border conditions lost, please usea a geometry with 2 variables per node (nvn=2)')
        Core.__init__(self, geometry)

    def elementMatrices(self) -> None:
        """Calculate the element matrices usign Guass Legendre quadrature.
        """
        # for e in tqdm(self.elements, unit='Element'):
        #     _x, _p = e.T(e.Z.T)
        #     jac, dpz = e.J(e.Z.T)
        #     detjac = np.linalg.det(jac)
        #     _j = np.linalg.inv(jac)
        #     dpx = _j @ dpz
        #     for i in range(e.n):
        #         for j in range(e.n):
        #             for k in range(len(e.Z)):
        #                 e.Ke[i, j] += (self.a(_x[k])*dpx[k][0][i]*dpx[k][0]
        #                                [j] + self.c(_x[k])*_p[k][i]*_p[k][j])*detjac[k]*e.W[k]
        #         for k in range(len(e.Z)):
        #             e.Fe[i][0] += self.f(_x[k])*_p[k][i]*detjac[k]*e.W[k]

    def postProcess(self) -> None:
        """Post process the solution. Shows graphs of displacement, rotation, shear and moment.
        """
        # X = []
        # U1 = []
        # U2 = []
        # fig = plt.figure()
        # ax1 = fig.add_subplot(1, 2, 1)
        # ax2 = fig.add_subplot(1, 2, 2)
        # for e in self.elements:
        #     _x, _u, du = e.giveSolution(True)
        #     X += _x.T[0].tolist()
        #     U1 += _u[0].tolist()
        #     U2 += (du[:, 0, 0]).tolist()
        # surf = ax1.plot(X, U1)
        # surf = ax2.plot(X, U2)
        # ax1.grid()
        # ax2.grid()
        # ax1.set_title(r'$U(x)$')
        # ax2.set_title(r'$\frac{dU}{dx}$')
