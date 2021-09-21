"""Euler Bernoulli Beam implementation [WIP]
"""


from .Elements.E1D.EulerBernoulliElement import EulerBernoulliElement
from .Core import Core, Geometry
from tqdm import tqdm
import numpy as np
import matplotlib.pyplot as plt


class EulerBernoulliBeam(Core):
    """Creates a Euler Bernoulli beam problem

    Args:
        geometry (Geometry): 1D 2 variables per node problem geometry. Geometry must have Euler Bernoulli elements.
        EI (float): Young's moduli multiplied by second moment of area (inertia).
        cf (float, optional): Soil coeficient. Defaults to 0.
    """

    def __init__(self, geometry: Geometry, EI: float, f: float = 0) -> None:
        """Creates a Euler Bernoulli beam problem

        Args:
            geometry (Geometry): 1D 2 variables per node problem geometry. Geometry must have Lineal elements.
            EI (float): Young's moduli multiplied by second moment of area (inertia).
            cf (float, optional): Soil coeficient. Defaults to 0.
        """
        self.a = EI
        self.f = f
        if isinstance(EI, float):
            self.a = lambda x: EI
        if isinstance(f, float):
            self.f = lambda x: f
        if geometry.nvn == 1:
            print(
                'Border conditions lost, please usea a geometry with 2 variables per node (nvn=2)')
        Core.__init__(self, geometry)
        for i in range(len(self.elements)):
            self.elements[i] = EulerBernoulliElement(
                self.elements[i].coords, self.elements[i].gdl)

    def elementMatrices(self) -> None:
        """Calculate the element matrices usign Guass Legendre quadrature.
        """
        for e in tqdm(self.elements, unit='Element'):
            _x, _p = e.T(e.Z.T)
            _h = e.hermit(e.Z.T)
            jac, dpz = e.J(e.Z.T)
            detjac = np.linalg.det(jac)
            # _j = np.linalg.inv(jac)
            # dpx = _j @ dpz
            _dh = e.dhermit(e.Z.T)
            for i in range(e.n):
                for j in range(e.n):
                    for k in range(len(e.Z)):
                        # + self.c(_x[k])*_p[k][i]*_p[k][j]
                        e.Ke[i, j] += (self.a(_x[k])*_dh[1][i][k]
                                       * _dh[1][j][k])*detjac[k]*e.W[k]
                for k in range(len(e.Z)):
                    e.Fe[i][0] += self.f(_x[k])*_h[k][i]*detjac[k]*e.W[k]

    def postProcess(self) -> None:
        """Post process the solution. Shows graphs of displacement, rotation, shear and moment.
        """
        X = []
        U1 = []
        U2 = []
        U3 = []
        U4 = []
        fig = plt.figure()
        ax1 = fig.add_subplot(2, 2, 1)
        ax2 = fig.add_subplot(2, 2, 2)
        ax3 = fig.add_subplot(2, 2, 3)
        ax4 = fig.add_subplot(2, 2, 4)
        for e in self.elements:
            _x, _u, du = e.giveSolution(True)
            X += _x.T[0].tolist()
            U1 += _u.tolist()
            U2 += (du[:, 0]).tolist()
            U3 += (du[:, 1]).tolist()
            U4 += (du[:, 2]).tolist()
        ax1.plot(X, U1)
        ax1.grid()
        ax2.plot(X, U2)
        ax2.grid()

        ax3.plot(X, U3)
        ax3.grid()

        ax4.plot(X, U4)
        ax4.grid()
        ax1.set_title(r'$U(x)$')
        ax2.set_title(r'$\frac{dU}{dx}$')
        ax3.set_title(r'$\frac{d^2U}{dx^2}$')
        ax4.set_title(r'$\frac{d^3U}{dx^3}$')
