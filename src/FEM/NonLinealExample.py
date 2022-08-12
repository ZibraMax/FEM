"""Solve a 1D 1 variable per node non lineal Finite Element problem
"""


from re import M
from .Core import Core, Geometry
from .Solvers import NoLineal
from tqdm import tqdm
import numpy as np
import matplotlib.pyplot as plt
from typing import Callable


class NonLinealSimpleEquation(Core):
    """Creates a nonlineal 1D equation with the form:

        .. math::
            -\\frac{d}{dx}\\left(a(x)u\\frac{du}{dx}\\right)=f(x)

        Args:
            geometry (Geometry): Input lineal geometry
            a (Callable): Function a
            f (Callable): Function f
        """

    def __init__(self, geometry: Geometry, a: Callable, f: Callable, **kargs) -> None:
        """Creates a nonlineal 1D equation with the form:

        .. math::
            -\\frac{d}{dx}\\left(a(x)u\\frac{du}{dx}\\right)=f(x)

        Args:
            geometry (Geometry): Input lineal geometry
            a (Callable): Function a
            f (Callable): Function f
        """

        self.a = a
        self.f = f
        Core.__init__(self, geometry, solver=NoLineal.Newton, **kargs)
        self.name = '1D non lineal sample equation'
        self.properties['a'] = None
        self.properties['f'] = None

    def elementMatrices(self) -> None:
        """Calculate the element matrices usign Reddy's non lineal finite element model. Element matrices and forces are calculated with Gauss-Legendre quadrature. Point number depends of element discretization.
        """

        for e in tqdm(self.elements, unit='Element'):
            e.Te = np.zeros(e.Ke.shape)
            e.Fe = np.zeros(e.Fe.shape)
            e.Ke = np.zeros(e.Ke.shape)
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
                        e.Ke[i, j] += \
                            (self.a(_x[k])*e.Ue[0]@_p[k]*dpx[k][0]
                                [i]*dpx[k][0][j])*detjac[k]*e.W[k]
                        e.Te[i, j] += \
                            (_p[k][j]*dpx[k][0]
                                [i]*e.Ue[0]@dpx[k][0])*detjac[k]*e.W[k]
                for k in range(len(e.Z)):  # Iterate over gauss points on domain
                    e.Fe[i][0] += self.f(_x[k])*_p[k][i]*detjac[k]*e.W[k]
            e.Te += e.Ke
            # e.Fe[:,0] = 2*self.G*self._phi*detjac@_p
            # e.Ke = (np.transpose(dpx,axes=[0,2,1]) @ dpx).T @ detjac

    def postProcess(self) -> None:
        """Generate graph of solution and solution derivative
        """

        X = []
        U1 = []
        U2 = []
        fig = plt.figure()
        ax1 = fig.add_subplot(1, 2, 1)
        ax2 = fig.add_subplot(1, 2, 2)
        for e in tqdm(self.elements, unit='Element'):
            _x, _u, du = e.giveSolution(True)
            X += _x.T[0].tolist()
            U1 += _u[0].tolist()
            U2 += (du[:, 0, 0]).tolist()
        ax1.plot(X, U1)
        ax2.plot(X, np.array(U2))
        ax1.grid()
        ax2.grid()
        ax1.set_title(r'$U(x)$')
        ax2.set_title(r'$\frac{dU}{dx}$')
