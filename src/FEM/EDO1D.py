"""Solve a 1D 1 variable per node Finite Element problem
"""


from .Core import Core, Geometry
from tqdm import tqdm
import numpy as np
import matplotlib.pyplot as plt
from typing import Callable


class EDO1D(Core):
    """Create a 1D 1 variable per node Finite Element problem

    The differential equation is:

    .. math::
        a(x)\\frac{d^2u}{dx^2}+c(x)u=f(x)

    Args:
        geometry (Geometry): 1D Geometry of the problem. Use the Mesh.Lineal class
        a (function): Function a, if a is constant you can use a = lambda x: [value]
        c (function): Function c, if c is constant you can use c = lambda x: [value]
        f (function): Function f, if f is constant you can use f = lambda x: [value]
    """

    def __init__(self, geometry: Geometry, a: Callable, c: Callable, f: Callable) -> None:
        """Create a 1D 1 variable per node Finite Element problem

        The differential equation is:

        .. math::
            a(x)\\frac{d^2u}{dx^2}+c(x)u=f(x)

        Args:
                geometry (Geometry): 1D Geometry of the problem. Use the Mesh.Lineal class
                a (function): Function a, if a is constant you can use a = lambda x: [value]
                c (function): Function c, if c is constant you can use c = lambda x: [value]
                f (function): Function f, if f is constant you can use f = lambda x: [value]
        """

        self.a = a
        self.c = c
        self.f = f
        Core.__init__(self, geometry)
        self.name = 'Generic 1D second order diferential equation'
        self.properties['a'] = None
        self.properties['b'] = None
        self.properties['c'] = None
        self.properties['WARNING'] = "It's not possible lo save callables"

    def elementMatrices(self) -> None:
        """Calculate the element matrices usign Reddy's (2005) finite element model. Element matrices and forces are calculated with Gauss-Legendre quadrature. Point number depends of element discretization.
        """

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
