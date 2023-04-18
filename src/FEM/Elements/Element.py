"""Element class
"""


import logging
import numpy as np
from typing import Callable

# TODO make list of avaliable Element atributes.


class Element():
    """Generates a generic element.

    Args:
        coords (np.ndarray): Vertical coordinates matrix
        _coords (np.ndarray): Vertical coordinates matrix for graphical interfaces
        gdl (np.ndarray): Degree of freedom matrix. Each row is a variable.
        border (bool): True if the element is part of the border domain of another element.
    """

    def __init__(self, coords: np.ndarray, _coords: np.ndarray, gdl: np.ndarray, border: bool = False, fast: bool = False) -> None:
        """Generates a generic element.

        Args:
            coords (np.ndarray): Vertical coordinates matrix
            _coords (np.ndarray): Vertical coordinates matrix for graphical interfaces
            gdl (np.ndarray): Degree of freedom matrix. Each row is a variable.
            border (bool): True if the element is part of the border domain of another element.
            fast (bool): If True, the element does not record Ke, Me, Fe, Qe. Defaults to False.
        """

        self.coords = coords
        self._coords = _coords
        self.border = border
        self.gdl = gdl
        self.fast = fast
        self.gdlm = []
        for i in range(len(self.gdl)):
            for j in range(len(self.gdl[i])):
                self.gdlm.append(self.gdl[i, j])
        self.n = int(len(self.gdl)*len(self.gdl[0]))
        # TODO this was only intended for 2D plane stress/strain elements
        self.properties = {'load_x': [], 'load_y': []}
        self.intBorders = False
        self._x, self._p = self.T(self.Z.T)
        self.jacs, self.dpz = self.J(self.Z.T)
        self._xcenter = self.T(self.center.T)[0].flatten()
        if not self.border:
            if not self.fast:
                self.Ke = np.zeros([self.n, self.n])
                self.Fe = np.zeros([self.n, 1])
                self.Qe = np.zeros([self.n, 1])

            # Specific transformations
            self.detjac = np.linalg.det(self.jacs)
            _j = np.linalg.inv(self.jacs)
            self.dpx = _j @ self.dpz
        self.Ue = np.zeros(self.gdl.shape)

    def restartMatrix(self) -> None:
        """Sets all element matrices and vectors to 0 state
        """
        if not self.border:
            self.Ke[:, :] = 0.0
            self.Fe[:, :] = 0.0
            self.Ue[:, :] = 0.0
            self.Qe[:, :] = 0.0

    def T(self, z: np.ndarray) -> np.ndarray:
        """Give the global coordinates of given natural coordiantes over element

        Args:
            z (np.ndarray): Natural coordinates matrix. Each row is a dimension, each column is a point.

        Returns:
            np.ndarray: Global coordinates matrix and shape functions
        """
        p = self.psis(z)
        return p@self.coords, p

    def TS(self, z):
        """Returns the transformation of a given set of points in the element. This method is used for border elements

        Args:
            z (np.ndarray): Natural coordinates matrix. Each row is a dimension, each column is a point.

        Returns:
            np.ndarray: Global coordinates matrix
        """
        return self.s0+self.dir*self.T(z)[0]

    def inverseMapping(self, x0: np.ndarray, n: int = 100) -> np.ndarray:
        """Give the natural coordinates of given global coordinates over elements using Newton's method

        Args:
            x0(np.ndarray): Global coordinates matrix
            n(int, optional): Máximun number of iterations. Defaults to 100.

        Returns:
            np.ndarray: Natural coordinates matrix
        """

        tol = 1*10**(-6)
        zi = np.zeros(x0.shape)+0.25
        for _ in range(n):
            xi = x0 - self.T(zi)[0].T
            _J = np.linalg.inv(self.J(zi)[0])
            xi = xi.T
            xi = xi.reshape(list(xi.shape)+[1])
            dz = _J@xi
            zi += dz[:, :, 0].T
            if np.max(np.abs(dz)) < tol:
                break
        return zi

    def J(self, z: np.ndarray) -> np.ndarray:
        """Calculate the jacobian matrix over a set of natural coordinates

        Args:
            z(np.ndarray): Natural coordinates matrix. Each row is a dimension, each column is a point.

        Returns:
            np.ndarray: Jacobian's matrices and shape function derivatives
        """
        dpsis = self.dpsis(z).T
        return dpsis @ self.coords, dpsis

    def giveSolution(self, SVSolution: bool = False, domain: str = 'domain') -> np.ndarray:
        """Calculate the interpolated solution over element domain

        Args:
            SVSolution(bool, optional): To calculate second variable solutions. Defaults to False.
            domain(str, optional): Where to give the solution ['domain' or 'gauss-points']. Defaults to 'domain'.

        Returns:
            np.ndarray: Arrays of coordinates, solutions and second variables solutions.
        """

        _z = self.domain
        if domain == 'gauss-points':
            _z = self.Z
        _x, _p = self.T(_z.T)
        if SVSolution:
            j, dpz = self.J(_z.T)
            dpx = np.linalg.inv(j) @ dpz
            # print((self.Ue @ np.transpose(dpx,axes=[0,2,1])).shape)
            return _x, self.Ue@_p.T, self.Ue @ np.transpose(dpx, axes=[0, 2, 1])
        return _x, self.Ue@_p.T

    def giveSolutionPoint(self, Z: np.ndarray, SVSolution: bool = False) -> np.ndarray:
        """Calculate the interpolated solution over given set of points

        Args:
            Z(np.ndarray): Natural coordintas to extract the solution
            SVSolution(bool, optional): To calculate second variable solution. Defaults to False.

        Returns:
            np.ndarray: Arrays of coordinates, solutions and second variables solutions.
        """

        _x, _p = self.T(Z)
        if SVSolution:
            j, dpz = self.J(Z)
            dpx = np.linalg.inv(j) @ dpz
            return _x, self.Ue@_p.T, self.Ue @ np.transpose(dpx, axes=[0, 2, 1])
        return _x, self.Ue@_p.T

    def setUe(self, U: np.ndarray) -> None:
        """Assing element local solution

        Args:
            U(np.ndarray): Global solution
        """

        for i in range(len(self.gdl)):
            self.Ue[i] = U[np.ix_(self.gdl[i])].flatten()
        n = len(self._coords)
        m = len(self.gdl)
        self._Ueg = self.Ue[np.ix_(np.linspace(
            0, m-1, m).astype(int), np.linspace(0, n-1, n).astype(int))]
        self._Ueg = np.array(self._Ueg.T.tolist()+[self._Ueg.T[0].tolist()])

    def integrate(self, f: Callable) -> float:
        """Calculate the integral of f function over element domain

        Args:
            f(function): Function to be integrated

        Returns:
            float: Integral value
        """

        integral = 0
        for w, z in zip(self.W, self.Z):
            integral += f(z)*w
        return integral
