"""Element class
"""


import numpy as np
from typing import Tuple, Callable


class Element():
    """Generates a generic element.

    Args:
        coords (np.ndarray): Vertical coordinates matrix
        _coords (np.ndarray): Vertical coordinates matrix for graphical interfaces
        gdl (np.ndarray): Degree of freedom matrix. Each row is a variable.
    """

    def __init__(self, coords: np.ndarray, _coords: np.ndarray, gdl: np.ndarray) -> None:
        """Generates a generic element.

        Args:
            coords (np.ndarray): Vertical coordinates matrix
            _coords (np.ndarray): Vertical coordinates matrix for graphical interfaces
            gdl (np.ndarray): Degree of freedom matrix. Each row is a variable.
        """

        self.coords = coords
        self._coords = _coords
        self.gdl = gdl
        self.gdlm = []
        for i in range(len(self.gdl)):
            for j in range(len(self.gdl[i])):
                self.gdlm.append(self.gdl[i, j])
        self.n = int(len(self.gdl)*len(self.gdl[0]))
        self.Ke = np.zeros([self.n, self.n])
        self.Fe = np.zeros([self.n, 1])
        self.Ue = np.zeros(self.gdl.shape)
        self.Qe = np.zeros([self.n, 1])

    def T(self, z: np.ndarray) -> np.ndarray:
        """Give the global coordinates of given natural coordiantes over element

        Args:
            z (np.ndarray): Natural coordinates matrix

        Returns:
            np.ndarray: Global coordinates matrix
        """

        p = self.psis(z)
        return p@self.coords, p

    def inverseMapping(self, x0: np.ndarray, n: int = 100) -> np.ndarray:
        """Give the natural coordinates of given global coordinates over elements using Newton's method

        Args:
            x0 (np.ndarray): Global coordinates matrix
            n (int, optional): Máximun number of iterations. Defaults to 100.

        Returns:
            np.ndarray: Natural coordinates matrix
        """

        tol = 1*10**(-6)
        zi = np.zeros(x0.shape)+0.1
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
            z (np.ndarray): Natural coordinates matrix

        Returns:
            np.ndarray: Jacobian's matrices
        """

        dpsis = self.dpsis(z).T
        return dpsis @ self.coords, dpsis

    def giveSolution(self, SVSolution: bool = False) -> np.ndarray:
        """Calculate the interpolated solution over element domain

        Args:
            SVSolution (bool, optional): To calculate second variable solutions. Defaults to False.

        Returns:
            np.ndarray: Arrays of coordinates, solutions and second variables solutions.
        """

        _z = self.domain
        _x, _p = self.T(_z.T)
        if SVSolution:
            j, dpz = self.J(_z.T)  # TODO Revisar con Reddy
            dpx = np.linalg.inv(j) @ dpz
            # print((self.Ue @ np.transpose(dpx,axes=[0,2,1])).shape)
            # TODO REVISAR VS
            return _x, self.Ue@_p.T, self.Ue @ np.transpose(dpx, axes=[0, 2, 1])
        return _x, self.Ue@_p.T

    def giveSolutionPoint(self, Z: np.ndarray, SVSolution: bool = False) -> np.ndarray:
        """Calculate the interpolated solution over given set of points

        Args:
            Z (np.ndarray): Natural coordintas to extract the solution
            SVSolution (bool, optional): To calculate second variable solution. Defaults to False.

        Returns:
            np.ndarray: Arrays of coordinates, solutions and second variables solutions.
        """
        _x, _p = self.T(Z)
        if SVSolution:
            j, dpz = self.J(Z)  # TODO Revisar con Reddy
            dpx = np.linalg.inv(j) @ dpz
            # TODO REVISAR VS
            return _x, self.Ue@_p.T, self.Ue @ np.transpose(dpx, axes=[0, 2, 1])
        return _x, self.Ue@_p.T

    def setUe(self, U: np.ndarray) -> None:
        """Assing element local solution

        Args:
            U (np.ndarray): Global solution
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
            f (function): Function to be integrated

        Returns:
            float: Integral value
        """

        integral = 0
        for w, z in zip(self.W, self.Z):
            integral += f(z)*w
        return integral
