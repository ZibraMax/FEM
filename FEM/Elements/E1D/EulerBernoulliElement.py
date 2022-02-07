"""Creates a 1D beam element
"""
from .LinealElement import LinealElement
import numpy as np


class EulerBernoulliElement(LinealElement):
    """Creates a 1D beam element
    """

    def __init__(self, coords, gdl, n=2, nvn=2) -> None:
        """Creates a 1D beam element

        Args:
            coords (list): Beam coordiantes
            gdl (list): Degrees of freedom
            n (int, optional): Number of gauss points. Defaults to 2.
            nvn (int, optional): Number of vairables per node. Defaults to 2.
        """
        gdlcopy = gdl.copy()
        gdl[0, 1] = gdlcopy[1, 0]
        gdl[1, 0] = gdlcopy[0, 1]
        if nvn > 2:
            gdl = gdlcopy.T.reshape([3, 2])

        LinealElement.__init__(self, coords, gdl, n=n)
        self.he = np.linalg.norm(self.coords[-1]-self.coords[0])
        self.Zr, self.Wr = np.polynomial.legendre.leggauss(n-1)
        self.n = 2*nvn

    def hermit(self, z: np.ndarray) -> np.ndarray:
        """
        Calculates the shape functions of the lineal element of a given natural coordinates

        Args:
            z (np.ndarray): Natural coordinates matrix

        Returns:
            np.ndarray: Shape function evaluated in Z points
        """
        he = self.he
        x = (he)/2*z+(he)/2
        return np.array([1-3*(x/he)**2+2*(x/he)**3, -x*(1-x/he)**2, 3*(x/he)**2-2*(x/he)**3, -x*((x/he)**2-x/he)]).T

    def dhermit(self, z: np.ndarray) -> np.ndarray:
        """Calculates the shape functions derivatives of the lineal element of a given natural coordinates

        Args:
            z (np.ndarray): Natural coordinates matrix

        Returns:
            np.ndarray: Shape function derivatives evaluated in Z points
        """

        h = self.he
        x = (h)/2*z+(h)/2
        return np.array([
            [-6/h*x/h*(1-x/h), -(1+3*(x/h)**2-4*x/h),
             6/h*x/h*(1-x/h), -x/h*(3*x/h-2)],
            [-6/h**2*(1-2*x/h), -2/h*(3*x/h-2), 6/h **
             2*(1-2*x/h), -2/h*(3*x/h-1)],
            [12/h**3+(x-x), -6/h**2+(x-x), -12/h**3+(x-x), -6/h**2+(x-x)]])

    def giveSolution(self, SVSolution: bool = False, domain: str = 'domain') -> np.ndarray:
        """Calculate the interpolated solution over element domain

        Args:
            SVSolution(bool, optional): To calculate second variable solutions. Defaults to False.

        Returns:
            np.ndarray: Arrays of coordinates, solutions and second variables solutions.
        """

        # TODO hacer una comprobación de frontera para evitar errores
        _z = self.domain
        if domain == 'gauss-points':
            _z = self.Z
        _x, _ = self.T(_z.T)
        _h = self.hermit(_z.T)
        if SVSolution:
            _dh = self.dhermit(_z.T)
            return _x, self.Ue.flatten()@_h.T, self.Ue.flatten() @ _dh.T
        return _x, self.Ue.flatten()@_h.T

    def giveSolutionPoint(self, Z: np.ndarray, SVSolution: bool = False) -> np.ndarray:
        """Calculate the interpolated solution over given set of points

        Args:
            Z(np.ndarray): Natural coordintas to extract the solution
            SVSolution(bool, optional): To calculate second variable solution. Defaults to False.

        Returns:
            np.ndarray: Arrays of coordinates, solutions and second variables solutions.
        """

        # TODO hacer una comprobación de frontera para evitar errores
        _x, _p = self.T(Z)
        if SVSolution:
            j, dpz = self.J(Z)  # TODO Revisar con Reddy
            dpx = np.linalg.inv(j) @ dpz
            # TODO REVISAR VS
            return _x, self.Ue@_p.T, self.Ue @ np.transpose(dpx, axes=[0, 2, 1])
        return _x, self.Ue@_p.T
