"""Lineal element definition
"""

from .Element1D import Element1D, np


class LinealElement(Element1D):
    """Create a lineal element

    Args:
        coords (np.ndarray): Element coordinates matrix
        gdl (np.ndarray): Element degree of freedom matrix
        n (int, optional): Number of gauss points. Defaults to 3.
    """

    def __init__(self, coords: np.ndarray, gdl: np.ndarray, n: int = 2, **kargs):
        """Create a lineal element

        Args:
            coords (np.ndarray): Element coordinates matrix
            gdl (np.ndarray): Element degree of freedom matrix
            n (int, optional): Number of gauss points. Defaults to 2.
        """
        Element1D.__init__(self, coords, gdl, n, **kargs)

    def psis(self, z: np.ndarray) -> np.ndarray:
        """Calculates the shape functions of the lineal element of a given natural coordinates

        Args:
            z (np.ndarray): Natural coordinates matrix

        Returns:
            np.ndarray: Shape function evaluated in Z points
        """

        return np.array([0.5*(1.0-z), 0.5*(1.0+z)]).T

    def dpsis(self, z: np.ndarray) -> np.ndarray:
        """Calculates the shape functions derivatives of the lineal element of a given natural coordinates

        Args:
            z (np.ndarray): Natural coordinates matrix

        Returns:
            np.ndarray: Shape function derivatives evaluated in Z points
        """
        kernel = 1+(z-z)
        return np.array([[-0.5*kernel], [0.5*kernel]])
