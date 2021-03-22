"""Defines a Lagrangian 1 order triangular element 
"""


from .Element2D import *
from .TriangularScheme import *


class LTriangular(Element2D, TriangularScheme):
    """Creates a lagrangian triangular element of order 1

    Args:
        coords (np.ndarray): Element coordinates matrix
        gdl (np.ndarray): Element gdl matrix
        n (int, optional): Number of Gauss Points. Defaults to 2.
    """

    def __init__(self, coords: np.ndarray, gdl: np.ndarray, n: int = 2) -> None:
        """Creates a lagrangian triangular element of order 1

        Args:
            coords (np.ndarray): Element coordinates matrix
            gdl (np.ndarray): Element gdl matrix
            n (int, optional): Number of Gauss Points. Defaults to 2.
        """

        coords = np.array(coords)
        Element2D.__init__(self, coords, coords, gdl)
        TriangularScheme.__init__(self, n)

    def psis(self, z: np.ndarray) -> np.ndarray:
        """Calculates the shape functions of a given natural coordinates

        Args:
            z (np.ndarray): Natural coordinates matrix

        Returns:
            np.ndarray: Shape function evaluated in Z points
        """

        return np.array([
            1.0-z[0]-z[1],
            z[0],
            z[1]]).T

    def dpsis(self, z: np.ndarray) -> np.ndarray:
        """Calculates the shape functions derivatives of a given natural coordinates

        Args:
            z (np.ndarray): Natural coordinates matrix

        Returns:
            np.ndarray: Shape function derivatives evaluated in Z points
        """

        kernell = (z[0]-z[0])
        return np.array([
            [-1.0*(1+kernell), -1.0*(1+kernell)],
            [1.0*(1+kernell), 0.0*(1+kernell)],
            [0.0*(1+kernell), 1.0*(1+kernell)]
        ])
