"""Defines a Lagrangian 1 order rectangular element 
"""


from .Element2D import Element2D, np
from .RectangularScheme import RectangularScheme
from ..E1D.LinealElement import LinealElement


class Quadrilateral(Element2D, RectangularScheme):
    """Creates a lagrangian rectangular element of order 1

    Args:
        coords (np.ndarray): Element coordinates matrix
        gdl (np.ndarray): Element gdl matrix
        n (int, optional): Number of Gauss Points. Defaults to 2.
    """

    def __init__(self, coords: np.ndarray, gdl: np.ndarray, n: int = 2) -> None:
        """Creates a lagrangian rectangular element of order 1

        Args:
            coords (np.ndarray): Element coordinates matrix
            gdl (np.ndarray): Element gdl matrix
            n (int, optional): Number of Gauss Points. Defaults to 2.
        """

        coords = np.array(coords)

        he1 = np.linalg.norm(coords[1]-coords[0])
        e1 = LinealElement(np.array([[0], [he1]]),
                           np.array([[-1, -1]]), border=True)

        he2 = np.linalg.norm(coords[2]-coords[1])
        e2 = LinealElement(np.array([[0], [he2]]),
                           np.array([[-1, -1]]), border=True)

        he3 = np.linalg.norm(coords[3]-coords[2])
        e3 = LinealElement(np.array([[0], [he3]]),
                           np.array([[-1, -1]]), border=True)

        he4 = np.linalg.norm(coords[0]-coords[3])
        e4 = LinealElement(np.array([[0], [he4]]),
                           np.array([[-1, -1]]), border=True)

        self.borders = [e1, e2, e3, e4]

        Element2D.__init__(self, coords, coords, gdl)
        RectangularScheme.__init__(self, n)

    def psis(self, z: np.ndarray) -> np.ndarray:
        """Calculates the shape functions of a given natural coordinates

        Args:
            z (np.ndarray): Natural coordinates matrix

        Returns:
            np.ndarray: Shape function evaluated in Z points
        """

        return np.array(
            [0.25*(1.0-z[0])*(1.0-z[1]),
             0.25*(1.0+z[0])*(1.0-z[1]),
             0.25*(1.0+z[0])*(1.0+z[1]),
             0.25*(1.0-z[0])*(1.0+z[1])]).T

    def dpsis(self, z: np.ndarray) -> np.ndarray:
        """Calculates the shape functions derivatives of a given natural coordinates

        Args:
            z (np.ndarray): Natural coordinates matrix

        Returns:
            np.ndarray: Shape function derivatives evaluated in Z points
        """

        return np.array(
            [[0.25*(z[1]-1.0), 0.25*(z[0]-1.0)],
             [-0.25*(z[1]-1.0), -0.25*(z[0]+1.0)],
             [0.25*(z[1]+1.0), 0.25*(1.0+z[0])],
             [-0.25*(1.0+z[1]), 0.25*(1.0-z[0])]])
