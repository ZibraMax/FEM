"""Defines a Lagrangian order 1 brick element 
"""


from ..E2D.Quadrilateral import Quadrilateral
from ..E2D.Serendipity import Serendipity
from .Element3D import Element3D, np
from .BrickScheme import BrickScheme


class Brick(Element3D, BrickScheme):
    """Creates a 3D brick element

    Args:
        coords (np.ndarray): Node coordinates matrix
        gdl (np.ndarray): Degrees of freedom matrix
        n (int, optional): Number of gauss points used for integration. Defaults to 3.
    """

    def __init__(self, coords: np.ndarray, gdl: np.ndarray, n: int = 3, **kargs) -> None:
        """Creates a 3D brick element

        Args:
            coords (np.ndarray): Node coordinates matrix
            gdl (np.ndarray): Degrees of freedom matrix
            n (int, optional): Number of gauss points used for integration. Defaults to 3.
        """

        coords = np.array(coords)
        self.faces = [
            [0, 1, 5, 4],
            [1, 2, 6, 5],
            [4, 5, 6, 7],
            [3, 2, 1, 0],
            [2, 3, 7, 6],
            [4, 7, 3, 0]]
        self.face_element = Quadrilateral

        Element3D.__init__(self, coords, coords, gdl, **kargs)
        BrickScheme.__init__(self, n, **kargs)

    def psis(self, _z: np.ndarray) -> np.ndarray:
        """Calculates the shape functions of a given natural coordinates

        Args:
            z (np.ndarray): Natural coordinates matrix

        Returns:
            np.ndarray: Shape function evaluated in Z points
        """

        z = _z[0]
        n = _z[1]
        g = _z[2]
        return 1/8*np.array(
            [(1-z)*(1-n)*(1-g),
             (1+z)*(1-n)*(1-g),
             (1+z)*(1+n)*(1-g),
             (1-z)*(1+n)*(1-g),
             (1-z)*(1-n)*(1+g),
             (1+z)*(1-n)*(1+g),
             (1+z)*(1+n)*(1+g),
             (1-z)*(1+n)*(1+g)]).T

    def dpsis(self, _z: np.ndarray) -> np.ndarray:
        """Calculates the shape functions derivatives of a given natural coordinates

        Args:
            z (np.ndarray): Natural coordinates matrix

        Returns:
            np.ndarray: Shape function derivatives evaluated in Z points
        """
        x = _z[0]
        y = _z[1]
        z = _z[2]
        return 1/8*np.array(
            [[(y-1.0)*(1.0-z), (x-1)*(1-z), -(1-x)*(1-y)],
             [(1-y)*(1-z), (-1-x)*(1-z), -(1+x)*(1-y)],
             [(1+y)*(1-z), (1+x)*(1-z), -(1+x)*(1+y)],
             [(-1.0-y)*(1-z), (1-x)*(1-z), -(1-x)*(1+y)],
             [(1-y)*(-1-z), -(1-x)*(1+z), (1-x)*(1-y)],
             [(1-y)*(1+z), -(1+x)*(1+z), (1+x)*(1-y)],
             [(1+y)*(1+z), (1+x)*(1+z), (1+x)*(1+y)],
             [-(1+y)*(1+z), (1-x)*(1+z), (1-x)*(1+y)]])


class BrickO2(Element3D, BrickScheme):
    """Creates a 3D second order brick element.

    Args:
        coords (np.ndarray): Node coordinates matrix
        gdl (np.ndarray): Degrees of freedom matrix
        n (int, optional): Number of gauss points used for integration. Defaults to 3.
    """

    def __init__(self, coords: np.ndarray, gdl: np.ndarray, n: int = 3, **kargs) -> None:
        """Creates a 3D second order brick element.

        Args:
            coords (np.ndarray): Node coordinates matrix
            gdl (np.ndarray): Degrees of freedom matrix
            n (int, optional): Number of gauss points used for integration. Defaults to 3.
        """

        coords = np.array(coords)
        self.faces = [
            [0, 1, 5, 4, 8, 13, 16, 12],
            [1, 2, 6, 5, 9, 14, 17, 13],
            [4, 5, 6, 7, 16, 17, 18, 19],
            [3, 2, 1, 0, 10, 9, 8, 11],
            [2, 3, 7, 6, 10, 15, 18, 14],
            [4, 7, 3, 0, 19, 15, 11, 12]]
        self.face_element = Serendipity

        Element3D.__init__(self, coords, coords, gdl, **kargs)
        BrickScheme.__init__(self, n, **kargs)

    def psis(self, _z: np.ndarray) -> np.ndarray:
        """Calculates the shape functions of a given natural coordinates

        Args:
            z (np.ndarray): Natural coordinates matrix

        Returns:
            np.ndarray: Shape function evaluated in Z points
        """

        x = _z[0]
        y = _z[1]
        z = _z[2]
        return 1/8*np.array([
            (1-x)*(1-y)*(1-z)*(-x-y-z-2),
            (1+x)*(1-y)*(1-z)*(x-y-z-2),
            (1+x)*(1+y)*(1-z)*(x+y-z-2),
            (1-x)*(1+y)*(1-z)*(-x+y-z-2),
            (1-x)*(1-y)*(1+z)*(-x-y+z-2),
            (1+x)*(1-y)*(1+z)*(x-y+z-2),
            (1+x)*(1+y)*(1+z)*(x+y+z-2),
            (1-x)*(1+y)*(1+z)*(-x+y+z-2),
            2*(1-x**2)*(1-y)*(1-z),
            2*(1+x)*(1-y**2)*(1-z),
            2*(1-x**2)*(1+y)*(1-z),
            2*(1-x)*(1-y**2)*(1-z),
            2*(1-x)*(1-y)*(1-z**2),
            2*(1+x)*(1-y)*(1-z**2),
            2*(1+x)*(1+y)*(1-z**2),
            2*(1-x)*(1+y)*(1-z**2),
            2*(1-x**2)*(1-y)*(1+z),
            2*(1+x)*(1-y**2)*(1+z),
            2*(1-x**2)*(1+y)*(1+z),
            2*(1-x)*(1-y**2)*(1+z)
        ]).T

    def dpsis(self, _z: np.ndarray) -> np.ndarray:
        """Calculates the shape functions derivatives of a given natural coordinates

        Args:
            z (np.ndarray): Natural coordinates matrix

        Returns:
            np.ndarray: Shape function derivatives evaluated in Z points
        """
        x = _z[0]
        y = _z[1]
        z = _z[2]  # I JUST WANT THIS TO WORK PROPERLY
        return 1/8*np.array([[-(1 - x)*(1 - y)*(1 - z) + (1 - z)*(y - 1)*(-x - y - z - 2), -(1 - x)*(1 - y)*(1 - z) + (1 - z)*(x - 1)*(-x - y - z - 2), -(1 - x)*(1 - y)*(1 - z) - (1 - x)*(1 - y)*(-x - y - z - 2)],
                             [(1 - y)*(1 - z)*(x + 1) + (1 - y)*(1 - z)*(x - y - z - 2), -(1 - y)*(1 - z)*(x + 1) + (
                                 1 - z)*(-x - 1)*(x - y - z - 2), -(1 - y)*(1 - z)*(x + 1) - (1 - y)*(x + 1)*(x - y - z - 2)],
                             [(1 - z)*(x + 1)*(y + 1) + (1 - z)*(y + 1)*(x + y - z - 2), (1 - z)*(x + 1)*(y + 1) + (
                                 1 - z)*(x + 1)*(x + y - z - 2), -(1 - z)*(x + 1)*(y + 1) - (x + 1)*(y + 1)*(x + y - z - 2)],
                             [-(1 - x)*(1 - z)*(y + 1) + (1 - z)*(-y - 1)*(-x + y - z - 2), (1 - x)*(1 - z)*(y + 1) + (1 - x)
                              * (1 - z)*(-x + y - z - 2), -(1 - x)*(1 - z)*(y + 1) - (1 - x)*(y + 1)*(-x + y - z - 2)],
                             [-(1 - x)*(1 - y)*(z + 1) + (1 - y)*(-z - 1)*(-x - y + z - 2), -(1 - x)*(1 - y)*(z + 1) -
                              (1 - x)*(z + 1)*(-x - y + z - 2), (1 - x)*(1 - y)*(z + 1) + (1 - x)*(1 - y)*(-x - y + z - 2)],
                             [(1 - y)*(x + 1)*(z + 1) + (1 - y)*(z + 1)*(x - y + z - 2), -(1 - y)*(x + 1)*(z + 1) -
                              (x + 1)*(z + 1)*(x - y + z - 2), (1 - y)*(x + 1)*(z + 1) + (1 - y)*(x + 1)*(x - y + z - 2)],
                             [(x + 1)*(y + 1)*(z + 1) + (y + 1)*(z + 1)*(x + y + z - 2), (x + 1)*(y + 1)*(z + 1) +
                              (x + 1)*(z + 1)*(x + y + z - 2), (x + 1)*(y + 1)*(z + 1) + (x + 1)*(y + 1)*(x + y + z - 2)],
                             [-(1 - x)*(y + 1)*(z + 1) - (y + 1)*(z + 1)*(-x + y + z - 2), (1 - x)*(y + 1)*(z + 1) + (1 - x)
                              * (z + 1)*(-x + y + z - 2), (1 - x)*(y + 1)*(z + 1) + (1 - x)*(y + 1)*(-x + y + z - 2)],
                             [-4*x*(1 - y)*(1 - z), (2 - 2*x**2) *
                              (z - 1), (2 - 2*x**2)*(y - 1)],
                             [2*(1 - y**2)*(1 - z), -2*y*(1 - z) *
                              (2*x + 2), (2*x + 2)*(y**2 - 1)],
                             [-4*x*(1 - z)*(y + 1), (1 - z) *
                              (2 - 2*x**2), (2 - 2*x**2)*(-y - 1)],
                             [-2*(1 - y**2)*(1 - z), -2*y*(1 - z) *
                              (2 - 2*x), (2 - 2*x)*(y**2 - 1)],
                             [-2*(1 - y)*(1 - z**2), (2 - 2*x) *
                              (z**2 - 1), -2*z*(1 - y)*(2 - 2*x)],
                             [2*(1 - y)*(1 - z**2), (2*x + 2) *
                                 (z**2 - 1), -2*z*(1 - y)*(2*x + 2)],
                             [2*(1 - z**2)*(y + 1), (1 - z**2) *
                                 (2*x + 2), -2*z*(2*x + 2)*(y + 1)],
                             [-2*(1 - z**2)*(y + 1), (1 - z**2) *
                              (2 - 2*x), -2*z*(2 - 2*x)*(y + 1)],
                             [-4*x*(1 - y)*(z + 1), (2 - 2*x**2) *
                              (-z - 1), (1 - y)*(2 - 2*x**2)],
                             [2*(1 - y**2)*(z + 1), -2*y*(2*x + 2)
                              * (z + 1), (1 - y**2)*(2*x + 2)],
                             [-4*x*(y + 1)*(z + 1), (2 - 2*x**2) *
                              (z + 1), (2 - 2*x**2)*(y + 1)],
                             [-2*(1 - y**2)*(z + 1), -2*y*(2 - 2*x)*(z + 1), (1 - y**2)*(2 - 2*x)]])