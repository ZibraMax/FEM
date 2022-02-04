"""Defines a Lagrangian 1 order rectangular element 
"""


from .Element3D import Element3D, np
from .BrickScheme import BrickScheme


class Brick(Element3D, BrickScheme):

    def __init__(self, coords: np.ndarray, gdl: np.ndarray, n: int = 3, **kargs) -> None:

        coords = np.array(coords)

        Element3D.__init__(self, coords, coords, gdl, **kargs)
        BrickScheme.__init__(self, n)

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
