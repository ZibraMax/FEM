"""Defines a Lagrangian 1 order rectangular element 
"""


from .Element3D import Element3D, np
from .BrickScheme import BrickScheme


class Quadrilateral(Element3D, BrickScheme):

    def __init__(self, coords: np.ndarray, gdl: np.ndarray, n: int = 2) -> None:

        coords = np.array(coords)

        Element3D.__init__(self, coords, coords, gdl)
        BrickScheme.__init__(self, n)

    def psis(self, z: np.ndarray) -> np.ndarray:
        """Calculates the shape functions of a given natural coordinates

        Args:
            z (np.ndarray): Natural coordinates matrix

        Returns:
            np.ndarray: Shape function evaluated in Z points
        """
        z = z[0]
        n = z[1]
        g = z[2]
        return 1/8*np.array(
            [(1-z)*(1-n)*(1-g),
             (1+z)*(1-n)*(1-g),
             (1+z)*(1+n)*(1-g),
             (1-z)*(1+n)*(1-g),
             (1-z)*(1-n)*(1+g),
             (1+z)*(1-n)*(1+g),
             (1+z)*(1+n)*(1+g),
             (1-z)*(1+n)*(1+g)]).T

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
