"""Define a o2 serendipity element
"""


from .Element2D import *
from .RectangularScheme import *


class Serendipity(Element2D, RectangularScheme):
    """Creates a Serendipity element

    Args:
        coords (np.ndarray): Coordinate matrix of element
        gdl (np.ndarray): Coordinate matrix of element GDL's
        n (int, optional): Number of gauss points. Defaults to 3.
    """

    def __init__(self, coords: np.ndarray, gdl: np.ndarray, n: int = 3) -> None:
        """Creates a Serendipity element

        Args:
            coords (np.ndarray): Coordinate matrix of element
            gdl (np.ndarray): Coordinate matrix of element GDL's
            n (int, optional): Number of gauss points. Defaults to 3.
        """

        _coords = np.array([coords[i] for i in range(4)])
        Element2D.__init__(self, np.array(coords), _coords, gdl)
        RectangularScheme.__init__(self, n)

    def psis(self, z: np.ndarray) -> np.ndarray:
        """Calculates the shape functions of a given natural coordinates

        Args:
            z (np.ndarray): Natural coordinates matrix

        Returns:
            np.ndarray: Shape function evaluated in Z points
        """

        return np.array([
            0.25*(1.0-z[0])*(1.0-z[1])*(-1.0-z[0]-z[1]),
            0.25*(1.0+z[0])*(1.0-z[1])*(-1.0+z[0]-z[1]),
            0.25*(1.0+z[0])*(1.0+z[1])*(-1.0+z[0]+z[1]),
            0.25*(1.0-z[0])*(1.0+z[1])*(-1.0-z[0]+z[1]),
            0.5*(1.0-z[0]**2.0)*(1.0-z[1]),
            0.5*(1.0+z[0])*(1.0-z[1]**2.0),
            0.5*(1.0-z[0]**2.0)*(1.0+z[1]),
            0.5*(1.0-z[0])*(1.0-z[1]**2.0)
        ]).T

    def dpsis(self, z: np.ndarray) -> np.ndarray:
        """Calculates the shape functions derivatives of a given natural coordinates

        Args:
            z (np.ndarray): Natural coordinates matrix

        Returns:
            np.ndarray: Shape function derivatives evaluated in Z points
        """

        return np.array(
            [[-0.25*(z[1]-1.0)*(2.0*z[0]+z[1]), -0.25*(z[0]-1.0)*(2.0*z[1]+z[0])],
             [-0.25*(z[1]-1.0)*(2.0*z[0]-z[1]),
              0.25*(z[0]+1.0)*(2.0*z[1]-z[0])],
             [0.25*(z[1]+1.0)*(2.0*z[0]+z[1]),
              0.25*(z[0]+1.0)*(2.0*z[1]+z[0])],
             [0.25*(z[1]+1.0)*(2.0*z[0]-z[1]), -
              0.25*(z[0]-1.0)*(2.0*z[1]-z[0])],
             [(z[1]-1.0)*z[0], 0.5*(z[0]**2.0-1.0)],
             [-0.5*(z[1]**2.0-1.0), -z[1]*(z[0]+1.0)],
             [-(z[1]+1.0)*z[0], -0.5*(z[0]**2.0-1.0)],
             [0.5*(z[1]**2.0-1.0), z[1]*(z[0]-1.0)]])
