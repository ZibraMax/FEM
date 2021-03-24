"""Define a o2 serendipity element
"""


from .Element2D import Element2D, np
from .RectangularScheme import RectangularScheme
from ..E1D.QuadraticElement import QuadraticElement


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
        coords = np.array(coords)

        he1 = np.linalg.norm(coords[1]-coords[0])
        e1 = QuadraticElement(np.array([[0], [he1*0.5], [he1]]),
                              np.array([[-1, -1, -1]]), border=True)

        he2 = np.linalg.norm(coords[2]-coords[1])
        e2 = QuadraticElement(np.array([[0], [he2*0.5], [he2]]),
                              np.array([[-1, -1, -1]]), border=True)

        he3 = np.linalg.norm(coords[3]-coords[2])
        e3 = QuadraticElement(np.array([[0], [he3*0.5], [he3]]),
                              np.array([[-1, -1, -1]]), border=True)

        he4 = np.linalg.norm(coords[0]-coords[3])
        e4 = QuadraticElement(np.array([[0], [he4*0.5], [he4]]),
                              np.array([[-1, -1, -1]]), border=True)

        self.borders = [e1, e2, e3, e4]
        Element2D.__init__(self, coords, _coords, gdl)
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
