"""Lineal element definition
"""


from .Element1D import Element1D, np


class QuadraticElement(Element1D):
    """Create a quadratic element

    Args:
        coords (np.ndarray): Element coordinates matrix
        gdl (np.ndarray): Element degree of freedom matrix
        n (int, optional): Number of gauss points. Defaults to 3.
    """

    def __init__(self, coords: np.ndarray, gdl: np.ndarray, n: int = 4, **kargs) -> None:
        """Create a quadratic element

        Args:
            coords (np.ndarray): Element coordinates matrix
            gdl (np.ndarray): Element degree of freedom matrix
            n (int, optional): Number of gauss points. Defaults to 3.
        """

        Element1D.__init__(self, coords, gdl, n, **kargs)

    def psis(self, z: np.ndarray) -> np.ndarray:
        """Calculates the shape functions of the lineal element of a given natural coordinates

        Args:
            z (np.ndarray): Natural coordinates matrix

        Returns:
            np.ndarray: Shape function evaluated in Z points
        """

        zm1 = z+1.0
        return np.array([1.0-3.0/2.0*zm1+zm1*zm1/2.0, 2.0*zm1*(1.0-zm1/2.0), z/2.0*zm1]).T

    def dpsis(self, z: np.ndarray) -> np.ndarray:
        """Calculates the shape functions derivatives of the lineal element of a given natural coordinates

        Args:
            z (np.ndarray): Natural coordinates matrix

        Returns:
            np.ndarray: Shape function derivatives evaluated in Z points
        """
        return np.array([[z-0.5], [-2.0*z], [z+0.5]])
