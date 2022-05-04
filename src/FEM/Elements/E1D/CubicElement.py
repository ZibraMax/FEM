"""Lineal element definition
"""


from .Element1D import Element1D, np


class CubicElement(Element1D):
    """Create a cubic element

    Args:
        coords (np.ndarray): Element coordinates matrix
        gdl (np.ndarray): Element degree of freedom matrix
        n (int, optional): Number of gauss points. Defaults to 4.
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

        return np.array([
            -9.0/16.0*(1.0-z)*(1.0/9.0-z*z),
            27.0/16.0*(1.0-z*z)*(1.0/3.0-z),
            27.0/16.0*(1.0-z*z)*(1.0/3.0+z),
            -9.0/16.0*(1.0+z)*(1.0/9.0-z*z)]).T

    def dpsis(self, z: np.ndarray) -> np.ndarray:
        """Calculates the shape functions derivatives of the lineal element of a given natural coordinates

        Args:
            z (np.ndarray): Natural coordinates matrix

        Returns:
            np.ndarray: Shape function derivatives evaluated in Z points
        """
        return np.array([
            [-27.0/16.0*z*z+9.0/8.0*z+1.0/16.0],
            [81.0/16.0*z*z-9.0/8.0*z-27.0/16.0],
            [-81.0/16.0*z*z-9.0/8.0*z+27.0/16.0],
            [27.0/16.0*z*z+9.0/8.0*z-1.0/16.0]])
