"""Defines a general 2D element"""


from ..Element import Element
import numpy as np


class Element3D(Element):
    """Create a 3D element

    Args:
        coords (np.ndarray): Element coordinate matrix
        _coords (np.ndarray): Element coordinate matrix for graphical interface purposes
        gdl (np.ndarray): Degree of freedom matrix
    """

    def __init__(self, coords: np.ndarray, _coords: np.ndarray, gdl: np.ndarray, **kargs) -> None:
        """Create a 2D element

        Args:
            coords (np.ndarray): Element coordinate matrix
            _coords (np.ndarray): Element coordinate matrix for graphical interface purposes
            gdl (np.ndarray): Degree of freedom matrix
        """

        Element.__init__(self, coords, _coords, gdl, **kargs)
        self._coordsg = np.array(
            self._coords.tolist()+[self._coords[0].tolist()])

    def draw(self) -> None:
        """Create a graph of element
        """
        pass

    def jacobianGraph(self) -> None:
        """Create the determinant jacobian graph
        """
        pass

    def isInside(self, x: np.ndarray) -> np.ndarray:
        """Test if a given points is inside element domain

        Args:
            x (np.ndarray): Point to be tested

        Returns:
            np.ndarray: Boolean array of test result
        """
        rest = [True]*len(x)
        i = -1
        for p in x:
            i += 1
            for f in self.faces:
                coords = self.coords[np.ix_(f)]
                gdl = np.array([[-1]*len(coords)])
                e = self.face_element(coords, gdl, fast=True, border=True)
                _x, _p = e.T(e.center.T)
                _j, _dp = e.J(e.center.T)
                normal = np.cross(_j[:, 0, :].T, _j[:, 1, :].T, axis=0)
                vector = p - _x[0]
                res = np.dot(vector, normal.flatten())
                if res > 0:
                    rest[i] = False
                    break
        return rest
