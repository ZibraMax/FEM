from .ContinumBase import ContinumBase
from ..E1D import LinealElement, QuadraticElement
import numpy as np


class BarBase(ContinumBase):

    def __init__(self, **kargs):
        ContinumBase.__init__(self, **kargs)

    def calculate_BNL(self, dpx) -> np.ndarray:
        """Calculate the BNL matrix for the element.

        This method should be implemented in derived classes.
        """
        return np.array([[dpx[0, 0], dpx[0, 1]]])

    def calculate_BL(self, dpx) -> np.ndarray:
        """Calculate the BL matrix for the element.

        This method should be implemented in derived classes.
        """
        dp1dx = dpx[0, 0]
        dp2dx = dpx[0, 1]
        dudx = np.sum(dpx * self.t_Ue.T)
        return np.array([[dp1dx+dudx*dp1dx, dp2dx+dudx*dp2dx]])

    def organize_S(self, S) -> tuple:
        """Organize the stress tensor S and the constitutive matrix C.

        This method should be implemented in derived classes.
        """
        return S.reshape([1, 1]), S.reshape([1, 1])

    def transformation_matrix(self, deformed=True) -> np.ndarray:
        """Calculate the transformation matrix for the element.

        This method should be implemented in derived classes.
        """
        R = self.rotation_matrix(deformed)
        o = np.zeros_like(R)
        # No se, no creo que sea tan facil tampoco...
        T = np.array([[*R, *o],
                      [*o, *R]])
        return T

    def rotation_matrix(self, deformed=True) -> np.ndarray:
        """Calculate the rotation matrix for the element.

        This method should be implemented in derived classes.
        """
        coords = self.coords + self.Ue.T*deformed
        RR = coords[-1] - coords[0]
        norm = np.linalg.norm(RR)
        R = RR / norm
        return R

    def get_local_jacobian(self, jac: np.ndarray, dni) -> np.ndarray:
        """Get the local Jacobian matrix for the element.

        Args:
            jac (np.ndarray): The Jacobian matrix from the integration points.

        Returns:
            np.ndarray: The local Jacobian matrix.
        """
        coords = self.coords
        RR = coords[-1] - coords[0]
        norm = np.linalg.norm(RR)
        return np.array([[norm/2]])


class BarLinear(BarBase, LinealElement):
    def __init__(self, coords: np.ndarray, gdl: np.ndarray, **kargs):
        LinealElement.__init__(self, coords, gdl, 1, **kargs)
        BarBase.__init__(self, **kargs)


class BarQuadratic(BarBase, QuadraticElement):
    def __init__(self, coords: np.ndarray, gdl: np.ndarray, **kargs):
        QuadraticElement.__init__(self, coords, gdl, 2, **kargs)
        BarBase.__init__(self, **kargs)
