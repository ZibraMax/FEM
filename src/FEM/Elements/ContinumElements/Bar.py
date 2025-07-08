from .ContinumBase import ContinumBase
import numpy as np


class Bar(ContinumBase):

    def __init__(self, coords: np.ndarray, _coords: np.ndarray, gdl: np.ndarray, **kargs):
        _coords = _coords.copy()
        ContinumBase.__init__(self, coords, _coords, gdl, **kargs)

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
        dudx = self.Ue.T @ dpx
        return np.array([[dp1dx+dudx*dp1dx, dp2dx+dudx*dp2dx]])

    def organize_S(self, S) -> tuple:
        """Organize the stress tensor S and the constitutive matrix C.

        This method should be implemented in derived classes.
        """
        return S.reshape([1, 1]), S.reshape([1, 1])

    def transformation_matrix(self) -> np.ndarray:
        """Calculate the transformation matrix for the element.

        This method should be implemented in derived classes.
        """
        RR = self.coords[-1] - self.coords[0]
        norm = np.linalg.norm(RR)
        R = RR / norm
        o = np.zeros_like(R)
        # No se, no creo que sea tan facil tampoco...
        T = np.array([[*R, *o],
                      [*o, *R]])
        return T

    def get_local_jacobian(self, jac: np.ndarray) -> np.ndarray:
        """Get the local Jacobian matrix for the element.

        Args:
            jac (np.ndarray): The Jacobian matrix from the integration points.

        Returns:
            np.ndarray: The local Jacobian matrix.
        """
        RR = self.coords[-1] - self.coords[0]
        norm = np.linalg.norm(RR)
        return np.array([[norm/2]])
