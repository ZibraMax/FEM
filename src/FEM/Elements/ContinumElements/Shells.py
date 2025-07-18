from .ContinumBase import ContinumBase
from ..E2D import LTriangular, Quadrilateral
import numpy as np


class ShellBase(ContinumBase):

    def __init__(self, **kargs):
        ContinumBase.__init__(self, **kargs)

    def calculate_BNL(self, dpx) -> np.ndarray:

        return

    def calculate_BL(self, dpx) -> np.ndarray:
        return

    def organize_S(self, S) -> tuple:
        return

    def transformation_matrix(self, deformed=True) -> np.ndarray:
        num_nodes = self.coords.shape[0]
        T_plane = self.rotation_matrix(deformed)
        T_e = np.zeros((3 * num_nodes, 2 * num_nodes))
        for i in range(num_nodes):
            row_start = 3 * i
            col_start = 2 * i
            T_e[row_start:row_start+3, col_start:col_start +
                2] = T_plane
        return T_e.T

    def rotation_matrix(self, deformed=True) -> np.ndarray:
        coords = self.coords + self.Ue.T*deformed
        # This should be in the class definition
        _, dpn = self.J(self.center.T)
        jac = dpn[0] @ coords
        e3 = np.cross(*jac, axis=0)
        e1 = jac[0].copy()
        e2 = jac[1].copy()

        # Ortonormalizar la cosa
        e1 /= np.linalg.norm(e1)
        e2 = np.cross(e3, e1)
        e2 /= np.linalg.norm(e2)
        T_plane = np.column_stack((e1, e2))
        return T_plane

    def get_local_jacobian(self, jac: np.ndarray, dni: np.array) -> np.ndarray:
        R = self.rotation_matrix(False)
        coords = self.coords.copy()
        t_coords = coords@R
        new_jac = dni @ t_coords
        return new_jac


class QuadShellLinear(ShellBase, Quadrilateral):
    def __init__(self, coords: np.ndarray, gdl: np.ndarray, **kargs):
        Quadrilateral.__init__(self, coords, gdl, 3, **kargs)
        ShellBase.__init__(self, **kargs)


class TriShelleLinear(ShellBase, LTriangular):
    def __init__(self, coords: np.ndarray, gdl: np.ndarray, **kargs):
        LTriangular.__init__(self, coords, gdl, 2, **kargs)
        ShellBase.__init__(self, **kargs)
