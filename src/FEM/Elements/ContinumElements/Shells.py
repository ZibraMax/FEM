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
        return

    def rotation_matrix(self, deformed=True) -> np.ndarray:
        return

    def get_local_jacobian(self, jac: np.ndarray, dni: np.array) -> np.ndarray:
        return


class QuadShellLinear(ShellBase, Quadrilateral):
    def __init__(self, coords: np.ndarray, gdl: np.ndarray, **kargs):
        Quadrilateral.__init__(self, coords, gdl, 3, **kargs)
        ShellBase.__init__(self, **kargs)


class TriShelleLinear(ShellBase, LTriangular):
    def __init__(self, coords: np.ndarray, gdl: np.ndarray, **kargs):
        LTriangular.__init__(self, coords, gdl, 2, **kargs)
        ShellBase.__init__(self, **kargs)
