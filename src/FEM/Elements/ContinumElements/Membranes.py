from .ContinumBase import ContinumBase
from ..E2D import LTriangular, Quadrilateral
import numpy as np


class MembraneBase(ContinumBase):

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

    def get_local_jacobian(self, jac: np.ndarray) -> np.ndarray:
        return


class QuadMembraneLinear(MembraneBase, Quadrilateral):
    def __init__(self, coords: np.ndarray, gdl: np.ndarray, **kargs):
        Quadrilateral.__init__(self, coords, gdl, 3, **kargs)
        MembraneBase.__init__(self, **kargs)


class TriMembraneLinear(MembraneBase, LTriangular):
    def __init__(self, coords: np.ndarray, gdl: np.ndarray, **kargs):
        LTriangular.__init__(self, coords, gdl, 2, **kargs)
        MembraneBase.__init__(self, **kargs)
