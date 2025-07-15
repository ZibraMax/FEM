from .ContinumBase import ContinumBase
from ..E2D import LTriangular, Quadrilateral
import numpy as np


class MembraneBase(ContinumBase):

    def __init__(self, **kargs):
        ContinumBase.__init__(self, **kargs)

    def _BL1(self, dNdx, dNdy):
        n_nodes = dNdx.shape[0]
        B = np.zeros((3, 2 * n_nodes))
        for i in range(n_nodes):
            B[0, 2 * i] = dNdx[i]
            B[1, 2 * i + 1] = dNdy[i]
            B[2, 2 * i] = dNdy[i]
            B[2, 2 * i + 1] = dNdx[i]
        return B

    def _BL2(self, dNdx, dNdy, Uet):
        n_nodes = dNdx.shape[0]
        BL2 = np.zeros((3, 2 * n_nodes))
        Ux = Uet[0]
        Uy = Uet[1]

        dudx = np.sum(dNdx*Ux)
        dvdx = np.sum(dNdx*Uy)
        dudy = np.sum(dNdy*Ux)
        dvdy = np.sum(dNdy*Uy)

        for i in range(n_nodes):
            BL2[0, 2 * i] = dudx*dNdx[i]
            BL2[1, 2 * i] = dudy*dNdy[i]
            BL2[2, 2 * i] = dudx*dNdy[i] + dudy*dNdx[i]

            BL2[0, 2 * i + 1] = dvdx*dNdx[i]
            BL2[1, 2 * i + 1] = dvdy*dNdy[i]
            BL2[2, 2 * i + 1] = dvdx*dNdy[i] + dvdy*dNdx[i]
        return BL2

    def calculate_BNL(self, dpx) -> np.ndarray:
        dNdx = dpx[0]
        dNdy = dpx[1]
        n_nodes = dNdx.shape[0]
        BNL = np.zeros((4, 2 * n_nodes))
        for i in range(n_nodes):
            BNL[0, 2 * i] = dNdx[i]
            BNL[1, 2 * i] = dNdy[i]

            BNL[2, 2 * i + 1] = dNdx[i]
            BNL[3, 2 * i + 1] = dNdy[i]

        return BNL

    def calculate_BL(self, dpx) -> np.ndarray:
        dNdx = dpx[0]
        dNdy = dpx[1]
        BL1 = self._BL1(dNdx, dNdy)
        BL2 = self._BL2(dNdx, dNdy, self.t_Ue)
        BL = BL1 + BL2
        return BL

    def organize_S(self, S) -> tuple:
        # S is 3x1 stress vector
        # s_stiff second piola kirkoff stress tensor
        # s_force second piola kirkoff stress vector
        S_stiff = np.zeros((4, 4))
        S_stiff[0:2, 0:2] = S
        S_stiff[2:4, 2:4] = S

        S_vect = np.zeros((3, 1))
        S_vect[0] = S[0, 0]
        S_vect[1] = S[1, 1]
        S_vect[2] = S[0, 1]
        return S_stiff, S_vect

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


class QuadMembraneLinear(MembraneBase, Quadrilateral):
    def __init__(self, coords: np.ndarray, gdl: np.ndarray, **kargs):
        Quadrilateral.__init__(self, coords, gdl, 3, **kargs)
        MembraneBase.__init__(self, **kargs)


class TriMembraneLinear(MembraneBase, LTriangular):
    def __init__(self, coords: np.ndarray, gdl: np.ndarray, **kargs):
        LTriangular.__init__(self, coords, gdl, 2, **kargs)
        MembraneBase.__init__(self, **kargs)
