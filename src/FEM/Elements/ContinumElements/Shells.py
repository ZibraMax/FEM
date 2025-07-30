from .ContinumBase import ContinumBase
from ..E2D import LTriangular, Quadrilateral
import numpy as np


class ShellBase(ContinumBase):

    def __init__(self, **kargs):
        ContinumBase.__init__(self, **kargs)

    def set_thickness(self, t: float, e3: np.ndarray = None) -> None:
        """Set the thickness of the shell element.

        Args:
            t (float): Thickness of the shell element.
        """

        # If thickness is not a list, create a list with n elements, n is the number of nodes
        if not isinstance(t, list):
            t = [t] * len(self.coords)
        if len(t) != len(self.coords):
            raise ValueError(
                "Thickness list must match the number of nodes in the element.")
        self.th = t
        # If e3 is not provided, calculate as the normal vector of the shell in the centroid
        if e3 is None:
            _, dpn = self.J(self.center.T)
            jac = dpn[0] @ self.coords
            e3 = np.cross(*jac, axis=0)
            e3 /= np.linalg.norm(e3)
        # If e3 is not a list, create a list with n elements, n is the number of nodes
        self.e3 = []
        if not isinstance(e3, list):
            e3 = np.array(e3)
            for i in range(len(self.coords)):
                self.e3.append(e3.copy())
        else:
            if len(e3) != len(self.coords):
                raise ValueError(
                    "e3 list must match the number of nodes in the element.")
            self.e3 = [np.array(i) for i in e3]  # Garantiza que sea np
        self.e1 = np.zeros_like(self.e3)
        self.e2 = np.zeros_like(self.e3)
        self.e1, self.e2 = self.calculate_e1_e2()

    def calculate_e1_e2(self, deformed=False):
        """Calculate the e1 and e2 vectors for the shell element."""
        E1 = []
        E2 = []
        for i in range(len(self.coords)):
            e3 = self.e3[i] + (self.Ue[3][i]*self.e1[i] -
                               self.Ue[4][i]*self.e2[i])*deformed
            e1 = np.cross([0, 1, 0], e3)
            e1 /= np.linalg.norm(e1)
            e2 = np.cross(e3, e1)
            e2 /= np.linalg.norm(e2)
            E1.append(e1)
            E2.append(e2)
        return E1, E2

    def calculate_dpxs(self, w, k) -> np.ndarray:
        DPXS = []
        _J = self.JS_inv[k]
        for psi, dpsiz, th in zip(self._p[k], self.dpz[k].T, self.th):
            dpz = np.zeros(6)
            dpz[0:2] = dpsiz
            dpz[3] = 1/2*w*th*dpsiz[0]
            dpz[4] = 1/2*w*th*dpsiz[1]
            dpz[5] = 1/2*th*psi

            dpx1 = _J @ dpz[:3]
            dpx2 = _J @ dpz[3:6]
            DPXS.append([*dpx1, *dpx2])
        return np.array(DPXS)

    def calculate_BNL(self, dpx) -> np.ndarray:
        n = len(self.coords)
        BNL = np.zeros((9, 5*n))
        _e1, _e2 = self.calculate_e1_e2(deformed=True)
        for i in range(n):
            e1 = _e1[i]
            e2 = _e2[i]
            __dpx = dpx[:, i]
            _dpx = __dpx[:3]
            _dpe = __dpx[3:6]

            BNL[0:3, 5*i] = _dpx
            BNL[3:6, 5*i+1] = _dpx
            BNL[6:9, 5*i+2] = _dpx

            BNL[0:3, 5*i+3] = _dpe*e1[0]
            BNL[3:6, 5*i+3] = _dpe*e1[1]
            BNL[6:9, 5*i+3] = _dpe*e1[2]

            BNL[0:3, 5*i+4] = -_dpe*e2[0]
            BNL[3:6, 5*i+4] = -_dpe*e2[1]
            BNL[6:9, 5*i+4] = -_dpe*e2[2]
        return BNL

    def calculate_BL(self, dpx) -> np.ndarray:
        return

    def organize_S(self, S) -> tuple:
        return

    def transformation_matrix(self, deformed=True) -> np.ndarray:
        return

    def rotation_matrix(self, deformed=True) -> np.ndarray:
        return np.eye(5)

    def get_local_jacobian(self, jac: np.ndarray, dni: np.array) -> np.ndarray:
        return

    def project_coords(self, deformed=True) -> None:
        R = self.rotation_matrix(deformed)[:3, :3]
        coords = self.coords
        self.t_coords = coords@R

    def project_disp(self, deformed=True) -> None:
        """Project the displacements to the original coordinates system."""
        R = self.rotation_matrix(deformed)
        self.t_Ue = self.Ue.T@R
        self.t_Ue = self.t_Ue.T

    def _get_spatial_derivatives(self, w, deformed=False):
        J = []
        coords = self.coords.copy() + self.Ue[:3].T*deformed

        for psi, dpsiz in zip(self._p, self.dpz):
            j = 0.0
            dx0dz = 0.0
            dx0dn = 0.0
            dx0dw = 0.0
            for i in range(len(coords)):
                e3 = self.e3[i] + (self.Ue[3][i]*self.e1[i] -
                                   self.Ue[4][i]*self.e2[i])*deformed
                dx0dz += dpsiz[0][i] * \
                    (coords[i] + w/2*self.th[i]*e3)
                dx0dn += dpsiz[1][i] * \
                    (coords[i] + w/2*self.th[i]*e3)
                dx0dw += 1/2*psi[i]*self.th[i]*e3
            j = np.zeros((3, 3))
            j[:, 0] = dx0dz
            j[:, 1] = dx0dn
            j[:, 2] = dx0dw
            J.append(j)
        return J

    def get_jacobians(self, w):
        self.JS = self._get_spatial_derivatives(w, deformed=False)
        self.JS_inv = [np.linalg.inv(j) for j in self.JS]

    def calculate_deformation_gradients(self, w):
        """Calculate the deformation gradients for the shell element."""
        self.FS = []
        JS2 = self._get_spatial_derivatives(w, deformed=True)
        for i in range(len(self.JS)):
            J_inv = self.JS_inv[i]  # Jacobian with undeformed coordinates
            J2 = JS2[i]  # Jacobian with deformed coordinates
            col1 = J_inv @ J2[:, 0]
            col2 = J_inv @ J2[:, 1]
            col3 = J_inv @ J2[:, 2]
            F = np.column_stack((col1, col2, col3))
            self.FS.append(F)


class QuadShellLinear(ShellBase, Quadrilateral):
    def __init__(self, coords: np.ndarray, gdl: np.ndarray, **kargs):
        Quadrilateral.__init__(self, coords, gdl, 3, **kargs)
        ShellBase.__init__(self, **kargs)


class TriShelleLinear(ShellBase, LTriangular):
    def __init__(self, coords: np.ndarray, gdl: np.ndarray, **kargs):
        LTriangular.__init__(self, coords, gdl, 2, **kargs)
        ShellBase.__init__(self, **kargs)
