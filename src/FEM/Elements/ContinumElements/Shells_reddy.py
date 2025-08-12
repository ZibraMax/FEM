from .ContinumBase import ContinumBase
from ..E2D import LTriangular, Quadrilateral
import numpy as np


class ShellBase(ContinumBase):

    def __init__(self, **kargs):
        ContinumBase.__init__(self, **kargs)

    def set_thickness(self, t: float, e3: np.ndarray = None, n_gauss_thickness: int = 2) -> None:
        """Set the thickness of the shell element.

        Args:
            t (float): Thickness of the shell element.
        """

        # If thickness is not a list, create a list with n elements, n is the number of nodes
        self.n_gauss_thickness = n_gauss_thickness
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
            __dpx = dpx[i]
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

    def calculate_A(self, F):
        dudx = F+np.eye(3)
        A = np.zeros((6, 9))
        A[0, 0] = dudx[0, 0]
        A[0, 3] = dudx[0, 1]
        A[0, 6] = dudx[0, 2]

        A[1, 1] = dudx[1, 0]
        A[1, 4] = dudx[1, 1]
        A[1, 7] = dudx[1, 2]

        A[2, 2] = dudx[2, 0]
        A[2, 5] = dudx[2, 1]
        A[2, 8] = dudx[2, 2]

        fila = np.array([dudx[0, 1], dudx[0, 0], 0, dudx[1, 1],
                         dudx[1, 0], 0, dudx[2, 1], dudx[2, 0], 0])
        A[3, :] = fila

        fila = np.array([dudx[0, 2], 0, dudx[0, 0], dudx[1, 2],
                        0, dudx[1, 0], dudx[2, 2], 0, dudx[2, 0]])
        A[4, :] = fila

        fila = np.array([0, dudx[0, 2], dudx[0, 1], 0,
                        dudx[1, 2], dudx[1, 1], 0, dudx[2, 2], dudx[2, 1]])
        A[5, :] = fila

        return 1/2*A

    def calculate_BL(self, dpx) -> np.ndarray:

        return

    def organize_S(self, S) -> tuple:
        S_stiff = np.zeros((9, 9))

        S_stiff[0:3, 0:3] = S
        S_stiff[3:6, 3:6] = S
        S_stiff[6:9, 6:9] = S

        # Create 6x1 Voigt stress vector
        S_vect = np.zeros((6, 1))
        S_vect[0] = S[0, 0]  # S11
        S_vect[1] = S[1, 1]  # S22
        S_vect[2] = S[2, 2]  # S33
        S_vect[3] = S[0, 1]  # S12
        S_vect[4] = S[1, 2]  # S23
        S_vect[5] = S[0, 2]  # S13

        return S_stiff, S_vect

    def transformation_matrix(self, deformed=True) -> np.ndarray:
        return 1

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
        return self.JS, self.JS_inv

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
        return self.FS

    def elementMatrices(self) -> None:
        """Calculate element matrices and vectors.

        This method should be implemented in derived classes.
        """

        """Calculate the lineal stiffness matrix.

        This method should be implemented in derived classes.
        """
        weights = self.W
        t_z, t_w = np.polynomial.legendre.leggauss(self.n_gauss_thickness)
        Ke = 0.0
        Fe = 0.0
        for z_t, w_t in zip(t_z, t_w):
            JS, _JS = self.get_jacobians(z_t)
            FS = self.calculate_deformation_gradients(z_t)

            for gi in range(len(weights)):
                detjac = np.linalg.det(JS[gi])
                wi = weights[gi] * w_t
                F = FS[gi]
                E = self.green_lagrange_strain(F)
                C, S = self.constitutive_model(E)
                const = 1
                S_stiff, S_force = self.organize_S(S)
                dpx = self.calculate_dpxs(z_t, gi)
                # Estos dpx pueden ser diferentes porque toca rotarlos. Esperemos un poco

                BNL = self.calculate_BNL(dpx)
                A = self.calculate_A(F)
                BL = A @ BNL

                Ke += const*(BL.T @ C @ BL + BNL.T @
                             S_stiff @ BNL) * detjac * wi
                Fe += const*(BL.T @ S_force) * detjac * wi

        # T1 = self.transformation_matrix(True)
        # T2 = self.transformation_matrix(True)
        # Ke = T1.T @ Ke @ T1
        # Fe = T2.T @ Fe
        return Ke, Fe


class QuadShellLinearReddy(ShellBase, Quadrilateral):
    def __init__(self, coords: np.ndarray, gdl: np.ndarray, **kargs):
        Quadrilateral.__init__(self, coords, gdl, 3, **kargs)
        ShellBase.__init__(self, **kargs)


class TriShelleLinearReddy(ShellBase, LTriangular):
    def __init__(self, coords: np.ndarray, gdl: np.ndarray, **kargs):
        LTriangular.__init__(self, coords, gdl, 2, **kargs)
        ShellBase.__init__(self, **kargs)
