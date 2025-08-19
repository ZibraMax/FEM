from .ContinumBase import ContinumBase
from ..E2D import LTriangular, Quadrilateral
import numpy as np


class ShellBase(ContinumBase):

    def __init__(self, **kargs):
        ContinumBase.__init__(self, **kargs)

    def set_thickness(self, t: float, e3: np.ndarray = None, n_gauss_thickness: int = 1) -> None:

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
        self.e1, self.e2, _ = self.calculate_e1_e2()

    def calculate_e1_e2(self, deformed=False):
        """Calculate the e1 and e2 vectors for the shell element."""
        E1 = []
        E2 = []
        E3 = []
        for i in range(len(self.coords)):
            e3 = self.e3[i] + (self.Ue[3][i]*self.e1[i] -
                               self.Ue[4][i]*self.e2[i])*deformed
            e1 = np.cross([0, 1, 0], e3)
            e1 /= np.linalg.norm(e1)
            e2 = np.cross(e3, e1)
            e2 /= np.linalg.norm(e2)
            e3 /= np.linalg.norm(e3)
            E1.append(e1)
            E2.append(e2)
            E3.append(e3)
        return E1, E2, E3

    def calculate_dpxs(self, w, k, J_INV) -> np.ndarray:
        # w is  gauss point in thickness
        # k is the index of the gausss point in s,r
        # J_INV is the inverse of the jacobian at the gauss point
        dpnodedx = []
        dpnodedy = []
        dpnodedz = []
        dpbnodedx = []
        dpbnodedy = []
        dpbnodedz = []

        for i in range(len(self.coords)):

            dpnodedx.append(J_INV[0, 0]*self.dpz[k][0]
                            [i] + J_INV[0, 1]*self.dpz[k][1][i])
            dpnodedy.append(J_INV[1, 0]*self.dpz[k][0]
                            [i] + J_INV[1, 1]*self.dpz[k][1][i])
            dpnodedz.append(J_INV[2, 0]*self.dpz[k][0]
                            [i] + J_INV[2, 1]*self.dpz[k][1][i])
            hk = self.th[i] / 2
            dpbnodedx.append(hk*(J_INV[0, 2]*self._p[k][i] + w*dpnodedx[-1]))
            dpbnodedy.append(hk*(J_INV[1, 2]*self._p[k][i] + w*dpnodedy[-1]))
            dpbnodedz.append(hk*(J_INV[2, 2]*self._p[k][i] + w*dpnodedz[-1]))
        dpxs = [dpnodedx, dpnodedy, dpnodedz,
                dpbnodedx, dpbnodedy, dpbnodedz]
        dpxs = np.array(dpxs)
        return dpxs

    def calculate_BNL(self, dpxs) -> np.ndarray:
        dpnodedx, dpnodedy, dpnodedz, dpbnodedx, dpbnodedy, dpbnodedz = dpxs
        BNL = np.zeros((9, len(self.coords)*5))
        e1s, e2s, _ = self.calculate_e1_e2(deformed=True)
        for i in range(len(self.coords)):
            e1 = e1s[i]
            e2 = e2s[i]
            BNL[0, i*5] = dpnodedx[i]
            BNL[0, i*5+3] = dpbnodedx[i]*e1[0]
            BNL[0, i*5+4] = dpbnodedx[i]*e2[0]

            BNL[1, i*5+1] = dpnodedy[i]
            BNL[1, i*5+3] = dpbnodedy[i]*e1[1]
            BNL[1, i*5+4] = dpbnodedy[i]*e2[1]

            BNL[2, i*5+2] = dpnodedz[i]
            BNL[2, i*5+3] = dpbnodedz[i]*e1[2]
            BNL[2, i*5+4] = dpbnodedz[i]*e2[2]

            BNL[3, i*5] = dpnodedy[i]
            BNL[3, i*5+3] = dpbnodedy[i]*e1[0]
            BNL[3, i*5+4] = dpbnodedy[i]*e2[0]

            BNL[4, i*5+1] = dpnodedx[i]
            BNL[4, i*5+3] = dpbnodedx[i]*e1[1]
            BNL[4, i*5+4] = dpbnodedx[i]*e2[1]

            BNL[5, i*5] = dpnodedz[i]
            BNL[5, i*5+3] = dpbnodedz[i]*e1[0]
            BNL[5, i*5+4] = dpbnodedz[i]*e2[0]

            BNL[6, i*5+2] = dpnodedx[i]
            BNL[6, i*5+3] = dpbnodedx[i]*e1[2]
            BNL[6, i*5+4] = dpbnodedx[i]*e2[2]

            BNL[7, i*5+1] = dpnodedz[i]
            BNL[7, i*5+3] = dpbnodedz[i]*e1[1]
            BNL[7, i*5+4] = dpbnodedz[i]*e2[1]

            BNL[8, i*5+2] = dpnodedy[i]
            BNL[8, i*5+3] = dpbnodedy[i]*e1[2]
            BNL[8, i*5+4] = dpbnodedy[i]*e2[2]

            BNL[:, i*5+4] *= -1

        return BNL

    def calculate_ls(self, F) -> np.ndarray:
        e1s, e2s, _ = self.calculate_e1_e2(deformed=True)

        lij = F
        lb1 = []
        lb2 = []
        for i in range(len(self.coords)):
            e1 = e1s[i]
            e2 = e2s[i]
            _lb1x = lij[0]@e1
            _lb1y = lij[1]@e1
            _lb1z = lij[2]@e1
            _lb1 = np.array([_lb1x, _lb1y, _lb1z])
            _lb2x = lij[0]@e2
            _lb2y = lij[1]@e2
            _lb2z = lij[2]@e2
            _lb2 = np.array([_lb2x, _lb2y, _lb2z])
            lb1.append(_lb1)
            lb2.append(_lb2)
        return lij, lb1, lb2

    def calculate_BL(self, dpxs, F) -> np.ndarray:
        dpnodedx, dpnodedy, dpnodedz, dpbnodedx, dpbnodedy, dpbnodedz = dpxs
        lij, lb1, lb2 = self.calculate_ls(F)
        n = len(self.coords)
        BL = np.zeros((6, 5*n))
        for i in range(n):
            _lb1 = lb1[i]
            _lb2 = lb2[i]

            BL[0, 5*i] = dpnodedx[i]*lij[0, 0]
            BL[0, 5*i+1] = dpnodedx[i]*lij[1, 0]
            BL[0, 5*i+2] = dpnodedx[i]*lij[2, 0]
            BL[0, 5*i+3] = dpbnodedx[i]*_lb1[0]
            BL[0, 5*i+4] = dpbnodedx[i]*_lb2[0]

            BL[1, 5*i] = dpnodedy[i]*lij[0, 1]
            BL[1, 5*i+1] = dpnodedy[i]*lij[1, 1]
            BL[1, 5*i+2] = dpnodedy[i]*lij[2, 1]
            BL[1, 5*i+3] = dpbnodedy[i]*_lb1[1]
            BL[1, 5*i+4] = dpbnodedy[i]*_lb2[1]

            BL[2, 5*i] = dpnodedz[i]*lij[0, 2]
            BL[2, 5*i+1] = dpnodedz[i]*lij[1, 2]
            BL[2, 5*i+2] = dpnodedz[i]*lij[2, 2]
            BL[2, 5*i+3] = dpbnodedz[i]*_lb1[2]
            BL[2, 5*i+4] = dpbnodedz[i]*_lb2[2]

            BL[3, 5*i] = dpnodedy[i]*lij[0, 0] + dpnodedx[i]*lij[0, 1]
            BL[3, 5*i+1] = dpnodedy[i]*lij[1, 0] + dpnodedx[i]*lij[1, 1]
            BL[3, 5*i+2] = dpnodedy[i]*lij[2, 0] + dpnodedx[i]*lij[2, 1]

            BL[3, 5*i+3] = dpbnodedy[i]*_lb1[0] + dpbnodedx[i]*_lb1[1]
            BL[3, 5*i+4] = dpbnodedy[i]*_lb2[0] + dpbnodedx[i]*_lb2[1]

            BL[4, 5*i] = dpnodedz[i]*lij[0, 0] + dpnodedx[i]*lij[0, 2]
            BL[4, 5*i+1] = dpnodedz[i]*lij[1, 0] + dpnodedx[i]*lij[1, 2]
            BL[4, 5*i+2] = dpnodedz[i]*lij[2, 0] + dpnodedx[i]*lij[2, 2]
            BL[4, 5*i+3] = dpbnodedz[i]*_lb1[0] + dpbnodedx[i]*_lb1[2]
            BL[4, 5*i+4] = dpbnodedz[i]*_lb2[0] + dpbnodedx[i]*_lb2[2]

            BL[5, 5*i] = dpnodedz[i]*lij[0, 1] + dpnodedy[i]*lij[0, 2]
            BL[5, 5*i+1] = dpnodedz[i]*lij[1, 1] + dpnodedy[i]*lij[1, 2]
            BL[5, 5*i+2] = dpnodedz[i]*lij[2, 1] + dpnodedy[i]*lij[2, 2]
            BL[5, 5*i+3] = dpbnodedz[i]*_lb1[1] + dpbnodedy[i]*_lb1[2]
            BL[5, 5*i+4] = dpbnodedz[i]*_lb2[1] + dpbnodedy[i]*_lb2[2]

            BL[:, i*5+4] *= -1

        return BL

    def organize_S(self, S) -> tuple:
        # S is 3x3

        S_force = np.zeros((6, 1))
        S_force[0, 0] = S[0, 0]  # S11
        S_force[1, 0] = S[1, 1]  # S22
        S_force[2, 0] = S[2, 2]  # S33
        S_force[3, 0] = S[0, 1]  # S[1, 0] = S[0, 1]  # S12
        S_force[4, 0] = S[0, 2]  # S[0, 2]  # S13
        S_force[5, 0] = S[1, 2]

        mS = np.zeros((9, 9))

        # Fill according to the pattern in the image
        mS[0, 0] = S[0, 0]
        mS[0, 3] = S[0, 1]
        mS[0, 5] = S[0, 2]

        mS[1, 1] = S[1, 1]
        mS[1, 4] = S[1, 0]
        mS[1, 7] = S[1, 2]

        mS[2, 2] = S[2, 2]
        mS[2, 6] = S[2, 0]
        mS[2, 8] = S[2, 1]

        mS[3, 3] = S[1, 1]
        mS[3, 5] = S[1, 2]

        mS[4, 4] = S[0, 0]
        mS[4, 7] = S[0, 2]

        mS[5, 5] = S[2, 2]

        mS[6, 6] = S[0, 0]
        mS[6, 8] = S[0, 1]

        mS[7, 7] = S[2, 2]

        mS[8, 8] = S[1, 1]

        # Make symmetric
        mS = mS + mS.T - np.diag(np.diag(mS))

        return mS, S_force

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
        FS = []
        JS2 = self._get_spatial_derivatives(w, deformed=True)
        for i in range(len(self.JS)):
            J_inv = self.JS_inv[i]  # Jacobian with undeformed coordinates
            J2 = JS2[i]  # Jacobian with deformed coordinates
            col1 = J_inv @ J2[:, 0]
            col2 = J_inv @ J2[:, 1]
            col3 = J_inv @ J2[:, 2]
            F = np.column_stack((col1, col2, col3))
            FS.append(F)
        return FS

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

                dpx = self.calculate_dpxs(z_t, gi, _JS[gi])
                BNL = self.calculate_BNL(dpx)
                F = FS[gi]
                E = self.green_lagrange_strain(F)
                C, S = self.constitutive_model(E)
                S_stiff, S_force = self.organize_S(S)

                BL = self.calculate_BL(dpx, F)
                Ke += (BL.T @ C @ BL + BNL.T @
                       S_stiff @ BNL) * detjac * wi
                Fe += (BL.T @ S_force) * detjac * wi

        # T1 = self.transformation_matrix(True)
        # T2 = self.transformation_matrix(True)
        # Ke = T1.T @ Ke @ T1
        # Fe = T2.T @ Fe
        return Ke, Fe

    def x(self, r, s, t, deformed=False):
        """Spatial transformation, from a gauss point (r,s,t) to an spatial point X,Y,Z
        """
        res = 0.0
        z = np.array([[r], [s]])
        _, phi = self.T(z)
        phi = phi[0]
        coords = self.coords + deformed*self.Ue[:3].T
        e3 = self.e3
        if deformed:
            _, _, e3 = self.calculate_e1_e2(deformed)

        for k in range(len(self.coords)):
            res += phi[k]*coords[k] + 1/2*self.th[k]*t*phi[k]*e3[k]
        return res

    def numerical_jacobian(self, r, s, t, deformed=False, h=0.0000001):
        dX = (self.x(r+h, s, t, deformed) - self.x(r-h, s, t, deformed))/(2*h)
        dY = (self.x(r, s+h, t, deformed) - self.x(r, s-h, t, deformed))/(2*h)
        dZ = (self.x(r, s, t+h, deformed) - self.x(r, s, t-h, deformed))/(2*h)

        return np.array([dX, dY, dZ]).T

    def u1(self, r, s, t):
        """Displacement at a gauss point (r,s,t)"""
        x0 = self.x(r, s, t, deformed=False)
        x1 = self.x(r, s, t, deformed=True)
        u = x1 - x0
        return u

    def u_incr(self, r, s, t, Ue):
        U1 = self.Ue.copy()
        U21 = Ue.copy()
        x_deformed = self.x(r, s, t, deformed=True)
        self.Ue = U1 + U21
        x_deformed_incr = self.x(r, s, t, deformed=True)
        self.Ue = U1
        return x_deformed_incr - x_deformed

    def u_incr_deriv(self, r, s, t, Ue, h=1e-8):
        du_dr = (self.u_incr(r+h, s, t, Ue) -
                 self.u_incr(r-h, s, t, Ue)) / (2*h)
        du_ds = (self.u_incr(r, s+h, t, Ue) -
                 self.u_incr(r, s-h, t, Ue)) / (2*h)
        du_dt = (self.u_incr(r, s, t+h, Ue) -
                 self.u_incr(r, s, t-h, Ue)) / (2*h)
        return np.array([du_dr, du_ds, du_dt]).T

    def u_deriv(self, r, s, t, h=1e-8):
        du_dr = (self.u1(r+h, s, t) - self.u1(r-h, s, t)) / (2*h)
        du_ds = (self.u1(r, s+h, t) - self.u1(r, s-h, t)) / (2*h)
        du_dt = (self.u1(r, s, t+h) - self.u1(r, s, t-h)) / (2*h)
        return np.array([du_dr, du_ds, du_dt]).T

    def u1_analytical(self, r, s, t):
        """Displacement at a gauss point (r,s,t) using analytical shape functions"""
        res = 0.0
        z = np.array([[r], [s]])
        _, phis = self.T(z)
        phis = phis[0]
        _, _, e3_0 = self.calculate_e1_e2(deformed=False)
        _, _, e3_1 = self.calculate_e1_e2(deformed=True)
        for i in range(len(self.coords)):
            res += phis[i]*self.Ue[:3, i] + 1/2*t * \
                self.th[i]*phis[i]*(e3_1[i]-e3_0[i])
        return res

    def du1_analytical_deriv(self, r, s, t, h=1e-8):
        du_dr = (self.u1_analytical(r+h, s, t) -
                 self.u1_analytical(r-h, s, t)) / (2*h)
        du_ds = (self.u1_analytical(r, s+h, t) -
                 self.u1_analytical(r, s-h, t)) / (2*h)
        du_dt = (self.u1_analytical(r, s, t+h) -
                 self.u1_analytical(r, s, t-h)) / (2*h)
        return np.array([du_dr, du_ds, du_dt]).T


class QuadShellLinear(ShellBase, Quadrilateral):
    def __init__(self, coords: np.ndarray, gdl: np.ndarray, **kargs):
        Quadrilateral.__init__(self, coords, gdl, 5, **kargs)
        ShellBase.__init__(self, **kargs)


class TriShelleLinear(ShellBase, LTriangular):
    def __init__(self, coords: np.ndarray, gdl: np.ndarray, **kargs):
        LTriangular.__init__(self, coords, gdl, 2, **kargs)
        ShellBase.__init__(self, **kargs)
