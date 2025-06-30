
from .Quadrilateral import Quadrilateral, np


# Esto debería ser un "template" en lugar de una clase concreta
class Membrane(Quadrilateral):

    def __init__(self, coords: np.ndarray, gdl: np.ndarray, n: int = 2, **kargs) -> None:
        Quadrilateral.__init__(self, coords, gdl, n, **kargs)

    def assign_properties(self, E: float, nu: float, t: float):
        self.E = E
        self.nu = nu
        self.t = t
        self.c11 = self.E / (1 - self.nu**2)
        self.c12 = self.nu * self.c11
        self.c22 = self.c11
        self.c66 = self.E / (2 * (1 + self.nu))
        self.C = np.array([
            [self.c11, self.c12, 0.0],
            [self.c12, self.c22, 0.0],
            [0.0, 0.0, self.c66]])

    def compute_B_matrix(self, dNdx, dNdy):
        """Builds membrane B-matrix (plane stress assumption)"""
        n_nodes = dNdx.shape[0]
        B = np.zeros((3, 2 * n_nodes))
        for i in range(n_nodes):
            B[0, 2 * i] = dNdx[i]
            B[1, 2 * i + 1] = dNdy[i]
            B[2, 2 * i] = dNdy[i]
            B[2, 2 * i + 1] = dNdx[i]
        return B

    def build_global_transform(self, T_plane, num_nodes):
        T_e = np.zeros((3 * num_nodes, 2 * num_nodes))
        for i in range(num_nodes):
            row_start = 3 * i
            col_start = 2 * i
            T_e[row_start:row_start+3, col_start:col_start +
                2] = T_plane  # insert 3×2 block

        return T_e

    def calculate_deformation_gradient(self, dpt, T_plane):
        X_0 = self.coords.copy()
        X_t = X_0 + self.Ue.T
        t_t = X_t @ T_plane
        t_0 = X_0 @ T_plane
        F = dpt @ t_t
        Uet = t_t - t_0
        return F, Uet

    # Considerar mover esta fuction a un lugar más general
    def calculate_strains(self, F):
        """Calculates Green-Lagrange amd Cauchy strains from deformation gradient"""
        E = 0.5 * (F.T @ F - np.eye(F.shape[0]))
        e = 0.5 * (np.trace(E) * np.eye(E.shape[0]) + E)
        return e, E

    def calculate_stress(self, e, E):
        """Calculates stress from strains using constitutive matrix"""
        e = e.flatten()[[0, 3, 1]]
        E = E.flatten()[[0, 3, 1]]
        S_vec = self.C @ E
        S = np.zeros((2, 2))
        S[0, 0] = S_vec[0]
        S[1, 1] = S_vec[1]
        S[0, 1] = S_vec[2]
        S[1, 0] = S_vec[2]

        s_vec = self.C @ e
        s = np.zeros((2, 2))
        s[0, 0] = s_vec[0]
        s[1, 1] = s_vec[1]
        s[0, 1] = s_vec[2]
        s[1, 0] = s_vec[2]

        return s, S

    def compute_BL2_matrix(self, dNdx, dNdy, Uet):
        """Builds membrane B-L2 matrix (plane stress assumption)"""
        n_nodes = dNdx.shape[0]
        BL2 = np.zeros((3, 2 * n_nodes))
        Ux = Uet.T[0]
        Uy = Uet.T[1]

        l11 = np.sum(dNdx*Ux)
        l21 = np.sum(dNdx*Uy)
        l12 = np.sum(dNdy*Ux)
        l22 = np.sum(dNdy*Uy)
        L = np.array([[l11, l12],
                      [l21, l22]])

        for i in range(n_nodes):
            BL2[0, 2 * i] = dNdx[i] * L[0, 0]
            BL2[0, 2 * i + 1] = dNdx[i] * L[1, 0]

            BL2[1, 2 * i] = dNdy[i] * L[0, 1]
            BL2[1, 2 * i + 1] = dNdy[i] * L[1, 1]

            BL2[2, 2 * i] = dNdy[i] * L[0, 0] + dNdx[i] * L[0, 1]
            BL2[2, 2 * i + 1] = dNdy[i] * L[1, 0] + dNdx[i] * L[1, 1]
        return BL2

    def compute_BNL_matrix(self, dNdx, dNdy):
        """Builds membrane B-NL matrix (plane stress assumption)"""
        n_nodes = dNdx.shape[0]
        BNL = np.zeros((4, 2 * n_nodes))
        for i in range(n_nodes):
            BNL[0, 2 * i] = dNdx[i]
            BNL[1, 2 * i] = dNdy[i]

            BNL[2, 2 * i + 1] = dNdx[i]
            BNL[3, 2 * i + 1] = dNdy[i]

        return BNL

    def calculate_orthonormal_basis(self, jac):
        e3 = np.cross(*jac, axis=0)
        detjac = np.linalg.norm(e3)
        e1 = jac[0].copy()
        e2 = jac[1].copy()

        # Ortonormalizar la cosa
        e1 /= np.linalg.norm(e1)
        e2 = np.cross(e3, e1)
        e2 /= np.linalg.norm(e2)

        # Matriz de proyeccion
        # 3x2 transformation to local plane
        T_plane = np.vstack((e1, e2)).T
        return detjac, e1, e2, e3, T_plane

    def stiffnessMatrix(self):
        Ke = 0.0
        Fe = 0.0
        _x, _p = self.T(self.Z.T)
        _j, _dp = self.J(self.Z.T)
        _j_t = _dp @ (self.coords + self.Ue.T)
        weights = self.W
        t = self.t
        center_jac, dp_center = self.J(self.center.T)
        _, _, _, _, T_plane = self.calculate_orthonormal_basis(center_jac[0])
        center_jac_t = dp_center @ (self.coords + self.Ue.T)
        _, _, _, _, T_plane_t = self.calculate_orthonormal_basis(
            center_jac_t[0])

        for x, jac, jac_t, wi, ni, dni in zip(_x, _j, _j_t, weights, _p, _dp):
            JP = np.linalg.pinv(jac)  # Puede que sea otra cosa
            detjac, e1, e2, e3, _ = self.calculate_orthonormal_basis(jac)
            _, _, _, _, _ = self.calculate_orthonormal_basis(jac_t)
            # En teoria estas son las derivadas de las funciones de forma pero en X, Y, Z
            dpx = JP @ dni
            # Ahora hay que proyectar esas derivadas en el plano usando la base ortonormal
            # Entonces ahoera estas derivadas estan en el plano de la menmbrana y se pueden usar para calcular rigidez
            dpt = dpx.T @ T_plane
            dpt = dpt.T  # Transpose porque si

            # Calculamos el gradiente de deformación
            F, Uet = self.calculate_deformation_gradient(dpt, T_plane)
            # Calculamos las deformaciones
            e, E = self.calculate_strains(F)
            # Calculamos los esfuersos
            s, S = self.calculate_stress(e, E)

            BL1 = self.compute_B_matrix(*dpt)
            BL2 = self.compute_BL2_matrix(*dpt, Uet)
            BNL = self.compute_BNL_matrix(*dpt)
            BL = BL1 + BL2
            T = self.build_global_transform(T_plane, len(self.gdl.T))
            T_t = self.build_global_transform(T_plane_t, len(self.gdl.T))
            Ke_l_local = t*(BL.T @ self.C @ BL)*detjac*wi
            S_nl_int = np.zeros((4, 4))
            S_nl_int[0:2, 0:2] = S
            S_nl_int[2:4, 2:4] = S

            # T = T_t

            Ke_nl_local = t*(BNL.T @ S_nl_int @ BNL)*detjac*wi

            S_vec = S.flatten()[[0, 3, 1]]

            Fe_local = t*(BL.T @ S_vec)*detjac*wi

            Ke_local = Ke_l_local + Ke_nl_local
            Ke += T @ Ke_local @ T.T
            Fe += T @ Fe_local
        return Ke, Fe
