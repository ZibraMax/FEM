
from .Quadrilateral import Quadrilateral, np


# Esto debería ser un "template" en lugar de una clase concreta
class Membrane(Quadrilateral):

    def __init__(self, coords: np.ndarray, gdl: np.ndarray, n: int = 2, **kargs) -> None:
        Quadrilateral.__init__(self, coords, gdl, n, **kargs)

    def assign_properties(self, E: float, nu: float, t: float):
        self.E = E
        self.nu = nu
        self.t = t

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

    def stiffnessMatrix(self):
        Ke = 0.0
        _x, _p = self.T(self.Z.T)
        _j, _dp = self.J(self.Z.T)
        weights = self.W
        t = self.t
        E = self.E
        nu = self.nu

        c11 = E / (1 - nu**2)  # Revisar
        c12 = nu * c11  # Revisar
        c22 = c11
        c66 = E / (2 * (1 + nu))
        C = np.array([
            [c11, c12, 0.0],
            [c12, c22, 0.0],
            [0.0, 0.0, c66]])
        for x, jac, wi, ni, dni in zip(_x, _j, weights, _p, _dp):
            JP = np.linalg.pinv(jac)
            e3 = np.cross(*jac, axis=0)
            detjac = np.linalg.norm(e3)
            e1 = jac[0].copy()  # This is important, creo
            e2 = jac[1].copy()  # This is important, creo

            # Ortonormalizar la cosa
            e1 /= np.linalg.norm(e1)
            e2 = np.cross(e3, e1)
            e2 /= np.linalg.norm(e2)

            # Matriz de proyeccion
            # 3x2 transformation to local plane
            T_plane = np.vstack((e1, e2)).T

            # En teoria estas son las derivadas de las funciones de forma pero en X, Y, Z
            dpx = JP @ dni
            # Ahora hay que proyectar esas derivadas en el plano usando la base ortonormal
            # Entonces ahoera estas derivadas estan en el plano de la menmbrana y se pueden usar para calcular rigidez
            dpt = dpx.T @ T_plane
            dpt = dpt.T  # Transpose porque si

            B = self.compute_B_matrix(dpt[0], dpt[1])
            T = self.build_global_transform(T_plane, len(self.gdl.T))
            Ke_local = t*(B.T@C@B)*detjac*wi
            Ke += T @ Ke_local @ T.T
        return Ke
