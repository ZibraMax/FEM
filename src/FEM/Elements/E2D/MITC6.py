
from .Element2D import Element2D, np
from .TriangularScheme import TriangularScheme
from ..E1D.QuadraticElement import QuadraticElement


def shape_functions_MITC6(r, s):
    N = np.array([
        (1 - r - s)*(1 - 2*r - 2*s),
        r*(2*r - 1),
        s*(2*s - 1),
        4*r*(1 - r - s),
        4*r*s,
        4*s*(1 - r - s)
    ])

    dNdr = np.array([
        -3 + 4*r + 4*s,
        4*r - 1,
        0,
        4*(1 - 2*r - s),
        4*s,
        -4*s
    ])

    dNds = np.array([
        -3 + 4*r + 4*s,
        0,
        4*s - 1,
        -4*r,
        4*r,
        4*(1 - r - 2*s)
    ])

    return N, dNdr, dNds


def assemble_B_membrane(dNdxy):
    B = np.zeros((3, 30))
    for i in range(6):
        B[0, i*5 + 0] = dNdxy[i, 0]  # du/dx
        B[1, i*5 + 1] = dNdxy[i, 1]  # dv/dy
        B[2, i*5 + 0] = dNdxy[i, 1]  # du/dy
        B[2, i*5 + 1] = dNdxy[i, 0]  # dv/dx
    return B


def assemble_B_bending(dNdxy):
    B = np.zeros((3, 30))
    for i in range(6):
        B[0, i*5 + 3] = dNdxy[i, 0]  # d(theta_x)/dx
        B[1, i*5 + 4] = dNdxy[i, 1]  # d(theta_y)/dy
        B[2, i*5 + 3] = dNdxy[i, 1]  # d(theta_x)/dy
        B[2, i*5 + 4] = dNdxy[i, 0]  # d(theta_y)/dx
    return B


def assemble_B_shear(dNdxy, N):
    B = np.zeros((2, 30))
    for i in range(6):
        # Shear strain gamma_xz depends on w,z and rotation theta_y
        B[0, i*5 + 2] = dNdxy[i, 0]    # dw/dx
        B[0, i*5 + 4] = -N[i]         # -theta_y

        # Shear strain gamma_yz depends on w,z and rotation theta_x
        B[1, i*5 + 2] = dNdxy[i, 1]    # dw/dy
        B[1, i*5 + 3] = -N[i]         # -theta_x
    return B


def compute_local_coordinates(X):
    v1 = X[1, :] - X[0, :]
    v2 = X[2, :] - X[0, :]
    e3 = np.cross(v1, v2)
    e3 /= np.linalg.norm(e3)
    e1 = v1 / np.linalg.norm(v1)
    e2 = np.cross(e3, e1)
    e2 /= np.linalg.norm(e2)
    basis = np.vstack((e1, e2, e3)).T
    X_local = np.zeros((6, 2))
    for i in range(6):
        dX = X[i, :] - X[0, :]
        X_local[i, 0] = np.dot(dX, e1)
        X_local[i, 1] = np.dot(dX, e2)
    return X_local, basis


def build_global_transformation(basis):
    n_nodes = 6
    n_dof = 5
    T = np.eye(n_nodes * n_dof)

    for i in range(n_nodes):
        idx = i * n_dof

        # Rotate translation part (u, v, w)
        T[idx:idx+3, idx:idx+3] = basis

        # Rotations (theta_x, theta_y) are left identity (local rotations)

    return T


def get_D_matrices(E=210e9, nu=0.3, t=0.01, kappa=5/6):
    C = E / (1 - nu**2)
    D_membrane = C * np.array([
        [1, nu, 0],
        [nu, 1, 0],
        [0, 0, (1-nu)/2]
    ])
    D_bending = (E * t**3 / (12*(1 - nu**2))) * np.array([
        [1, nu, 0],
        [nu, 1, 0],
        [0, 0, (1-nu)/2]
    ])
    G = E / (2*(1+nu))
    D_shear = kappa * G * t * np.array([
        [1, 0],
        [0, 1]
    ])
    return D_membrane, D_bending, D_shear


def assemble_Ktangent_Fint(x_current, u_disp, D_membrane, D_bending, D_shear):

    n_nodes = 6
    n_dof = 5
    K_tangent = np.zeros((n_nodes*n_dof, n_nodes*n_dof))
    F_int = np.zeros(n_nodes*n_dof)

    # Get local coordinates and basis (to compute Jacobian)
    X_local, basis = compute_local_coordinates(x_current)

    # Gauss integration points and weights
    gauss_points = np.array([
        [1/6, 1/6],
        [2/3, 1/6],
        [1/6, 2/3]
    ])
    weights = np.array([1/3, 1/3, 1/3])

    for gp, w in zip(gauss_points, weights):
        r, s = gp
        N, dNdr, dNds = shape_functions_MITC6(r, s)

        # Jacobian J
        J = np.zeros((2, 2))
        for i in range(n_nodes):
            J[0, 0] += dNdr[i] * X_local[i, 0]
            J[0, 1] += dNdr[i] * X_local[i, 1]
            J[1, 0] += dNds[i] * X_local[i, 0]
            J[1, 1] += dNds[i] * X_local[i, 1]

        detJ = np.linalg.det(J)
        invJ = np.linalg.inv(J)

        # Shape function derivatives in physical (x,y)
        dNdxy = np.zeros((n_nodes, 2))
        for i in range(n_nodes):
            dNdxy[i, :] = invJ @ np.array([dNdr[i], dNds[i]])

        # Assemble B-matrices
        Bm = assemble_B_membrane(dNdxy)
        Bb = assemble_B_bending(dNdxy)
        Bs = assemble_B_shear(dNdxy, N)

        # Extract displacements and rotations
        strain_membrane = Bm @ u_disp
        strain_bending = Bb @ u_disp
        strain_shear = Bs @ u_disp

        # Stress (linear elastic)
        stress_membrane = D_membrane @ strain_membrane
        stress_bending = D_bending @ strain_bending
        stress_shear = D_shear @ strain_shear

        # Internal force
        F_int += (Bm.T @ stress_membrane + Bb.T @
                  stress_bending + Bs.T @ stress_shear) * detJ * w

        # Material tangent stiffness
        K_material = (Bm.T @ D_membrane @ Bm +
                      Bb.T @ D_bending @ Bb +
                      Bs.T @ D_shear @ Bs) * detJ * w

        # Geometric stiffness (only membrane stresses contribute)
        K_geo = np.zeros((n_nodes*n_dof, n_nodes*n_dof))
        sigma_x = stress_membrane[0]
        sigma_y = stress_membrane[1]
        sigma_xy = stress_membrane[2]

        for i in range(n_nodes):
            for j in range(n_nodes):
                # Only translation DOFs (u,v) contribute to geometric stiffness
                value = (sigma_x * dNdxy[i, 0] * dNdxy[j, 0] +
                         sigma_y * dNdxy[i, 1] * dNdxy[j, 1] +
                         sigma_xy * (dNdxy[i, 0]*dNdxy[j, 1] + dNdxy[i, 1]*dNdxy[j, 0]))
                K_geo[i*n_dof+0, j*n_dof+0] += value * detJ * w  # u_x-u_x
                K_geo[i*n_dof+1, j*n_dof+1] += value * detJ * w  # u_y-u_y

        # Accumulate
        K_tangent += K_material + K_geo
    F_int = F_int.reshape([30, 1])
    return K_tangent, F_int


class MITC6(Element2D, TriangularScheme):

    def __init__(self, coords: np.ndarray, gdl: np.ndarray, n: int = 3, **kargs) -> None:
        self.internal_force = None
        coords = np.array(coords)
        _coords = np.array([coords[i] for i in range(3)])
        TriangularScheme.__init__(self, n, **kargs)
        self.boundaries = []
        Element2D.__init__(self, coords, _coords, gdl, **kargs)

    def set_material(self, E: float, nu: float, t: float) -> None:
        self.E = E
        self.nu = nu
        self.t = t
        self.kappa = 5/6
        self.D_membrane, self.D_bending, self.D_shear = get_D_matrices(
            E, nu, t, self.kappa)

    def psis(self, z: np.ndarray) -> np.ndarray:
        return np.array([
            2.0*(z[0]+z[1]-1.0)*(z[0]+z[1]-0.5),
            2.0*z[0]*(z[0]-0.5),
            2.0*z[1]*(z[1]-0.5),
            -4.0*(z[0]+z[1]-1.0)*(z[0]),
            4.0*z[0]*z[1],
            -4.0*z[1]*(z[0]+z[1]-1.0)]).T

    def dpsis(self, z: np.ndarray) -> np.ndarray:

        return np.array([
            [4.0*z[0]+4.0*z[1]-3.0, 4.0*z[1]+4.0*z[0]-3.0],
            [4.0*z[0]-1.0, 0*z[0]],
            [0*z[0], 4.0*z[1]-1.0],
            [-8.0*z[0]-4.0*(z[1]-1.0), -4.0*z[0]],
            [4.0*z[1], 4.0*z[0]],
            [-4.0*z[1], -8.0*z[1]-4.0*z[0]+4.0]
        ])

    def elementMatrixNonLineal(self):
        Ue = self.Ue.T.flatten().reshape([6, 5])
        coords = self.coords
        # Update coordinates with displacements
        x_current = coords  # + Ue[:, :3]
        u_disp = self.Ue.T.flatten()
        D_membrane, D_bending, D_shear = self.D_membrane, self.D_bending, self.D_shear
        return assemble_Ktangent_Fint(x_current, u_disp, D_membrane, D_bending, D_shear)
