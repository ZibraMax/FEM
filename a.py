import numpy as np


def shape_functions_Q4(xi, eta):
    """Q4 bilinear shape functions and derivatives w.r.t xi, eta"""
    N = 0.25 * np.array([
        (1 - xi) * (1 - eta),
        (1 + xi) * (1 - eta),
        (1 + xi) * (1 + eta),
        (1 - xi) * (1 + eta)
    ])

    dN_dxi = 0.25 * np.array([
        [-(1 - eta), -(1 - xi)],
        [(1 - eta), -(1 + xi)],
        [(1 + eta),  (1 + xi)],
        [-(1 + eta),  (1 - xi)]
    ])

    return N, dN_dxi


def compute_B_matrix(dNdx, dNdy):
    """Builds membrane B-matrix (plane stress assumption)"""
    n_nodes = dNdx.shape[0]
    B = np.zeros((3, 2 * n_nodes))
    for i in range(n_nodes):
        B[0, 2 * i] = dNdx[i]
        B[1, 2 * i + 1] = dNdy[i]
        B[2, 2 * i] = dNdy[i]
        B[2, 2 * i + 1] = dNdx[i]
    return B


def membrane_stiffness_Q4_3D(coords, thickness, E, nu):
    # 2x2 Gauss quadrature points and weights
    gauss_points = [(-1/np.sqrt(3), -1/np.sqrt(3)),
                    (1/np.sqrt(3), -1/np.sqrt(3)),
                    (1/np.sqrt(3),  1/np.sqrt(3)),
                    (-1/np.sqrt(3),  1/np.sqrt(3))]
    weights = [1, 1, 1, 1]

    D = (E / (1 - nu**2)) * np.array([
        [1, nu, 0],
        [nu, 1, 0],
        [0, 0, (1 - nu) / 2]
    ])

    K = np.zeros((8, 8))  # 4 nodes × 2 dof (in-plane)

    for (xi, eta), w in zip(gauss_points, weights):
        N, dN_dxi = shape_functions_Q4(xi, eta)

        # Compute Jacobian J (3x2): maps reference to 3D
        J = np.zeros((3, 2))
        for i in range(4):
            J[:, 0] += dN_dxi[i, 0] * coords[i]
            J[:, 1] += dN_dxi[i, 1] * coords[i]

        # ✅ Copy J before modifying it
        J_plane = J.copy()

        # Build local orthonormal basis using J_plane
        e1 = J_plane[:, 0]
        e2 = J_plane[:, 1]
        normal = np.cross(e1, e2)
        area_factor = np.linalg.norm(normal)
        n = normal / (area_factor + 1e-12)  # unit normal

        e1 /= np.linalg.norm(e1)
        e2 = np.cross(n, e1)
        e2 /= np.linalg.norm(e2)

        T_plane = np.vstack((e1, e2)).T  # 3x2 transformation to local plane

        # Now compute the projection matrix using (possibly modified) J
        JtJ_inv = np.linalg.pinv(J.T @ J)
        dXdxi = J @ JtJ_inv @ dN_dxi.T  # shape fn derivatives w.r.t. 3D coords
        dNdx_global = dXdxi.T

        # Project 3D derivatives to local 2D plane
        dN_local = dNdx_global @ T_plane  # shape: (4, 2)
        dNdx = dN_local[:, 0]
        dNdy = dN_local[:, 1]

        B = compute_B_matrix(dNdx, dNdy)

        detJ = area_factor
        K += B.T @ D @ B * detJ * w * thickness

    return K


# Define coordinates of the 4-node membrane in 3D
coords = np.array([
    [0.5, 0.0, 0.0],
    [1.0, 1.0, 1.0],
    [0.5, 1.0, 1.0],
    [0.0, 0.5, 0.0]
])

# Material properties
E = 210e9      # Young's modulus (Pa)
nu = 0.3       # Poisson's ratio
thickness = 0.01  # meters

# Compute stiffness matrix
K_membrane = membrane_stiffness_Q4_3D(coords, thickness, E, nu)
K_membrane
