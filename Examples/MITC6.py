import numpy as np
import matplotlib.pyplot as plt

# -------------------------
# Functions
# -------------------------


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
    """
    Builds the 30x30 rotation matrix T for 5-DOF shell element:
    - Translations (u,v,w) rotated with basis,
    - Rotations (theta_x, theta_y) remain local.

    Inputs:
    - basis: (3x3) matrix with e1, e2, e3 as columns.

    Output:
    - T: (30x30) global transformation matrix
    """
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
    """
    Assembles the tangent stiffness matrix and internal force vector
    for MITC6 5-DOF shell element (large displacement, updated geometry).

    Inputs:
    - x_current: (6x3) array of current nodal coordinates (deformed shape)
    - u_disp: (30,) displacement vector
    - D_membrane: (3x3) material matrix (membrane)
    - D_bending: (3x3) material matrix (bending)
    - D_shear: (2x2) material matrix (shear)

    Outputs:
    - K_tangent: (30x30) tangent stiffness matrix
    - F_int: (30,) internal force vector
    """

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

    return K_tangent, F_int


def generate_flat_panel():
    """
    Generate a simple flat triangular panel (ground plane) with 6-node MITC6 layout.
    """
    X = np.zeros((6, 3))

    # Corner nodes
    X[0, :] = [0.0, 0.0, 0.0]  # Node 1
    X[1, :] = [1.0, 0.0, 0.0]  # Node 2
    X[2, :] = [0.0, 1.0, 0.0]  # Node 3

    # Mid-side nodes (midpoints)
    X[3, :] = [(0.0 + 1.0)/2, (0.0 + 0.0)/2, 0.0]  # Node 4 (between Node 1-2)
    X[4, :] = [(1.0 + 0.0)/2, (0.0 + 1.0)/2, 0.0]  # Node 5 (between Node 2-3)
    X[5, :] = [(0.0 + 0.0)/2, (1.0 + 0.0)/2, 0.0]  # Node 6 (between Node 3-1)

    return X


# --- Geometry: Quarter-cylinder Shell ---
R = 1.0  # Radius 1 meter

X = generate_flat_panel()*10

n_nodes = 6
n_dof = 5
total_dofs = n_nodes * n_dof

# Initialize displacements
u = np.zeros(total_dofs)

# --- Material properties ---

D_membrane, D_bending, D_shear = get_D_matrices(
    E=210000, nu=0.3, t=1, kappa=5/6)

pressure_total = -0.01  # Negative = pressure downward (N/m^2)

# --- Newton loop ---
for i in range(5):
    # --- External Force: Uniform pressure ---
    pressure_total -= 10  # Negative = pressure downward (N/m^2)

    # Approximate element area (in undeformed shape) for one triangle
    area = 0.5  # m² (approximate for small element)

    # Distribute pressure load to nodes
    F_ext = np.zeros(total_dofs)
    F_ext[12] = pressure_total
    for i in range(n_nodes):
        F_ext[i*n_dof + 0] = 100.0  # Apply 100 N in x direction (tension)

    # --- Boundary Conditions: Clamp node 1 (full fixity) ---
    # u, v, w, theta_x, theta_y at node 1
    fixed_dofs = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
    free_dofs = np.setdiff1d(np.arange(total_dofs), fixed_dofs)
    # --- Newton-Raphson Parameters ---
    max_iter = 30
    tolerance = 1e-3
    for iteration in range(max_iter):
        print(f"--- Iteration {iteration} ---")

        # Only translation part for geometry update
        x_current = X  # + u.reshape((n_nodes, n_dof))[:, 0:3]
        # Assemble tangent stiffness matrix and internal force vector
        K_tangent, F_int = assemble_Ktangent_Fint(
            x_current, u, D_membrane, D_bending, D_shear)

        # Residual
        R = F_ext - F_int

        res_norm = np.linalg.norm(R[free_dofs])
        print(f"Residual norm = {res_norm:.4e}")

        if res_norm < tolerance:
            print("Converged successfully!")
            break

        # Solve for displacement increment
        K_ff = K_tangent[np.ix_(free_dofs, free_dofs)]
        R_f = R[free_dofs]

        delta_u_f = np.linalg.solve(K_ff, R_f)

        # Update displacements
        u[free_dofs] += delta_u_f

    else:
        print("Did not converge within maximum iterations.")
# --- Plot undeformed and deformed shapes ---

scale = 100  # Amplify deformation for visualization
X_deformed = X.copy()
for i in range(n_nodes):
    X_deformed[i, 0] += u[i*n_dof + 0] * scale
    X_deformed[i, 1] += u[i*n_dof + 1] * scale
    X_deformed[i, 2] += u[i*n_dof + 2] * scale

print(u)
fig = plt.figure(figsize=(10, 6))
ax = fig.add_subplot(111, projection='3d')

# Undeformed
ax.plot_trisurf(X[:, 0], X[:, 1], X[:, 2], color='cyan',
                alpha=0.5, label='Undeformed')
ax.plot(X[:, 0], X[:, 1],
        X[:, 2], "o", c='b', label='Deformed')
# Deformed
ax.plot_trisurf(X_deformed[:, 0], X_deformed[:, 1],
                X_deformed[:, 2], color='red', alpha=0.3, label='Deformed')
ax.plot(X_deformed[:, 0], X_deformed[:, 1],
        X_deformed[:, 2], 'o', c="r", label='Deformed')
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
ax.set_title('Curved MITC6 Element: Deformation')
plt.show()
