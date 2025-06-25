import FEM
import numpy as np
import matplotlib.pyplot as plt

coords = np.array([
    [0.5, 0.0, 0.0],
    [1.0, 1.0, 1.0],
    [0.5, 1.0, 1.0],
    [0.0, 0.5, 0.0]
])

coords = np.array([
    [0.0, 0.0, 0.0],  # Node 1
    [1.0, 0.0, 0.0],  # Node 2
    [1.0, 1.0, 0.0],  # Node 3
    [0.0, 1.0, 0.0]   # Node 4
])

coords = np.array([
    [0.0, 0.0, 0.0],  # Node 1
    [0.0, 1.0, 0.0],  # Node 2
    [0.0, 1.0, 1.0],  # Node 3
    [0.0, 0.0, 1.0]   # Node 4
])

gdl = np.array([[0, 1, 2],
                [3, 4, 5],
                [6, 7, 8],
                [9, 10, 11]]).T
m = gdl.shape[0]*gdl.shape[1]
E = 210e3
nu = 0.3

c11 = E / (1 - nu**2)  # Revisar
c12 = nu * c11  # Revisar
c22 = c11
c66 = E / (2 * (1 + nu))
C = np.array([
    [c11, c12, 0.0],
    [c12, c22, 0.0],
    [0.0, 0.0, c66]])
t = 0.05
e = FEM.Quadrilateral(coords, gdl, 2, boundary=True)
__x, __p = e.T(e.domain.T)
_x, _p = e.T(e.Z.T)
_j, _dp = e.J(e.Z.T)
weights = e.W
fig = plt.figure()
ax = fig.add_subplot(projection='3d')


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


def build_global_transform(T_plane, num_nodes):
    T_e = np.zeros((3 * num_nodes, 2 * num_nodes))
    for i in range(num_nodes):
        row_start = 3 * i
        col_start = 2 * i
        T_e[row_start:row_start+3, col_start:col_start +
            2] = T_plane  # insert 3×2 block

    return T_e


Ke = 0.0
# ax.plot_trisurf(*__x.T, alpha=0.3)
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
    T_plane = np.vstack((e1, e2)).T  # 3x2 transformation to local plane

    # En teoria estas son las derivadas de las funciones de forma pero en X, Y, Z
    dpx = JP @ dni
    # Ahora hay que proyectar esas derivadas en el plano usando la base ortonormal
    # Entonces ahoera estas derivadas estan en el plano de la menmbrana y se pueden usar para calcular rigidez
    dpt = dpx.T @ T_plane
    dpt = dpt.T  # Transpose porque si

    B = compute_B_matrix(dpt[0], dpt[1])
    T = build_global_transform(T_plane, len(e.gdl.T))
    Ke_local = t*(B.T@C@B)*detjac*wi
    Ke += T @ Ke_local @ T.T
    ax.plot(*x, 'o', c='k', zorder=200)
    ax.quiver(*x, *e3/3, color="green", zorder=200)
    ax.quiver(*x, *e1/3, color="blue", zorder=200)
    ax.quiver(*x, *e2/3, color="red", zorder=200)
plt.show()
