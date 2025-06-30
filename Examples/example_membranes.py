import FEM
import numpy as np
import matplotlib.pyplot as plt


def analytical_force(u, E, A, L, h):
    L0 = np.sqrt(L**2 + h**2)
    y = h - u
    l = np.sqrt(L**2 + y**2)
    strain = (l - L0) / L0
    N = E * A * strain
    vertical_component = (y / l)
    return 2 * N * vertical_component


a = 1
b = 1
L = (a**2 + b**2)**0.5
h = 0.5
t = 0.1
A = h*t

theta = np.arctan2(b, a)

coords = [
    [0, 0, 0],
    [a, 0, b],
    [2*a, 0, 0],
    [0, h, 0],
    [a, h, b],
    [2*a, h, 0],
]

conectivity = [
    [0, 1, 4, 3],
    [1, 2, 5, 4]]

geometry = FEM.Geometry3D(conectivity, coords, ["MEM"]*len(conectivity), 3)
geometry.cbe = [[0, 0],
                [1, 0],
                [2, 0],
                [3, 0],
                [4, 0],
                [6, 0],
                [7, 0],
                [8, 0],
                [9, 0],
                [11, 0],
                [12, 0],
                [15, 0],
                [17, 0]]
# [14, 1],
# [5, 1]]
P = 1000
geometry.cbn = [[14, -P], [5, -P]]
E = 210e3  # Young's modulus in MPa
nu = 0.3  # Poisson's ratio

O = FEM.ElasticityMembranes(geometry, E, nu, t, solver=FEM.MGDCM)
O.solver.set_increments(150)
O.solver.maxiter = 100
O.solver.momentum = False
O.solver.tol = 1e-1
O.solver.set_delta_lambda_bar(0.1)
O.solve()
O.elementMatrices()
K = np.array(O.K.todense())
R = K@O.U
stff = R[5]+R[14]
k = 2*E*h*t/L*(b/L)**2
print(f"Stiffness: {stff}, Expected: {k}")
O.exportJSON("membrane.json")
# Based on solutions create force displacement plot
displacements = []
load_factors = []
for i in range(len(O.solver.solutions)):
    O.solver.setSolution(i, elements=True)
    disp = [O.U[14], O.U[5], O.U[4], O.U[3]]

    displacements.append(disp)
    load_factors.append([O.solution_info['ld']])
displacements = np.array(displacements)
load_factors = np.array(load_factors)


us = np.linspace(0, -np.min(displacements), 300)
force = analytical_force(us, E, A, a, b)
plt.plot(us, -force/P/2,
         '-', c='k', lw=3, label="Analytical force")

plt.plot(-displacements[:, 0], load_factors,
         'o-', label="Displacement at node 14")
plt.plot(-displacements[:, 1], load_factors,
         'o-', label="Displacement at node 5")

plt.ylabel("Load factor")
plt.xlabel("Displacement")
plt.title("Force-Displacement Curve for Membrane Element")
plt.legend()
plt.grid()
plt.show()
