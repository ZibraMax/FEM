import FEM
import numpy as np
import matplotlib.pyplot as plt


a = 1
b = 1
L = (a**2 + b**2)**0.5
h = 0.5
t = 0.1

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
                [6, 0],
                [7, 0],
                [8, 0],
                [9, 0],
                [11, 0],
                [15, 0],
                [17, 0],
                [14, 1],
                [5, 1]]

E = 210e3  # Young's modulus in MPa
nu = 0.3  # Poisson's ratio

O = FEM.ElasticityMembranes(geometry, E, nu, t)

O.solve()
O.elementMatrices()
K = np.array(O.K.todense())
R = K@O.U
stff = R[5]+R[14]
k = 2*E*h*t/L*(b/L)**2
print(f"Stiffness: {stff}, Expected: {k}")
