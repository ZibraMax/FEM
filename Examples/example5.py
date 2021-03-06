import numpy as np
import matplotlib.pyplot as plt
from FEM.Torsion2D import Torsion2D
from FEM.Mesh.Delaunay import Delaunay

a = 1
b = 1
E = 200000
v = 0.2
G = E / (2 * (1 + v))
phi = 1
vertices = [[0, 0], [a, 0], [a, a], [0, a]]
params = Delaunay._strdelaunay(
    constrained=True, delaunay=True, a='0.0005', o=2)
geometria = Delaunay(vertices, params)
geometria.saveMesh('Mesh_tests/Square_torsion')
O = Torsion2D(geometria, G, phi)
O.solve()
plt.show()
