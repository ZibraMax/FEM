import numpy as np
import matplotlib.pyplot as plt
import FEM
from FEM import Mesh

a = 5
b = 5
E = 200000
v = 0.2

vertices = [[0,0],[a,0],[a,a],[0,a]]
params = Mesh.Delaunay._strdelaunay(constrained=True,delaunay=True,a='0.0005',o=2)
geometria = Mesh.Delaunay1V(vertices, params)
geometria.saveMesh('Mesh_tests/Square_torsion')
O = FEM.Torsion2D(geometria,G,phi)
O.solve()
plt.show()
