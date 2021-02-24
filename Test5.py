import numpy as np
import matplotlib.pyplot as plt
import FEM
from FEM import Mesh

a = 1
b = 1
E = 200000
v = 0.27
G = E/(2*(1+v))*0+1
phi =1

vertices = [[0,0],[a,0],[a,a],[0,a]]
params = Mesh.Delaunay._strdelaunay(constrained=True,delaunay=True,a='0.0003',o=2)
geometria = Mesh.Delaunay1V(vertices, params)
geometria.saveMesh('Mesh_tests/Square_torsion')
geometria.show()
plt.show()
O = FEM.Torsion2D(geometria,G,phi)
O.solve()
plt.show()
