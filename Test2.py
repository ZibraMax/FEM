import numpy as np
import matplotlib.pyplot as plt
import FEM
from FEM import Mesh

a = 0.3
b = 0.3
tw = 0.05
tf = 0.05
E = 200000
v = 0.27
G = E/(2*(1+v))
phi = 1

vertices = [[0,0],[a,0],[a,tf],[a/2+tw/2,tf],[a/2+tw/2,tf+b],[a,tf+b],[a,2*tf+b],[0,2*tf+b],[0,tf+b],[a/2-tw/2,tf+b],[a/2-tw/2,tf],[0,tf]]
params = Mesh.Delaunay._strdelaunay(constrained=True,delaunay=True,a='0.00003',o=2)
geometria = Mesh.Delaunay1V(vertices, params)
# geometria.saveMesh('I_test')
# geometria = Mesh.Geometry.loadmsh('Mesh_tests/I_test.msh')
print(len(geometria.elements))
O = FEM.Torsion2D(geometria,G,phi)
O.solve()
plt.show()
