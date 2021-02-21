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
phi =1

vertices = [[0,0],[a,0],[a,tf],[a/2+tw/2,tf],[a/2+tw/2,tf+b],[a,tf+b],[a,2*tf+b],[0,2*tf+b],[0,tf+b],[a/2-tw/2,tf+b],[a/2-tw/2,tf],[0,tf]]
params = Mesh.Delaunay._strdelaunay(constrained=True,delaunay=True,a='0.0005',o=1)
geometria = Mesh.Delaunay1V(vertices, params)
print(len(geometria.elements))
O = FEM.Torsion2D(geometria,G,phi)
print(O.solve())
print(O.elements[0].Ke)





# geometria.dibujarse()
plt.show()