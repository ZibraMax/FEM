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
geometria = Mesh.Geometry.loadGiDMsh('test.msh')
seg = [[366,3053],[3053,3795],[3795,1626],[1626,2471],[2471,1122]]
geometria.segments = seg
O = FEM.Torsion2D(geometria,G,phi)
O.solve()
plt.show()
